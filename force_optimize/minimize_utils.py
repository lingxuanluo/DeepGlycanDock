import os
from openff.toolkit import Molecule
from openmmforcefields.generators import SystemGenerator
from openmm import unit, LangevinIntegrator
from openmm.app import PDBFile, Simulation, element
from pdbfixer import PDBFixer
import traceback
from rdkit import Chem
from rdkit.Chem import AllChem
import torch
import numpy as np
from rdkit import Chem
import warnings
from openmm import unit, Platform, State
from joblib import wrap_non_picklable_objects
from joblib import delayed
import re
from openmm.app import Modeller
import sys
import loguru
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from cleaup import clean_structure,fix_pdb
from openmm.app.internal.pdbstructure import PdbStructure
import io
import subprocess
from loguru import logger
def run_command(command: str, cwd_path: str) -> None:
    r"""
    Create a child process and run the command in the cwd_path.
    It is more safe than os.system.
    """
    proc = subprocess.Popen(
        command,
        shell=True,
        cwd=cwd_path,
        executable="/bin/bash",
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    errorcode = proc.wait()
    if errorcode:
        path = cwd_path
        msg = (
            'Failed with command "{}" failed in '
            ""
            "{} with error code {}"
            "stdout: {}"
            "stderr: {}".format(command, path, errorcode, proc.stdout.read().decode(), proc.stderr.read().decode())
        )
        raise ValueError(msg)

def read_molecule(molecule_file, sanitize=False, calc_charges=False, remove_hs=False):
    if molecule_file.endswith('.mol2'):
        mol = Chem.MolFromMol2File(molecule_file, sanitize=False, removeHs=False)
    elif molecule_file.endswith('.sdf'):
        supplier = Chem.SDMolSupplier(molecule_file, sanitize=False, removeHs=False)
        mol = supplier[0]
    elif molecule_file.endswith('.pdbqt'):
        with open(molecule_file) as file:
            pdbqt_data = file.readlines()
        pdb_block = ''
        for line in pdbqt_data:
            pdb_block += '{}\n'.format(line[:66])
        mol = Chem.MolFromPDBBlock(pdb_block, sanitize=False, removeHs=False)
    elif molecule_file.endswith('.pdb'):
        mol = Chem.MolFromPDBFile(molecule_file, sanitize=False, removeHs=False)
    else:
        raise ValueError('Expect the format of the molecule_file to be '
                         'one of .mol2, .sdf, .pdbqt and .pdb, got {}'.format(molecule_file))

    try:
        if sanitize or calc_charges:
            Chem.SanitizeMol(mol)

        if calc_charges:
            # Compute Gasteiger charges on the molecule.
            try:
                AllChem.ComputeGasteigerCharges(mol)
            except:
                warnings.warn('Unable to compute charges for the molecule.')
        if remove_hs:
            mol = Chem.RemoveHs(mol, sanitize=sanitize)
    except Exception as e:
        logger.info(e)
        logger.info("RDKit was unable to read the molecule.")
        return None

    return mol
def read_abs_file_mol(file, remove_hs=False, sanitize=True):
    mol = read_molecule(file, remove_hs=remove_hs, sanitize=True)
    
    if file.endswith(".sdf") and mol is None:
        # mol = read_molecule(file, remove_hs=remove_hs, sanitize=True)
        if os.path.exists(file[:-4] + ".mol2"):
            logger.info('Using the .sdf file failed. We found a .mol2 file instead and are trying to use that.')
            mol = read_molecule(file[:-4] + ".mol2", remove_hs=remove_hs, sanitize=True)
    elif file.endswith(".mol2") and mol is None:
        if os.path.exists(file[:-4] + ".sdf"):
            logger.info('Using the .mol2 file failed. We found a .sdf file instead and are trying to use that.')
            mol = read_molecule(file[:-4] + ".sdf", remove_hs=remove_hs, sanitize=True)

    return mol
# from joblib.externals.loky import set_loky_pickler
def trySystem(system_generator,modeller,ligand_mol,lig_path):
   
    max_attempts = 100
    attempts = 0
    success = False
    while attempts < max_attempts and not success:
        try:
            system = system_generator.create_system(modeller.topology, molecules=ligand_mol)
            success = True  # Mark
        except Exception as e:
            # extract the error residue index from the error message
            logger.info(f'Try DELETE THIS ERROE {str(e)}!')
            match = re.search(r"residue (\d+)", str(e))
            if match:
                extracted_index = int(match.group(1)) - 1
                # located and record the residue to delete
                current_index = 0
                residue_to_delete = None
                for residue in modeller.topology.residues():
                    if current_index == extracted_index:
                        residue_to_delete = residue
                        break
                    current_index += 1

                modeller.delete([residue_to_delete])
        finally:
            attempts += 1
    if not success:
        logger.info("Try maximum times but cannot create system")
        return None
    else:
        logger.info(f"Try {attempts} times and system is created successfully")
        with open(os.path.join('/home/house/caoduanhua_tmp/DeepLearningForDock/DiffDockForScreen/diffScreen/Screen_dataset/create_system_pdbs',os.path.basename(lig_path).split('_')[0]+'_create_system.pdb'), "w") as f:
            PDBFile.writeFile(modeller.topology, modeller.positions, f)

    return modeller

def UpdatePose(lig_path,system_generator,modeller,protein_atoms,out_dir,device_num=0):
    try:
        # init save path
        out_base_dir = os.path.join(out_dir,lig_path.split('/')[-2])
        os.makedirs(out_base_dir,exist_ok=True)
        out_file = os.path.join(out_base_dir,os.path.splitext(os.path.basename(lig_path))[0] + '_minimized.sdf')
        if os.path.exists(out_file):
            return 0
        dockingpose = read_molecule(lig_path, remove_hs=True, sanitize=True)
        lig_mol = Molecule.from_rdkit(dockingpose,allow_undefined_stereo=True)
        lig_mol.assign_partial_charges(partial_charge_method='gasteiger')

        lig_top = lig_mol.to_topology()
        modeller.add(lig_top.to_openmm(), lig_top.get_positions().to_openmm())
        # create simulation system
        system=system_generator.create_system(modeller.topology,molecules=lig_mol)
        # keep protein atom static in smiulation 
        for atom in protein_atoms:
            system.setParticleMass(atom.index, 0.000*unit.dalton)
        # start simulation
        platform = GetPlatform()
        simulation = EnergyMinimized(modeller,system, platform,verbose=False,device_num=device_num)
        # get energy minimized conformer and modify the graph['ligand'].pos to scoring
        # use conformer mapping
        ligand_atoms = list(filter(lambda atom:  atom.residue.name == 'UNK',list(modeller.topology.atoms())))
        ligand_index = [atom.index for atom in ligand_atoms]
        new_coords = simulation.context.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.angstrom)[ligand_index]
        lig_mol = lig_mol.to_rdkit()
        conf = lig_mol.GetConformer()
        for i in range(lig_mol.GetNumAtoms()):
            x,y,z = new_coords.astype(np.double)[i]
            conf.SetAtomPosition(i,Point3D(x,y,z))
        try:
            writer = Chem.SDWriter(out_file)
            writer.write(lig_mol)
            writer.close()
        except:
            out_base_dir = os.path.join(out_dir,lig_path.split('/')[-2] + '_tmp')
            os.makedirs(out_base_dir,exist_ok=True)
            out_file = os.path.join(out_base_dir,os.path.splitext(os.path.basename(lig_path))[0] + '_minimized.sdf')
            if os.path.exists(out_file):
                return 0
            writer = Chem.SDWriter(out_file)
            writer.write(lig_mol)
            writer.close()
        return 0
        # return lig_mol
    except Exception as e:
        error_info = traceback.format_exc()
        logger.info(error_info)
        logger.warning(f' : {e}')
        with open('error_sdf.txt','a') as f:
            f.write(lig_path +': error by :' + error_info + '\n')
        return 1

def UpdateGrpah(graph,system_generator,modeller,protein_atoms,device_num=0):
    try:
        # raw_position = graph['ligand'].pos
        dockingpose = GetDockingPose(graph)
        lig_mol = Molecule.from_rdkit(dockingpose,allow_undefined_stereo=True)
        lig_mol.assign_partial_charges(partial_charge_method='gasteiger')
        # add ligand to modeller
        lig_top = lig_mol.to_topology()
        modeller.add(lig_top.to_openmm(), lig_top.get_positions().to_openmm())
        # create simulation system
        platform = GetPlatform()

        system = system_generator.create_system(modeller.topology,molecules=lig_mol)
        # keep protein atom static in smiulation 
        for atom in protein_atoms:
            system.setParticleMass(atom.index, 0.000*unit.dalton)
        # start simulation
        simulation = EnergyMinimized(modeller,system, platform,verbose=False,device_num=device_num)
        # get energy minimized conformer and modify the graph['ligand'].pos to scoring
        # conformer mapping
        ligand_atoms = list(filter(lambda atom:  atom.residue.name == 'UNK',list(modeller.topology.atoms())))
        ligand_index = [atom.index for atom in ligand_atoms]
        new_coords = simulation.context.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.angstrom)[ligand_index]

        new_coords -= graph.original_center.detach().cpu().numpy()
        lig_mol = lig_mol.to_rdkit()
        conf = lig_mol.GetConformer()
        for i in range(lig_mol.GetNumAtoms()):
            x,y,z = new_coords.astype(np.double)[i]
            conf.SetAtomPosition(i,Point3D(x,y,z))
        lig_mol = Chem.RemoveHs(lig_mol)
      
        graph['ligand'].pos = torch.from_numpy(lig_mol.GetConformer().GetPositions()).to(graph.original_center.device).float()
    
        return graph
    except Exception as e:
        error_info = traceback.format_exc()
        logger.info(error_info)
        warnings.warn(graph['name'][0]+f' : {e}')
        return 1

def DescribeState(state: State, name: str):
    """logger.info energy and force information about a simulation state."""
    max_force = max(np.linalg.norm([v.x, v.y, v.z]) for v in state.getForces())
    logger.info(f"{name} has energy {state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole):.2f} kJ/mol "
          f"with maximum force {max_force:.2f} kJ/(mol nm)")
def GetFFGenerator(protein_forcefield = 'amber/ff14SB.xml',water_forcefield = 'amber/tip3p_standard.xml',small_molecule_forcefield = 'openff-2.0.0',ignoreExternalBonds=False):
    """
    Get forcefield generator by different forcefield files
    """
    forcefield_kwargs = {'constraints': None, 'rigidWater': True, 'removeCMMotion': False, 'ignoreExternalBonds': ignoreExternalBonds, 'hydrogenMass': 4*unit.amu }
    # forcefield_kwargs = {'constraints': None, 'rigidWater': True, 'removeCMMotion': False, 'hydrogenMass': 4*unit.amu }
    system_generator = SystemGenerator(
                forcefields=[protein_forcefield, water_forcefield ],
                small_molecule_forcefield=small_molecule_forcefield,
                forcefield_kwargs=forcefield_kwargs)
    return system_generator
def GetfixedPDB(receptor_path):
    
    temp_fixd_pdbs = f'{os.path.dirname(os.path.dirname(os.path.abspath(__file__)))}/fixed_pdbs'
    os.makedirs(temp_fixd_pdbs,exist_ok=True)
    if not os.path.exists(os.path.join(temp_fixd_pdbs,os.path.basename(receptor_path).replace('.pdb','_fixer_processed_cleanup.pdb'))):
        alterations_info  = {}
        fixed_pdb = fix_pdb(receptor_path, alterations_info)
        fixed_pdb_file = io.StringIO(fixed_pdb)
        pdb_structure = PdbStructure(fixed_pdb_file)
        clean_structure(pdb_structure, alterations_info)
        fixer = PDBFile(pdb_structure)
        logger.info("Protein loaded with success!")
        PDBFile.writeFile(fixer.topology, fixer.positions, open(os.path.join(temp_fixd_pdbs,os.path.basename(receptor_path).replace('.pdb','_fixer_processed_cleanup.pdb')), 'w'))
        logger.info('Dont have processed by fixer try fix and save in disk')
    else:
        fixer = PDBFixer(os.path.join(temp_fixd_pdbs,os.path.basename(receptor_path).replace('.pdb','_fixer_processed_cleanup.pdb')))
        logger.info('There have a precessed pdb file use it!')
    return fixer

import copy
from rdkit.Geometry import Point3D
def GetDockingPose(graph):
    mol = copy.deepcopy(graph.mol[0] if type(graph.mol) == list else graph.mol)
    mol = Chem.RemoveHs(mol)
    docking_position = graph['ligand'].pos.detach().cpu().numpy() # without Hs and dont match with raw pocket
    docking_position = docking_position + graph.original_center.detach().cpu().numpy()
    conf = mol.GetConformer()
    for i in range(mol.GetNumAtoms()):
        x,y,z = docking_position.astype(np.double)[i]
        conf.SetAtomPosition(i,Point3D(x,y,z))
    return mol

@delayed
@wrap_non_picklable_objects
def GetPlatformPara():
    """Determine the best simulation platform available."""
    platform_name = os.getenv('PLATFORM')
    # properties = {'CudaDeviceIndex': '0'}
    if platform_name:
        platform = Platform.getPlatformByName(platform_name)
    else:
        platform = max((Platform.getPlatform(i) for i in range(Platform.getNumPlatforms())), key=lambda x: x.getSpeed())
    logger.info(f'Using platform {platform.getName()}')
    if platform.getName() in ['CUDA', 'OpenCL']:
        platform.setPropertyDefaultValue('Precision', 'mixed')
        logger.info(f'Set precision for platform {platform.getName()} to mixed')
    return platform
# @delayed
# @wrap_non_picklable_objects
def GetPlatform():
    """Determine the best simulation platform available."""
    platform_name = os.getenv('PLATFORM')
    # properties = {'CudaDeviceIndex': '0'}
    if platform_name:
        platform = Platform.getPlatformByName(platform_name)
    else:
        platform = max((Platform.getPlatform(i) for i in range(Platform.getNumPlatforms())), key=lambda x: x.getSpeed())
    logger.info(f'Using platform {platform.getName()}')
    if platform.getName() in ['CUDA', 'OpenCL']:
        platform.setPropertyDefaultValue('Precision', 'mixed')
        logger.info(f'Set precision for platform {platform.getName()} to mixed')
    return platform

def compute_distance(pos1, pos2, box_vectors=None):
    diff = pos1 - pos2
    if box_vectors:  # Handle periodic boundary conditions
        box = np.array([v.value_in_unit(unit.nanometers) for v in box_vectors])
        diff_nm = diff.value_in_unit(unit.nanometers)
        diff_nm -= np.round(diff_nm / box.diagonal()) * box.diagonal()
        diff = diff_nm * unit.nanometers
    return np.sqrt(np.sum(diff.value_in_unit(unit.nanometers)**2)) * unit.nanometers

# 计算残基的平均坐标（几何中心），忽略氢原子
def compute_residue_average_position(residue, positions):
    num_atoms = 0
    avg_pos = np.zeros(3) * unit.nanometers
    included_atoms = []
    excluded_hydrogens = []
    
    for atom in residue.atoms():
        if atom.element == element.hydrogen:
            excluded_hydrogens.append(atom.index)
            continue  # 忽略氢原子
        pos = positions[atom.index]  # Use Quantity directly
        avg_pos += pos
        num_atoms += 1
        included_atoms.append(atom.index)
    
    if num_atoms == 0:
        raise ValueError(f"Residue {residue.name} (ID: {residue.id}) has no non-hydrogen atoms.")
    
    # if excluded_hydrogens:
    #     print(f"Residue {residue.name} (ID: {residue.id}): Excluded {len(excluded_hydrogens)} hydrogen atoms: {excluded_hydrogens}")
    # print(f"Residue {residue.name} (ID: {residue.id}): Included {num_atoms} non-hydrogen atoms: {included_atoms}")
    
    return avg_pos / num_atoms

# 检查 UNK 残基与系统中其他原子的最小距离，忽略氢原子
def check_min_distance(unk_residue, topology, positions, min_distance=3.0 * unit.angstroms, box_vectors=None):
    unk_atoms = [atom for atom in unk_residue.atoms() if atom.element != element.hydrogen]
    if not unk_atoms:
        raise ValueError(f"UNK residue (ID: {unk_residue.id}) has no non-hydrogen atoms.")
    
    min_dist = float('inf') * unit.nanometers
    for residue in topology.residues():
        if residue.name != 'UNK':
            for atom in residue.atoms():
                if atom.element == element.hydrogen:
                    continue  # 忽略氢原子
                for unk_atom in unk_atoms:
                    dist = compute_distance(
                        positions[unk_atom.index],
                        positions[atom.index],
                        box_vectors
                    )
                    min_dist = min(min_dist, dist)
    return min_dist

# 移动 UNK 残基到远离其他残基平均坐标的方向
def move_unk_residue(modeller, system, min_distance=3.0 * unit.angstroms, step_size=0.1 * unit.angstroms, max_steps=1000):
    topology = modeller.topology
    positions = modeller.positions
    box_vectors = topology.getPeriodicBoxVectors() if topology.getPeriodicBoxVectors() else None

    # 找到 UNK 残基
    unk_residue = None
    for residue in topology.residues():
        if residue.name == 'UNK':
            unk_residue = residue
            break
    if unk_residue is None:
        raise ValueError("No residue named 'UNK' found in the topology.")

    # 调试：检查 UNK 残基的原子
    # print(f"Inspecting UNK residue (ID: {unk_residue.id})")
    for atom in unk_residue.atoms():
        is_hydrogen = atom.element == element.hydrogen
        pos = positions[atom.index]
        # print(f"  Atom {atom.index}: {atom.name}, Position: {pos.value_in_unit(unit.nanometers)} nm, Hydrogen: {is_hydrogen}")

    # 计算其他残基的平均坐标
    total_atoms = 0
    other_avg_pos = np.zeros(3) * unit.nanometers
    for residue in topology.residues():
        if residue.name != 'UNK':
            try:
                res_avg_pos = compute_residue_average_position(residue, positions)
                num_non_hydrogen_atoms = sum(1 for atom in residue.atoms() if atom.element != element.hydrogen)
                other_avg_pos += res_avg_pos * num_non_hydrogen_atoms
                total_atoms += num_non_hydrogen_atoms
            except ValueError as e:
                print(f"Skipping residue {residue.name} (ID: {residue.id}): {e}")
    
    if total_atoms == 0:
        raise ValueError("No valid residues with non-hydrogen atoms found to compute average position.")
    
    other_avg_pos /= total_atoms

    # 计算 UNK 残基的平均坐标
    unk_avg_pos = compute_residue_average_position(unk_residue, positions)

    # 计算移动方向（远离其他残基平均坐标）
    direction = (unk_avg_pos - other_avg_pos)
    norm = np.sqrt(np.sum(direction.value_in_unit(unit.nanometers)**2))
    if norm == 0:
        raise ValueError("Cannot compute direction: UNK and other residues have the same average position.")
    direction = direction / norm  # 归一化

    # 逐步移动 UNK 残基
    step = 0
    while step < max_steps:
        # 检查当前最小距离
        current_min_dist = check_min_distance(unk_residue, topology, positions, min_distance, box_vectors)
        if current_min_dist > min_distance and step != 0:
            print(f"Minimum distance achieved: {current_min_dist.value_in_unit(unit.angstroms):.2f} Å after {step} steps.")
            break

        # 移动 UNK 残基的所有原子（包括氢原子）
        for atom in unk_residue.atoms():
            pos = positions[atom.index]
            new_pos = pos + step_size * direction / unit.angstroms / unit.angstroms * unit.nanometers
            positions[atom.index] = new_pos

        step += 1

    if step == max_steps:
        print(f"Warning: Reached maximum steps ({max_steps}) without achieving minimum distance of {min_distance}.")

    # 更新 modeller 的位置
    modeller.positions = positions
    return modeller

def EnergyMinimized(modeller,system, platform,verbose=False,device_num = 0):
    modeller = move_unk_residue(modeller, system, min_distance=3.0 * unit.angstroms, step_size=0.01 * unit.angstroms)
    integrator = LangevinIntegrator(
    300 * unit.kelvin,
    1 / unit.picosecond,
    0.002 * unit.picoseconds,
    )
    properties = {'CudaDeviceIndex': f'{device_num}'}
    simulation = Simulation(modeller.topology, system = system, integrator = integrator, platform=platform,platformProperties=properties)
    simulation.context.setPositions(modeller.positions)
    if verbose:
        DescribeState(
            simulation.context.getState(
                getEnergy=True,
                getForces=True,
            ),
            "Original state",
        )


    # simulation.minimizeEnergy()
    if verbose:
        DescribeState(
            simulation.context.getState(
                getEnergy=True, 
                getForces=True),
            "Minimized state",
        )
    return simulation
