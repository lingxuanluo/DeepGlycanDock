import os
import subprocess
import warnings
from datetime import datetime
import signal
from contextlib import contextmanager
import numpy as np
import torch
import yaml
from rdkit import Chem
from rdkit.Chem import RemoveHs, MolToPDBFile
from torch_geometric.nn.data_parallel import DataParallel

from models.surface_score_model_glycan_pretrain import TensorProductScoreModel as SurfaceScoreModelV3 

from models.mdn_score_model_v6 import TensorProductScoreModelV6 as ConfidenceCGScoreModelV6
# from models.score_model_mdn_energy_v1 import TensorProductEnergyModel
from utils_glycan.diffusion_utils import get_timestep_embedding
from spyrmsd import rmsd, molecule
from Bio import PDB
from collections import defaultdict
def get_obrmsd(mol1_path, mol2_path, cache_name=None):
    cache_name = datetime.now().strftime('date%d-%m_time%H-%M-%S.%f') if cache_name is None else cache_name
    os.makedirs(".openbabel_cache", exist_ok=True)
    if not isinstance(mol1_path, str):
        MolToPDBFile(mol1_path, '.openbabel_cache/obrmsd_mol1_cache.pdb')
        mol1_path = '.openbabel_cache/obrmsd_mol1_cache.pdb'
    if not isinstance(mol2_path, str):
        MolToPDBFile(mol2_path, '.openbabel_cache/obrmsd_mol2_cache.pdb')
        mol2_path = '.openbabel_cache/obrmsd_mol2_cache.pdb'
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return_code = subprocess.run(f"obrms {mol1_path} {mol2_path} > .openbabel_cache/obrmsd_{cache_name}.rmsd",
                                     shell=True)
        print(return_code)
    obrms_output = read_strings_from_txt(f".openbabel_cache/obrmsd_{cache_name}.rmsd")
    rmsds = [line.split(" ")[-1] for line in obrms_output]
    return np.array(rmsds, dtype=np.float)


def remove_all_hs(mol,santize=None):
    params = Chem.RemoveHsParameters()
    params.removeAndTrackIsotopes = True
    params.removeDefiningBondStereo = True
    params.removeDegreeZero = True
    params.removeDummyNeighbors = True
    params.removeHigherDegrees = True
    params.removeHydrides = True
    params.removeInSGroups = True
    params.removeIsotopes = True
    params.removeMapped = True
    params.removeNonimplicit = True
    params.removeOnlyHNeighbors = True
    params.removeWithQuery = True
    params.removeWithWedgedBond = True
    if santize is not None:
        params.sanitize = santize
    return RemoveHs(mol, params)


def read_strings_from_txt(path):
    # every line will be one element of the returned list
    with open(path) as file:
        lines = file.readlines()
        return [line.rstrip() for line in lines]


def save_yaml_file(path, content):
    assert isinstance(path, str), f'path must be a string, got {path} which is a {type(path)}'
    content = yaml.dump(data=content)
    if '/' in path and os.path.dirname(path) and not os.path.exists(os.path.dirname(path)):
        os.makedirs(os.path.dirname(path),exist_ok=True)
    with open(path, 'w') as f:
        f.write(content)

# from accelerate.utils import DummyOptim,DummyScheduler
def get_optimizer_and_scheduler(args, model, accelerator,scheduler_mode='min'):
    optimizer_cls = (
        torch.optim.AdamW
        if accelerator.state.deepspeed_plugin is None
        or "optimizer" not in accelerator.state.deepspeed_plugin.deepspeed_config
        else None #DummyOptim
    )
    optimizer = optimizer_cls(filter(lambda p: p.requires_grad, model.parameters()), lr=args.lr, weight_decay=args.w_decay)
    
    if args.scheduler == 'plateau':
        # scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode=scheduler_mode, factor=0.7,
                                                            #    patience=args.scheduler_patience, min_lr=args.lr / 100)
        
        if (
            accelerator.state.deepspeed_plugin is None
            or "scheduler" not in accelerator.state.deepspeed_plugin.deepspeed_config
        ):
            scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode=scheduler_mode, factor=0.7,
                                                               patience=args.scheduler_patience, min_lr=args.lr / 100)
        # else:
        #     lr_scheduler = DummyScheduler(
        #         optimizer, total_num_steps=args.max_train_steps, warmup_num_steps=args.num_warmup_steps
        #     )
    else:
        print('No scheduler')
        scheduler = None
    return optimizer, scheduler

def get_model(args, device, t_to_sigma, no_parallel=False, model_type='score_model'):
    # ['score_model','mdn_model','energy_score_model']
    timestep_emb_func = get_timestep_embedding(
        embedding_type=args.embedding_type,
        embedding_dim=args.sigma_embed_dim,
        embedding_scale=args.embedding_scale)
    lm_embedding_type = None
    if args.esm_embeddings_path is not None: lm_embedding_type = 'esm'
    if model_type == 'mdn_model':
        model_class = ConfidenceCGScoreModelV6
        model = model_class(args,t_to_sigma=t_to_sigma,
                    device=device,
                    no_torsion=args.no_torsion,
                    timestep_emb_func=timestep_emb_func,
                    num_conv_layers=args.num_conv_layers,
                    lig_max_radius=args.max_radius,
                    scale_by_sigma=args.scale_by_sigma,
                    sigma_embed_dim=args.sigma_embed_dim,
                    ns=args.ns, nv=args.nv,
                    distance_embed_dim=args.distance_embed_dim,
                    cross_distance_embed_dim=args.cross_distance_embed_dim,
                    batch_norm=not args.no_batch_norm,
                    dropout=args.dropout,
                    use_second_order_repr=args.use_second_order_repr,
                    cross_max_distance=args.cross_max_distance,
                    dynamic_max_cross=args.dynamic_max_cross,
                    lm_embedding_type=lm_embedding_type,
                    mdn_dropout=args.mdn_dropout,n_gaussians = args.n_gaussians)
    
    elif model_type == 'surface_score_model':

        model_class = SurfaceScoreModelV3

        model = model_class(t_to_sigma=t_to_sigma,
                    device=device,
                    no_torsion=args.no_torsion,
                    timestep_emb_func=timestep_emb_func,
                    num_conv_layers=args.num_conv_layers,
                    lig_max_radius=args.max_radius,
                    scale_by_sigma=args.scale_by_sigma,
                    sigma_embed_dim=args.sigma_embed_dim,
                    ns=args.ns, nv=args.nv,
                    distance_embed_dim=args.distance_embed_dim,
                    cross_distance_embed_dim=args.cross_distance_embed_dim,
                    batch_norm=not args.no_batch_norm,
                    dropout=args.dropout,
                    use_second_order_repr=args.use_second_order_repr,
                    cross_max_distance=args.cross_max_distance,
                    dynamic_max_cross=args.dynamic_max_cross,
                    lm_embedding_type=lm_embedding_type,
                   )
    else:
        raise f'not support {model_type} type model setup'

    model.to(device)
    return model


def get_symmetry_rmsd(mol, coords1, coords2, mol2=None):
    with time_limit(10):
        mol = molecule.Molecule.from_rdkit(mol)
        mol2 = molecule.Molecule.from_rdkit(mol2) if mol2 is not None else mol2
        mol2_atomicnums = mol2.atomicnums if mol2 is not None else mol.atomicnums
        mol2_adjacency_matrix = mol2.adjacency_matrix if mol2 is not None else mol.adjacency_matrix
        RMSD = rmsd.symmrmsd(
            coords1,
            coords2,
            mol.atomicnums,
            mol2_atomicnums,
            mol.adjacency_matrix,
            mol2_adjacency_matrix,
        )
        return RMSD


class TimeoutException(Exception): pass


@contextmanager
def time_limit(seconds):
    def signal_handler(signum, frame):
        raise TimeoutException("Timed out!")

    signal.signal(signal.SIGALRM, signal_handler)
    signal.alarm(seconds)
    try:
        yield
    finally:
        signal.alarm(0)


class ExponentialMovingAverage:
    """ from https://github.com/yang-song/score_sde_pytorch/blob/main/models/ema.py
    Maintains (exponential) moving average of a set of parameters. """

    def __init__(self, parameters, decay, use_num_updates=True):
        """
        Args:
          parameters: Iterable of `torch.nn.Parameter`; usually the result of
            `model.parameters()`.
          decay: The exponential decay.
          use_num_updates: Whether to use number of updates when computing
            averages.
        """
        if decay < 0.0 or decay > 1.0:
            raise ValueError('Decay must be between 0 and 1')
        self.decay = decay
        self.num_updates = 0 if use_num_updates else None
        self.shadow_params = [p.clone().detach()
                              for p in parameters if p.requires_grad]
        self.collected_params = []

    def update(self, parameters):
        """
        Update currently maintained parameters.
        Call this every time the parameters are updated, such as the result of
        the `optimizer.step()` call.
        Args:
          parameters: Iterable of `torch.nn.Parameter`; usually the same set of
            parameters used to initialize this object.
        """
        decay = self.decay
        if self.num_updates is not None:
            self.num_updates += 1
            decay = min(decay, (1 + self.num_updates) / (10 + self.num_updates))
        one_minus_decay = 1.0 - decay
        with torch.no_grad():
            parameters = [p for p in parameters if p.requires_grad]
            for s_param, param in zip(self.shadow_params, parameters):
                s_param.sub_(one_minus_decay * (s_param - param))

    def copy_to(self, parameters):
        """
        Copy current parameters into given collection of parameters.
        Args:
          parameters: Iterable of `torch.nn.Parameter`; the parameters to be
            updated with the stored moving averages.
        """
        parameters = [p for p in parameters if p.requires_grad]
        for s_param, param in zip(self.shadow_params, parameters):
            if param.requires_grad:
                param.data.copy_(s_param.data)

    def store(self, parameters):
        """
        Save the current parameters for restoring later.
        Args:
          parameters: Iterable of `torch.nn.Parameter`; the parameters to be
            temporarily stored.
        """
        self.collected_params = [param.clone() for param in parameters]

    def restore(self, parameters):
        """
        Restore the parameters stored with the `store` method.
        Useful to validate the model with EMA parameters without affecting the
        original optimization process. Store the parameters before the
        `copy_to` method. After validation (or model saving), use this to
        restore the former parameters.
        Args:
          parameters: Iterable of `torch.nn.Parameter`; the parameters to be
            updated with the stored parameters.
        """
        for c_param, param in zip(self.collected_params, parameters):
            param.data.copy_(c_param.data)

    def state_dict(self):
        return dict(decay=self.decay, num_updates=self.num_updates,
                    shadow_params=self.shadow_params)

    def load_state_dict(self, state_dict, device):
        self.decay = state_dict['decay']
        self.num_updates = state_dict['num_updates']
        self.shadow_params = [tensor.to(device) for tensor in state_dict['shadow_params']]

def is_sugar_ring(ring, mol):
    if len(ring) not in [5, 6]:
        return False
    has_oxygen = any(mol.GetAtomWithIdx(idx).GetSymbol() == 'O' for idx in ring)
    if not has_oxygen:
        return False
    hydroxyl_count = 0
    for idx in ring:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetSymbol() == 'C':
            neighbors = atom.GetNeighbors()
            for neighbor in neighbors:
                if neighbor.GetSymbol() == 'O':
                    hydroxyl_count += 1
    return hydroxyl_count >= 1

def find_sugar_rings(mol):
    rings = mol.GetRingInfo().AtomRings()
    return [ring for ring in rings if is_sugar_ring(ring, mol)]

def analyze_sugar_rings(data_list, ligand_pred_list):
    pred_tensor = ligand_pred_list[-1]  # Shape: [2920, 1]
    return_list = defaultdict(list)
    for mol_idx, data in enumerate(data_list):
        mol = data['mol']
        if not isinstance(mol, Chem.Mol):
            print(f"Skipping molecule {mol_idx}: 'mol' is not an RDKit Mol object")
            continue

        # Verify heavy atom count
        mol_heavy_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1)

        # Extract predictions for this molecule
        start_idx = mol_idx * mol_heavy_atoms  # e.g., 0, 73, 146, ...
        end_idx = start_idx + mol_heavy_atoms  # e.g., 73, 146, 219, ...
        mol_pred_values = pred_tensor[start_idx:end_idx]  # Shape: [73, 1]

        # Find sugar rings
        sugar_rings = find_sugar_rings(mol)
        if not sugar_rings:
            print(f"Molecule {mol_idx}: No sugar rings found")
            continue

        for ring_idx, ring in enumerate(sugar_rings):
            # Get atom symbols for the ring
            ring_atoms = [mol.GetAtomWithIdx(idx).GetSymbol() for idx in ring]
            ring_composition = ''.join(ring_atoms)

            # Get predictions for atoms in this ring (adjust indices relative to molecule)
            ring_indices = torch.tensor(ring)
            # print(ring_indices)
            ring_pred_values = torch.sigmoid(mol_pred_values[ring_indices]).flatten().tolist()
            # print(ring_pred_values)
            avg_pred = sum(ring_pred_values) / len(ring_pred_values)

            # Print details
            return_list[mol_idx].append(f"Sugar Ring {ring_idx + 1}:")
            return_list[mol_idx].append(f"    Composition: {ring_composition}")
            return_list[mol_idx].append(f"    Atom Indices: {ring}")
            return_list[mol_idx].append(f"    Individual Predictions: {ring_pred_values}")
            return_list[mol_idx].append(f"    Average Prediction: {avg_pred:.4f}\n")
    return return_list

def read_pdb_heavy_atoms(pdb_path):
    """
    Read heavy atoms from a PDB file and group them by residue.
    Returns a list of (residue, atom_list) tuples.
    """
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_path)
    heavy_atoms_by_residue = []

    atom_idx = 0  # Global index for heavy atoms
    for model in structure:
        for chain in model:
            for residue in chain:
                atom_list = []
                for atom in residue:
                    if atom.element != 'H':  # Skip hydrogen
                        atom_list.append((atom_idx, atom.name, atom.element))
                        atom_idx += 1
                if atom_list:
                    heavy_atoms_by_residue.append((residue, atom_list))
    return heavy_atoms_by_residue

def analyze_pocket_residues(data_list, residue_pred_list):
    return_list = defaultdict(list)
    pred_tensor = residue_pred_list[-1]  # Shape: [3160, 1]
    
    for mol_idx, data in enumerate(data_list):
        # Get PDB path for this molecule
        pdb_path = data.get('pocket_path')
        if not pdb_path:
            print(f"Skipping molecule {mol_idx}: No 'pocket_path' found")
            continue

        # Read heavy atoms from PDB
        heavy_atoms_by_residue = read_pdb_heavy_atoms(pdb_path)
        num_residues = len(heavy_atoms_by_residue)
        
        if num_residues == 0:
            print(f"Skipping molecule {mol_idx}: No residues found in PDB")
            continue

        # Extract predictions for this molecule
        start_idx = sum(len(read_pdb_heavy_atoms(data_list[i]['pocket_path'])) 
                            for i in range(mol_idx)) if mol_idx > 0 else 0
        end_idx = start_idx + num_residues
        mol_pred_values = pred_tensor[start_idx:end_idx]

        for i, (residue, atom_list) in enumerate(heavy_atoms_by_residue):
            res_name = residue.resname
            res_id = residue.id[1]
            chain = residue.parent.id
            pred_value = torch.sigmoid(mol_pred_values[i]).item()

            return_list[mol_idx].append(f"  Residue {res_name} {res_id} (Chain {chain}):")
            return_list[mol_idx].append(f"    Prediction: {pred_value:.4f}\n")
    
    return return_list