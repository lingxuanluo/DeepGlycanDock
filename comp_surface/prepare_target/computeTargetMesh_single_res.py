import os
import sys
import numpy as np
import shutil
import glob
import pymesh
import Bio.PDB
from Bio.PDB import * 
from rdkit import Chem
import warnings
warnings.filterwarnings("ignore")
from IPython.utils import io
from sklearn.neighbors import KDTree
from scipy.spatial import distance

from default_config.masif_opts import masif_opts
from compute_normal import compute_normal
from computeAPBS import computeAPBS
from computeCharges import computeCharges, assignChargesToNewMesh
from computeHydrophobicity import computeHydrophobicity
from computeMSMS import computeMSMS
from fixmesh import fix_mesh
from save_ply import save_ply
from mol2graph import *

def parse_residue_list(residue_input):
    """Parse a residue input string like '1-50,60,62,63' into a list of residue numbers."""
    residue_numbers = set()
    parts = residue_input.split(',')
    for part in parts:
        if '-' in part:
            start, end = map(int, part.split('-'))
            residue_numbers.update(range(start, end + 1))
        else:
            residue_numbers.add(int(part))
    return residue_numbers

def compute_inp_surface(target_filename, residue_input, out_dir=None):
    sufix = f'_{residue_input.replace(",", "_").replace("-", "to")}.pdb'
    if out_dir is not None: 
        out_filename = os.path.join(out_dir, os.path.splitext(target_filename)[0].split('/')[-1])
        os.makedirs(out_filename, exist_ok=True)
        sufix = f'/{os.path.splitext(target_filename)[0].split("/")[-1]}_{residue_input.replace(",", "_").replace("-", "to")}.pdb'
    else:
        out_filename = os.path.splitext(target_filename)[0]
    
    if os.path.exists(out_filename + f"/{sufix.split('.')[0]}.ply"):
        print('Already done, skipping!')
        return 0
    
    input_filename = os.path.splitext(target_filename)[0]
    
    # Parse residue input
    selected_residues = parse_residue_list(residue_input)
    
    # Read protein and select specified residues
    parser = Bio.PDB.PDBParser(QUIET=True)
    structures = parser.get_structure('target', input_filename + '.pdb')
    structure = structures[0]

    class SelectNeighbors(Select):
        def accept_residue(self, residue):
            residue_id = residue.get_id()[1]  # Get residue number
            if residue_id in selected_residues:
                if all(a in [i.get_name() for i in residue.get_unpacked_list()] for a in ['N', 'CA', 'C', 'O']) or residue.resname == 'HOH':
                    return True
                else:
                    return False
            return False

    pdbio = PDBIO()
    pdbio.set_structure(structure)
    pdbio.save(out_filename + sufix, SelectNeighbors())
    
    # Compute MSMS of surface with hydrogens
    try:
        vertices1, faces1, normals1, names1, areas1 = computeMSMS(out_filename + sufix, protonate=True, one_cavity=None)
        
        # Since no ligand is used, select all vertices for the mesh
        faces_to_keep = list(range(len(faces1)))  # Keep all faces
        
        # Compute "charged" vertices
        if masif_opts['use_hbond']:
            vertex_hbond = computeCharges(input_filename, vertices1, names1)    
    
        # Assign hydrophobicity
        if masif_opts['use_hphob']:
            vertex_hphobicity = computeHydrophobicity(names1)    
        
        # If protonate = false, recompute MSMS without hydrogens
        vertices2 = vertices1
        faces2 = faces1
    
        # Fix the mesh
        mesh = pymesh.form_mesh(vertices2, faces2)
        mesh = pymesh.submesh(mesh, faces_to_keep, 0)
        with io.capture_output() as captured:
            regular_mesh = fix_mesh(mesh, masif_opts['mesh_res'])
        
    except Exception as e:
        print(f"Error in MSMS computation: {str(e)}")
        vertices1, faces1, normals1, names1, areas1 = computeMSMS(out_filename + sufix, protonate=True, one_cavity=None)
        
        faces_to_keep = list(range(len(faces1)))
        
        if masif_opts['use_hbond']:
            vertex_hbond = computeCharges(input_filename, vertices1, names1)    
        
        if masif_opts['use_hphob']:
            vertex_hphobicity = computeHydrophobicity(names1)    
        
        vertices2 = vertices1
        faces2 = faces1
    
        mesh = pymesh.form_mesh(vertices2, faces2)
        mesh = pymesh.submesh(mesh, faces_to_keep, 0)
        with io.capture_output() as captured:
            regular_mesh = fix_mesh(mesh, masif_opts['mesh_res'])
    
    # Compute the normals
    vertex_normal = compute_normal(regular_mesh.vertices, regular_mesh.faces)
    
    # Assign charges to new vertices
    if masif_opts['use_hbond']:
        vertex_hbond = assignChargesToNewMesh(regular_mesh.vertices, vertices1, vertex_hbond, masif_opts)
    
    if masif_opts['use_hphob']:
        vertex_hphobicity = assignChargesToNewMesh(regular_mesh.vertices, vertices1, vertex_hphobicity, masif_opts)
    
    if masif_opts['use_apbs']:
        vertex_charges = computeAPBS(regular_mesh.vertices, out_filename + sufix, out_filename + "_temp")
        
    # Compute the principal curvature components for the shape index
    regular_mesh.add_attribute("vertex_mean_curvature")
    H = regular_mesh.get_attribute("vertex_mean_curvature")
    regular_mesh.add_attribute("vertex_gaussian_curvature")
    K = regular_mesh.get_attribute("vertex_gaussian_curvature")
    elem = np.square(H) - K
    elem[elem < 0] = 1e-8
    k1 = H + np.sqrt(elem)
    k2 = H - np.sqrt(elem)
    si = (k1 + k2) / (k1 - k2)
    si = np.arctan(si) * (2 / np.pi)
    
    save_ply(out_filename + f"/{sufix.split('.')[0]}.ply", regular_mesh.vertices,
             regular_mesh.faces, normals=vertex_normal, charges=vertex_charges,
             normalize_charges=True, hbond=vertex_hbond, hphob=vertex_hphobicity,
             si=si)
    
    return 0

if __name__ == "__main__":
    from joblib import delayed, Parallel
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('--pdb_file', type=str, required=True, help='Path to the input PDB file')
    parser.add_argument('--residue_list', type=str, required=True, help='Comma-separated list of residue numbers or ranges (e.g., "1-50,60,62,63")')
    parser.add_argument('--out_dir', type=str, default='./surface_output', help='Output directory for surface files')
    args = parser.parse_args()
    
    # Convert to absolute paths
    args.pdb_file = os.path.abspath(args.pdb_file)
    args.out_dir = os.path.abspath(args.out_dir)
    
    # Ensure output directory exists
    os.makedirs(args.out_dir, exist_ok=True)
    os.chdir(args.out_dir)

    # Validate input file
    if not os.path.exists(args.pdb_file):
        raise FileNotFoundError(f"PDB file {args.pdb_file} does not exist")

    # Process the input
    try:
        res = compute_inp_surface(args.pdb_file, args.residue_list, args.out_dir)
        print(f"Processed: {args.pdb_file} with residues {args.residue_list}")
    except Exception as e:
        print(f"Error processing {args.pdb_file} with residues {args.residue_list}: {str(e)}")

    # Clean up temporary files
    files = glob.glob(os.path.join(args.out_dir, '*_temp*')) + glob.glob(os.path.join(args.out_dir, '*msms*'))
    for f in files:
        os.remove(f)
    print(f"Cleaned up {len(files)} temporary files")