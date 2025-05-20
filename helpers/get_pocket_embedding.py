# by caoduanhua email: caodh@zju.edu.cn
# cycle chain and residues for a single PDB file
import pickle
import os
from Bio.PDB import *
import warnings
warnings.filterwarnings('ignore')
import tqdm
from Bio.PDB import PDBParser
import torch
import pandas as pd
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('--pdb_file', type=str, required=True, help='Path to the input PDB file')
parser.add_argument('--pocket_file', type=str, required=True, help='Path to the pocket PDB file')
parser.add_argument('--embeddings_dir', type=str, default='~/esm_embedding/esm_embedding_output', help='Full protein embedding dir location')
parser.add_argument('--pocket_emb_save_dir', type=str, default='~/esm_embedding/esm_embedding_output_pocket_new', help='Pocket embedding save directory')
args = parser.parse_args()

# Define three-to-one amino acid code mapping
three_to_one = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G',
    'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'MSE': 'M', 'PHE': 'F', 'PRO': 'P',
    'PYL': 'O', 'SER': 'S', 'SEC': 'U', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V', 'ASX': 'B',
    'GLX': 'Z', 'XAA': 'X', 'XLE': 'J'
}

# Initialize
biopython_parser = PDBParser()
os.makedirs(args.pocket_emb_save_dir, exist_ok=True)
Assertion_list = []

# Extract protein name from PDB file
protein_name = os.path.splitext(os.path.basename(args.pdb_file))[0]

# Skip if output already exists
if os.path.exists(os.path.join(args.pocket_emb_save_dir, f'{protein_name}.pt')):
    print(f'{protein_name}.pt already exists, skipping!')
else:
    try:
        # Load structures
        pocket_structure = biopython_parser.get_structure(protein_name, args.pocket_file)[0]
        full_structure = biopython_parser.get_structure(protein_name, args.pdb_file)[0]
        pocket_embeddings = []
        pocket_infos_all = []

        # Iterate over chains in full protein
        for i, chain in enumerate(full_structure.get_chains()):
            chain_id = chain.get_id()
            try:
                pocket_chain = pocket_structure[chain_id]
            except KeyError:
                print(f'Chain {chain_id} not in {protein_name} pocket, skipping this chain')
                continue

            try:
                # Load embeddings for the chain
                embeddings_path_chain = os.path.join(args.embeddings_dir, f'{protein_name}.pdb_chain_{i}.pt')
                embeddings = torch.load(embeddings_path_chain)['representations'][33]
                assert len(list(chain.get_residues())) == len(embeddings), 'Embedding length must match residue count!'

            except AssertionError:
                # Clean up residues (remove HOH and incomplete residues)
                residue_list = list(chain.get_residues())
                for residue in residue_list:
                    if residue.get_resname() == 'HOH':
                        chain.detach_child(residue.get_id())
                        continue
                    c_alpha, n, c = None, None, None
                    for atom in residue:
                        if atom.name == 'CA':
                            c_alpha = list(atom.get_vector())
                        if atom.name == 'N':
                            n = list(atom.get_vector())
                        if atom.name == 'C':
                            c = list(atom.get_vector())
                    if c_alpha is None or n is None or c is None:
                        chain.detach_child(residue.get_id())

                assert len(list(chain.get_residues())) == len(embeddings), \
                    f'Embedding must equal residue count! {len(list(chain.get_residues()))},{len(embeddings)}'

            # Process pocket residues
            pocket_infos = []
            pocket_residue_list = list(pocket_chain.get_residues())
            for residue in pocket_residue_list:
                if residue.get_resname() == 'HOH':
                    continue
                c_alpha, n, c = None, None, None
                for atom in residue:
                    if atom.name == 'CA':
                        c_alpha = list(atom.get_vector())
                    if atom.name == 'N':
                        n = list(atom.get_vector())
                    if atom.name == 'C':
                        c = list(atom.get_vector())
                if c_alpha is not None and n is not None and c is not None:
                    pocket_infos.append(residue.get_id())
                else:
                    print(f'Skipping residue {residue.get_resname()} due to missing atoms')

            pocket_infos_all.extend(pocket_infos)
            # Identify pocket residues in full chain
            pocket_idx_list = [
                res_idx for res_idx, res in enumerate(chain.get_residues())
                if res.get_id() in pocket_infos
            ]
            pocket_embeddings.append(embeddings[pocket_idx_list])

        # Concatenate embeddings and save
        pocket_embeddings = torch.cat(pocket_embeddings, dim=0)
        assert len(pocket_embeddings) == len(pocket_infos_all), \
            f'Pocket embedding must match residue count! {len(pocket_embeddings)},{len(pocket_infos_all)}'
        torch.save(pocket_embeddings, os.path.join(args.pocket_emb_save_dir, f'{protein_name}.pt'))
        print(f'Successfully processed {protein_name}')

    except AssertionError as e:
        print(f'AssertionError: {e} for {protein_name}')
        Assertion_list.append(protein_name)
    except FileNotFoundError as e:
        print(f'FileNotFoundError: {e} for {protein_name}')
        Assertion_list.append(protein_name)
    except Exception as e:
        print(f'Error: {e} for {protein_name}')
        Assertion_list.append(protein_name)

print('Assertion_list:', Assertion_list)