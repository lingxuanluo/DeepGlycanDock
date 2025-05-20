import os
from argparse import ArgumentParser
from Bio.PDB import PDBParser
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from tqdm import tqdm
import pandas as pd

# Amino acid three-to-one letter code mapping
three_to_one = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'MSE': 'M', 'PHE': 'F',
    'PRO': 'P', 'PYL': 'O', 'SER': 'S', 'SEC': 'U', 'THR': 'T',
    'TRP': 'W', 'TYR': 'Y', 'VAL': 'V', 'ASX': 'B', 'GLX': 'Z',
    'XAA': 'X', 'XLE': 'J'
}

# Parse command-line arguments
parser = ArgumentParser(description="Convert PDB files to FASTA format")
parser.add_argument('--out_file', type=str, default="output.fasta", help="Output FASTA file path")
parser.add_argument('--protein_path', type=str, default=None, help="Path to a single PDB file")
parser.add_argument('--protein_ligand_csv', type=str, default=None, help="CSV file with protein paths")
args = parser.parse_args()

# Initialize PDB parser
biopython_parser = PDBParser(QUIET=True)

# Determine input PDB files
if args.protein_path:
    file_paths = [args.protein_path]
elif args.protein_ligand_csv:
    df = pd.read_csv(args.protein_ligand_csv)
    file_paths = list(set(df['protein_path'].tolist()))
else:
    raise ValueError("Must provide either --protein_path or --protein_ligand_csv")

# Process PDB files
sequences = []
ids = []
for file_path in tqdm(file_paths, desc="Processing PDB files"):
    try:
        structure = biopython_parser.get_structure('structure', file_path)
        structure = structure[0]  # Get first model
        for i, chain in enumerate(structure):
            seq = ''
            for residue in chain:
                if residue.get_resname() == 'HOH':  # Skip water molecules
                    continue
                # Check if residue is a valid amino acid by confirming presence of CA, N, and C atoms
                c_alpha = any(atom.name == 'CA' for atom in residue)
                n = any(atom.name == 'N' for atom in residue)
                c = any(atom.name == 'C' for atom in residue)
                if c_alpha and n and c:
                    res_name = residue.get_resname()
                    seq += three_to_one.get(res_name, '-')  # Use '-' for unknown residues
            if seq:  # Only add non-empty sequences
                sequences.append(seq)
                ids.append(f'{os.path.basename(file_path)}_chain_{i}')
    except Exception as e:
        print(f"Error processing {file_path}: {e}")

# Create SeqRecord objects
records = [
    SeqRecord(Seq(seq), id=str(index), description='')
    for index, seq in zip(ids, sequences) if seq
]

# Write to FASTA file
os.makedirs(os.path.dirname(args.out_file) or '.', exist_ok=True)
SeqIO.write(records, args.out_file, "fasta")
print(f"FASTA file written to {args.out_file}")