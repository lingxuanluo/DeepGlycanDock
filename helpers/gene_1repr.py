from unimol_tools import UniMolRepr
import numpy as np
import os, pickle
from rdkit import Chem
import sys

# Set up parameters
params = {
    'data_type': 'molecule',
    'remove_hs': False,
    'model_name': 'unimolv2',
    'model_size': '310m-40',
}

# Initialize UniMolRepr
clf = UniMolRepr(use_gpu='True', **params)

# Check if SDF file path is provided as command-line argument
if len(sys.argv) != 2:
    print("Usage: python script.py <sdf_file_path>")
    sys.exit(1)

sdf_file = sys.argv[1]

# Verify file exists and is an SDF file
if not os.path.exists(sdf_file) or not sdf_file.endswith('.sdf'):
    print("Error: Please provide a valid SDF file path")
    sys.exit(1)

# Generate output pickle filename
filename = os.path.basename(sdf_file)
file_key = os.path.splitext(filename)[0]
output_dir = os.path.dirname(sdf_file)
output_file = os.path.join(output_dir, f"{file_key}_unimol_repr.pkl")

# Check if output file already exists
if os.path.exists(output_file):
    print(f"Output file {output_file} already exists, skipping generation")
    sys.exit(0)

# Dictionary to store results
results_dict = {}

# Process the SDF file
# Read SDF file
supplier = Chem.SDMolSupplier(sdf_file)
mol = next(supplier)  # Get the first molecule from the SDF file

if mol is not None:
    # Convert to SMILES
    smiles = Chem.MolToSmiles(mol)
    
    # Get UniMol representation
    unimol_repr = clf.get_repr(data=[smiles], return_atomic_reprs=False)
    cls_repr = np.array(unimol_repr["cls_repr"])
    
    # Store in dictionary
    results_dict[file_key] = cls_repr
    
    # Save results to a pickle file
    with open(output_file, 'wb') as f:
        pickle.dump(results_dict, f)
    print(f"Saved representations to {output_file}")
else:
    print(f"Error: Could not process molecule from {sdf_file}")