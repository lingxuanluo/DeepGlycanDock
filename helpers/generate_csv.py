import pandas as pd
from pathlib import Path
import os
def generate_single_case_csv(protein_path, ligand_path, ref_ligand_path, output_csv):
    """
    Generate a CSV file with paths for a single protein, ligand, ref ligand, pocket, and surface.
    
    Args:
        protein_path (str): Path to the protein PDB file (e.g., single_case/8ufh_y75_prepared.pdb)
        ligand_path (str): Path to the ligand SDF file (e.g., single_case/ss-84.sdf)
        ref_ligand_path (str): Path to the reference ligand PDB file (e.g., single_case/ref_lig.pdb)
        output_csv (str): Path to save the output CSV file
    """
    # Convert input paths to Path objects
    protein_path = Path(protein_path)
    ligand_path = Path(ligand_path)
    ref_ligand_path = Path(ref_ligand_path)
    
    # Extract directory and protein name
    case_dir = protein_path.parent  # e.g., single_case
    protein_name = protein_path.stem  # e.g., 8ufh_y75_prepared
    
    # Construct pocket and surface paths based on the example
    pocket_path = case_dir / "surface" / case_dir.name / f"{protein_name}_ref_lig_8A.pdb"
    surface_path = case_dir / "surface" / case_dir.name / f"{protein_name}_ref_lig_8A.ply"
    
    # Check if all files exist
    if not (protein_path.exists() and 
            ligand_path.exists() and 
            ref_ligand_path.exists() and 
            pocket_path.exists() and 
            surface_path.exists()):
        raise FileNotFoundError("One or more required files are missing.")
    
    # Create DataFrame
    df = pd.DataFrame({
        'protein_path': [os.path.abspath(protein_path)],
        'ligand_path': [os.path.abspath(ligand_path)],
        'ref_ligand': [os.path.abspath(ref_ligand_path)],
        'pocket_path': [os.path.abspath(pocket_path)],
        'protein_surface': [os.path.abspath(surface_path)]
    })
    
    # Save to CSV
    df.to_csv(output_csv, index=False)
    print(f"CSV file saved to {output_csv}")

if __name__ == "__main__":
    # Example usage with provided paths
    parser = argparse.ArgumentParser(description="Generate CSV for protein-ligand data")
    parser.add_argument('--protein', required=True, help='Path to protein PDB file')
    parser.add_argument('--ligand', required=True, help='Path to ligand SDF file')
    parser.add_argument('--ref_ligand', required=True, help='Path to reference ligand PDB file')
    parser.add_argument('--output', required=True, help='Path to output CSV file')
    
    args = parser.parse_args()
    
    generate_single_case_csv(args.protein, args.ligand, args.ref_ligand, args.output)
