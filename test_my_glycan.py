import argparse
import os
import subprocess
from pathlib import Path
import pandas as pd
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def run_command(command, step_name):
    """Run a shell command and handle errors."""
    try:
        result = subprocess.run(command, shell=True, check=True, capture_output=True, text=True)
        logger.info(f"{step_name} completed successfully.")
        logger.debug(f"{step_name} output: {result.stdout}")
    except subprocess.CalledProcessError as e:
        logger.error(f"Error in {step_name}: {e.stderr}")
        raise

def generate_single_case_csv(protein_path, ligand_path, ref_ligand_path, pocket_residue, output_csv):
    """Generate a CSV file with paths for a single protein, ligand, ref ligand or pocket, and surface."""
    protein_path = Path(protein_path)
    ligand_path = Path(ligand_path)
    
    case_dir = protein_path.parent
    protein_name = protein_path.stem
    
    if pocket_residue:
        residue_str = pocket_residue.replace(",", "_").replace("-", "to")
        pocket_path = case_dir / "surface" / protein_name / f"{protein_name}_{residue_str}.pdb"
        surface_path = case_dir / "surface" / protein_name / f"{protein_name}_{residue_str}.ply"
    else:
        ref_ligand_path = Path(ref_ligand_path)
        ref_ligand_name = ref_ligand_path.stem
        pocket_path = case_dir / "surface" / case_dir.name / f"{protein_name}_{ref_ligand_name}_8A.pdb"
        surface_path = case_dir / "surface" / case_dir.name / f"{protein_name}_{ref_ligand_name}_8A.ply"
    
    if not (protein_path.exists() and ligand_path.exists()):
        raise FileNotFoundError("One or more input files are missing.")
    if ref_ligand_path and not Path(ref_ligand_path).exists():
        raise FileNotFoundError("Reference ligand file is missing.")
    
    df = pd.DataFrame({
        'protein_path': [str(protein_path)],
        'ligand_path': [str(ligand_path)],
        'ref_ligand': [str(ref_ligand_path)] if ref_ligand_path else None,
        'pocket_path': [str(pocket_path)],
        'protein_surface': [str(surface_path)]
    })
    
    output_csv = Path(output_csv)
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_csv, index=False)
    logger.info(f"CSV file saved to {output_csv}")
    return str(output_csv)

def run_docking_pipeline(protein_path, ligand_path, ref_ligand_path, pocket_residue, output_dir, large_ligand=False):
    """Run the full docking pipeline with the specified inputs."""
    protein_path = Path(protein_path)
    ligand_path = Path(ligand_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    case_dir = protein_path.parent
    protein_name = protein_path.stem
    ligand_name = ligand_path.stem
    
    if pocket_residue:
        residue_str = pocket_residue.replace(",", "_").replace("-", "to")
        pocket_embedding_path = case_dir / f"{protein_name}_{residue_str}_8A.pt"
    else:
        ref_ligand_path = Path(ref_ligand_path)
        ref_ligand_name = ref_ligand_path.stem
        pocket_embedding_path = case_dir / f"{protein_name}_{ref_ligand_name}_8A.pt"
    
    # Paths for intermediate and output files
    csv_path = case_dir / "output_paths.csv"
    fasta_path = case_dir / f"{protein_name}.fasta"
    esm_embedding_dir = case_dir / "esm_embedding"
    esm_embedding_save_dir = case_dir / "esm_embedding_save"
    ligand_embedding_path = case_dir / f"{ligand_name}_unimol_repr.pkl"
    surface_dir = case_dir / "surface"
    
    # Step 1: Generate ligand representation
    logger.info("Step 1: Generating ligand representation")
    cmd1 = f"python ./helpers/gene_1repr.py {ligand_path}"
    run_command(cmd1, "Ligand Representation")
    
    # Step 2: Generate protein surface
    logger.info("Step 2: Generating protein surface")
    surface_dir.mkdir(parents=True, exist_ok=True)
    if pocket_residue:
        cmd2 = f"python ./comp_surface/prepare_target/computeTargetMesh_single_res.py --pdb_file {protein_path} --residue_list \"{pocket_residue}\" --out_dir {surface_dir}"
    else:
        cmd2 = f"python ./comp_surface/prepare_target/computeTargetMesh_single.py --pdb_file {protein_path} --ligand_file {ref_ligand_path} --out_dir {surface_dir}"
    run_command(cmd2, "Protein Surface Generation")
    # print(cmd2)
    
    # Step 3: Generate ESM3 protein embeddings
    logger.info("Step 3: Generating ESM3 protein embeddings")
    esm_embedding_dir.mkdir(parents=True, exist_ok=True)
    cmd3a = f"python ./helpers/fasta_extract_pdb.py --protein_path {protein_path} --out_file {fasta_path}"
    cmd3b = f"python ./helpers/e3_test.py --fasta_path {fasta_path} --output_dir {esm_embedding_dir}"
    run_command(cmd3a, "FASTA Extraction")
    run_command(cmd3b, "ESM3 Embedding Generation")
    
    # Step 4: Generate CSV file
    logger.info("Step 4: Generating CSV file")
    generate_single_case_csv(protein_path, ligand_path, ref_ligand_path, pocket_residue, csv_path)
    
    # Step 5: Map pocket ESM embedding
    logger.info("Step 5: Mapping pocket ESM embedding")
    esm_embedding_save_dir.mkdir(parents=True, exist_ok=True)
    cmd5 = f"python ./datasets/get_pocket_embedding.py --protein_pocket_csv {csv_path} --embeddings_dir {esm_embedding_dir} --pocket_emb_save_dir {esm_embedding_save_dir}"
    run_command(cmd5, "Pocket ESM Embedding Mapping")
    
    # Step 6: Save pocket ESM embedding to single file
    logger.info("Step 6: Saving pocket ESM embedding to single file")
    cmd6 = f"python ./datasets/esm_pocket_embeddings_to_pt.py --esm_embeddings_path {esm_embedding_save_dir} --output_path {pocket_embedding_path}"
    run_command(cmd6, "Pocket ESM Embedding Saving")
    
    # Step 7: Final docking
    logger.info("Step 7: Running final docking")
    docking_out_dir = output_dir
    docking_out_dir.mkdir(parents=True, exist_ok=True)
    cmd7 = (
        f"python -m accelerate.commands.launch --main_process_port 29515 --num_processes 1 ./inference_accelerate_glycan{'_res' if pocket_residue else ''}.py "
        f"--data_csv {csv_path} "
        f"--model_dir ./model_weights/glycan_finetune "
        f"--ckpt epoch_best_model.pt "
        f"--confidence_model_dir ./model_weights/posepredict "
        f"--confidence_ckpt best_model.pt "
        f"--save_docking_result "
        f"--mdn_dist_threshold_test 3 "
        f"--esm_embeddings_path {pocket_embedding_path} "
        f"--run_name ./model_weights/posepredict_test_dist_3 "
        f"--project 'DeepGlycanDock'_eval_samples/carb_test "
        f"--out_dir {docking_out_dir} "
        f"--batch_size 40 "
        f"--batch_size_molecule 1 "
        f"--samples_per_complex 40 "
        f"--save_docking_result_number 40 "
        f"--head_index 0 "
        f"--tail_index 10000 "
        f"--inference_mode Screen "
        f"--ligand_embeddings_path {ligand_embedding_path}"
    )
    
    if large_ligand:
        cmd7 += " --keep_input_pose --force_optimize  "
    # print(cmd7)
    run_command(cmd7, "Final Docking")
    
    logger.info(f"Docking pipeline completed. Results saved in {docking_out_dir}")

def main():
    parser = argparse.ArgumentParser(description="Run docking pipeline with protein, ligand, and either reference ligand or pocket residue inputs.")
    parser.add_argument(
        "--protein",
        default="single_case/8ufh_y75_prepared.pdb",
        help="Path to protein PDB file (default: single_case/8ufh_y75_prepared.pdb)"
    )
    parser.add_argument(
        "--ligand",
        default="single_case/ss-84.sdf",
        help="Path to ligand SDF file (default: single_case/ss-84.sdf)"
    )

    parser.add_argument(
        "--ref_ligand",
        default=None,
        help="Path to reference ligand PDB file (default: None)"
    )
    parser.add_argument(
        "--pocket_residue",
        default=None,
        help="Comma-separated list of residue numbers or ranges for pocket (e.g., '1-50,60,61,62,63') (default: None)"
    )
    parser.add_argument(
        "--output_dir",
        default="./output2",
        help="Directory to save docking results (default: ./output2)"
    )
    parser.add_argument(
        "--large_ligand",
        action="store_true",
        help="If set, adds --keep_input_pose to the final docking command"
    )
    
    args = parser.parse_args()
    
    if args.pocket_residue and args.ref_ligand:
        raise ValueError("Cannot specify both --pocket_residue and --ref_ligand. Choose one.")
    if not args.pocket_residue and not args.ref_ligand:
        raise ValueError("Must specify either --pocket_residue or --ref_ligand.")
    
    try:
        run_docking_pipeline(args.protein, args.ligand, args.ref_ligand, args.pocket_residue, args.output_dir, args.large_ligand)
    except Exception as e:
        logger.error(f"Pipeline failed: {str(e)}")
        raise

if __name__ == "__main__":
    main()