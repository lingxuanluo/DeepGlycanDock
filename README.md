# DeepGlycanDock: A Surface-Based Diffusion Model for Glycan Binding Pose Prediction
DeepGlycanDock is a cutting-edge tool for accurate and reliable protein-glycan docking, developed by Xinheng He. Leveraging glycan pretraining, ESM3 protein embeddings, and glycan-ring-specific information, it delivers superior performance in predicting glycan binding poses.
This repository provides the code, instructions, and model weights to:

Generate accurate protein-ligand complexes, screen compounds using DeepGlycanDock, and evaluate the tool's performance.

For questions or support, please open an issue or contact us at he-xinheng@foxmail.com.
## Table of Contents
1. [Environment Setup](#Environment Setup)
2. [Running DeepGlycanDock on Your Complexes](#Running DeepGlycanDock on Your Complexes)
3. [License](#license)

## Environment Setup
Follow these steps to configure the DeepGlycanDock environment:

Create a Conda Environment
conda create -y -n DeepGlycanDock python=3.10
conda activate DeepGlycanDock


Install Mamba and Clean Cache
conda install -y -c conda-forge -c pytorch -c pyg mamba
conda clean -ya


Install PyTorch with CUDA Support
mamba install -y -c pytorch -c nvidia pytorch pytorch-cuda=12.1


Install Core Dependencies
mamba install -y -c conda-forge -c pytorch -c pyg numpy==1.20 scipy==1.8.1 pandas==2.1.2
conda clean -ya


Install Molecular Modeling Dependencies
mamba install -y -c conda-forge -c pytorch -c pyg openff-toolkit==0.15.2 openmm==8.1.1 openmmforcefields==0.12.0 pdbfixer==1.9
conda clean -ya


Install Additional Dependencies
mamba install -y -c conda-forge -c pytorch -c pyg babel==2.13.1 biopandas==0.4.1 openbabel==3.1.1 plyfile==1.0.1 prody==2.4.0 torch-ema==0.3 torchmetrics==1.2.1
conda clean -ya


Install PyTorch Geometric (PyG)
mamba install -y -c pyg pyg
pip install pyg_lib torch_scatter torch_sparse torch_cluster torch_spline_conv -f https://data.pyg.org/whl/torch-2.2.0+cu121.html


Install Remaining Dependencies
pip install -U --no-cache-dir spyrmsd scikit-learn==1.3.2 accelerate==0.15.0 biopython==1.79 e3nn==0.5.1 huggingface-hub==0.17.3 mdanalysis==2.4.0 posebusters==0.2.7 rdkit==2023.3.1 tokenizers==0.13.3 transformers==4.29.2 wandb==0.16.1


Install PyMesh
pip install pymesh
pip install https://github.com/nuvolos-cloud/PyMesh/releases/download/v0.3.1/pymesh2-0.3.1-cp310-cp310-linux_x86_64.whl


Install Additional Utilities
mamba install -y loguru
pip install dimorphite_dl prefetch_generator esm==3.1.1 --no-deps
pip install zstd==1.5.6.6 biotite==0.41.2 cloudpathlib==0.21.0 msgpack-numpy==0.4.8 tenacity==9.0.0 torchtext==0.17.2 --no-deps


Install MaSIF and Data Processing Dependencies
mamba install -y -c mx reduce
mamba install -y -c conda-forge openbabel


These scripts demonstrate how to use DeepGlycanDock for evaluation and compound screening.
## Running DeepGlycanDock on Your Complexes
To perform docking on your protein-glycan complexes, use the following commands:
Basic Docking
python test_my_glycan.py \
  --protein <path_to_protein> \
  --ligand <path_to_sdf_ligand> \
  --ref_ligand <path_to_pdb_ligand_define_site> \
  --output_dir ./output

Docking with Specified Pocket Residues
To specify pocket residues for docking:
python test_my_glycan_res.py \
  --protein single_case/8ufh_y75_prepared.pdb \
  --ligand single_case/ss-84.sdf \
  --pocket_residue "1-50,60,61,62,63" \
  --output_dir ./output

Replace the file paths and residue numbers with your own data as needed.
## License
This project is licensed under the MIT License.
