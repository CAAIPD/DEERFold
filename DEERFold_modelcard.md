# Model Card for DEERFold

DEERFold is a modified AlphaFold2 framework that integrates distance distributions from Double Electron-Electron Resonance (DEER) spectroscopy into protein structure prediction. It generates accurate and diverse protein conformations guided by experimental and simulated DEER data, outperforming base AlphaFold2 for proteins with limited homologs by producing multiple conformational states with fewer experimental constraints.

# Table of Contents

1. Model Details
2. Uses
3. Bias, Risks, and Limitations
4. How To Get Started With the Model
5. Training Details
6. Evaluation
7. Environmental Impact
8. Citation

## Model Details

### Model Description

DEERFold enhances protein structure prediction by incorporating DEER spectroscopy data into the AlphaFold2 framework. It leverages distance distributions to guide the generation of diverse and accurate protein conformations, particularly excelling for proteins with limited evolutionary homologs. By reducing the number of required experimental constraints, DEERFold increases throughput and extends its applicability to other probe-based methods, bridging experimental data and deep learning for efficient conformational modeling.

- **Developed by:** Tianqi Wu, Richard A. Stein, Te-Yu Kao, Benjamin Brown, Hassane S. Mchaourab
- **Shared by:** Center for Applied Artificial Intelligence in Protein Dynamics (CAAIPD)
- **Model type:** Modified AlphaFold2 for protein structure prediction with DEER data
- **License:** MIT License

### Model Sources

- **Repository:** https://github.com/CAAIPD/DEERFold
- **Model_weights & test data:** https://zenodo.org/records/15147304
- **Preprint:** https://doi.org/10.1101/2024.10.30.621127

## Uses

### Direct Use

DEERFold is designed for research purposes, enabling the prediction of protein structures using DEER spectroscopy data. It supports both constrained predictions (guided by DEER distance constraints) and unconstrained predictions (without DEER data), making it versatile for studying protein conformations.

### Out-of-Scope Use

DEERFold is tailored for protein structure prediction with DEER data and is not suitable for other biological sequences (e.g., DNA) or non-protein-related tasks such as natural language processing.

## Bias, Risks, and Limitations

DEERFold’s performance depends heavily on the quality of input DEER data. Compromised or low-quality DEER data may lead to erroneous models. The model is optimized for protein structures and may not generalize well to other data types or applications outside its intended scope.

## How to Get Started with the Model

To use DEERFold, follow these steps:

1. **Clone the repository:**

   ```
   git clone https://github.com/CAAIPD/DEERFold.git
   cd DEERFold
   ```

2. **Install DEERFold (takes \~20 minutes):**

   ```
   bash install_deerfold.sh
   ```

3. **Activate the environment:**

   ```
   export PATH="$HOME/miniforge3/bin:$PATH"
   conda activate ./env/deerfold_venv
   ```

4. **Download the model weights:**

   ```
   wget -O DEERFold.tar.gz "https://zenodo.org/records/15147304/files/DEERFold.tar.gz?download=1"
   tar xvzf DEERFold.tar.gz
   ```

5. **Run inference:**

   - **Unconstrained prediction:**

     ```
     python deerfold_inference.py <fasta_file> <out_dir> --model <model_weights_dir> --msa_dir <msa_dir> --neff <neff> --num <num>
     ```

     Example:

     ```
     python deerfold_inference.py examples/PfMATE/PfMATE.fasta out/PfMATE_unconstrained --model model/DEERFold.pt --msa_dir examples/msa --neff 5 --num 15
     ```
   - **Constrained prediction:**

     ```
     python deerfold_inference.py <fasta_file> <out_dir> --model <model_weights_dir> --splabel <csv_file> --msa_dir <msa_dir> --neff <neff> --num <num>
     ```

     Example:

     ```
     python deerfold_inference.py examples/PfMATE/PfMATE.fasta out/PfMATE_constrained --model model/DEERFold.pt --splabel examples/PfMATE/experiment.csv --msa_dir examples/msa --neff 5 --num 15
     ```

For detailed examples, refer to the examples/deerfold.ipynb notebook.

### Options

- `fasta_file`: Input sequence in FASTA format
- `out_dir`: Output directory for generated models
- `model`: Path to model weights (e.g., `model/DEERFold.pt`)
- `msa_dir`: Directory with multiple sequence alignments
- `neff`: MSA Neff value
- `num`: Number of models to generate
- `splabel`: CSV file with DEER constraints (for constrained prediction)
- `ref_pdbs`: Optional reference PDB files for RMSD/TM-score analysis

## Training Details

### Training Data

DEERFold was trained using:

- **Multiple Sequence Alignments (MSAs):** Sourced from the OpenFold dataset, containing 401,381 MSAs for 140,000 PDB chains and 16,000,000 UniClust30 clusters.
- **Protein Structures:** Extracted from the OpenFold dataset as mmCIF files.
- **DEER Data:** Available for download from [Zenodo](https://zenodo.org/records/1514730) or generated using provided scripts.

### Training Procedure

DEERFold was fine-tuned from AlphaFold2’s pretrained weights (`model_5_ptm`) without templates, following the OpenFold training protocol. For specifics, see the [OpenFold documentation](https://openfold.readthedocs.io/en/latest/).

## Evaluation

DEERFold’s performance is assessed by its ability to generate protein conformations consistent with input DEER constraints. Outputs are ranked by Earth Mover’s Distance (EMD) to the constraints, with top-ranking models best fitting the data. If reference PDB files are provided, RMSD and TM-score metrics can be computed. Detailed evaluation results are available in the preprint.

## Environmental Impact

- **Hardware Type:** 80GB NVIDIA A100 GPUs
- **Hours used:** 336 (14 days)
- **Platform:** Vanderbilt Data Science Institute - DGX A100

## Citation

If you use DEERFold or its results, please cite the preprint:

```
@article {Wu2024.10.30.621127,
	author = {Wu, Tianqi and Stein, Richard A. and Kao, Te-Yu and Brown, Benjamin and Mchaourab, Hassane S.},
	title = {Modeling Protein Conformations by Guiding AlphaFold2 with Distance Distributions. Application to Double Electron Electron Resonance (DEER) Spectroscopy},
	elocation-id = {2024.10.30.621127},
	year = {2024},
	doi = {10.1101/2024.10.30.621127},
	publisher = {Cold Spring Harbor Laboratory},
	eprint = {https://www.biorxiv.org/content/early/2024/11/01/2024.10.30.621127.full.pdf},
	journal = {bioRxiv}
}
```
