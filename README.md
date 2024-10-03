# DEERFold

### Description
In this work, we introduce DEERFold, a modified AlphaFold2 framework that incorporates distance constraint distributions from Double Electron-Electron Resonance (DEER) spectroscopy into protein structure prediction. DEERFold generates accurate and diverse protein conformations guided by experimental and simulated distance distributions, outperforming base AlphaFold2 for proteins with limited homologs to generate multiple conformational states. Critically, DEERFold reduces the number of required experimental constraints, increasing throughput and expanding capabilities beyond DEER to other probe-based methods. We envision DEERFold will advance protein structure prediction by bridging the gap between experimental data and deep learning techniques, enabling more accurate and efficient conformational modeling.

DEERFold is described in this [preprint](); if you use the code from this repository or the results, please cite the preprint.

<p align="center">
<img src="img/Fig1.png" width="1200" alt="Workflow of DEERFold">
</p>


<div style="display: flex; justify-content: center; gap: 20px;">
    <div>
        <img src="img/LmrP.gif" width="350" height="300" alt="LmrP" style="display: block;">
    </div>
    <div>
        <img src="img/PfMATE.gif" width="350" height="300" alt="PfMATE" style="display: block;">
    </div>
</div>

----

## Table of contents

- [DEERFold](#deerfold)
- [Installation](#installation)
    - [Datasets](#datasets)
    - [Loading pretrained models](#loading-pretrained-models)
- [Available models](#available-models)
- [Unconstrained prediction](#unconstrained-prediction) 
- [Constrained prediction](#constrained-prediction)
- [Analysis](#analysis-of-generations)

----

## Installation
To download our code, we recommend creating the conda environment with deerfold_environment.yml.
```
conda env create -f deerfold_environment.yml
```
To run DEERFold, you need to install the following tools: 
* [PyTorch v1.12](https://pytorch.org/get-started/previous-versions/)
```
pip install torch==1.12.1+cu113 torchvision==0.13.1+cu113 torchaudio==0.12.1 --extra-index-url https://download.pytorch.org/whl/cu113
```
* [Openfold 1.0.1](https://github.com/aqlaboratory/openfold/releases) As DEERfold is built on OpenFold, please refer to the [OpenFold GitHub](https://github.com/aqlaboratory/openfold#installation-linux) for installation instructions.  
* [chiLife 1.0.0.dev3](https://github.com/stollLab/chiLife/)
* [TM-score](https://zhanggroup.org/TM-score/)

We provide a notebook with installation guidance that can be found in [examples/deerfold.ipynb](https://github.com/CAAIPD/DEERFold/tree/main/examples/deerfold.ipynb). It also includes examples on how to generate predictions using our models. We recommend following this notebook if you would like to use our models to generate proteins.

**Available DEERFold models are:**
* ``` DEERFold_helix.pt: DEERFold model trained on exposed sites in helical regions```
* ``` DEERFold_strand.pt: DEERFold model trained on exposed sites in beta strand regions ```
* ``` DEERFold_com.pt: DEERFold model trained on exposed sites in both helical and beta strand regions  ```

Model weights can be downloaded [here]().

## Inference
### Constrained prediction
DEERFold can accept input DEER distance constraints in csv format, the format is shown as below:
```
18,95,0,0,0,0,2.00E-05,7.00E-05,0.00027,0.0009,0.0027,0.00708,0.01622,0,0,0,0,...
18,215,0,0,0,0,0,0,8.00E-05,0.000940009,0.007260073,0.035440354,0.109791098,0.215792158,0.269092691,0.212892129,0.106851069,...
18,240,0,0,1.00E-05,4.00E-05,0.000169995,0.00067998,0.00233993,0.006879794,0.017259482,0.036938892,0.067427977,0.10499685,0.139475816,0.158035259,0.152745418,...
```
The input distance constraints are in the format of distograms (shape LxLx100, each bin, except for the first and last, spans 1 Å and is centered at integer values from 1 to 99 Å. The first bin includes all values up to and including 1.5 Å, while the last bin captures all values above 99.5 Å).

To run inference on a sequence with the given DEER constraints, run the following script:
```
python predict_with_sp_batch.py <fasta_file> <msa_dir> <out_dir> --model <model_weights_dir> --sp <csv_file> --neff neff --num num --refs <ref1.pdb,ref2.pdb>

Options are:
--model   The directory which stores all the model weights .
--sp      Input file in csv format.
--neff    MSA neff value.
--num     number of models generated
--refs    reference model used for TMscore or RMSD analysis
```
* E.g. Generate 100 constrained models from DEERFold_helix based on the input DEER constraints(**--sp** examples/PfMATE/PfMATE_low.csv), set MSA Neff as 5
```
python predict_with_sp_batch.py examples//PfMATE/PfMATE.fasta examples/alignments out/PfMATE --model model --sp examples/PfMATE/PfMATE_low.csv --neff 5 --num 100
```
* E.g. If you have the reference models available, DEERFold provides RMSD and TM-score analysis of the prediction results compared to the reference PDB files(**--ref**).
```
python predict_with_sp_batch.py examples//PfMATE/PfMATE.fasta examples/alignments out/PfMATE --model model --sp examples/PfMATE/PfMATE_low.csv --neff 5 --num 100 --refs examples/PfMATE/6gwh.pdb,examples/PfMATE/6fhz.pdb
```

Outputs should be DEERFold predicted models in the format of PDB files, which are ranked by EMD distance between the prediction and the input distance constraints. **The top ranking models should be those most closely fitting the input distance constraints.**

An example of constrained model generation from DEERFold can be found in
[this notebook](https://github.com/CAAIPD/DEERFold/tree/main/examples/deerfold.ipynb).

### Unconstrained prediction

For comparison DEERFold can accept empty csv file to genreate unconstrained models. 

To unconditionally generate models from DEERFold, run the following script:

``` 
python predict_with_sp_batch.py <fasta_file> <msa_dir> <out_dir> --model <strand/helix/com> --neff neff --num num

E.g. Generate 100 unconstrained models from DEERFold_helix, set MSA Neff as 5
python predict_with_sp_batch.py examples/PfMATE.fasta examples/alignments out --model strand --neff 5 --num 100
```

An example of unconstrained model generation from DEERFold can be found in
[this notebook](https://github.com/CAAIPD/DEERFold/tree/main/examples/deerfold.ipynb).

## Training
### MSA Dataset
To train the model, you will first need to download the Multiple Sequence Alignments (MSAs) are from the [OpenFold dataset](https://registry.opendata.aws/openfold/), 
which contains 401,381 MSAs for 140,000 unique Protein Data Bank (PDB) chains and 16,000,000 UniClust30 clusters.

### Protein Structure Dataset
The correspond protein strucures are extracted from [OpenFold dataset](https://registry.opendata.aws/openfold/) and denoted as `mmcif_dir`.

### DEERFold Training Data
The filenames for train and validation DEERFold splits are saved in `data/train.json` and `data/val.json`

### DEER Data Generation
DEER data used for training can be downloaded [here](). Or you can generate the simulation data in csv file by chiLife for the protein structure in the format of mmCIF file, which the labeling sites are selected by the method [Kazmier et al.](https://www.sciencedirect.com/science/article/abs/pii/S1047847710003400?via%3Dihub). 

**Step 1** Generate the simulated DEER data in csv file from the protein structure in the format of mmCIF file.
```
python scripts/generate_DEER_data.py <mmcif_file> <chain_id> <out_dir>
```
**Step 2** From Step 1, you 'll get the simulated DEER data in csv file and a json file with sequence of corresponding `mmcif_file`. If your input sequence in fasta file doesn't align with the squence in `mmcif_file`, you can use the following script to generate the mapping file and output the reindexed csv file.
```
python scripts/fasta_align.py <json_file> <fasta_file> <out_dir>
```

### Loading Pretrained Alphafold2 Models
The model was fine-tuned using AlphaFold2's pretrained model weights without using templates ([model_5_ptm](https://storage.googleapis.com/alphafold/alphafold_params_2022-01-19.tar)), which was previously released versions of AlphaFold and AlphaFold-Multimer were trained using PDB structures with a release date before 2018-04-30, a cutoff date chosen to coincide with the start of the 2018 CASP13 assessment. 

### Finally, call the training sript:

The training procedure mainly follows the Openfold training procedure. For more details, please check [Openfold](https://github.com/aqlaboratory/openfold/blob/main/docs/source/original_readme.md) for instructions. 
```
python train_openfold.py <train_mmcif_dir> <train_alignment_dir>  <out_dir> --val_data_dir <val_mmcif_dir> --val_alignment_dir <val_alignment_dir> --checkpoint_best_val=True  --precision bf16 --gpus 4 --replace_sampler_ddp=True --seed 4242022 --obsolete_pdbs_file_path obsolete.dat  --train_mapping_path data/train.json --deepspeed_config_path deepspeed_config.json --wandb --wandb_project deerfold
```

