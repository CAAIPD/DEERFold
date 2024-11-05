import argparse
from copy import deepcopy
from datetime import date
import logging
import math
import numpy as np
import os

logging.basicConfig()
logger = logging.getLogger(__file__)
logger.setLevel(level=logging.INFO)

import pickle
import random
import sys
import time
import torch
import re
from Bio import SeqIO
import pandas as pd

torch_versions = torch.__version__.split(".")
torch_major_version = int(torch_versions[0])
torch_minor_version = int(torch_versions[1])
if(
    torch_major_version > 1 or 
    (torch_major_version == 1 and torch_minor_version >= 12)
):
    # Gives a large speedup on Ampere-class GPUs
    torch.set_float32_matmul_precision("high")

torch.set_grad_enabled(False)

from openfold.config import model_config
from openfold.data import feature_pipeline, data_pipeline
from openfold.model.model import AlphaFold
from openfold.np import residue_constants, protein
import openfold.np.relax.relax as relax
from openfold.utils.tensor_utils import (
    tensor_tree_map,
)
from openfold.data.msa_subsampling import subsample_msa_sequentially

# Global configuration class
class Config:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

# Global configuration instance
config = Config()

def precompute_alignments(tag, seq, alignment_dir):
    tmp_fasta_path = os.path.join(config.output_dir, f"tmp_{os.getpid()}.fasta")
    with open(tmp_fasta_path, "w") as fp:
        fp.write(f">{tag}\n{seq}")

    local_alignment_dir = os.path.join(alignment_dir, tag)

    if(config.use_precomputed_alignments is None and not os.path.isdir(local_alignment_dir)):
        logger.info(f"Generating alignments for {tag}...")
            
        os.makedirs(local_alignment_dir)

        alignment_runner = data_pipeline.AlignmentRunner(
            jackhmmer_binary_path=config.jackhmmer_binary_path,
            hhblits_binary_path=config.hhblits_binary_path,
            hhsearch_binary_path=config.hhsearch_binary_path,
            uniref90_database_path=config.uniref90_database_path,
            mgnify_database_path=config.mgnify_database_path,
            bfd_database_path=config.bfd_database_path,
            uniclust30_database_path=config.uniclust30_database_path,
            pdb70_database_path=config.pdb70_database_path,
            no_cpus=config.cpus,
        )
        alignment_runner.run(
            tmp_fasta_path, local_alignment_dir
        )
    else:
        logger.info(
            f"Using precomputed alignments for {tag} at {alignment_dir}..."
        )

    # Remove temporary FASTA file
    os.remove(tmp_fasta_path)

def run_model(model, batch, tag):
    with torch.no_grad(): 
        # Disable templates if there aren't any in the batch
        model.config.template.enabled = model.config.template.enabled and any([
            "template_" in k for k in batch
        ])

        logger.info(f"Running inference for {tag}...")
        t = time.perf_counter()
        out = model(batch)
        inference_time = time.perf_counter() - t
        logger.info(f"Inference time: {inference_time}")
   
    return out

def parse_fasta(data):
    data = re.sub('>$', '', data, flags=re.M)
    lines = [
        l.replace('\n', '')
        for prot in data.split('>') for l in prot.strip().split('\n', 1)
    ][1:]
    tags, seqs = lines[::2], lines[1::2]

    tags = [t.split()[0] for t in tags]

    return tags, seqs

def generate_feature_dict(
    tags,
    seqs,
    alignment_dir,
    data_processor
):
    tmp_fasta_path = os.path.join(config.output_dir, f"tmp_{os.getpid()}.fasta")
    if len(seqs) == 1:
        tag = tags[0]
        seq = seqs[0]
        with open(tmp_fasta_path, "w") as fp:
            fp.write(f">{tag}\n{seq}")

        local_alignment_dir = os.path.join(alignment_dir, tag)
        feature_dict = data_processor.process_fasta(
            fasta_path=tmp_fasta_path, alignment_dir=local_alignment_dir
        )
    else:
        with open(tmp_fasta_path, "w") as fp:
            fp.write(
                '\n'.join([f">{tag}\n{seq}" for tag, seq in zip(tags, seqs)])
            )
        feature_dict = data_processor.process_multiseq_fasta(
            fasta_path=tmp_fasta_path, super_alignment_dir=alignment_dir,
        )

    # Remove temporary FASTA file
    os.remove(tmp_fasta_path)

    return feature_dict

def get_model_basename(model_path):
    return os.path.splitext(
                os.path.basename(
                    os.path.normpath(model_path)
                )
            )[0]

def make_output_directory(output_dir, model_name, multiple_model_mode):
    if multiple_model_mode:
        prediction_dir = os.path.join(output_dir, model_name)
    else:
        prediction_dir = os.path.join(output_dir)
    os.makedirs(prediction_dir, exist_ok=True)
    return prediction_dir

def count_models_to_evaluate(openfold_checkpoint_path, jax_param_path):
    model_count = 0
    if openfold_checkpoint_path:
        model_count += len(openfold_checkpoint_path.split(","))
    if jax_param_path:
        model_count += len(jax_param_path.split(","))
    return model_count

def load_model():
    # Create the output directory
    model = AlphaFold(config.model_config)
    model = model.eval()
    checkpoint_basename = get_model_basename(config.checkpoint_path)

    sd = torch.load(config.checkpoint_path)['ema']['params']
    new_sd = {key.replace('module.model.', ''): value for key, value in sd.items()}

    model.load_state_dict(new_sd)
    
    model = model.to(config.model_device)
    logger.info(
        f"Loaded OpenFold parameters at {config.checkpoint_path}..."
    )
    return model

def list_files_with_extensions(dir, extensions):
    return [f for f in os.listdir(dir) if f.endswith(extensions)]


def fill_org_dist(fasta_file, sp_csv):
    record = SeqIO.read(fasta_file, 'fasta')
    n = len(record.seq)
    org_dist = np.zeros((n, n, 100))

    if os.path.getsize(sp_csv) != 0:
        df = pd.read_csv(sp_csv, header=None)

        for _, row in df.iterrows():
            i, j = int(row[0]), int(row[1])
            org_dist[i - 1][j - 1] = org_dist[j - 1][i - 1] = row[2:].to_numpy()
    
    return org_dist

def run_deerfold(**kwargs):
    """
    Main function to run the DEERFold pipeline.
    Updates the global config with provided kwargs and executes the pipeline.
    """
    # Update global configuration
    config.__dict__.update(kwargs)
    
    os.makedirs(config.output_dir, exist_ok=True)

    config.model_config = model_config('model_5_ptm')
    data_processor = data_pipeline.DataPipeline(
        template_featurizer=None,
    )

    random_seed = config.data_random_seed
    if random_seed is None:
        random_seed = random.randrange(2**32)
    
    np.random.seed(random_seed)
    torch.manual_seed(random_seed + 1)
    
    feature_processor = feature_pipeline.FeaturePipeline(config.model_config.data)

    if config.use_precomputed_alignments is None:
        alignment_dir = os.path.join(config.output_dir, "alignments")
    else:
        alignment_dir = config.use_precomputed_alignments


    with open(config.fasta_path, "r") as fp:
        data = fp.read()

    tag, seq = parse_fasta(data)
    tag = tag[0]
    seq = seq[0]

    model = load_model()
    output_directory = make_output_directory(config.output_dir, get_model_basename(config.checkpoint_path), False)

    output_name = f'{tag}_model_1'

    if config.features:
        feature_dict = pickle.load(open(config.features,'rb'))
    else:
        # Does nothing if the alignments have already been computed
        precompute_alignments(tag, seq, alignment_dir)
        feature_dict = generate_feature_dict(
            [tag],
            [seq],
            alignment_dir,
            data_processor,
        )

    if config.splabel.endswith('.csv'):
        sp = fill_org_dist(config.fasta_path, config.splabel)
        feature_dict['sp'] = sp
    else:
        logger.error("Input distance info need to be in CSV file")
        return

    # subsample MSAs to specified Neff
    msa = feature_dict['msa']

    if config.neff:
        logger.info(
            f"Subsampling MSA to Neff={config.neff}..."
        )
        indices = subsample_msa_sequentially(msa, neff=config.neff, eff_cutoff=0.62, cap_msa=False)
        feature_dict['msa'] = msa[indices]
        feature_dict['deletion_matrix_int'] = feature_dict['deletion_matrix_int'][indices]
#     print(feature_dict['msa'].shape, feature_dict['deletion_matrix_int'].shape)

    processed_feature_dict = feature_processor.process_features(
        feature_dict, mode='predict',
    )

    processed_feature_dict = {
        k:torch.as_tensor(v, device=config.model_device) 
        for k,v in processed_feature_dict.items()
    }

    out = run_model(model, processed_feature_dict, tag)

    # Toss out the recycling dimensions --- we don't need them anymore
    processed_feature_dict = tensor_tree_map(
        lambda x: np.array(x[..., -1].cpu()), 
        processed_feature_dict
    )
    out = tensor_tree_map(lambda x: np.array(x.cpu()), out)

    plddt = out["plddt"]

    plddt_b_factors = np.repeat(
        plddt[..., None], residue_constants.atom_type_num, axis=-1
    )

    unrelaxed_protein = protein.from_prediction(
        features=processed_feature_dict,
        result=out,
        b_factors=plddt_b_factors
    )

    unrelaxed_output_path = os.path.join(
        output_directory, f'{output_name}_unrelaxed.pdb'
    )


    with open(unrelaxed_output_path, 'w') as fp:
        fp.write(protein.to_pdb(unrelaxed_protein))

    logger.info(f"Output written to {unrelaxed_output_path}...")

    if not config.skip_relaxation:    
        amber_relaxer = relax.AmberRelaxation(
            **config.relax
        )
        
        relaxed_pdb_str, _, _ = amber_relaxer.process(prot=unrelaxed_protein)

        # Save the relaxed PDB.
        relaxed_output_path = os.path.join(
            output_directory, f'{output_name}_relaxed.pdb'
        )
        with open(relaxed_output_path, 'w') as fp:
            fp.write(relaxed_pdb_str)
        
        logger.info(f"Relaxed output written to {relaxed_output_path}...")

    if config.save_outputs:
        output_dict_path = os.path.join(
            output_directory, f'{output_name}_output_dict.pkl'
        )
        with open(output_dict_path, "wb") as fp:
            pickle.dump(feature_dict['msa'], fp, protocol=pickle.HIGHEST_PROTOCOL)

        logger.info(f"Model output written to {output_dict_path}...")
    
    return output_directory, output_name