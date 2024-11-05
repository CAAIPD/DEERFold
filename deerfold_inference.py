import sys
import os
import csv
import argparse
import numpy as np
import pandas as pd
import subprocess
from Bio import SeqIO
import MDAnalysis
import chilife as xl
from scipy.stats import wasserstein_distance
from deerfold_processor import run_deerfold

def cdf(arr):
    return np.cumsum(arr) / np.sum(arr)

def calculate_emd(pdb_file_1, input_csv, pairs):
    emd = []
    sigma = 2
    r = np.linspace(1, 100, 100)
    mdaload1 = MDAnalysis.Universe(pdb_file_1)
    for i, j in pairs:
        P = xl.distance_distribution(xl.SpinLabel.from_mmm('R1M', site=i, chain='A', protein=mdaload1), xl.SpinLabel.from_mmm('R1M', site=j, chain='A', protein=mdaload1), r=r)
        a = np.argmax(P)
        maxValue = r[a]
        gauss = 1/np.sqrt(2*np.pi*sigma**2)*np.exp(-(r-maxValue)**2/(2*sigma**2))
        dist1 = gauss/np.sum(gauss)
        peak_index = np.argmax(dist1)
        left_index = max(0, peak_index - 8)
        right_index = min(len(dist1) - 1, peak_index + 8)
        dist1 = dist1[left_index:right_index + 1]

        dist2 = input_csv[(i, j)]
        emd.append(wasserstein_distance(cdf(np.asarray(P)), cdf(dist2)))
    return np.mean(emd)

def get_rmsd(pdb1, pdb2, tmscore_path="./TMscore"):
    result = subprocess.run([tmscore_path, pdb1, pdb2], stdout=subprocess.PIPE)
    output = result.stdout.decode('utf-8')

    for line in output.split("\n"):
        if "RMSD of  the common residues" in line and "=" in line:
            return float(line.split("=")[1].split()[0])
    return None

def get_tmscore(pdb1, pdb2, tmscore_path="./TMscore"):
    result = subprocess.run([tmscore_path, pdb1, pdb2], stdout=subprocess.PIPE)
    output = result.stdout.decode('utf-8')

    for line in output.split("\n"):
        if "TM-score" in line and "=" in line:
            return float(line.split("=")[1].split()[0])
    return None

def create_empty_csv(file_path):
    # Ensure the directory exists
    os.makedirs(os.path.dirname(file_path), exist_ok=True)
    
    # Create an empty CSV file
    with open(file_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
    print(f"Empty CSV file created at: {file_path}")


def main(args):
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    model_files = args.models.split(',')
    pdbid, _ = os.path.splitext(os.path.basename(args.fasta))
    num_per_model = args.num_models // len(model_files)  # Divide equally among 3 models
    remainder = args.num_models % len(model_files)  # Handle any remainder

    pairs = []
    input_csv = {}
    if args.splabel is None or os.path.getsize(args.splabel) == 0:
        print(f"# of distance constraints: 0")
        args.splabel = f"{args.outdir}/empty.csv"
        create_empty_csv(args.splabel)
    else:
        df = pd.read_csv(args.splabel, header=None)
        count = 0
        residues = set()
        for _, row in df.iterrows():
            i, j = int(row[0]), int(row[1])
            input_csv[(i, j)] = np.array(row[2:])
            pairs.append((i,j))
            count += 1
            residues.add(i)
            residues.add(j)
        print(f"{args.splabel} # of distance constraints: {count}, residues are {residues}")

    model_num = 1
    for model_idx, model_file in enumerate(model_files):
        num_for_this_model = num_per_model + (1 if model_idx < remainder else 0)
        for i in range(num_for_this_model):
            output_path = f"{args.outdir}/{pdbid}_model_{model_num}.pdb"
            if not os.path.exists(output_path):
                output_dir, output_name = run_deerfold(
                    fasta_path=args.fasta,
                    splabel=args.splabel,
                    use_precomputed_alignments=args.msa_dir,
                    output_dir=f"{args.outdir}/{pdbid}",
                    checkpoint_path=model_file,
                    neff=args.neff,
                    # neff=None,
                    distograms=True,
                    model_device=args.model_device,
                    features=args.features,
                    data_random_seed=None,
                    skip_relaxation=True,
                    save_outputs=True
                )
                os.rename(f"{output_dir}/{output_name}_unrelaxed.pdb", output_path)
                os.system(f"rm -r {output_dir}")
                print(f"{output_path} .......Complete")
            else:
                print(f"{output_path} existed....Skip")
            model_num += 1

    if pairs:
        results2 = {}
        for i in range(args.num_models):
            emd_dist = calculate_emd(f"{args.outdir}/{pdbid}_model_{i + 1}.pdb", input_csv, pairs)
            results2[f"{pdbid}_model_{i + 1}.pdb"] = emd_dist
            print(f"Calculating emd scores for {pdbid}_model_{i + 1}.pdb....Done")

        sorted_result = sorted(results2.items(), key=lambda x: x[1])
        top_n = sorted_result[:args.num_models]
        # print(top_n)

        rename_operations = []
        for i, (old_name, _) in enumerate(top_n, start=1):
            temp_name = f"{args.outdir}/temp_{i}.pdb"
            final_name = f"{args.outdir}/{pdbid}_model_{i}.pdb"
            rename_operations.append((f"{args.outdir}/{old_name}", temp_name, final_name))

        for old_name, temp_name, final_name in rename_operations:
            os.rename(old_name, temp_name)

        for _, temp_name, final_name in rename_operations:
            os.rename(temp_name, final_name)
    
    print("\n")
    for ref_pdb in args.ref_pdbs.split(','):
        results2 = {}
        tms = {}
        for i in range(args.num_models):
            if not os.path.exists(f"{args.outdir}/{pdbid}_model_{i + 1}.pdb"):
                continue
            rmsd = get_rmsd(f"{args.outdir}/{pdbid}_model_{i + 1}.pdb", f"{ref_pdb}")
            results2[f"{pdbid}_model_{i + 1}.pdb"] = rmsd
            tmscore = get_tmscore(f"{args.outdir}/{pdbid}_model_{i + 1}.pdb", f"{ref_pdb}")
            tms[f"{pdbid}_model_{i + 1}.pdb"] = tmscore

        mean_value = sum(results2.values()) / len(results2)
        min_value = min(results2.values())
        min_key = [key for key, value in results2.items() if value == min_value]
        print(f"Average RMSD between prediction and {ref_pdb}: {mean_value}")
        print(f"Best RMSD between prediction({min_key}) and {ref_pdb}: {min_value}")
        print(results2, "\n")
        
        mean_value = sum(tms.values()) / len(tms)
        max_value = max(tms.values())
        max_key = [key for key, value in tms.items() if value == max_value]
        print(f"Average TMscore between prediction and {ref_pdb}: {mean_value}")
        print(f"Best TMscore between prediction({max_key}) and {ref_pdb}: {max_value}")
        print(tms, "\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="DEERFold processing script")
    parser.add_argument(
        "fasta", type=str,
        help="Path to FASTA file, one sequence per file"
    )
    parser.add_argument(
        "msa_dir", type=str,
        help="Directory containing precomputed MSA"
    )
    parser.add_argument(
        "outdir", type=str,
        help="Output directory for results"
    )
    parser.add_argument(
        "--splabel", type=str,
        help="Path to the spin label CSV file"
    )
    parser.add_argument(
        "--neff", type=float, default=5.0,
        help="Neff value for MSA subsampling"
    )
    parser.add_argument(
        "--models", type=str, default="models/DEERFold_helix.pt,models/DEERFold_strand.pt,models/DEERFold_com.pt",
        help="Model weights (default: path/to/default/model.pth)"
    )
    parser.add_argument(
        "--save_outputs", action="store_true", default=True,
        help="Whether to save all model outputs, including embeddings, etc."
    )
    parser.add_argument(
        "--features", type=str,
        help="Feature pickle"
    )
    parser.add_argument(
        "--model_device", type=str, default="cuda:0",
        help="""Name of the device on which to run the model. Any valid torch
             device name is accepted (e.g. "cpu", "cuda:0")"""
    )
    parser.add_argument(
        "--ref_pdbs", type=str,
        help="Comma-separated list of reference PDB files"
    )
    parser.add_argument(
        "--num_models", type=int, default=15,
        help="Total number of models to generate (default: 15)"
    )

    args = parser.parse_args()
    main(args)