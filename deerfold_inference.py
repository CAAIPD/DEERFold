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
from scripts.mmseqs_query_colabfold import run_mmseqs2, read_fasta_lists

def cdf(arr):
    return np.cumsum(arr) / np.sum(arr)

def calculate_emd(pdb_file_1, input_csv, pairs, nan_penalty_factor=1, plddt_penalty_factor=0.5):
    mdaload1 = MDAnalysis.Universe(pdb_file_1)
    atoms = mdaload1.select_atoms("chainID A")
    resids = np.unique(atoms.resids)

    pLDDT_dict = {}
    for resid in resids:
        residue_atoms = mdaload1.select_atoms(f"chainID A and resid {resid}")
        pLDDT = np.mean(residue_atoms.tempfactors)
        pLDDT_dict[resid] = pLDDT

    r = np.linspace(1, 100, 100)
    emd_values = []
    total_pairs = len(pairs)

    for i, j in pairs:
        P = xl.distance_distribution(
            xl.SpinLabel.from_mmm('R1M', site=i, chain='A', protein=mdaload1),
            xl.SpinLabel.from_mmm('R1M', site=j, chain='A', protein=mdaload1),
            r=r
        )
        dist2 = input_csv[(i, j)]
        emd_value = wasserstein_distance(cdf(np.asarray(P)), cdf(dist2))

        if not np.isnan(emd_value):
            emd_values.append(emd_value)

    if emd_values:
        mean_emd = np.mean(emd_values)
        nan_ratio = (total_pairs - len(emd_values)) / total_pairs
        score = mean_emd + nan_ratio * nan_penalty_factor
    else:
        score = float('inf')

    return score

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
    os.makedirs(os.path.dirname(file_path), exist_ok=True)
    with open(file_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
    print(f"Empty CSV file created at: {file_path}")

def generate_msa(fasta_path, outdir):
    fasta_fname = os.path.basename(fasta_path).split('.')[0]
    tmp_path = f"{outdir}_tmp_{fasta_fname}/"
    if not os.path.exists(tmp_path):
        os.makedirs(tmp_path)
    names, seqs = read_fasta_lists(fasta_path)
    msas = run_mmseqs2(seqs, prefix=tmp_path, user_agent="DEERFold/1.0")
    for name, msa in zip(names, msas):
        os.makedirs(f"{outdir}/{name}/", exist_ok=True)
        with open(f"{outdir}/{name}/{name}.a3m", "w") as f:
            f.write(msa)
    os.system(f"rm -r {tmp_path}")
    return outdir

def check_msa_exists(fasta_path, outdir):
    """Check if MSA files already exist in outdir based on FASTA headers."""
    names, _ = read_fasta_lists(fasta_path)
    for name in names:
        msa_file = f"{outdir}/{name}/{name}.a3m"
        if not os.path.exists(msa_file):
            return False
        print(f"MSA exists in: {outdir}/{name}/{name}.a3m")
    return True

def main(args):
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    # Handle MSA generation or use provided MSA directory
    if args.msa_dir and os.path.exists(args.msa_dir):
        print(f"Using precomputed MSA from: {args.msa_dir}")
    elif check_msa_exists(args.fasta, args.outdir):
        print(f"Using existing MSA in: {args.outdir}")
        args.msa_dir = args.outdir
    else:
        print("No MSA directory provided or found. Generating MSA...")
        args.msa_dir = generate_msa(args.fasta, args.outdir)

    model_files = args.models.split(',')
    pdbid, _ = os.path.splitext(os.path.basename(args.fasta))
    num_per_model = args.num_models // len(model_files)
    remainder = args.num_models % len(model_files)

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
    if args.ref_pdbs:
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
    parser.add_argument("fasta", type=str, help="Path to FASTA file, one sequence per file")
    parser.add_argument("outdir", type=str, help="Output directory for results")
    parser.add_argument("--splabel", type=str, help="Path to the spin label CSV file")
    parser.add_argument("--msa_dir", type=str, help="Directory containing precomputed MSA (optional)")
    parser.add_argument("--neff", type=float, help="Neff value for MSA subsampling")
    parser.add_argument("--models", type=str, default="model/DEERFold.pt", help="Model weights")
    parser.add_argument("--save_outputs", action="store_true", default=True, help="Whether to save all model outputs")
    parser.add_argument("--features", type=str, help="Feature pickle")
    parser.add_argument("--model_device", type=str, default="cuda:0", help="Device for model (e.g., 'cpu', 'cuda:0')")
    parser.add_argument("--ref_pdbs", type=str, help="Comma-separated list of reference PDB files")
    parser.add_argument("--num_models", type=int, default=15, help="Total number of models to generate")

    args = parser.parse_args()
    main(args)