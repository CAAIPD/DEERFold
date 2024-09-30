import sys
import os
import numpy as np
import pandas as pd
import subprocess
from Bio import SeqIO
import MDAnalysis
import chilife as xl
from scipy.stats import wasserstein_distance

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

def main(args):
    fasta_file, sp_csv, msa_dir, neff, ref_pdbs, outdir = args

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    pdbid, _ = os.path.splitext(os.path.basename(fasta_file))
    num = 5

    pairs = []
    input_csv = {}
    if os.path.getsize(sp_csv) == 0:
        print(f"# of distance constraints: 0")
    else:
        df = pd.read_csv(sp_csv, header=None)
        count = 0
        residues = set()
        for _, row in df.iterrows():
            i, j = int(row[0]), int(row[1])
            input_csv[(i, j)] = np.array(row[2:])
            pairs.append((i,j))
            count += 1
            residues.add(i)
            residues.add(j)
        print(f"{sp_csv} # of distance constraints: {count}, residues are {residues}")

    model_files = [
        "model/DEERFold_helix.pt",
        "model/DEERFold_strand.pt",
        "model/DEERFold_com.pt"
    ]

    for model_idx, model_file in enumerate(model_files):
        for i in range(num):
            model_num = model_idx * num + i + 1
            if not os.path.exists(f"{outdir}/{pdbid}_model_{model_num}.pdb"):
                os.system(f"python predict_with_sp.py {fasta_file} {sp_csv} --use_precomputed_alignments {msa_dir} --distograms --output_dir {outdir}/{pdbid} --checkpoint_path {model_file} --neff {neff} --skip_relaxation")
                os.system(f"mv {outdir}/{pdbid}/{pdbid}_model_1_unrelaxed.pdb {outdir}/{pdbid}_model_{model_num}.pdb")
                os.system(f"rm -r {outdir}/{pdbid}")
            else:
                print(f"{outdir}/{pdbid}_model_{model_num}.pdb existed....Skip")

    num = num * 3
    if pairs:
        results2 = {}
        for i in range(num):
            emd_dist = calculate_emd(f"{outdir}/{pdbid}_model_{i + 1}.pdb", input_csv, pairs)
            results2[f"{pdbid}_model_{i + 1}.pdb"] = emd_dist
            print(f"Calculating emd scores for {pdbid}_model_{i + 1}.pdb....Done")

        sorted_result = sorted(results2.items(), key=lambda x: x[1])
        top_n = sorted_result[:num]
        print(top_n)

        rename_operations = []
        for i, (old_name, _) in enumerate(top_n, start=1):
            temp_name = f"{outdir}/temp_{i}.pdb"
            final_name = f"{outdir}/{pdbid}_model_{i}.pdb"
            rename_operations.append((f"{outdir}/{old_name}", temp_name, final_name))

        for old_name, temp_name, final_name in rename_operations:
            os.rename(old_name, temp_name)

        for _, temp_name, final_name in rename_operations:
            os.rename(temp_name, final_name)

    for ref_pdb in ref_pdbs.split(','):
        results2 = {}
        tms = {}
        for i in range(num):
            if not os.path.exists(f"{outdir}/{pdbid}_model_{i + 1}.pdb"):
                continue
            rmsd = get_rmsd(f"{outdir}/{pdbid}_model_{i + 1}.pdb", f"{ref_pdb}")
            results2[f"{pdbid}_model_{i + 1}.pdb"] = rmsd
            tmscore = get_tmscore(f"{outdir}/{pdbid}_model_{i + 1}.pdb", f"{ref_pdb}")
            tms[f"{pdbid}_model_{i + 1}.pdb"] = tmscore

        mean_value = sum(results2.values()) / len(results2)
        min_value = min(results2.values())
        min_key = [key for key, value in results2.items() if value == min_value]
        print(results2)
        print(f"RMSD between prediction and {ref_pdb}: {mean_value}")
        print(f"Best RMSD between prediction({min_key}) and {ref_pdb}: {min_value}")
        
        mean_value = sum(tms.values()) / len(tms)
        max_value = max(tms.values())
        max_key = [key for key, value in tms.items() if value == max_value]
        print(tms)
        print(f"TMscore between prediction and {ref_pdb}: {mean_value}")
        print(f"Best TMscore between prediction({max_key}) and {ref_pdb}: {max_value}")

if __name__ == "__main__":
    main(sys.argv[1:])