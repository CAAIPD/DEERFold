# dssp_secondary_structure_analyzer.py

import re, os, json
import sys
from itertools import combinations
from collections import defaultdict
import math
import random
from geometry_utils import find_min_angle_pair
import numpy as np
from Bio.PDB import MMCIFParser, PDBIO, Model, MMCIFIO
from Bio.PDB.DSSP import DSSP
import chilife as xl
import MDAnalysis as mda
from io import StringIO
import pandas as pd
from matplotlib import pyplot as plt

def model_to_mdanalysis(filename):
    mdaload = mda.Universe(filename, format="PDB")
    return mdaload

def parse_dssp_file(mmcif_file, filename, chain_id, out_dir):
    m = MMCIFParser()
    structure = m.get_structure("protein", mmcif_file)
    new_structure = structure.__class__(structure.id)
    found = False
    
    for model in structure:
        if found:
            break
        if chain_id in [chain.id for chain in model]:
            new_model = Model.Model(0)
            for chain in model:
                if chain.id == chain_id:
                    if chain_id != 'X':
                        chains_to_remove = [chain for chain in model if chain.id == 'X']
                        if chains_to_remove:
                            model.detach_child('X')
                    # print(chain.id)
                    chain.id = 'X'
                    chain_id = 'X'
                    new_model.add(chain)
                    new_structure.add(new_model)
                    found = True
                    break
    
    if len(new_structure) == 0:
        raise ValueError(f"Chain {chain_id} not found in any model of the structure.")

    # Save the new structure as a PDB file
    io = PDBIO()
    output_pdb_file = f'{out_dir}/{filename}_{chain_id}.pdb'
    io.set_structure(new_structure)
    io.save(output_pdb_file)
    
    target_model = new_structure[0]
    dssp = DSSP(target_model, output_pdb_file, dssp='mkdssp')
    
    secondary_structure = ''
    seq = ''
    full_res = []
    residue_numbers = []
    accessibility = []
    coordinates = {}

    for key in dssp.keys():
        if key[0] == chain_id:
            res_num = key[1][1]
            seq += dssp[key][1]
            full_res.append(res_num)
            ss = dssp[key][2]
            # print(key, dssp[key])
            # Avoid acc 'NA' or residue number not in pdb
            if not isinstance(dssp[key][3], float) or ((' ', res_num, ' ') not in target_model[chain_id]):
                continue

            if ss == 'H':
                secondary_structure += 'H'
            elif ss == 'E':
                secondary_structure += 'E'
            else:
                secondary_structure += 'C'

            residue_numbers.append(res_num)
            accessibility.append(dssp[key][3])
            coordinates[res_num] = target_model[chain_id][res_num]['CA'].coord

    return output_pdb_file, seq, full_res, secondary_structure, residue_numbers, accessibility, coordinates

def find_secondary_structure_elements(ss_string, residue_numbers, accessibility):
    elements = []
    current_type = ss_string[0]
    start = residue_numbers[0]
    start_index = 0
    
    for i in range(1, len(ss_string)):
        if ss_string[i] != current_type:
            elements.append((current_type, start, residue_numbers[i-1], accessibility[start_index:i]))
            current_type = ss_string[i]
            start = residue_numbers[i]
            start_index = i
    
    elements.append((current_type, start, residue_numbers[-1], accessibility[start_index:]))
    return elements

def determine_direction(start, end):
    return "N to C" if start < end else "C to N"

def analyze_surface_exposure(accessibility, residue_start, threshold=20):
    exposed = [residue_start + i for i, acc in enumerate(accessibility) if acc > threshold]
    if not exposed:
        return "No exposure", []
    elif len(exposed) == len(accessibility):
        return "Fully exposed", exposed
    else:
        exposed_ranges = []
        start = exposed[0]
        for i in range(1, len(exposed)):
            if exposed[i] - exposed[i-1] > 1:
                exposed_ranges.append((start, exposed[i-1]))
                start = exposed[i]
        exposed_ranges.append((start, exposed[-1]))
        return ", ".join([f"{s}-{e}" for s, e in exposed_ranges]), exposed

def select_two_residues(exposed_residues):
    if len(exposed_residues) < 2:
        return exposed_residues
    first = random.choice(exposed_residues[:len(exposed_residues)//2])
    second = random.choice(exposed_residues[len(exposed_residues)//2:])
    return [first, second]

def calculate_distance(coord1, coord2):
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(coord1, coord2)))

def calculate_segment_distances(elements, coordinates, helix_cufoff = 5, strand_cutoff = 5):
    distances = defaultdict(list)
    pairs = []
    for ss_type in ['H', 'E']:
        same_type_elements = [e for e in elements if e[0] == ss_type]
        if ss_type == 'H':
            thre = helix_cufoff
        if ss_type == 'E':
            thre = strand_cutoff

        for (e1, e2) in combinations(same_type_elements, 2):
            if (e1[2] - e1[1] + 1 < thre) or (e2[2] - e2[1] + 1 < thre):
                continue
            e_index = {1:e1[1], 2:e1[2], 3:e2[1], 4:e2[2]}
            min_pair = find_min_angle_pair(np.array(coordinates[e1[1]]), np.array(coordinates[e1[2]]), np.array(coordinates[e2[1]]), np.array(coordinates[e2[2]]))
            result = [e_index[int(i)] for i in min_pair.split(',')]
            pairs.append((result[0], result[1]))
            pairs.append((result[2], result[3]))
            distance1 = calculate_distance(coordinates[result[0]], coordinates[result[1]])
            distance2 = calculate_distance(coordinates[result[2]], coordinates[result[3]])
            distances[f"{ss_type}: {e1[1]}-{e1[2]}"].append(f"{e2[1]}-{e2[2]}: {result[0]}-{result[1]}: {distance1:.2f}, {result[2]}-{result[3]}: {distance2:.2f}")
    return distances, pairs

def analyze_dssp_file(mmcif_file, chain_id, helix_cufoff, strand_cutoff, out_dir):
    filename = os.path.splitext(os.path.basename(mmcif_file))[0]
    target_model, seq, full_res, secondary_structure, residue_numbers, accessibility, coordinates = parse_dssp_file(mmcif_file, filename, chain_id, out_dir)
    # print(seq, full_res)
    seq_res = {
    "seq": seq,
    "res_num": full_res,
    }
    with open(f"{out_dir}/{filename}_{chain_id}.json", 'w') as fp:
        json.dump(seq_res, fp)

    elements = find_secondary_structure_elements(secondary_structure, residue_numbers, accessibility)

    for i, (ss_type, start, end, acc) in enumerate(elements):
        if ss_type in ['H', 'E']:
            direction = determine_direction(start, end)
            exposure, exposed_residues = analyze_surface_exposure(acc, start, 0.4)
            elements[i] = (ss_type, start, end, acc, exposed_residues)
            # print(f"{ss_type}: {start}-{end}, Direction: {direction}, Exposed residues: {exposure}")

    distances, pairs = calculate_segment_distances(elements, coordinates, int(helix_cufoff), int(strand_cutoff))
    # print("\nSegment Distances:")
    # for segments, distance in distances.items():
    #     print(f"{segments}: {distance}")
    print(pairs)
    mdaload = model_to_mdanalysis(target_model)
    data_list = []
    r = np.linspace(1, 100, 100)
    if target_model != out_dir+"/"+filename+"_"+chain_id+".pdb":
        new_chain_id = 'X'
    else:
        new_chain_id = chain_id

    for i, j in pairs:
        try:
            P = xl.distance_distribution(xl.SpinLabel.from_mmm('R1M', site=i, chain=new_chain_id, protein=mdaload), xl.SpinLabel.from_mmm('R1M', site=j, chain=new_chain_id, protein=mdaload), r=r)
            # plt.plot(r, P)
            # plt.savefig(f'{i}-{j}.png')
            # plt.close()
        except:
            continue

        a = np.argmax(P)
        maxValue = r[a]
        if maxValue > 15 and maxValue < 100:
            data_list.append([i, j] + P.tolist())
    output_df = pd.DataFrame(data_list)
    output_df.to_csv(f"{out_dir}/{filename}_{chain_id}.csv", index=None, header=False)
    if target_model != out_dir+"/"+filename+"_"+chain_id+".pdb":
        os.system(f"mv {out_dir}/{filename}_X.pdb {out_dir}/{filename}_{chain_id}.pdb")

if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage: python generate_DEER_data.py <mmcif_file_path> <chain_id> <helix_cufoff> <strand_cutoff> <out_dir>")
        sys.exit(1)
    
    mmcif_file_path = sys.argv[1]
    chain_id = sys.argv[2]
    helix_cufoff = sys.argv[3]
    strand_cutoff = sys.argv[4]
    out_dir = sys.argv[5]

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    try:
        analyze_dssp_file(mmcif_file_path, chain_id, helix_cufoff, strand_cutoff, out_dir)
    except Exception as e:
        print(f"An error occurred: {e}")
        sys.exit(1)