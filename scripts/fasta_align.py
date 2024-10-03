# File: pdb_fasta_align_and_map.py

import json, sys, os
from Bio import SeqIO
from Bio import pairwise2
from Bio.Seq import Seq

def read_json(file_path):
    with open(file_path, 'r') as f:
        data = json.load(f)
    return data['seq'], data['res_num']

def read_fasta(file_path):
    with open(file_path, 'r') as f:
        record = next(SeqIO.parse(f, 'fasta'))
    return str(record.seq)

def create_mapping(aligned_json_seq, aligned_fasta_seq, res_num):
    mapping = {}
    json_index = 0
    fasta_index = 0

    for i, (json_char, fasta_char) in enumerate(zip(aligned_json_seq, aligned_fasta_seq)):
        if json_char != '-' and fasta_char != '-':
            mapping[res_num[json_index]] = fasta_index
            json_index += 1
            fasta_index += 1
        elif json_char != '-':
            json_index += 1
        elif fasta_char != '-':
            fasta_index += 1
    
    return mapping

jsonfile = sys.argv[1] #'raw/pdb_seq/2pv7_A.json'
fasta = sys.argv[2] #'raw/fasta/2pv7_A.fasta'

id = os.path.splitext(os.path.basename(jsonfile))[0]

# Read input files
json_seq, res_num = read_json(jsonfile)
fasta_seq = read_fasta(fasta)

# Perform sequence alignment
alignments = pairwise2.align.globalms(fasta_seq, json_seq, 2, -1, -0.5, -0.1)
best_alignment = alignments[0]
aligned_fasta_seq, aligned_json_seq, score, start, end = best_alignment

# Create mapping
mapping = create_mapping(aligned_json_seq, aligned_fasta_seq, res_num)

# Print results
print(f"Aligned FASTA sequence: {aligned_fasta_seq}")
print(f"Aligned JSON sequence:  {aligned_json_seq}")
print(f"Alignment score: {score}")

print("\nMapping (res_num to FASTA position):")
for res, pos in mapping.items():
    print(f"Residue {res} -> FASTA position {pos + 1}")  # +1 for 1-based indexing

# Optionally, save the mapping to a file
with open(f'{id}_mapping.json', 'w') as f:
    json.dump({str(k): v + 1 for k, v in mapping.items()}, f, indent=2)  # +1 for 1-based indexing

print("\nMapping saved to 'residue_mapping.json'")

# Calculate and print statistics
total_res = len(res_num)
mapped_res = len(mapping)
print(f"\nTotal residues in JSON: {total_res}")
print(f"Mapped residues: {mapped_res}")
print(f"Mapping coverage: {mapped_res/total_res*100:.2f}%")

# Find unmapped residues
unmapped = set(res_num) - set(mapping.keys())
if unmapped:
    print("\nUnmapped residues:")
    print(sorted(unmapped))