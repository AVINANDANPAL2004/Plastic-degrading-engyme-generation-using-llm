from Bio import SeqIO
import gzip
import os

VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")

def is_valid(seq):
    return set(seq).issubset(VALID_AA)

# 🔥 Base directory (VERY IMPORTANT)
BASE = os.path.dirname(__file__)

# Correct file paths
petase_path = os.path.join(BASE, "petase_clean.fasta")
cutinase_path = os.path.join(BASE, "cutinase.fasta.gz")
output_path = os.path.join(BASE, "merged_clean_enzymes.fasta")


all_sequences = {}

# ---- Load PETase ----
for record in SeqIO.parse("data/petase_clean.fasta", "fasta"):
    seq = str(record.seq)
    if 200 <= len(seq) <= 400 and is_valid(seq):
        all_sequences[seq] = record

# ---- Load cutinase (.gz) ----
with gzip.open("data/cutinase.fasta.gz", "rt") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        seq = str(record.seq)
        if 200 <= len(seq) <= 400 and is_valid(seq):
            all_sequences[seq] = record

# ---- Save merged clean dataset ----
SeqIO.write(
    all_sequences.values(),
    "merged_clean_enzymes.fasta",
    "fasta"
)

print("Final cleaned dataset size:", len(all_sequences))