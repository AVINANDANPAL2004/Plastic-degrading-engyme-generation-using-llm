from Bio import SeqIO
import gzip
import os

BASE = os.path.dirname(__file__)

file_path = os.path.join(BASE, "cutinase.fasta.gz")

with gzip.open(file_path, "rt") as handle:
    seqs = list(SeqIO.parse(handle, "fasta"))


print("Total cutinase sequences:", len(seqs))

lengths = [len(s.seq) for s in seqs]
print("Min length:", min(lengths))
print("Max length:", max(lengths))

print("Sample sequence:", seqs[0].seq[:60])
