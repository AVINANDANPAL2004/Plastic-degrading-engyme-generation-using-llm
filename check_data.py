from Bio import SeqIO
import gzip

with gzip.open("data/uniprot_datasets.gz", "rt") as handle:
    seqs = list(SeqIO.parse(handle, "fasta"))

print("Total sequences:", len(seqs))
lengths = [len(s.seq) for s in seqs]
print("Min length:", min(lengths))
print("Max length:", max(lengths))
print("Sample sequence:", seqs[0].seq[:100])
