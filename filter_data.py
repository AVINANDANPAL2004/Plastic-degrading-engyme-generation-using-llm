from Bio import SeqIO
import gzip

filtered = []

with gzip.open("data/uniprot_datasets.gz", "rt") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        if 250 <= len(record.seq) <= 350:
            filtered.append(record)

SeqIO.write(filtered, "data/petase_clean.fasta", "fasta")

print("Filtered sequences:", len(filtered))
