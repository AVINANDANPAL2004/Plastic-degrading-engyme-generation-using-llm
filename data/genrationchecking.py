from Bio import SeqIO

# Load training sequences
train_seqs = [str(r.seq) for r in SeqIO.parse("data/merged_clean_enzymes.fasta", "fasta")]
train_set = set(train_seqs)

# Load generated sequences
generated = []
with open("generated_enzymes.fasta") as f:
    seq = ""
    for line in f:
        if line.startswith(">"):
            if seq:
                generated.append(seq)
                seq = ""
        else:
            seq += line.strip()
    if seq:
        generated.append(seq)

# Check duplicates
for i, seq in enumerate(generated):
    if seq in train_set:
        print(f"Generated sequence {i+1} is EXACT COPY of training data")
    else:
        print(f"Generated sequence {i+1} is NEW")