from Bio import SeqIO
from difflib import SequenceMatcher

# Load training sequences
train_seqs = [str(r.seq) for r in SeqIO.parse("data/merged_clean_enzymes.fasta", "fasta")]

# Load generated sequences
generated = [str(r.seq) for r in SeqIO.parse("generated_enzymes.fasta", "fasta")]

def similarity(a,b):
    return SequenceMatcher(None,a,b).ratio()

# Compare sequences
for i, g in enumerate(generated):

    max_sim = 0

    for t in train_seqs:
        s = similarity(g,t)

        if s > max_sim:
            max_sim = s

    print(f"Generated {i+1} max similarity: {max_sim:.2f}")