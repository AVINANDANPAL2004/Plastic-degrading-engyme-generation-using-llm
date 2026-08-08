from Bio import SeqIO
import os
# Step 1: Load clean protein sequences
# 🔥 Base path (IMPORTANT)
BASE = os.path.dirname(__file__)
file_path = os.path.join(BASE, "petase_clean.fasta")

# Step 1: Load clean protein sequences
sequences = [str(r.seq) for r in SeqIO.parse(file_path, "fasta")]




print("Total sequences loaded:", len(sequences))

# Step 2: Create vocabulary (unique amino acids)
vocab = sorted(set("".join(sequences)))
vocab_size = len(vocab)

print("Vocabulary size:", vocab_size)
print("Amino acids:", vocab)

# Step 3: Create mappings
stoi = {aa: idx for idx, aa in enumerate(vocab)}
itos = {idx: aa for aa, idx in stoi.items()}

# Step 4: Encode sequences
def encode_sequence(seq):
    return [stoi[aa] for aa in seq]

encoded_sequences = [encode_sequence(seq) for seq in sequences]

# Show example
print("\nExample:")
print("Original sequence:", sequences[0][:30])
print("Encoded sequence:", encoded_sequences[0][:30])
