import torch
import torch.nn as nn
from Bio import SeqIO

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

# Load dataset to rebuild vocabulary
seqs = [str(r.seq) for r in SeqIO.parse("data/merged_clean_enzymes.fasta", "fasta")]

vocab = sorted(set("".join(seqs)))
vocab.append("<PAD>")

stoi = {c:i for i,c in enumerate(vocab)}
itos = {i:c for c,i in stoi.items()}

vocab_size = len(vocab)
PAD_IDX = stoi["<PAD>"]

# Model configuration (same as training)
D_MODEL = 512
NUM_LAYERS = 6
N_HEAD = 8
FF_DIM = 2048

class ProteinLLM(nn.Module):
    def __init__(self):
        super().__init__()
        self.embed = nn.Embedding(vocab_size, D_MODEL)
        self.pos = nn.Embedding(512, D_MODEL)

        encoder_layer = nn.TransformerEncoderLayer(
            d_model=D_MODEL,
            nhead=N_HEAD,
            dim_feedforward=FF_DIM,
            dropout=0.1,
            batch_first=True
        )

        self.transformer = nn.TransformerEncoder(
            encoder_layer,
            num_layers=NUM_LAYERS
        )

        self.fc = nn.Linear(D_MODEL, vocab_size)

    def forward(self, x):
        B,T = x.shape
        pos_ids = torch.arange(T, device=x.device).unsqueeze(0).expand(B,T)
        x = self.embed(x) + self.pos(pos_ids)

        mask = torch.triu(torch.ones(T,T, device=x.device), diagonal=1).bool()
        x = self.transformer(x, mask=mask)

        return self.fc(x)

# Load trained model
model = ProteinLLM().to(device)
model.load_state_dict(torch.load("protein_llm_512_trainonly.pt", map_location=device))
model.eval()

# Sequence generation function
def generate_sequence(start="M", max_length=300, temperature=1.0):

    generated = [stoi[start]]

    for _ in range(max_length-1):

        x = torch.tensor(generated, dtype=torch.long).unsqueeze(0).to(device)

        with torch.no_grad():
            logits = model(x)

        next_logits = logits[0,-1] / temperature
        probs = torch.softmax(next_logits, dim=-1)

        next_token = torch.multinomial(probs,1).item()

        if next_token == PAD_IDX:
            break

        generated.append(next_token)

    return "".join([itos[i] for i in generated])


generated_sequences = []

for i in range(10):
    seq = generate_sequence(start="M", max_length=300)
    generated_sequences.append(seq)
    print(f"\nGenerated Protein {i+1}:\n{seq}")


with open("generated_enzymes.fasta","w") as f:
    for i,seq in enumerate(generated_sequences):
        f.write(f">generated_enzyme_{i+1}\n")
        f.write(seq+"\n")

print("\nGenerated sequences saved to generated_enzymes.fasta")