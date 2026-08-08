import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader, random_split
from torch.nn.utils.rnn import pad_sequence
from Bio import SeqIO
import json

# ===========================
# CONFIGURATION
# ===========================
D_MODEL = 512
NUM_LAYERS = 6
N_HEAD = 8
FF_DIM = 2048
BATCH_SIZE = 16
MAX_EPOCHS = 14
LR = 8e-4

# ===========================
# DEVICE
# ===========================
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print("Using device:", device)

# ===========================
# LOAD DATA
# ===========================
seqs = [str(r.seq) for r in SeqIO.parse("data/merged_clean_enzymes.fasta", "fasta")]
print("Total sequences:", len(seqs))

# ===========================
# VOCAB + PAD
# ===========================
vocab = sorted(set("".join(seqs)))
vocab.append("<PAD>")

stoi = {c: i for i, c in enumerate(vocab)}
itos = {i: c for c, i in stoi.items()}
vocab_size = len(vocab)
PAD_IDX = stoi["<PAD>"]

# ===========================
# DATASET
# ===========================
class ProteinDataset(Dataset):
    def __init__(self, sequences):
        self.data = []
        for s in sequences:
            encoded = [stoi[c] for c in s]
            self.data.append(torch.tensor(encoded, dtype=torch.long))

    def __len__(self):
        return len(self.data)

    def __getitem__(self, idx):
        s = self.data[idx]
        return s[:-1], s[1:]

dataset = ProteinDataset(seqs)

# ===========================
# 80% SPLIT (TRAIN ONLY)
# ===========================
train_size = int(0.8 * len(dataset))  # 13240
test_size = len(dataset) - train_size

train_dataset, _ = random_split(dataset, [train_size, test_size])

print("Training sequences:", len(train_dataset))

# ===========================
# COLLATE FUNCTION
# ===========================
def collate_fn(batch):
    inputs = [item[0] for item in batch]
    targets = [item[1] for item in batch]

    inputs_padded = pad_sequence(inputs, batch_first=True, padding_value=PAD_IDX)
    targets_padded = pad_sequence(targets, batch_first=True, padding_value=PAD_IDX)

    return inputs_padded, targets_padded

train_loader = DataLoader(
    train_dataset,
    batch_size=BATCH_SIZE,
    shuffle=True,
    collate_fn=collate_fn,
    pin_memory=True
)

# ===========================
# MODEL
# ===========================
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
        B, T = x.shape
        pos_ids = torch.arange(T, device=x.device).unsqueeze(0).expand(B, T)

        x = self.embed(x) + self.pos(pos_ids)

        mask = torch.triu(
            torch.ones(T, T, device=x.device),
            diagonal=1
        ).bool()

        x = self.transformer(x, mask=mask)

        return self.fc(x)

model = ProteinLLM().to(device)

# ===========================
# OPTIMIZER
# ===========================
loss_fn = nn.CrossEntropyLoss(ignore_index=PAD_IDX)

optimizer = torch.optim.AdamW(
    model.parameters(),
    lr=LR,
    weight_decay=1e-4
)

scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(
    optimizer,
    T_max=MAX_EPOCHS
)

scaler = torch.amp.GradScaler("cuda")

# ===========================
# TRAINING LOOP (TRAIN ONLY)
# ===========================
loss_history = []

for epoch in range(MAX_EPOCHS):
    model.train()
    total_loss = 0

    for inp, tgt in train_loader:
        inp = inp.to(device, non_blocking=True)
        tgt = tgt.to(device, non_blocking=True)

        optimizer.zero_grad()

        with torch.amp.autocast("cuda"):
            output = model(inp)
            loss = loss_fn(
                output.reshape(-1, vocab_size),
                tgt.reshape(-1)
            )

        scaler.scale(loss).backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        scaler.step(optimizer)
        scaler.update()

        total_loss += loss.item()

    scheduler.step()

    avg_loss = total_loss / len(train_loader)
    loss_history.append(avg_loss)

    print(f"Epoch {epoch+1} | Train Loss: {avg_loss:.4f}")

# ===========================
# SAVE MODEL
# ===========================
torch.save(model.state_dict(), "protein_llm_512_trainonly.pt")

with open("train_losses_512.json", "w") as f:
    json.dump(loss_history, f)

print("Training complete.")