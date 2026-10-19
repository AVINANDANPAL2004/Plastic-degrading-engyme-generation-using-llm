# Generating Plastic-Degrading Enzyme Sequences with a Protein Language Model

A research and learning project exploring how **Generative AI** can model protein sequences and generate candidate enzyme-like sequences related to plastic degradation.

The project uses a custom **Transformer-based protein language model implemented in PyTorch**. It is trained on cleaned amino-acid sequences, then sampled to produce candidate sequences. Supporting scripts perform basic sequence-similarity and exact-match checks against the training set.

> **Important scientific limitation:** generated sequences are computational candidates only. Valid amino-acid characters, novelty checks, or low sequence similarity do **not** establish that a sequence folds correctly, is a functional PET-degrading enzyme, or is safe. Experimental validation and appropriate bioinformatics review are required before drawing biological conclusions.

## Contents

- [Project goals](#project-goals)
- [How it works](#how-it-works)
- [Model configuration](#model-configuration)
- [Repository structure](#repository-structure)
- [Requirements](#requirements)
- [Setup](#setup)
- [Running the pipeline](#running-the-pipeline)
- [Outputs](#outputs)
- [Limitations and next steps](#limitations-and-next-steps)

## Project goals

- Explore generative modelling of protein amino-acid sequences.
- Prepare and combine PETase- and cutinase-related FASTA datasets.
- Train an autoregressive Transformer to predict the next amino-acid token.
- Generate candidate protein sequences and save them in FASTA format.
- Perform simple similarity and exact-copy checks against the training sequences.

## How it works

1. **Prepare the data.** FASTA records are filtered for valid standard amino-acid characters and a configured sequence-length range, then deduplicated by sequence when merging.
2. **Tokenize sequences.** Each amino-acid character is mapped to an integer ID. A padding token is added for batching.
3. **Train the language model.** The Transformer receives a sequence prefix and learns to predict the next token using cross-entropy loss and AdamW.
4. **Generate candidates.** Starting with methionine (`M`), the model samples subsequent amino-acid tokens until it reaches the configured maximum length or emits the padding token.
5. **Inspect novelty heuristically.** Separate scripts compare generated sequences with training sequences and check for exact duplicates.

## Model configuration

The current training script (`data/trainingllm.py`) uses the following configuration:

| Setting | Value |
|---|---:|
| Framework | PyTorch |
| Architecture | Transformer encoder with a causal attention mask |
| Embedding dimension | 512 |
| Transformer layers | 6 |
| Attention heads | 8 |
| Feed-forward dimension | 2,048 |
| Batch size | 16 |
| Maximum training epochs | 14 |
| Initial learning rate | 0.0008 |
| Optimizer | AdamW |
| Loss | Cross-entropy, ignoring padding tokens |
| Positional embedding limit | 512 tokens |
| Device | CUDA when available; otherwise CPU |

The repository includes `train_losses_512.json`, containing the recorded training-loss values, and `protein_llm_512_trainonly.pt`, a saved model state dictionary.

## Repository structure

```text
.
├── a_guidelines_and_req.txt
├── check_data.py
├── filter_data.py
├── merged_clean_enzymes.fasta
├── generated_enzymes.fasta
├── protein_llm_512_trainonly.pt
├── train_losses_512.json
└── data/
    ├── petase_clean.fasta
    ├── cutinase.fasta.gz
    ├── merged_clean_enzymes.fasta
    ├── uniprot_datasets.gz
    ├── uniprotkb_polyester_hydrolase_2026_02_05.fasta
    ├── cleaningandmerging.py
    ├── encode_sequences.py
    ├── trainingllm.py
    ├── generate.py
    ├── inferenceonly.py
    ├── run(with_proper_training).py
    ├── similaritycheck.py
    └── genrationchecking.py
```

Some files exist both at the repository root and inside `data/`. The scripts currently contain a mixture of root-relative and script-relative paths, so the working directory matters. The commands below are intended to be run from the **repository root**; inspect and adjust paths if you reorganize the files.

## Requirements

- Python 3.10 or newer is a reasonable starting point.
- PyTorch
- Biopython
- NumPy
- Matplotlib and pandas are listed in the project notes, though not every pipeline script uses them.

Training is resource-intensive. A CUDA-capable NVIDIA GPU is recommended; CPU execution is possible but may be significantly slower. The saved model checkpoint is approximately 77 MB.

## Setup

Clone the repository and enter its directory:

```bash
git clone https://github.com/AVINANDANPAL2004/Plastic-degrading-engyme-generation-using-llm.git
cd Plastic-degrading-engyme-generation-using-llm
```

Create and activate a virtual environment.

**Windows PowerShell**
```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
```

**macOS / Linux**
```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
```

Install the core dependencies:

```bash
pip install torch biopython numpy matplotlib pandas
```

For GPU training, install a PyTorch build compatible with your installed NVIDIA driver and CUDA environment. Follow the official PyTorch installation selector rather than assuming every machine supports CUDA.

## Running the pipeline

### Option A — Generate using the saved model

The repository includes a saved checkpoint and generated FASTA output. In principle, the generation script is:

```bash
python data/generate.py
```

This script expects the cleaned training file at `data/merged_clean_enzymes.fasta` and the checkpoint at the repository root as written in the current code. Run it from the repository root and ensure those paths exist. The script is configured to generate 10 sequences, starting with `M`, with a maximum length of 300 residues.

### Option B — Run the orchestration script

The repository includes an orchestration script intended to run data preparation, training (if no checkpoint is found), generation, and the two sequence checks:

```bash
python "data/run(with_proper_training).py"
```

The inference-only orchestration script is:

```bash
python data/inferenceonly.py
```

**Path note:** the orchestration scripts currently build paths relative to their own directory while some child scripts use root-relative paths. Therefore, the one-command pipeline may fail unless file locations and script paths are made consistent. If a command fails, inspect the referenced path and run the corresponding script from the directory it expects. The checkpoint-existence check in the orchestration scripts should also be aligned with the actual checkpoint location before relying on automatic training-skipping behavior.

### Useful individual scripts

| Script | Purpose |
|---|---|
| `data/cleaningandmerging.py` | Filter and merge PETase/cutinase FASTA sequences |
| `data/encode_sequences.py` | Demonstrate amino-acid-to-token encoding |
| `data/trainingllm.py` | Train the Transformer model |
| `data/generate.py` | Sample candidate sequences from the checkpoint |
| `data/similaritycheck.py` | Report the maximum `SequenceMatcher` similarity per generated sequence |
| `data/genrationchecking.py` | Check whether generated sequences exactly match a training sequence |
| `check_data.py` | Inspect the compressed dataset |
| `filter_data.py` | Filter sequences by length |

Because several scripts use relative paths, run them from the expected directory and verify input/output locations before executing a full run.

## Outputs

- **`generated_enzymes.fasta`** — generated candidate sequences in FASTA format.
- **`protein_llm_512_trainonly.pt`** — saved PyTorch model parameters.
- **`train_losses_512.json`** — recorded per-epoch training loss.
- **`merged_clean_enzymes.fasta`** — merged, filtered sequence dataset.

FASTA files can be opened with a text editor or analysed with tools such as Biopython. The checkpoint is a PyTorch state dictionary, not a standalone application.

## Limitations and next steps

- The training objective learns sequence-level statistical patterns; it does not directly optimize PET degradation activity.
- The current similarity script uses Python's `difflib.SequenceMatcher`. It is a lightweight heuristic, not a standard alignment-based protein identity or homology analysis.
- An exact-copy check only detects sequences identical to a training record; a sequence reported as new may still be highly similar to known proteins.
- The training script creates an 80/20 split but trains only on the 80% subset; it does not evaluate held-out loss or report a validation metric.
- The implementation should be tested for consistent file paths, reproducible data splits, padding/masking behavior, and checkpoint-vocabulary compatibility.
- Promising future directions include robust train/validation/test evaluation, alignment-based similarity searches, protein structure prediction, and experimental assays conducted by qualified researchers.

## Author

**Avinandan Pal**

GitHub: [@AVINANDANPAL2004](https://github.com/AVINANDANPAL2004)

---

This repository is an educational and computational research project. Generated protein sequences should not be described as functional plastic-degrading enzymes without appropriate validation.
