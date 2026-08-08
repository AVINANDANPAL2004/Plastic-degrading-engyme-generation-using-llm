#1.make sure that you install the dependencies mentioned in a_guidelines_and_req.txt.
#2.for training pls check your gpu specification and uncomment line no.23.

import os
import subprocess
import sys

BASE = os.path.dirname(os.path.abspath(__file__))

print("🚀 Starting Protein LLM System...")


# ------------------------------
# Step 2: Check model
# ------------------------------
model_path = os.path.join(BASE, "data", "protein_llm_512_trainonly.pt")

if not os.path.exists(model_path):
    print(" ")

    subprocess.run([sys.executable, os.path.join(BASE, "cleaningandmerging.py")])
    subprocess.run([sys.executable, os.path.join(BASE, "encode_sequences.py")])
    subprocess.run([sys.executable, os.path.join(BASE, "trainingllm.py")]) 

else:
    print("✅ Model found. Skipping training.")

# ------------------------------
# Step 3: Generate
# ------------------------------
print("🧬 Generating proteins...")
subprocess.run([sys.executable, os.path.join(BASE, "generate.py")])
print('Similarity evaluation')
subprocess.run([sys.executable, os.path.join(BASE, "similaritycheck.py")])
print('Novelty check')
subprocess.run([sys.executable, os.path.join(BASE, "genrationchecking.py")])