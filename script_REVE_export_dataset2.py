import scipy.io
import torch
from transformers import AutoModel
import os
from tqdm import tqdm
import glob
import mat73
import numpy as np

# ================== CONFIG ==================
DATA_FOLDER = "Dataset2_REVE_eeg"          # folder with your per-subject .mat files
OUTPUT_FOLDER = "Dataset2_REVE_features"  # where embeddings will be saved
BATCH_SIZE = 32               # safe number; increase if you have lots of RAM
# ===========================================

os.makedirs(OUTPUT_FOLDER, exist_ok=True)

# Load the two models ONCE (this is the heavy part)
print("Loading REVE position bank and model... (this takes ~30-60 seconds the first time)")
pos_bank = AutoModel.from_pretrained("brain-bzh/reve-positions", trust_remote_code=True)
model = AutoModel.from_pretrained("brain-bzh/reve-base", trust_remote_code=True)
model.eval()  # important: inference mode

# Get list of all subject .mat files
mat_files = glob.glob(os.path.join(DATA_FOLDER, "*_reve_eeg.mat"))  # adjust pattern if needed
print(f"Found {len(mat_files)} subject files to process.")

for mat_file in tqdm(mat_files, desc="Processing subjects"):
    # Load the .mat
    data_dict = mat73.loadmat(mat_file, use_attrdict=True)
    struct = data_dict['eeg_reve_struct']

    data = torch.from_numpy(struct['data']).float()          # (n_epochs, n_chans, 1000)
    
    # === ROBUST electrode_names handling ===
    electrode_raw = struct['electrode_names']
    
    electrode_list = []
    for i in electrode_raw: electrode_list.append(i[0])
    
    electrode_names = [str(ch).strip() for ch in electrode_list if str(ch).strip()]

    # print(f"FINAL Electrode names ({len(electrode_names)} channels): {electrode_names[:10]} ...")
    
    subj_id = struct['subject_id']
    label = struct['label']

    # Get positions (only needs to be done once per subject)
    positions = pos_bank(electrode_names)                    # (n_chans, 3)
    positions = positions.expand(data.size(0), -1, -1)       # (n_epochs, n_chans, 3)

    # Process in small batches to save memory
    all_embeddings = []
    for i in range(0, data.shape[0], BATCH_SIZE):
        batch_data = data[i:i+BATCH_SIZE]
        batch_pos = positions[i:i+BATCH_SIZE]

        with torch.no_grad():
            output = model(batch_data, batch_pos)
            # REVE returns last_hidden_state by default → we mean-pool it
            emb = output.mean(dim=1)   # mean-pool over the patch / temporal dimension
            all_embeddings.append(emb)

    # Concatenate all epochs of this subject
    embeddings = torch.cat(all_embeddings, dim=0)        # (n_epochs_total, hidden_dim)

    # Save as new .mat for MATLAB
    out_dict = {
        'embeddings': np.mean(embeddings.numpy(),axis=1),
        'subject_id': subj_id,
        'label': label,
        'n_epochs': embeddings.shape[0]
    }
    out_file = os.path.join(OUTPUT_FOLDER, f"{str(subj_id).zfill(3)}_embeddings.mat")
    scipy.io.savemat(out_file, out_dict)

    print(f"✓ Saved embeddings for subject {subj_id} → {out_file}")

print("🎉 ALL DONE! All embeddings are in the 'embeddings' folder.")