# feature_selection/04_mamba_ig_ranking.py
import os
import warnings
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from sklearn.preprocessing import StandardScaler
from captum.attr import IntegratedGradients
from mambapy.mamba import Mamba, MambaConfig
from tqdm import tqdm

warnings.filterwarnings("ignore")

# ==============================================================================
# GT-Mamba Feature Selection - Stage 4: Feature Importance Ranking via IG
# Purpose: Rank 5000 candidate CpGs based on their contribution to age prediction
#          using Integrated Gradients (XAI).
# ==============================================================================

# 1. Automatically locate the project root directory
CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DIR = os.path.dirname(CURRENT_DIR)
DATA_FILE = os.path.join(BASE_DIR, "results", "feature_selection", "top5000_merged_matrix.csv")
OUTPUT_FILE = os.path.join(BASE_DIR, "results", "feature_selection", "ig_5000_ranking.csv")

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")


def main():
    print("📥 Loading data for IG attribution...")
    df = pd.read_csv(DATA_FILE)

    # Feature and label preparation
    exclude_cols = ['age', 'dataset', 'split', 'id', 'sampleid']
    feature_cols = [c for c in df.columns if c.lower() not in exclude_cols]

    # Strictly extract the training partition to prevent data leakage during feature attribution
    train_df = df[df['split'].str.lower() == 'train'].reset_index(drop=True)

    X_raw = train_df[feature_cols].values.astype(np.float32)
    y_raw = train_df['age'].values.astype(np.float32).reshape(-1, 1)

    # Standardization for Deep Learning convergence
    X_std = StandardScaler().fit_transform(X_raw)
    y_std = StandardScaler().fit_transform(y_raw)

    X_t = torch.tensor(X_std).to(device)
    y_t = torch.tensor(y_std).to(device)

    # --- Probing Model Construction ---
    D = 10
    orig_dim = X_t.shape[1]

    # Zero-padding to ensure the feature dimension is strictly divisible by D
    pad_len = (D - (orig_dim % D)) % D
    X_padded = torch.cat([X_t, torch.zeros(X_t.shape[0], pad_len, device=device)], dim=1) if pad_len > 0 else X_t

    L = X_padded.shape[1] // D
    model = EnhancedMambaRegressor(X_padded.shape[1], L, D).to(device)

    # Rapid probing training to establish baseline gradients
    train_probing_model(model, X_padded, y_t)

    # --- Calculate Integrated Gradients ---
    print("\n🔍 Calculating Integrated Gradients (Attribution)...")
    model.eval()
    ig = IntegratedGradients(model)

    # Extract a subset of 500 representative samples to compute the average attribution efficiently
    subset_size = min(500, X_padded.shape[0])
    X_subset = X_padded[:subset_size]

    attr, _ = ig.attribute(X_subset, n_steps=50, return_convergence_delta=True)

    # Calculate the mean absolute importance across the subset and discard padding dimensions
    importance = np.abs(attr.detach().cpu().numpy()).mean(axis=0)[:orig_dim]

    # Save the feature ranking results
    ranking_df = pd.DataFrame({'feature': feature_cols, 'importance': importance})
    ranking_df = ranking_df.sort_values(by='importance', ascending=False)
    ranking_df.to_csv(OUTPUT_FILE, index=False)
    print(f"✅ Feature ranking saved to: {OUTPUT_FILE}")

# [Your EnhancedMambaRegressor class and train_probing_model definitions remain here]