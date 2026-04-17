# feature_selection/06_stepwise_feature_refinement.py
import os
import random
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, TensorDataset
from sklearn.metrics import mean_absolute_error, r2_score, mean_squared_error
import matplotlib.pyplot as plt
from mambapy.mamba import Mamba, MambaConfig
from tqdm import tqdm
import warnings

warnings.filterwarnings("ignore")

# ==============================================================================
# GT-Mamba Feature Selection - Stage 6: Stepwise Feature Refinement
# Purpose: Iteratively add candidate features to the baseline to find the
#          optimal feature count (e.g., 198) that minimizes MAE.
# Output: Final model weights (.pth) and the final optimized feature list.
# ==============================================================================

# 1. Automatically locate the project root and define paths
CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DIR = os.path.dirname(CURRENT_DIR)

# Configure adaptive paths
DATA_DIR = os.path.join(BASE_DIR, "results", "feature_selection")
FINAL_DATA_DIR = os.path.join(BASE_DIR, "data")
OUTPUT_DIR = os.path.join(BASE_DIR, "results", "final_model")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Input files
BASELINE_FILE = os.path.join(DATA_DIR, "mechanism_anchored_matrix_204.csv")  # Anchored baseline
UNIQUE_FILE = os.path.join(DATA_DIR, "top5000_merged_matrix.csv")  # Candidate pool
CORR_FILE = os.path.join(DATA_DIR, "ig_5000_ranking.csv")  # Ranking reference


def set_seed(seed=42):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True


set_seed(42)
device = "cuda" if torch.cuda.is_available() else "cpu"


class DynamicMambaRegressor(nn.Module):
    def __init__(self, input_dim, d_model=64):
        super().__init__()
        self.embedding = nn.Linear(input_dim, d_model)
        self.norm = nn.LayerNorm(d_model)
        self.dropout = nn.Dropout(0.2)
        config = MambaConfig(d_model=d_model, n_layers=2, d_state=8, d_conv=2, pscan=True)
        self.mamba = Mamba(config)
        self.head = nn.Sequential(
            nn.Linear(d_model, 32),
            nn.ReLU(),
            nn.Linear(32, 1)
        )

    def forward(self, x):
        x = self.embedding(x)
        x = self.norm(x)
        x = self.dropout(x)
        x = x.unsqueeze(1)
        x = self.mamba(x)
        x = x.mean(dim=1)
        return self.head(x).squeeze(-1)


def main():
    print(f"✅ Training on device: {device}")

    # --- Step 1: Prepare Data and Rank Features ---
    df_base = pd.read_csv(BASELINE_FILE, index_col=0)
    df_unique = pd.read_csv(UNIQUE_FILE, index_col=0)
    df_corr = pd.read_csv(CORR_FILE)

    common_idx = df_base.index.intersection(df_unique.index)
    df_base, df_unique = df_base.loc[common_idx], df_unique.loc[common_idx]

    y = df_base['Age'].values
    splits = df_base['Split'].values

    # Retrieve candidate sites based on importance ranking
    id_col = df_corr.columns[0]
    sorted_candidates = [c for c in df_corr[id_col] if c in df_unique.columns and c not in df_base.columns]
    base_feats = [c for c in df_base.columns if c.startswith('cg') or c.startswith('ch')]

    # --- Step 2: Stepwise Selection Loop ---
    results = []
    best_mae = float('inf')
    best_k = 0
    best_weights = None

    max_k = 200  # Based on the study design, scan the top 200 most promising candidates
    step = 10

    print(f"🚀 Starting optimization curve scanning (Step={step})...")

    for k in tqdm(range(0, max_k + 1, step)):
        X_curr = df_base[base_feats].values
        if k > 0:
            X_curr = np.hstack([X_curr, df_unique[sorted_candidates[:k]].values])

        X_t = torch.tensor(X_curr, dtype=torch.float32).to(device)
        y_t = torch.tensor(y, dtype=torch.float32).to(device)

        # Train and evaluate the dynamically sized model
        mae, r2, rmse, r_val, weights = train_evaluate(
            X_t[splits == 'train'], y_t[splits == 'train'],
            X_t[splits == 'val'], y_t[splits == 'val']
        )

        results.append({'Added_Count': k, 'MAE': mae, 'R2': r2})

        if mae < best_mae:
            best_mae, best_k, best_weights = mae, k, weights

    # --- Step 3: Save Final Outputs ---
    print(f"\n🏆 Optimization Finished! Best k: +{best_k} features | Best MAE: {best_mae:.4f}")

    # A. Save the optimal feature list (CRITICAL: This represents the final selected feature set)
    final_features = base_feats + sorted_candidates[:best_k]
    pd.DataFrame({'CpG': final_features}).to_csv(os.path.join(FINAL_DATA_DIR, "best_features_list.csv"), index=False)

    # B. Save the best model weights for downstream benchmarking
    torch.save(best_weights, os.path.join(BASE_DIR, "checkpoints", "best_model.pth"))

    # C. Save the optimization curve plot for visual validation
    res_df = pd.DataFrame(results)
    plt.figure(figsize=(10, 6))
    plt.plot(res_df['Added_Count'], res_df['MAE'], marker='o')
    plt.title('MAE Optimization Curve')
    plt.savefig(os.path.join(OUTPUT_DIR, "optimization_curve.png"))

    print(f"🎉 Final model and feature list are ready for benchmarking!")


def train_evaluate(X_train, y_train, X_val, y_val):
    # [Keep your train_evaluate logic unchanged]
    # ...
    return mae, r2, rmse, r_corr, best_model_state


if __name__ == "__main__":
    main()