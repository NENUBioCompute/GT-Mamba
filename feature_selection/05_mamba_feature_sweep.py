# feature_selection/05_mamba_feature_sweep.py
import os, gc, time
import numpy as np
import pandas as pd
import torch
from sklearn.metrics import r2_score, mean_absolute_error
from tqdm import tqdm

# ==============================================================================
# GT-Mamba Feature Selection - Stage 5: Feature Sweep (Step-wise Evaluation)
# Purpose: Iterate through different feature counts (500-3000) to identify
#          the optimal subset size for cross-platform robustness.
# ==============================================================================

# 1. Automatically locate the project root and define paths
CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DIR = os.path.dirname(CURRENT_DIR)
DATA_PATH = os.path.join(BASE_DIR, "results", "feature_selection", "top5000_merged_matrix.csv")
RANKING_PATH = os.path.join(BASE_DIR, "results", "feature_selection", "ig_5000_ranking.csv")
OUTPUT_CSV = os.path.join(BASE_DIR, "results", "feature_selection", "mamba_feature_sweep_results.csv")


def main():
    # ... Data Loading and Partitioning (Train/Val/Test) ...

    results = []
    # Sweep range: from 500 to 3000 features, with a step size of 100
    scan_range = range(500, 3001, 100)

    for feature_num in tqdm(scan_range, desc="Scanning Features"):
        # 1. Select the Top N features based on IG ranking
        sel_feats = top_features[:feature_num]

        # 2. Train the Mamba probing model
        # [Invoke your training logic, log metrics such as R2, MAE, etc.]

        # 3. Record the evaluation metrics for the current feature subset
        results.append({
            "features": feature_num,
            "R2": r2,
            "MAE": mae,
            "RMSE": rmse
        })

    # Save the comprehensive feature sweep report
    pd.DataFrame(results).to_csv(OUTPUT_CSV, index=False)
    print(f"🎉 Feature sweep complete. Report saved to: {OUTPUT_CSV}")


if __name__ == "__main__":
    main()