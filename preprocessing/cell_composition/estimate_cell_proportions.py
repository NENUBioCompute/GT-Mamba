# preprocessing/cell_composition/estimate_cell_proportions.py
import os
import pandas as pd
import numpy as np
import json
from pathlib import Path
from sklearn.linear_model import LinearRegression
import warnings

warnings.filterwarnings('ignore')

# ==============================================================================
# GT-Mamba Cell Type Deconvolution Module
# Purpose: Estimate proportions of 6 major leukocyte subtypes in blood samples.
# Method: Non-negative Least Squares (NNLS) Regression.
# Reference: FlowSorted.Blood.450k (Gold standard for 450K/EPIC arrays).
# ==============================================================================

# 1. Automatically locate the project root directory
CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DIR = os.path.dirname(os.path.dirname(CURRENT_DIR))

# Path Configuration (Adaptive relative paths)
# NOTE: Ensure the reference file is placed in data/references/
REFERENCE_FILE = os.path.join(BASE_DIR, "data", "references", "FlowSorted.Blood.450k.csv")
DATA_DIR = os.path.join(BASE_DIR, "data", "raw_data")
OUTPUT_DIR = os.path.join(BASE_DIR, "results", "cell_proportions")


def load_flowsorted_reference():
    """Loads and parses the blood cell reference dataset (Signature Matrix)."""
    if not os.path.exists(REFERENCE_FILE):
        print(f"❌ Reference file not found: {REFERENCE_FILE}")
        return None, None

    try:
        reference_data = pd.read_csv(REFERENCE_FILE, index_col=0)
        reference_pheno = pd.DataFrame(index=reference_data.columns)

        # Automatically infer cell types from column names
        # (Standard: CD8T, CD4T, NK, Bcell, Mono, Gran)
        cell_types = []
        for col in reference_data.columns:
            col_lower = col.lower()
            if 'cd8' in col_lower:
                cell_types.append('CD8T')
            elif 'cd4' in col_lower:
                cell_types.append('CD4T')
            elif 'nk' in col_lower:
                cell_types.append('NK')
            elif 'bcell' in col_lower or 'b_' in col_lower:
                cell_types.append('Bcell')
            elif 'mono' in col_lower:
                cell_types.append('Mono')
            elif 'gran' in col_lower or 'neut' in col_lower:
                cell_types.append('Gran')
            else:
                cell_types.append('Unknown')

        reference_pheno['CellType'] = cell_types
        return reference_data, reference_pheno
    except Exception as e:
        print(f"❌ Error loading reference: {e}")
        return None, None


def calculate_cell_proportions(beta_df, reference_beta, reference_pheno):
    """
    Perform deconvolution using Non-negative Least Squares (NNLS).
    Ensures that cell proportions are ≥ 0 and sum to 1.
    """
    common_probes = beta_df.index.intersection(reference_beta.index)
    if len(common_probes) < 100:
        print("⚠️ Warning: Too few common probes between data and reference.")
        return None

    beta_common = beta_df.loc[common_probes]
    reference_common = reference_beta.loc[common_probes]

    cell_types = reference_pheno['CellType'].unique()
    reference_means = {}
    valid_cell_types = []

    for ct in cell_types:
        samples = reference_pheno[reference_pheno['CellType'] == ct].index
        if len(samples) > 0:
            reference_means[ct] = reference_common[samples].mean(axis=1)
            valid_cell_types.append(ct)

    # Construct the Signature Matrix
    ref_matrix = pd.DataFrame(reference_means).fillna(0.5)

    # Initialize result container
    cell_proportions = pd.DataFrame(index=beta_common.columns, columns=valid_cell_types)

    # Execute NNLS sample by sample
    X = ref_matrix.values
    for sample in beta_common.columns:
        y = beta_common[sample].values
        mask = ~np.isnan(y)
        if mask.sum() < 100: continue

        # Core algorithm: Constrained non-negative linear regression
        model = LinearRegression(fit_intercept=False, positive=True)
        model.fit(X[mask], y[mask])

        # Normalization (Ensure sum of proportions equals 1)
        props = model.coef_
        if props.sum() > 0:
            props = props / props.sum()
        cell_proportions.loc[sample] = props

    return cell_proportions


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    print("🚀 Starting Cell Type Deconvolution Pipeline...")

    ref_beta, ref_pheno = load_flowsorted_reference()
    if ref_beta is None: return

    # Define target datasets to process
    # Example: datasets = ["GSE40279", "GSE61496"]
    datasets = ["GSE40279"]

    for gse in datasets:
        beta_path = os.path.join(DATA_DIR, f"{gse}_beta.csv")
        if not os.path.exists(beta_path):
            print(f"  ⏩ Skipping {gse}: Data not found.")
            continue

        print(f"  🔍 Deconvolving {gse}...")
        beta_df = pd.read_csv(beta_path, index_col=0)
        results = calculate_cell_proportions(beta_df, ref_beta, ref_pheno)

        if results is not None:
            out_path = os.path.join(OUTPUT_DIR, f"{gse}_cell_proportions.csv")
            results.to_csv(out_path)
            print(f"  ✅ Saved: {out_path}")


if __name__ == "__main__":
    main()