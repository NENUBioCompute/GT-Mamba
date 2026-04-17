# feature_selection/01_partial_correlation_screening.py
import os
import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression
import gc
from tqdm import tqdm
import warnings

warnings.filterwarnings("ignore")

# ==============================================================================
# GT-Mamba Feature Selection - Stage 1: Partial Correlation Screening
# Purpose: Identify Age-related CpGs while controlling for cell-type heterogeneity.
# Method: Residual-based Partial Correlation (Regressing out 6 major cell types)
# ==============================================================================

# 1. Automatically locate the project root directory
CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DIR = os.path.dirname(CURRENT_DIR)

# Configure paths
DATA_DIR = os.path.join(BASE_DIR, "data", "processed_data")
OUTPUT_DIR = os.path.join(BASE_DIR, "results", "feature_selection")
REF_FILE = os.path.join(BASE_DIR, "data", "final_merged_cpg140_matrix.csv")

os.makedirs(OUTPUT_DIR, exist_ok=True)

# List of datasets (10 discovery cohorts)
GEO_LIST = [
    'GSE40279', 'GSE53128', 'GSE55763', 'GSE60132', 'GSE64495',
    'GSE65638', 'GSE72773', 'GSE72775', 'GSE73103', 'GSE87571'
]

# 6 major cell types used for regression adjustment (Deconvolution results)
USE_CELL_TYPES = ['cell_CD4T', 'cell_CD8T', 'cell_Bcell', 'cell_Mono', 'cell_NK', 'cell_Eosinophils']


def main():
    # --- Step 1: Load train/test split reference (Prevent data leakage) ---
    print(f"📥 Loading train/test split reference: {os.path.basename(REF_FILE)}")
    ref_df = pd.read_csv(REF_FILE)

    # Automatically identify ID and Split columns
    id_col = next((col for col in ref_df.columns if col.lower() in ['id', 'sampleid']), ref_df.columns[0])
    split_col = next((col for col in ref_df.columns if col.lower() in ['split', 'group']), None)

    ref_df = ref_df.rename(columns={id_col: 'ID', split_col: 'Split'})
    ref_df['ID'] = ref_df['ID'].astype(str).str.lower().str.strip()
    id_map = dict(zip(ref_df['ID'], ref_df['Split']))

    # --- Step 2: Integrate multi-cohort training data (Global Training Pool) ---
    train_fragments = []
    valid_cpg_sets = []

    print("📥 Integrating multi-cohort training data...")
    for gse in tqdm(GEO_LIST):
        beta_path = os.path.join(DATA_DIR, f"{gse}_beta.csv")
        pheno_path = os.path.join(DATA_DIR, f"{gse}_pheno_healthy_with_celltypes.csv")

        if not os.path.exists(beta_path) or not os.path.exists(pheno_path):
            continue

        # Load only samples marked as 'train'.
        # STRICTLY PROHIBIT test set data from participating in the screening process!
        pheno = pd.read_csv(pheno_path, index_col=0)
        pheno.index = pheno.index.str.lower().str.strip()
        pheno_train = pheno[pheno.index.map(id_map) == 'train'].copy()

        if pheno_train.empty: continue

        # Standardize Age unit to Years
        if 'Age_unit' in pheno_train.columns and "month" in str(pheno_train['Age_unit'].iloc[0]).lower():
            pheno_train['Age'] = pheno_train['Age'] / 12

        # Read Beta matrix (Rows=Samples, Columns=CpGs)
        beta = pd.read_csv(beta_path, index_col=0).T
        beta.index = beta.index.str.lower().str.strip()
        beta_train = beta.loc[pheno_train.index].astype(np.float32)

        # Record common CpG sites
        valid_cpg_sets.append(set(beta_train.columns))

        # Merge current cohort fragments
        merged = pd.concat([pheno_train[['Age'] + USE_CELL_TYPES], beta_train], axis=1)
        train_fragments.append(merged)
        del beta, pheno, beta_train;
        gc.collect()

    # Take the intersection of CpGs across all cohorts
    common_cpgs = sorted(list(set.intersection(*valid_cpg_sets)))
    global_train_df = pd.concat([frag[['Age'] + USE_CELL_TYPES + common_cpgs] for frag in train_fragments], axis=0)
    del train_fragments;
    gc.collect()

    print(f"✅ Global matrix ready: {global_train_df.shape[0]} samples, {len(common_cpgs)} CpGs.")

    # --- Step 3: Block-matrix partial correlation calculation (Core Algorithm) ---
    result_df = batch_partial_correlation(global_train_df, common_cpgs, USE_CELL_TYPES)

    # --- Step 4: Save Top 5000 candidate CpGs ---
    result_df = result_df.sort_values('abs_corr', ascending=False).reset_index(drop=True)
    save_path = os.path.join(OUTPUT_DIR, "top_5000_partial_corr_cpgs.csv")
    result_df.head(5000).to_csv(save_path, index=False)
    print(f"🎉 Stage 1 complete. Top 5000 candidates saved to: {save_path}")


def batch_partial_correlation(df, cpg_list, cell_list, batch_size=20000):
    """Efficiently compute the partial correlation matrix in batches."""
    X_cells = df[cell_list].values
    y_age = df['Age'].values.reshape(-1, 1)

    # Calculate Age residuals
    res_age = y_age - LinearRegression().fit(X_cells, y_age).predict(X_cells)
    res_age -= res_age.mean()
    ss_age = np.sum(res_age ** 2)

    all_corrs = []
    for i in tqdm(range(0, len(cpg_list), batch_size), desc="Calculating"):
        batch = cpg_list[i:i + batch_size]
        Y = df[batch].values

        # Calculate CpG residuals
        res_cpg = Y - LinearRegression().fit(X_cells, Y).predict(X_cells)
        res_cpg -= res_cpg.mean(axis=0)

        # Pearson correlation coefficient (Residual version)
        num = np.dot(res_age.T, res_cpg).flatten()
        den = np.sqrt(ss_age * np.sum(res_cpg ** 2, axis=0))
        all_corrs.append(pd.DataFrame({'cpg': batch, 'partial_corr': num / den, 'abs_corr': np.abs(num / den)}))

    return pd.concat(all_corrs)


if __name__ == "__main__":
    main()