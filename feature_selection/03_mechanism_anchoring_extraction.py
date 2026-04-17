# feature_selection/03_mechanism_anchoring_extraction.py
import os
import pandas as pd
import numpy as np
from tqdm import tqdm
import warnings, time

warnings.filterwarnings("ignore")

# ==============================================================================
# GT-Mamba Feature Selection - Stage 3: Mechanism-Anchoring (Biological Priors)
# Purpose: Extract CpG sites within the Promoter Region (TSS -2000bp to +500bp)
#          of aging-related proteins (e.g., ProtAge20).
# Rationale: Anchoring the model with biological priors ensures robustness
#            and interpretability across different aging cohorts.
# ==============================================================================

# 1. Automatically locate the project root directory
CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DIR = os.path.dirname(CURRENT_DIR)

# Configure adaptive paths
DATA_DIR = os.path.join(BASE_DIR, "data")
RESULT_DIR = os.path.join(BASE_DIR, "results", "feature_selection")
os.makedirs(RESULT_DIR, exist_ok=True)

# Dependency file paths (Ensure these files are placed in the data/ directory)
ANNO_FILE = os.path.join(DATA_DIR, "450k_annotation.csv")
REFSEQ_FILE = os.path.join(DATA_DIR, "ncbiRefSeq_hg19.csv")
GENES_204_FILE = os.path.join(DATA_DIR, "204aps.csv")
GENES_20_FILE = os.path.join(DATA_DIR, "20aps.csv")

# Input matrix (Derived from Stage 2)
MERGED_FILE = os.path.join(RESULT_DIR, "top5000_merged_matrix.csv")

# Core biological parameters: TSS window range (Promoter Region)
WINDOW = (-2000, 500)


def extract_window_cpg(genes, window_range, anno_df, ref_df, desc):
    """
    Extract matching CpG probes based on a target gene list and TSS window.
    """
    ref_subset = ref_df[ref_df['name2'].isin(genes)].copy()

    # Determine the exact TSS position based on the forward (+) or reverse (-) strand
    ref_subset["tss"] = np.where(ref_subset["strand"] == "+", ref_subset["txStart"], ref_subset["txEnd"])

    selected_cpg = set()
    start_time = time.time()

    for _, row in tqdm(ref_subset.iterrows(), total=ref_subset.shape[0], desc=f"🔍 Scanning {desc}"):
        chr_name, tss = row["chrom"], row["tss"]

        # Locate probes falling within the specified window range in the annotation file
        region = anno_df[(anno_df["chr"] == chr_name) &
                         (anno_df["pos"].between(tss + window_range[0], tss + window_range[1]))]
        selected_cpg.update(region["Name"].values)

    duration = time.time() - start_time
    print(f"✅ {desc} Extraction Complete: {len(selected_cpg)} CpGs found | Time: {duration:.2f}s")
    return list(selected_cpg)


def main():
    # --- Step 1: Load Biological Prior Information ---
    print("📖 Loading annotation and gene-protein mapping files...")
    anno_df = pd.read_csv(ANNO_FILE)
    ref_df = pd.read_csv(REFSEQ_FILE)
    genes_204 = pd.read_csv(GENES_204_FILE)["Gene_name"].unique()
    genes_20 = pd.read_csv(GENES_20_FILE)["Gene_name"].unique()

    # --- Step 2: Extract Anchoring Sites ---
    cpg_204 = extract_window_cpg(genes_204, WINDOW, anno_df, ref_df, "ProtAge-204")
    cpg_20 = extract_window_cpg(genes_20, WINDOW, anno_df, ref_df, "ProtAge-20")

    # Save the extracted list for mandatory retention in subsequent Tier-based screening
    pd.DataFrame({"CpG": cpg_204}).to_csv(os.path.join(RESULT_DIR, "cpg_204_prior_list.csv"), index=False)

    # --- Step 3: Matrix Slicing and Alignment ---
    if os.path.exists(MERGED_FILE):
        print(f"\n📥 Loading candidate matrix: {os.path.basename(MERGED_FILE)}")
        df = pd.read_csv(MERGED_FILE)

        # Automatically identify and standardize foundational metadata columns
        cols_lower = {c.lower(): c for c in df.columns}
        base_cols = [cols_lower[w] for w in ["id", "age", "dataset", "split"] if w in cols_lower]

        # Generate the mechanism-anchored matrix specific to the 204 gene set
        valid_cpgs = [c for c in cpg_204 if c in df.columns]
        subset_df = df[base_cols + valid_cpgs]

        save_path = os.path.join(RESULT_DIR, "mechanism_anchored_matrix_204.csv")
        subset_df.to_csv(save_path, index=False)
        print(f"🎉 Final Anchored Matrix Saved: {save_path} | Shape: {subset_df.shape}")
    else:
        print(f"⚠️ Merge file not found: {MERGED_FILE}. Only CpG lists generated.")


if __name__ == "__main__":
    main()