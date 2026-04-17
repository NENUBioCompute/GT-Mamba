# preprocessing/raw_data_pipeline/run_pipeline.py
import os
import pandas as pd
import subprocess
import logging

# ==============================================================================
# GT-Mamba: Official Data Preprocessing Pipeline
# Purpose: Execute end-to-end preprocessing including Imputation (R-side),
#          Normalization, and Phenotype-Genotype alignment.
# ==============================================================================

# 1. Automatically locate the project root directory
CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
# Navigate up twice to reach the GT-Mamba root from preprocessing/raw_data_pipeline/
BASE_DIR = os.path.dirname(os.path.dirname(CURRENT_DIR))

# Path configuration (pointing to the 'data' folder in the root)
RAW_BETA_DIR = os.path.join(BASE_DIR, "data", "raw_beta")
RAW_PHENO_DIR = os.path.join(BASE_DIR, "data", "raw_pheno")
OUTPUT_DIR = os.path.join(BASE_DIR, "data", "processed_data")

# R script path (located in the same directory as this Python script)
R_SCRIPT_PATH = os.path.join(CURRENT_DIR, "impute_methylation.R")

# 2. Dataset Configuration (10 Discovery Cohorts + 5 External Validation Cohorts = 15 total)
DEV_COHORTS = [
    "GSE40279", "GSE42861", "GSE87571", "GSE55763", "GSE61107",
    "GSE51032", "GSE73103", "GSE105018", "GSE111629", "GSE125105"
]

EXT_COHORTS = [
    "GSE72777", "GSE61496", "GSE77445", "GSE110554",
    "GSE132203"  # Includes the newly added EPIC 850k dataset
]

# Logging configuration
logging.basicConfig(level=logging.INFO, format='%(asctime)s - [%(levelname)s] - %(message)s')
logger = logging.getLogger('GT-Mamba-Prep')

def run_pipeline():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    all_cohorts = DEV_COHORTS + EXT_COHORTS
    logger.info("="*60)
    logger.info(f"🚀 Starting GT-Mamba preprocessing for {len(all_cohorts)} cohorts...")
    logger.info("="*60)

    for gse_id in all_cohorts:
        # Define input paths
        input_beta = os.path.join(RAW_BETA_DIR, f"{gse_id}_beta.csv.gz")
        input_pheno = os.path.join(RAW_PHENO_DIR, f"{gse_id}_pheno.csv")

        # Define intermediate and final output paths
        temp_beta = os.path.join(OUTPUT_DIR, f"{gse_id}_temp.csv.gz")
        final_beta = os.path.join(OUTPUT_DIR, f"{gse_id}_beta.csv.gz")
        final_pheno = os.path.join(OUTPUT_DIR, f"{gse_id}_pheno.csv")

        if not os.path.exists(input_beta):
            logger.warning(f"⚠️ Dataset missing: {gse_id}, skipping.")
            continue

        # ----------------------------------------------------
        # Step 1: Invoke R script (methyLImp2 imputation + BMIQ normalization)
        # ----------------------------------------------------
        logger.info(f"[{gse_id}] Running R script (Imputation + BMIQ)...")
        cmd = ["Rscript", R_SCRIPT_PATH, input_beta, temp_beta]
        try:
            subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError:
            logger.error(f"❌ [{gse_id}] R script failed.")
            continue

        # ----------------------------------------------------
        # Step 2: Phenotype filtering and sample alignment (Python)
        # ----------------------------------------------------
        logger.info(f"[{gse_id}] Filtering samples and aligning phenotypes...")
        try:
            beta_df = pd.read_csv(temp_beta, index_col=0)
            pheno_df = pd.read_csv(input_pheno)

            if 'Age' not in pheno_df.columns:
                logger.error(f"❌ [{gse_id}] 'Age' column missing in phenotype file.")
                continue

            # Filter samples with valid age and align Sample IDs across Beta and Pheno matrices
            valid_pheno = pheno_df.dropna(subset=['Age'])
            valid_sample_ids = valid_pheno['SampleID'].astype(str).tolist()
            common_samples = sorted(list(set(beta_df.columns.astype(str)) & set(valid_sample_ids)))

            if not common_samples:
                logger.error(f"❌ [{gse_id}] No matching samples found between Beta and Phenotype matrices.")
                continue

            # Extract aligned subsets and save results
            beta_final = beta_df[common_samples]
            pheno_final = valid_pheno[valid_pheno['SampleID'].astype(str).isin(common_samples)]
            pheno_final = pheno_final.set_index('SampleID').loc[common_samples].reset_index()

            # Save processed data in compressed format to save storage
            beta_final.to_csv(final_beta, compression='gzip')
            pheno_final.to_csv(final_pheno, index=False)

            # Cleanup intermediate temporary files
            if os.path.exists(temp_beta):
                os.remove(temp_beta)

            logger.info(f"✅ [{gse_id}] Preprocessing complete. Samples aligned: {len(common_samples)}")

        except Exception as e:
            logger.error(f"💥 [{gse_id}] Post-processing filtering failed: {e}")

    logger.info("🎉 Multi-cohort Preprocessing Pipeline Completed Successfully.")

if __name__ == "__main__":
    run_pipeline()