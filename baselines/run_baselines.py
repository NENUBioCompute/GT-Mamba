# baselines/run_baselines.py
import pandas as pd
import numpy as np
import os
import pyaging as pya
import anndata as ad
import logging
import time
import torch
from sklearn.metrics import mean_absolute_error, r2_score
from scipy.stats import pearsonr
import warnings

warnings.filterwarnings('ignore')


# --- Logging Setup ---
def setup_logging(save_dir, log_name="SOTA_Benchmark_Final.log"):
    os.makedirs(save_dir, exist_ok=True)
    log_path = os.path.join(save_dir, log_name)
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s [%(levelname)s] %(message)s',
        handlers=[
            logging.FileHandler(log_path, mode='w', encoding='utf-8'),
            logging.StreamHandler()
        ]
    )
    return logging.getLogger(__name__)


if __name__ == "__main__":
    # 1. Automatically locate the project root directory (navigate up one level from 'baselines/')
    BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    # 2. Configure relative paths (assuming the user places downloaded raw data in data/raw/)
    DATA_DIR = os.path.join(BASE_DIR, "data", "raw")
    DATA_850K_DIR = os.path.join(DATA_DIR, "850K")
    SAVE_DIR = os.path.join(BASE_DIR, "results", "baselines")

    os.makedirs(SAVE_DIR, exist_ok=True)
    logger = setup_logging(SAVE_DIR)

    GEO_LIST = ["GSE40279", "GSE61496", "GSE72777", "GSE77445", "GSE132203"]

    SOTA_CLOCKS = [
        "horvath2013", "hannum",
        "dnamphenoage", "grimage",
        "pchorvath2013", "pchannum", "pcphenoage",
        "altumage",
        "grimage2", "dunedinpace"
    ]

    logger.info("=" * 60)
    logger.info("⚠️ NOTE: Running this benchmark requires full whole-genome methylation matrices.")
    if torch.cuda.is_available():
        logger.info(f"🔥 GPU: {torch.cuda.get_device_name(0)}")
    else:
        logger.warning("⚠️ Currently using CPU for computation.")
    logger.info("=" * 60)

    results = []

    for gse in GEO_LIST:
        logger.info(f"\n▶️ Attempting to process cohort: {gse}")

        try:
            start_time = time.time()
            # --- Data Loading & Validation ---
            if gse == "GSE132203":
                beta_path = os.path.join(DATA_850K_DIR, "GSE132203_Geo_Submission_GTPEpic.csv.gz")
                pheno_path = os.path.join(DATA_850K_DIR, "GSE132203_pheno.csv")
            else:
                beta_path = os.path.join(DATA_DIR, f"{gse}_beta.csv")
                pheno_path = os.path.join(DATA_DIR, f"{gse}_pheno.csv")

            # Gracefully skip if raw data is missing, preventing hard crashes
            if not os.path.exists(beta_path):
                logger.warning(f"⏩ Raw data {beta_path} not found. Skipping. Ensure full datasets are downloaded from GEO.")
                continue

            # --- Load Data ---
            if gse == "GSE132203":
                beta = pd.read_csv(beta_path, index_col=0, compression='gzip')
                real_samples = [c for c in beta.columns if "_" in str(c) and "PVal" not in str(c)]
                beta = beta[real_samples].T
            else:
                beta = pd.read_csv(beta_path, index_col=0).T

            pheno = pd.read_csv(pheno_path, index_col=0)

            # --- Alignment (Beta & Phenotype) ---
            common = beta.index.intersection(pheno.index)
            if len(common) == 0:
                logger.error(f"❌ Sample alignment failed for {gse}")
                continue

            beta_sub = beta.loc[common]
            y_true = pheno.loc[common, 'Age'].values

            # --- 🌟 Standardized Imputation Strategy (Reviewer 2 Alignment) 🌟 ---
            logger.info("🩹 Applying standardized cohort-level mean imputation and zero-padding...")
            # Step 1: Cohort-level mean imputation for columns with partial missingness
            beta_sub = beta_sub.fillna(beta_sub.mean())
            # Step 2: Zero-padding for columns/probes that are entirely NaN (absent from the platform)
            beta_sub = beta_sub.fillna(0.0)

            clean_cpgs = [str(c).strip().replace('"', '').replace("'", "") for c in beta_sub.columns]

            adata = ad.AnnData(
                X=beta_sub.values.astype(np.float32),
                obs=pd.DataFrame(index=beta_sub.index),
                var=pd.DataFrame(index=clean_cpgs)
            )

            logger.info(f"🧠 Initiating age prediction using {len(SOTA_CLOCKS)} Epigenetic Clocks...")
            pya.pred.predict_age(adata, SOTA_CLOCKS)

            # --- Save Sample-level Predictions ---
            adata.obs['True_Age'] = y_true
            detail_path = os.path.join(SAVE_DIR, f"{gse}_SOTA_predictions.csv")
            adata.obs.to_csv(detail_path)

            metrics_row = {'Dataset': gse, 'Samples': len(common)}

            # --- Clock-wise Evaluation ---
            for clock in SOTA_CLOCKS:
                if clock in adata.obs.columns:
                    y_pred = adata.obs[clock].values
                    mask = ~np.isnan(y_true) & ~np.isnan(y_pred)

                    if mask.sum() > 5:
                        # RAW PREDICTIONS
                        mae_raw = mean_absolute_error(y_true[mask], y_pred[mask])
                        r_val, _ = pearsonr(y_true[mask], y_pred[mask])
                        r2_val = r2_score(y_true[mask], y_pred[mask])

                        # LINEAR CALIBRATION (Slope & Intercept correction)
                        slope, intercept = np.polyfit(y_pred[mask], y_true[mask], 1)
                        y_pred_calib = slope * y_pred[mask] + intercept
                        mae_calib = mean_absolute_error(y_true[mask], y_pred_calib)

                        # Store metrics
                        metrics_row[f'{clock}_MAE_raw'] = mae_raw
                        metrics_row[f'{clock}_MAE_calib'] = mae_calib
                        metrics_row[f'{clock}_R'] = r_val
                        metrics_row[f'{clock}_R2'] = r2_val

                        logger.info(f"    ∟ {clock:<12} | RAW={mae_raw:.2f} | CALIB={mae_calib:.2f}")
                    else:
                        logger.warning(f"    ∟ {clock} has insufficient valid data (Too many NaNs)")
                else:
                    logger.warning(f"    ∟ {clock} model did not output results")

            results.append(metrics_row)
            logger.info(f"✅ {gse} evaluation complete (Time: {time.time() - start_time:.1f}s)")

        except Exception as e:
            logger.error(f"💥 {gse} crashed: {e}")

    # --- Save Final Summary Report ---
    if results:
        results_df = pd.DataFrame(results).round(4)
        save_path = os.path.join(SAVE_DIR, "SOTA_raw_vs_calibrated.csv")
        results_df.to_csv(save_path, index=False)
        logger.info("=" * 60)
        logger.info(f"🏁 Benchmarking complete! Summary saved to: {save_path}")
        logger.info("=" * 60)
    else:
        logger.warning("⚠️ No datasets were successfully evaluated. Please check if raw data exists.")
