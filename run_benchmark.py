import os
import torch
import pandas as pd
import numpy as np
import warnings
from sklearn.metrics import mean_absolute_error

# --- Modular Imports ---
from models.gt_mamba import GraphMambaRegressor
from utils.logger import setup_logging
from utils.metrics import calculate_metrics, print_metrics_report
from configs.config import CONFIG

warnings.filterwarnings('ignore')


def main():
    # 1. Environment Setup
    BASE_DIR = os.path.dirname(os.path.abspath(__file__))
    DATA_DIR = os.path.join(BASE_DIR, "data")
    FEAT_PATH = os.path.join(DATA_DIR, "feature_list_198.csv")
    MODEL_PATH = os.path.join(BASE_DIR, "checkpoints", "best_model.pth")
    SAVE_DIR = os.path.join(BASE_DIR, "results", "sota")
    os.makedirs(SAVE_DIR, exist_ok=True)

    logger = setup_logging(SAVE_DIR, log_name="benchmark_log.log")
    logger.info("🚀 GT-Mamba Full Benchmark Execution Started...")

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    logger.info(f"🖥️ Execution Device: {device}")

    # 2. Load Core Feature List (198-CpG list)
    try:
        df_feat = pd.read_csv(FEAT_PATH, header=None)
        if len(df_feat) == 199: df_feat = pd.read_csv(FEAT_PATH, header=0)
        target_cpgs = [str(c).strip().replace('"', '') for c in df_feat.iloc[:198, 0].tolist()]
        logger.info(f"✅ Successfully loaded {len(target_cpgs)} topological features.")
    except Exception as e:
        logger.error(f"❌ Failed to load feature list: {e}")
        return

    # 3. Initialize Model
    try:
        model = GraphMambaRegressor(num_nodes=CONFIG['NUM_NODES']).to(device)
        model.load_state_dict(torch.load(MODEL_PATH, map_location=device))
        model.eval()
        logger.info(f"✅ Weights loaded: {os.path.basename(MODEL_PATH)}")
    except Exception as e:
        logger.error(f"❌ Model initialization failed: {e}")
        return

    # 4. Benchmarking Independent GEO Cohorts
    GEO_LIST = ["GSE40279", "GSE61496", "GSE72777", "GSE77445", "GSE132203"]
    summary_results = []

    for gse in GEO_LIST:
        logger.info("-" * 45)
        logger.info(f"▶️ Processing Cohort: {gse}")

        beta_path = os.path.join(DATA_DIR, f"{gse}_beta.csv")
        pheno_path = os.path.join(DATA_DIR, f"{gse}_pheno.csv")

        # Specific handling for the 850K EPIC array compressed file
        if gse == "GSE132203":
            beta_path = os.path.join(DATA_DIR, "GSE132203_Geo_Submission_GTPEpic.csv.gz")

        if not os.path.exists(beta_path) or not os.path.exists(pheno_path):
            logger.warning(f"⚠️ {gse} data/pheno file not found. Skipping...")
            continue

        try:
            # Data Loading (with Chunk-based strategy for 850K arrays)
            if gse == "GSE132203" and beta_path.endswith('.gz'):
                logger.info(f"⏳ Loading compressed 850K array for {gse}...")
                found_chunks = []
                reader = pd.read_csv(beta_path, index_col=0, chunksize=100000, compression='gzip')
                for chunk in reader:
                    chunk.index = [str(i).strip().replace('"', '') for i in chunk.index]
                    matched = chunk.index.intersection(target_cpgs)
                    if len(matched) > 0: found_chunks.append(chunk.loc[matched])
                beta_sub = pd.concat(found_chunks).T
            else:
                beta_sub = pd.read_csv(beta_path, index_col=0).T

            pheno = pd.read_csv(pheno_path, index_col=0)
            common = beta_sub.index.intersection(pheno.index)

            if len(common) == 0:
                logger.warning(f"⚠️ No overlapping samples found for {gse}. Skipping...")
                continue

            # --- FEATURE ALIGNMENT & SAMPLE-WISE MEAN IMPUTATION ---
            logger.info(f"🔄 Applying MARE Alignment & Sample-wise Imputation for {len(common)} samples...")
            x_df = beta_sub.loc[common].reindex(columns=target_cpgs)
            x_df = x_df.apply(lambda row: row.fillna(x_df.mean(axis=1)[row.name]), axis=1).fillna(0.5).clip(0, 1)
            y_true = pheno.loc[common, 'Age'].values

            # --- BATCH INFERENCE TO PREVENT OOM ---
            logger.info("⚡ Running Graph-Mamba Inference in batches...")
            X_tensor = torch.tensor(x_df.values, dtype=torch.float32).to(device)  # Fixed dimension here
            batch_size = 16  # Anti-OOM safe batch size
            all_preds = []

            with torch.no_grad():
                for i in range(0, len(X_tensor), batch_size):
                    batch_x = X_tensor[i: i + batch_size]
                    batch_preds = model(batch_x).cpu().numpy().flatten()
                    all_preds.append(batch_preds)

            y_pred = np.concatenate(all_preds)

            # Metric Evaluation via utils.metrics
            metrics_raw = calculate_metrics(y_true, y_pred)

            # Linear Calibration (Post-hoc for robust evaluation)
            slope, intercept = np.polyfit(y_pred, y_true, 1)
            y_pred_calib = slope * y_pred + intercept
            metrics_calib = calculate_metrics(y_true, y_pred_calib)

            logger.info(f"🏆 {gse} Results - RAW MAE: {metrics_raw['MAE']:.4f} | CALIB MAE: {metrics_calib['MAE']:.4f}")

            summary_results.append({
                'Dataset': gse,
                'N': len(common),
                'MAE_Raw': metrics_raw['MAE'],
                'MAE_Calib': metrics_calib['MAE'],
                'Pearson_r': metrics_raw['Pearson_r']
            })

            # Save per-cohort predictions for detailed analysis
            pred_df = pd.DataFrame({
                'Sample_ID': common,
                'True_Age': y_true,
                'Pred_Age_Raw': y_pred,
                'Pred_Age_Calib': y_pred_calib
            })
            pred_df.to_csv(os.path.join(SAVE_DIR, f"{gse}_predictions.csv"), index=False)

        except Exception as e:
            logger.error(f"💥 Fatal error processing {gse}: {e}")

    # Save final benchmark summary
    if summary_results:
        summary_path = os.path.join(SAVE_DIR, "benchmark_summary.csv")
        pd.DataFrame(summary_results).to_csv(summary_path, index=False)
        logger.info(f"✅ Benchmark completely finished. Final report saved to: {summary_path}")


if __name__ == "__main__":
    main()