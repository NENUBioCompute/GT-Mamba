import os
import torch
import pandas as pd
import warnings
import numpy as np

# --- Modular Imports ---
from models.gt_mamba import GraphMambaRegressor
from utils.metrics import calculate_metrics, print_metrics_report
from configs.config import CONFIG

warnings.filterwarnings('ignore')


def main():
    print("\n" + "=" * 65)
    print("🚀 GT-Mamba: Quick Inference Demo for Reviewers")
    print("=" * 65)

    # 1. Environment and Path Configuration
    BASE_DIR = os.path.dirname(os.path.abspath(__file__))
    DEMO_DATA_PATH = os.path.join(BASE_DIR, "data", "methylation_matrix_198.csv")
    FEAT_PATH = os.path.join(BASE_DIR, "data", "feature_list_198.csv")
    MODEL_PATH = os.path.join(BASE_DIR, "checkpoints", "best_model.pth")
    SAVE_DIR = os.path.join(BASE_DIR, "results", "sota")
    os.makedirs(SAVE_DIR, exist_ok=True)

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"[*] 🖥️  Execution Device: {device}")

    # 2. Load Core 198-CpG Feature List
    try:
        df_feat = pd.read_csv(FEAT_PATH, header=None)
        if len(df_feat) == 199:
            df_feat = pd.read_csv(FEAT_PATH, header=0)
        target_cpgs = [str(c).strip().replace('"', '') for c in df_feat.iloc[:198, 0].tolist()]
        print(f"[*] ✅ Successfully loaded {len(target_cpgs)} topological features.")
    except Exception as e:
        print(f"[!] Error reading feature list: {e}")
        return

    # 3. Initialize Model and Load Pre-trained Weights
    print("[*] 🧠 Loading pre-trained weights...")
    try:
        # NUM_NODES is pulled from configs/config.py (should be 198)
        model = GraphMambaRegressor(num_nodes=CONFIG['NUM_NODES']).to(device)
        model.load_state_dict(torch.load(MODEL_PATH, map_location=device))
        model.eval()
        print(f"[*] ✅ Weights loaded from: {os.path.basename(MODEL_PATH)}")
    except Exception as e:
        print(f"[!] Model initialization failed: {e}")
        return

    # 4. Load Demonstration Dataset
    if not os.path.exists(DEMO_DATA_PATH):
        print(f"[!] Demo dataset not found at {DEMO_DATA_PATH}")
        return

    print(f"[*] 📂 Processing cohort: {os.path.basename(DEMO_DATA_PATH)}")
    demo_df = pd.read_csv(DEMO_DATA_PATH, index_col=0)

    # Extract ground truth age if available in columns (case-insensitive check)
    age_col = next((c for c in demo_df.columns if c.lower() == 'age'), None)
    true_ages = demo_df[age_col].values if age_col else None

    # 5. Feature Alignment & Sample-wise Mean Imputation
    print("[*] 🔄 Aligning topological nodes and performing Sample-wise Mean Imputation...")
    x_df = demo_df.reindex(columns=target_cpgs)

    # Fill missing probes with the mean of available probes for EACH specific sample
    # This is the core MARE robust alignment logic
    x_df = x_df.apply(lambda row: row.fillna(x_df.mean(axis=1)[row.name]), axis=1) \
        .fillna(0.5) \
        .clip(0, 1)

    # 6. Model Inference (Graph-based State Space Reasoning)
    print("[*] ⚡ Executing GT-Mamba inference engine...")

    # [FIXED] Removed .unsqueeze(-1) because model.forward() handles it internally
    X_tensor = torch.tensor(x_df.values, dtype=torch.float32).to(device)

    with torch.no_grad():
        preds = model(X_tensor).cpu().numpy().flatten()

    # 7. Metrics Evaluation and Reporting
    sample_ids = demo_df.index.tolist()

    if true_ages is not None:
        # Use the unified metrics module
        metrics = calculate_metrics(true_ages, preds)
        print_metrics_report(metrics, cohort_name="Demo Matrix (198-CpG)")
    else:
        print("\n[*] 📊 Prediction Summary (Top 10 Samples):")
        print(f"{'Sample ID':<25} | {'Predicted Age':<15}")
        print("-" * 45)
        for i in range(min(10, len(preds))):
            print(f"{str(sample_ids[i]):<25} | {preds[i]:<15.2f}")

    # 8. Save Detailed Predictions
    out_path = os.path.join(SAVE_DIR, "demo_predictions.csv")
    out_df = pd.DataFrame({'Sample_ID': sample_ids, 'Predicted_Age': preds})
    if true_ages is not None:
        out_df.insert(1, 'True_Age', true_ages)
        out_df['AgeAccel'] = out_df['Predicted_Age'] - out_df['True_Age']

    out_df.to_csv(out_path, index=False)

    print(f"\n✅ Demo task completed. Results exported to: {out_path}")
    print("=" * 65 + "\n")


if __name__ == "__main__":
    main()