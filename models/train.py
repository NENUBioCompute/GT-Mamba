import os
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset
import numpy as np
import pandas as pd
import warnings

# --- Modular Imports ---
from models.gt_mamba import GraphMambaRegressor
from utils.logger import setup_logging
from utils.reproducibility import set_seed
from utils.metrics import calculate_metrics, print_metrics_report
from configs.config import CONFIG

warnings.filterwarnings("ignore")


def main():
    # 1. Path and Module Initialization
    BASE_DIR = os.path.dirname(os.path.abspath(__file__))
    DATA_PATH = os.path.join(BASE_DIR, 'data', 'methylation_matrix_198.csv')
    SAVE_DIR = os.path.join(BASE_DIR, 'checkpoints')
    RESULT_DIR = os.path.join(BASE_DIR, 'results', 'train_outputs')

    os.makedirs(SAVE_DIR, exist_ok=True)
    os.makedirs(RESULT_DIR, exist_ok=True)

    # Apply global seed for deterministic results
    set_seed(CONFIG['RANDOM_SEED'])
    logger = setup_logging(RESULT_DIR, log_name="training_log.log")
    logger.info("🔥 GT-Mamba Training Pipeline Initiated.")

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    logger.info(f"🖥️ Execution Device: {device}")

    # 2. Data Preparation and Sample-wise Imputation
    try:
        df = pd.read_csv(DATA_PATH)
        # Separate meta columns from features
        meta_cols = ['age', 'split', 'dataset', 'id', 'Age', 'Split', 'Dataset', 'ID']
        feat_cols = [c for c in df.columns if c not in meta_cols]

        logger.info("🔄 Applying Sample-wise Mean Imputation...")
        df[feat_cols] = df[feat_cols].apply(lambda row: row.fillna(row.mean()), axis=1).fillna(0.5)

        # Dataset Partitioning (Based on 'Split' column in the CSV)
        split = df['Split'].values
        X_all = df[feat_cols].values
        y_all = df['Age'].values

        def create_loader(X_sub, y_sub, shuffle=False):
            # [FIXED] Removed .unsqueeze(-1) to prevent 4D tensor error.
            # The model's forward() handles the feature dimension internally.
            X_t = torch.tensor(X_sub, dtype=torch.float32)
            y_t = torch.tensor(y_sub, dtype=torch.float32)
            return DataLoader(TensorDataset(X_t, y_t), batch_size=CONFIG['BATCH_SIZE'], shuffle=shuffle)

        train_loader = create_loader(X_all[split == 'train'], y_all[split == 'train'], shuffle=True)
        val_loader = create_loader(X_all[split == 'val'], y_all[split == 'val'])
        test_loader = create_loader(X_all[split == 'test'], y_all[split == 'test'])

        logger.info(f"📊 Partition Sizes: Train={len(train_loader.dataset)}, Val={len(val_loader.dataset)}")
    except Exception as e:
        logger.error(f"❌ Data loading/preprocessing failed: {e}")
        return

    # 3. Model Architecture and Optimization Setup
    model = GraphMambaRegressor(num_nodes=CONFIG['NUM_NODES']).to(device)
    optimizer = optim.AdamW(model.parameters(), lr=CONFIG['LR'], weight_decay=1e-4)
    scheduler = optim.lr_scheduler.ReduceLROnPlateau(optimizer, 'min', patience=10, factor=0.5)
    criterion = nn.MSELoss()

    # 4. Main Training Loop
    best_val_loss = float('inf')
    logger.info(f"🚀 Training for {CONFIG['EPOCHS']} epochs (Patience={CONFIG['PATIENCE']})...")

    for epoch in range(CONFIG['EPOCHS']):
        model.train()
        epoch_loss = 0
        for X_batch, y_batch in train_loader:
            X_batch, y_batch = X_batch.to(device), y_batch.to(device)
            optimizer.zero_grad()
            outputs = model(X_batch)
            loss = criterion(outputs, y_batch)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            optimizer.step()
            epoch_loss += loss.item()

        # Validation Phase
        model.eval()
        val_loss = 0
        with torch.no_grad():
            for X_batch, y_batch in val_loader:
                X_batch, y_batch = X_batch.to(device), y_batch.to(device)
                val_loss += criterion(model(X_batch), y_batch).item()

        avg_val_loss = val_loss / len(val_loader)
        scheduler.step(avg_val_loss)

        if (epoch + 1) % 10 == 0:
            logger.info(
                f"Epoch [{epoch + 1:03d}/{CONFIG['EPOCHS']}] | Train MSE: {epoch_loss / len(train_loader):.4f} | Val MSE: {avg_val_loss:.4f}")

        # Persistence: Save only the best performing weights
        if avg_val_loss < best_val_loss:
            best_val_loss = avg_val_loss
            torch.save(model.state_dict(), os.path.join(SAVE_DIR, 'best_model.pth'))
            logger.info(f"🌟 Checkpoint: Best model saved (Validation MSE: {best_val_loss:.4f})")

    # 5. Final Evaluation on Independent Test Partition
    logger.info("🔍 Final Evaluation on Held-out Test Set...")
    model.load_state_dict(torch.load(os.path.join(SAVE_DIR, 'best_model.pth')))
    model.eval()

    y_true_list, y_pred_list = [], []
    with torch.no_grad():
        for X_batch, y_batch in test_loader:
            y_pred_list.extend(model(X_batch.to(device)).cpu().numpy())
            y_true_list.extend(y_batch.numpy())

    # Generate Performance Report
    final_metrics = calculate_metrics(y_true_list, y_pred_list)
    print_metrics_report(final_metrics, "Internal Test Set")
    logger.info(f"🏆 Final Performance Summary: MAE={final_metrics['MAE']:.4f}, R2={final_metrics['R2']:.4f}")
    logger.info(f"✅ Training session completed. Weights and logs are stored in {SAVE_DIR}.")


if __name__ == "__main__":
    main()