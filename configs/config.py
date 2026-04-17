# --- 1. configs/config.py ---

CONFIG = {
    # Architecture
    'D_MODEL': 64,
    'D_STATE': 16,
    'D_CONV': 4,
    'NUM_NODES': 198,
    'NUM_HEADS': 4,
    'DROPOUT': 0.1,

    # Reproducibility
    'RANDOM_SEED': 42,

    # Training Hyperparameters
    'LR': 1e-4,
    'BATCH_SIZE': 8,
    'EPOCHS': 400,
    'PATIENCE': 30
}