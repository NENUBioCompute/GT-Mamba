import torch
import numpy as np
import random
import os


def set_seed(seed=None):
    """
    Freeze all random seeds to ensure reproducible results.
    If seed is None, it will try to use CONFIG['RANDOM_SEED'].
    """
    if seed is None:
        # Avoid circular import: import inside function
        from configs.config import CONFIG
        seed = CONFIG.get('RANDOM_SEED', 42)

    random.seed(seed)
    os.environ['PYTHONHASHSEED'] = str(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)  # if you are using multi-GPU

    # Critical for CuDNN determinism
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False

    print(f"✅ Reproducibility: Global seed set to {seed}")