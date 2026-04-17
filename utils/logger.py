# utils/logger.py
import os
import logging


def setup_logging(save_dir, log_name="GT_Mamba_Benchmark.log"):
    os.makedirs(save_dir, exist_ok=True)
    log_path = os.path.join(save_dir, log_name)

    # 防止重复添加 handler
    logger = logging.getLogger("GT_Mamba")
    if not logger.handlers:
        logger.setLevel(logging.INFO)
        formatter = logging.Formatter('%(asctime)s [%(levelname)s] %(message)s')

        # 文件输出
        fh = logging.FileHandler(log_path, mode='w', encoding='utf-8')
        fh.setFormatter(formatter)
        logger.addHandler(fh)

        # 控制台输出
        ch = logging.StreamHandler()
        ch.setFormatter(formatter)
        logger.addHandler(ch)

    return logger