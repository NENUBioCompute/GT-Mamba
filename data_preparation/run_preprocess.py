import os
import pandas as pd
import subprocess
import logging

# ==========================================
# 1. 数据集配置 (您的14个数据集)
# ==========================================
# 开发集 (10个)
DEV_COHORTS = [
    "GSE40279", "GSE42861", "GSE87571", "GSE55763", "GSE61107",
    "GSE51032", "GSE73103", "GSE105018", "GSE111629", "GSE125105"
]
# 外部验证集 (4个)
EXT_COHORTS = [
    "GSE72777", "GSE61496", "GSE77445", "GSE110554"
]

# 路径设置 (请根据实际情况修改)
RAW_BETA_DIR = "../data/raw_beta"
RAW_PHENO_DIR = "../data/raw_pheno"
OUTPUT_DIR = "../data/processed_data"
R_SCRIPT_PATH = "./impute_methylation.R"

# 日志配置
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger('GT-Mamba-Prep')


def run_pipeline():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    all_cohorts = DEV_COHORTS + EXT_COHORTS
    logger.info(f"Starting preprocessing for {len(all_cohorts)} cohorts...")

    for gse_id in all_cohorts:
        # 1. 定义输入路径
        input_beta = os.path.join(RAW_BETA_DIR, f"{gse_id}_beta.csv.gz")
        input_pheno = os.path.join(RAW_PHENO_DIR, f"{gse_id}_pheno.csv")

        # 2. 定义输出路径 (已去除 _final)
        # 中间临时文件
        temp_beta = os.path.join(OUTPUT_DIR, f"{gse_id}_temp.csv.gz")

        # ===> 修改处：生成的文件名不带 _final <===
        final_beta = os.path.join(OUTPUT_DIR, f"{gse_id}_beta.csv.gz")
        final_pheno = os.path.join(OUTPUT_DIR, f"{gse_id}_pheno.csv")

        if not os.path.exists(input_beta):
            logger.warning(f"Dataset missing: {gse_id}, skipping.")
            continue

        # ----------------------------------------------------
        # Step 1: 调用 R 脚本 (methyLImp2 插补 + BMIQ 归一化)
        # ----------------------------------------------------
        logger.info(f"[{gse_id}] Running R script (Imputation + BMIQ)...")
        # 将结果输出到 temp_beta
        cmd = ["Rscript", R_SCRIPT_PATH, input_beta, temp_beta]
        try:
            subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError:
            logger.error(f"[{gse_id}] R script failed.")
            continue

        # ----------------------------------------------------
        # Step 2: Python 端进行表型过滤 (Phenotype Filtering)
        # ----------------------------------------------------
        logger.info(f"[{gse_id}] Filtering samples by Age...")

        try:
            # 读取 R 处理好的 Beta 矩阵
            beta_df = pd.read_csv(temp_beta, index_col=0)
            # 读取原始表型
            pheno_df = pd.read_csv(input_pheno)

            # 确保有 Age 列
            if 'Age' not in pheno_df.columns:
                logger.error(f"[{gse_id}] 'Age' column missing in phenotype.")
                continue

            # 找出有 Age 的样本
            valid_pheno = pheno_df.dropna(subset=['Age'])
            valid_sample_ids = valid_pheno['SampleID'].tolist()

            # 取交集 (Beta矩阵列名 vs 表型SampleID)
            common_samples = sorted(list(set(beta_df.columns) & set(valid_sample_ids)))

            if not common_samples:
                logger.error(f"[{gse_id}] No matching samples found.")
                continue

            # 切片保存
            beta_final = beta_df[common_samples]
            pheno_final = valid_pheno[valid_pheno['SampleID'].isin(common_samples)]

            # 保持顺序一致
            pheno_final = pheno_final.set_index('SampleID').loc[common_samples].reset_index()

            # 保存到最终文件名 (无 _final)
            beta_final.to_csv(final_beta, compression='gzip')
            pheno_final.to_csv(final_pheno, index=False)

            # 删除临时文件
            if os.path.exists(temp_beta):
                os.remove(temp_beta)

            logger.info(f"[{gse_id}] Done. Samples: {len(common_samples)}")

        except Exception as e:
            logger.error(f"[{gse_id}] Filtering failed: {e}")

    logger.info("Pipeline Completed.")


if __name__ == "__main__":
    run_pipeline()