import os
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from tqdm import tqdm

data_dir = "/datapool/home/info_wang/wanghui/data"
output_dir = "/datapool/home/info_wang/wanghui/file"
os.makedirs(output_dir, exist_ok=True)

geo_list = [
    'GSE50660', 'GSE53128', 'GSE55763', 'GSE60132', 'GSE64495',
    'GSE65638', 'GSE72773', 'GSE72775', 'GSE73103', 'GSE87571'
]

all_beta, all_pheno = [], []

print("📥 开始加载各个数据集...")
for gse in tqdm(geo_list, desc="加载数据集"):
    beta_file = os.path.join(data_dir, f"{gse}_beta.csv")
    pheno_file = os.path.join(data_dir, f"{gse}_pheno_healthy_with_celltypes.csv")

    if not os.path.exists(beta_file) or not os.path.exists(pheno_file):
        print(f"⚠️ 文件缺失: {gse}, 跳过")
        continue

    beta = pd.read_csv(beta_file, index_col=0).T
    pheno = pd.read_csv(pheno_file, index_col=0)

    # 统一 Age 单位
    if 'Age_unit' in pheno.columns:
        unit = pheno['Age_unit'].iloc[0].lower()
        if unit.startswith("month"):
            pheno['Age'] = pheno['Age'] / 12
        elif unit.startswith("week"):
            pheno['Age'] = pheno['Age'] / 52
        elif unit.startswith("day"):
            pheno['Age'] = pheno['Age'] / 365

    pheno['Dataset'] = gse

    # 保留样本名列，改为 ID
    pheno['ID'] = pheno.index

    # 忽略大小写和空格匹配样本
    beta.index = beta.index.str.lower().str.strip()
    pheno.index = pheno.index.str.lower().str.strip()
    common_samples = beta.index.intersection(pheno.index)
    beta = beta.loc[common_samples]
    pheno = pheno.loc[common_samples, ['ID', 'Age', 'Dataset'] + [c for c in pheno.columns if c not in ['ID','Age','Dataset']]]

    # 清理列名：小写 + 去掉首尾空格 + 去掉中间多余空格
    beta.columns = beta.columns.str.lower().str.strip().str.replace(r"\s+", "", regex=True)

    all_beta.append(beta)
    all_pheno.append(pheno)

print("✅ 全部数据加载完成")

# 取所有数据集共有且无缺失的 CpG
print("🔍 筛选所有数据集共有且无缺失的 CpG 位点...")
common_cpgs = set(all_beta[0].columns)
for beta in tqdm(all_beta[1:], desc="筛选公共 CpG"):
    beta_cols = beta.columns
    common_cpgs = {c for c in common_cpgs if c in beta_cols and beta[c].notna().all()}

common_cpgs = sorted(list(common_cpgs))
print(f"✅ 共 {len(common_cpgs)} 个 CpG 位点可用")

# 合并数据
beta_all = pd.concat([b[common_cpgs] for b in all_beta], axis=0)
pheno_all = pd.concat(all_pheno, axis=0)
full_df = pd.concat([pheno_all.reset_index(drop=True), beta_all.reset_index(drop=True)], axis=1)

# 划分 7:2:1
train_val_ratio = 0.9
val_ratio = 2 / 9  # 验证集占训练+验证

train_val_df, test_df = train_test_split(full_df, test_size=0.1, random_state=42, shuffle=True)
train_df, val_df = train_test_split(train_val_df, test_size=val_ratio, random_state=42, shuffle=True)

# 保存 CSV
train_csv = os.path.join(output_dir, "train.csv")
val_csv = os.path.join(output_dir, "val.csv")
test_csv = os.path.join(output_dir, "test.csv")

train_df.to_csv(train_csv, index=False)
val_df.to_csv(val_csv, index=False)
test_df.to_csv(test_csv, index=False)

print(f"✅ 划分完成 | 训练集: {train_df.shape[0]} | 验证集: {val_df.shape[0]} | 测试集: {test_df.shape[0]}")
print(f"✅ CSV 文件保存到 {output_dir}")
