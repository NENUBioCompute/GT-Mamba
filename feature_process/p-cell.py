# import os
# import pandas as pd
# import numpy as np
# from sklearn.linear_model import LinearRegression
# from sklearn.model_selection import train_test_split
# from tqdm import tqdm
# import warnings
# warnings.filterwarnings("ignore")
#
# data_dir = "/datapool/home/info_wang/wanghui/data"
# output_dir = "/datapool/home/info_wang/wanghui/file"
# os.makedirs(output_dir, exist_ok=True)
#
# geo_list = [
#     'GSE50660', 'GSE53128', 'GSE55763', 'GSE60132', 'GSE64495',
#     'GSE65638', 'GSE72773', 'GSE72775', 'GSE73103', 'GSE87571'
# ]
#
# all_beta, all_pheno = [], []
#
# print("📥 开始加载各个数据集...")
# for gse in tqdm(geo_list, desc="加载数据集"):
#     beta_file = os.path.join(data_dir, f"{gse}_beta.csv")
#     pheno_file = os.path.join(data_dir, f"{gse}_pheno_healthy_with_celltypes.csv")
#     if not os.path.exists(beta_file) or not os.path.exists(pheno_file):
#         print(f"⚠️ 文件缺失: {gse}, 跳过")
#         continue
#
#     beta = pd.read_csv(beta_file, index_col=0).T
#     pheno = pd.read_csv(pheno_file, index_col=0)
#
#     # 统一 Age 单位
#     if 'Age_unit' in pheno.columns:
#         unit = pheno['Age_unit'].iloc[0].lower()
#         if unit.startswith("month"):
#             pheno['Age'] = pheno['Age'] / 12
#         elif unit.startswith("week"):
#             pheno['Age'] = pheno['Age'] / 52
#         elif unit.startswith("day"):
#             pheno['Age'] = pheno['Age'] / 365
#
#     pheno['Dataset'] = gse
#     pheno['ID'] = pheno.index
#
#     # 对齐样本
#     beta.index = beta.index.str.lower().str.strip()
#     pheno.index = pheno.index.str.lower().str.strip()
#     common_samples = beta.index.intersection(pheno.index)
#     beta = beta.loc[common_samples]
#     pheno = pheno.loc[common_samples, ['ID', 'Age', 'Dataset'] + [c for c in pheno.columns if c not in ['ID','Age','Dataset']]]
#
#     beta.columns = beta.columns.str.lower().str.strip().str.replace(r"\s+", "", regex=True)
#
#     all_beta.append(beta)
#     all_pheno.append(pheno)
#
# print("✅ 全部数据加载完成")
#
# # 取所有数据集共有且无缺失的 CpG
# print("🔍 筛选所有数据集共有且无缺失的 CpG 位点...")
# common_cpgs = set(all_beta[0].columns)
# for beta in tqdm(all_beta[1:], desc="筛选公共 CpG"):
#     beta_cols = beta.columns
#     common_cpgs = {c for c in common_cpgs if c in beta_cols and beta[c].notna().all()}
# common_cpgs = sorted(list(common_cpgs))
# print(f"✅ 共 {len(common_cpgs)} 个 CpG 位点可用")
#
# # 合并数据
# beta_all = pd.concat([b[common_cpgs] for b in all_beta], axis=0)
# pheno_all = pd.concat(all_pheno, axis=0)
# full_df = pd.concat([pheno_all.reset_index(drop=True), beta_all.reset_index(drop=True)], axis=1)
#
# # 划分训练/验证/测试
# train_val_ratio = 0.9
# val_ratio = 2 / 9  # 验证集占训练+验证
# train_val_df, test_df = train_test_split(full_df, test_size=0.1, random_state=42, shuffle=True)
# train_df, val_df = train_test_split(train_val_df, test_size=val_ratio, random_state=42, shuffle=True)
#
# # 保存划分文件
# train_df.to_csv(os.path.join(output_dir, "train.csv"), index=False)
# val_df.to_csv(os.path.join(output_dir, "val.csv"), index=False)
# test_df.to_csv(os.path.join(output_dir, "test.csv"), index=False)
#
# print(f"✅ 划分完成 | 训练集: {train_df.shape[0]} | 验证集: {val_df.shape[0]} | 测试集: {test_df.shape[0]}")
#
# # ================== 偏相关计算 ==================
# TOP_K = 5000
# CELL_TYPES = [c for c in train_df.columns if c.startswith('cell_')]
# cpg_cols = [c for c in common_cpgs]
#
# def calculate_partial_corr(df, cpg_cols, cell_types):
#     results = []
#     X_cells = df[cell_types].values
#     y_age = df['Age'].values
#     for probe in tqdm(cpg_cols, desc="计算偏相关"):
#         X_cpg = df[probe].values
#         # 残差法
#         age_resid = y_age - LinearRegression().fit(X_cells, y_age).predict(X_cells)
#         cpg_resid = X_cpg - LinearRegression().fit(X_cells, X_cpg).predict(X_cells)
#         if np.std(age_resid) == 0 or np.std(cpg_resid) == 0:
#             corr = 0.0
#         else:
#             corr = np.corrcoef(age_resid, cpg_resid)[0,1]
#         results.append({'cpg': probe, 'partial_corr': corr})
#     return pd.DataFrame(results)
#
# # 只在训练集上计算偏相关
# print("🔹 在训练集上计算 CpG 偏相关...")
# partial_corr_df = calculate_partial_corr(train_df, cpg_cols, CELL_TYPES)
# partial_corr_df['abs_corr'] = partial_corr_df['partial_corr'].abs()
# partial_corr_df = partial_corr_df.sort_values('abs_corr', ascending=False).reset_index(drop=True)
# partial_corr_df['rank'] = partial_corr_df.index + 1
# partial_corr_df.to_csv(os.path.join(output_dir, "all_cpgs_partial_corr.csv"), index=False)
# print(f"✅ 全部 CpG 偏相关结果已保存")
#
# # 取 Top K
# top_cpgs = partial_corr_df.head(TOP_K)['cpg'].tolist()
# print(f"🔹 Top {TOP_K} CpG 已选择")
#
# # 构建最终矩阵（合并训练/验证/测试集）
# final_rows = []
# for split_name, df in zip(['train','val','test'], [train_df, val_df, test_df]):
#     subset = df[['ID','Age','Dataset'] + top_cpgs].copy()
#     subset['Split'] = split_name
#     final_rows.append(subset)
#
# final_matrix = pd.concat(final_rows, axis=0).reset_index(drop=True)
# final_matrix.to_csv(os.path.join(output_dir, f"top{TOP_K}_cpgs_matrix.csv"), index=False)
# print(f"✅ 最终矩阵已保存: top{TOP_K}_cpgs_matrix.csv | 形状: {final_matrix.shape}")
import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression
import os
import gc
from tqdm import tqdm
import warnings

# 忽略不必要的 Pandas 警告
warnings.filterwarnings("ignore")

# ================= 1. 配置路径与参数 =================
BASE_DIR = "/datapool/home/info_wang/wanghui"
DATA_DIR = os.path.join(BASE_DIR, "data")  # 原始 beta/pheno 文件目录
OUTPUT_DIR = os.path.join(BASE_DIR, "file")  # 结果输出目录
REF_FILE = os.path.join(OUTPUT_DIR, "final_merged_cpg140_matrix.csv")  # 参考划分文件

os.makedirs(OUTPUT_DIR, exist_ok=True)

# 数据集列表
GEO_LIST = [
    'GSE50660', 'GSE53128', 'GSE55763', 'GSE60132', 'GSE64495',
    'GSE65638', 'GSE72773', 'GSE72775', 'GSE73103', 'GSE87571'
]

# 【关键配置】用于回归校正的细胞类型
# 根据您的图片，保留以下 6 种，排除 'cell_Gran' 和 'cell_Unkno'
USE_CELL_TYPES = [
    'cell_CD4T',
    'cell_CD8T',
    'cell_Bcell',
    'cell_Mono',
    'cell_NK',
    'cell_Eosinophils'  # 注意：请确认您的 pheno 文件中列名是否完全一致 (如 cell_Eosin vs cell_Eosinophils)
]

# ================= 2. 加载参考划分 (Master Key) =================
# ================= 2. 加载参考划分 (Master Key) [修复版] =================
print(f"📥 [Step 1] 加载参考划分文件: {REF_FILE}")
if not os.path.exists(REF_FILE):
    raise FileNotFoundError(f"❌ 参考文件未找到: {REF_FILE}")

# 1. 读取文件（先不限制 usecols，以防列名不匹配）
ref_df = pd.read_csv(REF_FILE)
print(f"ℹ️  参考文件列名预览: {ref_df.columns.tolist()}")

# 2. 智能识别 ID 列
# 尝试常见的 ID 列名变体
possible_id_names = ['ID', 'id', 'Id', 'sampleID', 'SampleID', 'geo_accession']
id_col_name = next((col for col in ref_df.columns if col in possible_id_names), None)

if id_col_name is None:
    # 如果找不到，默认第一列是 ID
    id_col_name = ref_df.columns[0]
    print(f"⚠️  未找到标准 ID 列名，默认使用第一列: {id_col_name}")

# 3. 智能识别 Split 列
possible_split_names = ['Split', 'split', 'SPLIT', 'Group', 'group']
split_col_name = next((col for col in ref_df.columns if col in possible_split_names), None)

if split_col_name is None:
    raise KeyError(f"❌ 在参考文件中找不到 Split 列，请检查列名是否为: {possible_split_names}")

# 4. 统一重命名为标准格式
ref_df = ref_df.rename(columns={id_col_name: 'ID', split_col_name: 'Split'})

# 5. 仅保留需要的列
ref_df = ref_df[['ID', 'Split']]

# 6. 确保 ID 格式统一 (小写去空格)
ref_df['ID'] = ref_df['ID'].astype(str).str.lower().str.strip()

# 构建映射字典: ID -> Split
id_map = dict(zip(ref_df['ID'], ref_df['Split']))

print(f"✅ 参考集加载完成，包含 {len(id_map)} 个样本。")
print(f"   (ID列来源: '{id_col_name}' -> 'ID', Split列来源: '{split_col_name}' -> 'Split')")
# ================= 3. 加载原始数据并构建 Global Train =================
print("\n📥 [Step 2] 开始加载原始数据并筛选训练集...")

train_fragments = []
valid_cpg_sets = []

for gse in tqdm(GEO_LIST, desc="Loading Datasets"):
    beta_path = os.path.join(DATA_DIR, f"{gse}_beta.csv")
    pheno_path = os.path.join(DATA_DIR, f"{gse}_pheno_healthy_with_celltypes.csv")

    if not os.path.exists(beta_path) or not os.path.exists(pheno_path):
        print(f"⚠️  文件缺失: {gse}，跳过")
        continue

    # --- 3.1 读取 Pheno ---
    pheno = pd.read_csv(pheno_path, index_col=0)
    pheno.index = pheno.index.str.lower().str.strip()

    # 筛选：只保留在 Reference 中存在的样本
    valid_indices = [idx for idx in pheno.index if idx in id_map]
    if not valid_indices:
        continue
    pheno = pheno.loc[valid_indices]

    # 标记 Split
    pheno['Split'] = pheno.index.map(id_map)

    # 【关键】只保留 Split == 'train' 的样本！
    # 这一步能大幅减少内存占用，并且防止验证/测试集数据泄露进入特征筛选环节
    pheno_train = pheno[pheno['Split'] == 'train'].copy()

    if pheno_train.empty:
        continue

    # --- 3.2 检查细胞列 ---
    # 检查当前数据集是否包含所有需要的细胞列
    missing_cells = [c for c in USE_CELL_TYPES if c not in pheno_train.columns]
    if missing_cells:
        print(f"⚠️  {gse} 缺失细胞列 {missing_cells}，跳过该数据集")
        continue

    # --- 3.3 统一 Age 单位 ---
    if 'Age_unit' in pheno_train.columns:
        unit = str(pheno_train['Age_unit'].iloc[0]).lower()
        if "month" in unit:
            pheno_train['Age'] = pheno_train['Age'] / 12
        elif "week" in unit:
            pheno_train['Age'] = pheno_train['Age'] / 52
        elif "day" in unit:
            pheno_train['Age'] = pheno_train['Age'] / 365

    # --- 3.4 读取 Beta ---
    # 读取全部 beta，转置 (行=样本)
    beta = pd.read_csv(beta_path, index_col=0).T
    beta.index = beta.index.str.lower().str.strip()

    # 只取训练集样本
    beta_train = beta.loc[pheno_train.index]

    # 简单清洗列名
    beta_train.columns = beta_train.columns.str.lower().str.strip()
    # 移除全NA列
    beta_train = beta_train.dropna(axis=1, how='all')

    # 记录该数据集有效的 CpG
    valid_cpg_sets.append(set(beta_train.columns))

    # --- 3.5 合并 Beta 和 Pheno (仅必要列) ---
    # 转为 float32 节省内存
    beta_train = beta_train.astype(np.float32)
    pheno_train['Age'] = pheno_train['Age'].astype(np.float32)
    for cell in USE_CELL_TYPES:
        pheno_train[cell] = pheno_train[cell].astype(np.float32)

    # 合并
    merged = pd.concat([pheno_train[['Age'] + USE_CELL_TYPES], beta_train], axis=1)
    train_fragments.append(merged)

    # 释放内存
    del beta, pheno, beta_train, pheno_train, merged
    gc.collect()

# ================= 4. 整合 Global Train Matrix =================
print("\n🔄 [Step 3] 整合全局训练集...")

if not valid_cpg_sets:
    raise ValueError("❌ 未加载到任何有效数据，请检查路径和参考文件 ID 匹配情况。")

# 取所有数据集共有的 CpG
common_cpgs = sorted(list(set.intersection(*valid_cpg_sets)))
print(f"✅ 筛选出 {len(common_cpgs)} 个公共 CpG 位点 (Intersection)。")

# 拼接所有碎片 (仅保留公共 CpG)
final_cols = ['Age'] + USE_CELL_TYPES + common_cpgs
global_train_df = pd.concat(
    [frag[final_cols] for frag in train_fragments],
    axis=0
).reset_index(drop=True)

# 再次清理内存
del train_fragments
gc.collect()

print(f"✅ 全局训练集构建完成 | 样本数: {global_train_df.shape[0]} | 特征数: {len(common_cpgs)}")
print(f"ℹ️  使用的细胞去混杂协变量: {USE_CELL_TYPES}")


# ================= 5. 极速偏相关计算 (核心算法) =================

def batch_partial_correlation(df, cpg_list, cell_list, batch_size=20000):
    """
    分块矩阵化计算残差偏相关:
    1. Regress Age ~ Cells -> get Age_Residuals
    2. Regress CpGs ~ Cells -> get CpG_Residuals (Batch by Batch)
    3. Corr(Age_Resid, CpG_Resid)
    """
    # 准备矩阵
    X_cells = df[cell_list].values  # (N, 6)
    y_age = df['Age'].values.reshape(-1, 1)  # (N, 1)

    total_cpgs = len(cpg_list)
    all_results = []

    print(f"\n⚡ [Step 4] 开始计算偏相关 (Batch Size = {batch_size})...")

    # --- 1. 计算 Age 残差 (只需计算一次) ---
    print("   ↳ 计算 Age 残差...")
    model_age = LinearRegression().fit(X_cells, y_age)
    age_resid = y_age - model_age.predict(X_cells)
    age_resid = age_resid - age_resid.mean()  # 中心化
    ss_age = np.sum(age_resid ** 2)  # Age 残差平方和 (标量)

    # --- 2. 循环计算 CpG 残差 (分块进行) ---
    for start_idx in tqdm(range(0, total_cpgs, batch_size), desc="   ↳ Batch Processing"):
        end_idx = min(start_idx + batch_size, total_cpgs)
        batch_cols = cpg_list[start_idx:end_idx]

        # 提取当前批次 CpG 矩阵
        Y_batch = df[batch_cols].values

        # 回归校正 CpG
        model_cpg = LinearRegression().fit(X_cells, Y_batch)
        cpg_resid = Y_batch - model_cpg.predict(X_cells)
        cpg_resid = cpg_resid - cpg_resid.mean(axis=0)  # 中心化

        # 计算 Pearson 相关系数
        # 分子: dot(age_res, cpg_res)
        numerator = np.dot(age_resid.T, cpg_resid).flatten()

        # 分母: sqrt(ss_age * ss_cpg)
        ss_cpg = np.sum(cpg_resid ** 2, axis=0)
        denominator = np.sqrt(ss_age * ss_cpg)

        # 安全除法 (避免除以0)
        with np.errstate(divide='ignore', invalid='ignore'):
            corrs = numerator / denominator
            corrs[denominator == 0] = 0.0

        # 收集结果
        batch_res = pd.DataFrame({
            'cpg': batch_cols,
            'partial_corr': corrs,
            'abs_corr': np.abs(corrs)
        })
        all_results.append(batch_res)

        # 显式清理临时变量
        del Y_batch, cpg_resid, numerator, denominator, model_cpg

    return pd.concat(all_results, ignore_index=True)


# 执行计算
result_df = batch_partial_correlation(global_train_df, common_cpgs, USE_CELL_TYPES)

# ================= 6. 保存结果 =================
print("\n📊 [Step 5] 正在排序并保存结果...")

# 按绝对值排序
result_df = result_df.sort_values('abs_corr', ascending=False).reset_index(drop=True)
result_df['rank'] = result_df.index + 1

# 保存完整相关性表
save_path_full = os.path.join(OUTPUT_DIR, "all_cpgs_partial_corr.csv")
result_df.to_csv(save_path_full, index=False)
print(f"✅ 完整偏相关结果已保存至: {save_path_full}")

# 保存 Top 5000 列表 (供后续图模型构建使用)
save_path_top = os.path.join(OUTPUT_DIR, "top_5000_cpg_list.csv")
top_5000 = result_df.head(5000)['cpg']
top_5000.to_csv(save_path_top, index=False, header=False)
print(f"✅ Top 5000 CpG 列表已保存至: {save_path_top}")

# 打印预览
print("\n🏆 Top 5 最显著位点:")
print(result_df.head())

print("\n🎉 全部流程执行完毕！")