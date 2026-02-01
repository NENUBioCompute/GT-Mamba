# import os
# import pandas as pd
# import numpy as np
# from tqdm import tqdm
#
# # ==================== 配置 ====================
# data_dir = '/datapool/home/info_wang/wanghui/data'
# selected_features_path = '/datapool/home/info_wang/wanghui/file/top_5000_cpg_list.csv'
# output_matrix_path = '/datapool/home/info_wang/wanghui/file/top5000_cpgs_matrix.csv'
#
# geo_list = [
#     'GSE50660', 'GSE53128', 'GSE55763', 'GSE60132', 'GSE64495',
#     'GSE65638', 'GSE72773', 'GSE72775', 'GSE73103', 'GSE87571'
# ]
#
# # ==================== 读取选定位点 ====================
# selected_features = pd.read_csv(selected_features_path, header=None).iloc[:, 0].tolist()
#
# # ==================== 合并矩阵 ====================
# all_dfs = []
#
# for gse in tqdm(geo_list, desc="合并 GEO 数据集"):
#     beta_file = os.path.join(data_dir, f"{gse}_beta.csv")
#     pheno_file = os.path.join(data_dir, f"{gse}_pheno.csv")
#
#     if not os.path.exists(beta_file) or not os.path.exists(pheno_file):
#         print(f"⚠️ {gse} 文件缺失，跳过")
#         continue
#
#     beta_df = pd.read_csv(beta_file, index_col=0).T  # 行样本，列 CpG
#     pheno_df = pd.read_csv(pheno_file, index_col=0)
#
#     # 对齐样本
#     common_samples = beta_df.index.intersection(pheno_df.index)
#     if len(common_samples) == 0:
#         print(f"⚠️ {gse} 没有共同样本，跳过")
#         continue
#
#     dropped_samples = set(beta_df.index) - set(common_samples)
#     if dropped_samples:
#         print(f"⚠️ {gse} 丢弃样本: {dropped_samples}")
#
#     beta_df = beta_df.loc[common_samples]
#     pheno_df = pheno_df.loc[common_samples]
#
#     # 统一年龄单位为年
#     unit = pheno_df["Age_unit"].iloc[0].lower()
#     age = pheno_df["Age"].astype(float).copy()
#     if unit == "month":
#         age /= 12
#     elif unit == "week":
#         age /= 52
#     elif unit == "day":
#         age /= 365
#     beta_df['Age'] = age.values
#
#     beta_df['Dataset'] = gse
#
#     # 只保留 selected_features 中存在的列
#     features_in_data = [c for c in selected_features if c in beta_df.columns]
#     beta_df = beta_df[features_in_data + ['Age', 'Dataset']]
#
#     # 删除含缺失值的行
#     beta_df.dropna(inplace=True)
#
#     all_dfs.append(beta_df)
#
# # 合并所有数据集
# final_matrix = pd.concat(all_dfs, axis=0)
# final_matrix.index.name = 'Sample'
#
# # 保存
# final_matrix.to_csv(output_matrix_path)
# print(f"✅ 最终矩阵已保存: {output_matrix_path}")
# print(f"矩阵形状: {final_matrix.shape}")
import os
import pandas as pd
import numpy as np
from tqdm import tqdm
import warnings

warnings.filterwarnings('ignore')

# ==================== 配置路径 ====================
data_dir = '/datapool/home/info_wang/wanghui/data'
selected_features_path = '/datapool/home/info_wang/wanghui/file/top_5000_cpg_list.csv'
output_matrix_path = '/datapool/home/info_wang/wanghui/file/top5000_cpgs_matrix.csv'

# 参考划分文件路径 (包含 ID 和 Split)
reference_split_file = '/datapool/home/info_wang/wanghui/file/final_merged_cpg140_matrix.csv'

geo_list = [
    'GSE50660', 'GSE53128', 'GSE55763', 'GSE60132', 'GSE64495',
    'GSE65638', 'GSE72773', 'GSE72775', 'GSE73103', 'GSE87571'
]

# ==================== 1. 加载 Split 参考信息 ====================
print(f"📥 正在加载参考划分信息: {reference_split_file}")
if not os.path.exists(reference_split_file):
    raise FileNotFoundError(f"❌ 参考文件未找到: {reference_split_file}")

ref_df = pd.read_csv(reference_split_file)

# 智能识别 ID 和 Split 列名
id_col = next((c for c in ref_df.columns if c.lower() in ['id', 'sampleid', 'geo_accession']), ref_df.columns[0])
split_col = next((c for c in ref_df.columns if c.lower() in ['split', 'group']), None)

if split_col is None:
    raise KeyError("❌ 参考文件中未找到 'Split' 列")

# 统一 ID 格式 (小写去空格)
ref_df[id_col] = ref_df[id_col].astype(str).str.lower().str.strip()

# 构建映射字典: ID -> Split
id_split_map = dict(zip(ref_df[id_col], ref_df[split_col]))
print(f"✅ 参考信息加载完成，包含 {len(id_split_map)} 个样本的划分数据")

# ==================== 2. 读取选定位点 ====================
if not os.path.exists(selected_features_path):
    raise FileNotFoundError(f"❌ 特征列表文件未找到: {selected_features_path}")

# 尝试读取特征列表 (兼容有头无头)
df_tmp = pd.read_csv(selected_features_path, header=None)
if isinstance(df_tmp.iloc[0, 0], str) and df_tmp.iloc[0, 0].startswith('cg'):
    selected_features = df_tmp.iloc[:, 0].tolist()
else:
    # 可能是之前保存的 Series 格式，第二列才是特征名，或者是第一列
    # 这里假设第一列就是 CpG 名称
    selected_features = df_tmp.iloc[:, 0].tolist()

print(f"✅ 待提取特征数: {len(selected_features)}")

# ==================== 3. 合并矩阵 ====================
all_dfs = []

for gse in tqdm(geo_list, desc="合并 GEO 数据集"):
    beta_file = os.path.join(data_dir, f"{gse}_beta.csv")
    pheno_file = os.path.join(data_dir, f"{gse}_pheno.csv")

    if not os.path.exists(beta_file) or not os.path.exists(pheno_file):
        print(f"⚠️ {gse} 文件缺失，跳过")
        continue

    # 读取数据
    beta_df = pd.read_csv(beta_file, index_col=0).T  # 行样本，列 CpG
    pheno_df = pd.read_csv(pheno_file, index_col=0)

    # 索引标准化 (确保 ID 可以匹配)
    beta_df.index = beta_df.index.str.lower().str.strip()
    pheno_df.index = pheno_df.index.str.lower().str.strip()

    # 设置 Index 名称，方便后续 reset_index
    beta_df.index.name = 'ID'

    # 对齐样本 (Beta ∩ Pheno ∩ Split_Map)
    common_samples = beta_df.index.intersection(pheno_df.index)
    valid_samples = [s for s in common_samples if s in id_split_map]

    if len(valid_samples) == 0:
        continue

    # 切片数据
    beta_df = beta_df.loc[valid_samples]
    pheno_df = pheno_df.loc[valid_samples]

    # 添加元数据列
    beta_df['Split'] = beta_df.index.map(id_split_map)
    beta_df['Dataset'] = gse

    # 处理年龄
    if "Age_unit" in pheno_df.columns:
        unit = str(pheno_df["Age_unit"].iloc[0]).lower()
        age = pheno_df["Age"].astype(float).copy()
        if "month" in unit:
            age /= 12
        elif "week" in unit:
            age /= 52
        elif "day" in unit:
            age /= 365
        beta_df['Age'] = age.values
    else:
        beta_df['Age'] = pheno_df['Age'].values

    # 特征筛选
    features_in_data = [c for c in selected_features if c in beta_df.columns]

    # 只保留需要的列
    final_cols = ['Age', 'Dataset', 'Split'] + features_in_data
    beta_df = beta_df[final_cols]

    # 删除含缺失值的行
    beta_df.dropna(inplace=True)

    all_dfs.append(beta_df)

# ==================== 4. 保存结果 ====================
if all_dfs:
    final_matrix = pd.concat(all_dfs, axis=0)

    # 【关键修改】将 Index (ID) 重置为普通列
    final_matrix.reset_index(inplace=True)
    # 现在的列名里应该有了 'ID' (因为之前设置了 index.name='ID')
    # 如果之前没有设置成功，这就叫 'index'，我们重命名一下保险
    if 'ID' not in final_matrix.columns:
        final_matrix.rename(columns={'index': 'ID'}, inplace=True)

    # 【关键修改】调整列顺序，确保 ID 在第一列
    # 期望顺序: ID, Age, Dataset, Split, [CpGs...]
    meta_cols = ['ID', 'Age', 'Dataset', 'Split']
    feature_cols = [c for c in final_matrix.columns if c not in meta_cols]

    # 重新排列
    final_matrix = final_matrix[meta_cols + feature_cols]

    print(f"✅ 正在保存...")
    # index=False, 因为 ID 已经是一列数据了
    final_matrix.to_csv(output_matrix_path, index=False)

    print(f"✅ 最终矩阵已保存: {output_matrix_path}")
    print(f"   矩阵形状: {final_matrix.shape}")
    print(f"   前5列: {final_matrix.columns[:5].tolist()}")
else:
    print("❌ 未能合并任何数据，请检查输入文件。")