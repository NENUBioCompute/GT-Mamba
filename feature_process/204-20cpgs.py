import os
import pandas as pd
import numpy as np
from tqdm import tqdm
import warnings, time

warnings.filterwarnings("ignore")

# ===================== 配置 =====================
data_dir = "/datapool/home/info_wang/wanghui/file"
output_dir = data_dir
os.makedirs(output_dir, exist_ok=True)

# 注释与基因文件路径
anno_file = os.path.join(data_dir, "450k_annotation.csv")
refseq_file = os.path.join(data_dir, "ncbiRefSeq_hg19.csv")
genes_204_file = os.path.join(data_dir, "204aps.csv")
genes_20_file = os.path.join(data_dir, "20aps.csv")

# 已经生成的合并矩阵
merged_file = os.path.join(data_dir, "final_window_cpg_matrix.csv")

# 窗口范围
window = (-2000, 500)

# ===================== 加载注释和基因信息 =====================
print("📖 加载注释与基因信息...")
anno_df = pd.read_csv(anno_file)
ref_df = pd.read_csv(refseq_file)
genes_204 = pd.read_csv(genes_204_file)["Gene_name"].unique()
genes_20 = pd.read_csv(genes_20_file)["Gene_name"].unique()
print(f"✅ 注释加载完成：anno({len(anno_df)}), refseq({len(ref_df)}), 基因数 204={len(genes_204)} / 20={len(genes_20)})")


# ===================== 提取窗口 CpG 函数 =====================
def extract_window_cpg(genes, window, desc):
    ref_subset = ref_df[ref_df['name2'].isin(genes)].copy()
    ref_subset["tss"] = np.where(ref_subset["strand"] == "+", ref_subset["txStart"], ref_subset["txEnd"])
    selected_cpg = set()

    start_time = time.time()
    for _, row in tqdm(ref_subset.iterrows(), total=ref_subset.shape[0], desc=f"🔍 {desc}"):
        chr_name, tss = row["chrom"], row["tss"]
        region = anno_df[(anno_df["chr"] == chr_name) &
                         (anno_df["pos"].between(tss + window[0], tss + window[1]))]
        selected_cpg.update(region["Name"].values)
    duration = time.time() - start_time
    print(f"✅ {desc} 提取完成：{len(selected_cpg)} 个 CpG | 用时 {duration:.2f}s")
    return list(selected_cpg)


# ===================== 提取两个集合的 CpG =====================
cpg_204 = extract_window_cpg(genes_204, window, "204基因窗口CpG")
cpg_20 = extract_window_cpg(genes_20, window, "20基因窗口CpG")

# 保存列表
pd.DataFrame({"CpG": cpg_204}).to_csv(os.path.join(data_dir, "cpg_204_list.csv"), index=False)
pd.DataFrame({"CpG": cpg_20}).to_csv(os.path.join(data_dir, "cpg_20_list.csv"), index=False)
print(f"💾 已保存 CpG 列表：cpg_204_list.csv ({len(cpg_204)} 个)，cpg_20_list.csv ({len(cpg_20)} 个)")

# ===================== 加载合并矩阵 =====================
print("\n📥 加载 final_window_cpg_matrix.csv ...")
load_start = time.time()
df = pd.read_csv(merged_file)
print(f"✅ 加载完成: {df.shape} | 用时 {time.time() - load_start:.2f}s")

# ===================== 不区分大小写匹配列名 =====================
print("\n⚙️ 自动识别基础列（忽略大小写）...")
cols_lower = {c.lower(): c for c in df.columns}
wanted = ["id", "age", "dataset", "split"]

base_cols = []
for w in wanted:
    if w in cols_lower:
        base_cols.append(cols_lower[w])
print(f"✅ 识别到基础列: {base_cols}")

# ===================== 拆分矩阵 =====================
split_tasks = [
    ("204gene", cpg_204),
    ("20gene", cpg_20)
]

for name, cpg_list in split_tasks:
    print(f"\n🧩 处理 {name} 集合 ...")
    t0 = time.time()

    # ✅ 只保留存在的 CpG 列
    valid_cpgs = [c for c in cpg_list if c in df.columns]
    missing_cpgs = len(cpg_list) - len(valid_cpgs)

    subset_cols = base_cols + valid_cpgs
    subset_df = df[subset_cols]
    save_path = os.path.join(output_dir, f"final_window_cpg_matrix_{name}.csv")
    subset_df.to_csv(save_path, index=False)

    print(f"✅ 保存 {name} 完成: {subset_df.shape} | 缺失 {missing_cpgs} 个 CpG | 用时 {time.time() - t0:.2f}s | 文件: {save_path}")

print("\n🎉 全部完成！")
