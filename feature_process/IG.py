import os
import warnings
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, TensorDataset
from sklearn.preprocessing import StandardScaler
from captum.attr import IntegratedGradients
from mambapy.mamba import Mamba, MambaConfig
from tqdm import tqdm

# ===================== 0. 全局配置 =====================
warnings.filterwarnings("ignore")


# 固定随机种子
def seed_everything(seed=42):
    import random
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)


seed_everything(42)
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"✅ 使用设备: {device}")

# 路径配置
DATA_FILE = "/datapool/home/info_wang/wanghui/file/top5000_cpgs_matrix.csv"
OUTPUT_FILE = "/datapool/home/info_wang/wanghui/file/ig_5000_ranking.csv"

# ===================== 1. 数据加载与预处理 =====================
print("\n📥 [1/4] 正在加载数据...")
if not os.path.exists(DATA_FILE):
    raise FileNotFoundError(f"❌ 文件未找到: {DATA_FILE}")

# 读取数据
df = pd.read_csv(DATA_FILE)

# --- 列名清洗 (去空格、转小写、处理标点) ---
df.columns = df.columns.str.strip().str.lower().str.replace(r"\s+", "", regex=True)

# --- 智能修复关键列名 ---
rename_map = {}
for c in df.columns:
    if c == 'age': rename_map[c] = 'age'
    if c == 'split': rename_map[c] = 'split'
    if c == 'dataset': rename_map[c] = 'dataset'
df.rename(columns=rename_map, inplace=True)

if 'split' not in df.columns or 'age' not in df.columns:
    print(f"当前列名: {df.columns.tolist()[:10]}")
    raise KeyError("❌ 数据缺失 'split' 或 'age' 列")

# --- 提取特征列 ---
exclude_cols = ['age', 'dataset', 'split', 'id', 'sampleid']
feature_cols = [c for c in df.columns if c not in exclude_cols]
print(f"   特征数量: {len(feature_cols)}")

# --- 只使用训练集进行探测 ---
# 我们需要解释的是“模型在训练集上学到了什么规律”
train_df = df[df['split'] == 'train'].reset_index(drop=True)

# 准备 X 和 y
X_raw = train_df[feature_cols].values.astype(np.float32)
y_raw = train_df['age'].values.astype(np.float32).reshape(-1, 1)

# --- 【关键步骤】标准化 (Standardization) ---
# 1. 特征标准化
x_scaler = StandardScaler()
X_std = x_scaler.fit_transform(X_raw)

# 2. 标签标准化 (这是解决 Loss 不下降的关键!)
y_scaler = StandardScaler()
y_std = y_scaler.fit_transform(y_raw)

# 转为 Tensor
X_t = torch.tensor(X_std).to(device)
y_t = torch.tensor(y_std).to(device)

print(f"   数据形状: X={X_t.shape}, y={y_t.shape}")

# ===================== 2. 模型构建与维度对齐 =====================
D = 10  # Mamba 嵌入维度
orig_dim = X_t.shape[1]
pad_len = (D - (orig_dim % D)) % D

# --- 维度补齐 (Padding) ---
# 如果特征数不是 D 的倍数，补 0
if pad_len > 0:
    print(f"   执行维度补齐: 末尾填充 {pad_len} 个 0")
    padding = torch.zeros(X_t.shape[0], pad_len, device=device)
    X_train_padded = torch.cat([X_t, padding], dim=1)
else:
    X_train_padded = X_t

input_dim_padded = X_train_padded.shape[1]
L = input_dim_padded // D


# --- 模型定义 ---
class EnhancedMambaRegressor(nn.Module):
    def __init__(self, input_dim, L, D):
        super().__init__()
        self.linear = nn.Linear(input_dim, L * D)
        self.norm = nn.LayerNorm(L * D)
        self.dropout = nn.Dropout(0.2)
        # Mamba 核心参数
        self.mamba = Mamba(MambaConfig(d_model=D, n_layers=2, d_state=8, d_conv=2, pscan=True))
        self.pool = nn.AdaptiveAvgPool1d(1)
        # 输出层
        self.fc = nn.Sequential(
            nn.Linear(D, 64),
            nn.GELU(),
            nn.Linear(64, 1)
        )
        self.L, self.D = L, D

    def forward(self, x):
        x = self.linear(x)
        x = self.norm(x)
        x = self.dropout(x)
        x = x.view(x.size(0), self.L, self.D)
        x = self.mamba(x)
        x = x.transpose(1, 2)
        x = self.pool(x).squeeze(-1)
        return self.fc(x).squeeze(-1)


model = EnhancedMambaRegressor(input_dim_padded, L, D).to(device)

# ===================== 3. 快速训练 (Probing) =====================
print("\n🚀 [2/4] 正在训练探测模型...")

# 优化器
optimizer = torch.optim.AdamW(model.parameters(), lr=1e-3, weight_decay=1e-4)
loss_fn = nn.MSELoss()

# DataLoader
batch_size = 64
dataset = TensorDataset(X_train_padded, y_t)
loader = DataLoader(dataset, batch_size=batch_size, shuffle=True)

model.train()
EPOCHS = 50  # 训练 50 轮以确保收敛

for epoch in range(EPOCHS):
    total_loss = 0
    for bx, by in loader:
        optimizer.zero_grad()
        pred = model(bx)  # 预测的是标准化后的年龄
        loss = loss_fn(pred, by.squeeze(-1))
        loss.backward()
        optimizer.step()
        total_loss += loss.item()

    avg_mse = total_loss / len(loader)

    # --- 实时还原真实误差 (岁) ---
    # RMSE_real = RMSE_scaled * std_age
    current_rmse_scaled = np.sqrt(avg_mse)
    real_error_years = current_rmse_scaled * y_scaler.scale_[0]

    # 每 10 轮打印一次
    if (epoch + 1) % 10 == 0 or epoch == 0:
        print(f"   Epoch {epoch + 1}/{EPOCHS} | Scaled MSE: {avg_mse:.4f} | 真实误差: {real_error_years:.2f} 岁")

print("   ✅ 探测模型训练完成！")

# ===================== 4. 计算 IG (Integrated Gradients) =====================
print("\n🔍 [3/4] 正在计算 Integrated Gradients...")

model.eval()
ig = IntegratedGradients(model)

# 选取计算样本 (500个足够代表整体分布，全量计算太慢)
subset_size = min(500, X_train_padded.shape[0])
# 为了结果稳定，不随机打乱，直接取前 N 个
X_subset = X_train_padded[:subset_size]

attr_list = []
ig_batch_size = 16

# 禁用梯度计算以节省显存 (IG 内部会处理梯度)
model.zero_grad()

for i in tqdm(range(0, subset_size, ig_batch_size), desc="Calculating IG"):
    batch = X_subset[i:i + ig_batch_size]
    # 计算归因
    batch_attr, _ = ig.attribute(batch, n_steps=50, return_convergence_delta=True)
    attr_list.append(batch_attr.detach().cpu())

# 合并结果
attr_all = torch.cat(attr_list, dim=0).numpy()

# 计算平均绝对重要性
feature_importance_padded = np.abs(attr_all).mean(axis=0)

# ===================== 5. 保存结果 =====================
print("\n💾 [4/4] 保存特征排名...")

# --- 关键：去除 Padding ---
# 只保留前 len(feature_cols) 个权重，后面的是 padding 的 0
feature_importance = feature_importance_padded[:len(feature_cols)]

# 构建 DataFrame
ranking_df = pd.DataFrame({
    'feature': feature_cols,
    'importance': feature_importance
})

# 降序排列 (最重要的在前)
ranking_df = ranking_df.sort_values(by='importance', ascending=False)

# 保存
ranking_df.to_csv(OUTPUT_FILE, index=False)
print(f"✅ 排名已保存至: {OUTPUT_FILE}")
print(f"   Top 5 特征: {ranking_df['feature'].head(5).tolist()}")
print(f"   Top 5 得分: {ranking_df['importance'].head(5).tolist()}")

print("\n🎉 流程结束！现在可以运行 3000 特征扫描脚本了。")