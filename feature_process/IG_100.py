import os, gc, time, random
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, TensorDataset
from torch.cuda.amp import autocast, GradScaler
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score
from mambapy.mamba import Mamba, MambaConfig
from tqdm import tqdm
import warnings

# ===============================
# 1. 全局配置
# ===============================
warnings.simplefilter(action='ignore', category=FutureWarning)

# 路径配置
BASE_DIR = "/datapool/home/info_wang/wanghui/file"
DATA_PATH = os.path.join(BASE_DIR, "top5000_cpgs_matrix.csv")
RANKING_PATH = os.path.join(BASE_DIR, "ig_5000_ranking.csv")
OUTPUT_DIR = BASE_DIR
os.makedirs(OUTPUT_DIR, exist_ok=True)


# 固定随机种子
def seed_everything(seed=42):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
    os.environ['PYTHONHASHSEED'] = str(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False


seed_everything(42)
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"✅ 使用设备: {device}")

# ===============================
# 2. 鲁棒的数据加载
# ===============================
print("📥 正在加载数据...")
if not os.path.exists(DATA_PATH):
    raise FileNotFoundError(f"❌ 数据文件未找到: {DATA_PATH}")

df = pd.read_csv(DATA_PATH)

# --- 列名清洗 ---
df.columns = df.columns.str.strip().str.lower().str.replace(r"\s+", "", regex=True)
rename_map = {}
for c in df.columns:
    if c == 'age': rename_map[c] = 'Age'
    if c == 'split': rename_map[c] = 'Split'
    if c == 'dataset': rename_map[c] = 'Dataset'
df.rename(columns=rename_map, inplace=True)

if 'Split' not in df.columns or 'Age' not in df.columns:
    raise KeyError(f"❌ 数据缺失关键列。当前列: {df.columns.tolist()[:10]}")

# --- 加载排名 ---
if not os.path.exists(RANKING_PATH):
    raise FileNotFoundError(f"❌ 排名文件未找到: {RANKING_PATH}")
rank_df = pd.read_csv(RANKING_PATH)
# 假设第一列是特征名
feature_col = rank_df.columns[0]
top_features_raw = rank_df[feature_col].tolist()

# --- 特征对齐 ---
# 必须清洗排名里的特征名，以防格式不一致
clean_features = [f.strip().lower() for f in top_features_raw]
# 只保留矩阵中存在的特征
matrix_cols = set(df.columns)
top_features = [f for f in clean_features if f in matrix_cols]

print(f"✅ 数据加载成功 | 样本数: {df.shape[0]}")
print(f"   排名特征数: {len(top_features_raw)} -> 有效特征数: {len(top_features)}")

# ===============================
# 3. 数据集划分
# ===============================
df_train = df[df['Split'] == 'train'].reset_index(drop=True)
df_val = df[df['Split'] == 'val'].reset_index(drop=True)
df_test = df[df['Split'] == 'test'].reset_index(drop=True)

print(f"📊 划分: Train={len(df_train)} | Val={len(df_val)} | Test={len(df_test)}")

# 释放内存
del df
gc.collect()


# ===============================
# 4. 模型定义
# ===============================
class EnhancedMambaRegressor(nn.Module):
    def __init__(self, input_dim, L, D):
        super().__init__()
        self.L, self.D = L, D
        self.linear = nn.Linear(input_dim, L * D)
        self.norm = nn.LayerNorm(L * D)
        self.dropout = nn.Dropout(0.2)
        # Mamba 块参数
        self.mamba = Mamba(MambaConfig(d_model=D, n_layers=2, d_state=8, d_conv=2, pscan=True))
        self.pool = nn.AdaptiveAvgPool1d(1)
        self.fc = nn.Sequential(
            nn.Linear(D, 64),
            nn.GELU(),
            nn.Dropout(0.2),
            nn.Linear(64, 1)
        )

    def forward(self, x):
        x = self.linear(x)
        x = self.norm(x)
        x = self.dropout(x)
        x = x.view(x.size(0), self.L, self.D)
        x = self.mamba(x)
        x = x.transpose(1, 2)
        x = self.pool(x).squeeze(-1)
        return self.fc(x).squeeze(-1)


# ===============================
# 5. 训练函数 (带目标值标准化)
# ===============================
def run_training(model, train_loader, val_loader, optimizer, num_epochs=100):
    loss_fn = nn.MSELoss()
    scaler = GradScaler()
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode='min', patience=5, factor=0.5)

    best_val_loss = float('inf')
    best_state = None
    patience = 0
    max_patience = 8  # 早停轮数

    for epoch in range(num_epochs):
        # --- 训练 ---
        model.train()
        for X, y in train_loader:
            X, y = X.to(device), y.to(device)
            optimizer.zero_grad()
            with autocast():
                pred = model(X)
                loss = loss_fn(pred, y)  # 注意：这里 y 是标准化后的
            scaler.scale(loss).backward()
            scaler.step(optimizer)
            scaler.update()

        # --- 验证 ---
        model.eval()
        val_loss = 0
        with torch.no_grad():
            for X, y in val_loader:
                X, y = X.to(device), y.to(device)
                pred = model(X)
                val_loss += loss_fn(pred, y).item()
        val_loss /= len(val_loader)

        scheduler.step(val_loss)

        if val_loss < best_val_loss:
            best_val_loss = val_loss
            best_state = model.state_dict()
            patience = 0
        else:
            patience += 1
            if patience >= max_patience:
                break

    if best_state:
        model.load_state_dict(best_state)
    return model


# ===============================
# 6. 特征扫描主循环
# ===============================
results = []
best_overall_r2 = -float('inf')
best_model_state = None
best_feature_num = None

# 扫描范围：500 到 3000，步长 100
scan_range = range(500, 3001, 100)
# 确保不超过实际有效特征数
max_valid = len(top_features)
scan_range = [n for n in scan_range if n <= max_valid]

print(f"\n🚀 开始特征扫描 (Range: {scan_range[0]} - {scan_range[-1]})...")

for feature_num in tqdm(scan_range, desc="Scanning"):
    t0 = time.time()

    # 1. 准备特征数据
    sel_feats = top_features[:feature_num]

    # 转换为 float32
    X_train_raw = df_train[sel_feats].values.astype(np.float32)
    y_train_raw = df_train['Age'].values.astype(np.float32).reshape(-1, 1)

    X_val_raw = df_val[sel_feats].values.astype(np.float32)
    y_val_raw = df_val['Age'].values.astype(np.float32).reshape(-1, 1)

    X_test_raw = df_test[sel_feats].values.astype(np.float32)
    y_test_raw = df_test['Age'].values.astype(np.float32).reshape(-1, 1)

    # 2. 【关键】标准化 (Standardization)
    # 特征标准化
    x_scaler = StandardScaler()
    X_train = x_scaler.fit_transform(X_train_raw)
    X_val = x_scaler.transform(X_val_raw)
    X_test = x_scaler.transform(X_test_raw)

    # 标签标准化 (解决 Loss 不下降的核心)
    y_scaler = StandardScaler()
    y_train = y_scaler.fit_transform(y_train_raw).ravel()  # ravel 转为 1D
    y_val = y_scaler.transform(y_val_raw).ravel()
    # 测试集标签保持原始值用于最终评估，不需要 transform，后面我们会反标准化预测值

    # 3. 维度调整
    D = 16  # 嵌入维度
    L = feature_num // D
    input_dim = L * D

    X_train = X_train[:, :input_dim]
    X_val = X_val[:, :input_dim]
    X_test = X_test[:, :input_dim]

    # 4. DataLoader
    batch_size = 16
    train_loader = DataLoader(TensorDataset(torch.tensor(X_train), torch.tensor(y_train)),
                              batch_size=batch_size, shuffle=True)
    val_loader = DataLoader(TensorDataset(torch.tensor(X_val), torch.tensor(y_val)),
                            batch_size=batch_size, shuffle=False)
    # 测试集 loader 不需要 y，因为我们要反归一化预测值后再对比
    test_loader_x = DataLoader(torch.tensor(X_test), batch_size=batch_size, shuffle=False)

    # 5. 训练模型
    model = EnhancedMambaRegressor(input_dim, L, D).to(device)
    optimizer = torch.optim.AdamW(model.parameters(), lr=1e-3, weight_decay=1e-4)

    model = run_training(model, train_loader, val_loader, optimizer, num_epochs=100)

    # 6. 测试集评估
    model.eval()
    preds_scaled = []
    with torch.no_grad():
        for batch_x in test_loader_x:
            batch_x = batch_x.to(device)
            p = model(batch_x).cpu().numpy()
            preds_scaled.extend(p)

    preds_scaled = np.array(preds_scaled).reshape(-1, 1)

    # 【关键】反标准化预测值 (还原为真实年龄)
    y_pred_real = y_scaler.inverse_transform(preds_scaled).ravel()
    y_true_real = y_test_raw.ravel()  # 原始真实年龄

    # 计算真实指标
    r2 = r2_score(y_true_real, y_pred_real)
    pearson_r = np.corrcoef(y_true_real, y_pred_real)[0, 1]
    mae = mean_absolute_error(y_true_real, y_pred_real)
    rmse = np.sqrt(mean_squared_error(y_true_real, y_pred_real))

    # 记录
    results.append({
        "features": feature_num,
        "R2": r2,
        "Pearson_R": pearson_r,
        "MAE": mae,
        "RMSE": rmse
    })

    tqdm.write(f"Feat: {feature_num} | R2: {r2:.3f} | MAE: {mae:.2f} | RMSE: {rmse:.2f}")

    # 保存最佳状态
    if r2 > best_overall_r2:
        best_overall_r2 = r2
        best_feature_num = feature_num
        best_model_state = model.state_dict().copy()

    # 清理显存
    del model, optimizer, X_train, X_val, X_test
    torch.cuda.empty_cache()

# ===============================
# 7. 保存结果
# ===============================
# 保存 CSV
res_df = pd.DataFrame(results)
res_path = os.path.join(OUTPUT_DIR, "mamba_3000_feature_sweep.csv")
res_df.to_csv(res_path, index=False)
print(f"\n✅ 测试指标已保存: {res_path}")

# 保存最佳模型
if best_model_state is not None:
    model_path = os.path.join(OUTPUT_DIR, f"best_mamba_model_{best_feature_num}.pth")
    torch.save(best_model_state, model_path)
    print(f"🏆 全局最优模型已保存: {model_path}")
    print(f"   最优特征数: {best_feature_num}")
    print(f"   最优 R²: {best_overall_r2:.4f}")

print("\n🎉 脚本运行结束")