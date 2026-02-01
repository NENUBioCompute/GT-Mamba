import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset
from torch.optim.lr_scheduler import ReduceLROnPlateau
from mambapy.mamba import Mamba, MambaConfig
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score, mean_absolute_error, mean_squared_error
import warnings
import os
import random

warnings.filterwarnings("ignore")

# ===========================
# 0. 全局配置
# ===========================
CONFIG = {
    'GRAPH_PATH': '/datapool/home/info_wang/wanghui/file/graph/methylation_graph(5).pt',
    'DATA_PATH': "/datapool/home/info_wang/wanghui/file/CPG/methylation_matrix_198.csv",
    'SAVE_DIR': '/datapool/home/info_wang/wanghui/file/graph',
    'DEVICE': torch.device("cuda" if torch.cuda.is_available() else "cpu"),
    'LR': 1e-4,
    'BATCH_SIZE': 8,
    'EPOCHS': 400,
    'PATIENCE': 30,
    'D_MODEL': 64,
    'D_STATE': 16,
    'D_CONV': 4
}

print(f"✅ Using device: {CONFIG['DEVICE']}")


# ===========================
# 1. 随机种子设置
# ===========================
def set_seed(seed=42):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = False


set_seed(42)
print("✅ 随机种子已固定: 42")


# ===========================
# 2. 数据加载与预处理
# ===========================
def load_data():
    print("⏳ Loading data...")
    if not os.path.exists(CONFIG['GRAPH_PATH']):
        raise FileNotFoundError(f"找不到图文件: {CONFIG['GRAPH_PATH']}")
    graph_data = torch.load(CONFIG['GRAPH_PATH'])

    if not os.path.exists(CONFIG['DATA_PATH']):
        raise FileNotFoundError(f"找不到数据文件: {CONFIG['DATA_PATH']}")
    df = pd.read_csv(CONFIG['DATA_PATH'])

    if df.isnull().values.any():
        print("⚠️ 警告: 数据包含 NaN，正在填充为 0")
        df = df.fillna(0)

    meta_cols = ['age', 'split', 'dataset', 'id', 'Age', 'Split', 'Dataset', 'ID']
    feat_cols = [c for c in df.columns if c not in meta_cols]

    X = df[feat_cols].values
    y = df['Age'].values
    split = df['Split'].values

    X_train = torch.tensor(X[split == 'train'], dtype=torch.float32)
    y_train = torch.tensor(y[split == 'train'], dtype=torch.float32)
    X_val = torch.tensor(X[split == 'val'], dtype=torch.float32)
    y_val = torch.tensor(y[split == 'val'], dtype=torch.float32)
    X_test = torch.tensor(X[split == 'test'], dtype=torch.float32)
    y_test = torch.tensor(y[split == 'test'], dtype=torch.float32)

    print(f"  特征数量 (Input Dim): {X_train.shape[1]}")
    return X_train, X_val, X_test, y_train, y_val, y_test, graph_data


# ===========================
# 3. 绘图函数 (新增：论文级散点图)
# ===========================
def plot_results(y_true, y_pred, save_dir):
    """绘制专业散点图并保存"""
    r2 = r2_score(y_true, y_pred)
    mae = mean_absolute_error(y_true, y_pred)
    rmse = np.sqrt(mean_squared_error(y_true, y_pred))
    r = np.corrcoef(y_true, y_pred)[0, 1]

    plt.figure(figsize=(6, 6))
    ax = plt.gca()

    # 绘制散点 (深蓝色)
    plt.scatter(y_true, y_pred, alpha=0.6, color='#1f77b4', edgecolors='white', linewidth=0.6, s=35)

    # 绘制对角线
    data_min = min(y_true.min(), y_pred.min())
    data_max = max(y_true.max(), y_pred.max())
    plot_min = data_min - 2
    plot_max = data_max + 2
    plt.plot([plot_min, plot_max], [plot_min, plot_max], 'r--', lw=1.5, label='Identity Line')

    plt.xlabel('True Age (Years)', fontsize=12)
    plt.ylabel('Predicted Age (Years)', fontsize=12)
    plt.title('Chronological vs. Predicted Age', fontsize=13)

    # 添加指标文本框
    metrics_str = (f'$R^2 = {r2:.4f}$\n'
                   f'$MAE = {mae:.4f}$\n'
                   f'$RMSE = {rmse:.4f}$\n'
                   f'$r = {r:.4f}$')

    plt.text(0.05, 0.95, metrics_str, transform=ax.transAxes, fontsize=11,
             verticalalignment='top', horizontalalignment='left',
             bbox=dict(boxstyle='round,pad=0.4', facecolor='white', alpha=0.9, edgecolor='gray'))

    plt.xlim(plot_min, plot_max)
    plt.ylim(plot_min, plot_max)
    plt.grid(True, linestyle=':', alpha=0.5)
    plt.legend(loc='lower right', fontsize=10, frameon=True)

    save_path = os.path.join(save_dir, 'scatter_result(5).png')
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"📈 散点图已保存: {save_path}")


# ===========================
# 4. 模型定义 (GAT + Mamba)
# ===========================
class GraphMaskedAttention(nn.Module):
    def __init__(self, input_dim, num_heads=4, dropout=0.1):#修改dropout0.1-0.25
        super().__init__()
        self.num_heads = num_heads
        self.head_dim = input_dim // num_heads
        self.scale = self.head_dim ** -0.5

        self.q_proj = nn.Linear(input_dim, input_dim)
        self.k_proj = nn.Linear(input_dim, input_dim)
        self.v_proj = nn.Linear(input_dim, input_dim)
        self.out = nn.Linear(input_dim, input_dim)
        self.dropout = nn.Dropout(dropout)
        self.norm = nn.LayerNorm(input_dim)

    def forward(self, x, adj_mask):
        B, N, C = x.shape
        residual = x
        q = self.q_proj(x).view(B, N, self.num_heads, self.head_dim).transpose(1, 2)
        k = self.k_proj(x).view(B, N, self.num_heads, self.head_dim).transpose(1, 2)
        v = self.v_proj(x).view(B, N, self.num_heads, self.head_dim).transpose(1, 2)

        scores = (q @ k.transpose(-2, -1)) * self.scale
        if adj_mask is not None:
            scores = scores + adj_mask.unsqueeze(0).unsqueeze(0)

        attn = torch.softmax(scores, dim=-1)
        attn = self.dropout(attn)
        out = (attn @ v).transpose(1, 2).reshape(B, N, C)
        out = self.out(out)
        out = self.dropout(out)
        return self.norm(residual + out)


class GraphMambaRegressor(nn.Module):
    def __init__(self, num_nodes, graph_edge_index=None):
        super().__init__()
        d_model = CONFIG['D_MODEL']

        mask = torch.full((num_nodes, num_nodes), float('-inf'))
        mask.fill_diagonal_(0.0)
        if graph_edge_index is not None:
            edge_index = graph_edge_index.long()
            mask[edge_index[0], edge_index[1]] = 0.0
        self.register_buffer('adj_mask', mask)

        self.node_embedding = nn.Linear(1, d_model)
        self.norm_in = nn.LayerNorm(d_model)
        self.gat = GraphMaskedAttention(d_model, num_heads=4)

        config = MambaConfig(d_model=d_model, n_layers=2, d_state=CONFIG['D_STATE'], d_conv=CONFIG['D_CONV'],
                             pscan=True)
        self.mamba = Mamba(config)

        self.head = nn.Sequential(
            nn.LayerNorm(d_model),
            nn.Dropout(0.1),#修改dropout0.1-0.3
            nn.Linear(d_model, 32),
            nn.GELU(),
            nn.Linear(32, 1)
        )

    def forward(self, x):
        x = x.unsqueeze(-1)
        x = self.node_embedding(x)
        x = self.norm_in(x)
        x = self.gat(x, self.adj_mask)
        x = self.mamba(x)
        x = x.mean(dim=1)
        return self.head(x).squeeze(-1)


# ===========================
# 5. 训练与评估流程
# ===========================
def train_and_evaluate(model, train_loader, val_loader, test_loader):
    model = model.to(CONFIG['DEVICE'])
    optimizer = optim.AdamW(model.parameters(), lr=CONFIG['LR'], weight_decay=1e-4)
    scheduler = ReduceLROnPlateau(optimizer, 'min', patience=10, factor=0.5, verbose=True)
    criterion = nn.MSELoss()

    best_val_loss = float('inf')
    train_losses, val_losses = [], []
    patience_counter = 0

    print(f"\n🚀 开始训练 (Epochs={CONFIG['EPOCHS']})...")

    for epoch in range(CONFIG['EPOCHS']):
        model.train()
        batch_losses = []
        for X_b, y_b in train_loader:
            X_b, y_b = X_b.to(CONFIG['DEVICE']), y_b.to(CONFIG['DEVICE'])
            optimizer.zero_grad()
            pred = model(X_b)
            loss = criterion(pred, y_b)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            optimizer.step()
            batch_losses.append(loss.item())

        train_loss = np.mean(batch_losses)
        train_losses.append(train_loss)

        model.eval()
        val_batch_losses = []
        with torch.no_grad():
            for X_b, y_b in val_loader:
                X_b, y_b = X_b.to(CONFIG['DEVICE']), y_b.to(CONFIG['DEVICE'])
                pred = model(X_b)
                val_batch_losses.append(criterion(pred, y_b).item())

        val_loss = np.mean(val_batch_losses)
        val_losses.append(val_loss)
        scheduler.step(val_loss)

        if (epoch + 1) % 10 == 0:
            print(f"Epoch {epoch + 1:3d} | Train Loss: {train_loss:.4f} | Val Loss: {val_loss:.4f}")

        if val_loss < best_val_loss:
            best_val_loss = val_loss
            torch.save(model.state_dict(), os.path.join(CONFIG['SAVE_DIR'], 'best_model(5).pth'))
            patience_counter = 0
        else:
            patience_counter += 1
            if patience_counter >= CONFIG['PATIENCE']:
                print(f"🛑 早停触发 (Epoch {epoch + 1})")
                break

    # 绘制 Loss 曲线
    plt.figure(figsize=(10, 6))
    plt.plot(train_losses, label='Train Loss')
    plt.plot(val_losses, label='Val Loss')
    plt.title('Training Dynamics')
    plt.xlabel('Epochs')
    plt.ylabel('MSE Loss')
    plt.legend()
    plt.savefig(os.path.join(CONFIG['SAVE_DIR'], 'loss_curve(5).png'))
    print("📈 Loss 曲线已保存")

    # --- 最终测试 ---
    print("\n🔍 正在评估测试集...")
    model.load_state_dict(torch.load(os.path.join(CONFIG['SAVE_DIR'], 'best_model(5).pth')))
    model.eval()

    y_true, y_pred = [], []
    with torch.no_grad():
        for X_b, y_b in test_loader:
            X_b = X_b.to(CONFIG['DEVICE'])
            pred = model(X_b)
            y_true.extend(y_b.cpu().numpy())
            y_pred.extend(pred.cpu().numpy())

    y_true = np.array(y_true)
    y_pred = np.array(y_pred)

    # 指标计算
    mae = mean_absolute_error(y_true, y_pred)
    r2 = r2_score(y_true, y_pred)
    rmse = np.sqrt(mean_squared_error(y_true, y_pred))
    r = np.corrcoef(y_true, y_pred)[0, 1]

    print("=" * 40)
    print(f"最终测试结果:")
    print(f"  MAE : {mae:.4f}")
    print(f"  R²  : {r2:.4f}")
    print(f"  RMSE: {rmse:.4f}")
    print(f"  R   : {r:.4f}")
    print("=" * 40)

    # [关键] 调用新的绘图函数
    plot_results(y_true, y_pred, CONFIG['SAVE_DIR'])


# ===========================
# 6. 主程序
# ===========================
def main():
    X_train, X_val, X_test, y_train, y_val, y_test, graph_data = load_data()
    train_loader = DataLoader(TensorDataset(X_train, y_train), batch_size=CONFIG['BATCH_SIZE'], shuffle=True)
    val_loader = DataLoader(TensorDataset(X_val, y_val), batch_size=CONFIG['BATCH_SIZE'])
    test_loader = DataLoader(TensorDataset(X_test, y_test), batch_size=CONFIG['BATCH_SIZE'])

    try:
        edge_index = graph_data.edge_index
        print(f"🔗 图结构加载成功: {edge_index.shape[1]} 条边")
    except AttributeError:
        print("⚠️ 警告: 图文件中未找到 edge_index")
        edge_index = None

    model = GraphMambaRegressor(num_nodes=X_train.shape[1], graph_edge_index=edge_index)
    train_and_evaluate(model, train_loader, val_loader, test_loader)


if __name__ == "__main__":
    main()
