# import os
# import random
# import numpy as np
# import pandas as pd
# import torch
# import torch.nn as nn
# from torch.utils.data import DataLoader, TensorDataset
# from sklearn.metrics import mean_absolute_error, r2_score, mean_squared_error
# import matplotlib.pyplot as plt
# from tqdm import tqdm
# from mambapy.mamba import Mamba, MambaConfig
# from datetime import datetime
# import warnings
#
# warnings.filterwarnings("ignore")  # 忽略所有警告
#
#
# # ===========================
# # 1. 随机种子
# # ===========================
# def set_seed(seed=42):
#     random.seed(seed)
#     np.random.seed(seed)
#     torch.manual_seed(seed)
#     torch.cuda.manual_seed_all(seed)
#     torch.backends.cudnn.deterministic = True
#     torch.backends.cudnn.benchmark = False
#
#
# set_seed(42)
#
# # ===========================
# # 2. 设备
# # ===========================
# device = "cuda" if torch.cuda.is_available() else "cpu"
# print(f"✅ Using device: {device}")
#
#
# # ===========================
# # 3. 模型结构
# # ===========================
# class EnhancedMambaRegressor(nn.Module):
#     def __init__(self, input_dim, D=10):
#         super().__init__()
#         L = input_dim // D  # 动态计算L
#         self.linear = nn.Linear(input_dim, L * D)  # 线性层，将输入映射到L*D维度
#         self.norm = nn.LayerNorm(L * D)  # 层归一化
#         self.dropout = nn.Dropout(0.2)  # Dropout正则化
#         self.mamba = Mamba(MambaConfig(d_model=D, n_layers=2, d_state=8, d_conv=2, pscan=True))  # Mamba模块
#         self.pool = nn.AdaptiveAvgPool1d(1)  # 自适应池化
#         self.fc = nn.Sequential(
#             nn.Linear(D, 32),
#             nn.ReLU(),
#             nn.Dropout(0.2),
#             nn.Linear(32, 1)
#         )
#
#     def forward(self, x):
#         x = self.linear(x)  # 线性映射
#         x = self.norm(x)  # 层归一化
#         x = self.dropout(x)  # Dropout
#         L = x.size(1) // self.fc[0].in_features  # 获取L维度
#         x = x.view(x.size(0), L, self.fc[0].in_features)  # 重新reshape为 (batch_size, L, D)
#         x = self.mamba(x)  # Mamba模块
#         x = x.transpose(1, 2)  # 转置
#         x = self.pool(x).squeeze(-1)  # 池化并去掉多余的维度
#         return self.fc(x).squeeze(-1)  # 输出
#
#
# # ===========================
# # 4. 数据加载
# # ===========================
# def load_cpg_matrix_by_split(baseline_file, unique_file, cpg_list):
#     df_base = pd.read_csv(baseline_file)
#     df_unique = pd.read_csv(unique_file)
#     if "split" not in df_base.columns:
#         raise KeyError("❌ baseline 文件中没有 split 列")
#
#     cpg_list = [c for c in cpg_list if c.startswith("cg")]
#     X_list = []
#     for cpg in cpg_list:
#         if cpg in df_base.columns:
#             X_list.append(df_base[cpg])
#         elif cpg in df_unique.columns:
#             X_list.append(df_unique[cpg])
#         else:
#             raise ValueError(f"CpG {cpg} 不在 baseline 或 unique 文件中！")
#
#     X_all = pd.concat(X_list, axis=1).values
#     y_all = df_base["Age"].values
#     split_all = df_base["split"].values
#     train_idx = split_all == "train"
#     val_idx = split_all == "val"
#
#     X_train = torch.tensor(X_all[train_idx], dtype=torch.float32).to(device)
#     y_train = torch.tensor(y_all[train_idx], dtype=torch.float32).to(device)
#     X_val = torch.tensor(X_all[val_idx], dtype=torch.float32).to(device)
#     y_val = torch.tensor(y_all[val_idx], dtype=torch.float32).to(device)
#     return X_train, y_train, X_val, y_val
#
#
# # ===========================
# # 5. 训练+验证 (混合精度与早停)
# # ===========================
# def train_and_eval(X_train, y_train, X_val, y_val, D=10, epochs=100, patience=15):
#     input_dim = X_train.shape[1]  # 输入特征的维度
#     model = EnhancedMambaRegressor(input_dim, D).to(device)
#     optimizer = torch.optim.AdamW(model.parameters(), lr=0.001)
#     loss_fn = nn.MSELoss()
#     loader = DataLoader(TensorDataset(X_train, y_train), batch_size=8, shuffle=True)
#     scaler = torch.cuda.amp.GradScaler()
#
#     best_val_mae = float('inf')
#     best_state = None
#     best_r2 = None
#     best_r = None
#     best_rmse = None
#
#     val_loss_patience = 0  # 记录验证损失没有下降的轮次
#     for epoch in range(epochs):
#         model.train()
#         for xb, yb in loader:
#             optimizer.zero_grad()
#             with torch.cuda.amp.autocast():
#                 pred = model(xb)
#                 loss = loss_fn(pred, yb)
#             scaler.scale(loss).backward()
#             scaler.step(optimizer)
#             scaler.update()
#
#         # 每个训练周期结束后输出一次日志
#         model.eval()
#         with torch.no_grad():
#             with torch.cuda.amp.autocast():
#                 pred_val = model(X_val).cpu().numpy()
#             y_true = y_val.cpu().numpy()
#
#             val_mae = mean_absolute_error(y_true, pred_val)
#             val_r2 = r2_score(y_true, pred_val)
#             val_r = np.corrcoef(y_true, pred_val)[0, 1]  # 计算相关系数
#             val_rmse = np.sqrt(mean_squared_error(y_true, pred_val))  # 计算RMSE
#
#             # 更新最优模型
#             if val_mae < best_val_mae:
#                 best_val_mae = val_mae
#                 best_r2 = val_r2
#                 best_r = val_r
#                 best_rmse = val_rmse
#                 best_state = model.state_dict()
#
#             # 早停机制
#             val_loss_patience += 1
#             if val_loss_patience >= patience:  # 早停
#                 print(f"Early stopping after {epoch + 1} epochs")
#                 break
#
#     return best_val_mae, best_r2, best_r, best_rmse, best_state
#
#
# # ===========================
# # 6. 只基于年龄相关性选择位点
# # ===========================
# def select_top_cpgs_based_on_correlation(cpg_info_file, baseline):
#     cpg_df = pd.read_csv(cpg_info_file)
#     ranked = cpg_df[~cpg_df["CpG_ID"].isin(baseline)]  # 选择不在baseline中的CpGs
#     ranked_sorted = ranked.sort_values(by="correlation", ascending=False)
#     return ranked_sorted["CpG_ID"].tolist()
#
#
# # ===========================
# # 7. 保存最佳甲基化矩阵
# # ===========================
# def save_best_methylation_matrix(baseline_file, unique_file, best_cpgs, output_file):
#     # 读取数据
#     df_base = pd.read_csv(baseline_file)
#     df_unique = pd.read_csv(unique_file)
#
#     # 清理列名，去除前后的空格，并统一转换为小写
#     df_base.columns = df_base.columns.str.strip().str.lower()
#     df_unique.columns = df_unique.columns.str.strip().str.lower()
#
#     # 打印列名，检查是否存在目标CpG位点
#     print("Columns in df_base:", df_base.columns)
#     print("Columns in df_unique:", df_unique.columns)
#
#     # 确保所有位点都在列中
#     missing_cpgs = [cpg.lower() for cpg in best_cpgs if
#                     cpg.lower() not in df_base.columns and cpg.lower() not in df_unique.columns]
#     if missing_cpgs:
#         print(f"❌ Missing CpGs: {missing_cpgs}")
#         return
#
#     # 创建包含最佳位点的矩阵
#     selected_columns = ["id", "age", "split", "dataset"] + [cpg.lower() for cpg in best_cpgs]
#
#     # 合并数据
#     df_best = pd.concat([df_base[selected_columns], df_unique[selected_columns]], axis=0)
#
#     # 保存最终文件
#     df_best.to_csv(output_file, index=False)
#     print(f"✅ Best methylation matrix saved to {output_file}")
#
#
# # ===========================
# # 8. 主流程
# # ===========================
# def main():
#     global baseline_file, unique_file, baseline, cpg_info_file, output_dir
#     baseline_file = "/datapool/home/info_wang/wanghui/file/matrix_138.csv"
#     unique_file = "/datapool/home/info_wang/wanghui/file/matrix_combined_unique.csv"
#     cpg_info_file = "/datapool/home/info_wang/wanghui/file/all_cpgs_partial_corr.csv"
#     output_dir = "/datapool/home/info_wang/wanghui/file/cpg_weighted_expansion/"
#
#     os.makedirs(output_dir, exist_ok=True)
#
#     base_cols = pd.read_csv(baseline_file, nrows=0).columns.tolist()
#     baseline = [c for c in base_cols if c.startswith("cg")]
#
#     # 根据年龄相关性选择Top CpG位点
#     top_cpgs = select_top_cpgs_based_on_correlation(cpg_info_file, baseline)
#
#     # 定义保存结果的列表
#     mae_list = []
#     r2_list = []
#     r_list = []
#     rmse_list = []
#     best_mae = float('inf')
#     best_state = None
#     best_cpgs = []
#
#     # 循环增加位点，按步长10增加
#     for i in tqdm(range(40), desc="Progress", unit="step"):
#         k = (i + 1) * 10  # 每步增加10个位点
#         current_cpgs = top_cpgs[:k]
#         print(f"\nTesting {k} CpG sites...")
#
#         # 加载数据并训练模型
#         X_train, y_train, X_val, y_val = load_cpg_matrix_by_split(baseline_file, unique_file, baseline + current_cpgs)
#         mae, r2, r, rmse, best_state = train_and_eval(X_train, y_train, X_val, y_val, D=10)
#
#         # 保存结果
#         mae_list.append(mae)
#         r2_list.append(r2)
#         r_list.append(r)
#         rmse_list.append(rmse)
#
#         # 打印四个指标，保留四位小数
#         print(f"MAE: {mae:.4f}, R²: {r2:.4f}, R: {r:.4f}, RMSE: {rmse:.4f}")
#
#         # 保存当前最优的模型（如果当前结果最好）
#         if mae < best_mae:
#             best_mae = mae
#             best_state = best_state
#             best_cpgs = current_cpgs
#
#     # =======================
#     # 绘制评估结果曲线图
#     # =======================
#     plt.figure(figsize=(10, 6))
#     plt.plot(range(1, 41), mae_list, label='MAE', marker='o')
#     plt.plot(range(1, 41), r2_list, label='R²', marker='o')
#     plt.xlabel('Number of CpG Sites')
#     plt.ylabel('Metric Value')
#     plt.title('MAE and R² over different CpG Sites')
#     plt.legend()
#     plt.grid(True)
#     plt.savefig(os.path.join(output_dir, 'R2&MAE_curve.png'))
#     plt.show()
#
#     # =======================
#     # 保存评估结果
#     # =======================
#     results_df = pd.DataFrame({
#         'CpG_Count': list(range(10, 401, 10)),
#         'MAE': mae_list,
#         'R²': r2_list,
#         'R': r_list,
#         'RMSE': rmse_list
#     })
#     results_df.to_csv(os.path.join(output_dir, 'shap10_results.csv'), index=False)
#     print(f"✅ Results saved to shap10_results.csv")
#
#     # =======================
#     # 保存最佳甲基化矩阵
#     # =======================
#     save_best_methylation_matrix(baseline_file, unique_file, best_cpgs,
#                                  os.path.join(output_dir, "best_methylation_matrix.csv"))
#     print(f"✅ Best methylation matrix saved to {output_dir}best_methylation_matrix.csv")
#
#     # =======================
#     # 保存最终模型
#     # =======================
#     torch.save(best_state, os.path.join(output_dir, "best_model.pth"))
#     print(f"✅ Best model saved to best_model.pth")
#
#
# if __name__ == "__main__":
#     main()
import os
import random
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, TensorDataset
from sklearn.metrics import mean_absolute_error, r2_score, mean_squared_error
import matplotlib.pyplot as plt
from mambapy.mamba import Mamba, MambaConfig
from tqdm import tqdm
import warnings

warnings.filterwarnings("ignore")


# ===========================
# 1. 基础配置
# ===========================
def set_seed(seed=42):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True


set_seed(42)
device = "cuda" if torch.cuda.is_available() else "cpu"
print(f"✅ Using device: {device}")


# ===========================
# 2. 动态 Mamba 模型
# ===========================
class DynamicMambaRegressor(nn.Module):
    def __init__(self, input_dim, d_model=64):
        super().__init__()
        self.embedding = nn.Linear(input_dim, d_model)
        self.norm = nn.LayerNorm(d_model)
        self.dropout = nn.Dropout(0.2)
        config = MambaConfig(d_model=d_model, n_layers=2, d_state=8, d_conv=2, pscan=True)
        self.mamba = Mamba(config)
        self.head = nn.Sequential(
            nn.Linear(d_model, 32),
            nn.ReLU(),
            nn.Linear(32, 1)
        )

    def forward(self, x):
        x = self.embedding(x)
        x = self.norm(x)
        x = self.dropout(x)
        x = x.unsqueeze(1)
        x = self.mamba(x)
        x = x.mean(dim=1)
        return self.head(x).squeeze(-1)


# ===========================
# 3. 数据准备与排序
# ===========================
def prepare_data_with_ranking(baseline_file, unique_file, corr_file):
    print("⏳ 读取数据并排序...")
    df_base = pd.read_csv(baseline_file, index_col=0)
    df_unique = pd.read_csv(unique_file, index_col=0)

    common_idx = df_base.index.intersection(df_unique.index)
    if len(common_idx) == 0:
        raise ValueError("❌ 错误：两个矩阵的样本ID没有重合！")

    df_base = df_base.loc[common_idx]
    df_unique = df_unique.loc[common_idx]

    y = df_base['Age'].values
    splits = df_base['Split'].values

    print(f"📖 读取相关性文件: {os.path.basename(corr_file)}")
    df_corr = pd.read_csv(corr_file)

    cols = df_corr.columns
    id_col = next((c for c in cols if c.lower() in ['cpg_id', 'cpg', 'id', 'feature']), df_corr.columns[0])
    numeric_cols = df_corr.select_dtypes(include=[np.number]).columns
    corr_col = numeric_cols[0] if len(numeric_cols) > 0 else None

    if corr_col is None:
        raise ValueError("❌ 在相关性文件中找不到数值列！")

    df_corr['abs_corr'] = df_corr[corr_col].abs()
    df_sorted = df_corr.sort_values(by='abs_corr', ascending=False)

    unique_cols = set(df_unique.columns)
    sorted_candidates = [c for c in df_sorted[id_col] if c in unique_cols]
    base_feats = [c for c in df_base.columns if c.startswith('cg') or c.startswith('ch')]

    print(f"✅ 排序完成: 基准={len(base_feats)}, 候选={len(sorted_candidates)}")
    return df_base, df_unique, base_feats, sorted_candidates, y, splits


# ===========================
# 4. 训练评估函数 (修改：返回模型状态)
# ===========================
def train_evaluate(X_train, y_train, X_val, y_val, epochs=150, patience=15):
    model = DynamicMambaRegressor(input_dim=X_train.shape[1]).to(device)
    optimizer = torch.optim.AdamW(model.parameters(), lr=1e-4)
    criterion = nn.MSELoss()

    loader = DataLoader(TensorDataset(X_train, y_train), batch_size=16, shuffle=True)

    # 早停相关变量
    best_val_mae = float('inf')
    best_metrics = (0, 0, 0, 0)  # (mae, r2, rmse, r)
    best_model_state = None
    patience_counter = 0

    for epoch in range(epochs):
        # --- 训练阶段 ---
        model.train()
        for xb, yb in loader:
            optimizer.zero_grad()
            pred = model(xb)
            loss = criterion(pred, yb)
            loss.backward()
            optimizer.step()

        # --- 验证阶段 ---
        model.eval()
        with torch.no_grad():
            val_pred_t = model(X_val)
            val_pred = val_pred_t.cpu().numpy()
            y_true = y_val.cpu().numpy()

            # 计算当前轮次的指标
            mae = mean_absolute_error(y_true, val_pred)
            r2 = r2_score(y_true, val_pred)
            rmse = np.sqrt(mean_squared_error(y_true, val_pred))
            r_corr = np.corrcoef(y_true, val_pred)[0, 1]

        # --- 早停逻辑 ---
        if mae < best_val_mae:
            best_val_mae = mae
            # 暂存当前最佳的所有指标和模型参数
            best_metrics = (mae, r2, rmse, r_corr)
            best_model_state = model.state_dict()
            patience_counter = 0  # 重置计数器
        else:
            patience_counter += 1

        if patience_counter >= patience:
            # 既然已经连续 N 轮没提升了，直接停止，无需打印
            break

    # 如果训练完了还没有 best_model_state (极端情况), 用最后的
    if best_model_state is None:
        best_model_state = model.state_dict()
        best_metrics = (mae, r2, rmse, r_corr)

    # 返回最佳时刻的指标和参数
    return *best_metrics, best_model_state


# ===========================
# 5. 主流程
# ===========================
def main():
    baseline_file = "/datapool/home/info_wang/wanghui/file/matrix_138.csv"
    unique_file = "/datapool/home/info_wang/wanghui/file/matrix_unique_combined.csv"
    corr_file = "/datapool/home/info_wang/wanghui/file/all_cpgs_partial_corr.csv"

    output_dir = "/datapool/home/info_wang/wanghui/file/CPG"
    if not os.path.exists(output_dir): os.makedirs(output_dir)

    df_base, df_unique, base_feats, sorted_candidates, y, splits = \
        prepare_data_with_ranking(baseline_file, unique_file, corr_file)

    data_dict = {
        'base': df_base[base_feats].values,
        'cand_df': df_unique,
        'y': torch.tensor(y, dtype=torch.float32).to(device),
        'train_mask': splits == 'train',
        'val_mask': splits == 'val'
    }

    results = []
    best_mae = float('inf')
    best_k = 0
    best_model_weights = None  # 用于暂存最佳模型参数

    max_k = 200
    step = 10
    steps = list(range(0, max_k + 1, step))

    print(f"\n🚀 开始筛选 (Max={max_k}, Step={step})...")

    for k in tqdm(steps):
        if k == 0:
            X_curr = data_dict['base']
        else:
            top_k_cols = sorted_candidates[:k]
            X_cand = data_dict['cand_df'][top_k_cols].values
            X_curr = np.hstack([data_dict['base'], X_cand])

        X_tensor = torch.tensor(X_curr, dtype=torch.float32).to(device)

        # ... Inside main loop ...

        # 训练并获取指标 (epochs 可以设大一点，比如 100，利用早停来控制)
        mae, r2, rmse, r_val, model_weights = train_evaluate(
            X_tensor[data_dict['train_mask']],
            data_dict['y'][data_dict['train_mask']],
            X_tensor[data_dict['val_mask']],
            data_dict['y'][data_dict['val_mask']],
            epochs=150,  # 增加最大轮数
            patience=15  # 设置耐心值
        )

        results.append({
            'Added_Count': k,
            'Total_Features': len(base_feats) + k,
            'MAE': mae,
            'R2': r2,
            'RMSE': rmse,
            'R': r_val
        })

        # 如果是最佳结果，保存权重到内存
        if mae < best_mae:
            best_mae = mae
            best_k = k
            best_model_weights = model_weights  # 暂存参数

        if k % 50 == 0:
            tqdm.write(f"   [+{k}] MAE:{mae:.4f} | R2:{r2:.4f} | Best k={best_k}")

    # 1. 保存指标文件
    res_df = pd.DataFrame(results)
    metrics_file = os.path.join(output_dir, "stepwise_metrics_full.csv")
    res_df.to_csv(metrics_file, index=False)
    print(f"\n✅ 指标记录: {metrics_file}")

    print(f"\n🏆 筛选结束！")
    print(f"   最佳结果: 基准 + Top {best_k} 个候选")
    print(f"   最佳MAE: {best_mae:.4f}")

    # ===========================
    # 结果生成
    # ===========================
    print("\n💾 正在生成最终文件...")

    final_candidate_cols = sorted_candidates[:best_k]
    all_best_features = base_feats + final_candidate_cols

    # A. 保存最佳位点列表 (CSV)
    features_df = pd.DataFrame({'CpG': all_best_features})
    features_list_path = os.path.join(output_dir, "best_features_list.csv")
    features_df.to_csv(features_list_path, index=False)
    print(f"   - 位点列表: {features_list_path}")

    # B. 保存最佳甲基化矩阵
    meta_cols = [c for c in df_base.columns if c in ['Age', 'Split', 'Dataset']]
    df_final = df_base[meta_cols].join(df_base[base_feats])
    if best_k > 0:
        df_cand_best = df_unique[final_candidate_cols]
        df_final = df_final.join(df_cand_best)

    matrix_path = os.path.join(output_dir, "best_performance_matrix.csv")
    df_final.to_csv(matrix_path)
    print(f"   - 最佳矩阵: {matrix_path}")

    # C. 保存最佳模型 (新增)
    if best_model_weights is not None:
        model_path = os.path.join(output_dir, "best_model.pth")
        torch.save(best_model_weights, model_path)
        print(f"   - 最佳模型: {model_path}")
        print("     (注意: 加载时模型 input_dim 需设置为: {})".format(len(all_best_features)))

    # D. 绘图
    plt.figure(figsize=(10, 6))
    plt.plot(res_df['Added_Count'], res_df['MAE'], marker='o', label='MAE')
    plt.axvline(x=best_k, color='r', linestyle='--', label=f'Best k={best_k}')
    plt.title('MAE Optimization Curve')
    plt.xlabel('Added Features')
    plt.ylabel('MAE')
    plt.legend()
    plt.savefig(os.path.join(output_dir, "optimization_curve.png"))
    print("   - 曲线图已保存")


if __name__ == "__main__":
    main()