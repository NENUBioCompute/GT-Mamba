import pandas as pd
import os
import numpy as np
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
import seaborn as sns

# --- 🌟 Core Modification 1: Move theme settings globally to ensure all plots have grids ---
sns.set_style("whitegrid")

# --- Configuration Paths ---
SAVE_DIR = "/datapool/home/info_wang/wanghui/file/sota"
GEO_LIST = ["GSE40279", "GSE61496", "GSE72777", "GSE77445", "GSE132203"]


# --- 🌟 Core Modification 2: Calculate AgeAccel using the top-tier journal standard residual method ---
def calculate_residuals(df, pred_col, true_col="True_Age"):
    """
    Calculate residuals via simple linear regression (Predicted ~ True).
    This eliminates the underlying interference of chronological age on prediction bias.
    """
    valid_idx = df[[pred_col, true_col]].dropna().index
    if len(valid_idx) < 2:
        return pd.Series(index=df.index, dtype=float)

    y = df.loc[valid_idx, pred_col].values
    x = df.loc[valid_idx, true_col].values

    # 1D linear fit using Ordinary Least Squares (OLS): y = mx + c
    A = np.vstack([x, np.ones(len(x))]).T
    m, c = np.linalg.lstsq(A, y, rcond=None)[0]

    # Calculate residuals (Actual predicted value - Fitted expected value)
    residuals = y - (m * x + c)

    # Create a NaN Series of the same length as the original DataFrame, then fill in valid residuals
    res_series = pd.Series(index=df.index, dtype=float)
    res_series.loc[valid_idx] = residuals
    return res_series


# --- Plotting Function: Generate a 1x3 correlation Panel ---
def plot_supplementary_corr(df, gse_name, save_dir):
    plot_configs = [
        ('PhenoAge Accel', 'AgeAccel', 'PhenoAge_AgeAccel'),
        ('GrimAge2 Accel', 'AgeAccel', 'GrimAge2_AgeAccel'),
        ('DunedinPACE', 'AgeAccel', 'dunedinpace')
    ]

    # Check which columns actually exist in the dataframe
    active_configs = [p for p in plot_configs if p[1] in df.columns and p[2] in df.columns]
    if not active_configs: return None

    fig, axes = plt.subplots(1, len(active_configs), figsize=(6 * len(active_configs), 5.5))
    if len(active_configs) == 1: axes = [axes]

    # Note: Local sns.set_style("whitegrid") has been removed, unified via global settings above.

    for i, (label, col1, col2) in enumerate(active_configs):
        ax = axes[i]
        # Drop NaNs to prevent calculation errors
        valid = df[[col1, col2]].dropna()
        if len(valid) < 5: continue

        x, y = valid[col1], valid[col2]
        r, p = pearsonr(x, y)

        n_val = len(valid)

        # Plot: Regression line + 95% confidence interval
        sns.regplot(x=x, y=y, ax=ax,
                    scatter_kws={'alpha': 0.4, 'color': '#2c3e50', 's': 30},
                    line_kws={'color': '#e74c3c', 'lw': 2})

        # Annotate statistical metrics
        p_str = f"P < 0.001" if p < 0.001 else f"P = {p:.3f}"
        stats_text = f"Pearson $r$ = {r:.3f}\n{p_str}\n$N$ = {n_val}"

        ax.text(0.05, 0.95, stats_text, transform=ax.transAxes,
                fontsize=11, fontweight='bold', verticalalignment='top',
                bbox=dict(facecolor='white', alpha=0.8, edgecolor='gray', boxstyle='round,pad=0.5'))

        ax.set_title(f"{gse_name}: GT-Mamba vs {label}", fontsize=13, fontweight='bold')
        ax.set_xlabel("GT-Mamba Age Acceleration", fontsize=11)
        ax.set_ylabel(f"{label}", fontsize=11)

    plt.tight_layout()
    img_path = os.path.join(save_dir, f"Supp_Fig_Corr_{gse_name}.png")
    plt.savefig(img_path, dpi=300)
    plt.close()
    return img_path


# --- Main Execution Logic ---
print("=" * 85)
print("🚀 Calculating AgeAccel correlations (rigorous residual method) and generating HD plots with grids...")
print("=" * 85)

results = []

for gse in GEO_LIST:
    mamba_path = os.path.join(SAVE_DIR, f"{gse}_GT_Mamba_predictions.csv")
    sota_path = os.path.join(SAVE_DIR, f"{gse}_SOTA_predictions.csv")

    if not os.path.exists(mamba_path) or not os.path.exists(sota_path):
        print(f"⚠️ Skipping {gse}: Missing prediction files")
        continue

    # 1. Read data
    df_mamba = pd.read_csv(mamba_path, index_col=0)
    df_sota = pd.read_csv(sota_path, index_col=0)

    # 2. Extract SOTA metrics
    cols_to_keep = ['dnamphenoage', 'grimage2', 'dunedinpace', 'True_Age']
    available_cols = [c for c in cols_to_keep if c in df_sota.columns]
    df_sota_sub = df_sota[available_cols]

    # 3. Merge and align samples
    df_merged = df_mamba.join(df_sota_sub, rsuffix='_sota', how='inner')

    # 4. Core Application: Calculate AgeAccel for SOTA clocks using the residual method
    if 'dnamphenoage' in df_merged.columns:
        df_merged['PhenoAge_AgeAccel'] = calculate_residuals(df_merged, 'dnamphenoage', 'True_Age')
    if 'grimage2' in df_merged.columns:
        df_merged['GrimAge2_AgeAccel'] = calculate_residuals(df_merged, 'grimage2', 'True_Age')

    # 5. Statistical calculations
    res_dict = {'Dataset': gse, 'Samples': len(df_merged)}


    def get_corr_metrics(df, c1, c2):
        if c1 in df.columns and c2 in df.columns:
            valid = df[[c1, c2]].dropna()
            if len(valid) > 5:
                r, p = pearsonr(valid[c1], valid[c2])
                return round(r, 4), f"{p:.2e}"
        return np.nan, np.nan


    res_dict['R_vs_PhenoAge'], res_dict['P_vs_PhenoAge'] = get_corr_metrics(df_merged, 'AgeAccel', 'PhenoAge_AgeAccel')
    res_dict['R_vs_GrimAge2'], res_dict['P_vs_GrimAge2'] = get_corr_metrics(df_merged, 'AgeAccel', 'GrimAge2_AgeAccel')
    res_dict['R_vs_DunedinPACE'], res_dict['P_vs_DunedinPACE'] = get_corr_metrics(df_merged, 'AgeAccel', 'dunedinpace')

    results.append(res_dict)

    # 6. Generate plots
    img_path = plot_supplementary_corr(df_merged, gse, SAVE_DIR)
    if img_path:
        print(f"✅ {gse} Success: {len(df_merged)} samples processed, correlation plot saved")

# --- Summary of Results ---
final_df = pd.DataFrame(results)
print("\n" + "-" * 85)
print("📊 Summary Statistics (Pearson R, Residual-corrected):")
cols_show = ['Dataset', 'Samples', 'R_vs_PhenoAge', 'R_vs_GrimAge2', 'R_vs_DunedinPACE']
print(final_df[cols_show].to_string(index=False))
print("-" * 85)

save_csv = os.path.join(SAVE_DIR, "Comprehensive_Prognostic_Results.csv")
final_df.to_csv(save_csv, index=False)
print(f"🚀 All tasks completed! Plots and CSV have been saved to: {SAVE_DIR}")