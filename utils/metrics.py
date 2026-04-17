# utils/metrics.py (Academic Standard Version)

import numpy as np
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from scipy.stats import pearsonr


def calculate_metrics(y_true, y_pred):
    """
    Calculate a comprehensive set of regression metrics for epigenetic age prediction.

    Returns:
        dict: A dictionary containing MAE, RMSE, R2, and Pearson Correlation (r).
    """
    # Ensure inputs are numpy arrays
    y_true = np.array(y_true).flatten()
    y_pred = np.array(y_pred).flatten()

    # 1. Mean Absolute Error (MAE) - The primary metric for aging clocks
    mae = mean_absolute_error(y_true, y_pred)

    # 2. Root Mean Squared Error (RMSE)
    rmse = np.sqrt(mean_squared_error(y_true, y_pred))

    # 3. Coefficient of Determination (R²)
    r2 = r2_score(y_true, y_pred)

    # 4. Pearson Correlation Coefficient (r)
    # pearsonr returns (correlation, p-value)
    r_corr, p_val = pearsonr(y_true, y_pred)

    return {
        'MAE': mae,
        'RMSE': rmse,
        'R2': r2,
        'Pearson_r': r_corr,
        'p_value': p_val
    }


def print_metrics_report(metrics_dict, cohort_name="Target Cohort"):
    """
    Print a formatted performance report to the console.
    """
    print(f"\n" + "-" * 40)
    print(f"📊 Performance Report: {cohort_name}")
    print("-" * 40)
    print(f"MAE      : {metrics_dict['MAE']:.4f} years")
    print(f"RMSE     : {metrics_dict['RMSE']:.4f} years")
    print(f"R²       : {metrics_dict['R2']:.4f}")
    print(f"Pearson r: {metrics_dict['Pearson_r']:.4f} (p={metrics_dict['p_value']:.2e})")
    print("-" * 40 + "\n")

# --- Example Usage ---
# from utils.metrics import calculate_metrics, print_metrics_report
# metrics = calculate_metrics(y_true, y_pred)
# print_metrics_report(metrics, "GSE132203")