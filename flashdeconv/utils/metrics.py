"""
Evaluation metrics for deconvolution results.

This module provides functions to evaluate the accuracy of
cell type proportion estimates.
"""

import numpy as np
from typing import Optional, Tuple


def compute_rmse(
    pred: np.ndarray,
    true: np.ndarray,
    per_cell_type: bool = False,
) -> np.ndarray:
    """
    Compute Root Mean Squared Error between predicted and true proportions.

    Parameters
    ----------
    pred : ndarray of shape (n_spots, n_cell_types)
        Predicted proportions.
    true : ndarray of shape (n_spots, n_cell_types)
        True proportions.
    per_cell_type : bool, default=False
        If True, return RMSE per cell type.

    Returns
    -------
    rmse : float or ndarray
        RMSE value(s).
    """
    squared_errors = (pred - true) ** 2

    if per_cell_type:
        return np.sqrt(np.mean(squared_errors, axis=0))
    else:
        return np.sqrt(np.mean(squared_errors))


def compute_mae(
    pred: np.ndarray,
    true: np.ndarray,
    per_cell_type: bool = False,
) -> np.ndarray:
    """
    Compute Mean Absolute Error.

    Parameters
    ----------
    pred : ndarray of shape (n_spots, n_cell_types)
        Predicted proportions.
    true : ndarray of shape (n_spots, n_cell_types)
        True proportions.
    per_cell_type : bool, default=False
        If True, return MAE per cell type.

    Returns
    -------
    mae : float or ndarray
        MAE value(s).
    """
    abs_errors = np.abs(pred - true)

    if per_cell_type:
        return np.mean(abs_errors, axis=0)
    else:
        return np.mean(abs_errors)


def compute_correlation(
    pred: np.ndarray,
    true: np.ndarray,
    method: str = "pearson",
    per_cell_type: bool = False,
) -> np.ndarray:
    """
    Compute correlation between predicted and true proportions.

    Parameters
    ----------
    pred : ndarray of shape (n_spots, n_cell_types)
        Predicted proportions.
    true : ndarray of shape (n_spots, n_cell_types)
        True proportions.
    method : str, default="pearson"
        Correlation method: "pearson" or "spearman".
    per_cell_type : bool, default=False
        If True, return correlation per cell type.

    Returns
    -------
    corr : float or ndarray
        Correlation value(s).
    """
    if method == "spearman":
        from scipy.stats import spearmanr

        def corr_func(x, y):
            return spearmanr(x, y)[0]
    else:
        def corr_func(x, y):
            return np.corrcoef(x, y)[0, 1]

    if per_cell_type:
        n_cell_types = pred.shape[1]
        corrs = np.zeros(n_cell_types)
        for k in range(n_cell_types):
            corrs[k] = corr_func(pred[:, k], true[:, k])
        return corrs
    else:
        return corr_func(pred.flatten(), true.flatten())


def compute_jsd(
    pred: np.ndarray,
    true: np.ndarray,
    epsilon: float = 1e-10,
) -> np.ndarray:
    """
    Compute Jensen-Shannon Divergence per spot.

    Parameters
    ----------
    pred : ndarray of shape (n_spots, n_cell_types)
        Predicted proportions (should sum to 1).
    true : ndarray of shape (n_spots, n_cell_types)
        True proportions (should sum to 1).
    epsilon : float, default=1e-10
        Small value to avoid log(0).

    Returns
    -------
    jsd : ndarray of shape (n_spots,)
        JSD per spot.
    """
    # Ensure valid probability distributions
    pred = np.clip(pred, epsilon, 1 - epsilon)
    true = np.clip(true, epsilon, 1 - epsilon)

    # Normalize
    pred = pred / pred.sum(axis=1, keepdims=True)
    true = true / true.sum(axis=1, keepdims=True)

    # Compute mixture
    m = 0.5 * (pred + true)

    # KL divergences
    kl_pred = np.sum(pred * np.log(pred / m), axis=1)
    kl_true = np.sum(true * np.log(true / m), axis=1)

    # JSD
    jsd = 0.5 * (kl_pred + kl_true)

    return jsd


def evaluate_deconvolution(
    pred: np.ndarray,
    true: np.ndarray,
    cell_type_names: Optional[np.ndarray] = None,
) -> dict:
    """
    Comprehensive evaluation of deconvolution results.

    Parameters
    ----------
    pred : ndarray of shape (n_spots, n_cell_types)
        Predicted proportions.
    true : ndarray of shape (n_spots, n_cell_types)
        True proportions.
    cell_type_names : ndarray, optional
        Names of cell types for reporting.

    Returns
    -------
    metrics : dict
        Dictionary containing various evaluation metrics.
    """
    n_cell_types = pred.shape[1]

    if cell_type_names is None:
        cell_type_names = [f"CellType_{i}" for i in range(n_cell_types)]

    metrics = {
        "overall": {
            "rmse": float(compute_rmse(pred, true)),
            "mae": float(compute_mae(pred, true)),
            "pearson": float(compute_correlation(pred, true, "pearson")),
            "spearman": float(compute_correlation(pred, true, "spearman")),
            "mean_jsd": float(np.mean(compute_jsd(pred, true))),
        },
        "per_cell_type": {},
    }

    # Per-cell-type metrics
    rmse_per = compute_rmse(pred, true, per_cell_type=True)
    mae_per = compute_mae(pred, true, per_cell_type=True)
    pearson_per = compute_correlation(pred, true, "pearson", per_cell_type=True)
    spearman_per = compute_correlation(pred, true, "spearman", per_cell_type=True)

    for k, name in enumerate(cell_type_names):
        metrics["per_cell_type"][name] = {
            "rmse": float(rmse_per[k]),
            "mae": float(mae_per[k]),
            "pearson": float(pearson_per[k]),
            "spearman": float(spearman_per[k]),
            "mean_proportion_true": float(np.mean(true[:, k])),
            "mean_proportion_pred": float(np.mean(pred[:, k])),
        }

    return metrics


def compute_rare_cell_detection(
    pred: np.ndarray,
    true: np.ndarray,
    threshold: float = 0.05,
) -> Tuple[float, float, float]:
    """
    Evaluate detection of rare cell types.

    Parameters
    ----------
    pred : ndarray of shape (n_spots, n_cell_types)
        Predicted proportions.
    true : ndarray of shape (n_spots, n_cell_types)
        True proportions.
    threshold : float, default=0.05
        Threshold for "present" vs "absent".

    Returns
    -------
    precision : float
        Precision for rare cell detection.
    recall : float
        Recall for rare cell detection.
    f1 : float
        F1 score.
    """
    # Identify "rare" entries (low but non-zero)
    rare_mask = (true > 0) & (true < threshold)

    if not np.any(rare_mask):
        return np.nan, np.nan, np.nan

    # Binarize predictions
    pred_present = pred > (threshold / 2)  # More lenient threshold for predictions

    # Compute metrics
    true_positive = np.sum(pred_present & rare_mask)
    false_positive = np.sum(pred_present & ~rare_mask & (true == 0))
    false_negative = np.sum(~pred_present & rare_mask)

    precision = true_positive / (true_positive + false_positive + 1e-10)
    recall = true_positive / (true_positive + false_negative + 1e-10)
    f1 = 2 * precision * recall / (precision + recall + 1e-10)

    return precision, recall, f1
