"""
Uncertainty quantification for FlashDeconv.

Three tiers of uncertainty estimation:

Tier 0 (always-on, zero cost):
    - Prediction entropy: Shannon entropy of proportion vector
    - Reconstruction residual: normalized per-spot fitting error

Tier 1 (analytical, cheap):
    - Hessian-diagonal Laplace approximation for per-type variance
    - Uses the BCD denominator (XtX[k,k] + lambda * n_neighbors)
    - Produces approximate confidence intervals via delta method

Tier 2 (bootstrap, optional):
    - Poisson parametric bootstrap on original counts
    - Warm-start BCD for fast convergence
    - Produces empirical percentile confidence intervals
"""

import numpy as np
from scipy import sparse
from typing import Dict, Optional, Tuple


def compute_entropy(proportions: np.ndarray) -> np.ndarray:
    """
    Shannon entropy of per-spot proportion vectors.

    Low entropy = confident (dominated by few types).
    High entropy = uncertain (spread across many types).
    Maximum = log(K) when all types have equal proportion.

    Parameters
    ----------
    proportions : ndarray of shape (n_spots, n_cell_types)

    Returns
    -------
    entropy : ndarray of shape (n_spots,)
        Shannon entropy in nats.
    """
    p = np.clip(proportions, 1e-10, 1.0)
    return -np.sum(p * np.log(p), axis=1)


def compute_reconstruction_residual(
    Y_sketch: np.ndarray,
    X_sketch: np.ndarray,
    beta: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Per-spot reconstruction error in sketch space.

    Parameters
    ----------
    Y_sketch : ndarray of shape (n_spots, sketch_dim)
    X_sketch : ndarray of shape (n_cell_types, sketch_dim)
    beta : ndarray of shape (n_spots, n_cell_types)

    Returns
    -------
    residual_ss : ndarray of shape (n_spots,)
        Sum of squared residuals per spot.
    residual_norm : ndarray of shape (n_spots,)
        Normalized residual: ||e_i||^2 / ||y_i||^2.
    """
    # Predicted: beta @ X_sketch -> (n_spots, sketch_dim)
    Y_hat = beta @ X_sketch
    residual = Y_sketch - Y_hat
    residual_ss = np.sum(residual ** 2, axis=1)

    y_ss = np.sum(Y_sketch ** 2, axis=1)
    y_ss = np.maximum(y_ss, 1e-10)
    residual_norm = residual_ss / y_ss

    return residual_ss, residual_norm


def compute_hessian_variance(
    beta: np.ndarray,
    XtX: np.ndarray,
    residual_ss: np.ndarray,
    sketch_dim: int,
    lambda_: float,
    n_neighbors: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Laplace approximation for per-type variance using Hessian diagonal.

    For each spot i and active cell type k (beta_ik > 0):
        Var(beta_ik) ≈ sigma2_i / (XtX[k,k] + lambda * n_neighbors_i)

    where sigma2_i is the per-spot noise variance estimated from residuals.

    For proportions (delta method):
        Var(p_ik) ≈ Var(beta_ik) / (sum_j beta_ij)^2

    Parameters
    ----------
    beta : ndarray of shape (n_spots, n_cell_types)
        Converged raw abundances.
    XtX : ndarray of shape (n_cell_types, n_cell_types)
        Gram matrix of sketched reference.
    residual_ss : ndarray of shape (n_spots,)
        Per-spot sum of squared residuals.
    sketch_dim : int
        Dimensionality of sketch space.
    lambda_ : float
        Spatial regularization strength used in the solve.
    n_neighbors : ndarray of shape (n_spots,)
        Number of spatial neighbors per spot.

    Returns
    -------
    var_beta : ndarray of shape (n_spots, n_cell_types)
        Approximate variance of raw abundances.
    var_prop : ndarray of shape (n_spots, n_cell_types)
        Approximate variance of proportions (via delta method).
    """
    n_spots, n_types = beta.shape
    active = beta > 0  # (n_spots, n_types)
    n_active = active.sum(axis=1)  # (n_spots,)

    # Per-spot noise variance: sigma2_i = RSS_i / (d - K_active_i)
    dof = np.maximum(sketch_dim - n_active, 1.0)
    sigma2 = residual_ss / dof  # (n_spots,)

    # Hessian diagonal per type: h_kk = XtX[k,k] + lambda * n_neighbors_i
    xtx_diag = np.diag(XtX)  # (n_types,)
    h_diag = xtx_diag[np.newaxis, :] + lambda_ * n_neighbors[:, np.newaxis]  # (n_spots, n_types)

    # Var(beta_ik) = sigma2_i / h_kk (only for active types)
    var_beta = np.zeros((n_spots, n_types), dtype=np.float64)
    var_beta[active] = (sigma2[:, np.newaxis] / h_diag)[active]

    # Delta method for proportions: Var(p_k) ≈ Var(beta_k) / (sum beta)^2
    beta_sum = np.sum(beta, axis=1, keepdims=True)
    beta_sum_sq = np.maximum(beta_sum ** 2, 1e-20)
    var_prop = var_beta / beta_sum_sq

    return var_beta, var_prop


def compute_confidence_scores(
    proportions: np.ndarray,
    var_prop: np.ndarray,
    alpha: float = 0.05,
) -> Dict[str, np.ndarray]:
    """
    Compute per-spot, per-type confidence metrics from variance estimates.

    Parameters
    ----------
    proportions : ndarray of shape (n_spots, n_cell_types)
    var_prop : ndarray of shape (n_spots, n_cell_types)
    alpha : float
        Significance level for confidence intervals. Default 0.05 (95% CI).

    Returns
    -------
    dict with keys:
        'ci_half_width': ndarray (n_spots, n_cell_types)
            Half-width of approximate confidence interval.
        'ci_lower': ndarray (n_spots, n_cell_types)
            Lower bound, clipped to [0, 1].
        'ci_upper': ndarray (n_spots, n_cell_types)
            Upper bound, clipped to [0, 1].
        'cv': ndarray (n_spots, n_cell_types)
            Coefficient of variation (std / mean).
        'detection_confident': ndarray (n_spots, n_cell_types) bool
            Whether lower CI bound > 0.01 (confidently detected).
    """
    from scipy.stats import norm
    z = norm.ppf(1 - alpha / 2)

    std_prop = np.sqrt(np.maximum(var_prop, 0))
    ci_half = z * std_prop

    ci_lower = np.clip(proportions - ci_half, 0, 1)
    ci_upper = np.clip(proportions + ci_half, 0, 1)

    # Coefficient of variation
    cv = np.zeros_like(proportions)
    nonzero = proportions > 1e-6
    cv[nonzero] = std_prop[nonzero] / proportions[nonzero]

    detection_confident = ci_lower > 0.01

    return {
        'ci_half_width': ci_half,
        'ci_lower': ci_lower,
        'ci_upper': ci_upper,
        'cv': cv,
        'detection_confident': detection_confident,
    }


def poisson_bootstrap(
    Y: np.ndarray,
    X: np.ndarray,
    coords: np.ndarray,
    model_params: Dict,
    beta_init: np.ndarray,
    n_bootstrap: int = 100,
    max_iter_boot: int = 20,
    seed: int = 42,
    verbose: bool = False,
) -> Dict[str, np.ndarray]:
    """
    Poisson parametric bootstrap for empirical confidence intervals.

    For each bootstrap replicate:
    1. Resample counts: Y* ~ Poisson(Y) in original gene space
    2. Preprocess and sketch using the same pipeline
    3. Solve BCD with warm start from converged solution
    4. Collect proportion estimates

    Parameters
    ----------
    Y : ndarray of shape (n_spots, n_genes)
        Original count matrix (before preprocessing).
    X : ndarray of shape (n_cell_types, n_genes)
        Reference signatures (before preprocessing).
    coords : ndarray of shape (n_spots, 2)
        Spatial coordinates.
    model_params : dict
        Parameters to reconstruct the FlashDeconv pipeline:
        'gene_idx', 'sketch_dim', 'preprocess', 'lambda_', 'rho',
        'leverage_scores', 'random_state', 'adjacency', 'sketch_matrix'.
    beta_init : ndarray of shape (n_spots, n_cell_types)
        Converged solution for warm start.
    n_bootstrap : int
        Number of bootstrap replicates.
    max_iter_boot : int
        Max BCD iterations per bootstrap (warm start converges fast).
    seed : int
        Random seed for reproducibility.
    verbose : bool

    Returns
    -------
    dict with keys:
        'boot_mean': ndarray (n_spots, n_cell_types)
        'boot_std': ndarray (n_spots, n_cell_types)
        'boot_ci_lower': ndarray (n_spots, n_cell_types) - 2.5th percentile
        'boot_ci_upper': ndarray (n_spots, n_cell_types) - 97.5th percentile
        'boot_cv': ndarray (n_spots, n_cell_types) - coefficient of variation
        'n_bootstrap': int
    """
    from .solver import bcd_solve, normalize_proportions
    from .sketching import sketch_data
    from .deconv import FlashDeconv

    rng = np.random.default_rng(seed)

    gene_idx = model_params['gene_idx']
    sketch_dim = model_params['sketch_dim']
    preprocess = model_params['preprocess']
    lambda_ = model_params['lambda_']
    rho = model_params['rho']
    leverage_scores = model_params['leverage_scores']
    random_state = model_params['random_state']
    A = model_params['adjacency']

    n_spots, n_types = beta_init.shape

    # Pre-subset reference (same for all bootstraps)
    X_subset = X[:, gene_idx]

    # We need a dummy FlashDeconv instance for preprocessing
    dummy = FlashDeconv(preprocess=preprocess)

    # Collect bootstrap proportions
    boot_props = np.zeros((n_bootstrap, n_spots, n_types), dtype=np.float32)

    for b in range(n_bootstrap):
        if verbose and (b % 10 == 0):
            print(f"  Bootstrap {b+1}/{n_bootstrap}...")

        # Poisson resample original counts
        Y_dense = Y.toarray() if sparse.issparse(Y) else np.asarray(Y)
        Y_boot = rng.poisson(np.maximum(Y_dense, 0).astype(np.float64))
        Y_boot = Y_boot.astype(np.float32)

        # Subset to same genes
        Y_boot_subset = Y_boot[:, gene_idx]

        # Preprocess
        Y_tilde, X_tilde = dummy._preprocess_data(
            Y_boot_subset, X_subset, preprocess
        )

        # Sketch with same random state (same hash functions)
        Y_sketch, X_sketch, _ = sketch_data(
            Y_tilde, X_tilde,
            sketch_dim=sketch_dim,
            leverage_scores=leverage_scores,
            random_state=random_state,
        )

        # Warm-start BCD
        beta_boot, _ = bcd_solve(
            Y_sketch, X_sketch, A,
            lambda_=lambda_,
            rho=rho,
            max_iter=max_iter_boot,
            tol=1e-3,  # Looser tolerance for speed
            verbose=False,
        )

        boot_props[b] = normalize_proportions(beta_boot).astype(np.float32)

    # Compute statistics
    boot_mean = boot_props.mean(axis=0)
    boot_std = boot_props.std(axis=0)
    boot_ci_lower = np.percentile(boot_props, 2.5, axis=0)
    boot_ci_upper = np.percentile(boot_props, 97.5, axis=0)

    boot_cv = np.zeros_like(boot_mean)
    nonzero = boot_mean > 1e-6
    boot_cv[nonzero] = boot_std[nonzero] / boot_mean[nonzero]

    return {
        'boot_mean': boot_mean,
        'boot_std': boot_std,
        'boot_ci_lower': boot_ci_lower,
        'boot_ci_upper': boot_ci_upper,
        'boot_cv': boot_cv,
        'n_bootstrap': n_bootstrap,
    }
