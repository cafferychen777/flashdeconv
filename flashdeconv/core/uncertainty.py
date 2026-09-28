"""
Uncertainty quantification for FlashDeconv.

All quantities in this module are *model-based*: they describe sampling
variability of the estimate under the fitted model (the given reference,
the log-space linear mixture, the selected genes and the fitted penalties).
They do not include error from reference-tissue mismatch, cell types missing
from the reference, or the log-space mixture approximation. In benchmarks
against true cell-type fractions these sources dominate the error, so the
intervals below are not calibrated confidence intervals for true
proportions. They are useful for comparing the stability of estimates across
spots and cell types.

Tier 0 (cheap diagnostics):
    - Prediction entropy of each proportion vector.
    - Reconstruction residual in the solver's gene representation.

Tier 1 (analytical, opt-in):
    - Per-spot active-set sandwich (HC0 robust) variance with the full
      inverse Hessian of the penalized objective, delta method to
      proportions, and boundary-aware one-sided upper bounds for cell types
      estimated at zero (``compute_active_set_uncertainty``).
    - One-sided Wald z-statistics for "abundance > 0" and Benjamini-Hochberg
      detection calls (``bh_detection``).
    - ``compute_hessian_variance`` (the pre-0.2.0 diagonal-Hessian Laplace
      approximation) is deprecated: it ignores collinearity between cell
      types and is anti-conservative even when the model holds.

Tier 2 (bootstrap, opt-in):
    - Poisson bootstrap of the observed counts with refits of the full
      estimator warm-started at the fitted solution
      (``poisson_bootstrap``).
"""

import warnings
from typing import Dict, Optional, Tuple

import numpy as np
from numba import njit, prange
from scipy import sparse
from scipy.stats import norm


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
    Per-spot reconstruction error in the solver's gene representation.

    Parameters
    ----------
    Y_sketch : ndarray or sparse matrix of shape (n_spots, p)
        Spatial data in the gene representation used by the solver.
    X_sketch : ndarray of shape (n_cell_types, p)
    beta : ndarray of shape (n_spots, n_cell_types)

    Returns
    -------
    residual_ss : ndarray of shape (n_spots,)
        Sum of squared residuals per spot.
    residual_norm : ndarray of shape (n_spots,)
        Normalized residual: ||e_i||^2 / ||y_i||^2.
    """
    if sparse.issparse(Y_sketch):
        # Row blocks keep the dense temporaries bounded.
        Y_csr = Y_sketch.tocsr()
        n_spots, p = Y_csr.shape
        residual_ss = np.empty(n_spots, dtype=np.float64)
        y_ss = np.empty(n_spots, dtype=np.float64)
        block = max(1, (1 << 22) // max(p, 1))
        for s in range(0, n_spots, block):
            e = min(s + block, n_spots)
            Y_blk = Y_csr[s:e].toarray()
            res_blk = Y_blk - beta[s:e] @ X_sketch
            residual_ss[s:e] = np.sum(res_blk ** 2, axis=1)
            y_ss[s:e] = np.sum(Y_blk ** 2, axis=1)
        y_ss = np.maximum(y_ss, 1e-10)
        return residual_ss, residual_ss / y_ss

    # Predicted: beta @ X_sketch -> (n_spots, p)
    Y_hat = beta @ X_sketch
    residual = Y_sketch - Y_hat
    residual_ss = np.sum(residual ** 2, axis=1)

    y_ss = np.sum(Y_sketch ** 2, axis=1)
    y_ss = np.maximum(y_ss, 1e-10)
    residual_norm = residual_ss / y_ss

    return residual_ss, residual_norm


# ---------------------------------------------------------------------------
# Tier 1: active-set sandwich variance
# ---------------------------------------------------------------------------

_Z_CAP = 50.0  # z-statistic assigned when a standard error is exactly zero


@njit(cache=True, fastmath=False)
def _spot_active_set(i, y, Xs, Xs2, G, b, deg_i, nsum, lam, rho, z, sandwich,
                     lo, hi, se_p, zst):
    """Active-set intervals for one spot (writes row i of the outputs).

    Objective (per spot, neighbours held fixed):
        0.5 ||y - b Xs||^2 + 0.5 lam sum_j ||b - b_j||^2 + rho ||b||_1, b >= 0,
    whose gradient is -Xs e + lam (deg_i b - nsum) + rho with
    nsum = sum_j b_j (the neighbours' abundances).
    Hessian on the active set A: H_A = G_AA + lam * deg_i * I.
    Variance of b_A: H_A^-1 M_AA H_A^-1 with M = Xs diag(e^2) Xs^T (HC0),
    or s^2 H_A^-1 (model-based). Inactive k: one Newton step for b_k from 0
    on the augmented set A + {k}, Wald upper bound truncated at b_k >= 0.
    """
    K, p = Xs.shape
    act = np.where(b > 0.0)[0]
    na = act.shape[0]

    # Residual e = y - b_A Xs_A (only active rows contribute)
    e = y.copy()
    for u in range(na):
        a = act[u]
        ba = b[a]
        for g in range(p):
            e[g] -= ba * Xs[a, g]
    ee = 0.0
    for g in range(p):
        ee += e[g] * e[g]
    w = e * e

    S = 0.0
    for u in range(na):
        S += b[act[u]]
    if S <= 0.0:
        # No fitted signal: the spot carries no information on composition.
        for k in range(K):
            lo[i, k] = 0.0
            hi[i, k] = 1.0
            se_p[i, k] = np.nan
            zst[i, k] = 0.0
        return

    # Gradient of the data term g = Xs e (needed for inactive types only)
    grad = np.zeros(K)
    # M rows for active types against all types, and diag(M)
    MA = np.zeros((na, K))
    Mdiag = np.zeros(K)
    for k in range(K):
        if b[k] > 0.0:
            continue
        acc = 0.0
        for g in range(p):
            acc += Xs[k, g] * e[g]
        grad[k] = acc
    if sandwich:
        for k in range(K):
            acc = 0.0
            for g in range(p):
                acc += Xs2[k, g] * w[g]
            Mdiag[k] = acc
        for u in range(na):
            a = act[u]
            for k in range(K):
                acc = 0.0
                for g in range(p):
                    acc += Xs[a, g] * Xs[k, g] * w[g]
                MA[u, k] = acc

    pen = lam * deg_i
    HA = np.empty((na, na))
    for u in range(na):
        for v in range(na):
            HA[u, v] = G[act[u], act[v]]
        HA[u, u] += pen
    HAi = np.linalg.pinv(HA)

    if sandwich:
        MAA = np.empty((na, na))
        for u in range(na):
            for v in range(na):
                MAA[u, v] = MA[u, act[v]]
        VA = HAi @ MAA @ HAi
    else:
        VA = (ee / max(p - na, 1)) * HAi

    # Active types: delta method to proportions p_u = b_u / S
    rowsum = np.zeros(na)
    tot = 0.0
    for u in range(na):
        for v in range(na):
            rowsum[u] += VA[u, v]
        tot += rowsum[u]
    for u in range(na):
        k = act[u]
        pk = b[k] / S
        var = (VA[u, u] - 2.0 * pk * rowsum[u] + pk * pk * tot) / (S * S)
        sd = np.sqrt(max(var, 0.0))
        se_p[i, k] = sd
        lo[i, k] = min(max(pk - z * sd, 0.0), 1.0)
        hi[i, k] = min(max(pk + z * sd, 0.0), 1.0)
        sb = np.sqrt(max(VA[u, u], 0.0))
        zst[i, k] = b[k] / sb if sb > 0.0 else _Z_CAP

    # Inactive types: Schur complement of the augmented Hessian
    h = np.empty(na)
    for k in range(K):
        if b[k] > 0.0:
            continue
        for u in range(na):
            h[u] = G[act[u], k]
        v = HAi @ h
        s = G[k, k] + pen
        for u in range(na):
            s -= h[u] * v[u]
        lo[i, k] = 0.0
        if s <= 1e-12 * (G[k, k] + pen + 1e-300):
            # Type k is (numerically) in the span of the active set.
            hi[i, k] = 1.0
            se_p[i, k] = np.nan
            zst[i, k] = 0.0
            continue
        # KKT: the penalized gradient on the active coordinates is zero.
        # At b_k = 0 the gradient in k is -(Xs e)_k - lam * nsum_k + rho, so
        # the Newton step on A + {k} gives b_k = (grad_k + lam nsum_k - rho) / s.
        bk = (grad[k] + lam * nsum[k] - rho) / s
        if sandwich:
            # c = [-v/s, 1/s];  Var = c^T M_BB c
            q = Mdiag[k]
            for u in range(na):
                q -= 2.0 * v[u] * MA[u, k]
                for t in range(na):
                    q += v[u] * v[t] * MA[u, act[t]]
            vk = q / (s * s)
        else:
            vk = (ee / max(p - na - 1, 1)) / s
        sk = np.sqrt(max(vk, 0.0))
        up = max(bk + z * sk, 0.0)
        hi[i, k] = min(up / (S + up), 1.0)
        se_p[i, k] = sk / S
        if sk > 0.0:
            zst[i, k] = bk / sk
        else:
            zst[i, k] = -_Z_CAP if bk <= 0.0 else _Z_CAP


@njit(parallel=True, cache=True)
def _active_set_csr(indptr, indices, data, Xs, Xs2, G, beta, deg, nsum, lam,
                    rho, z, sandwich):
    n, K = beta.shape
    p = Xs.shape[1]
    lo = np.zeros((n, K))
    hi = np.zeros((n, K))
    se_p = np.zeros((n, K))
    zst = np.zeros((n, K))
    for i in prange(n):
        y = np.zeros(p)
        for q in range(indptr[i], indptr[i + 1]):
            y[indices[q]] = data[q]
        _spot_active_set(i, y, Xs, Xs2, G, beta[i].copy(), deg[i], nsum[i],
                         lam, rho, z, sandwich, lo, hi, se_p, zst)
    return lo, hi, se_p, zst


@njit(parallel=True, cache=True)
def _active_set_dense(Y, Xs, Xs2, G, beta, deg, nsum, lam, rho, z, sandwich):
    n, K = beta.shape
    lo = np.zeros((n, K))
    hi = np.zeros((n, K))
    se_p = np.zeros((n, K))
    zst = np.zeros((n, K))
    for i in prange(n):
        _spot_active_set(i, Y[i].copy(), Xs, Xs2, G, beta[i].copy(), deg[i],
                         nsum[i], lam, rho, z, sandwich, lo, hi, se_p, zst)
    return lo, hi, se_p, zst


def compute_active_set_uncertainty(
    Y_sketch,
    X_sketch: np.ndarray,
    beta: np.ndarray,
    lambda_: float,
    rho_abs: float,
    n_neighbors: np.ndarray,
    alpha: float = 0.05,
    variance: str = "sandwich",
    neighbor_sum: Optional[np.ndarray] = None,
) -> Dict[str, np.ndarray]:
    """
    Per-spot active-set variance with boundary-aware intervals.

    For spot ``i`` with active set ``A = {k : beta_ik > 0}`` and neighbours
    held fixed, the Hessian of the penalized objective on ``A`` is
    ``H_A = G_AA + lambda * n_neighbors_i * I`` (``G = X_sketch X_sketch^T``).
    The abundance covariance is the sandwich ``H_A^-1 M_AA H_A^-1`` with
    ``M = X_sketch diag(e_i^2) X_sketch^T`` (HC0, robust to heteroscedastic
    residuals) or, with ``variance="model"``, ``s_i^2 H_A^-1``. The full
    inverse is used, so collinearity between cell types widens the
    intervals. Proportion variances follow from the delta method for
    ``p = beta / sum(beta)``.

    For a type estimated at zero, a single Newton step from 0 on the
    augmented set ``A + {k}`` (using the KKT condition on ``A``) gives
    ``b_k = (x_k^T e_i + lambda * sum_j beta_jk - rho) / s`` (``s`` the Schur
    complement of ``H_A`` in the augmented Hessian; the neighbour sum is the
    gradient of the spatial penalty ``0.5 lambda Tr(beta^T L beta)`` at
    ``beta_ik = 0``) and a standard error ``s_k``; the interval is
    ``[0, u / (S + u)]`` with ``u = max(b_k + z s_k, 0)`` and
    ``S = sum(beta_i)``, i.e. a Wald bound truncated to the parameter space.

    All quantities are conditional on the fitted model (reference, genes,
    penalties); see the module docstring.

    Parameters
    ----------
    Y_sketch : ndarray or sparse matrix of shape (n_spots, p)
        Spatial data in the solver's representation.
    X_sketch : ndarray of shape (n_cell_types, p)
    beta : ndarray of shape (n_spots, n_cell_types)
        Fitted abundances.
    lambda_ : float
        Spatial regularization used in the fit.
    rho_abs : float
        L1 weight in the solver's scale (``rho_sparsity * mean(diag(G))``).
    n_neighbors : ndarray of shape (n_spots,)
    alpha : float, default=0.05
        Two-sided level for active types; the same ``z_{1-alpha/2}`` is used
        for the upper bound of inactive types.
    variance : {"sandwich", "model"}, default="sandwich"
    neighbor_sum : ndarray of shape (n_spots, n_cell_types), optional
        ``A @ beta``: per-spot sum of the neighbours' abundances. Required for
        correct upper bounds of zero-estimated types when ``lambda_ > 0``;
        if omitted it is taken as zero (with a warning when it matters).

    Returns
    -------
    dict with ``ci_lower``, ``ci_upper``, ``se_prop`` (standard error on the
    proportion scale; NaN where undefined) and ``z_score`` (one-sided Wald
    statistic for ``beta_ik > 0``), each of shape (n_spots, n_cell_types).
    """
    if variance not in ("sandwich", "model"):
        raise ValueError(f"variance must be 'sandwich' or 'model', got {variance!r}")
    z = float(norm.ppf(1 - alpha / 2))
    Xs = np.ascontiguousarray(X_sketch, dtype=np.float64)
    Xs2 = Xs * Xs
    G = Xs @ Xs.T
    B = np.ascontiguousarray(beta, dtype=np.float64)
    deg = np.ascontiguousarray(n_neighbors, dtype=np.float64)
    if neighbor_sum is None:
        if lambda_ > 0 and np.any(deg > 0):
            warnings.warn(
                "neighbor_sum not given: upper bounds of zero-estimated types "
                "ignore the spatial penalty; pass neighbor_sum=A @ beta.",
                RuntimeWarning,
                stacklevel=2,
            )
        nsum = np.zeros_like(B)
    else:
        nsum = np.ascontiguousarray(neighbor_sum, dtype=np.float64)
        if nsum.shape != B.shape:
            raise ValueError(
                f"neighbor_sum must have shape {B.shape}, got {nsum.shape}"
            )
    sandwich = variance == "sandwich"
    if sparse.issparse(Y_sketch):
        Yc = Y_sketch.tocsr()
        lo, hi, se_p, zst = _active_set_csr(
            Yc.indptr.astype(np.int64), Yc.indices.astype(np.int64),
            Yc.data.astype(np.float64), Xs, Xs2, G, B, deg, nsum,
            float(lambda_), float(rho_abs), z, sandwich,
        )
    else:
        Yd = np.ascontiguousarray(Y_sketch, dtype=np.float64)
        lo, hi, se_p, zst = _active_set_dense(
            Yd, Xs, Xs2, G, B, deg, nsum, float(lambda_), float(rho_abs), z,
            sandwich,
        )
    return {"ci_lower": lo, "ci_upper": hi, "se_prop": se_p, "z_score": zst}


def bh_detection(pvalues: np.ndarray, q: float = 0.1) -> np.ndarray:
    """
    Benjamini-Hochberg step-up procedure over all entries.

    Parameters
    ----------
    pvalues : ndarray (any shape)
        One-sided p-values; NaN entries are never rejected.
    q : float, default=0.1
        Target false discovery rate.

    Returns
    -------
    reject : bool ndarray of the same shape as ``pvalues``.
    """
    if not 0 < q < 1:
        raise ValueError(f"q must be in (0, 1), got {q}")
    pv = np.asarray(pvalues, dtype=np.float64)
    flat = pv.ravel()
    ok = np.isfinite(flat)
    m = int(ok.sum())
    reject = np.zeros(flat.shape, dtype=bool)
    if m == 0:
        return reject.reshape(pv.shape)
    idx = np.where(ok)[0]
    order = idx[np.argsort(flat[idx], kind="stable")]
    thr = q * np.arange(1, m + 1) / m
    passed = np.nonzero(flat[order] <= thr)[0]
    if passed.size:
        reject[order[: passed[-1] + 1]] = True
    return reject.reshape(pv.shape)


# ---------------------------------------------------------------------------
# Deprecated diagonal-Hessian Laplace approximation (pre-0.2.0 Tier 1)
# ---------------------------------------------------------------------------

def compute_hessian_variance(
    beta: np.ndarray,
    XtX: np.ndarray,
    residual_ss: np.ndarray,
    sketch_dim: int,
    lambda_: float,
    n_neighbors: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Deprecated diagonal-Hessian Laplace approximation.

    Uses ``sigma2_i / (XtX[k,k] + lambda * n_neighbors_i)``, i.e. ``1/diag(H)``
    instead of ``diag(H^-1)``. This ignores collinearity between cell types
    and under-covers even when the data follow the model; types estimated at
    zero get variance 0. Kept only for ``method="laplace_diag"``; use
    :func:`compute_active_set_uncertainty` instead.

    Parameters
    ----------
    beta : ndarray of shape (n_spots, n_cell_types)
    XtX : ndarray of shape (n_cell_types, n_cell_types)
    residual_ss : ndarray of shape (n_spots,)
    sketch_dim : int
        Dimension p of the solver's representation.
    lambda_ : float
    n_neighbors : ndarray of shape (n_spots,)

    Returns
    -------
    var_beta, var_prop : ndarray of shape (n_spots, n_cell_types)
    """
    n_spots, n_types = beta.shape
    active = beta > 0
    n_active = active.sum(axis=1)

    dof = np.maximum(sketch_dim - n_active, 1.0)
    sigma2 = residual_ss / dof

    xtx_diag = np.diag(XtX)
    h_diag = xtx_diag[np.newaxis, :] + lambda_ * n_neighbors[:, np.newaxis]

    var_beta = np.zeros((n_spots, n_types), dtype=np.float64)
    var_beta[active] = (sigma2[:, np.newaxis] / h_diag)[active]

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
    Symmetric Wald intervals from a proportion variance (used by the
    deprecated ``method="laplace_diag"``).

    Returns
    -------
    dict with 'ci_half_width', 'ci_lower', 'ci_upper' (clipped to [0, 1]),
    'cv' and 'detection_confident' (lower bound > 0.01; realised false
    discovery rates of 0.19-0.38 were observed in benchmarks, prefer the
    BH-based ``detected`` of the default method).
    """
    z = norm.ppf(1 - alpha / 2)

    std_prop = np.sqrt(np.maximum(var_prop, 0))
    ci_half = z * std_prop

    ci_lower = np.clip(proportions - ci_half, 0, 1)
    ci_upper = np.clip(proportions + ci_half, 0, 1)

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


# ---------------------------------------------------------------------------
# Tier 2: Poisson bootstrap
# ---------------------------------------------------------------------------

def poisson_bootstrap(
    Y,
    X: np.ndarray,
    coords: Optional[np.ndarray],
    model_params: Dict,
    beta_init: np.ndarray,
    n_bootstrap: int = 100,
    max_iter_boot: int = 1000,
    tol: float = 1e-4,
    seed: int = 42,
    verbose: bool = False,
    alpha: float = 0.05,
) -> Dict[str, np.ndarray]:
    """
    Poisson bootstrap of the observed counts with warm-started refits.

    For each replicate:
    1. Resample counts on the selected genes: Y* ~ Poisson(Y).
    2. Re-apply the same preprocessing and the fitted gene representation
       (fixed gene weights, or the same CountSketch matrix).
    3. Re-solve with the fitted lambda and rho, starting from the fitted
       ``beta_init``, with the same iteration limit and tolerance as the
       main fit.

    Genes, gene weights, lambda and rho are held at their fitted values,
    and the reference is fixed, so the spread reflects count sampling
    under the fitted model only (see the module docstring).

    Parameters
    ----------
    Y : ndarray or sparse matrix of shape (n_spots, n_genes)
        Original count matrix (before preprocessing).
    X : ndarray of shape (n_cell_types, n_genes)
        Reference signatures (before preprocessing).
    coords : ignored
        Kept for backward compatibility (the adjacency is passed in
        ``model_params``).
    model_params : dict
        'gene_idx', 'sketch_dim', 'preprocess', 'lambda_', 'rho',
        'leverage_scores', 'random_state', 'adjacency', and optionally
        'gene_weighting' ("expected" or "countsketch"; default
        "countsketch"), 'gene_weights' (required for "expected") and
        'sketch_matrix' (the fitted CountSketch matrix; if absent it is
        rebuilt from 'random_state', which reproduces the fit only for an
        integer seed).
    beta_init : ndarray of shape (n_spots, n_cell_types)
        Fitted abundances; every refit starts here.
    n_bootstrap : int
    max_iter_boot : int, default=1000
        Iteration limit per refit (use the main fit's ``max_iter``).
    tol : float, default=1e-4
        Convergence tolerance per refit (use the main fit's ``tol``).
    seed : int
    verbose : bool
    alpha : float, default=0.05
        Percentile interval level.

    Returns
    -------
    dict with 'boot_mean', 'boot_std', 'boot_ci_lower', 'boot_ci_upper'
    (percentile interval), 'boot_cv', 'n_bootstrap', 'n_converged' (refits
    that met ``tol``), 'n_iterations' and 'final_change' (relative change
    at the last iteration, per refit) and 'warm_start' (True).
    """
    from .preprocessing import preprocess_expression
    from .sketching import apply_gene_weights, project_to_sketch, sketch_data
    from .solver import bcd_solve, normalize_proportions

    if n_bootstrap < 1:
        raise ValueError(f"n_bootstrap must be at least 1, got {n_bootstrap}")
    if not 0.0 < alpha < 1.0:
        raise ValueError(f"alpha must be in (0, 1), got {alpha}")

    rng = np.random.default_rng(seed)

    gene_idx = model_params['gene_idx']
    sketch_dim = model_params['sketch_dim']
    preprocess = model_params['preprocess']
    lambda_ = model_params['lambda_']
    rho = model_params['rho']
    leverage_scores = model_params['leverage_scores']
    random_state = model_params['random_state']
    A = model_params['adjacency']
    gene_weighting = model_params.get('gene_weighting', 'countsketch')
    gene_weights = model_params.get('gene_weights')
    sketch_matrix = model_params.get('sketch_matrix')
    if gene_weighting == 'expected' and gene_weights is None:
        raise ValueError("model_params['gene_weights'] is required for "
                         "gene_weighting='expected'")

    n_spots, n_types = beta_init.shape

    # Preprocessing uses only the selected genes, so resampling them alone
    # is equivalent to resampling the full matrix.
    X_subset = X[:, gene_idx]
    if sparse.issparse(Y):
        Y_sel = sparse.csr_matrix(Y[:, gene_idx], dtype=np.float64)
        Y_sel.data = np.maximum(Y_sel.data, 0)
    else:
        Y_sel = np.maximum(np.asarray(Y)[:, gene_idx], 0).astype(np.float64)

    boot_props = np.zeros((n_bootstrap, n_spots, n_types), dtype=np.float32)
    n_iterations = np.zeros(n_bootstrap, dtype=np.int64)
    converged = np.zeros(n_bootstrap, dtype=bool)
    final_change = np.zeros(n_bootstrap, dtype=np.float64)

    for b in range(n_bootstrap):
        if verbose and (b % 10 == 0):
            print(f"  Bootstrap {b+1}/{n_bootstrap}...")

        if sparse.issparse(Y_sel):
            Y_boot = Y_sel.copy()
            Y_boot.data = rng.poisson(Y_boot.data).astype(np.float64)
            Y_boot.eliminate_zeros()
        else:
            Y_boot = rng.poisson(Y_sel).astype(np.float64)

        Y_tilde, X_tilde = preprocess_expression(Y_boot, X_subset, preprocess)
        if sparse.issparse(Y_tilde):
            Y_tilde = sparse.csr_matrix(Y_tilde)

        if gene_weighting == 'expected':
            Y_sketch, X_sketch = apply_gene_weights(Y_tilde, X_tilde, gene_weights)
        elif sketch_matrix is not None:
            Y_sketch, X_sketch = project_to_sketch(Y_tilde, X_tilde, sketch_matrix)
        else:
            Y_sketch, X_sketch, _ = sketch_data(
                Y_tilde, X_tilde,
                sketch_dim=sketch_dim,
                leverage_scores=leverage_scores,
                random_state=random_state,
            )

        beta_boot, info = bcd_solve(
            Y_sketch, X_sketch, A,
            lambda_=lambda_,
            rho=rho,
            max_iter=max_iter_boot,
            tol=tol,
            verbose=False,
            beta_init=beta_init,
        )
        n_iterations[b] = info['n_iterations']
        converged[b] = info['converged']
        final_change[b] = info['final_change']

        boot_props[b] = normalize_proportions(beta_boot).astype(np.float32)

    n_conv = int(converged.sum())
    if n_conv < n_bootstrap:
        warnings.warn(
            f"{n_bootstrap - n_conv} of {n_bootstrap} bootstrap refits did not "
            f"reach tol={tol} within {max_iter_boot} iterations; increase "
            f"max_iter_boot (and check that the main fit converged).",
            RuntimeWarning,
            stacklevel=3,
        )

    boot_mean = boot_props.mean(axis=0)
    boot_std = boot_props.std(axis=0)
    boot_ci_lower = np.percentile(boot_props, 100 * alpha / 2, axis=0)
    boot_ci_upper = np.percentile(boot_props, 100 * (1 - alpha / 2), axis=0)

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
        'n_converged': n_conv,
        'n_iterations': n_iterations,
        'final_change': final_change,
        'warm_start': True,
    }
