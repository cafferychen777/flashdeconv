"""
Block Coordinate Descent solver for FlashDeconv.

This module implements a Numba-accelerated BCD algorithm for solving
the spatial-regularized non-negative least squares problem.

Optimizations:
- Precompute Gram matrix G = X @ X^T (avoids O(K*K*d) per iteration)
- Precompute H = X @ Y^T (avoids O(T*N*K*d) recomputation across iterations)
"""

import numpy as np
from numba import jit, prange
from typing import Optional, Tuple, List
from scipy import sparse


@jit(nopython=True, cache=True)
def soft_threshold(x: float, threshold: float) -> float:
    """Soft thresholding operator for L1 regularization."""
    if x > threshold:
        return x - threshold
    elif x < -threshold:
        return x + threshold
    else:
        return 0.0


@jit(nopython=True, cache=True)
def update_spot_with_Xty(
    Xty_i: np.ndarray,
    XtX: np.ndarray,
    beta_i: np.ndarray,
    neighbor_sum: np.ndarray,
    n_neighbors: int,
    lambda_: float,
    rho: float,
) -> np.ndarray:
    """
    Update cell type abundances for a single spot using precomputed Xty.

    Uses maintained residual r = XtX @ beta_i to avoid recomputing the
    full dot product per coordinate.  When beta_i[k] does not change
    (common under non-negativity + L1), the O(K) residual update is skipped.

    Parameters
    ----------
    Xty_i : ndarray of shape (n_cell_types,)
        Precomputed X_sketch @ y_i for spot i.
    XtX : ndarray of shape (n_cell_types, n_cell_types)
        Precomputed X_sketch @ X_sketch.T (Gram matrix)
    beta_i : ndarray of shape (n_cell_types,)
        Current estimates (modified in place).
    neighbor_sum : ndarray of shape (n_cell_types,)
        Sum of beta over neighbors.
    n_neighbors : int
        Number of neighbors for this spot.
    lambda_ : float
        Spatial regularization strength.
    rho : float
        Sparsity regularization strength.

    Returns
    -------
    beta_i : ndarray of shape (n_cell_types,)
        Updated estimates.
    """
    n_cell_types = beta_i.shape[0]

    # Initialize maintained residual with a dense matrix-vector product.
    # This is equivalent to the nested loops below, but faster in numba.
    r = XtX @ beta_i

    # Coordinate descent over cell types
    for k in range(n_cell_types):
        old_k = beta_i[k]

        # Partial residual: exclude self-contribution from r
        residual_k = Xty_i[k] - r[k] + XtX[k, k] * old_k

        # Add spatial regularization term
        if n_neighbors > 0:
            residual_k += lambda_ * neighbor_sum[k]

        # Compute denominator
        denom = XtX[k, k] + lambda_ * n_neighbors

        # Apply soft thresholding for L1 and project to non-negative
        if denom > 1e-10:
            beta_k_new = soft_threshold(residual_k, rho) / denom
            beta_i[k] = max(0.0, beta_k_new)
        else:
            beta_i[k] = 0.0

        # Update maintained residual incrementally
        delta = beta_i[k] - old_k
        if delta != 0.0:
            for kk in range(n_cell_types):
                r[kk] += delta * XtX[kk, k]

    return beta_i


@jit(nopython=True, parallel=True, cache=True)
def _bcd_iteration_fused(
    H: np.ndarray,
    XtX: np.ndarray,
    beta_in: np.ndarray,
    beta_out: np.ndarray,
    neighbor_indices: np.ndarray,
    neighbor_indptr: np.ndarray,
    lambda_: float,
    rho: float,
    spot_diffs: np.ndarray,
    spot_abs: np.ndarray,
) -> None:
    """
    Single BCD iteration with double-buffered I/O and fused convergence stats.

    Reads from beta_in, writes to beta_out, and computes per-spot max absolute
    change (spot_diffs) and max absolute old value (spot_abs) in the same pass.

    Parameters
    ----------
    H : ndarray of shape (n_cell_types, n_spots)
        Precomputed X_sketch @ Y_sketch.T matrix.
    XtX : ndarray of shape (n_cell_types, n_cell_types)
        Precomputed Gram matrix.
    beta_in : ndarray of shape (n_spots, n_cell_types)
        Current (read-only) cell type abundances.
    beta_out : ndarray of shape (n_spots, n_cell_types)
        Output buffer for updated abundances.
    neighbor_indices : ndarray
        Flattened neighbor indices (CSR format).
    neighbor_indptr : ndarray
        Index pointers for neighbors (CSR format).
    lambda_ : float
        Spatial regularization.
    rho : float
        Sparsity regularization.
    spot_diffs : ndarray of shape (n_spots,)
        Output: per-spot max |beta_new - beta_old|.
    spot_abs : ndarray of shape (n_spots,)
        Output: per-spot max |beta_old|.
    """
    n_spots = beta_in.shape[0]
    n_cell_types = beta_in.shape[1]

    for i in prange(n_spots):
        Xty_i = H[:, i]

        # Copy old values into output buffer (workspace for coordinate descent)
        for k in range(n_cell_types):
            beta_out[i, k] = beta_in[i, k]

        # Get neighbors
        start = neighbor_indptr[i]
        end = neighbor_indptr[i + 1]
        n_neighbors = end - start

        # Compute neighbor sum (reads from beta_in: Jacobi pattern)
        neighbor_sum = np.zeros(n_cell_types)
        for idx in range(start, end):
            j = neighbor_indices[idx]
            for k in range(n_cell_types):
                neighbor_sum[k] += beta_in[j, k]

        # Update this spot (modifies beta_out[i] in place via CD)
        update_spot_with_Xty(
            Xty_i, XtX, beta_out[i], neighbor_sum, n_neighbors, lambda_, rho
        )

        # Fused convergence statistics
        local_diff = 0.0
        local_abs = 0.0
        for k in range(n_cell_types):
            d = abs(beta_out[i, k] - beta_in[i, k])
            if d > local_diff:
                local_diff = d
            a = abs(beta_in[i, k])
            if a > local_abs:
                local_abs = a
        spot_diffs[i] = local_diff
        spot_abs[i] = local_abs


def precompute_gram_matrix(X_sketch: np.ndarray) -> np.ndarray:
    """
    Precompute the Gram matrix X_sketch @ X_sketch.T.

    Parameters
    ----------
    X_sketch : ndarray of shape (n_cell_types, sketch_dim)
        Sketched reference.

    Returns
    -------
    XtX : ndarray of shape (n_cell_types, n_cell_types)
        Gram matrix.
    """
    return X_sketch @ X_sketch.T


def precompute_XtY(X_sketch: np.ndarray, Y_sketch: np.ndarray) -> np.ndarray:
    """
    Precompute H = X_sketch @ Y_sketch.T.

    This avoids redundant computation of X @ y_i across iterations.
    Complexity: O(N * K * d) once, instead of O(T * N * K * d) total.

    Parameters
    ----------
    X_sketch : ndarray of shape (n_cell_types, sketch_dim)
        Sketched reference signatures.
    Y_sketch : ndarray of shape (n_spots, sketch_dim)
        Sketched spatial data.

    Returns
    -------
    H : ndarray of shape (n_cell_types, n_spots)
        Precomputed X_sketch @ Y_sketch.T matrix.
    """
    return X_sketch @ Y_sketch.T


def compute_objective(
    beta: np.ndarray,
    H: np.ndarray,
    XtX: np.ndarray,
    YtY: float,
    L: sparse.spmatrix,
    lambda_: float,
    rho: float,
) -> float:
    """
    Compute objective using precomputed matrices (no N*d temporaries).

    Uses the algebraic expansion:
        0.5*||Y - bX||^2 = 0.5*(||Y||^2 - 2*Tr(Y^T b X) + Tr(b^T b X X^T))

    Parameters
    ----------
    beta : ndarray of shape (n_spots, n_cell_types)
        Current solution.
    H : ndarray of shape (n_cell_types, n_spots)
        Precomputed X @ Y^T.
    XtX : ndarray of shape (n_cell_types, n_cell_types)
        Precomputed X @ X^T (Gram matrix).
    YtY : float
        Precomputed ||Y||_F^2.
    L : sparse matrix
        Graph Laplacian.
    lambda_, rho : float
        Regularization parameters.

    Returns
    -------
    obj : float
        Objective value.
    """
    # Fidelity: 0.5*(YtY - 2*cross + quad) via precomputed matrices
    # cross = Tr(Y^T beta X) = sum_{ik} beta_{ik} H_{ki}
    cross = np.sum(beta * H.T)
    # quad = ||beta X||_F^2 = Tr(beta^T beta @ XtX)
    BtB = beta.T @ beta  # (K, K)
    quad = np.sum(BtB * XtX)
    fidelity = 0.5 * (YtY - 2.0 * cross + quad)

    # Spatial smoothing term
    Lbeta = L @ beta
    spatial = lambda_ * np.sum(beta * Lbeta)

    # Sparsity term
    sparsity = rho * np.sum(np.abs(beta))

    return fidelity + spatial + sparsity


def bcd_solve(
    Y_sketch: np.ndarray,
    X_sketch: np.ndarray,
    A: sparse.spmatrix,
    lambda_: float = 0.1,
    rho: float = 0.01,
    max_iter: int = 100,
    tol: float = 1e-4,
    verbose: bool = False,
) -> Tuple[np.ndarray, dict]:
    """
    Solve the spatial-regularized deconvolution problem via BCD.

    Minimize:
    0.5 * ||Y_sketch - beta @ X_sketch||_F^2 + lambda * Tr(beta^T L beta) + rho * ||beta||_1
    subject to: beta >= 0

    Parameters
    ----------
    Y_sketch : ndarray of shape (n_spots, sketch_dim)
        Sketched spatial data.
    X_sketch : ndarray of shape (n_cell_types, sketch_dim)
        Sketched reference signatures.
    A : sparse matrix of shape (n_spots, n_spots)
        Spatial adjacency matrix.
    lambda_ : float, default=0.1
        Spatial regularization strength.
    rho : float, default=0.01
        Sparsity regularization strength (L1).
    max_iter : int, default=100
        Maximum number of iterations.
    tol : float, default=1e-4
        Convergence tolerance (relative change in beta).
    verbose : bool, default=False
        Whether to print progress.

    Returns
    -------
    beta : ndarray of shape (n_spots, n_cell_types)
        Estimated cell type abundances.
    info : dict
        Optimization information including convergence status.
    """
    n_spots, sketch_dim = Y_sketch.shape
    n_cell_types = X_sketch.shape[0]

    # Handle empty input
    if n_spots == 0 or n_cell_types == 0:
        beta = np.empty((n_spots, n_cell_types), dtype=np.float64)
        info = {
            'converged': True,
            'n_iterations': 0,
            'final_objective': 0.0,
            'objectives': [],
            'final_change': 0.0,
        }
        return beta, info

    # Precompute matrices for efficiency
    XtX = precompute_gram_matrix(X_sketch)   # (K, K)
    H = precompute_XtY(X_sketch, Y_sketch)   # (K, N)
    YtY = float(np.sum(Y_sketch ** 2))       # scalar, for fast objective

    # Scale rho so that the user-facing parameter is a dimensionless fraction.
    # The partial residual r_ik ~ O(diag(G)); without scaling, rho ~ 0.01
    # is negligible vs r_ik ~ 7000 and soft-thresholding has no effect.
    #
    # Note: the non-negativity constraint beta >= 0 already provides the
    # primary sparsity (typically zeroing 50-77% of entries).  L1 is a
    # secondary refinement that eliminates weakly-positive coefficients
    # in the range (0, rho * G_bar].  Its marginal benefit increases with
    # K / (true types per spot), i.e., fine-grained cell type taxonomies.
    gram_diag_mean = np.mean(np.diag(XtX))
    rho = rho * gram_diag_mean

    # Convert adjacency to CSR for efficient neighbor access
    A_csr = A.tocsr()
    neighbor_indices = A_csr.indices.astype(np.int64)
    neighbor_indptr = A_csr.indptr.astype(np.int64)

    # Compute Laplacian for objective
    from flashdeconv.core.spatial import compute_laplacian
    L = compute_laplacian(A)

    # Double-buffered iteration: swap roles each iteration, zero allocation
    beta_a = np.ones((n_spots, n_cell_types), dtype=np.float64) / n_cell_types
    beta_b = np.empty_like(beta_a)

    # Pre-allocate convergence stat buffers (N floats each, reused every iter)
    spot_diffs = np.empty(n_spots, dtype=np.float64)
    spot_abs = np.empty(n_spots, dtype=np.float64)

    # Track convergence
    objectives = []
    converged = False
    iteration = -1
    rel_change = 0.0

    for iteration in range(max_iter):
        # BCD iteration: read beta_a, write beta_b, fuse convergence stats
        _bcd_iteration_fused(
            H, XtX, beta_a, beta_b,
            neighbor_indices, neighbor_indptr,
            lambda_, rho,
            spot_diffs, spot_abs,
        )

        # Convergence from fused stats (no extra full-matrix scan)
        max_diff = spot_diffs.max()
        max_abs_old = spot_abs.max()
        rel_change = max_diff / (max_abs_old + 1e-10)

        if verbose and (iteration % 10 == 0 or iteration == max_iter - 1):
            obj = compute_objective(
                beta_b, H, XtX, YtY, L, lambda_, rho
            )
            objectives.append(obj)
            print(f"Iteration {iteration}: objective = {obj:.6f}, rel_change = {rel_change:.6e}")

        # Swap buffers (pointer swap, zero copy)
        beta_a, beta_b = beta_b, beta_a

        if rel_change < tol:
            converged = True
            if verbose:
                print(f"Converged at iteration {iteration}")
            break

    # After swap, beta_a always holds the latest result
    final_obj = compute_objective(
        beta_a, H, XtX, YtY, L, lambda_, rho
    )

    info = {
        'converged': converged,
        'n_iterations': iteration + 1,
        'final_objective': final_obj,
        'objectives': objectives if verbose else [],
        'final_change': rel_change,
    }

    return beta_a, info


def normalize_proportions(beta: np.ndarray) -> np.ndarray:
    """
    Normalize beta to sum to 1 per spot (cell type proportions).

    Parameters
    ----------
    beta : ndarray of shape (n_spots, n_cell_types)
        Raw abundances.

    Returns
    -------
    proportions : ndarray of shape (n_spots, n_cell_types)
        Normalized proportions (sum to 1 per row).
    """
    row_sums = np.sum(beta, axis=1, keepdims=True)
    zero_mask = (row_sums == 0).ravel()
    row_sums = np.maximum(row_sums, 1e-10)
    proportions = beta / row_sums
    # All-zero rows get uniform distribution (maximum entropy)
    if np.any(zero_mask):
        proportions[zero_mask] = 1.0 / beta.shape[1]
    return proportions
