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

    This is the optimized version that takes precomputed X @ y_i directly,
    avoiding redundant computation across iterations.

    Solves:
    beta_i = argmin 0.5 * ||y_i - X_sketch^T beta_i||^2
             + lambda * ||beta_i - mean(beta_neighbors)||^2
             + rho * ||beta_i||_1
    subject to beta_i >= 0

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

    # Coordinate descent over cell types
    for k in range(n_cell_types):
        # Compute gradient contribution from other cell types
        residual_k = Xty_i[k]
        for j in range(n_cell_types):
            if j != k:
                residual_k -= XtX[k, j] * beta_i[j]

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

    return beta_i


@jit(nopython=True, parallel=True, cache=True)
def bcd_iteration(
    H: np.ndarray,
    XtX: np.ndarray,
    beta: np.ndarray,
    neighbor_indices: np.ndarray,
    neighbor_indptr: np.ndarray,
    lambda_: float,
    rho: float,
) -> np.ndarray:
    """
    Single BCD iteration over all spots (parallelized).

    Uses precomputed H = X_sketch @ Y_sketch.T for efficiency.

    Parameters
    ----------
    H : ndarray of shape (n_cell_types, n_spots)
        Precomputed X_sketch @ Y_sketch.T matrix.
    XtX : ndarray of shape (n_cell_types, n_cell_types)
        Precomputed Gram matrix.
    beta : ndarray of shape (n_spots, n_cell_types)
        Current cell type abundances.
    neighbor_indices : ndarray
        Flattened neighbor indices (CSR format).
    neighbor_indptr : ndarray
        Index pointers for neighbors (CSR format).
    lambda_ : float
        Spatial regularization.
    rho : float
        Sparsity regularization.

    Returns
    -------
    beta_new : ndarray of shape (n_spots, n_cell_types)
        Updated abundances.
    """
    n_spots = beta.shape[0]
    n_cell_types = beta.shape[1]

    beta_new = beta.copy()

    for i in prange(n_spots):
        # Get precomputed X @ y_i from H matrix (O(K) lookup instead of O(K*d) computation)
        Xty_i = H[:, i]
        beta_i = beta_new[i].copy()

        # Get neighbors
        start = neighbor_indptr[i]
        end = neighbor_indptr[i + 1]
        n_neighbors = end - start

        # Compute neighbor sum
        neighbor_sum = np.zeros(n_cell_types)
        for idx in range(start, end):
            j = neighbor_indices[idx]
            for k in range(n_cell_types):
                neighbor_sum[k] += beta[j, k]

        # Update this spot
        beta_new[i] = update_spot_with_Xty(
            Xty_i, XtX, beta_i, neighbor_sum, n_neighbors, lambda_, rho
        )

    return beta_new


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
    Y_sketch: np.ndarray,
    X_sketch: np.ndarray,
    beta: np.ndarray,
    L: sparse.spmatrix,
    lambda_: float,
    rho: float,
) -> float:
    """
    Compute the objective function value.

    L(beta) = 0.5 * ||Y - beta @ X||_F^2 + lambda * Tr(beta^T L beta) + rho * ||beta||_1

    Parameters
    ----------
    Y_sketch, X_sketch, beta : arrays
        Data and solution.
    L : sparse matrix
        Graph Laplacian.
    lambda_, rho : float
        Regularization parameters.

    Returns
    -------
    obj : float
        Objective value.
    """
    # Fidelity term
    residual = Y_sketch - beta @ X_sketch
    fidelity = 0.5 * np.sum(residual ** 2)

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

    # Initialize beta (uniform proportions)
    beta = np.ones((n_spots, n_cell_types), dtype=np.float64) / n_cell_types

    # Precompute matrices for efficiency
    # G = X @ X^T: Gram matrix (K x K)
    XtX = precompute_gram_matrix(X_sketch)
    # H = X @ Y^T: Avoids recomputing X @ y_i in each iteration (K x N)
    H = precompute_XtY(X_sketch, Y_sketch)

    # Convert adjacency to CSR for efficient neighbor access
    A_csr = A.tocsr()
    neighbor_indices = A_csr.indices.astype(np.int64)
    neighbor_indptr = A_csr.indptr.astype(np.int64)

    # Compute Laplacian for objective
    from flashdeconv.core.spatial import compute_laplacian
    L = compute_laplacian(A)

    # Track convergence
    objectives = []
    converged = False

    for iteration in range(max_iter):
        beta_old = beta.copy()

        # BCD iteration (using precomputed H)
        beta = bcd_iteration(
            H, XtX, beta,
            neighbor_indices, neighbor_indptr,
            lambda_, rho
        )

        # Check convergence
        change = np.max(np.abs(beta - beta_old))
        rel_change = change / (np.max(np.abs(beta_old)) + 1e-10)

        if verbose and (iteration % 10 == 0 or iteration == max_iter - 1):
            obj = compute_objective(Y_sketch, X_sketch, beta, L, lambda_, rho)
            objectives.append(obj)
            print(f"Iteration {iteration}: objective = {obj:.6f}, rel_change = {rel_change:.6e}")

        if rel_change < tol:
            converged = True
            if verbose:
                print(f"Converged at iteration {iteration}")
            break

    # Final objective
    final_obj = compute_objective(Y_sketch, X_sketch, beta, L, lambda_, rho)

    info = {
        'converged': converged,
        'n_iterations': iteration + 1,
        'final_objective': final_obj,
        'objectives': objectives if verbose else [],
        'final_change': rel_change,
    }

    return beta, info


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
    row_sums = np.maximum(row_sums, 1e-10)
    return beta / row_sums
