"""
Spatial graph Laplacian regularization for FlashDeconv.

This module implements the graph Laplacian for spatial smoothing,
which encourages neighboring spots to have similar cell type compositions.
"""

import numpy as np
from scipy import sparse
from typing import Union, Optional, Tuple, List

ArrayLike = Union[np.ndarray, sparse.spmatrix]


def compute_degree_matrix(A: sparse.spmatrix) -> sparse.dia_matrix:
    """
    Compute the degree matrix from adjacency matrix.

    D_ii = sum_j A_ij (number of neighbors for node i)

    Parameters
    ----------
    A : sparse matrix of shape (n, n)
        Adjacency matrix.

    Returns
    -------
    D : sparse.dia_matrix of shape (n, n)
        Diagonal degree matrix.
    """
    degrees = np.asarray(A.sum(axis=1)).flatten()
    return sparse.diags(degrees, format='dia')


def compute_laplacian(
    A: sparse.spmatrix,
    normalized: bool = False,
) -> sparse.csr_matrix:
    """
    Compute the graph Laplacian matrix.

    L = D - A (unnormalized)
    L_norm = I - D^{-1/2} A D^{-1/2} (normalized)

    Parameters
    ----------
    A : sparse matrix of shape (n_spots, n_spots)
        Adjacency matrix.
    normalized : bool, default=False
        Whether to compute normalized Laplacian.

    Returns
    -------
    L : sparse.csr_matrix of shape (n_spots, n_spots)
        Graph Laplacian matrix.
    """
    n = A.shape[0]
    D = compute_degree_matrix(A)

    if normalized:
        # L_norm = I - D^{-1/2} A D^{-1/2}
        degrees = np.asarray(A.sum(axis=1)).flatten()
        degrees_inv_sqrt = np.zeros_like(degrees)
        nonzero = degrees > 0
        degrees_inv_sqrt[nonzero] = 1.0 / np.sqrt(degrees[nonzero])

        D_inv_sqrt = sparse.diags(degrees_inv_sqrt, format='dia')
        L = sparse.eye(n) - D_inv_sqrt @ A @ D_inv_sqrt
    else:
        # L = D - A
        L = D - A

    return L.tocsr()


def get_neighbor_indices(A: sparse.spmatrix) -> List[np.ndarray]:
    """
    Extract neighbor indices for each spot from adjacency matrix.

    Parameters
    ----------
    A : sparse matrix of shape (n_spots, n_spots)
        Adjacency matrix.

    Returns
    -------
    neighbors : list of ndarray
        neighbors[i] contains indices of spots neighboring spot i.
    """
    A_csr = A.tocsr()
    n_spots = A_csr.shape[0]

    neighbors = []
    for i in range(n_spots):
        start, end = A_csr.indptr[i], A_csr.indptr[i + 1]
        neighbors.append(A_csr.indices[start:end].copy())

    return neighbors


def get_neighbor_counts(A: sparse.spmatrix) -> np.ndarray:
    """
    Get number of neighbors for each spot.

    Parameters
    ----------
    A : sparse matrix of shape (n_spots, n_spots)
        Adjacency matrix.

    Returns
    -------
    counts : ndarray of shape (n_spots,)
        Number of neighbors per spot.
    """
    return np.asarray(A.sum(axis=1)).flatten().astype(np.int32)


def compute_laplacian_quadratic(
    beta: np.ndarray,
    L: sparse.spmatrix,
) -> float:
    """
    Compute the Laplacian regularization term: Tr(beta^T L beta).

    This equals sum_{(i,j) in E} ||beta_i - beta_j||^2 for unnormalized L.

    Parameters
    ----------
    beta : ndarray of shape (n_spots, n_cell_types)
        Cell type abundance matrix.
    L : sparse matrix of shape (n_spots, n_spots)
        Graph Laplacian.

    Returns
    -------
    value : float
        Value of the quadratic form.
    """
    # Tr(beta^T L beta) = sum_k (beta[:, k]^T L beta[:, k])
    Lbeta = L @ beta
    return np.sum(beta * Lbeta)


def compute_neighbor_sum(
    beta: np.ndarray,
    neighbors: List[np.ndarray],
) -> np.ndarray:
    """
    Compute sum of beta over neighbors for each spot.

    Parameters
    ----------
    beta : ndarray of shape (n_spots, n_cell_types)
        Cell type abundances.
    neighbors : list of ndarray
        Neighbor indices for each spot.

    Returns
    -------
    neighbor_sums : ndarray of shape (n_spots, n_cell_types)
        Sum of neighbor betas for each spot.
    """
    n_spots, n_cell_types = beta.shape
    neighbor_sums = np.zeros((n_spots, n_cell_types), dtype=beta.dtype)

    for i in range(n_spots):
        if len(neighbors[i]) > 0:
            neighbor_sums[i] = np.sum(beta[neighbors[i]], axis=0)

    return neighbor_sums


def auto_tune_lambda(
    Y_sketch: np.ndarray,
    X_sketch: np.ndarray,
    A: sparse.spmatrix,
    alpha: float = 0.005,
) -> float:
    """
    Auto-tune spatial regularization parameter lambda.

    The key insight is that lambda must be scaled relative to XtX to have
    any effect on the optimization. In the BCD update:
        denom = XtX[k,k] + lambda * n_neighbors
    For lambda to contribute meaningfully, we need:
        lambda * n_neighbors ~ alpha * XtX[k,k]

    Parameters
    ----------
    Y_sketch : ndarray of shape (n_spots, sketch_dim)
        Sketched spatial data.
    X_sketch : ndarray of shape (n_cell_types, sketch_dim)
        Sketched reference.
    A : sparse matrix of shape (n_spots, n_spots)
        Adjacency matrix.
    alpha : float, default=0.005
        Relative strength of spatial regularization (0-1).
        Controls the contribution of spatial term to the Hessian diagonal:
        - alpha=0: no spatial regularization
        - alpha=0.005: spatial term contributes ~0.5% (recommended)
        - alpha=0.1: spatial term contributes ~10%
        Higher values enforce more spatial smoothing.

    Returns
    -------
    lambda_ : float
        Estimated regularization parameter scaled to data.
    """
    # Compute data scale from XtX
    XtX = X_sketch @ X_sketch.T
    avg_XtX_diag = np.mean(np.diag(XtX))

    # Average number of neighbors
    avg_neighbors = np.mean(np.asarray(A.sum(axis=1)).flatten())

    # Set lambda to balance terms
    # When lambda * n_neighbors = alpha * XtX[k,k],
    # the spatial term contributes alpha/(1+alpha) of the Hessian diagonal
    lambda_ = alpha * avg_XtX_diag / max(avg_neighbors, 1.0)

    return float(lambda_)
