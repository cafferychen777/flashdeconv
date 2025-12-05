"""
Gene selection utilities for FlashDeconv.

This module implements feature selection for structure-preserving sketching:
- Highly variable gene (HVG) selection
- Cell-type specific marker gene identification
- Leverage score computation for importance weighting
"""

import numpy as np
from scipy import sparse
from scipy.sparse import diags
from typing import Union, Optional, Tuple, List

ArrayLike = Union[np.ndarray, sparse.spmatrix]


def select_hvg(
    Y: ArrayLike,
    n_top: int = 2000,
    min_mean: float = 0.0125,
    max_mean: float = 3.0,
    min_disp: float = 0.5,
) -> np.ndarray:
    """
    Select highly variable genes using the Seurat v3 method.

    Supports both sparse and dense input matrices. For sparse matrices,
    operations are performed without converting to dense, preserving
    memory efficiency. Key insight: log1p(0) = 0, so sparsity is preserved.

    Parameters
    ----------
    Y : array-like of shape (n_spots, n_genes)
        Raw count matrix (sparse or dense).
    n_top : int, default=2000
        Number of top HVGs to select.
    min_mean : float, default=0.0125
        Minimum mean expression.
    max_mean : float, default=3.0
        Maximum mean expression.
    min_disp : float, default=0.5
        Minimum dispersion.

    Returns
    -------
    hvg_idx : ndarray of shape (n_hvg,)
        Indices of highly variable genes.
    """
    N, n_genes = Y.shape

    if sparse.issparse(Y):
        # Sparse-friendly implementation (memory efficient)
        # Key insight: log1p(0) = 0, so sparsity is preserved!

        # Step 1: Row normalize using diagonal matrix multiplication
        lib_size = np.array(Y.sum(axis=1)).flatten()
        lib_size = np.maximum(lib_size, 1.0)
        D = diags(10000.0 / lib_size)
        Y_norm = D @ Y  # Still sparse!

        # Step 2: Log1p transform in-place - zeros stay zeros!
        # No copy needed: directly modify .data array
        Y_norm.data = np.log1p(Y_norm.data)

        # Step 3: Compute mean and variance per gene (memory-efficient)
        # Use sum instead of mean to avoid scipy's internal dense conversion
        gene_sums = np.array(Y_norm.sum(axis=0)).flatten()
        gene_means = gene_sums / N

        # Compute E[X^2] without creating power(2) matrix copy
        # bincount efficiently computes column-wise sum of squares
        data_sq = Y_norm.data ** 2  # O(nnz) temporary
        col_sum_sq = np.bincount(Y_norm.indices, weights=data_sq, minlength=n_genes)
        mean_sq = col_sum_sq / N

        # Sample variance (ddof=1): Var = N/(N-1) * (E[X^2] - E[X]^2)
        gene_vars = N / (N - 1) * (mean_sq - gene_means ** 2)
        gene_vars = np.maximum(gene_vars, 0)  # numerical stability

    else:
        # Dense implementation (original)
        Y_dense = np.asarray(Y)

        # Normalize by total counts
        total_counts = np.sum(Y_dense, axis=1, keepdims=True)
        total_counts = np.maximum(total_counts, 1)
        Y_norm = Y_dense / total_counts * 10000

        # Log transform
        Y_log = np.log1p(Y_norm)

        # Compute mean and variance
        gene_means = np.mean(Y_log, axis=0)
        gene_vars = np.var(Y_log, axis=0, ddof=1)

    # Compute standardized dispersion
    # Bin genes by mean expression
    n_bins = 20
    bins = np.percentile(gene_means[gene_means > 0], np.linspace(0, 100, n_bins + 1))
    bins = np.unique(bins)

    # Assign genes to bins
    gene_bins = np.digitize(gene_means, bins) - 1
    gene_bins = np.clip(gene_bins, 0, len(bins) - 2)

    # Compute normalized dispersion within each bin
    normalized_dispersion = np.zeros(n_genes)

    for i in range(len(bins) - 1):
        mask = gene_bins == i
        if np.sum(mask) > 1:
            bin_vars = gene_vars[mask]
            bin_mean = np.mean(bin_vars)
            bin_std = np.std(bin_vars) + 1e-10
            normalized_dispersion[mask] = (bin_vars - bin_mean) / bin_std

    # Filter by mean expression
    mean_mask = (gene_means >= min_mean) & (gene_means <= max_mean)
    disp_mask = normalized_dispersion >= min_disp

    valid_mask = mean_mask & disp_mask
    valid_idx = np.where(valid_mask)[0]

    if len(valid_idx) < n_top:
        # If not enough genes pass filters, take top by dispersion
        sorted_idx = np.argsort(normalized_dispersion)[::-1]
        hvg_idx = sorted_idx[:n_top]
    else:
        # Take top n_top by normalized dispersion
        valid_disp = normalized_dispersion[valid_idx]
        top_idx = np.argsort(valid_disp)[::-1][:n_top]
        hvg_idx = valid_idx[top_idx]

    return np.sort(hvg_idx)


def select_markers(
    X: np.ndarray,
    n_markers: int = 50,
    method: str = "diff",
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Select cell-type specific marker genes.

    Parameters
    ----------
    X : ndarray of shape (n_cell_types, n_genes)
        Reference signature matrix.
    n_markers : int, default=50
        Number of markers per cell type.
    method : str, default="diff"
        Selection method:
        - "diff": Difference between max and second max
        - "ratio": Ratio of expression to others
        - "specificity": Specificity score

    Returns
    -------
    marker_idx : ndarray
        Indices of marker genes (union across all cell types).
    marker_assignments : ndarray of shape (n_markers * n_cell_types,)
        Cell type assignment for each marker.
    """
    n_cell_types, n_genes = X.shape

    # Normalize X row-wise
    X_norm = X / (np.sum(X, axis=1, keepdims=True) + 1e-10)

    if method == "diff":
        # For each gene, compute difference between top and second expression
        sorted_expr = np.sort(X_norm, axis=0)[::-1]
        specificity = sorted_expr[0] - sorted_expr[1]  # Max - second max

    elif method == "ratio":
        # Ratio of max to mean of others
        max_expr = np.max(X_norm, axis=0)
        mean_others = (np.sum(X_norm, axis=0) - max_expr) / (n_cell_types - 1 + 1e-10)
        specificity = max_expr / (mean_others + 1e-10)

    elif method == "specificity":
        # Tau specificity score
        max_expr = np.max(X_norm, axis=0)
        specificity = np.sum(1 - X_norm / (max_expr + 1e-10), axis=0) / (n_cell_types - 1)

    else:
        raise ValueError(f"Unknown method: {method}")

    # For each cell type, find top markers
    all_markers = []
    marker_assignments = []

    for k in range(n_cell_types):
        # Genes where this cell type has highest expression
        is_top = np.argmax(X_norm, axis=0) == k
        top_genes = np.where(is_top)[0]

        if len(top_genes) > 0:
            # Rank by specificity
            gene_spec = specificity[top_genes]
            sorted_idx = np.argsort(gene_spec)[::-1][:n_markers]
            markers_k = top_genes[sorted_idx]
        else:
            # Fallback: genes with highest expression in this cell type
            sorted_idx = np.argsort(X_norm[k])[::-1][:n_markers]
            markers_k = sorted_idx

        all_markers.extend(markers_k)
        marker_assignments.extend([k] * len(markers_k))

    marker_idx = np.unique(all_markers)

    return marker_idx, np.array(marker_assignments)


def compute_leverage_scores(
    X: np.ndarray,
    regularization: float = 1e-6,
) -> np.ndarray:
    """
    Compute leverage scores for each gene.

    Leverage scores measure the importance of each gene in the
    column space of the reference matrix. Genes with higher
    leverage scores are more important for distinguishing cell types.

    Parameters
    ----------
    X : ndarray of shape (n_cell_types, n_genes)
        Reference signature matrix.
    regularization : float, default=1e-6
        Regularization for numerical stability.

    Returns
    -------
    leverage : ndarray of shape (n_genes,)
        Leverage scores for each gene.
    """
    n_cell_types, n_genes = X.shape

    # Center the data
    X_centered = X - np.mean(X, axis=0, keepdims=True)

    # Compute SVD of X.T (genes x cell_types)
    # X.T = U @ diag(s) @ Vt
    # U has shape (n_genes, k), each row is a gene's projection onto PCs
    try:
        U, s, Vt = np.linalg.svd(X_centered.T, full_matrices=False)
    except np.linalg.LinAlgError:
        # Fallback to simple variance-based scores
        var_scores = np.var(X, axis=0)
        return var_scores / (var_scores.sum() + regularization)

    # Leverage scores are the row norms of U (weighted by singular values)
    # U has shape (n_genes, k) where k = min(n_genes, n_cell_types)
    k = min(n_cell_types, n_genes, len(s))

    # Weight by singular values squared (importance of each PC)
    weights = s[:k] ** 2 / (s[:k] ** 2 + regularization)

    # Compute leverage: sum of squared loadings weighted by PC importance
    # U[:, :k] has shape (n_genes, k), result has shape (n_genes,)
    leverage = np.sum((U[:, :k] ** 2) * weights, axis=1)

    # Normalize to sum to 1
    leverage = leverage / (np.sum(leverage) + regularization)

    return leverage


def select_informative_genes(
    Y: ArrayLike,
    X: np.ndarray,
    n_hvg: int = 2000,
    n_markers_per_type: int = 50,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Select informative genes for sketching.

    Combines highly variable genes from spatial data with
    cell-type specific markers from reference.

    Parameters
    ----------
    Y : array-like of shape (n_spots, n_genes)
        Raw spatial count matrix.
    X : ndarray of shape (n_cell_types, n_genes)
        Reference signature matrix.
    n_hvg : int, default=2000
        Number of highly variable genes.
    n_markers_per_type : int, default=50
        Number of markers per cell type.

    Returns
    -------
    gene_idx : ndarray
        Indices of selected informative genes.
    leverage_scores : ndarray
        Leverage scores for selected genes.
    """
    # Select HVGs
    hvg_idx = select_hvg(Y, n_top=n_hvg)

    # Select markers
    marker_idx, _ = select_markers(X, n_markers=n_markers_per_type)

    # Union of HVGs and markers
    gene_idx = np.union1d(hvg_idx, marker_idx)

    # Compute leverage scores for selected genes
    X_subset = X[:, gene_idx]
    leverage_scores = compute_leverage_scores(X_subset)

    return gene_idx, leverage_scores
