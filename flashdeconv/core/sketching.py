"""
Structure-preserving randomized sketching for FlashDeconv.

This module implements the dimensionality reduction via sparse
CountSketch matrices with importance sampling to preserve
signals from rare cell types.
"""

import numpy as np
from scipy import sparse
from typing import Union, Optional, Tuple

ArrayLike = Union[np.ndarray, sparse.spmatrix]


def build_countsketch_matrix(
    n_genes: int,
    sketch_dim: int,
    leverage_scores: Optional[np.ndarray] = None,
    random_state: Optional[int] = None,
) -> sparse.csr_matrix:
    """
    Build a sparse CountSketch matrix with importance sampling.

    CountSketch assigns each row (gene) to exactly one column (sketch dimension)
    with a random sign (+1 or -1). Importance sampling weights the assignment
    probabilities by leverage scores.

    Parameters
    ----------
    n_genes : int
        Number of genes (rows of the sketch matrix).
    sketch_dim : int
        Sketch dimension d (columns of the sketch matrix).
    leverage_scores : ndarray of shape (n_genes,), optional
        Importance weights for each gene. Uniform if not provided.
    random_state : int, optional
        Random seed for reproducibility.

    Returns
    -------
    Omega : sparse.csr_matrix of shape (n_genes, sketch_dim)
        Sparse CountSketch matrix.
    """
    if random_state is not None:
        np.random.seed(random_state)

    # Default to uniform weights
    if leverage_scores is None:
        leverage_scores = np.ones(n_genes) / n_genes
    else:
        # Normalize to probability distribution
        leverage_scores = leverage_scores / (np.sum(leverage_scores) + 1e-10)

    # Assign each gene to a sketch dimension
    # Higher leverage = higher probability of spreading across multiple dimensions

    # For CountSketch: each gene maps to exactly one bucket
    # But we use leverage scores to influence which bucket

    # Standard CountSketch: uniform random assignment
    # We modify to use weighted sampling but ensure each gene is assigned once

    # Compute bucket assignments
    # Use leverage scores to create a "soft" assignment that respects importance

    # Simple approach: hash-based assignment with leverage-weighted repetition
    bucket_assignments = np.random.randint(0, sketch_dim, size=n_genes)
    signs = np.random.choice([-1, 1], size=n_genes)

    # Scale by sqrt(leverage) to preserve more signal from important genes
    # This is the "importance sampling" twist
    scale_factors = np.sqrt(leverage_scores * n_genes + 1e-10)
    scale_factors = np.clip(scale_factors, 0.1, 10.0)  # Prevent extreme values

    # Build sparse matrix
    row_idx = np.arange(n_genes)
    col_idx = bucket_assignments
    data = signs * scale_factors

    Omega = sparse.csr_matrix(
        (data, (row_idx, col_idx)),
        shape=(n_genes, sketch_dim),
        dtype=np.float64,
    )

    # Normalize columns for numerical stability
    col_norms = np.sqrt(np.asarray(Omega.power(2).sum(axis=0)).flatten())
    col_norms = np.maximum(col_norms, 1e-10)

    # Scale by sqrt(n_genes / sketch_dim) to preserve norms approximately
    scale = np.sqrt(n_genes / sketch_dim)
    Omega = Omega.multiply(scale / col_norms)

    return Omega.tocsr()


def build_sparse_rademacher_matrix(
    n_genes: int,
    sketch_dim: int,
    sparsity: float = 0.1,
    leverage_scores: Optional[np.ndarray] = None,
    random_state: Optional[int] = None,
) -> sparse.csr_matrix:
    """
    Build a sparse Rademacher (random sign) matrix.

    Alternative to CountSketch with controllable sparsity.
    Each entry is 0 with probability (1-sparsity), or ±1/sqrt(sparsity)
    with probability sparsity/2 each.

    Parameters
    ----------
    n_genes : int
        Number of genes.
    sketch_dim : int
        Sketch dimension.
    sparsity : float, default=0.1
        Fraction of non-zero entries per column.
    leverage_scores : ndarray, optional
        Importance weights (higher = more likely to be non-zero).
    random_state : int, optional
        Random seed.

    Returns
    -------
    Omega : sparse.csr_matrix of shape (n_genes, sketch_dim)
        Sparse Rademacher matrix.
    """
    if random_state is not None:
        np.random.seed(random_state)

    if leverage_scores is None:
        leverage_scores = np.ones(n_genes) / n_genes
    else:
        leverage_scores = leverage_scores / (np.sum(leverage_scores) + 1e-10)

    # Compute per-gene sparsity based on leverage
    # Higher leverage = higher probability of being sampled
    gene_probs = sparsity * (1 + leverage_scores * n_genes)
    gene_probs = np.clip(gene_probs, 0.01, 1.0)

    rows, cols, data = [], [], []

    scale = 1.0 / np.sqrt(sparsity * n_genes / sketch_dim)

    for j in range(sketch_dim):
        # Sample genes for this column
        mask = np.random.random(n_genes) < gene_probs
        selected_genes = np.where(mask)[0]

        if len(selected_genes) == 0:
            # Ensure at least one gene per column
            selected_genes = np.array([np.random.randint(n_genes)])

        # Random signs
        signs = np.random.choice([-1, 1], size=len(selected_genes))

        rows.extend(selected_genes)
        cols.extend([j] * len(selected_genes))
        data.extend(signs * scale)

    Omega = sparse.csr_matrix(
        (data, (rows, cols)),
        shape=(n_genes, sketch_dim),
        dtype=np.float64,
    )

    return Omega


def project_to_sketch(
    Y_tilde: Union[np.ndarray, sparse.spmatrix],
    X_tilde: np.ndarray,
    Omega: sparse.spmatrix,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Project data into low-dimensional sketch space.

    Y_sketch = Y_tilde @ Omega  (N x d)
    X_sketch = X_tilde @ Omega  (K x d)

    Supports sparse Y_tilde for memory efficiency. The output is always
    dense since the sketch dimension d is small.

    Parameters
    ----------
    Y_tilde : array-like of shape (n_spots, n_genes)
        Transformed spatial data (sparse or dense).
    X_tilde : ndarray of shape (n_cell_types, n_genes)
        Transformed reference signatures.
    Omega : sparse matrix of shape (n_genes, sketch_dim)
        Sketching matrix.

    Returns
    -------
    Y_sketch : ndarray of shape (n_spots, sketch_dim)
        Projected spatial data (always dense).
    X_sketch : ndarray of shape (n_cell_types, sketch_dim)
        Projected reference signatures.
    """
    # Ensure Omega is in efficient format
    if sparse.issparse(Omega):
        Omega = Omega.tocsr()

    # Project spatial data (sparse @ sparse works efficiently)
    Y_sketch = Y_tilde @ Omega

    # Ensure output is dense (needed for downstream BCD solver)
    if sparse.issparse(Y_sketch):
        Y_sketch = Y_sketch.toarray()

    # Project reference
    X_sketch = X_tilde @ Omega
    if sparse.issparse(X_sketch):
        X_sketch = X_sketch.toarray()

    return Y_sketch, X_sketch


def sketch_data(
    Y_tilde: Union[np.ndarray, sparse.spmatrix],
    X_tilde: np.ndarray,
    sketch_dim: int = 512,
    leverage_scores: Optional[np.ndarray] = None,
    method: str = "countsketch",
    random_state: Optional[int] = None,
) -> Tuple[np.ndarray, np.ndarray, sparse.spmatrix]:
    """
    Full sketching pipeline.

    Parameters
    ----------
    Y_tilde : ndarray of shape (n_spots, n_genes)
        Transformed spatial data.
    X_tilde : ndarray of shape (n_cell_types, n_genes)
        Transformed reference.
    sketch_dim : int, default=512
        Target dimension.
    leverage_scores : ndarray, optional
        Gene importance weights.
    method : str, default="countsketch"
        Sketching method ("countsketch" or "rademacher").
    random_state : int, optional
        Random seed.

    Returns
    -------
    Y_sketch : ndarray of shape (n_spots, sketch_dim)
        Sketched spatial data.
    X_sketch : ndarray of shape (n_cell_types, sketch_dim)
        Sketched reference.
    Omega : sparse matrix
        The sketching matrix used.
    """
    n_genes = Y_tilde.shape[1]

    if method == "countsketch":
        Omega = build_countsketch_matrix(
            n_genes, sketch_dim, leverage_scores, random_state
        )
    elif method == "rademacher":
        Omega = build_sparse_rademacher_matrix(
            n_genes, sketch_dim, leverage_scores=leverage_scores,
            random_state=random_state
        )
    else:
        raise ValueError(f"Unknown sketching method: {method}")

    Y_sketch, X_sketch = project_to_sketch(Y_tilde, X_tilde, Omega)

    return Y_sketch, X_sketch, Omega
