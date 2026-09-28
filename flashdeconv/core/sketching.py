"""
Structure-preserving gene representations for FlashDeconv.

This module implements two gene-space representations:

- Expected leverage-weighted CountSketch weights (default). Each selected
  gene is kept (no projection) and scaled by the exact expectation, over the
  random bucket assignment, of its weight in the column-normalized
  leverage-weighted CountSketch. The result is deterministic.
- The randomized leverage-weighted CountSketch projection (legacy), which
  hashes genes into ``sketch_dim`` buckets with random signs.
"""

import numpy as np
from scipy import sparse
from typing import Union, Optional, Tuple

from flashdeconv.utils.random import check_random_state

ArrayLike = Union[np.ndarray, sparse.spmatrix]


def build_countsketch_matrix(
    n_genes: int,
    sketch_dim: int,
    leverage_scores: Optional[np.ndarray] = None,
    random_state: Optional[int] = None,
) -> sparse.csr_matrix:
    """
    Build a sparse CountSketch matrix with leverage-based amplitude scaling.

    CountSketch assigns each row (gene) to exactly one column (sketch dimension)
    via uniform random hashing, with a random sign (+1 or -1). Leverage scores
    scale the entry amplitude so that high-leverage genes contribute more to the
    sketch, preserving signal from rare cell types.

    Parameters
    ----------
    n_genes : int
        Number of genes (rows of the sketch matrix).
    sketch_dim : int
        Sketch dimension d (columns of the sketch matrix).
    leverage_scores : ndarray of shape (n_genes,), optional
        Importance weights for amplitude scaling. Uniform if not provided.
    random_state : int, optional
        Random seed for reproducibility.

    Returns
    -------
    Omega : sparse.csr_matrix of shape (n_genes, sketch_dim)
        Sparse CountSketch matrix.
    """
    rng = check_random_state(random_state)

    # Default to uniform weights
    if leverage_scores is None:
        leverage_scores = np.ones(n_genes) / n_genes
    else:
        # Normalize to probability distribution
        leverage_scores = leverage_scores / (np.sum(leverage_scores) + 1e-10)

    # Standard CountSketch: uniform random bucket assignment + random sign
    bucket_assignments = rng.randint(0, sketch_dim, size=n_genes)
    signs = rng.choice([-1, 1], size=n_genes)

    # Scale amplitude by sqrt(leverage) so high-leverage genes contribute more
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


def expected_countsketch_weights(
    leverage_scores: Optional[np.ndarray],
    n_genes: int,
    sketch_dim: int = 512,
    n_quad: int = 3000,
    max_block_elements: int = 2 ** 24,
) -> np.ndarray:
    """
    Exact expected per-gene weight of the leverage-weighted CountSketch.

    In :func:`build_countsketch_matrix`, gene ``j`` receives the amplitude
    ``a_j = clip(sqrt(g * l_j), 0.1, 10)`` (``l_j`` the normalized leverage,
    ``g`` the number of genes) and is hashed to one of ``d`` buckets; each
    bucket (column) is then normalized to unit norm. The absolute weight of
    gene ``j`` is therefore ``a_j / sqrt(a_j^2 + S_j)``, where
    ``S_j = sum_{i != j} a_i^2 B_i`` and ``B_i ~ Bernoulli(1/d)`` indicates
    that gene ``i`` shares the bucket of gene ``j``. This function returns
    its exact expectation under random hashing (the global
    ``sqrt(g / d)`` factor is omitted; it does not change the solution with
    automatic ``lambda_spatial``):

        w_j = a_j * pi^(-1/2) * int_0^inf t^(-1/2) exp(-t a_j^2)
              * prod_{i != j} [1 - (1 - exp(-t a_i^2)) / d] dt,

    using ``x^(-1/2) = pi^(-1/2) int_0^inf t^(-1/2) exp(-t x) dt``. The
    integral is evaluated with the trapezoid rule in ``log t`` on
    ``n_quad`` log-spaced nodes in ``[1e-7, 1e5]``; the product is computed
    in log space. The result is deterministic (no hashing).

    Parameters
    ----------
    leverage_scores : ndarray of shape (n_genes,) or None
        Gene leverage scores (normalized internally). Uniform if None.
    n_genes : int
        Number of genes g.
    sketch_dim : int, default=512
        Number of CountSketch buckets d.
    n_quad : int, default=3000
        Number of quadrature nodes.
    max_block_elements : int, default=2**24
        Maximum size of the (genes x nodes) work array; larger problems are
        processed in blocks of quadrature nodes to bound memory.

    Returns
    -------
    weights : ndarray of shape (n_genes,)
        Positive, finite per-gene weights.
    """
    if sketch_dim <= 0:
        raise ValueError(f"sketch_dim must be positive, got {sketch_dim}")
    if leverage_scores is None:
        p = np.ones(n_genes) / n_genes
    else:
        leverage_scores = np.asarray(leverage_scores, dtype=np.float64)
        if leverage_scores.shape != (n_genes,):
            raise ValueError(
                f"leverage_scores must have shape ({n_genes},), "
                f"got {leverage_scores.shape}"
            )
        p = leverage_scores / (np.sum(leverage_scores) + 1e-10)
    if n_genes == 0:
        return np.empty(0, dtype=np.float64)
    a2 = np.clip(np.sqrt(p * n_genes + 1e-10), 0.1, 10.0) ** 2

    if sketch_dim == 1:
        # Every gene shares the single bucket: the weight is deterministic.
        return np.sqrt(a2) / np.sqrt(np.sum(a2))

    lt = np.linspace(np.log(1e-7), np.log(1e5), n_quad)
    t = np.exp(lt)
    trapz = getattr(np, "trapezoid", getattr(np, "trapz", None))

    def _integrand(tb):
        # log of the per-gene factor 1 - (1 - exp(-t a_i^2)) / d
        neg = -np.outer(a2, tb)
        logf = np.log1p(-(1.0 - np.exp(neg)) / sketch_dim)
        tot = logf.sum(axis=0)
        # t^(1/2) accounts for t^(-1/2) dt = t^(1/2) d(log t); the product
        # over i != j is exp(tot - logf_j).
        return np.sqrt(tb)[None, :] * np.exp(neg + tot[None, :] - logf)

    block = max(2, int(max_block_elements // max(n_genes, 1)))
    if block >= n_quad:
        val = trapz(_integrand(t), lt, axis=1)
    else:
        # Trapezoid over consecutive node blocks sharing their end nodes.
        val = np.zeros(n_genes, dtype=np.float64)
        start = 0
        while start < n_quad - 1:
            stop = min(start + block, n_quad)
            val += trapz(_integrand(t[start:stop]), lt[start:stop], axis=1)
            start = stop - 1
    return np.sqrt(a2) * (val / np.sqrt(np.pi))


def apply_gene_weights(
    Y_tilde: Union[np.ndarray, sparse.spmatrix],
    X_tilde: np.ndarray,
    weights: np.ndarray,
) -> Tuple[Union[np.ndarray, sparse.csr_matrix], np.ndarray]:
    """
    Scale each gene column of Y_tilde and X_tilde by its weight.

    Parameters
    ----------
    Y_tilde : array-like of shape (n_spots, n_genes)
        Transformed spatial data (sparse or dense).
    X_tilde : ndarray of shape (n_cell_types, n_genes)
        Transformed reference signatures.
    weights : ndarray of shape (n_genes,)
        Per-gene weights.

    Returns
    -------
    Y_w : array-like of shape (n_spots, n_genes)
        Weighted spatial data (CSR if the input was sparse, else dense).
    X_w : ndarray of shape (n_cell_types, n_genes)
        Weighted reference signatures.
    """
    if sparse.issparse(Y_tilde):
        Y_w = (Y_tilde @ sparse.diags(weights)).tocsr()
    else:
        Y_w = Y_tilde * weights[None, :]
    X_w = np.asarray(X_tilde) * weights[None, :]
    return Y_w, X_w


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
    rng = check_random_state(random_state)

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
        mask = rng.random(n_genes) < gene_probs
        selected_genes = np.where(mask)[0]

        if len(selected_genes) == 0:
            # Ensure at least one gene per column
            selected_genes = np.array([rng.randint(n_genes)])

        # Random signs
        signs = rng.choice([-1, 1], size=len(selected_genes))

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
