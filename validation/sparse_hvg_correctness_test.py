#!/usr/bin/env python3
"""
Test correctness of sparse HVG selection.

Key insight: log1p(0) = 0, so sparsity IS preserved after log transform!
The simplified version in the demo was wrong because it skipped log transform.

This script implements a numerically equivalent sparse version and verifies
100% overlap with the dense version.
"""

import numpy as np
from scipy import sparse
from scipy.sparse import csr_matrix, diags
import time


def hvg_dense_original(Y_dense: np.ndarray, n_top: int = 2000) -> np.ndarray:
    """
    Original dense HVG selection (from genes.py).
    This is the reference implementation.
    """
    n_genes = Y_dense.shape[1]

    # Normalize by total counts
    total_counts = np.sum(Y_dense, axis=1, keepdims=True)
    total_counts = np.maximum(total_counts, 1)
    Y_norm = Y_dense / total_counts * 10000  # Scale to 10k counts

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

    # Take top by normalized dispersion (simplified - skip filtering for test)
    hvg_idx = np.argsort(normalized_dispersion)[::-1][:n_top]

    return np.sort(hvg_idx), gene_means, gene_vars, normalized_dispersion


def hvg_sparse_equivalent(Y_sparse: csr_matrix, n_top: int = 2000) -> np.ndarray:
    """
    Sparse HVG selection - NUMERICALLY EQUIVALENT to dense version.

    Key insight: log1p(0) = 0, so sparsity is preserved!
    We can do:
    1. Row normalize using diagonal matrix multiplication (stays sparse)
    2. log1p on .data (stays sparse, because log1p(0)=0)
    3. Compute mean/var from sparse matrix
    4. Bin-based normalization (on small gene-level arrays)
    """
    N, G = Y_sparse.shape

    # Step 1: Row normalize (CPM-like)
    # Y_norm = Y / total_counts * 10000
    lib_size = np.array(Y_sparse.sum(axis=1)).flatten()
    lib_size = np.maximum(lib_size, 1.0)

    # D @ Y where D = diag(10000 / lib_size)
    D = diags(10000.0 / lib_size)
    Y_norm = D @ Y_sparse  # Still sparse!

    # Step 2: Log1p transform
    # CRITICAL: log1p(0) = 0, so zeros stay zeros!
    Y_log = Y_norm.copy()
    Y_log.data = np.log1p(Y_log.data)

    # Step 3: Compute mean and variance per gene
    # For sparse: E[X] = sum(X) / N
    gene_means = np.array(Y_log.mean(axis=0)).flatten()

    # For variance: Var(X) = E[X^2] - E[X]^2
    # But we need sample variance (ddof=1): Var = sum((X - mean)^2) / (N-1)
    # = (sum(X^2) - N*mean^2) / (N-1)
    # = (N * E[X^2] - N * mean^2) / (N-1)
    # = N / (N-1) * (E[X^2] - mean^2)

    mean_sq = np.array(Y_log.power(2).mean(axis=0)).flatten()
    gene_vars = N / (N - 1) * (mean_sq - gene_means ** 2)
    gene_vars = np.maximum(gene_vars, 0)  # numerical stability

    # Step 4: Bin-based normalization (same as dense, operates on small arrays)
    n_bins = 20
    valid_means = gene_means[gene_means > 0]
    if len(valid_means) > 0:
        bins = np.percentile(valid_means, np.linspace(0, 100, n_bins + 1))
        bins = np.unique(bins)
    else:
        bins = np.array([0, 1])

    gene_bins = np.digitize(gene_means, bins) - 1
    gene_bins = np.clip(gene_bins, 0, max(len(bins) - 2, 0))

    normalized_dispersion = np.zeros(G)

    for i in range(len(bins) - 1):
        mask = gene_bins == i
        if np.sum(mask) > 1:
            bin_vars = gene_vars[mask]
            bin_mean = np.mean(bin_vars)
            bin_std = np.std(bin_vars) + 1e-10
            normalized_dispersion[mask] = (bin_vars - bin_mean) / bin_std

    # Select top genes
    hvg_idx = np.argsort(normalized_dispersion)[::-1][:n_top]

    return np.sort(hvg_idx), gene_means, gene_vars, normalized_dispersion


def test_equivalence():
    """Test that sparse and dense versions produce identical results."""
    print("=" * 60)
    print("Testing Sparse HVG Correctness")
    print("=" * 60)

    # Create test data with realistic sparsity
    N, G = 5000, 10000
    density = 0.15
    n_top = 2000

    print(f"\nTest parameters: N={N}, G={G}, density={density}, n_top={n_top}")

    np.random.seed(42)

    # Create sparse matrix with Poisson-like counts
    nnz = int(N * G * density)
    data = np.random.exponential(3, nnz).astype(np.float64)
    rows = np.random.randint(0, N, nnz)
    cols = np.random.randint(0, G, nnz)
    Y_sparse = csr_matrix((data, (rows, cols)), shape=(N, G))

    # Convert to dense for comparison
    Y_dense = Y_sparse.toarray()

    print(f"Sparse matrix: {Y_sparse.data.nbytes / 1e6:.1f} MB")
    print(f"Dense matrix: {Y_dense.nbytes / 1e6:.1f} MB")

    # Run dense version
    print("\n--- Dense Version ---")
    t0 = time.time()
    hvg_dense, means_dense, vars_dense, disp_dense = hvg_dense_original(Y_dense, n_top)
    t_dense = time.time() - t0
    print(f"Time: {t_dense:.3f}s")

    # Run sparse version
    print("\n--- Sparse Version ---")
    t0 = time.time()
    hvg_sparse, means_sparse, vars_sparse, disp_sparse = hvg_sparse_equivalent(Y_sparse, n_top)
    t_sparse = time.time() - t0
    print(f"Time: {t_sparse:.3f}s")

    # Compare intermediate results
    print("\n--- Numerical Comparison ---")

    mean_diff = np.max(np.abs(means_dense - means_sparse))
    print(f"Gene means max diff: {mean_diff:.2e}")

    var_diff = np.max(np.abs(vars_dense - vars_sparse))
    print(f"Gene vars max diff: {var_diff:.2e}")

    disp_diff = np.max(np.abs(disp_dense - disp_sparse))
    print(f"Norm dispersion max diff: {disp_diff:.2e}")

    # Compare HVG selection
    print("\n--- HVG Selection Comparison ---")
    overlap = len(set(hvg_dense) & set(hvg_sparse))
    print(f"HVG overlap: {overlap}/{n_top} ({100*overlap/n_top:.1f}%)")

    # Check if they are identical
    is_identical = np.array_equal(hvg_dense, hvg_sparse)
    print(f"Identical HVG sets: {is_identical}")

    if not is_identical:
        # Find differences
        only_dense = set(hvg_dense) - set(hvg_sparse)
        only_sparse = set(hvg_sparse) - set(hvg_dense)
        print(f"Only in dense: {len(only_dense)} genes")
        print(f"Only in sparse: {len(only_sparse)} genes")

        # Check if differences are due to ties
        if len(only_dense) > 0:
            example = list(only_dense)[0]
            dense_rank = np.where(np.argsort(disp_dense)[::-1] == example)[0][0]
            sparse_rank = np.where(np.argsort(disp_sparse)[::-1] == example)[0][0]
            print(f"Example gene {example}: dense_rank={dense_rank}, sparse_rank={sparse_rank}")
            print(f"  disp_dense={disp_dense[example]:.6f}, disp_sparse={disp_sparse[example]:.6f}")

    # Verify numerical precision is the issue, not algorithm difference
    print("\n--- Precision Analysis ---")
    rtol = 1e-10
    means_close = np.allclose(means_dense, means_sparse, rtol=rtol)
    vars_close = np.allclose(vars_dense, vars_sparse, rtol=rtol)
    disp_close = np.allclose(disp_dense, disp_sparse, rtol=rtol)

    print(f"Means close (rtol={rtol}): {means_close}")
    print(f"Vars close (rtol={rtol}): {vars_close}")
    print(f"Disp close (rtol={rtol}): {disp_close}")

    return overlap == n_top


def test_memory_efficiency():
    """Test memory efficiency of sparse version."""
    import tracemalloc

    print("\n" + "=" * 60)
    print("Testing Memory Efficiency")
    print("=" * 60)

    N, G = 20000, 20000
    density = 0.1
    n_top = 2000

    print(f"\nLarger test: N={N}, G={G}, density={density}")

    np.random.seed(42)
    nnz = int(N * G * density)
    data = np.random.exponential(3, nnz).astype(np.float64)
    rows = np.random.randint(0, N, nnz)
    cols = np.random.randint(0, G, nnz)
    Y_sparse = csr_matrix((data, (rows, cols)), shape=(N, G))

    sparse_size = Y_sparse.data.nbytes / 1e6
    dense_size = N * G * 8 / 1e6
    print(f"Sparse input: {sparse_size:.1f} MB")
    print(f"Dense equivalent: {dense_size:.1f} MB")

    # Sparse version
    print("\n--- Sparse Version ---")
    tracemalloc.start()
    t0 = time.time()
    hvg_sparse, _, _, _ = hvg_sparse_equivalent(Y_sparse, n_top)
    t_sparse = time.time() - t0
    current, peak_sparse = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    print(f"Time: {t_sparse:.2f}s")
    print(f"Peak memory: {peak_sparse / 1e6:.1f} MB")

    # Dense version
    print("\n--- Dense Version ---")
    Y_dense = Y_sparse.toarray()
    tracemalloc.start()
    t0 = time.time()
    hvg_dense, _, _, _ = hvg_dense_original(Y_dense, n_top)
    t_dense = time.time() - t0
    current, peak_dense = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    print(f"Time: {t_dense:.2f}s")
    print(f"Peak memory: {peak_dense / 1e6:.1f} MB")

    print(f"\n--- Summary ---")
    print(f"Memory reduction: {peak_dense / peak_sparse:.1f}x")
    print(f"Speed improvement: {t_dense / t_sparse:.1f}x")


if __name__ == "__main__":
    success = test_equivalence()
    test_memory_efficiency()

    print("\n" + "=" * 60)
    if success:
        print("SUCCESS: Sparse version is numerically equivalent!")
    else:
        print("WARNING: Results differ - need investigation")
    print("=" * 60)
