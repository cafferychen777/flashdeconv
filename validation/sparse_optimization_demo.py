#!/usr/bin/env python3
"""
Sparse Matrix Optimization Verification Scripts

This script demonstrates and verifies each optimization point in the
sparse matrix optimization plan. Run each demo independently to verify
correctness and memory efficiency.

Usage:
    python sparse_optimization_demo.py --demo hvg
    python sparse_optimization_demo.py --demo logcpm
    python sparse_optimization_demo.py --demo sketch
    python sparse_optimization_demo.py --demo all
"""

import numpy as np
from scipy import sparse
from scipy.sparse import csr_matrix, diags
import time
import argparse
import tracemalloc
from typing import Tuple


# =============================================================================
# Demo 1: Sparse HVG Selection
# =============================================================================

def demo_hvg_dense(Y_dense: np.ndarray, n_top: int = 2000) -> np.ndarray:
    """Original dense HVG selection (for comparison)."""
    # Normalize by total counts
    total_counts = np.sum(Y_dense, axis=1, keepdims=True)
    total_counts = np.maximum(total_counts, 1)
    Y_norm = Y_dense / total_counts * 10000

    # Log transform
    Y_log = np.log1p(Y_norm)

    # Compute mean and variance
    gene_means = np.mean(Y_log, axis=0)
    gene_vars = np.var(Y_log, axis=0, ddof=1)

    # Simple dispersion-based selection
    mean_safe = np.maximum(gene_means, 1e-12)
    dispersion = gene_vars / mean_safe

    # Normalized dispersion
    disp_mean = np.mean(dispersion[dispersion > 0])
    disp_std = np.std(dispersion[dispersion > 0]) + 1e-10
    norm_disp = (dispersion - disp_mean) / disp_std

    hvg_idx = np.argsort(norm_disp)[-n_top:]
    return np.sort(hvg_idx)


def demo_hvg_sparse(Y_sparse: csr_matrix, n_top: int = 2000) -> np.ndarray:
    """
    Sparse-friendly HVG selection.

    Key operations that work on sparse without densifying:
    - Y.mean(axis=0): computes mean efficiently
    - Y.power(2): squares only non-zero elements
    - Y.sum(axis=1): computes row sums

    Variance formula: Var(Y) = E[Y^2] - E[Y]^2
    """
    N, G = Y_sparse.shape

    # Step 1: Compute library sizes (row sums) - sparse friendly
    lib_size = np.array(Y_sparse.sum(axis=1)).flatten()
    lib_size = np.maximum(lib_size, 1.0)

    # Step 2: For log-normalized mean and variance, we need to think carefully.
    # Direct approach: compute mean and E[Y^2] on raw counts, then adjust.
    #
    # Actually, for HVG we care about dispersion = var/mean on log-normalized data.
    # But computing log on sparse is tricky (log(0) = -inf).
    #
    # Simplified approach: Use raw count statistics as proxy
    # This is actually what many methods do (Seurat v3 uses raw counts for HVG)

    # Mean per gene (on raw counts, normalized by total)
    # mean_i = sum_j (Y_ij / lib_j) / N
    # This can be computed as: (1/N) * sum_j (Y_ij / lib_j)

    # Create diagonal matrix D = diag(1/lib_size)
    D_inv = diags(1.0 / lib_size)

    # Y_norm = D_inv @ Y (row-normalized, still sparse)
    Y_norm = D_inv @ Y_sparse

    # Mean per gene
    mean = np.array(Y_norm.mean(axis=0)).flatten()

    # E[Y^2] per gene
    mean_sq = np.array(Y_norm.power(2).mean(axis=0)).flatten()

    # Variance = E[Y^2] - E[Y]^2
    var = mean_sq - mean ** 2
    var = np.maximum(var, 0)  # numerical stability

    # Dispersion = var / mean
    mean_safe = np.maximum(mean, 1e-12)
    dispersion = var / mean_safe

    # Normalized dispersion
    valid_disp = dispersion[dispersion > 0]
    if len(valid_disp) > 0:
        disp_mean = np.mean(valid_disp)
        disp_std = np.std(valid_disp) + 1e-10
    else:
        disp_mean, disp_std = 0, 1

    norm_disp = (dispersion - disp_mean) / disp_std

    # Select top genes
    hvg_idx = np.argsort(norm_disp)[-n_top:]
    return np.sort(hvg_idx)


def run_hvg_demo():
    """Demonstrate sparse vs dense HVG selection."""
    print("=" * 60)
    print("Demo 1: Sparse HVG Selection")
    print("=" * 60)

    # Create test data
    N, G = 10000, 20000  # 10K spots, 20K genes
    density = 0.1  # 10% non-zero

    print(f"\nCreating test data: {N} spots x {G} genes, density={density}")

    # Create sparse matrix
    np.random.seed(42)
    data = np.random.exponential(2, int(N * G * density))
    rows = np.random.randint(0, N, len(data))
    cols = np.random.randint(0, G, len(data))
    Y_sparse = csr_matrix((data, (rows, cols)), shape=(N, G))

    print(f"Sparse matrix size: {Y_sparse.data.nbytes / 1e6:.1f} MB")
    print(f"Dense equivalent: {N * G * 8 / 1e6:.1f} MB")

    # Run sparse version
    print("\n--- Sparse HVG Selection ---")
    tracemalloc.start()
    t0 = time.time()
    hvg_sparse = demo_hvg_sparse(Y_sparse, n_top=2000)
    t_sparse = time.time() - t0
    current, peak_sparse = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    print(f"Time: {t_sparse:.2f}s, Peak memory: {peak_sparse / 1e6:.1f} MB")

    # Run dense version (for comparison, on smaller data)
    print("\n--- Dense HVG Selection (for comparison) ---")
    tracemalloc.start()
    t0 = time.time()
    Y_dense = Y_sparse.toarray()
    hvg_dense = demo_hvg_dense(Y_dense, n_top=2000)
    t_dense = time.time() - t0
    current, peak_dense = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    print(f"Time: {t_dense:.2f}s, Peak memory: {peak_dense / 1e6:.1f} MB")

    # Compare results
    overlap = len(set(hvg_sparse) & set(hvg_dense))
    print(f"\n--- Comparison ---")
    print(f"Overlap in top 2000 HVGs: {overlap}/2000 ({100*overlap/2000:.1f}%)")
    print(f"Memory reduction: {peak_dense / peak_sparse:.1f}x")

    return hvg_sparse, hvg_dense


# =============================================================================
# Demo 2: Sparse Log-CPM Normalization
# =============================================================================

def demo_logcpm_dense(Y_dense: np.ndarray) -> np.ndarray:
    """Original dense Log-CPM (for comparison)."""
    lib_size = Y_dense.sum(axis=1, keepdims=True)
    lib_size = np.maximum(lib_size, 1.0)
    Y_cpm = Y_dense / lib_size * 1e4
    return np.log1p(Y_cpm)


def demo_logcpm_sparse(Y_sparse: csr_matrix) -> csr_matrix:
    """
    Sparse-friendly Log-CPM normalization.

    Key insight: log1p(0) = 0, so zeros stay zeros after log1p.
    We can operate on .data directly.
    """
    N, G = Y_sparse.shape

    # Step 1: Compute library sizes
    lib_size = np.array(Y_sparse.sum(axis=1)).flatten()
    lib_size[lib_size == 0] = 1.0

    # Step 2: Build diagonal scaling matrix
    scale_factors = 1e4 / lib_size
    D = diags(scale_factors)

    # Step 3: CPM = D @ Y (sparse @ sparse = sparse)
    Y_cpm = D @ Y_sparse

    # Step 4: Log1p on .data (preserves sparsity!)
    Y_log = Y_cpm.copy()
    Y_log.data = np.log1p(Y_log.data)

    return Y_log


def run_logcpm_demo():
    """Demonstrate sparse vs dense Log-CPM."""
    print("\n" + "=" * 60)
    print("Demo 2: Sparse Log-CPM Normalization")
    print("=" * 60)

    # Create test data
    N, G = 10000, 4000  # After gene selection
    density = 0.15

    print(f"\nCreating test data: {N} spots x {G} genes, density={density}")

    np.random.seed(42)
    data = np.random.exponential(5, int(N * G * density))
    rows = np.random.randint(0, N, len(data))
    cols = np.random.randint(0, G, len(data))
    Y_sparse = csr_matrix((data, (rows, cols)), shape=(N, G))

    print(f"Sparse matrix size: {Y_sparse.data.nbytes / 1e6:.1f} MB")

    # Sparse version
    print("\n--- Sparse Log-CPM ---")
    tracemalloc.start()
    t0 = time.time()
    Y_log_sparse = demo_logcpm_sparse(Y_sparse)
    t_sparse = time.time() - t0
    current, peak_sparse = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    print(f"Time: {t_sparse:.3f}s, Peak memory: {peak_sparse / 1e6:.1f} MB")
    print(f"Output is sparse: {sparse.issparse(Y_log_sparse)}")
    print(f"Output sparsity: {Y_log_sparse.nnz / (N * G) * 100:.1f}%")

    # Dense version
    print("\n--- Dense Log-CPM ---")
    Y_dense = Y_sparse.toarray()
    tracemalloc.start()
    t0 = time.time()
    Y_log_dense = demo_logcpm_dense(Y_dense)
    t_dense = time.time() - t0
    current, peak_dense = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    print(f"Time: {t_dense:.3f}s, Peak memory: {peak_dense / 1e6:.1f} MB")

    # Verify correctness
    Y_log_sparse_dense = Y_log_sparse.toarray()
    max_diff = np.max(np.abs(Y_log_sparse_dense - Y_log_dense))
    print(f"\n--- Correctness Check ---")
    print(f"Max absolute difference: {max_diff:.2e}")
    print(f"Numerically equivalent: {max_diff < 1e-10}")

    return Y_log_sparse, Y_log_dense


# =============================================================================
# Demo 3: Sparse @ Dense Sketching
# =============================================================================

def demo_sketch_sparse(Y_sparse: csr_matrix, sketch_dim: int = 512) -> np.ndarray:
    """
    Sparse @ Dense matrix multiplication for sketching.

    scipy.sparse handles this efficiently without converting Y to dense.
    """
    G = Y_sparse.shape[1]

    # Build a simple random sketch matrix (dense, small: G x d)
    np.random.seed(42)
    Omega = np.random.randn(G, sketch_dim) / np.sqrt(sketch_dim)

    # Sparse @ Dense = Dense (efficiently computed)
    Y_sketch = Y_sparse @ Omega

    return Y_sketch


def demo_sketch_dense(Y_dense: np.ndarray, sketch_dim: int = 512) -> np.ndarray:
    """Dense sketching (for comparison)."""
    G = Y_dense.shape[1]
    np.random.seed(42)
    Omega = np.random.randn(G, sketch_dim) / np.sqrt(sketch_dim)
    return Y_dense @ Omega


def run_sketch_demo():
    """Demonstrate sparse @ dense sketching."""
    print("\n" + "=" * 60)
    print("Demo 3: Sparse @ Dense Sketching")
    print("=" * 60)

    # Create test data
    N, G = 50000, 4000
    density = 0.15
    sketch_dim = 512

    print(f"\nCreating test data: {N} spots x {G} genes, density={density}")

    np.random.seed(42)
    data = np.random.exponential(5, int(N * G * density))
    rows = np.random.randint(0, N, len(data))
    cols = np.random.randint(0, G, len(data))
    Y_sparse = csr_matrix((data, (rows, cols)), shape=(N, G))

    sparse_mb = Y_sparse.data.nbytes / 1e6
    dense_mb = N * G * 8 / 1e6
    print(f"Sparse matrix: {sparse_mb:.1f} MB")
    print(f"Dense equivalent: {dense_mb:.1f} MB")

    # Sparse version
    print(f"\n--- Sparse @ Dense Sketching (to {sketch_dim} dims) ---")
    tracemalloc.start()
    t0 = time.time()
    Y_sketch_sparse = demo_sketch_sparse(Y_sparse, sketch_dim)
    t_sparse = time.time() - t0
    current, peak_sparse = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    print(f"Time: {t_sparse:.3f}s, Peak memory: {peak_sparse / 1e6:.1f} MB")
    print(f"Output shape: {Y_sketch_sparse.shape}, dtype: {Y_sketch_sparse.dtype}")

    # Dense version
    print("\n--- Dense @ Dense Sketching ---")
    Y_dense = Y_sparse.toarray()
    tracemalloc.start()
    t0 = time.time()
    Y_sketch_dense = demo_sketch_dense(Y_dense, sketch_dim)
    t_dense = time.time() - t0
    current, peak_dense = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    print(f"Time: {t_dense:.3f}s, Peak memory: {peak_dense / 1e6:.1f} MB")

    # Verify correctness
    max_diff = np.max(np.abs(Y_sketch_sparse - Y_sketch_dense))
    print(f"\n--- Correctness Check ---")
    print(f"Max absolute difference: {max_diff:.2e}")
    print(f"Numerically equivalent: {max_diff < 1e-10}")

    return Y_sketch_sparse, Y_sketch_dense


# =============================================================================
# Demo 4: scipy.sparse Key Operations Reference
# =============================================================================

def run_sparse_ops_reference():
    """
    Reference guide for scipy.sparse operations.

    This demonstrates which operations preserve sparsity and which don't.
    """
    print("\n" + "=" * 60)
    print("Demo 4: scipy.sparse Operations Reference")
    print("=" * 60)

    # Create small test matrix
    Y = csr_matrix([[1, 0, 3], [0, 5, 0], [7, 0, 9]])
    print(f"\nTest matrix Y (3x3):\n{Y.toarray()}")

    operations = [
        ("Y.sum(axis=0)", lambda: Y.sum(axis=0), "Column sums"),
        ("Y.sum(axis=1)", lambda: Y.sum(axis=1), "Row sums"),
        ("Y.mean(axis=0)", lambda: Y.mean(axis=0), "Column means"),
        ("Y.power(2)", lambda: Y.power(2), "Element-wise square"),
        ("Y.multiply(2)", lambda: Y.multiply(2), "Scalar multiply"),
        ("diags([1,2,3]) @ Y", lambda: diags([1, 2, 3]) @ Y, "Diag @ Sparse"),
        ("Y @ np.ones((3,2))", lambda: Y @ np.ones((3, 2)), "Sparse @ Dense"),
        ("Y / Y.sum(axis=1)", lambda: Y / Y.sum(axis=1), "Division (DANGER!)"),
    ]

    print("\n| Operation | Output Type | Sparse? | Notes |")
    print("|-----------|-------------|---------|-------|")

    for name, op, notes in operations:
        try:
            result = op()
            is_sparse = sparse.issparse(result)
            result_type = type(result).__name__
            sparse_str = "Yes" if is_sparse else "NO!"
            print(f"| {name:30} | {result_type:15} | {sparse_str:7} | {notes} |")
        except Exception as e:
            print(f"| {name:30} | ERROR | - | {str(e)[:30]} |")

    print("\n" + "-" * 60)
    print("KEY TAKEAWAYS:")
    print("-" * 60)
    print("1. Y / scalar -> preserves sparsity")
    print("2. Y / array (broadcasting) -> CONVERTS TO DENSE!")
    print("3. diags() @ Y -> preserves sparsity (use for row normalization)")
    print("4. Y.power(n) -> preserves sparsity")
    print("5. Y @ dense -> returns dense (expected for sketching)")
    print("6. Y.data = np.log1p(Y.data) -> in-place, preserves structure")


# =============================================================================
# Main
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description="Sparse optimization demos")
    parser.add_argument(
        "--demo",
        choices=["hvg", "logcpm", "sketch", "ops", "all"],
        default="all",
        help="Which demo to run"
    )
    args = parser.parse_args()

    if args.demo in ["hvg", "all"]:
        run_hvg_demo()

    if args.demo in ["logcpm", "all"]:
        run_logcpm_demo()

    if args.demo in ["sketch", "all"]:
        run_sketch_demo()

    if args.demo in ["ops", "all"]:
        run_sparse_ops_reference()

    print("\n" + "=" * 60)
    print("All demos completed!")
    print("=" * 60)


if __name__ == "__main__":
    main()
