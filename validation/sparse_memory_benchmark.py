#!/usr/bin/env python3
"""
Memory Benchmark for FlashDeconv Sparse Optimization

This script benchmarks memory usage at different scales to verify
that the sparse optimization is working correctly.

Usage:
    python sparse_memory_benchmark.py --scale small   # Quick test
    python sparse_memory_benchmark.py --scale medium  # Xenium-scale
    python sparse_memory_benchmark.py --scale large   # Atlas-scale (needs ~16GB RAM)
"""

import numpy as np
from scipy import sparse
from scipy.sparse import csr_matrix, diags
import time
import argparse
import tracemalloc
import gc
from dataclasses import dataclass
from typing import Dict, Any


@dataclass
class BenchmarkResult:
    """Results from a benchmark run."""
    stage: str
    time_sec: float
    peak_mb: float
    current_mb: float
    notes: str = ""


class MemoryTracker:
    """Context manager for tracking memory usage."""

    def __init__(self, stage_name: str):
        self.stage_name = stage_name
        self.start_time = None
        self.result = None

    def __enter__(self):
        gc.collect()
        tracemalloc.start()
        self.start_time = time.time()
        return self

    def __exit__(self, *args):
        elapsed = time.time() - self.start_time
        current, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        self.result = BenchmarkResult(
            stage=self.stage_name,
            time_sec=elapsed,
            peak_mb=peak / 1e6,
            current_mb=current / 1e6
        )


# =============================================================================
# Simulated Pipeline Stages (Current Implementation)
# =============================================================================

def current_select_hvg(Y_sparse: csr_matrix, n_top: int = 2000) -> np.ndarray:
    """Current implementation: converts to dense immediately."""
    # This is what the current code does - MEMORY BOMB!
    Y_dense = Y_sparse.toarray()

    # Rest of HVG selection...
    total_counts = np.sum(Y_dense, axis=1, keepdims=True)
    total_counts = np.maximum(total_counts, 1)
    Y_norm = Y_dense / total_counts * 10000
    Y_log = np.log1p(Y_norm)

    gene_means = np.mean(Y_log, axis=0)
    gene_vars = np.var(Y_log, axis=0, ddof=1)
    mean_safe = np.maximum(gene_means, 1e-12)
    dispersion = gene_vars / mean_safe

    hvg_idx = np.argsort(dispersion)[-n_top:]
    return np.sort(hvg_idx)


def current_preprocess(Y_subset: np.ndarray) -> np.ndarray:
    """Current implementation: expects dense input."""
    lib_size = Y_subset.sum(axis=1, keepdims=True)
    lib_size = np.maximum(lib_size, 1.0)
    Y_cpm = Y_subset / lib_size * 1e4
    return np.log1p(Y_cpm)


# =============================================================================
# Optimized Pipeline Stages (Proposed Implementation)
# =============================================================================

def optimized_select_hvg(Y_sparse: csr_matrix, n_top: int = 2000) -> np.ndarray:
    """Optimized: stays sparse throughout."""
    N, G = Y_sparse.shape

    # Normalize by library size using diagonal matrix
    lib_size = np.array(Y_sparse.sum(axis=1)).flatten()
    lib_size = np.maximum(lib_size, 1.0)
    D_inv = diags(1.0 / lib_size)
    Y_norm = D_inv @ Y_sparse  # Still sparse!

    # Compute statistics on sparse matrix
    mean = np.array(Y_norm.mean(axis=0)).flatten()
    mean_sq = np.array(Y_norm.power(2).mean(axis=0)).flatten()
    var = np.maximum(mean_sq - mean ** 2, 0)

    mean_safe = np.maximum(mean, 1e-12)
    dispersion = var / mean_safe

    hvg_idx = np.argsort(dispersion)[-n_top:]
    return np.sort(hvg_idx)


def optimized_preprocess(Y_sparse: csr_matrix) -> csr_matrix:
    """Optimized: keeps matrix sparse."""
    lib_size = np.array(Y_sparse.sum(axis=1)).flatten()
    lib_size[lib_size == 0] = 1.0

    D = diags(1e4 / lib_size)
    Y_cpm = D @ Y_sparse

    Y_log = Y_cpm.copy()
    Y_log.data = np.log1p(Y_log.data)
    return Y_log


def sketch_sparse(Y_sparse: csr_matrix, sketch_dim: int = 512) -> np.ndarray:
    """Sketch from sparse matrix."""
    G = Y_sparse.shape[1]
    np.random.seed(42)
    Omega = np.random.randn(G, sketch_dim) / np.sqrt(sketch_dim)
    return Y_sparse @ Omega  # Sparse @ Dense = Dense


def sketch_dense(Y_dense: np.ndarray, sketch_dim: int = 512) -> np.ndarray:
    """Sketch from dense matrix."""
    G = Y_dense.shape[1]
    np.random.seed(42)
    Omega = np.random.randn(G, sketch_dim) / np.sqrt(sketch_dim)
    return Y_dense @ Omega


# =============================================================================
# Full Pipeline Benchmarks
# =============================================================================

def benchmark_current_pipeline(N: int, G: int, G_selected: int,
                               density: float) -> Dict[str, BenchmarkResult]:
    """Benchmark the current (non-optimized) pipeline."""
    results = {}

    print(f"\n{'='*60}")
    print(f"CURRENT PIPELINE (N={N:,}, G={G:,})")
    print(f"{'='*60}")

    # Create sparse input
    print("\n[1/5] Creating sparse input data...")
    np.random.seed(42)
    nnz = int(N * G * density)
    data = np.random.exponential(5, nnz).astype(np.float32)
    rows = np.random.randint(0, N, nnz)
    cols = np.random.randint(0, G, nnz)
    Y_sparse = csr_matrix((data, (rows, cols)), shape=(N, G))

    input_mb = Y_sparse.data.nbytes / 1e6
    print(f"    Sparse input size: {input_mb:.1f} MB")
    print(f"    Dense equivalent would be: {N * G * 4 / 1e6:.1f} MB")

    # Stage 1: HVG selection (THE BOTTLENECK!)
    print("\n[2/5] HVG Selection (current: converts to dense)...")
    try:
        with MemoryTracker("hvg_selection") as tracker:
            hvg_idx = current_select_hvg(Y_sparse, n_top=G_selected)
        results["hvg"] = tracker.result
        print(f"    Time: {tracker.result.time_sec:.2f}s")
        print(f"    Peak memory: {tracker.result.peak_mb:.1f} MB")
    except MemoryError:
        print("    ❌ OUT OF MEMORY!")
        results["hvg"] = BenchmarkResult("hvg_selection", -1, -1, -1, "OOM")
        return results

    # Stage 2: Subset genes
    print("\n[3/5] Gene subsetting (current: toarray())...")
    with MemoryTracker("subset") as tracker:
        Y_subset = Y_sparse[:, hvg_idx].toarray()
    results["subset"] = tracker.result
    print(f"    Output shape: {Y_subset.shape}")
    print(f"    Peak memory: {tracker.result.peak_mb:.1f} MB")

    # Stage 3: Preprocessing
    print("\n[4/5] Log-CPM preprocessing...")
    with MemoryTracker("preprocess") as tracker:
        Y_tilde = current_preprocess(Y_subset)
    results["preprocess"] = tracker.result
    print(f"    Peak memory: {tracker.result.peak_mb:.1f} MB")

    # Stage 4: Sketching
    print("\n[5/5] Sketching to 512 dims...")
    with MemoryTracker("sketch") as tracker:
        Y_sketch = sketch_dense(Y_tilde, sketch_dim=512)
    results["sketch"] = tracker.result
    print(f"    Output shape: {Y_sketch.shape}")
    print(f"    Peak memory: {tracker.result.peak_mb:.1f} MB")

    return results


def benchmark_optimized_pipeline(N: int, G: int, G_selected: int,
                                  density: float) -> Dict[str, BenchmarkResult]:
    """Benchmark the optimized pipeline."""
    results = {}

    print(f"\n{'='*60}")
    print(f"OPTIMIZED PIPELINE (N={N:,}, G={G:,})")
    print(f"{'='*60}")

    # Create sparse input
    print("\n[1/5] Creating sparse input data...")
    np.random.seed(42)
    nnz = int(N * G * density)
    data = np.random.exponential(5, nnz).astype(np.float32)
    rows = np.random.randint(0, N, nnz)
    cols = np.random.randint(0, G, nnz)
    Y_sparse = csr_matrix((data, (rows, cols)), shape=(N, G))

    input_mb = Y_sparse.data.nbytes / 1e6
    print(f"    Sparse input size: {input_mb:.1f} MB")

    # Stage 1: HVG selection (OPTIMIZED!)
    print("\n[2/5] HVG Selection (optimized: stays sparse)...")
    with MemoryTracker("hvg_selection") as tracker:
        hvg_idx = optimized_select_hvg(Y_sparse, n_top=G_selected)
    results["hvg"] = tracker.result
    print(f"    Time: {tracker.result.time_sec:.2f}s")
    print(f"    Peak memory: {tracker.result.peak_mb:.1f} MB")

    # Stage 2: Subset genes (OPTIMIZED: keep sparse!)
    print("\n[3/5] Gene subsetting (optimized: stays sparse)...")
    with MemoryTracker("subset") as tracker:
        Y_subset = Y_sparse[:, hvg_idx]  # Still sparse!
    results["subset"] = tracker.result
    print(f"    Output type: {type(Y_subset).__name__}")
    print(f"    Peak memory: {tracker.result.peak_mb:.1f} MB")

    # Stage 3: Preprocessing (OPTIMIZED!)
    print("\n[4/5] Log-CPM preprocessing (optimized: stays sparse)...")
    with MemoryTracker("preprocess") as tracker:
        Y_tilde = optimized_preprocess(Y_subset)
    results["preprocess"] = tracker.result
    print(f"    Output type: {type(Y_tilde).__name__}")
    print(f"    Peak memory: {tracker.result.peak_mb:.1f} MB")

    # Stage 4: Sketching
    print("\n[5/5] Sketching to 512 dims (sparse @ dense)...")
    with MemoryTracker("sketch") as tracker:
        Y_sketch = sketch_sparse(Y_tilde, sketch_dim=512)
    results["sketch"] = tracker.result
    print(f"    Output shape: {Y_sketch.shape}")
    print(f"    Peak memory: {tracker.result.peak_mb:.1f} MB")

    return results


def print_comparison(current: Dict, optimized: Dict):
    """Print comparison table."""
    print(f"\n{'='*60}")
    print("COMPARISON SUMMARY")
    print(f"{'='*60}")

    print("\n| Stage | Current (MB) | Optimized (MB) | Reduction |")
    print("|-------|--------------|----------------|-----------|")

    total_current = 0
    total_optimized = 0

    for stage in ["hvg", "subset", "preprocess", "sketch"]:
        if stage in current and stage in optimized:
            c = current[stage].peak_mb
            o = optimized[stage].peak_mb

            if c > 0 and o > 0:
                reduction = c / o
                print(f"| {stage:12} | {c:12.1f} | {o:14.1f} | {reduction:7.1f}x |")
                total_current += c
                total_optimized += o
            elif c < 0:
                print(f"| {stage:12} | {'OOM':>12} | {o:14.1f} | {'∞':>9} |")

    if total_current > 0 and total_optimized > 0:
        print(f"| {'TOTAL':12} | {total_current:12.1f} | {total_optimized:14.1f} | {total_current/total_optimized:7.1f}x |")


# =============================================================================
# Main
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description="Memory benchmark for sparse optimization")
    parser.add_argument(
        "--scale",
        choices=["tiny", "small", "medium", "large", "atlas"],
        default="small",
        help="Scale of benchmark"
    )
    parser.add_argument(
        "--optimized-only",
        action="store_true",
        help="Only run optimized pipeline (for large scale tests)"
    )
    args = parser.parse_args()

    # Define test scales
    scales = {
        "tiny":   {"N": 1000,    "G": 5000,  "G_selected": 1000, "density": 0.15},
        "small":  {"N": 5000,    "G": 20000, "G_selected": 2000, "density": 0.10},
        "medium": {"N": 50000,   "G": 20000, "G_selected": 3000, "density": 0.10},
        "large":  {"N": 200000,  "G": 20000, "G_selected": 4000, "density": 0.08},
        "atlas":  {"N": 1000000, "G": 20000, "G_selected": 4000, "density": 0.05},
    }

    params = scales[args.scale]
    print(f"\n{'#'*60}")
    print(f"# FlashDeconv Sparse Optimization Benchmark")
    print(f"# Scale: {args.scale}")
    print(f"# N={params['N']:,}, G={params['G']:,}, density={params['density']}")
    print(f"{'#'*60}")

    # Estimate memory requirements
    full_dense_gb = params['N'] * params['G'] * 8 / 1e9
    subset_dense_gb = params['N'] * params['G_selected'] * 8 / 1e9
    print(f"\nExpected memory (current pipeline):")
    print(f"  - HVG stage (full dense): {full_dense_gb:.1f} GB")
    print(f"  - After gene selection: {subset_dense_gb:.1f} GB")

    # Run benchmarks
    if not args.optimized_only:
        if args.scale in ["large", "atlas"]:
            print(f"\n⚠️  Warning: Current pipeline will likely OOM at {args.scale} scale!")
            print("    Use --optimized-only to skip current pipeline benchmark.")

        try:
            current_results = benchmark_current_pipeline(**params)
        except MemoryError:
            print("\n❌ Current pipeline: OUT OF MEMORY")
            current_results = {}
    else:
        current_results = {}

    optimized_results = benchmark_optimized_pipeline(**params)

    if current_results:
        print_comparison(current_results, optimized_results)

    print("\n" + "=" * 60)
    print("Benchmark complete!")
    print("=" * 60)


if __name__ == "__main__":
    main()
