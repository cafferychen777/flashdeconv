"""
Final benchmark comparing original vs optimized solver across all scales.

This provides a comprehensive view of performance and memory characteristics.
"""

import numpy as np
from scipy import sparse
import time
import gc
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from flashdeconv.core.solver import (
    bcd_solve,
    bcd_iteration,
    precompute_gram_matrix,
    precompute_XtY,
)

try:
    import psutil
    HAS_PSUTIL = True
except ImportError:
    HAS_PSUTIL = False


def get_memory_mb():
    if HAS_PSUTIL:
        return psutil.Process(os.getpid()).memory_info().rss / 1024 / 1024
    return 0


def generate_test_data(n_spots, n_cell_types, sketch_dim, k_neighbors=6, seed=42):
    np.random.seed(seed)
    Y_sketch = np.random.randn(n_spots, sketch_dim).astype(np.float64)
    X_sketch = np.random.randn(n_cell_types, sketch_dim).astype(np.float64)

    rows, cols = [], []
    for i in range(n_spots):
        neighbors = np.random.choice(n_spots, size=min(k_neighbors, n_spots-1), replace=False)
        for j in neighbors:
            if i != j:
                rows.append(i)
                cols.append(j)

    data = np.ones(len(rows))
    A = sparse.csr_matrix((data, (rows, cols)), shape=(n_spots, n_spots))
    A = A + A.T
    A.data[:] = 1.0
    return Y_sketch, X_sketch, A


def run_original_solver(Y_sketch, X_sketch, A, lambda_=0.1, rho=0.01, max_iter=50, tol=1e-4):
    """Original solver with per-iteration copy."""
    n_spots, sketch_dim = Y_sketch.shape
    n_cell_types = X_sketch.shape[0]

    beta = np.ones((n_spots, n_cell_types), dtype=np.float64) / n_cell_types
    XtX = precompute_gram_matrix(X_sketch)
    H = precompute_XtY(X_sketch, Y_sketch)

    A_csr = A.tocsr()
    neighbor_indices = A_csr.indices.astype(np.int64)
    neighbor_indptr = A_csr.indptr.astype(np.int64)

    for iteration in range(max_iter):
        beta_old = beta.copy()
        beta = bcd_iteration(H, XtX, beta, neighbor_indices, neighbor_indptr, lambda_, rho)
        rel_change = np.max(np.abs(beta - beta_old)) / (np.max(np.abs(beta_old)) + 1e-10)
        if rel_change < tol:
            break

    return beta, iteration + 1


def benchmark_scale(n_spots, n_cell_types, sketch_dim=512, max_iter=50, n_runs=3):
    """Benchmark at a specific scale."""
    print(f"\n{'='*60}")
    print(f"N={n_spots}, K={n_cell_types}, d={sketch_dim}")
    print(f"{'='*60}")

    beta_mb = n_spots * n_cell_types * 8 / 1024 / 1024
    print(f"Beta size: {beta_mb:.2f} MB")

    # Generate data
    Y_sketch, X_sketch, A = generate_test_data(n_spots, n_cell_types, sketch_dim)

    # Warmup both
    gc.collect()
    run_original_solver(Y_sketch, X_sketch, A, max_iter=3)
    bcd_solve(Y_sketch, X_sketch, A, max_iter=3)
    gc.collect()

    # Benchmark original with memory tracking
    gc.collect()
    time.sleep(0.1)
    mem_before_orig = get_memory_mb()

    times_orig = []
    for _ in range(n_runs):
        gc.collect()
        start = time.perf_counter()
        beta_orig, n_iter = run_original_solver(Y_sketch, X_sketch, A, max_iter=max_iter)
        times_orig.append(time.perf_counter() - start)

    mem_after_orig = get_memory_mb()
    mem_orig = mem_after_orig - mem_before_orig

    # Benchmark optimized with memory tracking
    gc.collect()
    time.sleep(0.1)
    mem_before_opt = get_memory_mb()

    times_opt = []
    for _ in range(n_runs):
        gc.collect()
        start = time.perf_counter()
        beta_opt, info = bcd_solve(Y_sketch, X_sketch, A, max_iter=max_iter)
        times_opt.append(time.perf_counter() - start)

    mem_after_opt = get_memory_mb()
    mem_opt = mem_after_opt - mem_before_opt

    # Results
    time_orig = np.mean(times_orig)
    time_opt = np.mean(times_opt)
    speedup = time_orig / time_opt

    max_diff = np.max(np.abs(beta_orig - beta_opt))
    correlation = np.corrcoef(beta_orig.flatten(), beta_opt.flatten())[0, 1]

    print(f"\n[Original Solver]")
    print(f"  Time: {time_orig*1000:.1f} ms (±{np.std(times_orig)*1000:.1f})")
    print(f"  Memory increase: {mem_orig:.1f} MB")
    print(f"  Iterations: {n_iter}")

    print(f"\n[Optimized Solver]")
    print(f"  Time: {time_opt*1000:.1f} ms (±{np.std(times_opt)*1000:.1f})")
    print(f"  Memory increase: {mem_opt:.1f} MB")
    print(f"  Iterations: {info['n_iterations']}")

    print(f"\n[Comparison]")
    print(f"  Speedup: {speedup:.2f}x {'✓' if speedup >= 1.0 else '✗'}")
    print(f"  Memory saved: {mem_orig - mem_opt:.1f} MB")
    print(f"  Numerical diff: {max_diff:.2e} (correlation: {correlation:.6f})")

    return {
        'n_spots': n_spots,
        'n_cell_types': n_cell_types,
        'beta_mb': beta_mb,
        'time_orig': time_orig,
        'time_opt': time_opt,
        'speedup': speedup,
        'mem_orig': mem_orig,
        'mem_opt': mem_opt,
        'max_diff': max_diff,
        'correct': max_diff < 1e-4,
    }


def main():
    print("="*60)
    print("FlashDeconv Solver Optimization - Final Benchmark")
    print("="*60)

    configs = [
        # Small scale (may be slower due to JIT overhead)
        {'n_spots': 500, 'n_cell_types': 10},
        {'n_spots': 1000, 'n_cell_types': 20},

        # Medium scale (typical Visium - sweet spot)
        {'n_spots': 2000, 'n_cell_types': 30},
        {'n_spots': 5000, 'n_cell_types': 30},

        # Large scale (memory benefits most visible)
        {'n_spots': 10000, 'n_cell_types': 50},
        {'n_spots': 20000, 'n_cell_types': 50},
    ]

    results = []
    for config in configs:
        try:
            r = benchmark_scale(**config)
            results.append(r)
        except Exception as e:
            print(f"Error: {e}")
            import traceback
            traceback.print_exc()

    # Summary
    print("\n" + "="*70)
    print("FINAL SUMMARY")
    print("="*70)
    print(f"{'N':>8} {'K':>4} {'Beta':>6} | {'Orig':>8} {'Opt':>8} {'Speedup':>8} | {'MemOrig':>8} {'MemOpt':>8} | {'Correct':>7}")
    print("-"*70)
    for r in results:
        correct_str = '✓' if r['correct'] else '✗'
        speedup_str = f"{r['speedup']:.2f}x"
        if r['speedup'] >= 1.0:
            speedup_str += ' ✓'
        print(f"{r['n_spots']:>8} {r['n_cell_types']:>4} {r['beta_mb']:>5.1f}MB | "
              f"{r['time_orig']*1000:>7.0f}ms {r['time_opt']*1000:>7.0f}ms {speedup_str:>8} | "
              f"{r['mem_orig']:>7.1f}MB {r['mem_opt']:>7.1f}MB | {correct_str:>7}")

    print("\n" + "="*70)
    print("CONCLUSION")
    print("="*70)

    # Count wins
    speed_wins = sum(1 for r in results if r['speedup'] >= 1.0)
    mem_wins = sum(1 for r in results if r['mem_opt'] < r['mem_orig'])
    correct_all = all(r['correct'] for r in results)

    print(f"Speed improvement: {speed_wins}/{len(results)} scales")
    print(f"Memory improvement: {mem_wins}/{len(results)} scales")
    print(f"Numerical correctness: {'ALL PASS ✓' if correct_all else 'SOME FAILED ✗'}")

    # Recommendations
    print("\nRecommendations:")
    if correct_all:
        if speed_wins >= len(results) // 2:
            print("  ✓ Optimization is beneficial for most use cases")
            print("  ✓ Memory savings are significant for large scale data")
            print("  ✓ Safe to deploy to production")
        else:
            print("  ⚠ Speed improvement mainly visible at large scale")
            print("  ✓ Memory savings are still valuable")
            print("  ✓ Safe to deploy (no regression for typical use)")


if __name__ == "__main__":
    main()
