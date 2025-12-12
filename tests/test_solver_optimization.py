"""
Integration tests for solver optimization.

Verifies that the optimized double-buffer implementation produces
numerically equivalent results to the original implementation.
"""

import numpy as np
from scipy import sparse
import pytest
import sys
import os

# Add project root to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from flashdeconv.core.solver import (
    bcd_solve,
    bcd_iteration,
    _bcd_iteration_inplace,
    precompute_gram_matrix,
    precompute_XtY,
    normalize_proportions,
)


def generate_test_data(n_spots, n_cell_types, sketch_dim, k_neighbors=6, seed=42):
    """Generate synthetic test data."""
    np.random.seed(seed)

    Y_sketch = np.random.randn(n_spots, sketch_dim).astype(np.float64)
    X_sketch = np.random.randn(n_cell_types, sketch_dim).astype(np.float64)

    # Random spatial adjacency (k-NN style)
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
    """
    Run the original (non-optimized) solver implementation.

    This manually implements the old behavior with beta.copy() per iteration.
    """
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

        # Use the legacy bcd_iteration function (which does beta.copy())
        beta = bcd_iteration(
            H, XtX, beta,
            neighbor_indices, neighbor_indptr,
            lambda_, rho
        )

        rel_change = np.max(np.abs(beta - beta_old)) / (np.max(np.abs(beta_old)) + 1e-10)
        if rel_change < tol:
            break

    return beta, iteration + 1


class TestSolverOptimization:
    """Test suite for solver optimization correctness."""

    @pytest.mark.parametrize("n_spots,n_cell_types", [
        (100, 5),
        (500, 10),
        (1000, 20),
        (2000, 30),
    ])
    def test_numerical_equivalence(self, n_spots, n_cell_types):
        """Test that optimized solver produces same results as original."""
        sketch_dim = 256
        Y_sketch, X_sketch, A = generate_test_data(n_spots, n_cell_types, sketch_dim)

        # Run original
        beta_orig, n_iter_orig = run_original_solver(
            Y_sketch, X_sketch, A, max_iter=30
        )

        # Run optimized (current bcd_solve)
        beta_opt, info_opt = bcd_solve(
            Y_sketch, X_sketch, A, max_iter=30
        )

        # Check numerical equivalence
        max_diff = np.max(np.abs(beta_orig - beta_opt))
        mean_diff = np.mean(np.abs(beta_orig - beta_opt))
        correlation = np.corrcoef(beta_orig.flatten(), beta_opt.flatten())[0, 1]

        print(f"\nN={n_spots}, K={n_cell_types}:")
        print(f"  Max diff: {max_diff:.2e}")
        print(f"  Mean diff: {mean_diff:.2e}")
        print(f"  Correlation: {correlation:.6f}")
        print(f"  Iterations: orig={n_iter_orig}, opt={info_opt['n_iterations']}")

        # Assertions
        assert max_diff < 1e-4, f"Max diff {max_diff} exceeds threshold"
        assert correlation > 0.9999, f"Correlation {correlation} too low"

    def test_convergence_behavior(self):
        """Test that convergence behavior is preserved."""
        Y_sketch, X_sketch, A = generate_test_data(1000, 20, 256)

        # Run with verbose to check convergence
        beta, info = bcd_solve(
            Y_sketch, X_sketch, A,
            max_iter=100,
            tol=1e-6,
            verbose=False
        )

        assert info['converged'] or info['n_iterations'] == 100
        assert info['final_change'] >= 0
        assert info['final_objective'] >= 0

    def test_iteration_function_equivalence(self):
        """Test that _bcd_iteration_inplace produces same result as bcd_iteration."""
        n_spots, n_cell_types, sketch_dim = 500, 15, 128
        Y_sketch, X_sketch, A = generate_test_data(n_spots, n_cell_types, sketch_dim)

        XtX = precompute_gram_matrix(X_sketch)
        H = precompute_XtY(X_sketch, Y_sketch)

        A_csr = A.tocsr()
        neighbor_indices = A_csr.indices.astype(np.int64)
        neighbor_indptr = A_csr.indptr.astype(np.int64)

        # Initial beta
        np.random.seed(123)
        beta_init = np.random.rand(n_spots, n_cell_types).astype(np.float64)
        beta_init = beta_init / beta_init.sum(axis=1, keepdims=True)

        # Run original bcd_iteration
        beta_orig = bcd_iteration(
            H, XtX, beta_init.copy(),
            neighbor_indices, neighbor_indptr,
            0.1, 0.01
        )

        # Run optimized _bcd_iteration_inplace
        beta_src = beta_init.copy()
        beta_dst = np.empty_like(beta_src)
        _bcd_iteration_inplace(
            H, XtX, beta_src, beta_dst,
            neighbor_indices, neighbor_indptr,
            0.1, 0.01
        )

        max_diff = np.max(np.abs(beta_orig - beta_dst))
        print(f"\nSingle iteration max diff: {max_diff:.2e}")

        assert max_diff < 1e-10, f"Single iteration diff {max_diff} too large"

    def test_normalize_proportions(self):
        """Test that normalization produces valid proportions."""
        beta = np.array([[1.0, 2.0, 3.0], [0.5, 0.5, 0.0]])
        props = normalize_proportions(beta)

        # Check rows sum to 1
        row_sums = props.sum(axis=1)
        np.testing.assert_allclose(row_sums, 1.0, rtol=1e-10)

        # Check proportions are non-negative
        assert np.all(props >= 0)

    def test_edge_cases(self):
        """Test edge cases."""
        # Very small problem
        Y_sketch, X_sketch, A = generate_test_data(10, 3, 32)
        beta, info = bcd_solve(Y_sketch, X_sketch, A, max_iter=10)
        assert beta.shape == (10, 3)

        # Single cell type
        Y_sketch, X_sketch, A = generate_test_data(50, 1, 32)
        beta, info = bcd_solve(Y_sketch, X_sketch, A, max_iter=10)
        assert beta.shape == (50, 1)

    @pytest.mark.parametrize("lambda_,rho", [
        (0.0, 0.0),
        (0.1, 0.01),
        (1.0, 0.1),
        (10.0, 0.5),
    ])
    def test_regularization_parameters(self, lambda_, rho):
        """Test with various regularization parameters."""
        Y_sketch, X_sketch, A = generate_test_data(200, 10, 64)

        beta_orig, _ = run_original_solver(
            Y_sketch, X_sketch, A,
            lambda_=lambda_, rho=rho, max_iter=20
        )

        beta_opt, _ = bcd_solve(
            Y_sketch, X_sketch, A,
            lambda_=lambda_, rho=rho, max_iter=20
        )

        max_diff = np.max(np.abs(beta_orig - beta_opt))
        print(f"\nlambda={lambda_}, rho={rho}: max_diff={max_diff:.2e}")

        assert max_diff < 1e-4


class TestPerformanceImprovement:
    """Test that optimized version is actually faster."""

    def test_performance_medium_scale(self):
        """Test performance on medium scale data."""
        import time

        n_spots, n_cell_types = 2000, 25
        Y_sketch, X_sketch, A = generate_test_data(n_spots, n_cell_types, 256)

        # Warmup
        bcd_solve(Y_sketch, X_sketch, A, max_iter=5)
        run_original_solver(Y_sketch, X_sketch, A, max_iter=5)

        # Time original
        start = time.perf_counter()
        for _ in range(3):
            run_original_solver(Y_sketch, X_sketch, A, max_iter=30)
        time_orig = (time.perf_counter() - start) / 3

        # Time optimized
        start = time.perf_counter()
        for _ in range(3):
            bcd_solve(Y_sketch, X_sketch, A, max_iter=30)
        time_opt = (time.perf_counter() - start) / 3

        speedup = time_orig / time_opt
        print(f"\nN={n_spots}, K={n_cell_types}:")
        print(f"  Original: {time_orig*1000:.1f} ms")
        print(f"  Optimized: {time_opt*1000:.1f} ms")
        print(f"  Speedup: {speedup:.2f}x")

        # We expect at least some speedup (or at worst, not much slower)
        assert speedup > 0.8, f"Optimized version is too slow: {speedup:.2f}x"


def run_all_tests():
    """Run all tests and print summary."""
    print("="*60)
    print("FlashDeconv Solver Optimization Tests")
    print("="*60)

    test_suite = TestSolverOptimization()
    perf_suite = TestPerformanceImprovement()

    # Numerical equivalence tests
    print("\n[1] Numerical Equivalence Tests")
    print("-"*40)
    for n_spots, n_cell_types in [(100, 5), (500, 10), (1000, 20), (2000, 30)]:
        try:
            test_suite.test_numerical_equivalence(n_spots, n_cell_types)
            print(f"  PASSED: N={n_spots}, K={n_cell_types}")
        except AssertionError as e:
            print(f"  FAILED: N={n_spots}, K={n_cell_types}: {e}")

    # Single iteration test
    print("\n[2] Single Iteration Equivalence")
    print("-"*40)
    try:
        test_suite.test_iteration_function_equivalence()
        print("  PASSED")
    except AssertionError as e:
        print(f"  FAILED: {e}")

    # Convergence test
    print("\n[3] Convergence Behavior")
    print("-"*40)
    try:
        test_suite.test_convergence_behavior()
        print("  PASSED")
    except AssertionError as e:
        print(f"  FAILED: {e}")

    # Edge cases
    print("\n[4] Edge Cases")
    print("-"*40)
    try:
        test_suite.test_edge_cases()
        print("  PASSED")
    except AssertionError as e:
        print(f"  FAILED: {e}")

    # Regularization parameters
    print("\n[5] Regularization Parameters")
    print("-"*40)
    for lambda_, rho in [(0.0, 0.0), (0.1, 0.01), (1.0, 0.1), (10.0, 0.5)]:
        try:
            test_suite.test_regularization_parameters(lambda_, rho)
            print(f"  PASSED: lambda={lambda_}, rho={rho}")
        except AssertionError as e:
            print(f"  FAILED: lambda={lambda_}, rho={rho}: {e}")

    # Performance test
    print("\n[6] Performance Test")
    print("-"*40)
    try:
        perf_suite.test_performance_medium_scale()
        print("  PASSED")
    except AssertionError as e:
        print(f"  FAILED: {e}")

    print("\n" + "="*60)
    print("All tests completed!")
    print("="*60)


if __name__ == "__main__":
    run_all_tests()
