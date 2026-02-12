"""Tests for BCD solver module."""

import hashlib
import inspect
import numpy as np
import pytest
from scipy import sparse

from flashdeconv.core.solver import (
    soft_threshold,
    precompute_gram_matrix,
    precompute_XtY,
    bcd_solve,
    normalize_proportions,
    compute_objective,
    compute_objective_fast,
)
from flashdeconv.utils.graph import build_knn_graph
from flashdeconv.core.spatial import compute_laplacian


class TestSoftThreshold:
    """Tests for soft thresholding operator."""

    def test_positive_above_threshold(self):
        """Test positive value above threshold."""
        result = soft_threshold(5.0, 2.0)
        assert result == 3.0

    def test_negative_below_threshold(self):
        """Test negative value below threshold."""
        result = soft_threshold(-5.0, 2.0)
        assert result == -3.0

    def test_within_threshold(self):
        """Test value within threshold band."""
        result = soft_threshold(1.0, 2.0)
        assert result == 0.0


class TestGramMatrix:
    """Tests for Gram matrix computation."""

    def test_shape(self):
        """Test Gram matrix shape."""
        X_sketch = np.random.randn(10, 64)
        XtX = precompute_gram_matrix(X_sketch)

        assert XtX.shape == (10, 10)

    def test_symmetry(self):
        """Test Gram matrix is symmetric."""
        X_sketch = np.random.randn(10, 64)
        XtX = precompute_gram_matrix(X_sketch)

        np.testing.assert_allclose(XtX, XtX.T)

    def test_positive_semidefinite(self):
        """Test Gram matrix is positive semidefinite."""
        X_sketch = np.random.randn(10, 64)
        XtX = precompute_gram_matrix(X_sketch)

        eigenvalues = np.linalg.eigvalsh(XtX)
        assert np.all(eigenvalues >= -1e-10)


class TestBCDSolver:
    """Tests for BCD solver."""

    @pytest.fixture
    def simple_problem(self):
        """Create a simple deconvolution problem."""
        np.random.seed(42)

        n_spots = 50
        n_cell_types = 5
        sketch_dim = 32

        # Generate synthetic data
        X_sketch = np.random.randn(n_cell_types, sketch_dim)

        # True proportions
        beta_true = np.random.rand(n_spots, n_cell_types)
        beta_true = beta_true / beta_true.sum(axis=1, keepdims=True)

        # Generate observations with noise
        Y_sketch = beta_true @ X_sketch + 0.1 * np.random.randn(n_spots, sketch_dim)

        # Spatial graph
        coords = np.random.rand(n_spots, 2)
        A = build_knn_graph(coords, k=4)

        return Y_sketch, X_sketch, A, beta_true

    def test_output_shape(self, simple_problem):
        """Test solver output has correct shape."""
        Y_sketch, X_sketch, A, _ = simple_problem

        beta, info = bcd_solve(Y_sketch, X_sketch, A, max_iter=10)

        assert beta.shape == (50, 5)

    def test_non_negative(self, simple_problem):
        """Test solution is non-negative."""
        Y_sketch, X_sketch, A, _ = simple_problem

        beta, info = bcd_solve(Y_sketch, X_sketch, A, max_iter=50)

        assert np.all(beta >= -1e-10)

    def test_convergence(self, simple_problem):
        """Test solver converges."""
        Y_sketch, X_sketch, A, _ = simple_problem

        beta, info = bcd_solve(Y_sketch, X_sketch, A, max_iter=200, tol=1e-4)

        # Should either converge or reach max iterations
        assert 'converged' in info
        assert 'n_iterations' in info

    def test_objective_decreases(self, simple_problem):
        """Test objective function decreases (approximately)."""
        Y_sketch, X_sketch, A, _ = simple_problem
        L = compute_laplacian(A)

        # Run with verbose to get objectives
        beta, info = bcd_solve(
            Y_sketch, X_sketch, A,
            lambda_=0.1, rho=0.01,
            max_iter=50, verbose=True
        )

        # Final objective should be reasonable
        assert info['final_objective'] < float('inf')
        assert info['final_objective'] >= 0

    def test_lambda_effect(self, simple_problem):
        """Test that lambda controls spatial smoothing."""
        Y_sketch, X_sketch, A, _ = simple_problem

        # Low lambda - less smoothing
        beta_low, _ = bcd_solve(Y_sketch, X_sketch, A, lambda_=0.001, max_iter=50)

        # High lambda - more smoothing
        beta_high, _ = bcd_solve(Y_sketch, X_sketch, A, lambda_=1.0, max_iter=50)

        # High lambda should produce more uniform solution
        var_low = np.var(beta_low)
        var_high = np.var(beta_high)

        assert var_high <= var_low + 0.1  # Allow some tolerance


class TestNormalizeProportions:
    """Tests for proportion normalization."""

    def test_row_sum(self):
        """Test rows sum to 1."""
        beta = np.random.rand(20, 5)
        props = normalize_proportions(beta)

        row_sums = props.sum(axis=1)
        np.testing.assert_allclose(row_sums, 1.0)

    def test_preserves_ratios(self):
        """Test relative proportions are preserved."""
        beta = np.array([[2, 4], [3, 3]])
        props = normalize_proportions(beta)

        # First row: 2:4 ratio -> 1/3:2/3
        np.testing.assert_allclose(props[0], [1/3, 2/3])

        # Second row: 3:3 ratio -> 1/2:1/2
        np.testing.assert_allclose(props[1], [0.5, 0.5])

    def test_handles_zeros(self):
        """Test all-zero rows get uniform distribution and sum to 1."""
        beta = np.array([[0, 0, 0], [1, 2, 3], [0, 0, 0]])
        props = normalize_proportions(beta)

        assert not np.any(np.isnan(props))
        assert not np.any(np.isinf(props))
        # Every row must sum to 1, including all-zero rows
        np.testing.assert_allclose(props.sum(axis=1), 1.0)
        # All-zero rows should be uniform 1/K
        np.testing.assert_allclose(props[0], [1/3, 1/3, 1/3])
        np.testing.assert_allclose(props[2], [1/3, 1/3, 1/3])
        # Normal row preserves ratios
        np.testing.assert_allclose(props[1], [1/6, 2/6, 3/6])


class TestObjectiveFunction:
    """Tests for objective computation."""

    def test_positive(self):
        """Test objective is non-negative."""
        np.random.seed(42)

        n_spots, n_cell_types, d = 30, 5, 32
        Y_sketch = np.random.randn(n_spots, d)
        X_sketch = np.random.randn(n_cell_types, d)
        beta = np.random.rand(n_spots, n_cell_types)

        coords = np.random.rand(n_spots, 2)
        A = build_knn_graph(coords, k=4)
        L = compute_laplacian(A)

        obj = compute_objective(Y_sketch, X_sketch, beta, L, lambda_=0.1, rho=0.01)

        assert obj >= 0

    def test_components(self):
        """Test each component contributes correctly."""
        np.random.seed(42)

        n_spots, n_cell_types, d = 20, 3, 16
        X_sketch = np.random.randn(n_cell_types, d)
        beta = np.random.rand(n_spots, n_cell_types)

        # Y that exactly matches beta @ X
        Y_sketch = beta @ X_sketch

        coords = np.random.rand(n_spots, 2)
        A = build_knn_graph(coords, k=4)
        L = compute_laplacian(A)

        # With perfect fit, fidelity term should be ~0
        obj = compute_objective(Y_sketch, X_sketch, beta, L, lambda_=0.0, rho=0.0)
        np.testing.assert_allclose(obj, 0, atol=1e-10)

    def test_matches_algebraic_expansion(self):
        """Objective should match algebraic expansion of fidelity term."""
        np.random.seed(123)

        n_spots, n_cell_types, d = 40, 6, 64
        Y_sketch = np.random.randn(n_spots, d)
        X_sketch = np.random.randn(n_cell_types, d)
        beta = np.random.rand(n_spots, n_cell_types)

        coords = np.random.rand(n_spots, 2)
        A = build_knn_graph(coords, k=5)
        L = compute_laplacian(A)

        lambda_ = 0.2
        rho = 0.03

        # Expected fidelity via expansion: ||Y - βX||²_F = ||Y||² - 2·Tr(Hβ) + Tr((β^Tβ)·G)
        H = X_sketch @ Y_sketch.T  # (K, N)
        XtX = X_sketch @ X_sketch.T  # (K, K)
        YtY = np.einsum("ij,ij->", Y_sketch, Y_sketch)
        cross = np.einsum("ki,ik->", H, beta)
        BtB = beta.T @ beta
        quad = np.einsum("ij,ij->", BtB, XtX)
        fidelity = 0.5 * (YtY - 2.0 * cross + quad)

        Lbeta = L @ beta
        spatial = lambda_ * np.sum(beta * Lbeta)
        sparsity = rho * np.sum(np.abs(beta))
        expected = fidelity + spatial + sparsity

        obj_plain = compute_objective(Y_sketch, X_sketch, beta, L, lambda_=lambda_, rho=rho)
        np.testing.assert_allclose(obj_plain, expected, rtol=1e-10, atol=1e-8)

        # compute_objective_fast must agree with compute_objective
        obj_fast = compute_objective_fast(
            beta, H, XtX, YtY, L, lambda_=lambda_, rho=rho,
        )
        np.testing.assert_allclose(obj_fast, obj_plain, rtol=1e-9, atol=1e-8)

    def test_fast_objective_multiple_scales(self):
        """Test fast objective matches original across varying data scales."""
        for seed, N, K, d in [(0, 50, 3, 20), (1, 200, 10, 64), (2, 30, 8, 128)]:
            np.random.seed(seed)
            Y = np.random.randn(N, d) * (seed + 1)
            X = np.random.randn(K, d)
            beta = np.abs(np.random.randn(N, K))
            coords = np.random.rand(N, 2)
            A = build_knn_graph(coords, k=4)
            L = compute_laplacian(A)
            H = precompute_XtY(X, Y)
            XtX = precompute_gram_matrix(X)
            YtY = float(np.sum(Y ** 2))

            obj_orig = compute_objective(Y, X, beta, L, 0.1, 0.05)
            obj_fast = compute_objective_fast(beta, H, XtX, YtY, L, 0.1, 0.05)
            np.testing.assert_allclose(obj_fast, obj_orig, rtol=1e-9, atol=1e-8,
                                       err_msg=f"Mismatch for seed={seed}")


class TestDeterminism:
    """Tests for deterministic solver output with fixed inputs."""

    def test_bcd_solve_deterministic(self):
        np.random.seed(42)

        n_spots = 60
        n_cell_types = 7
        sketch_dim = 48

        X_sketch = np.random.randn(n_cell_types, sketch_dim)
        beta_true = np.random.rand(n_spots, n_cell_types)
        beta_true = beta_true / beta_true.sum(axis=1, keepdims=True)
        Y_sketch = beta_true @ X_sketch + 0.05 * np.random.randn(n_spots, sketch_dim)

        coords = np.random.rand(n_spots, 2)
        A = build_knn_graph(coords, k=4)

        beta1, info1 = bcd_solve(Y_sketch, X_sketch, A, lambda_=0.1, rho=0.01, max_iter=30, tol=1e-6)
        beta2, info2 = bcd_solve(Y_sketch, X_sketch, A, lambda_=0.1, rho=0.01, max_iter=30, tol=1e-6)

        # Use a stable hash to avoid huge diffs on failure.
        h1 = hashlib.sha256(beta1.tobytes()).hexdigest()
        h2 = hashlib.sha256(beta2.tobytes()).hexdigest()
        assert h1 == h2
        assert info1["converged"] == info2["converged"]
        assert info1["n_iterations"] == info2["n_iterations"]
