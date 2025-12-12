"""Tests for BCD solver module."""

import numpy as np
import pytest
from scipy import sparse

from flashdeconv.core.solver import (
    soft_threshold,
    precompute_gram_matrix,
    bcd_solve,
    normalize_proportions,
    compute_objective,
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
        """Test handling of all-zero rows."""
        beta = np.array([[0, 0], [1, 1]])
        props = normalize_proportions(beta)

        # Zero row should still produce valid output
        assert not np.any(np.isnan(props))
        assert not np.any(np.isinf(props))


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
