"""Tests for sketching module."""

import numpy as np
import pytest
from scipy import sparse

from flashdeconv.core.sketching import (
    build_countsketch_matrix,
    build_sparse_rademacher_matrix,
    project_to_sketch,
    sketch_data,
)


class TestCountSketchMatrix:
    """Tests for CountSketch matrix construction."""

    def test_output_shape(self):
        """Test matrix has correct shape."""
        n_genes = 1000
        d = 128

        Omega = build_countsketch_matrix(n_genes, d)

        assert Omega.shape == (n_genes, d)

    def test_sparsity(self):
        """Test matrix is sparse with one nonzero per row."""
        n_genes = 100
        d = 32

        Omega = build_countsketch_matrix(n_genes, d)

        # Each row should have exactly one nonzero
        nnz_per_row = np.diff(Omega.tocsr().indptr)
        assert np.all(nnz_per_row == 1)

    def test_reproducibility(self):
        """Test random state reproducibility."""
        Omega1 = build_countsketch_matrix(100, 32, random_state=42)
        Omega2 = build_countsketch_matrix(100, 32, random_state=42)

        np.testing.assert_array_equal(Omega1.toarray(), Omega2.toarray())

    def test_with_leverage_scores(self):
        """Test with leverage score weighting."""
        n_genes = 100
        leverage = np.random.random(n_genes)
        leverage = leverage / leverage.sum()

        Omega = build_countsketch_matrix(n_genes, 32, leverage_scores=leverage)

        assert Omega.shape == (n_genes, 32)


class TestSparseRademacherMatrix:
    """Tests for sparse Rademacher matrix."""

    def test_output_shape(self):
        """Test correct shape."""
        Omega = build_sparse_rademacher_matrix(1000, 128, sparsity=0.1)

        assert Omega.shape == (1000, 128)

    def test_sparsity_level(self):
        """Test approximate sparsity."""
        n_genes = 1000
        d = 100
        target_sparsity = 0.1

        Omega = build_sparse_rademacher_matrix(n_genes, d, sparsity=target_sparsity)

        actual_density = Omega.nnz / (n_genes * d)
        # Allow some tolerance due to randomness
        assert 0.05 < actual_density < 0.3


class TestProjectToSketch:
    """Tests for sketch projection."""

    def test_output_shapes(self):
        """Test projected matrices have correct shapes."""
        n_spots = 50
        n_genes = 200
        n_cell_types = 10
        d = 64

        Y_tilde = np.random.randn(n_spots, n_genes)
        X_tilde = np.random.randn(n_cell_types, n_genes)
        Omega = build_countsketch_matrix(n_genes, d)

        Y_sketch, X_sketch = project_to_sketch(Y_tilde, X_tilde, Omega)

        assert Y_sketch.shape == (n_spots, d)
        assert X_sketch.shape == (n_cell_types, d)

    def test_linearity(self):
        """Test projection is linear."""
        n_genes = 100
        d = 32

        Omega = build_countsketch_matrix(n_genes, d)

        Y1 = np.random.randn(10, n_genes)
        Y2 = np.random.randn(10, n_genes)
        X = np.random.randn(5, n_genes)

        Y1_sketch, _ = project_to_sketch(Y1, X, Omega)
        Y2_sketch, _ = project_to_sketch(Y2, X, Omega)
        Y_combined_sketch, _ = project_to_sketch(Y1 + Y2, X, Omega)

        np.testing.assert_allclose(Y_combined_sketch, Y1_sketch + Y2_sketch)


class TestSketchData:
    """Tests for full sketching pipeline."""

    def test_pipeline(self):
        """Test full sketching pipeline."""
        n_spots = 100
        n_genes = 500
        n_cell_types = 15
        d = 128

        Y_tilde = np.random.randn(n_spots, n_genes)
        X_tilde = np.random.randn(n_cell_types, n_genes)

        Y_sketch, X_sketch, Omega = sketch_data(Y_tilde, X_tilde, sketch_dim=d)

        assert Y_sketch.shape == (n_spots, d)
        assert X_sketch.shape == (n_cell_types, d)
        assert Omega.shape == (n_genes, d)

    def test_different_methods(self):
        """Test different sketching methods."""
        Y_tilde = np.random.randn(50, 200)
        X_tilde = np.random.randn(10, 200)

        for method in ["countsketch", "rademacher"]:
            Y_sk, X_sk, Omega = sketch_data(
                Y_tilde, X_tilde, sketch_dim=64, method=method
            )
            assert Y_sk.shape == (50, 64)
