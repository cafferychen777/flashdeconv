"""Tests for preprocessing module."""

import numpy as np
import pytest
from scipy import sparse

from flashdeconv.core.preprocessing import (
    compute_size_factors,
    compute_gene_means,
    estimate_overdispersion,
    pearson_residuals,
    preprocess_spatial_data,
)


class TestSizeFactors:
    """Tests for size factor computation."""

    def test_dense_input(self):
        """Test with dense matrix."""
        Y = np.array([[10, 20, 30], [5, 10, 15], [100, 200, 300]])
        sf = compute_size_factors(Y)

        assert sf.shape == (3,)
        np.testing.assert_array_equal(sf, [60, 30, 600])

    def test_sparse_input(self):
        """Test with sparse matrix."""
        Y = sparse.csr_matrix([[10, 20, 30], [5, 10, 15], [100, 200, 300]])
        sf = compute_size_factors(Y)

        assert sf.shape == (3,)
        np.testing.assert_array_equal(sf, [60, 30, 600])

    def test_zero_row(self):
        """Test handling of zero-count rows."""
        Y = np.array([[10, 20], [0, 0], [5, 5]])
        sf = compute_size_factors(Y)

        # Zero rows should get sf=1 (minimum)
        assert sf[1] == 1.0


class TestGeneMeans:
    """Tests for gene mean computation."""

    def test_basic(self):
        """Test basic gene mean computation."""
        Y = np.array([[10, 20], [10, 20], [10, 20]])
        sf = compute_size_factors(Y)
        gm = compute_gene_means(Y, sf)

        assert gm.shape == (2,)
        # All spots same, so gene means = gene sums / total counts
        np.testing.assert_allclose(gm, [30 / 90, 60 / 90])


class TestOverdispersion:
    """Tests for overdispersion estimation."""

    def test_estimate_range(self):
        """Test that theta is in valid range."""
        np.random.seed(42)
        Y = np.random.negative_binomial(n=5, p=0.5, size=(100, 50))

        theta = estimate_overdispersion(Y)

        assert 0.01 <= theta <= 100.0

    def test_poisson_like(self):
        """Test with Poisson-like data (low overdispersion)."""
        np.random.seed(42)
        # Poisson has variance = mean, so theta should be high
        Y = np.random.poisson(lam=10, size=(100, 50))

        theta = estimate_overdispersion(Y)

        # For Poisson, theta should be large (>10)
        assert theta > 5.0


class TestPearsonResiduals:
    """Tests for Pearson residuals transformation."""

    def test_output_shape(self):
        """Test output shape matches input."""
        Y = np.random.poisson(10, size=(50, 100))
        residuals = pearson_residuals(Y)

        assert residuals.shape == Y.shape

    def test_dense_output(self):
        """Test that output is always dense."""
        Y = sparse.random(50, 100, density=0.3, format='csr')
        Y.data[:] = np.abs(Y.data) * 10

        residuals = pearson_residuals(Y)

        assert isinstance(residuals, np.ndarray)
        assert not sparse.issparse(residuals)

    def test_clipping(self):
        """Test that extreme values are clipped."""
        Y = np.zeros((10, 10))
        Y[0, 0] = 10000  # Extreme value

        residuals = pearson_residuals(Y, clip_value=10.0)

        assert np.abs(residuals).max() <= 10.0


class TestPreprocessPipeline:
    """Tests for full preprocessing pipeline."""

    def test_pipeline_output(self):
        """Test pipeline returns all expected outputs."""
        np.random.seed(42)
        Y = np.random.negative_binomial(n=5, p=0.3, size=(100, 200))

        Y_tilde, sf, theta, gm = preprocess_spatial_data(Y)

        assert Y_tilde.shape == Y.shape
        assert sf.shape == (100,)
        assert isinstance(theta, float)
        assert gm.shape == (200,)

    def test_sparse_input(self):
        """Test pipeline with sparse input."""
        np.random.seed(42)
        Y_dense = np.random.negative_binomial(n=5, p=0.3, size=(100, 200))
        Y_sparse = sparse.csr_matrix(Y_dense)

        Y_tilde, sf, theta, gm = preprocess_spatial_data(Y_sparse)

        assert Y_tilde.shape == Y_dense.shape
        assert not sparse.issparse(Y_tilde)
