"""Integration tests for FlashDeconv."""

import numpy as np
import pytest
from scipy import sparse

from flashdeconv import FlashDeconv


def generate_synthetic_data(
    n_spots: int = 100,
    n_genes: int = 500,
    n_cell_types: int = 5,
    noise_level: float = 0.1,
    random_state: int = 42,
):
    """
    Generate synthetic spatial transcriptomics data for testing.

    Parameters
    ----------
    n_spots : int
        Number of spatial spots.
    n_genes : int
        Number of genes.
    n_cell_types : int
        Number of cell types.
    noise_level : float
        Noise level for Poisson sampling.
    random_state : int
        Random seed.

    Returns
    -------
    Y : ndarray
        Spatial count matrix.
    X : ndarray
        Reference signature matrix.
    coords : ndarray
        Spatial coordinates.
    beta_true : ndarray
        True cell type proportions.
    """
    np.random.seed(random_state)

    # Generate reference signatures (log-normal expression)
    X = np.exp(np.random.randn(n_cell_types, n_genes) * 0.5 + 1)

    # Add cell-type specific markers
    for k in range(n_cell_types):
        marker_genes = np.random.choice(n_genes, size=20, replace=False)
        X[k, marker_genes] *= 5  # Upregulate markers

    # Generate spatial coordinates (grid-like)
    side = int(np.ceil(np.sqrt(n_spots)))
    x = np.tile(np.arange(side), side)[:n_spots]
    y = np.repeat(np.arange(side), side)[:n_spots]
    coords = np.column_stack([x, y]).astype(float)

    # Add small jitter
    coords += np.random.randn(n_spots, 2) * 0.1

    # Generate true proportions with spatial smoothness
    beta_true = np.zeros((n_spots, n_cell_types))
    for k in range(n_cell_types):
        # Create spatial pattern
        center = np.random.rand(2) * side
        distances = np.sqrt(np.sum((coords - center) ** 2, axis=1))
        beta_true[:, k] = np.exp(-distances / (side / 2))

    # Normalize to sum to 1
    beta_true = beta_true / beta_true.sum(axis=1, keepdims=True)

    # Generate expected expression
    expected = beta_true @ X

    # Scale by sequencing depth
    depth = np.random.gamma(shape=5, scale=1000, size=n_spots)
    expected = expected * depth[:, np.newaxis]

    # Sample counts (Negative Binomial approximation via Poisson)
    Y = np.random.poisson(expected * (1 + noise_level * np.random.rand(*expected.shape)))

    return Y, X, coords, beta_true


class TestFlashDeconvIntegration:
    """Integration tests for FlashDeconv class."""

    @pytest.fixture
    def synthetic_data(self):
        """Generate synthetic test data."""
        return generate_synthetic_data(
            n_spots=100,
            n_genes=500,
            n_cell_types=5,
            noise_level=0.1,
            random_state=42,
        )

    def test_fit_basic(self, synthetic_data):
        """Test basic fitting works."""
        Y, X, coords, _ = synthetic_data

        model = FlashDeconv(
            sketch_dim=64,
            lambda_spatial=0.1,
            max_iter=20,
        )

        model.fit(Y, X, coords)

        assert model._fitted
        assert model.beta_ is not None
        assert model.proportions_ is not None

    def test_fit_transform(self, synthetic_data):
        """Test fit_transform returns proportions."""
        Y, X, coords, _ = synthetic_data

        model = FlashDeconv(sketch_dim=64, max_iter=20)
        proportions = model.fit_transform(Y, X, coords)

        assert proportions.shape == (100, 5)
        np.testing.assert_allclose(proportions.sum(axis=1), 1.0)

    def test_proportions_sum_to_one(self, synthetic_data):
        """Test that proportions sum to 1."""
        Y, X, coords, _ = synthetic_data

        model = FlashDeconv(sketch_dim=64, max_iter=30)
        model.fit(Y, X, coords)

        row_sums = model.proportions_.sum(axis=1)
        np.testing.assert_allclose(row_sums, 1.0)

    def test_non_negative_proportions(self, synthetic_data):
        """Test that proportions are non-negative."""
        Y, X, coords, _ = synthetic_data

        model = FlashDeconv(sketch_dim=64, max_iter=30)
        model.fit(Y, X, coords)

        assert np.all(model.proportions_ >= 0)

    def test_dominant_cell_type(self, synthetic_data):
        """Test get_dominant_cell_type method."""
        Y, X, coords, _ = synthetic_data

        model = FlashDeconv(sketch_dim=64, max_iter=20)
        model.fit(Y, X, coords)

        dominant = model.get_dominant_cell_type()

        assert dominant.shape == (100,)
        assert np.all((dominant >= 0) & (dominant < 5))

    def test_auto_lambda(self, synthetic_data):
        """Test automatic lambda tuning."""
        Y, X, coords, _ = synthetic_data

        model = FlashDeconv(
            sketch_dim=64,
            lambda_spatial="auto",
            max_iter=20,
        )
        model.fit(Y, X, coords)

        assert hasattr(model, 'lambda_used_')
        assert model.lambda_used_ > 0

    def test_summary(self, synthetic_data):
        """Test summary method."""
        Y, X, coords, _ = synthetic_data

        model = FlashDeconv(sketch_dim=64, max_iter=20)
        model.fit(Y, X, coords)

        summary = model.summary()

        assert summary['fitted'] is True
        assert summary['n_spots'] == 100
        assert summary['n_cell_types'] == 5
        assert 'converged' in summary

    def test_sparse_input(self, synthetic_data):
        """Test with sparse input matrix."""
        Y, X, coords, _ = synthetic_data
        Y_sparse = sparse.csr_matrix(Y)

        model = FlashDeconv(sketch_dim=64, max_iter=20)
        model.fit(Y_sparse, X, coords)

        assert model._fitted
        assert model.proportions_.shape == (100, 5)

    def test_reproducibility(self, synthetic_data):
        """Test results are reproducible with same random state."""
        Y, X, coords, _ = synthetic_data

        model1 = FlashDeconv(sketch_dim=64, max_iter=30, random_state=123)
        prop1 = model1.fit_transform(Y, X, coords)

        model2 = FlashDeconv(sketch_dim=64, max_iter=30, random_state=123)
        prop2 = model2.fit_transform(Y, X, coords)

        np.testing.assert_allclose(prop1, prop2)

    def test_different_sketch_dims(self, synthetic_data):
        """Test with different sketch dimensions."""
        Y, X, coords, _ = synthetic_data

        for d in [32, 64, 128]:
            model = FlashDeconv(sketch_dim=d, max_iter=20)
            model.fit(Y, X, coords)

            assert model.proportions_.shape == (100, 5)

    def test_recovery_accuracy(self, synthetic_data):
        """Test that estimates correlate with true proportions."""
        Y, X, coords, beta_true = synthetic_data

        model = FlashDeconv(
            sketch_dim=128,
            lambda_spatial=0.1,
            max_iter=100,
            random_state=42,
        )
        proportions = model.fit_transform(Y, X, coords)

        # Compute correlation
        corr = np.corrcoef(proportions.flatten(), beta_true.flatten())[0, 1]

        # Should have positive correlation with ground truth
        assert corr > 0.3, f"Correlation too low: {corr}"


class TestFlashDeconvErrors:
    """Test error handling."""

    def test_unfitted_error(self):
        """Test error when accessing results before fitting."""
        model = FlashDeconv()

        with pytest.raises(RuntimeError, match="not been fitted"):
            model.get_cell_type_proportions()

    def test_dimension_mismatch(self):
        """Test error with mismatched dimensions."""
        Y = np.random.rand(50, 100)  # 100 genes
        X = np.random.rand(5, 200)   # 200 genes (mismatch!)
        coords = np.random.rand(50, 2)

        model = FlashDeconv(sketch_dim=32, max_iter=10)

        # This should raise an error due to gene dimension mismatch
        with pytest.raises(Exception):
            model.fit(Y, X, coords)
