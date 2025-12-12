"""Tests for spatial graph and Laplacian modules."""

import numpy as np
import pytest
from scipy import sparse

from flashdeconv.utils.graph import (
    build_knn_graph,
    build_radius_graph,
    coords_to_adjacency,
)
from flashdeconv.core.spatial import (
    compute_laplacian,
    get_neighbor_indices,
    get_neighbor_counts,
    compute_laplacian_quadratic,
)


class TestKNNGraph:
    """Tests for KNN graph construction."""

    def test_output_shape(self):
        """Test adjacency matrix shape."""
        coords = np.random.rand(50, 2)
        A = build_knn_graph(coords, k=6)

        assert A.shape == (50, 50)

    def test_symmetry(self):
        """Test graph is symmetric (undirected)."""
        coords = np.random.rand(30, 2)
        A = build_knn_graph(coords, k=4)

        diff = A - A.T
        assert diff.nnz == 0 or np.allclose(diff.data, 0)

    def test_no_self_loops_default(self):
        """Test no self-loops by default."""
        coords = np.random.rand(20, 2)
        A = build_knn_graph(coords, k=3, include_self=False)

        diag = A.diagonal()
        assert np.all(diag == 0)

    def test_with_self_loops(self):
        """Test self-loops when requested."""
        coords = np.random.rand(20, 2)
        A = build_knn_graph(coords, k=3, include_self=True)

        diag = A.diagonal()
        assert np.all(diag > 0)


class TestRadiusGraph:
    """Tests for radius-based graph."""

    def test_grid_pattern(self):
        """Test on regular grid."""
        # Create 3x3 grid
        x = np.arange(3)
        y = np.arange(3)
        xx, yy = np.meshgrid(x, y)
        coords = np.column_stack([xx.ravel(), yy.ravel()])

        # With radius=1.5, center connects to all 8 neighbors (including diagonals at sqrt(2)≈1.414)
        A = build_radius_graph(coords, radius=1.5)
        center_idx = 4  # Index of (1,1)
        n_neighbors = int(A[center_idx].sum())
        assert n_neighbors == 8  # All 8 surrounding cells

        # With radius=1.1, center connects only to 4 direct neighbors (not diagonals)
        A_small = build_radius_graph(coords, radius=1.1)
        n_neighbors_small = int(A_small[center_idx].sum())
        assert n_neighbors_small == 4


class TestCoordsToAdjacency:
    """Tests for coords_to_adjacency function."""

    def test_knn_method(self):
        """Test KNN method."""
        coords = np.random.rand(30, 2)
        A = coords_to_adjacency(coords, method="knn", k=5)

        assert A.shape == (30, 30)

    def test_radius_method(self):
        """Test radius method."""
        coords = np.random.rand(30, 2)
        A = coords_to_adjacency(coords, method="radius", radius=0.3)

        assert A.shape == (30, 30)

    def test_grid_method(self):
        """Test grid method."""
        # Create regular grid
        x = np.arange(5)
        y = np.arange(5)
        xx, yy = np.meshgrid(x, y)
        coords = np.column_stack([xx.ravel(), yy.ravel()])

        A = coords_to_adjacency(coords, method="grid")

        assert A.shape == (25, 25)


class TestLaplacian:
    """Tests for graph Laplacian."""

    def test_laplacian_shape(self):
        """Test Laplacian has correct shape."""
        A = sparse.random(50, 50, density=0.1, format='csr')
        A = A + A.T  # Make symmetric

        L = compute_laplacian(A)

        assert L.shape == (50, 50)

    def test_laplacian_row_sum(self):
        """Test unnormalized Laplacian rows sum to zero."""
        A = sparse.random(30, 30, density=0.2, format='csr')
        A = A + A.T
        A.setdiag(0)  # Remove self-loops

        L = compute_laplacian(A, normalized=False)

        row_sums = np.asarray(L.sum(axis=1)).flatten()
        np.testing.assert_allclose(row_sums, 0, atol=1e-10)

    def test_normalized_laplacian(self):
        """Test normalized Laplacian properties."""
        A = sparse.random(30, 30, density=0.2, format='csr')
        A = A + A.T
        A.setdiag(0)

        L = compute_laplacian(A, normalized=True)

        # Diagonal should be <= 1
        diag = L.diagonal()
        assert np.all(diag <= 1 + 1e-10)


class TestNeighborHelpers:
    """Tests for neighbor helper functions."""

    def test_neighbor_indices(self):
        """Test neighbor index extraction."""
        coords = np.random.rand(20, 2)
        A = build_knn_graph(coords, k=4)

        neighbors = get_neighbor_indices(A)

        assert len(neighbors) == 20
        for i, nb in enumerate(neighbors):
            assert isinstance(nb, np.ndarray)

    def test_neighbor_counts(self):
        """Test neighbor count computation."""
        coords = np.random.rand(20, 2)
        A = build_knn_graph(coords, k=4)

        counts = get_neighbor_counts(A)

        assert counts.shape == (20,)
        assert np.all(counts >= 0)


class TestLaplacianQuadratic:
    """Tests for Laplacian quadratic form."""

    def test_positive_semidefinite(self):
        """Test Tr(beta^T L beta) >= 0."""
        coords = np.random.rand(30, 2)
        A = build_knn_graph(coords, k=5)
        L = compute_laplacian(A)

        beta = np.random.rand(30, 5)

        value = compute_laplacian_quadratic(beta, L)

        assert value >= -1e-10  # Allow small numerical errors

    def test_zero_for_constant(self):
        """Test Tr(beta^T L beta) = 0 for constant beta."""
        coords = np.random.rand(20, 2)
        A = build_knn_graph(coords, k=4)
        L = compute_laplacian(A)

        # Constant beta
        beta = np.ones((20, 3)) * 0.5

        value = compute_laplacian_quadratic(beta, L)

        np.testing.assert_allclose(value, 0, atol=1e-10)
