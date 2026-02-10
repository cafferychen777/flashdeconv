"""
Spatial graph construction utilities for FlashDeconv.

This module implements various methods to construct spatial neighbor
graphs from spot coordinates.
"""

import numpy as np
from scipy import sparse
from scipy.spatial import cKDTree
from typing import Union, Optional, Tuple, Literal

ArrayLike = Union[np.ndarray, sparse.spmatrix]


def build_knn_graph(
    coords: np.ndarray,
    k: int = 6,
    include_self: bool = False,
) -> sparse.csr_matrix:
    """
    Build k-nearest neighbor spatial graph.

    Parameters
    ----------
    coords : ndarray of shape (n_spots, 2) or (n_spots, 3)
        Spatial coordinates of spots.
    k : int, default=6
        Number of nearest neighbors.
    include_self : bool, default=False
        Whether to include self-loops.

    Returns
    -------
    A : sparse.csr_matrix of shape (n_spots, n_spots)
        Binary adjacency matrix.
    """
    n_spots = coords.shape[0]

    # Clamp k to valid range: cannot query more neighbors than exist
    k_actual = min(k, n_spots - 1)

    if k_actual <= 0:
        # Single spot or empty: no neighbors possible
        if include_self and n_spots > 0:
            return sparse.eye(n_spots, dtype=np.float64, format="csr")
        return sparse.csr_matrix((n_spots, n_spots), dtype=np.float64)

    # Build KD-tree for efficient neighbor search
    tree = cKDTree(coords)

    # Query k_actual+1 neighbors (includes self)
    distances, indices = tree.query(coords, k=k_actual + 1)

    # Build adjacency matrix (vectorized)
    # indices has shape (n_spots, k_actual+1), includes self
    row_idx = np.repeat(np.arange(n_spots), k_actual + 1)
    col_idx = indices.ravel()

    if not include_self:
        # Remove self-loops
        mask = row_idx != col_idx
        row_idx = row_idx[mask]
        col_idx = col_idx[mask]

    data = np.ones(len(row_idx), dtype=np.float64)
    A = sparse.csr_matrix((data, (row_idx, col_idx)), shape=(n_spots, n_spots))

    # Make symmetric (undirected graph)
    A = A + A.T
    A.data[:] = 1.0  # Binary adjacency

    return A


def build_radius_graph(
    coords: np.ndarray,
    radius: float,
    include_self: bool = False,
) -> sparse.csr_matrix:
    """
    Build radius-based neighbor graph.

    Parameters
    ----------
    coords : ndarray of shape (n_spots, 2) or (n_spots, 3)
        Spatial coordinates of spots.
    radius : float
        Maximum distance for two spots to be neighbors.
    include_self : bool, default=False
        Whether to include self-loops.

    Returns
    -------
    A : sparse.csr_matrix of shape (n_spots, n_spots)
        Binary adjacency matrix.
    """
    n_spots = coords.shape[0]

    # Build KD-tree
    tree = cKDTree(coords)

    # Query all pairs within radius
    pairs = tree.query_pairs(r=radius, output_type='ndarray')

    if len(pairs) == 0:
        # No neighbors found - return empty matrix
        return sparse.csr_matrix((n_spots, n_spots), dtype=np.float64)

    # Build symmetric adjacency
    rows = np.concatenate([pairs[:, 0], pairs[:, 1]])
    cols = np.concatenate([pairs[:, 1], pairs[:, 0]])
    data = np.ones(len(rows), dtype=np.float64)

    A = sparse.csr_matrix((data, (rows, cols)), shape=(n_spots, n_spots))

    if include_self:
        A = A + sparse.eye(n_spots, dtype=np.float64)

    return A


def build_grid_graph(
    coords: np.ndarray,
    grid_spacing: Optional[float] = None,
) -> sparse.csr_matrix:
    """
    Build graph assuming regular grid structure (e.g., Visium).

    Connects spots to their 6 hexagonal neighbors or 4/8 grid neighbors.

    Parameters
    ----------
    coords : ndarray of shape (n_spots, 2)
        Spatial coordinates of spots.
    grid_spacing : float, optional
        Expected spacing between grid points. Auto-detected if not provided.

    Returns
    -------
    A : sparse.csr_matrix of shape (n_spots, n_spots)
        Binary adjacency matrix.
    """
    n_spots = coords.shape[0]

    if grid_spacing is None:
        # Auto-detect grid spacing from nearest neighbor distances
        tree = cKDTree(coords)
        distances, _ = tree.query(coords, k=2)
        grid_spacing = np.median(distances[:, 1])

    # Use slightly larger radius to account for hexagonal grids
    radius = grid_spacing * 1.5

    return build_radius_graph(coords, radius)


def coords_to_adjacency(
    coords: np.ndarray,
    method: Literal["knn", "radius", "grid"] = "knn",
    k: int = 6,
    radius: Optional[float] = None,
) -> sparse.csr_matrix:
    """
    Convert spatial coordinates to adjacency matrix.

    Parameters
    ----------
    coords : ndarray of shape (n_spots, 2) or (n_spots, 3)
        Spatial coordinates.
    method : str, default="knn"
        Graph construction method:
        - "knn": k-nearest neighbors
        - "radius": fixed radius
        - "grid": regular grid (auto-detect spacing)
    k : int, default=6
        Number of neighbors for KNN method.
    radius : float, optional
        Radius for radius-based method.

    Returns
    -------
    A : sparse.csr_matrix
        Adjacency matrix.
    """
    if method == "knn":
        return build_knn_graph(coords, k=k)
    elif method == "radius":
        if radius is None:
            raise ValueError("radius must be specified for radius method")
        return build_radius_graph(coords, radius=radius)
    elif method == "grid":
        return build_grid_graph(coords)
    else:
        raise ValueError(f"Unknown method: {method}")


def estimate_optimal_k(
    coords: np.ndarray,
    min_k: int = 4,
    max_k: int = 20,
) -> int:
    """
    Estimate optimal k for KNN graph based on local density.

    Parameters
    ----------
    coords : ndarray of shape (n_spots, 2)
        Spatial coordinates.
    min_k : int, default=4
        Minimum k value.
    max_k : int, default=20
        Maximum k value.

    Returns
    -------
    k : int
        Estimated optimal k.
    """
    n_spots = coords.shape[0]

    # Parameter validation
    if min_k < 0:
        raise ValueError(f"min_k must be non-negative, got {min_k}")
    if max_k < min_k:
        raise ValueError(f"max_k must be >= min_k, got max_k={max_k} < min_k={min_k}")

    # Cannot have neighbors with fewer than 2 points
    if n_spots <= 1:
        return 0

    # Clamp to valid range: k can never exceed n_spots-1
    upper_bound = n_spots - 1
    max_k_eff = min(max_k, upper_bound)
    min_k_eff = min(min_k, upper_bound)

    # For small datasets, use smaller k
    if n_spots < 100:
        return min_k_eff

    # Estimate based on local density variation
    tree = cKDTree(coords)
    distances, _ = tree.query(coords, k=max_k_eff + 1)

    # Find k where distance growth stabilizes
    mean_distances = np.mean(distances, axis=0)[1:]  # Exclude self

    # Look for "elbow" in distance growth
    growth_rates = np.diff(mean_distances)
    normalized_growth = growth_rates / (mean_distances[:-1] + 1e-10)

    # Find where growth rate drops significantly
    threshold = np.median(normalized_growth) * 0.5
    stable_k = np.argmax(normalized_growth < threshold) + 1

    k = max(min_k_eff, min(stable_k + min_k_eff, max_k_eff))

    return k
