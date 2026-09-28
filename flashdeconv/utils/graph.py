"""
Spatial graph construction utilities for FlashDeconv.

This module implements various methods to construct spatial neighbor
graphs from spot coordinates.
"""

import numpy as np
from scipy import sparse
from scipy.spatial import cKDTree
from typing import Union, Optional, Literal

ArrayLike = Union[np.ndarray, sparse.spmatrix]


def _validate_coords(coords: np.ndarray) -> None:
    """Check that coords is a 2D array with at least 1 coordinate dimension."""
    if coords.ndim != 2 or coords.shape[1] == 0:
        raise ValueError(
            f"coords must be 2D with at least 1 coordinate dimension, "
            f"got shape {coords.shape}"
        )


def _knn_indices_tiebreak(coords: np.ndarray, k: int) -> np.ndarray:
    """
    Indices of the ``k`` nearest other spots, ties broken by spot index.

    The order in which ``cKDTree.query`` returns equidistant points depends
    on the tree traversal and differs across platforms (e.g. macOS arm64 vs
    Linux x86-64). On regular grids (Visium array coordinates, simulated
    grids) the k-th nearest distance is typically shared by several spots,
    so the graph itself would differ. Here neighbors are ranked by
    (distance, spot index) over *all* candidates tied with the k-th
    distance, which makes the graph a function of the coordinates only.

    Returns an array of shape (n_spots, k); requires 1 <= k <= n_spots - 1.
    """
    n = coords.shape[0]
    tree = cKDTree(coords)
    spot_idx = np.arange(n)
    out = np.empty((n, k), dtype=np.intp)
    todo = spot_idx
    m = min(k + 2, n)
    while todo.size:
        dist, ind = tree.query(coords[todo], k=m)
        dist = dist.reshape(len(todo), m)
        ind = ind.reshape(len(todo), m)
        # Drop the spot itself (it may be absent when it has > m-1 duplicates);
        # sort remaining candidates by (distance, index).
        not_self = ind != todo[:, None]
        dist = np.where(not_self, dist, np.inf)
        order = np.lexsort((ind, dist), axis=-1)
        dist = np.take_along_axis(dist, order, axis=1)
        ind = np.take_along_axis(ind, order, axis=1)
        # All candidates tied with the k-th neighbor were retrieved iff the
        # farthest retrieved other spot is strictly farther (or all spots
        # were retrieved).
        n_other = m - 1  # at least m - 1 non-self candidates per row
        complete = (m == n) | (dist[:, n_other - 1] > dist[:, k - 1])
        out[todo[complete]] = ind[complete, :k]
        todo = todo[~complete]
        m = min(2 * m, n)
    return out


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
    _validate_coords(coords)
    n_spots = coords.shape[0]

    # Clamp k to valid range: cannot query more neighbors than exist
    k_actual = min(k, n_spots - 1)

    if k_actual <= 0:
        # Single spot or empty: no neighbors possible
        if include_self and n_spots > 0:
            return sparse.eye(n_spots, dtype=np.float64, format="csr")
        return sparse.csr_matrix((n_spots, n_spots), dtype=np.float64)

    spot_idx = np.arange(n_spots)
    col_idx = _knn_indices_tiebreak(coords, k_actual).ravel()
    row_idx = np.repeat(spot_idx, k_actual)
    if include_self:
        row_idx = np.concatenate([row_idx, spot_idx])
        col_idx = np.concatenate([col_idx, spot_idx])

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
    _validate_coords(coords)
    n_spots = coords.shape[0]

    # Build KD-tree
    tree = cKDTree(coords)

    # Query all pairs within radius
    pairs = tree.query_pairs(r=radius, output_type='ndarray')

    if len(pairs) == 0:
        # No neighbors found
        if include_self and n_spots > 0:
            return sparse.eye(n_spots, dtype=np.float64, format="csr")
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
    _validate_coords(coords)
    n_spots = coords.shape[0]

    if n_spots <= 1:
        return sparse.csr_matrix((n_spots, n_spots), dtype=np.float64)

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
