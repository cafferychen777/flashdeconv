"""Core algorithms for FlashDeconv."""

from flashdeconv.core.deconv import FlashDeconv
from flashdeconv.core.sketching import (
    build_countsketch_matrix,
    project_to_sketch,
)
from flashdeconv.core.spatial import (
    compute_laplacian,
    get_neighbor_indices,
)
from flashdeconv.core.solver import bcd_solve

__all__ = [
    "FlashDeconv",
    "build_countsketch_matrix",
    "project_to_sketch",
    "compute_laplacian",
    "get_neighbor_indices",
    "bcd_solve",
]
