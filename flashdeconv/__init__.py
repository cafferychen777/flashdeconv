"""
FlashDeconv: Fast Linear Algebra for Scalable Hybrid Deconvolution

A high-performance spatial transcriptomics deconvolution method that combines:
- Log-normalization of spatial and reference expression
- Structure-preserving gene representation (deterministic expected
  leverage-weighted CountSketch weights; randomized CountSketch as legacy option)
- Spatial graph Laplacian regularization
- Numba-accelerated Block Coordinate Descent solver
- Diagnostics for cell types missing from the reference

Example
-------
Scanpy-style API (recommended):

>>> import flashdeconv as fd
>>> fd.tl.deconvolve(adata_st, adata_ref, cell_type_key="celltype")
>>> adata_st.obsm['flashdeconv']  # cell type proportions

NumPy API (for more control):

>>> from flashdeconv import FlashDeconv
>>> model = FlashDeconv(sketch_dim=512)
>>> proportions = model.fit_transform(Y, X, coords)
"""

__version__ = "0.2.0"
__author__ = "FlashDeconv Team"

from flashdeconv.core.deconv import FlashDeconv
from flashdeconv.core.refcheck import (
    reference_fit_scores,
    suggest_missing_types,
    unexplained_genes,
)
from flashdeconv import tl

__all__ = [
    "FlashDeconv",
    "tl",
    "reference_fit_scores",
    "unexplained_genes",
    "suggest_missing_types",
    "__version__",
]
