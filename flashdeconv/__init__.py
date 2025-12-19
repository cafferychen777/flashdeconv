"""
FlashDeconv: Fast Linear Algebra for Scalable Hybrid Deconvolution

A high-performance spatial transcriptomics deconvolution method that combines:
- Variance-stabilizing transformation with platform effect correction
- Structure-preserving randomized sketching
- Spatial graph Laplacian regularization
- Numba-accelerated Block Coordinate Descent solver

Example
-------
Scanpy-style API (recommended):

>>> import flashdeconv as fd
>>> fd.tl.deconvolve(adata_st, adata_ref, cell_type_key="celltype")
>>> adata_st.obsm['flashdeconv']  # cell type proportions

NumPy API (for more control):

>>> from flashdeconv import FlashDeconv
>>> model = FlashDeconv(sketch_dim=512, lambda_spatial=5000)
>>> proportions = model.fit_transform(Y, X, coords)
"""

__version__ = "0.1.1"
__author__ = "FlashDeconv Team"

from flashdeconv.core.deconv import FlashDeconv
from flashdeconv import tl

__all__ = ["FlashDeconv", "tl", "__version__"]
