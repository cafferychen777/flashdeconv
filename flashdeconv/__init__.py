"""
FlashDeconv: Fast Linear Algebra for Scalable Hybrid Deconvolution

A high-performance spatial transcriptomics deconvolution method that combines:
- Variance-stabilizing transformation with platform effect correction
- Structure-preserving randomized sketching
- Spatial graph Laplacian regularization
- Numba-accelerated Block Coordinate Descent solver

Example
-------
>>> from flashdeconv import FlashDeconv
>>> model = FlashDeconv(sketch_dim=512, lambda_spatial=5000)
>>> beta = model.fit_transform(adata_st, adata_ref)
"""

__version__ = "0.1.0"
__author__ = "FlashDeconv Team"

from flashdeconv.core.deconv import FlashDeconv

__all__ = ["FlashDeconv", "__version__"]
