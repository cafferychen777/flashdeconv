"""
Expression preprocessing for FlashDeconv.

Transforms the spatial counts ``Y`` and the reference signatures ``X`` (both
restricted to the selected genes) into the space in which the regression is
solved. Shared by :class:`flashdeconv.FlashDeconv` and the Poisson bootstrap.
"""

from typing import Tuple, Union

import numpy as np
from scipy import sparse

ArrayLike = Union[np.ndarray, sparse.spmatrix]

#: Accepted values of the ``preprocess`` parameter.
PREPROCESS_METHODS = ("log_cpm", "pearson", "raw")

# Scale factor of the library-size normalization (counts per 10,000).
LIBRARY_SCALE = 1e4

# Fixed negative-binomial dispersion of the uncentered Pearson residuals.
PEARSON_THETA = 100.0


def preprocess_expression(
    Y: ArrayLike,
    X: np.ndarray,
    method: str = "log_cpm",
) -> Tuple[ArrayLike, np.ndarray]:
    """
    Preprocess spatial counts and reference signatures.

    Sparse ``Y`` stays sparse (``log1p(0) = 0``).

    Parameters
    ----------
    Y : array-like of shape (n_spots, n_genes)
        Spatial count matrix (sparse or dense).
    X : ndarray of shape (n_cell_types, n_genes)
        Reference signature matrix.
    method : {"log_cpm", "pearson", "raw"}, default="log_cpm"
        - "log_cpm": ``log1p`` of library-size-normalized expression scaled
          to 10,000 counts per spot / cell type (CP10k; the name is kept for
          backward compatibility). Applied to each row of ``Y`` and ``X``.
        - "pearson": uncentered Pearson residuals ``Y / sigma`` with
          ``sigma^2 = mu + mu^2 / 100`` and ``mu`` the per-gene mean
          (computed separately for ``Y`` and ``X``); values stay
          non-negative.
        - "raw": no transformation (cast to float64).

    Returns
    -------
    Y_norm : array-like of shape (n_spots, n_genes)
        Preprocessed spatial data (sparse if ``Y`` was sparse).
    X_norm : ndarray of shape (n_cell_types, n_genes)
        Preprocessed reference.
    """
    if method == "log_cpm":
        if sparse.issparse(Y):
            lib_size = np.array(Y.sum(axis=1)).flatten()
            lib_size[lib_size == 0] = 1.0
            Y_norm = sparse.diags(LIBRARY_SCALE / lib_size) @ Y
            # log1p on the stored values only: zeros stay zeros.
            Y_norm.data = np.log1p(Y_norm.data)
        else:
            Y_cpm = Y / (Y.sum(axis=1, keepdims=True) + 1e-10) * LIBRARY_SCALE
            Y_norm = np.log1p(Y_cpm)

        X_cpm = X / (X.sum(axis=1, keepdims=True) + 1e-10) * LIBRARY_SCALE
        X_norm = np.log1p(X_cpm)
        return Y_norm, X_norm

    if method == "pearson":
        # Divide by sigma only (no centering) so that NNLS inputs stay >= 0.
        if sparse.issparse(Y):
            Y_mean = np.asarray(Y.mean(axis=0)).flatten() + 1e-6
            Y_sigma = np.sqrt(Y_mean + Y_mean**2 / PEARSON_THETA)
            Y_norm = Y.multiply(1.0 / Y_sigma)
        else:
            Y_mean = Y.mean(axis=0, keepdims=True) + 1e-6
            Y_sigma = np.sqrt(Y_mean + Y_mean**2 / PEARSON_THETA)
            Y_norm = Y / Y_sigma

        X_mean = X.mean(axis=0, keepdims=True) + 1e-6
        X_sigma = np.sqrt(X_mean + X_mean**2 / PEARSON_THETA)
        X_norm = X / X_sigma
        return Y_norm, X_norm

    if method == "raw":
        return Y.astype(np.float64, copy=False), X.astype(np.float64, copy=False)

    raise ValueError(
        f"Unknown preprocess method: {method!r}. "
        f"Choose from {', '.join(repr(m) for m in PREPROCESS_METHODS)}."
    )
