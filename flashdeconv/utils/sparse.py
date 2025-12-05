"""
Sparse matrix utilities for FlashDeconv.

This module provides helper functions for efficient sparse matrix operations.
"""

import numpy as np
from scipy import sparse
from typing import Union, Tuple

ArrayLike = Union[np.ndarray, sparse.spmatrix]


def ensure_sparse_csr(A: ArrayLike) -> sparse.csr_matrix:
    """
    Ensure matrix is in CSR sparse format.

    Parameters
    ----------
    A : array-like or sparse matrix
        Input matrix.

    Returns
    -------
    A_csr : sparse.csr_matrix
        Matrix in CSR format.
    """
    if sparse.issparse(A):
        return A.tocsr()
    else:
        return sparse.csr_matrix(A)


def ensure_dense(A: ArrayLike) -> np.ndarray:
    """
    Ensure matrix is dense numpy array.

    Parameters
    ----------
    A : array-like or sparse matrix
        Input matrix.

    Returns
    -------
    A_dense : ndarray
        Dense numpy array.
    """
    if sparse.issparse(A):
        return A.toarray()
    else:
        return np.asarray(A)


def sparse_row_norms(A: sparse.spmatrix, ord: int = 2) -> np.ndarray:
    """
    Compute row norms of sparse matrix efficiently.

    Parameters
    ----------
    A : sparse matrix
        Input sparse matrix.
    ord : int, default=2
        Norm order (1 for L1, 2 for L2).

    Returns
    -------
    norms : ndarray of shape (n_rows,)
        Row norms.
    """
    A_csr = A.tocsr()

    if ord == 1:
        return np.asarray(np.abs(A_csr).sum(axis=1)).flatten()
    elif ord == 2:
        return np.sqrt(np.asarray(A_csr.power(2).sum(axis=1)).flatten())
    else:
        raise ValueError(f"Unsupported norm order: {ord}")


def sparse_col_norms(A: sparse.spmatrix, ord: int = 2) -> np.ndarray:
    """
    Compute column norms of sparse matrix efficiently.

    Parameters
    ----------
    A : sparse matrix
        Input sparse matrix.
    ord : int, default=2
        Norm order.

    Returns
    -------
    norms : ndarray of shape (n_cols,)
        Column norms.
    """
    A_csc = A.tocsc()

    if ord == 1:
        return np.asarray(np.abs(A_csc).sum(axis=0)).flatten()
    elif ord == 2:
        return np.sqrt(np.asarray(A_csc.power(2).sum(axis=0)).flatten())
    else:
        raise ValueError(f"Unsupported norm order: {ord}")


def sparse_normalize_rows(A: sparse.spmatrix, norm: str = "l2") -> sparse.csr_matrix:
    """
    Normalize rows of sparse matrix.

    Parameters
    ----------
    A : sparse matrix
        Input matrix.
    norm : str, default="l2"
        Normalization type: "l1", "l2", or "max".

    Returns
    -------
    A_norm : sparse.csr_matrix
        Row-normalized matrix.
    """
    A_csr = A.tocsr().astype(np.float64)

    if norm == "l1":
        norms = sparse_row_norms(A_csr, ord=1)
    elif norm == "l2":
        norms = sparse_row_norms(A_csr, ord=2)
    elif norm == "max":
        norms = np.asarray(np.abs(A_csr).max(axis=1).todense()).flatten()
    else:
        raise ValueError(f"Unsupported norm: {norm}")

    norms = np.maximum(norms, 1e-10)
    D_inv = sparse.diags(1.0 / norms, format='csr')

    return D_inv @ A_csr


def sparse_column_subset(
    A: sparse.spmatrix,
    indices: np.ndarray,
) -> sparse.csr_matrix:
    """
    Efficiently subset columns of sparse matrix.

    Parameters
    ----------
    A : sparse matrix
        Input matrix.
    indices : ndarray
        Column indices to keep.

    Returns
    -------
    A_subset : sparse.csr_matrix
        Subsetted matrix.
    """
    return A.tocsc()[:, indices].tocsr()


def sparse_row_subset(
    A: sparse.spmatrix,
    indices: np.ndarray,
) -> sparse.csr_matrix:
    """
    Efficiently subset rows of sparse matrix.

    Parameters
    ----------
    A : sparse matrix
        Input matrix.
    indices : ndarray
        Row indices to keep.

    Returns
    -------
    A_subset : sparse.csr_matrix
        Subsetted matrix.
    """
    return A.tocsr()[indices, :]


def create_sparse_identity(n: int) -> sparse.csr_matrix:
    """
    Create sparse identity matrix.

    Parameters
    ----------
    n : int
        Size of matrix.

    Returns
    -------
    I : sparse.csr_matrix
        Identity matrix.
    """
    return sparse.eye(n, format='csr', dtype=np.float64)


def sparse_block_diag(matrices: list) -> sparse.csr_matrix:
    """
    Create block diagonal sparse matrix.

    Parameters
    ----------
    matrices : list of sparse matrices
        Diagonal blocks.

    Returns
    -------
    block : sparse.csr_matrix
        Block diagonal matrix.
    """
    return sparse.block_diag(matrices, format='csr')


def compute_sparse_nnz_stats(A: sparse.spmatrix) -> dict:
    """
    Compute statistics about sparsity of matrix.

    Parameters
    ----------
    A : sparse matrix
        Input matrix.

    Returns
    -------
    stats : dict
        Sparsity statistics.
    """
    A_csr = A.tocsr()
    n_rows, n_cols = A_csr.shape
    nnz = A_csr.nnz
    total = n_rows * n_cols

    return {
        "shape": (n_rows, n_cols),
        "nnz": nnz,
        "density": nnz / total if total > 0 else 0.0,
        "sparsity": 1 - (nnz / total) if total > 0 else 1.0,
        "nnz_per_row_mean": nnz / n_rows if n_rows > 0 else 0.0,
        "nnz_per_row_std": np.std(np.diff(A_csr.indptr)),
        "memory_dense_mb": total * 8 / 1e6,
        "memory_sparse_mb": (A_csr.data.nbytes + A_csr.indices.nbytes +
                            A_csr.indptr.nbytes) / 1e6,
    }
