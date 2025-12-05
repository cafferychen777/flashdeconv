"""
Data loading utilities for FlashDeconv.

This module provides functions to load spatial transcriptomics data
and single-cell references from AnnData objects.
"""

import numpy as np
from scipy import sparse
from typing import Union, Optional, Tuple, Any

ArrayLike = Union[np.ndarray, sparse.spmatrix]


def load_spatial_data(
    adata: Any,
    layer: Optional[str] = None,
    coord_key: str = "spatial",
) -> Tuple[ArrayLike, np.ndarray, np.ndarray]:
    """
    Extract spatial data from AnnData object.

    Parameters
    ----------
    adata : AnnData
        Spatial transcriptomics AnnData object.
    layer : str, optional
        Layer to use for counts. Uses .X if not specified.
    coord_key : str, default="spatial"
        Key in adata.obsm for spatial coordinates.

    Returns
    -------
    Y : array-like of shape (n_spots, n_genes)
        Count matrix.
    coords : ndarray of shape (n_spots, 2)
        Spatial coordinates.
    gene_names : ndarray of shape (n_genes,)
        Gene names.
    """
    # Get expression matrix
    if layer is not None:
        Y = adata.layers[layer]
    else:
        Y = adata.X

    # Get coordinates
    if coord_key in adata.obsm:
        coords = np.array(adata.obsm[coord_key])
    elif "X_spatial" in adata.obsm:
        coords = np.array(adata.obsm["X_spatial"])
    else:
        # Try to get from obs columns
        if "x" in adata.obs and "y" in adata.obs:
            coords = np.column_stack([adata.obs["x"], adata.obs["y"]])
        elif "array_row" in adata.obs and "array_col" in adata.obs:
            coords = np.column_stack([
                adata.obs["array_row"],
                adata.obs["array_col"]
            ])
        else:
            raise ValueError(
                f"Could not find spatial coordinates. "
                f"Expected key '{coord_key}' in adata.obsm or 'x'/'y' in adata.obs"
            )

    # Get gene names
    gene_names = np.array(adata.var_names)

    return Y, coords, gene_names


def load_reference(
    adata_ref: Any,
    cell_type_key: str = "cell_type",
    layer: Optional[str] = None,
    method: str = "mean",
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Extract reference signature matrix from single-cell AnnData.

    Parameters
    ----------
    adata_ref : AnnData
        Single-cell reference AnnData object.
    cell_type_key : str, default="cell_type"
        Key in adata_ref.obs for cell type annotations.
    layer : str, optional
        Layer to use. Uses .X if not specified.
    method : str, default="mean"
        Aggregation method: "mean" or "sum".

    Returns
    -------
    X : ndarray of shape (n_cell_types, n_genes)
        Reference signature matrix.
    cell_type_names : ndarray of shape (n_cell_types,)
        Cell type names.
    gene_names : ndarray of shape (n_genes,)
        Gene names.
    """
    # Get expression matrix
    if layer is not None:
        expr = adata_ref.layers[layer]
    else:
        expr = adata_ref.X

    if sparse.issparse(expr):
        expr = expr.toarray()
    else:
        expr = np.array(expr)

    # Get cell type labels
    if cell_type_key not in adata_ref.obs:
        raise ValueError(f"Cell type key '{cell_type_key}' not found in adata.obs")

    cell_types = np.array(adata_ref.obs[cell_type_key])
    unique_types = np.unique(cell_types)
    n_cell_types = len(unique_types)
    n_genes = expr.shape[1]

    # Aggregate by cell type
    X = np.zeros((n_cell_types, n_genes), dtype=np.float64)

    for i, ct in enumerate(unique_types):
        mask = cell_types == ct
        if method == "mean":
            X[i] = np.mean(expr[mask], axis=0)
        elif method == "sum":
            X[i] = np.sum(expr[mask], axis=0)
        else:
            raise ValueError(f"Unknown aggregation method: {method}")

    gene_names = np.array(adata_ref.var_names)

    return X, unique_types, gene_names


def align_genes(
    Y: ArrayLike,
    X: np.ndarray,
    genes_spatial: np.ndarray,
    genes_ref: np.ndarray,
) -> Tuple[ArrayLike, np.ndarray, np.ndarray]:
    """
    Align genes between spatial data and reference.

    Parameters
    ----------
    Y : array-like of shape (n_spots, n_genes_spatial)
        Spatial count matrix.
    X : ndarray of shape (n_cell_types, n_genes_ref)
        Reference signatures.
    genes_spatial : ndarray
        Gene names in spatial data.
    genes_ref : ndarray
        Gene names in reference.

    Returns
    -------
    Y_aligned : array-like of shape (n_spots, n_common)
        Aligned spatial data.
    X_aligned : ndarray of shape (n_cell_types, n_common)
        Aligned reference.
    common_genes : ndarray
        Names of common genes.
    """
    # Find common genes
    common_genes = np.intersect1d(genes_spatial, genes_ref)

    if len(common_genes) == 0:
        raise ValueError("No common genes found between spatial data and reference")

    # Get indices
    spatial_idx = np.array([
        np.where(genes_spatial == g)[0][0] for g in common_genes
    ])
    ref_idx = np.array([
        np.where(genes_ref == g)[0][0] for g in common_genes
    ])

    # Subset
    if sparse.issparse(Y):
        Y_aligned = Y[:, spatial_idx]
    else:
        Y_aligned = Y[:, spatial_idx]

    X_aligned = X[:, ref_idx]

    return Y_aligned, X_aligned, common_genes


def result_to_anndata(
    beta: np.ndarray,
    adata: Any,
    cell_type_names: Optional[np.ndarray] = None,
    key_added: str = "flashdeconv",
) -> Any:
    """
    Store deconvolution results in AnnData object.

    Parameters
    ----------
    beta : ndarray of shape (n_spots, n_cell_types)
        Cell type proportions.
    adata : AnnData
        Spatial AnnData object.
    cell_type_names : ndarray, optional
        Names of cell types for columns.
    key_added : str, default="flashdeconv"
        Key to store results under in adata.obsm.

    Returns
    -------
    adata : AnnData
        Modified AnnData with results in .obsm[key_added].
    """
    import pandas as pd

    # Create DataFrame with results
    if cell_type_names is not None:
        columns = cell_type_names
    else:
        columns = [f"CellType_{i}" for i in range(beta.shape[1])]

    df = pd.DataFrame(
        beta,
        index=adata.obs_names,
        columns=columns,
    )

    # Store in obsm
    adata.obsm[key_added] = df

    # Also store proportions in obs for easy access
    for col in columns:
        adata.obs[f"{key_added}_{col}"] = df[col].values

    # Store dominant cell type
    adata.obs[f"{key_added}_dominant"] = columns[np.argmax(beta, axis=1)]

    return adata


def prepare_data(
    adata_st: Any,
    adata_ref: Any,
    cell_type_key: str = "cell_type",
    spatial_coord_key: str = "spatial",
    layer_st: Optional[str] = None,
    layer_ref: Optional[str] = None,
) -> Tuple[ArrayLike, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Convenience function to prepare data for FlashDeconv.

    Parameters
    ----------
    adata_st : AnnData
        Spatial transcriptomics data.
    adata_ref : AnnData
        Single-cell reference data.
    cell_type_key : str, default="cell_type"
        Key for cell type annotations in reference.
    spatial_coord_key : str, default="spatial"
        Key for spatial coordinates.
    layer_st : str, optional
        Layer for spatial data.
    layer_ref : str, optional
        Layer for reference data.

    Returns
    -------
    Y : array-like
        Aligned spatial counts.
    X : ndarray
        Aligned reference signatures.
    coords : ndarray
        Spatial coordinates.
    cell_type_names : ndarray
        Cell type names.
    gene_names : ndarray
        Common gene names.
    """
    # Load data
    Y, coords, genes_st = load_spatial_data(
        adata_st, layer=layer_st, coord_key=spatial_coord_key
    )
    X, cell_type_names, genes_ref = load_reference(
        adata_ref, cell_type_key=cell_type_key, layer=layer_ref
    )

    # Align genes
    Y, X, gene_names = align_genes(Y, X, genes_st, genes_ref)

    return Y, X, coords, cell_type_names, gene_names
