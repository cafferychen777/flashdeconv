"""Spatial deconvolution tools for FlashDeconv."""

from typing import Union, Optional, Any

import numpy as np


def deconvolve(
    adata_st: Any,
    adata_ref: Any,
    cell_type_key: str = "cell_type",
    *,
    sketch_dim: int = 512,
    lambda_spatial: Union[float, str] = 5000.0,
    n_hvg: int = 2000,
    n_markers_per_type: int = 50,
    layer_st: Optional[str] = None,
    layer_ref: Optional[str] = None,
    spatial_key: str = "spatial",
    key_added: str = "flashdeconv",
    copy: bool = False,
) -> Optional[Any]:
    """
    Run FlashDeconv spatial deconvolution.

    Estimates cell type proportions for each spatial location by decomposing
    the observed gene expression as a mixture of reference cell type signatures.

    Parameters
    ----------
    adata_st
        Spatial transcriptomics AnnData with coordinates in ``.obsm[spatial_key]``.
    adata_ref
        Single-cell reference AnnData with cell type labels in ``.obs[cell_type_key]``.
    cell_type_key
        Column in ``adata_ref.obs`` containing cell type annotations.
    sketch_dim
        Dimension of the sketched space. Higher values preserve more information
        but increase computation. Default: 512.
    lambda_spatial
        Spatial regularization strength. Use ``"auto"`` for automatic tuning.
        Higher values encourage spatially smooth cell type distributions.
    n_hvg
        Number of highly variable genes to select. Default: 2000.
    n_markers_per_type
        Number of marker genes per cell type. Default: 50.
    layer_st
        Layer in ``adata_st`` to use for counts. Uses ``.X`` if None.
    layer_ref
        Layer in ``adata_ref`` to use for counts. Uses ``.X`` if None.
    spatial_key
        Key in ``adata_st.obsm`` for spatial coordinates. Default: "spatial".
    key_added
        Key under which results are stored. Default: "flashdeconv".
    copy
        If True, return a copy instead of modifying ``adata_st`` in-place.

    Returns
    -------
    If ``copy=True``, returns modified AnnData. Otherwise, modifies ``adata_st``
    in-place and returns ``None``.

    Adds fields to ``adata_st``:

    - ``.obsm['{key_added}']`` : DataFrame of cell type proportions (n_spots x n_types)
    - ``.obs['{key_added}_{celltype}']`` : Proportion column for each cell type
    - ``.obs['{key_added}_dominant']`` : Dominant cell type per spot
    - ``.uns['{key_added}_params']`` : Parameters used for deconvolution

    Examples
    --------
    Basic usage:

    >>> import scanpy as sc
    >>> import flashdeconv as fd
    >>> adata_st = sc.read_h5ad("visium.h5ad")
    >>> adata_ref = sc.read_h5ad("reference.h5ad")
    >>> fd.tl.deconvolve(adata_st, adata_ref, cell_type_key="celltype")
    >>> adata_st.obsm['flashdeconv']  # cell type proportions

    Visualize with scanpy:

    >>> sc.pl.spatial(adata_st, color='flashdeconv_dominant')
    >>> sc.pl.spatial(adata_st, color=['flashdeconv_Neuron', 'flashdeconv_Astrocyte'])

    Use deconvolution for downstream clustering:

    >>> sc.pp.neighbors(adata_st, use_rep='flashdeconv')
    >>> sc.tl.leiden(adata_st)
    """
    from flashdeconv.core.deconv import FlashDeconv
    from flashdeconv.io import prepare_data, result_to_anndata

    adata = adata_st.copy() if copy else adata_st

    # Prepare data using existing IO functions
    Y, X, coords, cell_type_names, gene_names = prepare_data(
        adata,
        adata_ref,
        cell_type_key=cell_type_key,
        layer_st=layer_st,
        layer_ref=layer_ref,
        spatial_coord_key=spatial_key,
    )

    # Run deconvolution using existing core
    model = FlashDeconv(
        sketch_dim=sketch_dim,
        lambda_spatial=lambda_spatial,
        n_hvg=n_hvg,
        n_markers_per_type=n_markers_per_type,
        verbose=False,
    )
    proportions = model.fit_transform(Y, X, coords, cell_type_names=cell_type_names)

    # Store results using existing IO function
    result_to_anndata(proportions, adata, cell_type_names, key_added=key_added)

    # Store parameters as metadata
    adata.uns[f"{key_added}_params"] = {
        "sketch_dim": sketch_dim,
        "lambda_spatial": float(model.lambda_used_),
        "n_hvg": n_hvg,
        "n_markers_per_type": n_markers_per_type,
        "n_genes_used": len(gene_names),
        "n_cell_types": len(cell_type_names),
        "cell_type_names": list(cell_type_names),
        "converged": model.info_.get("converged", False),
        "n_iterations": model.info_.get("n_iterations", 0),
    }

    return adata if copy else None
