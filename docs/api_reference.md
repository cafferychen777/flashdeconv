# API Reference

[Back to README](../README.md)

### `flashdeconv.FlashDeconv`

Main class for spatial deconvolution.

```python
from flashdeconv import FlashDeconv

model = FlashDeconv(sketch_dim=512, lambda_spatial="auto")
```

#### Constructor parameters

| Parameter | Type | Default | Description |
|:----------|:-----|:--------|:------------|
| `sketch_dim` | `int` | 512 | Number of CountSketch buckets d. Sets the expected gene weights (default mode) or the randomized sketch dimension (`gene_weighting="countsketch"`). |
| `lambda_spatial` | `float` or `"auto"` | `"auto"` | Spatial regularization strength. `"auto"` tunes based on data scale. |
| `rho_sparsity` | `float` | 0.01 | L1 sparsity penalty (dimensionless fraction, internally scaled). |
| `n_hvg` | `int` | 2000 | Number of highly variable genes to select. |
| `n_markers_per_type` | `int` | 50 | Number of marker genes per cell type. |
| `spatial_method` | `str` | `"knn"` | Graph construction: `"knn"`, `"radius"`, or `"grid"`. |
| `k_neighbors` | `int` | 6 | Number of neighbors for KNN graph. |
| `radius` | `float` or `None` | `None` | Radius for radius-based graph (required when `spatial_method="radius"`). |
| `max_iter` | `int` | 1000 | Maximum BCD solver iterations; the solver stops earlier once the relative change in the abundances falls below `tol`. |
| `tol` | `float` | 1e-4 | Convergence tolerance: largest absolute change of any coefficient relative to the largest coefficient. |
| `preprocess` | `str` | `"log_cpm"` | Preprocessing: `"log_cpm"` (`log1p` of expression normalized to 10,000 counts per spot / cell type; the name is historical), `"pearson"` (uncentered Pearson residuals), or `"raw"` (no transformation). |
| `random_state` | `int`, `RandomState` or `None` | 0 | Random seed for the legacy CountSketch projection. Not used by the deterministic default. |
| `verbose` | `bool` | `False` | Whether to print progress. |
| `gene_weighting` | `str` | `"expected"` | Gene representation. `"expected"`: keep the selected genes and scale gene j by w_j, the exact expectation (over random bucket assignment) of its weight in the column-normalized leverage-weighted CountSketch with `sketch_dim` buckets; deterministic. `"countsketch"`: legacy randomized leverage-weighted CountSketch projection (identical to versions ≤ 0.1.6). A fixed numeric `lambda_spatial` acts on the Gram matrix of the chosen representation, whose scale differs between the two options; `"auto"` adapts to either. |

#### Methods

**`fit(Y, X, coords, cell_type_names=None)`**

Fit the deconvolution model.

| Parameter | Type | Description |
|:----------|:-----|:------------|
| `Y` | ndarray or sparse (N, G) | Spatial transcriptomics count matrix (non-negative, finite). Sparse input is never densified; formats other than CSR/CSC are converted to CSR. |
| `X` | ndarray (K, G) | Reference cell type signature matrix (non-negative, finite; same genes and order as `Y`). |
| `coords` | ndarray (N, 2) or (N, 3) | Spatial coordinates. |
| `cell_type_names` | ndarray (K,), optional | Cell type names. |

Returns `self`. Raises `ValueError` for mismatched shapes, empty inputs, or NaN, infinite or negative values.

**`fit_transform(Y, X, coords, **kwargs)`**

Fit and return cell type proportions. Same parameters as `fit()`. Returns ndarray of shape (N, K).

**`get_cell_type_proportions()`** — Return normalized proportions (N, K).

**`get_abundances()`** — Return unnormalized regression coefficients (N, K).

**`get_dominant_cell_type()`** — Return index of dominant cell type per spot (N,).

**`summary()`** — Return dict with model parameters and fit statistics.

**`compute_uncertainty(alpha=0.05, method="sandwich", fdr_q=0.1)`**

Opt-in, model-based uncertainty (not computed by `fit`). Intervals quantify sampling uncertainty of the estimate under the fitted model. They do not include error from reference–tissue mismatch, cell types missing from the reference, or the log-space mixture approximation, which in benchmarks dominate the error against true cell-type fractions. Do not read them as calibrated confidence intervals for true proportions; use them to compare the stability of estimates across spots and cell types.

- `method="sandwich"` (default): per-spot active-set sandwich (HC0) variance with the full inverse Hessian `(G_AA + λ·n_neighbors·I)⁻¹`, neighbours held fixed, delta method to proportions. Types estimated at zero get a one-sided interval `[0, upper]` from a Newton step on the augmented active set whose gradient includes the spatial penalty (the pull `λ·Σ_neighbours β_jk` of the neighbouring spots), truncated at zero. `method="model"` uses the homoscedastic variance `s²·H_A⁻¹`. `method="laplace_diag"` is the deprecated pre-0.2.0 diagonal-Hessian approximation (under-covers even when the model holds; zero-width intervals for types at zero).
- `detected`: Benjamini–Hochberg calls at level `fdr_q` over all spot × type entries, from one-sided Wald p-values for abundance > 0. FDR control against true cell presence holds only when the reference matches the tissue (realised FDR ≈ 0.03 at q = 0.1 on colorectal Xenium bins with a matched reference; 0.18–0.25 under reference mismatch).
- Cost: about 20 µs per spot on 16 cores for K = 38 types and 422 genes (5.5 s for 2.4 × 10⁵ bins, 21 s for 9.7 × 10⁵ bins).

Returns dict with keys: `entropy`, `residual_ss`, `residual_norm`, `se_prop`, `var_prop`, `ci_lower`, `ci_upper`, `ci_half_width`, `cv`, `z_score`, `p_value`, `detected`, `detection_confident` (alias of `detected`; for `laplace_diag`, the old rule lower bound > 0.01), `mean_ci_width`, `method`, `alpha`, `fdr_q`.

**`bootstrap_uncertainty(n_bootstrap=100, max_iter_boot=None, seed=42, verbose=False, tol=None, alpha=0.05)`**

Poisson bootstrap of the observed counts, `Y* ~ Poisson(Y)`, refitting with the fitted genes, gene weights (or CountSketch matrix in legacy mode), λ and ρ, warm-started at the fitted abundances with the model's `max_iter` and `tol` (defaults). Model-based in the same sense as above. Warns if any refit does not converge. Returns dict with keys: `boot_mean`, `boot_std`, `boot_ci_lower`, `boot_ci_upper` (percentile), `boot_cv`, `n_bootstrap`, `n_converged`, `n_iterations`, `final_change`, `warm_start`.

#### Attributes (after fitting)

| Attribute | Type | Description |
|:----------|:-----|:------------|
| `proportions_` | ndarray (N, K) | Cell type proportions (sum to 1 per spot). |
| `beta_` | ndarray (N, K) | Non-negative regression coefficients before row normalization; not absolute cell counts. |
| `gene_idx_` | ndarray | Indices of genes used for deconvolution. |
| `gene_weights_` | ndarray or `None` | Per-gene weights of the selected genes (`gene_weighting="expected"`); `None` in legacy mode. |
| `lambda_used_` | float | Actual lambda value used (relevant when `lambda_spatial="auto"`). |
| `adjacency_` | sparse (N, N) | Binary spatial adjacency matrix. |
| `Y_repr_`, `X_repr_` | (N, p), (K, p) | Spatial data and reference in the solver's gene representation (p = number of selected genes by default; `sketch_dim` in legacy mode). `Y_sketch_` / `X_sketch_` are deprecated aliases. |
| `info_` | dict | Optimization info: `converged`, `n_iterations`, `final_objective`, `final_change`, `XtX` (Gram matrix), `n_neighbors` (per-spot degree). |

### `flashdeconv.tl.deconvolve`

Scanpy-style entry point. Runs deconvolution and stores results in `adata_st`.

```python
fd.tl.deconvolve(
    adata_st, adata_ref,
    cell_type_key="cell_type",
    sketch_dim=512, lambda_spatial="auto", rho_sparsity=0.01,
    n_hvg=2000, n_markers_per_type=50,
    spatial_method="knn", k_neighbors=6, radius=None,
    max_iter=1000, tol=1e-4, preprocess="log_cpm",
    layer_st=None, layer_ref=None,
    spatial_key="spatial", key_added="flashdeconv",
    random_state=0, gene_weighting="expected", copy=False,
)
```

| Parameter | Type | Default | Description |
|:----------|:-----|:--------|:------------|
| `adata_st` | AnnData | — | Spatial transcriptomics data with coordinates in `.obsm[spatial_key]`. |
| `adata_ref` | AnnData | — | Single-cell reference with cell type labels in `.obs[cell_type_key]`. |
| `cell_type_key` | `str` | `"cell_type"` | Column in `adata_ref.obs` for cell type annotations. |
| `layer_st` | `str` or `None` | `None` | Layer in `adata_st` to use. Uses `.X` if `None`. |
| `layer_ref` | `str` or `None` | `None` | Layer in `adata_ref` to use. Uses `.X` if `None`. |
| `spatial_key` | `str` | `"spatial"` | Key in `adata_st.obsm` for spatial coordinates. |
| `key_added` | `str` | `"flashdeconv"` | Key for storing results. |
| `copy` | `bool` | `False` | If `True`, return a copy instead of modifying in-place. |

All other parameters (`sketch_dim`, `lambda_spatial`, `gene_weighting`, etc.) are forwarded to `FlashDeconv` — see [constructor parameters](#constructor-parameters).

**Stores in `adata_st`:**
- `.obsm[key_added]` — DataFrame of cell type proportions (N x K)
- `.obs[f"{key_added}_dominant"]` — Dominant cell type per spot (Categorical)
- `.uns[f"{key_added}_params"]` — Parameters used for deconvolution

### Reference diagnostics

Opt-in checks for spatial regions that the reference cannot explain, typically because a cell lineage is missing from the reference. They detect missing lineages whose expression cannot be reproduced by a mixture of the reference types; closely related subtypes that the remaining types approximate are not detected.

**`flashdeconv.reference_fit_scores(model, eta=0.01, max_em_iter=10, em_tol=1e-4, pool=True, alpha=0.05, null="auto")`**

Per-bin score of how poorly the fitted reference mixture explains the observed counts, adjusted for UMI depth and composition, with an optional spatially pooled version. Scores are calibrated against the bulk of the section, so bins are flagged relative to the rest of the tissue (`score > 1.645` at `alpha=0.05`).

`null` sets that calibration: `"left_half"` uses the median and the spread of the lower half (the previous default); `"central"` fits the peak of the score distribution (central matching); `"auto"` (default) uses whichever gives the smaller spread. The choice only moves the flag threshold, so bins are ranked identically. `"auto"` matters in sections where many bins are explained *better* than expected, for example deep, homogeneous tumour epithelium. These bins form a long lower tail, which inflates the lower-half spread and suppresses almost all flags. The estimator used, with its centre and scale, is returned in `scores["null"]`.

**`flashdeconv.unexplained_genes(model, bins, gene_names=None, min_count=20, pseudocount=0.5)`**

Genes over-represented in a set of bins (for example, the flagged bins) relative to what the fitted reference mixture predicts.

**`flashdeconv.suggest_missing_types(genes, atlas_profiles, atlas_genes, atlas_names, n_markers=50)`**

Ranks the cell types of a broader atlas by how strongly their markers are among the unexplained genes.

### `flashdeconv.io`

I/O utilities for loading data from AnnData objects.

**`load_spatial_data(adata, layer=None, coord_key="spatial")`**

Extract count matrix, coordinates, and gene names from a spatial AnnData object. Looks for coordinates in `adata.obsm[coord_key]`, then `adata.obsm["X_spatial"]`, then `adata.obs[["x", "y"]]`.

| Parameter | Type | Default | Description |
|:----------|:-----|:--------|:------------|
| `adata` | AnnData | — | Spatial transcriptomics AnnData. |
| `layer` | `str` or `None` | `None` | Layer to use for counts. Uses `.X` if `None`. |
| `coord_key` | `str` | `"spatial"` | Key in `adata.obsm` for coordinates. |

Returns `(Y, coords, gene_names)`.

**`load_reference(adata_ref, cell_type_key="cell_type", layer=None, method="mean")`**

Aggregate single-cell reference into cell type signatures. Raises `ValueError` if any cell lacks a label (filter such cells first).

| Parameter | Type | Default | Description |
|:----------|:-----|:--------|:------------|
| `adata_ref` | AnnData | — | Single-cell reference AnnData. |
| `cell_type_key` | `str` | `"cell_type"` | Column in `adata_ref.obs` for cell type labels. |
| `layer` | `str` or `None` | `None` | Layer to use. Uses `.X` if `None`. |
| `method` | `str` | `"mean"` | Aggregation method: `"mean"` or `"sum"`. |

Returns `(X, cell_type_names, gene_names)`.

**`align_genes(Y, X, genes_spatial, genes_ref)`**

Intersect and align genes between spatial and reference data. Returns `(Y_aligned, X_aligned, common_genes)`.

**`prepare_data(adata_st, adata_ref, cell_type_key="cell_type", spatial_coord_key="spatial", layer_st=None, layer_ref=None)`**

Convenience wrapper combining `load_spatial_data`, `load_reference`, and `align_genes`. Returns `(Y, X, coords, cell_type_names, gene_names)`.

**`result_to_anndata(beta, adata, cell_type_names=None, key_added="flashdeconv")`**

Store deconvolution results in AnnData. Adds `.obsm[key_added]` (DataFrame) and `.obs[f"{key_added}_dominant"]` (Categorical).

### `flashdeconv.utils`

Graph construction and evaluation metrics.

#### Graph construction

**`build_knn_graph(coords, k=6, include_self=False)`**

Build k-nearest neighbor spatial graph from coordinates. Each spot is linked to its `k` nearest other spots (also when coordinates are duplicated) and the graph is symmetrized, so degrees can exceed `k`.

| Parameter | Type | Default | Description |
|:----------|:-----|:--------|:------------|
| `coords` | ndarray (N, 2) or (N, 3) | — | Spatial coordinates. |
| `k` | `int` | 6 | Number of nearest neighbors. |
| `include_self` | `bool` | `False` | Whether to include self-loops. |

Returns `scipy.sparse.csr_matrix` (N, N) binary adjacency matrix.

**`build_radius_graph(coords, radius, include_self=False)`**

Build radius-based neighbor graph. Parameters same as `build_knn_graph` except `radius: float` replaces `k`.

**`coords_to_adjacency(coords, method="knn", k=6, radius=None)`**

Convert coordinates to adjacency matrix. Dispatches to `build_knn_graph`, `build_radius_graph`, or grid-based construction depending on `method`.

#### Evaluation metrics

All evaluation functions take `pred` and `true` as ndarray of shape (N, K).

**`compute_rmse(pred, true, per_cell_type=False)`** — Root mean squared error. Returns float or ndarray (K,) if `per_cell_type=True`.

**`compute_mae(pred, true, per_cell_type=False)`** — Mean absolute error. Returns float or ndarray (K,).

**`compute_correlation(pred, true, method="pearson", per_cell_type=False)`** — Pearson or Spearman correlation. Returns float or ndarray (K,).

**`compute_jsd(pred, true, epsilon=1e-10)`** — Jensen-Shannon divergence per spot. Returns ndarray (N,).

**`evaluate_deconvolution(pred, true, cell_type_names=None)`** — Comprehensive evaluation returning a dict with `overall` metrics (RMSE, MAE, Pearson, Spearman, mean JSD) and `per_cell_type` breakdown.
