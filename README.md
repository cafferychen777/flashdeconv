# FlashDeconv

**Fast Linear Algebra for Scalable Hybrid Deconvolution**

FlashDeconv is a high-performance spatial transcriptomics deconvolution method that estimates cell type proportions from spatial gene expression data using single-cell reference signatures.

## Key Features

- **Ultra-fast**: Process 1 million spots in ~3 minutes on CPU
- **Memory-efficient**: O(N) linear scaling via structure-preserving sketching
- **No GPU required**: Runs on commodity hardware (32GB RAM sufficient for 1M spots)
- **Statistically rigorous**: Log-CPM normalization with leverage-weighted gene selection
- **Spatially-aware**: Sparse graph Laplacian regularization for spatial coherence
- **Rare cell detection**: Leverage scores prioritize marker genes over high-variance genes

## Installation

```bash
# From source
git clone https://github.com/cafferychen777/flashdeconv.git
cd flashdeconv
pip install -e .

# With development dependencies
pip install -e ".[dev]"

# With scanpy integration
pip install -e ".[scanpy]"
```

## Quick Start

```python
from flashdeconv import FlashDeconv

# Initialize model
model = FlashDeconv(
    sketch_dim=512,       # Sketch dimension (default: 512)
    lambda_spatial="auto", # Spatial regularization (auto-tuned)
    rho_sparsity=0.01,    # Sparsity regularization
)

# Fit and get cell type proportions
proportions = model.fit_transform(Y, X, coords)

# Y: spatial count matrix (n_spots x n_genes)
# X: reference signatures (n_cell_types x n_genes)
# coords: spatial coordinates (n_spots x 2)
```

## With AnnData

```python
from flashdeconv import FlashDeconv
from flashdeconv.io import prepare_data, result_to_anndata

# Prepare data from AnnData objects
Y, X, coords, cell_type_names, gene_names = prepare_data(
    adata_st,           # Spatial AnnData
    adata_ref,          # Single-cell reference AnnData
    cell_type_key="cell_type",
)

# Run deconvolution
model = FlashDeconv(verbose=True)
proportions = model.fit_transform(Y, X, coords, cell_type_names=cell_type_names)

# Store results back in AnnData
adata_st = result_to_anndata(proportions, adata_st, cell_type_names)

# Access results
adata_st.obsm["flashdeconv"]  # Cell type proportions DataFrame
adata_st.obs["flashdeconv_dominant"]  # Dominant cell type per spot
```

## Method Overview

FlashDeconv introduces a three-stage framework:

### 1. Gene Selection & Preprocessing
- **Leverage-weighted Gene Selection**: Selects informative genes using leverage scores that prioritize cell-type-specific markers over high-variance genes, enabling accurate detection of rare cell populations.
- **Log-CPM Normalization**: Default preprocessing that normalizes counts per million and applies log1p transformation for variance stabilization.

### 2. Structure-Preserving Sketching
- **CountSketch with Importance Sampling**: Compresses gene dimension (~20,000) to sketch space (~512) using sparse random projections weighted by leverage scores.
- **Theoretical Guarantees**: Preserves distance relationships via Johnson-Lindenstrauss lemma, ensuring rare cell type markers are retained with high probability.

### 3. Spatial Graph Regularized Optimization
- **Sparse Graph Laplacian**: O(N) memory complexity via k-NN graph construction, enabling million-scale analysis without dense covariance matrices.
- **Block Coordinate Descent (BCD)**: Numba-accelerated solver with closed-form updates and non-negativity constraints for extreme speed.

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `sketch_dim` | 512 | Dimension of sketch space |
| `lambda_spatial` | 5000.0 | Spatial regularization strength (use "auto" for automatic tuning) |
| `rho_sparsity` | 0.01 | L1 sparsity regularization |
| `n_hvg` | 2000 | Number of highly variable genes |
| `n_markers_per_type` | 50 | Markers per cell type |
| `k_neighbors` | 6 | Neighbors for spatial graph |
| `max_iter` | 100 | Maximum BCD iterations |
| `tol` | 1e-4 | Convergence tolerance |
| `preprocess` | "log_cpm" | Preprocessing method: "log_cpm", "pearson", or "raw" |

## Benchmarks

| Dataset Size | FlashDeconv Runtime | Memory |
|--------------|---------------------|--------|
| 10K spots | < 1 sec | < 1 GB |
| 100K spots | ~4 sec | ~2 GB |
| 1M spots | ~3 min | ~21 GB |

*Benchmarks on Apple MacBook Pro M2 Max with 32GB unified memory (no GPU required). FlashDeconv exhibits O(N) linear scaling for both time and memory.*

## API Reference

### FlashDeconv Class

```python
class FlashDeconv:
    def __init__(
        self,
        sketch_dim=512,
        lambda_spatial=5000.0,   # or "auto" for automatic tuning
        rho_sparsity=0.01,
        n_hvg=2000,
        n_markers_per_type=50,
        spatial_method="knn",
        k_neighbors=6,
        max_iter=100,
        tol=1e-4,
        preprocess="log_cpm",    # "log_cpm", "pearson", or "raw"
        random_state=None,
        verbose=False,
    ): ...

    def fit(self, Y, X, coords, gene_names=None, cell_type_names=None) -> self
    def fit_transform(self, Y, X, coords, **kwargs) -> np.ndarray
    def get_cell_type_proportions(self) -> np.ndarray
    def get_abundances(self) -> np.ndarray
    def get_dominant_cell_type(self) -> np.ndarray
    def summary(self) -> dict
```

### Attributes (after fitting)

- `beta_`: Raw cell type abundances (n_spots, n_cell_types)
- `proportions_`: Normalized proportions that sum to 1 (n_spots, n_cell_types)
- `gene_idx_`: Indices of genes used for deconvolution
- `lambda_used_`: Actual spatial regularization value used
- `info_`: Optimization information (converged, n_iterations, final_objective)
- `cell_type_names_`: Names of cell types (if provided)

## Citation

If you use FlashDeconv in your research, please cite:

```bibtex
@article{flashdeconv2024,
  title={FlashDeconv: Fast Linear Algebra for Scalable Hybrid Deconvolution},
  author={FlashDeconv Team},
  journal={bioRxiv},
  year={2024}
}
```

## License

GPL-3.0 License
