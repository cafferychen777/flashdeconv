# FlashDeconv

**Fast Linear Algebra for Scalable Hybrid Deconvolution**

FlashDeconv is a high-performance spatial transcriptomics deconvolution method that estimates cell type proportions from spatial gene expression data using single-cell reference signatures.

## Key Features

- **Ultra-fast**: Process 1 million spots in ~10 minutes on CPU
- **Memory-efficient**: Uses randomized sketching to reduce memory by orders of magnitude
- **No GPU required**: Runs efficiently on standard laptops
- **Statistically rigorous**: Handles count data properly via Pearson residuals
- **Spatially-aware**: Graph Laplacian regularization for spatial coherence
- **Platform-robust**: Automatic correction for platform effects between scRNA-seq and ST

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

### 1. Variance Stabilization & Platform Correction
- **Pearson Residuals**: Transforms raw counts to Pearson residuals for variance stabilization, handling the Poisson-Gamma mixture distribution of ST data.
- **Platform Correction**: Estimates and corrects platform-specific gene capture efficiency using a "Pseudo-bulk Ratio" strategy, ensuring compatibility between scRNA-seq references and spatial data.

### 2. Structure-Preserving Sketching
- **CountSketch with Leverage Scores**: Compresses gene dimension (~20,000) to sketch space (~512) using sparse CountSketch.
- **Importance Sampling**: Weighs genes by their leverage scores (biological importance) to ensure rare cell type markers are preserved with high probability, offering theoretical guarantees via the Johnson-Lindenstrauss lemma.

### 3. Spatial Graph Regularized Optimization
- **Graph Laplacian Regularization**: Explicitly models spatial autocorrelation, forcing neighboring spots to have similar cell type compositions to denoise "drop-out" events.
- **Block Coordinate Descent (BCD)**: A Numba-accelerated solver that exploits pre-computed Gram matrices and closed-form updates with non-negativity constraints for extreme speed.

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

## Benchmarks

| Dataset Size | Cell2location | CARD | FlashDeconv |
|--------------|---------------|------|-------------|
| 10K spots | ~1 hour (GPU) | ~30 min | **~30 sec** |
| 100K spots | ~10 hours (GPU) | OOM | **~3 min** |
| 1M spots | >24 hours (GPU) | OOM | **~10 min** |

*Benchmarks on MacBook Air M1 with 8GB RAM (FlashDeconv) vs. NVIDIA A100 GPU (Cell2location)*

## API Reference

### FlashDeconv Class

```python
class FlashDeconv:
    def __init__(
        self,
        sketch_dim=512,
        lambda_spatial="auto",
        rho_sparsity=0.01,
        n_hvg=2000,
        n_markers_per_type=50,
        spatial_method="knn",
        k_neighbors=6,
        max_iter=100,
        tol=1e-4,
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

- `beta_`: Raw cell type abundances
- `proportions_`: Normalized proportions (sum to 1)
- `gene_idx_`: Indices of genes used
- `lambda_used_`: Actual lambda value used
- `theta_`: Estimated overdispersion
- `info_`: Optimization information

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
