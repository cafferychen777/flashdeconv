# FlashDeconv

[![PyPI version](https://img.shields.io/pypi/v/flashdeconv.svg)](https://pypi.org/project/flashdeconv/)
[![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![DOI](https://zenodo.org/badge/1114934837.svg)](https://doi.org/10.5281/zenodo.18109003)

**Spatial deconvolution with linear scalability for atlas-scale data.**

FlashDeconv estimates cell type proportions from spatial transcriptomics data (Visium, Visium HD, Stereo-seq). It is designed for large-scale analyses where computational efficiency is essential, while maintaining attention to low-abundance cell populations through leverage-score-based feature weighting.

> **Paper:** Yang, C., Chen, J. & Zhang, X. FlashDeconv enables atlas-scale, multi-resolution spatial deconvolution via structure-preserving sketching. *bioRxiv* (2025). [DOI: 10.64898/2025.12.22.696108](https://doi.org/10.64898/2025.12.22.696108)

---

## Installation

```bash
pip install flashdeconv
```

For development or additional I/O support, see [Installation Options](#installation-options).

---

## Quick Start

```python
import scanpy as sc
import flashdeconv as fd

# Load data
adata_st = sc.read_h5ad("spatial.h5ad")
adata_ref = sc.read_h5ad("reference.h5ad")

# Deconvolve
fd.tl.deconvolve(adata_st, adata_ref, cell_type_key="cell_type")

# Results stored in adata_st.obsm["flashdeconv"]
sc.pl.spatial(adata_st, color="flashdeconv_Hepatocyte")
```

---

## Overview

Spatial deconvolution methods offer different trade-offs. Probabilistic approaches like Cell2Location and RCTD provide rigorous uncertainty quantification; methods like CARD incorporate spatial structure through dense kernel matrices. FlashDeconv takes a complementary approach, prioritizing computational efficiency for million-scale datasets.

### Design Principles

1. **Linear complexity** — O(N) time and memory through randomized sketching and sparse graph regularization.

2. **Leverage-based feature weighting** — Variance-based selection (PCA, HVG) can underweight markers of low-abundance populations. We use leverage scores from the reference SVD to identify genes that define distinct transcriptomic directions, regardless of expression magnitude.

3. **Sparse spatial regularization** — Graph Laplacian smoothing with O(N) complexity, avoiding the O(N²) cost of dense kernel methods.

---

## Performance

### Scalability

| Spots | Time | Memory |
|:------|:-----|:-------|
| 10,000 | < 1 sec | < 1 GB |
| 100,000 | ~4 sec | ~2 GB |
| 1,000,000 | ~3 min | ~21 GB |

Benchmarked on MacBook Pro M2 Max (32GB unified memory), CPU-only.

### Accuracy

On the [Spotless benchmark](https://github.com/saeyslab/spotless-benchmark):

| Metric | FlashDeconv | RCTD | Cell2Location |
|:-------|:------------|:-----|:--------------|
| Pearson (56 datasets) | 0.944 | 0.905 | 0.895 |

Performance varies by tissue type and experimental conditions. We recommend evaluating on data similar to your use case.

---

## Algorithm

FlashDeconv solves a graph-regularized non-negative least squares problem:

```
minimize  ½‖Y - βX‖²_F + λ·Tr(βᵀLβ) + ρ‖β‖₁,  subject to β ≥ 0
```

where Y is spatial expression, X is reference signatures, L is the graph Laplacian, and β represents cell type abundances.

![FlashDeconv Framework](https://raw.githubusercontent.com/cafferychen777/flashdeconv/main/figures/figure1.jpeg)

**Pipeline:**
1. Select informative genes (HVG ∪ markers) and compute leverage scores
2. Compress gene space via weighted CountSketch (G → 512 dimensions)
3. Construct sparse k-NN spatial graph
4. Solve via block coordinate descent with spatial smoothing

---

## API

### Scanpy-style

```python
fd.tl.deconvolve(
    adata_st,                    # Spatial AnnData
    adata_ref,                   # Reference AnnData
    cell_type_key="cell_type",   # Column in adata_ref.obs
    key_added="flashdeconv",     # Key for results
)
```

### NumPy

```python
from flashdeconv import FlashDeconv

model = FlashDeconv(
    sketch_dim=512,
    lambda_spatial="auto",
    n_hvg=2000,
    k_neighbors=6,
    random_state=0,
)
proportions = model.fit_transform(Y, X, coords)
```

### Parameters

| Parameter | Default | Description |
|:----------|:--------|:------------|
| `sketch_dim` | 512 | Sketch dimension |
| `lambda_spatial` | "auto" | Spatial regularization (auto-tuned) |
| `n_hvg` | 2000 | Highly variable genes |
| `k_neighbors` | 6 | Spatial graph neighbors |
| `preprocess` | "log_cpm" | Normalization method |
| `random_state` | 0 | Random seed for reproducibility |

### Output

| Attribute | Description |
|:----------|:------------|
| `proportions_` | Cell type proportions (N × K), sum to 1 |
| `beta_` | Raw abundances (N × K) |
| `info_` | Convergence statistics |

---

## Input Formats

- **Spatial data:** AnnData, NumPy array (N × G), or SciPy sparse matrix
- **Reference:** AnnData (aggregated by cell type) or NumPy array (K × G)
- **Coordinates:** Extracted from `adata.obsm["spatial"]` or NumPy array (N × 2)

---

## Reference Quality

Deconvolution accuracy depends on reference quality:

| Requirement | Guideline |
|:------------|:----------|
| Cells per type | ≥ 500 recommended |
| Marker fold-change | ≥ 5× for distinguishability |
| Signature correlation | < 0.95 between types |

See [Reference Data Guide](docs/reference_data_guide.md) for details.

---

## Installation Options

```bash
# Standard
pip install flashdeconv

# With AnnData support
pip install flashdeconv[io]

# Development
git clone https://github.com/cafferychen777/flashdeconv.git
cd flashdeconv && pip install -e ".[dev]"
```

**Requirements:** Python ≥ 3.9, numpy, scipy, numba. Optional: scanpy, anndata.

---

## Citation

```bibtex
@article{yang2025flashdeconv,
  title={FlashDeconv enables atlas-scale, multi-resolution spatial deconvolution
         via structure-preserving sketching},
  author={Yang, Chen and Chen, Jun and Zhang, Xianyang},
  journal={bioRxiv},
  year={2025},
  doi={10.64898/2025.12.22.696108}
}
```

---

## Resources

- [Paper reproducibility code](https://github.com/cafferychen777/flashdeconv-reproducibility)
- [Stereo-seq guide](docs/stereo_seq_guide.md) — Platform-specific considerations
- [GitHub Issues](https://github.com/cafferychen777/flashdeconv/issues)
- [BSD-3-Clause License](LICENSE)

---

## Acknowledgments

We thank the developers of [Spotless](https://github.com/saeyslab/spotless-benchmark), [Cell2Location](https://github.com/BayraktarLab/cell2location), [RCTD](https://github.com/dmcable/spacexr), [CARD](https://github.com/YingMa0107/CARD), and other deconvolution methods whose work contributed to this field.
