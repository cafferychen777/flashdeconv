# FlashDeconv

[![PyPI version](https://img.shields.io/pypi/v/flashdeconv.svg)](https://pypi.org/project/flashdeconv/)
[![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![DOI](https://zenodo.org/badge/1114934837.svg)](https://doi.org/10.5281/zenodo.18109003)

**Atlas-scale spatial deconvolution that preserves rare cell signals.**

FlashDeconv estimates cell type proportions from spatial transcriptomics data (Visium, Visium HD, Stereo-seq). It processes **1 million spots in ~3 minutes** on a laptop while detecting rare cell types that variance-based methods systematically miss.

> **Paper:** Yang, C., Chen, J. & Zhang, X. FlashDeconv enables atlas-scale, multi-resolution spatial deconvolution via structure-preserving sketching. *bioRxiv* (2025). [DOI: 10.64898/2025.12.22.696108](https://doi.org/10.64898/2025.12.22.696108)

---

## Why FlashDeconv?

Traditional deconvolution methods face a fundamental trade-off:

| Method | Speed | Rare Cell Detection | Spatial Modeling |
|:-------|:------|:--------------------|:-----------------|
| Cell2Location | Hours | Good | Implicit |
| RCTD | 20+ min/10K spots | Good | None |
| CARD | Minutes | Moderate | O(N²) dense kernel |
| **FlashDeconv** | **~3 min/1M spots** | **Best (+197% vs HVG)** | **O(N) sparse graph** |

**The core problem:** Variance-based feature selection (PCA, HVG) conflates abundance with biological importance. A gene expressed in 30% of cells dominates variance calculations, while markers of rare but critical populations (cancer stem cells, Tuft cells, endothelial cells) are systematically discarded.

**Our solution:** FlashDeconv uses *leverage scores*—a geometric measure from the reference SVD—to identify genes that define transcriptomically distinct directions, regardless of expression magnitude. This preserves rare cell signals during the 40× gene-space compression.

---

## Quick Start

```bash
pip install flashdeconv
```

```python
import scanpy as sc
import flashdeconv as fd

# Load data
adata_st = sc.read_h5ad("visium_hd.h5ad")
adata_ref = sc.read_h5ad("reference.h5ad")

# Deconvolve (one line)
fd.tl.deconvolve(adata_st, adata_ref, cell_type_key="cell_type")

# Results
adata_st.obsm["flashdeconv"]  # Cell type proportions (DataFrame)
sc.pl.spatial(adata_st, color="flashdeconv_Hepatocyte")
```

---

## Features

- **Linear Scalability** — O(N) time and memory. 1M spots in ~3 min, 100K spots in ~4 sec.
- **Rare Cell Detection** — Leverage-score weighting improves detection by 124–197% over PCA/HVG.
- **Spatial Coherence** — Sparse graph Laplacian regularization without O(N²) cost.
- **Hardware Friendly** — CPU-only. Runs on laptops with 32GB RAM (no GPU required).
- **Visium HD Ready** — Optimized for subcellular-resolution platforms (2–16μm bins).

---

## Algorithm

FlashDeconv solves a graph-regularized non-negative least squares problem in a compressed sketch space:

```
minimize  ½||Y - βX||²_F  +  λ·Tr(β^T L β)  +  ρ||β||₁
   β≥0
```

where Y is spatial expression, X is reference signatures, L is the spatial graph Laplacian, and β is cell type abundances.

![FlashDeconv Framework](https://raw.githubusercontent.com/cafferychen777/flashdeconv/main/figures/figure1.jpeg)

### Three-Stage Pipeline

**1. Gene Selection & Leverage Scoring**
   - Select HVGs from spatial data + markers from reference
   - Compute leverage scores: ℓ_g = Σ_j (V_gj)² from reference SVD
   - Leverage captures *transcriptomic distinctiveness*, not abundance

**2. Structure-Preserving Sketching**
   - Build weighted CountSketch matrix Ω (G × 512)
   - Project: Y_sketch = Y·Ω, X_sketch = X·Ω
   - Leverage weighting ensures rare cell markers survive hash collisions

**3. Spatial-Regularized Optimization**
   - Build k-NN graph from coordinates (default k=6)
   - Solve via Numba-accelerated Block Coordinate Descent
   - Auto-tune λ to balance data fidelity and spatial smoothness

---

## Benchmarks

### Runtime & Memory

| Spots | Time | Memory | Hardware |
|:------|:-----|:-------|:---------|
| 10K | < 1s | < 1GB | MacBook Pro M2 Max |
| 100K | ~4s | ~2GB | 32GB unified memory |
| 1M | ~3min | ~21GB | No GPU required |

### Accuracy (Spotless Benchmark)

| Metric | FlashDeconv | RCTD | Cell2Location |
|:-------|:------------|:-----|:--------------|
| Pearson (56 datasets) | **0.944** | 0.905 | 0.895 |
| Rare cell AUPR | **0.960** | 0.95 | 0.95 |
| Reference stability | **#1** | #2 | — |

---

## Scientific Insights from the Paper

### The Resolution Horizon (8–16μm)

At 8μm resolution, 61.5% of Visium HD bins contain a single dominant cell type. This collapses to 13.3% at 16μm—a 78% reduction in signal purity. Beyond this threshold, cell-cell correlations undergo *sign inversion*: spatially segregated populations (Paneth vs Goblet cells) appear spuriously co-localized due to geometric mixing.

### GOLD vs NOISE Genes

Genes partition into four quadrants by variance and leverage:

| Quadrant | Variance | Leverage | Biological Meaning |
|:---------|:---------|:---------|:-------------------|
| **GOLD** | Low | High | Rare cell markers (vascular, stem) |
| NOISE | High | Low | 6× more unannotated transcripts |

GOLD genes reconstruct clear anatomical structures; NOISE genes show speckle-like patterns.

### Tuft-Stem Chemosensory Niche

At 8μm resolution, FlashDeconv detected focal Tuft cell niches (up to 61% purity) with 16.8-fold enrichment for intestinal stem cells—a spatial relationship invisible at conventional Visium resolution (55μm) and systematically missed by HVG selection.

---

## API Reference

### High-Level (Scanpy-style)

```python
fd.tl.deconvolve(
    adata_st,                    # Spatial AnnData
    adata_ref,                   # Reference AnnData
    cell_type_key="cell_type",   # Column in adata_ref.obs
    key_added="flashdeconv",     # Key for results
    random_state=0,              # Reproducibility
)
```

### Low-Level (NumPy)

```python
from flashdeconv import FlashDeconv

model = FlashDeconv(
    sketch_dim=512,          # Reduced dimension (default: 512)
    lambda_spatial="auto",   # Spatial regularization (auto-tuned)
    rho_sparsity=0.01,       # L1 sparsity penalty
    n_hvg=2000,              # Highly variable genes
    n_markers_per_type=50,   # Markers per cell type
    k_neighbors=6,           # Spatial graph neighbors
    preprocess="log_cpm",    # "log_cpm", "pearson", or "raw"
    random_state=0,
)

proportions = model.fit_transform(Y, X, coords)  # (n_spots, n_cell_types)
```

### Key Parameters

| Parameter | Default | Description |
|:----------|:--------|:------------|
| `sketch_dim` | 512 | Sketch space dimension. Higher = more accurate, slower |
| `lambda_spatial` | "auto" | Spatial regularization strength (auto-tuned to data scale) |
| `rho_sparsity` | 0.01 | L1 sparsity penalty for sparse solutions |
| `n_hvg` | 2000 | Number of highly variable genes |
| `k_neighbors` | 6 | Neighbors for spatial graph (matches Visium hexagonal geometry) |
| `preprocess` | "log_cpm" | Normalization: log_cpm (recommended), pearson, or raw |

### Output Attributes

| Attribute | Shape | Description |
|:----------|:------|:------------|
| `proportions_` | (N, K) | Normalized proportions (sum to 1) |
| `beta_` | (N, K) | Raw abundances |
| `gene_idx_` | (G',) | Indices of genes used (HVG ∪ markers) |
| `info_` | dict | Convergence info (converged, n_iterations, final_objective) |

---

## Input Formats

**Spatial Data:** NumPy array (N×G), SciPy sparse matrix, or AnnData

**Reference:** NumPy array (K×G) or AnnData (auto-aggregated by cell type)

**Coordinates:** NumPy array (N×2) or extracted from `adata.obsm["spatial"]`

---

## Reference Data Quality

Deconvolution accuracy depends critically on reference quality:

| Requirement | Threshold | Consequence if Violated |
|:------------|:----------|:------------------------|
| Cells per type | ≥500 | Unstable signatures |
| Marker fold-change | ≥5× | Cannot distinguish similar types |
| Signature correlation | <0.95 | Algorithm cannot separate |

See [Reference Data Guide](docs/reference_data_guide.md) for detailed recommendations.

---

## Installation

```bash
# From PyPI (recommended)
pip install flashdeconv

# With AnnData integration
pip install flashdeconv[io]

# Development version
git clone https://github.com/cafferychen777/flashdeconv.git
cd flashdeconv && pip install -e ".[dev]"
```

**Requirements:** Python ≥3.9, numpy, scipy, numba. Optional: scanpy, anndata.

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

- **Reproducibility:** [flashdeconv-reproducibility](https://github.com/cafferychen777/flashdeconv-reproducibility)
- **Issues:** [GitHub Issues](https://github.com/cafferychen777/flashdeconv/issues)
- **License:** [BSD-3-Clause](LICENSE)

---

## Acknowledgments

We thank the developers of [Spotless](https://github.com/OmicsML/Spotless-Benchmark), [Cell2Location](https://github.com/BayraktarLab/cell2location), and [RCTD](https://github.com/dmcable/spacexr) for their contributions to spatial transcriptomics benchmarking.
