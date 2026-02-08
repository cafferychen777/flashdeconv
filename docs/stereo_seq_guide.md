# Stereo-seq Platform Guide

*Understanding how Stereo-seq's technical characteristics affect deconvolution.*

---

## Why Stereo-seq Needs Attention

Stereo-seq differs from Visium in ways that matter for deconvolution:

| Property | Stereo-seq | Visium | Implication |
|:---------|:-----------|:-------|:------------|
| **Resolution** | 500 nm (DNB) | 55 μm (spot) | Bins may contain partial cells |
| **Sparsity** | ~95% zeros | ~85% zeros | Need robust normalization |
| **Layout** | Random DNB | Hexagonal grid | k-NN graph, not fixed neighbors |
| **Capture** | DNA nanoball | Oligo array | Different UMI distributions |

> **First Principle:** *Stereo-seq's high sparsity means most gene-spot entries are zero. Normalization must handle zeros gracefully without amplifying noise.*

---

## The Sparsity Problem

With 95% sparsity, Pearson residual normalization can be unstable:

```
Pearson residual: r = (y - μ) / σ

When y = 0 (95% of data):
  r = -μ/σ ≈ small negative value

When y > 0 (5% of data):
  r = (y - μ)/σ  → can be large if σ is small
```

This creates a bimodal distribution that may not suit NNLS regression.

**Log-CPM** avoids this issue: `log1p(0) = 0`, preserving sparsity structure while stabilizing variance for non-zero values.

---

## Recommended Settings

Based on validation with MOSTA adult mouse brain (38,746 spots, 94.9% sparsity):

```python
fd.tl.deconvolve(
    adata_st,
    adata_ref,
    cell_type_key="cell_type",
    preprocess="log_cpm",   # Handles sparsity well
    n_hvg=2000,             # Sufficient gene coverage
)
```

Other parameters (`lambda_spatial`, `k_neighbors`, `sketch_dim`) have **low sensitivity** on Stereo-seq data. Defaults work well.

---

## Validation Strategy

Without ground truth cell type labels, validate against known biology:

### 1. Anatomical Consistency

Cell type distributions should match known tissue structure:

| Region | Expected Dominant Types |
|:-------|:------------------------|
| Cortex layers | Excitatory neurons |
| White matter / fiber tracts | Oligodendrocytes |
| Meninges | VLMCs, endothelial |
| Hippocampus CA1-3 | Excitatory neurons |

```python
# If you have region annotations
results_df = pd.DataFrame(
    adata_st.obsm['flashdeconv'],
    columns=cell_type_names
)
results_df['region'] = adata_st.obs['annotation'].values

# Check mean proportions by region
results_df.groupby('region')[cell_type_names].mean()
```

### 2. Marker Gene Correlation

Cell type proportions should correlate with known marker expression:

```python
markers = {'Excitatory': 'Slc17a7', 'Astrocytes': 'Aqp4', 'Oligodendrocytes': 'Mbp'}

for ct, gene in markers.items():
    if gene in adata_st.var_names:
        r = np.corrcoef(
            adata_st[:, gene].X.toarray().flatten(),
            adata_st.obsm['flashdeconv'][ct]
        )[0, 1]
        print(f"{ct} vs {gene}: r = {r:.3f}")
```

### 3. Spatial Coherence

Adjacent spots should have similar cell type compositions (biological expectation):

```python
from scipy.spatial.distance import pdist, squareform

# Compute proportion similarity between neighbors
coords = adata_st.obsm['spatial']
props = adata_st.obsm['flashdeconv']

# This should show that nearby spots have correlated proportions
```

---

## Troubleshooting

### Results differ from RCTD

This is expected. Different methods use different:
- **Normalization** (log-CPM vs raw counts vs Pearson)
- **Regression models** (NNLS vs Poisson-lognormal likelihood)
- **Spatial priors** (graph Laplacian vs none)

Compare against biology, not other algorithms. If both methods produce anatomically consistent results, the difference reflects modeling assumptions, not errors.

### One cell type dominates

Check reference balance:

```python
print(adata_ref.obs['cell_type'].value_counts())
```

If one type has 10× more cells, its signature is over-represented. Consider downsampling to ~2000 cells per type.

### Rare cell type not detected

Increase `n_hvg` to include more marker genes:

```python
fd.tl.deconvolve(..., n_hvg=3000)
```

Or verify the reference has sufficient cells (≥500) for that type.

### Unknown cells absorb most proportions

If results show one generic type (e.g., "Unknown", "Unassigned") dominating, this is a reference quality issue, not a platform issue. See [The Unknown Cell Problem](reference_data_guide.md#the-unknown-cell-problem) for the mathematical explanation and solution.

---

## Benchmark Summary

MOSTA adult mouse brain validation:

| Metric | Value |
|:-------|:------|
| Runtime | 1.3–1.6 s |
| Memory | < 2 GB |
| Convergence | 100% |
| Anatomical consistency | High (see validation) |

Parameter sensitivity:

| Parameter | Sensitivity |
|:----------|:------------|
| `preprocess` | High (log_cpm vs pearson: r ≈ 0.72–0.94) |
| `n_hvg` | High (affects gene coverage) |
| `lambda_spatial` | Low (auto-tuning sufficient) |
| `k_neighbors` | Low (4–12 all produce similar results) |

---

## See Also

- [Reference Data Guide](reference_data_guide.md) — Building quality reference signatures
- [README](../README.md) — API reference and quick start
