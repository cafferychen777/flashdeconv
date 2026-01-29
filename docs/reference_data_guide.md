# Building High-Quality Reference Data for Spatial Deconvolution

*A practical guide to maximize deconvolution accuracy through better reference data.*

---

## Why Reference Quality Matters

Spatial deconvolution is fundamentally a **regression problem**: we estimate cell type proportions by decomposing spatial gene expression into a linear combination of reference signatures. The quality of these signatures directly determines the upper bound of deconvolution accuracy.

> **First Principle:** *A deconvolution algorithm can only distinguish cell types that are distinguishable in the reference data.*

This guide distills empirical lessons from extensive benchmarking into actionable recommendations.

---

## Quick Checklist

Before running deconvolution, verify your reference data passes these checks:

| Check | Threshold | Command |
|:------|:----------|:--------|
| **No Unknown cells** | 0 cells | `~obs['cell_type'].str.contains('Unknown\|UND\|Unassigned')` |
| **Cells per type** | ≥500 (ideal: ≥2000) | `adata.obs['cell_type'].value_counts()` |
| **Marker expression** | ≥80% cells express core markers | `(adata[mask, marker].X > 0).mean()` |
| **Marker fold-change** | ≥5× vs other types | `type_mean / other_mean` |
| **Signature correlation** | <0.95 between types | `np.corrcoef(signatures)` |

If any check fails, your deconvolution results for that cell type will be unreliable.

---

## The Unknown Cell Problem

> **First Principle:** *Deconvolution assumes each cell type has a distinct signature. Unknown cells violate this assumption—their signature approximates the average of other types, making them match any mixed spot equally well.*

### Why Unknown Absorbs Proportions: A Causal Chain

**Step 1: Unknown ≈ Linear Combination of Others**

Unknown cells are defined by exclusion—they didn't cluster with any specific type. Their mean expression is therefore a weighted average:

```
sig_Unknown ≈ α₁·sig_Type1 + α₂·sig_Type2 + ... + αₖ·sig_TypeK
```

In zebrafish embryo data: R² = 0.96 when fitting Unknown from other types.

**Step 2: This Creates Equivalent Solutions**

Deconvolution solves `min ||Y - Xβ||²` subject to `β ≥ 0`. When one column is a linear combination of others, infinitely many β achieve the same minimum residual:

```
β = [0.30, 0.25, 0.25, 0.20, 0.00]  →  residual = r
β = [0.15, 0.12, 0.12, 0.10, 0.51]  →  residual = r  (same!)
```

**Step 3: NNLS Prefers Unknown**

Among equivalent solutions, NNLS favors Unknown because:
- Unknown signature correlates highest with average spot expression (0.76 vs 0.73 for specific types)
- Spatial spots are mixtures; Unknown is also a mixture
- The algorithm finds a "shortcut": one universal signature instead of decomposing into specific types

### The Residual Paradox

This is the key insight that reveals the problem:

| Condition | UND Proportion | Residual |
|:----------|:---------------|:---------|
| With Unknown | 68% | 53.1 |
| Without Unknown | 0% | 53.3 |
| **Difference** | **−68%** | **+0.3%** |

**Unknown absorbs 68% of proportions but improves residual by only 0.3%.**

This proves Unknown provides no additional explanatory power—it's a mathematical shortcut, not a biological signal. Removing Unknown forces NNLS to find the biologically meaningful solution with nearly identical fit quality.

### Solution

```python
# Always filter before deconvolution
ref = ref[~ref.obs['cell_type'].isin([
    'Unknown', 'UND', 'Unassigned', 'Ambiguous', 'Mixed', 'Doublet'
])].copy()
```

This affects **all regression-based methods** (RCTD, Cell2Location, CARD, SPOTlight, CIBERSORT, etc.)—it's a property of the mathematics, not any specific implementation.

---

## Cell Count Requirements

Signature stability depends directly on sample size:

```
┌─────────────────────────────────────────────────────────────────┐
│  Cells per Type    │  Signature Quality    │  Recommendation    │
├─────────────────────────────────────────────────────────────────┤
│  > 2000            │  Excellent            │  ✓ Use directly    │
│  500 - 2000        │  Good                 │  ✓ Acceptable      │
│  200 - 500         │  Marginal             │  ⚠ Verify markers  │
│  < 200             │  Unstable             │  ✗ Merge or remove │
└─────────────────────────────────────────────────────────────────┘
```

**Why this matters:** With <200 cells, random sampling noise dominates the mean expression profile. The resulting signature captures noise rather than biology.

---

## The Dual-Marker Annotation Strategy

Standard annotation using only positive markers is insufficient. High-quality annotation requires **both positive and negative markers**:

```python
def annotate_cell_type(adata, positive_markers, negative_markers):
    """
    Strict annotation with dual-marker validation.

    A cell is assigned to a type only if it:
    1. Expresses positive markers (identity)
    2. Does NOT express negative markers (exclusion)
    """
    # Score positive markers
    sc.tl.score_genes(adata, positive_markers, score_name='pos_score')

    # Score negative markers (should be low)
    sc.tl.score_genes(adata, negative_markers, score_name='neg_score')

    # Combined score: reward positive, penalize negative
    adata.obs['final_score'] = adata.obs['pos_score'] - 0.5 * adata.obs['neg_score']

    return adata.obs['final_score']
```

### Recommended Markers by Cell Type

| Cell Type | Positive Markers | Negative Markers |
|:----------|:-----------------|:-----------------|
| **Astrocytes** | Aqp4, Aldh1l1, Slc1a2, Gfap, Glul | Mbp, Cx3cr1, Slc17a7 |
| **Oligodendrocytes** | Mbp, Plp1, Mog, Mag, Mobp | Aqp4, Slc17a7, Gad1 |
| **Microglia** | Cx3cr1, P2ry12, Tmem119, Csf1r | Mbp, Aqp4, Slc17a7 |
| **OPCs** | Pdgfra, Cspg4, Sox10, Olig1 | Mbp, Mog, Aqp4 |
| **Excitatory neurons** | Slc17a7, Camk2a, Satb2 | Gad1, Gad2, Mbp |
| **Inhibitory neurons** | Gad1, Gad2, Slc32a1 | Slc17a7, Mbp |
| **Endothelial** | Pecam1, Cldn5, Vwf | Slc17a7, Gad1, Mbp |

---

## Signature Quality Control

### Step 1: Compute Signature Correlation Matrix

High inter-signature correlation (>0.95) indicates the algorithm cannot distinguish those cell types:

```python
import numpy as np
from scipy.stats import pearsonr

def check_signature_separability(adata, cell_type_key):
    """Compute pairwise signature correlations."""
    cell_types = adata.obs[cell_type_key].unique()

    # Compute mean expression per cell type
    signatures = {}
    for ct in cell_types:
        mask = adata.obs[cell_type_key] == ct
        signatures[ct] = np.asarray(adata[mask].X.mean(axis=0)).flatten()

    # Correlation matrix
    print("Signature correlations (flag if > 0.95):")
    for ct1 in cell_types:
        for ct2 in cell_types:
            if ct1 < ct2:
                r, _ = pearsonr(signatures[ct1], signatures[ct2])
                flag = "⚠️" if r > 0.95 else ""
                print(f"  {ct1} vs {ct2}: r={r:.3f} {flag}")
```

### Step 2: Validate Marker Expression

For each cell type, verify its markers are actually expressed:

```python
def validate_markers(adata, cell_type, markers, cell_type_key):
    """Check marker expression in annotated cells."""
    mask = adata.obs[cell_type_key] == cell_type

    print(f"\n{cell_type} marker validation:")
    for marker in markers:
        if marker not in adata.var_names:
            print(f"  {marker}: NOT IN DATA")
            continue

        expr = adata[mask, marker].X
        if hasattr(expr, 'toarray'):
            expr = expr.toarray()

        pct_expressing = (expr > 0).mean() * 100
        mean_in_type = expr.mean()
        mean_in_others = adata[~mask, marker].X.mean()
        fold_change = mean_in_type / (mean_in_others + 1e-6)

        status = "✓" if pct_expressing > 80 and fold_change > 5 else "⚠️"
        print(f"  {status} {marker}: {pct_expressing:.0f}% expressing, FC={fold_change:.1f}x")
```

---

## Handling Problematic Cell Types

### Problem: High Signature Correlation

**Symptom:** Two cell types have correlation >0.95
**Solution:** Merge into a single category or use only marker genes

```python
# Option 1: Merge similar types
adata.obs['cell_type_merged'] = adata.obs['cell_type'].replace({
    'CD4_naive': 'CD4_T',
    'CD4_memory': 'CD4_T',
})

# Option 2: Marker-only signature (reduces correlation)
marker_genes = ['CD4', 'CD8A', 'CD19', 'CD14', ...]  # curated list
adata_markers = adata[:, marker_genes]
```

### Problem: Low Cell Count

**Symptom:** Cell type has <200 cells
**Solution:**
1. **Merge** with related type (e.g., T cell subtypes → T cells)
2. **Exclude** from analysis (better than unreliable estimates)
3. **Augment** from external datasets (same tissue, same species)

### Problem: Marker Genes Not Expressed

**Symptom:** Core markers show <50% expression or FC <2×
**Solution:** Re-annotate using alternative markers or cluster-based refinement

---

## Complete QC Pipeline

```python
def prepare_reference(adata_raw, cell_type_key, marker_dict):
    """
    Full QC pipeline for reference data preparation.

    Parameters
    ----------
    adata_raw : AnnData
        Raw single-cell reference data
    cell_type_key : str
        Column in .obs with cell type annotations
    marker_dict : dict
        {cell_type: {'positive': [...], 'negative': [...]}}

    Returns
    -------
    adata_clean : AnnData
        QC-passed reference data
    """
    adata = adata_raw.copy()

    # Step 1: Filter low-count cell types
    counts = adata.obs[cell_type_key].value_counts()
    valid_types = counts[counts >= 200].index
    adata = adata[adata.obs[cell_type_key].isin(valid_types)]
    print(f"Step 1: Kept {len(valid_types)} cell types with ≥200 cells")

    # Step 2: Re-annotate with strict markers
    scores = {}
    for ct, markers in marker_dict.items():
        pos = [g for g in markers.get('positive', []) if g in adata.var_names]
        neg = [g for g in markers.get('negative', []) if g in adata.var_names]

        if pos:
            sc.tl.score_genes(adata, pos, score_name=f'{ct}_pos')
            scores[ct] = adata.obs[f'{ct}_pos'].values

            if neg:
                sc.tl.score_genes(adata, neg, score_name=f'{ct}_neg')
                scores[ct] = scores[ct] - 0.5 * adata.obs[f'{ct}_neg'].values

    # Step 3: Assign to highest-scoring type (with margin)
    score_df = pd.DataFrame(scores, index=adata.obs_names)
    max_score = score_df.max(axis=1)
    second_score = score_df.apply(lambda x: x.nlargest(2).iloc[1], axis=1)

    confident_mask = (max_score > second_score + 0.5) & (max_score > 0)
    adata.obs['cell_type_strict'] = 'Ambiguous'
    adata.obs.loc[confident_mask, 'cell_type_strict'] = score_df.loc[confident_mask].idxmax(axis=1)

    # Step 4: Remove ambiguous cells
    adata_clean = adata[adata.obs['cell_type_strict'] != 'Ambiguous']
    print(f"Step 4: Kept {adata_clean.n_obs}/{adata.n_obs} cells after filtering ambiguous")

    # Step 5: Verify signature separability
    check_signature_separability(adata_clean, 'cell_type_strict')

    return adata_clean
```

---

## Common Pitfalls

| Pitfall | Consequence | Prevention |
|:--------|:------------|:-----------|
| **Keeping Unknown/Unassigned cells** | Absorbs 60%+ of proportions | Filter before deconvolution |
| **Automatic clustering annotation** | Includes mislabeled cells | Validate with markers |
| **Using different species** | Gene names/expression don't match | Use same species |
| **Ignoring batch effects** | Technical variation in signature | Batch-correct first |
| **Too many subtypes** | Highly correlated signatures | Merge related types |
| **Only positive markers** | Cross-contamination | Add negative markers |

---

## Recommended Reference Data Sources

| Tissue | Recommended Source | Notes |
|:-------|:-------------------|:------|
| **Brain** | Allen Brain Atlas | Well-curated, multiple regions |
| **Liver** | MacParland et al., Nat Commun 2018 | Hepatocyte zonation |
| **Kidney** | Kidney Cell Atlas (Stewart et al.) | Comprehensive |
| **Immune** | Tabula Sapiens / Tabula Muris | Pan-tissue immune |
| **Tumor** | Tissue-matched normal + tumor | Avoid cell line artifacts |

> **General rule:** Prefer published, peer-reviewed references from the same tissue and species as your spatial data.

---

## Summary

**Reference data quality determines deconvolution accuracy.** The algorithm is only as good as the signatures it's given.

Key principles:
1. **Sufficient cells** (≥500 per type) for stable signatures
2. **Dual markers** (positive + negative) for clean annotation
3. **Low correlation** (<0.95) between signatures for distinguishability
4. **Validated expression** (≥80% expressing, FC ≥5×) for each marker

> *"Garbage in, garbage out. But gold in, gold out."*

---

## See Also

- [FlashDeconv README](../README.md) — Quick start and API reference
- [Stereo-seq Guide](stereo_seq_guide.md) — Platform-specific considerations and parameter sensitivity
