"""
Experiment 2 v2: Unbiased Enrichment Analysis - Improved Approach

Key insight: Instead of Fisher's exact test on fixed marker sets,
directly measure which cell types each gene is most specific to.

Method:
1. For each gene, compute specificity to each cell type
2. Assign each gene to its "best" cell type (highest specificity)
3. Compare: Do GOLD genes preferentially mark rare cell types?
4. Compare: What about NOISE genes?
"""

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

from flashdeconv.utils.genes import compute_leverage_scores

OUTPUT_DIR = Path("validation/results/unbiased_enrichment")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def load_data():
    """Load mouse brain scRNA-seq reference."""
    print("Loading data...")
    adata = sc.read_h5ad("validation/mouse_brain/scrna_reference.h5ad")
    annotations = pd.read_csv("validation/mouse_brain/cell_annotation.csv", index_col=0)

    common_cells = adata.obs_names.intersection(annotations.index)
    adata = adata[common_cells].copy()
    adata.obs['cell_type'] = annotations.loc[common_cells, 'annotation_1'].values

    ensembl_to_symbol = dict(zip(adata.var_names, adata.var['SYMBOL']))

    print(f"  Cells: {adata.n_obs}")
    print(f"  Genes: {adata.n_vars}")
    print(f"  Cell types: {adata.obs['cell_type'].nunique()}")

    return adata, ensembl_to_symbol


def compute_gene_specificity(adata):
    """
    Compute cell type specificity for each gene.
    Returns: DataFrame with gene, best_cell_type, specificity_score
    """
    print("\nComputing gene specificity scores...")

    # Build mean expression per cell type
    cell_types = list(adata.obs['cell_type'].unique())
    n_genes = adata.n_vars

    X_mean = np.zeros((len(cell_types), n_genes))
    for i, ct in enumerate(cell_types):
        mask = adata.obs['cell_type'] == ct
        X_mean[i, :] = np.asarray(adata[mask].X.mean(axis=0)).flatten()

    # Normalize to get expression fraction per cell type
    X_norm = X_mean / (X_mean.sum(axis=0, keepdims=True) + 1e-10)

    # For each gene, find best cell type and compute specificity
    best_ct_idx = np.argmax(X_norm, axis=0)
    best_ct = [cell_types[i] for i in best_ct_idx]

    # Specificity = max fraction / mean of other fractions
    max_frac = np.max(X_norm, axis=0)
    mean_others = (X_norm.sum(axis=0) - max_frac) / (len(cell_types) - 1 + 1e-10)
    specificity = max_frac / (mean_others + 1e-10)

    gene_spec_df = pd.DataFrame({
        'ensembl': adata.var_names.tolist(),
        'best_cell_type': best_ct,
        'specificity_score': specificity,
        'max_fraction': max_frac,
    })

    return gene_spec_df, cell_types


def compute_gene_quadrants(adata, ensembl_to_symbol):
    """Compute variance-leverage quadrants for all genes."""
    print("\nComputing variance-leverage quadrants...")

    # Build reference matrix
    cell_types = adata.obs['cell_type'].unique()
    X = np.zeros((len(cell_types), adata.n_vars))
    for i, ct in enumerate(cell_types):
        mask = adata.obs['cell_type'] == ct
        X[i, :] = np.asarray(adata[mask].X.mean(axis=0)).flatten()

    # Compute variance
    Y = adata.X.toarray() if hasattr(adata.X, 'toarray') else np.asarray(adata.X)
    gene_var = np.var(Y, axis=0)

    # Compute leverage
    leverage = compute_leverage_scores(X)

    # Log transform
    log_var = np.log10(gene_var + 1e-10)
    log_lev = np.log10(leverage + 1e-10)

    # Define quadrants
    var_thresh = np.median(log_var)
    lev_thresh = np.median(log_lev)

    gene_names = adata.var_names.tolist()
    gene_symbols = [ensembl_to_symbol.get(g, g) for g in gene_names]

    quadrant_df = pd.DataFrame({
        'ensembl': gene_names,
        'symbol': gene_symbols,
        'variance': gene_var,
        'leverage': leverage,
        'log_variance': log_var,
        'log_leverage': log_lev,
    })

    # Assign quadrants
    quadrant_df['quadrant'] = 'Other'
    quadrant_df.loc[(log_var < var_thresh) & (log_lev >= lev_thresh), 'quadrant'] = 'GOLD'
    quadrant_df.loc[(log_var >= var_thresh) & (log_lev < lev_thresh), 'quadrant'] = 'NOISE'
    quadrant_df.loc[(log_var >= var_thresh) & (log_lev >= lev_thresh), 'quadrant'] = 'High Var / High Lev'
    quadrant_df.loc[(log_var < var_thresh) & (log_lev < lev_thresh), 'quadrant'] = 'Low Var / Low Lev'

    return quadrant_df


def analyze_quadrant_cell_types(quadrant_df, gene_spec_df, adata):
    """
    Analyze which cell types each quadrant's genes are specific to.
    """
    print("\nAnalyzing cell type specificity by quadrant...")

    # Get cell type abundances
    ct_counts = adata.obs['cell_type'].value_counts()
    ct_fractions = ct_counts / ct_counts.sum()

    # Merge quadrant and specificity info
    merged = quadrant_df.merge(gene_spec_df, on='ensembl')

    results = []

    for quadrant in ['GOLD', 'NOISE', 'High Var / High Lev', 'Low Var / Low Lev']:
        subset = merged[merged['quadrant'] == quadrant]

        if len(subset) == 0:
            continue

        # Count genes per cell type
        ct_gene_counts = subset['best_cell_type'].value_counts()

        for ct in ct_gene_counts.index:
            n_genes = ct_gene_counts[ct]
            abundance = ct_fractions.get(ct, 0)

            # Expected genes if uniform distribution
            expected = len(subset) / len(ct_fractions)
            enrichment = n_genes / expected if expected > 0 else 0

            results.append({
                'quadrant': quadrant,
                'cell_type': ct,
                'n_genes': n_genes,
                'pct_of_quadrant': 100 * n_genes / len(subset),
                'cell_type_abundance': abundance * 100,
                'enrichment_ratio': enrichment,
            })

    return pd.DataFrame(results)


def create_visualization(results_df, quadrant_df, gene_spec_df, adata):
    """Create visualization showing cell type specificity by quadrant."""
    import matplotlib.pyplot as plt

    print("\nCreating visualization...")

    # Get cell type abundances
    ct_counts = adata.obs['cell_type'].value_counts()
    ct_fractions = ct_counts / ct_counts.sum()

    # Merge data
    merged = quadrant_df.merge(gene_spec_df, on='ensembl')

    fig = plt.figure(figsize=(16, 10))

    # Panel A: Cell type abundance distribution for GOLD vs NOISE genes
    ax1 = fig.add_subplot(221)

    gold_genes = merged[merged['quadrant'] == 'GOLD']
    noise_genes = merged[merged['quadrant'] == 'NOISE']

    # Get abundances of cell types that GOLD genes mark
    gold_ct_abundances = [ct_fractions.get(ct, 0) * 100 for ct in gold_genes['best_cell_type']]
    noise_ct_abundances = [ct_fractions.get(ct, 0) * 100 for ct in noise_genes['best_cell_type']]

    # Box plot
    data = [gold_ct_abundances, noise_ct_abundances]
    bp = ax1.boxplot(data, labels=['GOLD\n(Low Var / High Lev)', 'NOISE\n(High Var / Low Lev)'],
                     patch_artist=True)
    bp['boxes'][0].set_facecolor('#f39c12')
    bp['boxes'][1].set_facecolor('#95a5a6')

    # Statistical test
    stat, pval = stats.mannwhitneyu(gold_ct_abundances, noise_ct_abundances, alternative='less')
    ax1.set_ylabel('Cell Type Abundance (%)', fontsize=11)
    ax1.set_title(f'Cell Type Abundance of Marked Populations\n(Mann-Whitney p = {pval:.2e})',
                  fontsize=12, fontweight='bold')
    ax1.set_yscale('log')

    # Panel B: Histogram comparison
    ax2 = fig.add_subplot(222)

    bins = np.logspace(-2, 2, 30)
    ax2.hist(gold_ct_abundances, bins=bins, alpha=0.7, label=f'GOLD (n={len(gold_genes)})',
             color='#f39c12', density=True)
    ax2.hist(noise_ct_abundances, bins=bins, alpha=0.7, label=f'NOISE (n={len(noise_genes)})',
             color='#95a5a6', density=True)

    ax2.axvline(np.median(gold_ct_abundances), color='#d68910', linestyle='--',
                label=f'GOLD median: {np.median(gold_ct_abundances):.2f}%')
    ax2.axvline(np.median(noise_ct_abundances), color='#7f8c8d', linestyle='--',
                label=f'NOISE median: {np.median(noise_ct_abundances):.2f}%')

    ax2.set_xlabel('Cell Type Abundance (%)', fontsize=11)
    ax2.set_ylabel('Density', fontsize=11)
    ax2.set_xscale('log')
    ax2.legend(fontsize=9)
    ax2.set_title('Distribution of Cell Type Abundances', fontsize=12, fontweight='bold')

    # Panel C: Top cell types for GOLD genes
    ax3 = fig.add_subplot(223)

    gold_ct_counts = gold_genes['best_cell_type'].value_counts().head(15)
    y_pos = range(len(gold_ct_counts))

    colors = ['#f39c12' if ct_fractions.get(ct, 0) < 0.02 else '#d4a056'
              for ct in gold_ct_counts.index]

    bars = ax3.barh(y_pos, gold_ct_counts.values, color=colors)
    ax3.set_yticks(y_pos)

    labels = [f"{ct[:18]}{'...' if len(ct) > 18 else ''} ({ct_fractions.get(ct, 0)*100:.1f}%)"
              for ct in gold_ct_counts.index]
    ax3.set_yticklabels(labels, fontsize=9)
    ax3.set_xlabel('Number of GOLD Genes', fontsize=11)
    ax3.set_title('Top Cell Types Marked by GOLD Genes\n(Dark = rare <2%)',
                  fontsize=12, fontweight='bold')
    ax3.invert_yaxis()

    # Panel D: Top cell types for NOISE genes
    ax4 = fig.add_subplot(224)

    noise_ct_counts = noise_genes['best_cell_type'].value_counts().head(15)
    y_pos = range(len(noise_ct_counts))

    colors = ['#5d6d7e' if ct_fractions.get(ct, 0) >= 0.02 else '#95a5a6'
              for ct in noise_ct_counts.index]

    bars = ax4.barh(y_pos, noise_ct_counts.values, color=colors)
    ax4.set_yticks(y_pos)

    labels = [f"{ct[:18]}{'...' if len(ct) > 18 else ''} ({ct_fractions.get(ct, 0)*100:.1f}%)"
              for ct in noise_ct_counts.index]
    ax4.set_yticklabels(labels, fontsize=9)
    ax4.set_xlabel('Number of NOISE Genes', fontsize=11)
    ax4.set_title('Top Cell Types Marked by NOISE Genes\n(Dark = abundant ≥2%)',
                  fontsize=12, fontweight='bold')
    ax4.invert_yaxis()

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "unbiased_enrichment_v2.png", dpi=300,
                bbox_inches='tight', facecolor='white')
    plt.savefig(OUTPUT_DIR / "unbiased_enrichment_v2.pdf",
                bbox_inches='tight', facecolor='white')
    plt.close()

    print(f"  Saved: {OUTPUT_DIR / 'unbiased_enrichment_v2.png'}")

    # Return statistics
    return {
        'gold_median_abundance': np.median(gold_ct_abundances),
        'noise_median_abundance': np.median(noise_ct_abundances),
        'mannwhitney_pvalue': pval,
        'gold_n_rare_ct': sum(1 for a in gold_ct_abundances if a < 2),
        'noise_n_rare_ct': sum(1 for a in noise_ct_abundances if a < 2),
        'gold_pct_rare': 100 * sum(1 for a in gold_ct_abundances if a < 2) / len(gold_ct_abundances),
        'noise_pct_rare': 100 * sum(1 for a in noise_ct_abundances if a < 2) / len(noise_ct_abundances),
    }


def print_summary(stats):
    """Print summary of analysis."""
    print("\n" + "=" * 60)
    print("SUMMARY: Cell Type Specificity Analysis")
    print("=" * 60)

    print(f"\n  GOLD genes (Low Var / High Lev):")
    print(f"    Median cell type abundance: {stats['gold_median_abundance']:.2f}%")
    print(f"    Genes marking rare cells (<2%): {stats['gold_pct_rare']:.1f}%")

    print(f"\n  NOISE genes (High Var / Low Lev):")
    print(f"    Median cell type abundance: {stats['noise_median_abundance']:.2f}%")
    print(f"    Genes marking rare cells (<2%): {stats['noise_pct_rare']:.1f}%")

    print(f"\n  Statistical Test:")
    print(f"    Mann-Whitney U test p-value: {stats['mannwhitney_pvalue']:.2e}")

    print("\n" + "=" * 60)
    print("KEY INSIGHT")
    print("=" * 60)

    if stats['mannwhitney_pvalue'] < 0.05:
        print("\n  ✓ GOLD genes mark RARER cell types than NOISE genes!")
        print(f"    GOLD median: {stats['gold_median_abundance']:.2f}% abundance")
        print(f"    NOISE median: {stats['noise_median_abundance']:.2f}% abundance")
        print("\n  This proves that low-variance/high-leverage genes are")
        print("  specific to rare cell populations, validating the biological")
        print("  relevance of the GOLD gene class.")
    else:
        print("\n  ~ No significant difference in cell type abundance")


def main():
    print("=" * 60)
    print("EXPERIMENT 2 v2: Cell Type Specificity Analysis")
    print("=" * 60)

    # Load data
    adata, ensembl_to_symbol = load_data()

    # Compute gene specificity
    gene_spec_df, cell_types = compute_gene_specificity(adata)

    # Compute gene quadrants
    quadrant_df = compute_gene_quadrants(adata, ensembl_to_symbol)

    # Analyze cell type specificity by quadrant
    results_df = analyze_quadrant_cell_types(quadrant_df, gene_spec_df, adata)
    results_df.to_csv(OUTPUT_DIR / "cell_type_specificity.csv", index=False)

    # Create visualization
    stats = create_visualization(results_df, quadrant_df, gene_spec_df, adata)

    # Print summary
    print_summary(stats)

    # Save stats
    pd.DataFrame([stats]).to_csv(OUTPUT_DIR / "specificity_stats.csv", index=False)

    print(f"\n  All results saved to: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
