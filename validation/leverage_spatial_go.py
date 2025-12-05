"""
Deep Dive Part 2: Spatial Visualization & GO Enrichment Analysis
================================================================

This script provides visual and functional evidence that:
1. GOLD genes (Low Var / High Leverage) show clear anatomical patterns in spatial data
2. NOISE genes (High Var / Low Leverage) show random/artifact patterns
3. GOLD genes are enriched in specific biological pathways (angiogenesis, etc.)
4. NOISE genes have no clear functional enrichment

Requirements:
- scanpy, matplotlib, seaborn
- gseapy (for GO enrichment)
"""

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Output directory
OUTPUT_DIR = Path("validation/results/leverage_deep_dive")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# =============================================================================
# Load Pre-computed Results
# =============================================================================
def load_gene_lists():
    """Load GOLD and NOISE gene lists from previous analysis."""
    gold_df = pd.read_csv(OUTPUT_DIR / "experiment2_gold_genes_symbols.csv")
    noise_df = pd.read_csv(OUTPUT_DIR / "experiment2_noise_genes_symbols.csv")

    print(f"Loaded {len(gold_df)} GOLD genes and {len(noise_df)} NOISE genes")
    print(f"\nTop 10 GOLD genes: {gold_df['symbol'].head(10).tolist()}")
    print(f"\nTop 10 NOISE genes: {noise_df['symbol'].head(10).tolist()}")

    return gold_df, noise_df

def load_visium_data():
    """Load Visium spatial transcriptomics data."""
    print("\nLoading Visium mouse brain data...")
    adata = sc.read_h5ad("validation/mouse_brain/visium_mouse_brain.h5ad")
    adata.var_names_make_unique()
    print(f"  Shape: {adata.shape}")
    print(f"  Has spatial: {'spatial' in adata.obsm}")
    return adata

# =============================================================================
# Experiment: Spatial Visualization
# =============================================================================
def spatial_visualization(visium_adata, gold_df, noise_df, n_genes=3):
    """
    Create spatial expression plots comparing GOLD vs NOISE genes.

    Expected outcome:
    - GOLD genes should show clear anatomical structures (e.g., blood vessels)
    - NOISE genes should show random/salt-pepper noise patterns
    """
    print("\n" + "="*70)
    print("SPATIAL VISUALIZATION: GOLD vs NOISE Genes")
    print("="*70)

    # Get available genes in Visium data
    available_genes = set(visium_adata.var_names)

    # Select top GOLD genes that are available
    gold_genes = []
    for gene in gold_df['symbol']:
        if gene in available_genes:
            gold_genes.append(gene)
        if len(gold_genes) >= n_genes:
            break

    # Select top NOISE genes that are available
    noise_genes = []
    for gene in noise_df['symbol']:
        if gene in available_genes:
            noise_genes.append(gene)
        if len(noise_genes) >= n_genes:
            break

    print(f"\nGOLD genes for visualization: {gold_genes}")
    print(f"NOISE genes for visualization: {noise_genes}")

    if len(gold_genes) == 0 or len(noise_genes) == 0:
        print("ERROR: No matching genes found in Visium data!")
        return None

    # Normalize data for visualization if not already done
    if 'log1p' not in visium_adata.uns:
        sc.pp.normalize_total(visium_adata, target_sum=1e4)
        sc.pp.log1p(visium_adata)

    # Create figure with side-by-side comparison
    n_cols = max(len(gold_genes), len(noise_genes))
    fig, axes = plt.subplots(2, n_cols, figsize=(5*n_cols, 10))

    # Plot GOLD genes (top row)
    for i, gene in enumerate(gold_genes):
        ax = axes[0, i] if n_cols > 1 else axes[0]
        sc.pl.spatial(visium_adata, color=gene, ax=ax, show=False,
                     title=f"GOLD: {gene}", cmap='plasma',
                     colorbar_loc=None, frameon=False)

    # Plot NOISE genes (bottom row)
    for i, gene in enumerate(noise_genes):
        ax = axes[1, i] if n_cols > 1 else axes[1]
        sc.pl.spatial(visium_adata, color=gene, ax=ax, show=False,
                     title=f"NOISE: {gene}", cmap='plasma',
                     colorbar_loc=None, frameon=False)

    # Fill empty axes if any
    for i in range(len(gold_genes), n_cols):
        ax = axes[0, i] if n_cols > 1 else axes[0]
        ax.axis('off')
    for i in range(len(noise_genes), n_cols):
        ax = axes[1, i] if n_cols > 1 else axes[1]
        ax.axis('off')

    # Add row labels
    fig.text(0.02, 0.75, 'GOLD Genes\n(Low Var / High Lev)',
             fontsize=14, fontweight='bold', va='center', rotation=90,
             color='#f39c12')
    fig.text(0.02, 0.25, 'NOISE Genes\n(High Var / Low Lev)',
             fontsize=14, fontweight='bold', va='center', rotation=90,
             color='#e74c3c')

    plt.suptitle('Spatial Expression: GOLD vs NOISE Genes\n'
                 'GOLD genes should show anatomical structure; NOISE genes should look random',
                 fontsize=14, y=1.02)
    plt.tight_layout()

    plt.savefig(OUTPUT_DIR / "spatial_gold_vs_noise.png", dpi=200, bbox_inches='tight')
    plt.savefig(OUTPUT_DIR / "spatial_gold_vs_noise.pdf", bbox_inches='tight')
    print(f"\nSaved: {OUTPUT_DIR / 'spatial_gold_vs_noise.png'}")

    # Also create individual high-resolution plots
    for gene in gold_genes + noise_genes:
        fig, ax = plt.subplots(figsize=(8, 8))
        category = "GOLD" if gene in gold_genes else "NOISE"
        sc.pl.spatial(visium_adata, color=gene, ax=ax, show=False,
                     title=f"{category}: {gene}", cmap='plasma', frameon=False)
        plt.savefig(OUTPUT_DIR / f"spatial_{category.lower()}_{gene}.png",
                   dpi=200, bbox_inches='tight')
        plt.close()

    print(f"Saved individual gene plots")

    return {'gold_genes': gold_genes, 'noise_genes': noise_genes}

# =============================================================================
# Experiment: GO Enrichment Analysis
# =============================================================================
def go_enrichment_analysis(gold_df, noise_df, n_genes=100):
    """
    Perform GO enrichment analysis on GOLD vs NOISE gene sets.

    Expected outcome:
    - GOLD genes: Enriched in specific biological processes (angiogenesis, etc.)
    - NOISE genes: No significant enrichment or housekeeping functions
    """
    print("\n" + "="*70)
    print("GO ENRICHMENT ANALYSIS: GOLD vs NOISE Genes")
    print("="*70)

    try:
        import gseapy as gp
        has_gseapy = True
    except ImportError:
        print("WARNING: gseapy not installed. Saving gene lists for manual analysis.")
        has_gseapy = False

    # Get gene lists (use symbol column)
    gold_genes = gold_df['symbol'].head(n_genes).tolist()
    noise_genes = noise_df['symbol'].head(n_genes).tolist()

    # Filter out Gm genes for GO analysis (they won't have annotations)
    gold_genes_filtered = [g for g in gold_genes if not g.startswith('Gm') and not g.startswith('ENSMUSG')]
    noise_genes_filtered = [g for g in noise_genes if not g.startswith('Gm') and not g.startswith('ENSMUSG')]

    print(f"\nGOLD genes for GO analysis: {len(gold_genes_filtered)} (filtered from {len(gold_genes)})")
    print(f"NOISE genes for GO analysis: {len(noise_genes_filtered)} (filtered from {len(noise_genes)})")

    # Save gene lists
    pd.DataFrame({'gene': gold_genes_filtered}).to_csv(
        OUTPUT_DIR / "go_gold_genes.csv", index=False)
    pd.DataFrame({'gene': noise_genes_filtered}).to_csv(
        OUTPUT_DIR / "go_noise_genes.csv", index=False)
    print(f"\nGene lists saved for manual GO analysis (Enrichr/DAVID)")

    if not has_gseapy:
        print("\nTo perform GO enrichment:")
        print("  1. Go to https://maayanlab.cloud/Enrichr/")
        print("  2. Upload go_gold_genes.csv and go_noise_genes.csv")
        print("  3. Compare GO Biological Process results")
        return None

    # Run Enrichr analysis
    results = {}

    for name, genes in [('GOLD', gold_genes_filtered), ('NOISE', noise_genes_filtered)]:
        print(f"\nRunning Enrichr for {name} genes...")

        if len(genes) < 5:
            print(f"  Skipping: only {len(genes)} genes")
            continue

        try:
            enr = gp.enrichr(
                gene_list=genes,
                gene_sets=['GO_Biological_Process_2021', 'GO_Molecular_Function_2021',
                          'GO_Cellular_Component_2021', 'KEGG_2021_Human'],
                organism='mouse',
                outdir=str(OUTPUT_DIR / f"enrichr_{name.lower()}"),
                cutoff=0.05
            )
            results[name] = enr.results

            # Save results
            enr.results.to_csv(OUTPUT_DIR / f"go_enrichment_{name.lower()}.csv", index=False)
            print(f"  Saved: go_enrichment_{name.lower()}.csv")

            # Show top results
            top_bp = enr.results[enr.results['Gene_set'].str.contains('Biological_Process')]
            if len(top_bp) > 0:
                print(f"\n  Top GO Biological Process terms for {name}:")
                for _, row in top_bp.head(5).iterrows():
                    print(f"    - {row['Term'][:50]}... (p={row['Adjusted P-value']:.2e})")

        except Exception as e:
            print(f"  Error: {e}")
            results[name] = None

    # Create comparison visualization
    if 'GOLD' in results and results['GOLD'] is not None:
        plot_go_comparison(results)

    return results

def plot_go_comparison(results):
    """Create visualization comparing GO enrichment between GOLD and NOISE genes."""

    fig, axes = plt.subplots(1, 2, figsize=(14, 8))

    for idx, (name, color) in enumerate([('GOLD', '#f39c12'), ('NOISE', '#e74c3c')]):
        ax = axes[idx]

        if name not in results or results[name] is None:
            ax.text(0.5, 0.5, f'No results for {name}', ha='center', va='center')
            ax.set_title(f'{name} Genes')
            continue

        # Get top GO BP terms
        df = results[name]
        bp_df = df[df['Gene_set'].str.contains('Biological_Process')].head(10)

        if len(bp_df) == 0:
            ax.text(0.5, 0.5, 'No significant GO terms', ha='center', va='center')
            ax.set_title(f'{name} Genes')
            continue

        # Plot
        y_pos = range(len(bp_df))
        bars = ax.barh(y_pos, -np.log10(bp_df['Adjusted P-value']), color=color, alpha=0.7)

        # Labels
        terms = [t.split('(')[0][:40] for t in bp_df['Term']]
        ax.set_yticks(y_pos)
        ax.set_yticklabels(terms, fontsize=9)
        ax.set_xlabel('-log10(Adjusted P-value)', fontsize=11)
        ax.set_title(f'{name} Genes: GO Biological Process', fontsize=12, fontweight='bold')
        ax.axvline(x=-np.log10(0.05), color='red', linestyle='--', alpha=0.5, label='p=0.05')
        ax.legend()
        ax.invert_yaxis()

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "go_enrichment_comparison.png", dpi=150, bbox_inches='tight')
    plt.savefig(OUTPUT_DIR / "go_enrichment_comparison.pdf", bbox_inches='tight')
    print(f"\nSaved: {OUTPUT_DIR / 'go_enrichment_comparison.png'}")

# =============================================================================
# Additional Analysis: Expression Pattern Correlation
# =============================================================================
def analyze_expression_patterns(visium_adata, gold_df, noise_df, n_genes=10):
    """
    Analyze spatial autocorrelation of GOLD vs NOISE genes.
    Higher autocorrelation = more structured spatial pattern.
    """
    print("\n" + "="*70)
    print("SPATIAL AUTOCORRELATION ANALYSIS")
    print("="*70)

    available_genes = set(visium_adata.var_names)

    # Get genes
    gold_genes = [g for g in gold_df['symbol'].head(n_genes*2) if g in available_genes][:n_genes]
    noise_genes = [g for g in noise_df['symbol'].head(n_genes*2) if g in available_genes][:n_genes]

    if len(gold_genes) < 3 or len(noise_genes) < 3:
        print("Not enough genes found in Visium data")
        return None

    # Calculate coefficient of variation in spatial neighborhoods
    # Higher CV in neighborhoods = more structured pattern

    from scipy.spatial.distance import cdist

    coords = visium_adata.obsm['spatial']

    # Build spatial neighbor matrix (k nearest neighbors)
    k = 6
    dist_matrix = cdist(coords, coords)

    results = []

    for gene, category in [(g, 'GOLD') for g in gold_genes] + [(g, 'NOISE') for g in noise_genes]:
        if gene not in available_genes:
            continue

        expr = visium_adata[:, gene].X.toarray().flatten() if hasattr(visium_adata[:, gene].X, 'toarray') else visium_adata[:, gene].X.flatten()

        # For each spot, compute local variance
        local_vars = []
        for i in range(len(coords)):
            neighbors = np.argsort(dist_matrix[i])[:k+1]  # Include self
            local_expr = expr[neighbors]
            local_vars.append(np.var(local_expr))

        # Global variance
        global_var = np.var(expr)

        # Ratio: structured patterns have low local variance relative to global
        # (similar values locally, different globally)
        mean_local_var = np.mean(local_vars)
        structure_score = global_var / (mean_local_var + 1e-10)

        results.append({
            'gene': gene,
            'category': category,
            'global_variance': global_var,
            'mean_local_variance': mean_local_var,
            'structure_score': structure_score,
            'mean_expression': np.mean(expr),
        })

    df = pd.DataFrame(results)
    df.to_csv(OUTPUT_DIR / "spatial_structure_scores.csv", index=False)

    # Plot comparison
    fig, ax = plt.subplots(figsize=(8, 6))

    gold_scores = df[df['category'] == 'GOLD']['structure_score']
    noise_scores = df[df['category'] == 'NOISE']['structure_score']

    positions = [1, 2]
    bp = ax.boxplot([gold_scores, noise_scores], positions=positions, widths=0.6,
                    patch_artist=True)

    colors = ['#f39c12', '#e74c3c']
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)

    ax.set_xticklabels(['GOLD Genes\n(Low Var / High Lev)', 'NOISE Genes\n(High Var / Low Lev)'])
    ax.set_ylabel('Spatial Structure Score\n(Global Var / Local Var)', fontsize=11)
    ax.set_title('Spatial Pattern Structure: GOLD vs NOISE Genes\n'
                 'Higher score = more organized spatial pattern', fontsize=12)

    # Add individual points
    for i, (scores, color) in enumerate(zip([gold_scores, noise_scores], colors)):
        x = np.random.normal(positions[i], 0.04, size=len(scores))
        ax.scatter(x, scores, alpha=0.6, color=color, s=50, edgecolor='white')

    # Statistical test
    from scipy.stats import mannwhitneyu
    stat, pval = mannwhitneyu(gold_scores, noise_scores, alternative='greater')
    ax.text(0.95, 0.95, f'Mann-Whitney p = {pval:.2e}', transform=ax.transAxes,
            ha='right', va='top', fontsize=10,
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "spatial_structure_comparison.png", dpi=150, bbox_inches='tight')
    plt.savefig(OUTPUT_DIR / "spatial_structure_comparison.pdf", bbox_inches='tight')
    print(f"\nSaved: {OUTPUT_DIR / 'spatial_structure_comparison.png'}")

    print(f"\nSpatial Structure Score (mean ± std):")
    print(f"  GOLD: {gold_scores.mean():.2f} ± {gold_scores.std():.2f}")
    print(f"  NOISE: {noise_scores.mean():.2f} ± {noise_scores.std():.2f}")
    print(f"  Mann-Whitney p-value: {pval:.2e}")

    return df

# =============================================================================
# Main
# =============================================================================
def main():
    print("="*70)
    print("LEVERAGE DEEP DIVE PART 2: SPATIAL & FUNCTIONAL ANALYSIS")
    print("="*70)

    # Load pre-computed gene lists
    gold_df, noise_df = load_gene_lists()

    # Load Visium data
    visium_adata = load_visium_data()

    # 1. Spatial Visualization
    spatial_results = spatial_visualization(visium_adata, gold_df, noise_df, n_genes=3)

    # 2. Spatial Structure Analysis
    structure_df = analyze_expression_patterns(visium_adata, gold_df, noise_df, n_genes=15)

    # 3. GO Enrichment Analysis
    go_results = go_enrichment_analysis(gold_df, noise_df, n_genes=100)

    print("\n" + "="*70)
    print("ANALYSIS COMPLETE")
    print("="*70)
    print(f"\nResults saved to: {OUTPUT_DIR}")
    print("\nKey outputs:")
    print("  - spatial_gold_vs_noise.png: Side-by-side spatial comparison")
    print("  - spatial_structure_comparison.png: Quantitative structure analysis")
    print("  - go_gold_genes.csv / go_noise_genes.csv: Gene lists for GO analysis")
    print("  - go_enrichment_*.csv: GO enrichment results (if gseapy available)")

if __name__ == "__main__":
    main()
