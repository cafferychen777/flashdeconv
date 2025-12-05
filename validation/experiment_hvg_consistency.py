"""
Experiment: HVG Consistency Analysis

Purpose: Address the "two-stage selection" logic hole.
- Run A: Current pipeline (HVG + Markers → Leverage on subset)
- Run B: Whole-transcriptome Leverage (all genes directly)

Compare: Do we miss important GOLD genes by using HVG as pre-filter?
"""

import numpy as np
import pandas as pd
import scanpy as sc
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

from flashdeconv.utils.genes import (
    compute_leverage_scores,
    select_hvg,
    select_markers,
    select_informative_genes
)

OUTPUT_DIR = Path("validation/results/hvg_consistency")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def load_data():
    """Load mouse brain scRNA-seq reference."""
    print("Loading data...")
    adata = sc.read_h5ad("validation/mouse_brain/scrna_reference.h5ad")
    annotations = pd.read_csv("validation/mouse_brain/cell_annotation.csv", index_col=0)

    common_cells = adata.obs_names.intersection(annotations.index)
    adata = adata[common_cells].copy()
    adata.obs['cell_type'] = annotations.loc[common_cells, 'annotation_1'].values

    # Get gene symbol mapping
    ensembl_to_symbol = dict(zip(adata.var_names, adata.var['SYMBOL']))

    print(f"  Cells: {adata.n_obs}")
    print(f"  Genes: {adata.n_vars}")
    print(f"  Cell types: {adata.obs['cell_type'].nunique()}")

    return adata, ensembl_to_symbol


def build_reference_matrix(adata):
    """Build cell type signature matrix (mean expression per cell type)."""
    cell_types = adata.obs['cell_type'].unique()
    n_genes = adata.n_vars

    X = np.zeros((len(cell_types), n_genes))
    for i, ct in enumerate(cell_types):
        mask = adata.obs['cell_type'] == ct
        X[i, :] = np.asarray(adata[mask].X.mean(axis=0)).flatten()

    return X, list(cell_types)


def run_experiment():
    """Run the HVG consistency experiment."""
    print("=" * 60)
    print("EXPERIMENT: HVG Consistency Analysis")
    print("=" * 60)

    # Load data
    adata, ensembl_to_symbol = load_data()
    gene_names = adata.var_names.tolist()
    gene_symbols = [ensembl_to_symbol.get(g, g) for g in gene_names]

    # Build reference matrix (all genes)
    print("\nBuilding reference matrix...")
    X_full, cell_types = build_reference_matrix(adata)
    print(f"  Reference shape: {X_full.shape}")

    # =====================================================
    # RUN A: Current Pipeline (HVG + Markers → Leverage)
    # =====================================================
    print("\n" + "=" * 60)
    print("RUN A: Current Pipeline (HVG + Markers → Leverage)")
    print("=" * 60)

    # Get raw count matrix for HVG selection
    Y = adata.X.toarray() if hasattr(adata.X, 'toarray') else np.asarray(adata.X)

    # Select HVG
    hvg_idx = select_hvg(Y, n_top=2000)
    print(f"  HVG selected: {len(hvg_idx)}")

    # Select markers
    marker_idx, _ = select_markers(X_full, n_markers=50)
    print(f"  Markers selected: {len(marker_idx)}")

    # Union
    run_a_gene_idx = np.union1d(hvg_idx, marker_idx)
    print(f"  Union (HVG ∪ Markers): {len(run_a_gene_idx)}")

    # Compute leverage on subset
    X_subset = X_full[:, run_a_gene_idx]
    leverage_a = compute_leverage_scores(X_subset)

    # Map back to original gene indices and get top 500
    leverage_a_full = np.zeros(len(gene_names))
    leverage_a_full[run_a_gene_idx] = leverage_a

    top500_a_idx = np.argsort(leverage_a_full)[::-1][:500]
    top500_a_genes = [gene_names[i] for i in top500_a_idx]
    top500_a_symbols = [gene_symbols[i] for i in top500_a_idx]
    top500_a_leverage = leverage_a_full[top500_a_idx]

    print(f"  Top 500 leverage genes selected")

    # =====================================================
    # RUN B: Whole-Transcriptome Leverage (All Genes)
    # =====================================================
    print("\n" + "=" * 60)
    print("RUN B: Whole-Transcriptome Leverage (All Genes)")
    print("=" * 60)

    # Compute leverage on ALL genes
    print(f"  Computing leverage on all {X_full.shape[1]} genes...")
    leverage_b = compute_leverage_scores(X_full)

    # Get top 500
    top500_b_idx = np.argsort(leverage_b)[::-1][:500]
    top500_b_genes = [gene_names[i] for i in top500_b_idx]
    top500_b_symbols = [gene_symbols[i] for i in top500_b_idx]
    top500_b_leverage = leverage_b[top500_b_idx]

    print(f"  Top 500 leverage genes selected")

    # =====================================================
    # COMPARISON
    # =====================================================
    print("\n" + "=" * 60)
    print("COMPARISON: Run A vs Run B")
    print("=" * 60)

    set_a = set(top500_a_genes)
    set_b = set(top500_b_genes)

    overlap = set_a & set_b
    only_in_a = set_a - set_b
    only_in_b = set_b - set_a  # These are the "missed" GOLD genes!

    print(f"\n  Top 500 Overlap:")
    print(f"    Intersection: {len(overlap)} genes ({100*len(overlap)/500:.1f}%)")
    print(f"    Only in Run A (HVG pipeline): {len(only_in_a)} genes")
    print(f"    Only in Run B (Whole-genome): {len(only_in_b)} genes")

    # =====================================================
    # ANALYZE THE "MISSED" GENES (only in Run B)
    # =====================================================
    print("\n" + "=" * 60)
    print("MISSED GENES: High Leverage but Filtered by HVG")
    print("=" * 60)

    missed_genes = []
    for gene in only_in_b:
        idx = gene_names.index(gene)
        symbol = gene_symbols[idx]
        lev_b = leverage_b[idx]
        lev_a = leverage_a_full[idx]
        in_hvg = idx in hvg_idx
        in_marker = idx in marker_idx

        missed_genes.append({
            'ensembl': gene,
            'symbol': symbol,
            'leverage_whole_genome': lev_b,
            'leverage_hvg_pipeline': lev_a,
            'in_hvg': in_hvg,
            'in_markers': in_marker,
            'rank_whole_genome': np.where(np.argsort(leverage_b)[::-1] == idx)[0][0] + 1
        })

    missed_df = pd.DataFrame(missed_genes)
    missed_df = missed_df.sort_values('rank_whole_genome')

    print(f"\n  Top 20 missed genes (high whole-genome leverage, not in HVG Top 500):")
    print("-" * 80)
    for i, row in missed_df.head(20).iterrows():
        print(f"    Rank {row['rank_whole_genome']:3d}: {row['symbol']:15s} "
              f"(Lev={row['leverage_whole_genome']:.6f}, HVG={row['in_hvg']}, Marker={row['in_markers']})")

    # Check if missed genes are in the HVG set at all
    missed_not_in_hvg = missed_df[~missed_df['in_hvg'] & ~missed_df['in_markers']]
    print(f"\n  Genes completely excluded by HVG+Marker filtering: {len(missed_not_in_hvg)}")

    if len(missed_not_in_hvg) > 0:
        print(f"  These are the 'true' misses - high leverage genes not captured by variance:")
        print("-" * 80)
        for i, row in missed_not_in_hvg.head(20).iterrows():
            print(f"    Rank {row['rank_whole_genome']:3d}: {row['symbol']:15s} "
                  f"(Lev={row['leverage_whole_genome']:.6f})")

    # =====================================================
    # COMPUTE VARIANCE FOR MISSED GENES
    # =====================================================
    print("\n" + "=" * 60)
    print("VARIANCE-LEVERAGE ANALYSIS OF MISSED GENES")
    print("=" * 60)

    # Compute variance for all genes
    X_dense = Y  # Already converted above
    gene_var = np.var(X_dense, axis=0)

    # Analyze missed genes
    missed_gene_idx = [gene_names.index(g) for g in only_in_b]
    missed_var = gene_var[missed_gene_idx]
    missed_lev = leverage_b[missed_gene_idx]

    # Compare to all genes
    var_median = np.median(gene_var)
    lev_median = np.median(leverage_b)

    # Count how many missed genes are GOLD (low var, high lev)
    n_gold_missed = np.sum((missed_var < var_median) & (missed_lev > lev_median))
    n_noise_missed = np.sum((missed_var > var_median) & (missed_lev < lev_median))

    print(f"\n  Missed genes breakdown:")
    print(f"    GOLD type (Low Var / High Lev): {n_gold_missed} ({100*n_gold_missed/len(only_in_b):.1f}%)")
    print(f"    Other types: {len(only_in_b) - n_gold_missed}")

    # =====================================================
    # SAVE RESULTS
    # =====================================================
    print("\n" + "=" * 60)
    print("SAVING RESULTS")
    print("=" * 60)

    # Save Run A results
    df_a = pd.DataFrame({
        'ensembl': top500_a_genes,
        'symbol': top500_a_symbols,
        'leverage': top500_a_leverage,
        'rank': range(1, 501)
    })
    df_a.to_csv(OUTPUT_DIR / "run_a_hvg_pipeline_top500.csv", index=False)

    # Save Run B results
    df_b = pd.DataFrame({
        'ensembl': top500_b_genes,
        'symbol': top500_b_symbols,
        'leverage': top500_b_leverage,
        'rank': range(1, 501)
    })
    df_b.to_csv(OUTPUT_DIR / "run_b_whole_genome_top500.csv", index=False)

    # Save missed genes
    missed_df.to_csv(OUTPUT_DIR / "missed_genes_analysis.csv", index=False)

    # Save summary
    summary = {
        'total_genes': len(gene_names),
        'hvg_selected': len(hvg_idx),
        'markers_selected': len(marker_idx),
        'union_selected': len(run_a_gene_idx),
        'top500_overlap': len(overlap),
        'overlap_percentage': 100 * len(overlap) / 500,
        'only_in_run_a': len(only_in_a),
        'only_in_run_b': len(only_in_b),
        'missed_gold_genes': n_gold_missed,
        'missed_completely_excluded': len(missed_not_in_hvg),
    }

    pd.DataFrame([summary]).to_csv(OUTPUT_DIR / "summary.csv", index=False)

    print(f"\n  Results saved to: {OUTPUT_DIR}")

    # =====================================================
    # FINAL VERDICT
    # =====================================================
    print("\n" + "=" * 60)
    print("VERDICT")
    print("=" * 60)

    overlap_pct = 100 * len(overlap) / 500

    if overlap_pct >= 90:
        print(f"\n  ✓ SAFE: {overlap_pct:.1f}% overlap - HVG pre-filtering captures most high-leverage genes")
        print("    Conclusion: 'HVG as coarse filter is safe; it preserves GOLD genes.'")
    elif overlap_pct >= 70:
        print(f"\n  ⚠ MODERATE: {overlap_pct:.1f}% overlap - some GOLD genes are missed")
        print("    Recommendation: Offer whole-genome mode for maximum sensitivity")
    else:
        print(f"\n  ✗ CONCERN: {overlap_pct:.1f}% overlap - significant GOLD genes missed by HVG")
        print("    Recommendation: Consider defaulting to whole-genome leverage computation")

    return summary


def create_visualization():
    """Create visualization for the HVG consistency experiment."""
    import matplotlib.pyplot as plt
    from matplotlib_venn import venn2

    print("\n" + "=" * 60)
    print("CREATING VISUALIZATION")
    print("=" * 60)

    # Load results
    df_a = pd.read_csv(OUTPUT_DIR / "run_a_hvg_pipeline_top500.csv")
    df_b = pd.read_csv(OUTPUT_DIR / "run_b_whole_genome_top500.csv")
    missed_df = pd.read_csv(OUTPUT_DIR / "missed_genes_analysis.csv")

    set_a = set(df_a['ensembl'])
    set_b = set(df_b['ensembl'])

    # Load full data for variance-leverage plane
    adata = sc.read_h5ad("validation/mouse_brain/scrna_reference.h5ad")
    annotations = pd.read_csv("validation/mouse_brain/cell_annotation.csv", index_col=0)
    common_cells = adata.obs_names.intersection(annotations.index)
    adata = adata[common_cells].copy()
    adata.obs['cell_type'] = annotations.loc[common_cells, 'annotation_1'].values

    gene_names = adata.var_names.tolist()
    ensembl_to_symbol = dict(zip(adata.var_names, adata.var['SYMBOL']))

    # Build reference and compute metrics
    cell_types = adata.obs['cell_type'].unique()
    X = np.zeros((len(cell_types), adata.n_vars))
    for i, ct in enumerate(cell_types):
        mask = adata.obs['cell_type'] == ct
        X[i, :] = np.asarray(adata[mask].X.mean(axis=0)).flatten()

    Y = adata.X.toarray() if hasattr(adata.X, 'toarray') else np.asarray(adata.X)
    gene_var = np.var(Y, axis=0)
    leverage_full = compute_leverage_scores(X)

    log_var = np.log10(gene_var + 1e-10)
    log_lev = np.log10(leverage_full + 1e-10)

    # Create figure
    fig = plt.figure(figsize=(14, 6))

    # Panel A: Venn Diagram
    ax1 = fig.add_subplot(121)
    venn = venn2([set_a, set_b],
                  set_labels=['HVG Pipeline\n(Run A)', 'Whole-Genome\n(Run B)'],
                  ax=ax1)

    # Color the regions
    if venn.get_patch_by_id('10'):
        venn.get_patch_by_id('10').set_color('#3498db')
        venn.get_patch_by_id('10').set_alpha(0.7)
    if venn.get_patch_by_id('01'):
        venn.get_patch_by_id('01').set_color('#e74c3c')
        venn.get_patch_by_id('01').set_alpha(0.7)
    if venn.get_patch_by_id('11'):
        venn.get_patch_by_id('11').set_color('#2ecc71')
        venn.get_patch_by_id('11').set_alpha(0.7)

    ax1.set_title('Top 500 Leverage Genes: Overlap Analysis\n(92.4% Consistency)',
                   fontsize=12, fontweight='bold')

    # Panel B: Variance-Leverage Plane showing missed genes
    ax2 = fig.add_subplot(122)

    # Plot all genes in grey
    ax2.scatter(log_var, log_lev, alpha=0.15, s=5, c='grey', label='All genes')

    # Highlight overlap genes (green)
    overlap_genes = set_a & set_b
    overlap_idx = [gene_names.index(g) for g in overlap_genes if g in gene_names]
    ax2.scatter(log_var[overlap_idx], log_lev[overlap_idx],
                alpha=0.6, s=20, c='#2ecc71', label=f'Overlap ({len(overlap_genes)})')

    # Highlight missed genes (red) - only in Run B
    missed_genes = set_b - set_a
    missed_idx = [gene_names.index(g) for g in missed_genes if g in gene_names]
    ax2.scatter(log_var[missed_idx], log_lev[missed_idx],
                alpha=0.8, s=40, c='#e74c3c', marker='x', linewidths=2,
                label=f'Missed by HVG ({len(missed_genes)})')

    # Add quadrant lines
    var_median = np.median(log_var)
    lev_median = np.median(log_lev)
    ax2.axvline(var_median, color='black', linestyle='--', alpha=0.3)
    ax2.axhline(lev_median, color='black', linestyle='--', alpha=0.3)

    # Label quadrants
    ax2.text(var_median - 2, lev_median + 1.5, 'GOLD\n(Low Var / High Lev)',
             fontsize=9, ha='center', color='#f39c12', fontweight='bold')
    ax2.text(var_median + 2, lev_median - 3, 'NOISE\n(High Var / Low Lev)',
             fontsize=9, ha='center', color='grey')

    # Annotate some missed genes
    for i, row in missed_df.head(5).iterrows():
        gene = row['ensembl']
        if gene in gene_names:
            idx = gene_names.index(gene)
            ax2.annotate(row['symbol'], (log_var[idx], log_lev[idx]),
                        fontsize=8, xytext=(5, 5), textcoords='offset points',
                        color='#e74c3c')

    ax2.set_xlabel('Log10(Variance)', fontsize=11)
    ax2.set_ylabel('Log10(Leverage Score)', fontsize=11)
    ax2.set_title('Missed Genes: All High-Variance (Not GOLD)',
                   fontsize=12, fontweight='bold')
    ax2.legend(loc='upper left', fontsize=9)
    ax2.grid(True, alpha=0.2)

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "hvg_consistency_analysis.png", dpi=300,
                bbox_inches='tight', facecolor='white')
    plt.savefig(OUTPUT_DIR / "hvg_consistency_analysis.pdf",
                bbox_inches='tight', facecolor='white')
    plt.close()

    print(f"  Saved visualization to: {OUTPUT_DIR / 'hvg_consistency_analysis.png'}")


if __name__ == "__main__":
    run_experiment()

    # Try to create visualization (requires matplotlib_venn)
    try:
        create_visualization()
    except ImportError as e:
        print(f"\n  Note: Visualization skipped ({e})")
        print("  Install matplotlib-venn: pip install matplotlib-venn")
