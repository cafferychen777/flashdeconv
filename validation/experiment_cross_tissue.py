"""
Cross-Tissue Validation: Proving GOLD/NOISE Distinction is Universal

Goal: Demonstrate that the variance-leverage quadrant separation
generalizes across diverse tissues (brain, kidney), not just one dataset.

Key hypothesis: GOLD genes (low variance, high leverage) should mark
rare, structurally important cell types in ANY tissue.

Expected results:
- Brain: Endothelial markers (Cldn5, Rgs5) - already validated
- Kidney: Podocyte markers (Nphs1, Nphs2, Podxl) - rare glomerular cells
"""

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats
from pathlib import Path
import sys
import warnings
warnings.filterwarnings('ignore')

sys.path.insert(0, str(Path(__file__).parent.parent))
from flashdeconv.utils.genes import compute_leverage_scores


OUTPUT_DIR = Path("validation/results/cross_tissue")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def load_kidney_reference():
    """Load kidney scRNA-seq reference from converted files."""
    print("Loading kidney scRNA-seq reference...")

    from scipy.io import mmread

    # Load converted data
    counts_path = "validation/benchmark_data/converted/kidney_counts.mtx"
    genes_path = "validation/benchmark_data/converted/kidney_genes.csv"
    meta_path = "validation/benchmark_data/converted/kidney_metadata.csv"

    # Read counts matrix
    counts = mmread(counts_path).tocsr()
    genes_df = pd.read_csv(genes_path)
    meta_df = pd.read_csv(meta_path)

    print(f"  Loaded: {counts.shape[0]} genes x {counts.shape[1]} cells")
    print(f"  Cell types: {meta_df['cell_type'].nunique()}")

    # Print cell type distribution
    ct_counts = meta_df['cell_type'].value_counts()
    print(f"\n  Cell type distribution:")
    for ct, count in ct_counts.items():
        pct = count / len(meta_df) * 100
        rare_mark = " ← RARE" if pct < 1 else ""
        print(f"    {ct}: {count} ({pct:.1f}%){rare_mark}")

    # Create AnnData object
    adata = sc.AnnData(X=counts.T.tocsr())  # cells x genes
    adata.var_names = genes_df['gene'].values
    adata.obs['cell_type'] = meta_df['cell_type'].values

    # Basic preprocessing
    sc.pp.filter_genes(adata, min_cells=10)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    print(f"\n  After filtering: {adata.n_vars} genes")

    return adata


def load_brain_reference():
    """Load mouse brain scRNA-seq reference (already processed)."""
    print("\nLoading mouse brain scRNA-seq reference...")

    adata = sc.read_h5ad("validation/mouse_brain/scrna_reference.h5ad")
    annotations = pd.read_csv("validation/mouse_brain/cell_annotation.csv", index_col=0)

    common_cells = adata.obs_names.intersection(annotations.index)
    adata = adata[common_cells].copy()
    adata.obs['cell_type'] = annotations.loc[common_cells, 'annotation_1'].values

    print(f"  Loaded: {adata.n_vars} genes x {adata.n_obs} cells")
    print(f"  Cell types: {adata.obs['cell_type'].nunique()}")

    return adata


def compute_quadrants(adata, tissue_name):
    """Compute variance-leverage quadrants for a tissue."""
    print(f"\nComputing quadrants for {tissue_name}...")

    # Build reference matrix (mean expression per cell type)
    cell_types = list(adata.obs['cell_type'].unique())
    n_genes = adata.n_vars

    X_ref = np.zeros((len(cell_types), n_genes))
    for i, ct in enumerate(cell_types):
        mask = adata.obs['cell_type'] == ct
        X_ref[i, :] = np.asarray(adata[mask].X.mean(axis=0)).flatten()

    # Compute variance (across all cells)
    Y = adata.X.toarray() if hasattr(adata.X, 'toarray') else np.asarray(adata.X)
    gene_var = np.var(Y, axis=0)

    # Compute leverage scores
    leverage = compute_leverage_scores(X_ref)

    # Log transform
    log_var = np.log10(gene_var + 1e-10)
    log_lev = np.log10(leverage + 1e-10)

    # Define thresholds (median-based)
    var_thresh = np.median(log_var)
    lev_thresh = np.median(log_lev)

    # Create dataframe
    quadrant_df = pd.DataFrame({
        'gene': adata.var_names.tolist(),
        'variance': gene_var,
        'leverage': leverage,
        'log_variance': log_var,
        'log_leverage': log_lev,
    })

    # Add symbol column if available
    if 'SYMBOL' in adata.var.columns:
        quadrant_df['symbol'] = adata.var['SYMBOL'].values

    # Assign quadrants
    quadrant_df['quadrant'] = 'Other'
    quadrant_df.loc[(log_var < var_thresh) & (log_lev >= lev_thresh), 'quadrant'] = 'GOLD'
    quadrant_df.loc[(log_var >= var_thresh) & (log_lev < lev_thresh), 'quadrant'] = 'NOISE'
    quadrant_df.loc[(log_var >= var_thresh) & (log_lev >= lev_thresh), 'quadrant'] = 'High Var / High Lev'
    quadrant_df.loc[(log_var < var_thresh) & (log_lev < lev_thresh), 'quadrant'] = 'Low Var / Low Lev'

    # Count quadrants
    print(f"  Quadrant counts:")
    for q in ['GOLD', 'NOISE', 'High Var / High Lev', 'Low Var / Low Lev']:
        count = (quadrant_df['quadrant'] == q).sum()
        print(f"    {q}: {count}")

    return quadrant_df, var_thresh, lev_thresh


def check_known_markers(quadrant_df, tissue_name, known_markers, symbol_col=None):
    """Check if known markers appear in GOLD quadrant."""
    print(f"\n  Checking known markers for {tissue_name}:")

    results = []
    for marker, description in known_markers.items():
        # Try different gene name formats
        gene_row = None

        # First check symbol column if available
        if symbol_col is not None and symbol_col in quadrant_df.columns:
            matches = quadrant_df[quadrant_df[symbol_col].str.upper() == marker.upper()]
            if len(matches) > 0:
                gene_row = matches.iloc[0]

        # Otherwise try gene column
        if gene_row is None:
            for gene_name in [marker, marker.upper(), marker.lower(), marker.capitalize()]:
                matches = quadrant_df[quadrant_df['gene'].str.upper() == gene_name.upper()]
                if len(matches) > 0:
                    gene_row = matches.iloc[0]
                    break

        if gene_row is not None:
            quadrant = gene_row['quadrant']
            leverage_rank = (quadrant_df['leverage'] >= gene_row['leverage']).sum()
            variance_rank = (quadrant_df['variance'] >= gene_row['variance']).sum()

            status = "✓" if quadrant == 'GOLD' else "○" if quadrant == 'High Var / High Lev' else "✗"
            print(f"    {status} {marker} ({description}): {quadrant}")
            print(f"       Leverage rank: {leverage_rank}/{len(quadrant_df)}, Variance rank: {variance_rank}/{len(quadrant_df)}")

            results.append({
                'marker': marker,
                'description': description,
                'quadrant': quadrant,
                'leverage': gene_row['leverage'],
                'variance': gene_row['variance'],
                'in_gold': quadrant == 'GOLD'
            })
        else:
            print(f"    ? {marker} ({description}): NOT FOUND in dataset")
            results.append({
                'marker': marker,
                'description': description,
                'quadrant': 'NOT FOUND',
                'leverage': np.nan,
                'variance': np.nan,
                'in_gold': False
            })

    return pd.DataFrame(results)


def run_go_enrichment(quadrant_df, tissue_name):
    """Run GO enrichment on GOLD genes."""
    print(f"\n  Running GO enrichment for {tissue_name} GOLD genes...")

    try:
        from gprofiler import GProfiler

        gold_genes = quadrant_df[quadrant_df['quadrant'] == 'GOLD']['gene'].tolist()

        # Determine organism
        organism = 'mmusculus'  # Mouse for both brain and kidney

        gp = GProfiler(return_dataframe=True)
        results = gp.profile(
            organism=organism,
            query=gold_genes,
            sources=['GO:BP'],
            significance_threshold_method='fdr'
        )

        if len(results) > 0:
            results = results[results['p_value'] < 0.05].head(15)
            print(f"    Top GO terms (p < 0.05):")
            for _, row in results.head(5).iterrows():
                print(f"      - {row['name']}: p={row['p_value']:.2e}")
            return results
        else:
            print("    No significant GO terms found")
            return pd.DataFrame()

    except Exception as e:
        print(f"    GO enrichment failed: {e}")
        return pd.DataFrame()


def create_comparison_figure(brain_df, kidney_df, brain_markers_df, kidney_markers_df,
                             brain_go, kidney_go):
    """Create a clean comparison figure suitable for Nature-style journals."""
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec

    print("\nCreating comparison figure...")

    fig = plt.figure(figsize=(16, 10))
    gs = gridspec.GridSpec(2, 2, hspace=0.35, wspace=0.45)

    colors = {'GOLD': '#f39c12', 'NOISE': '#95a5a6',
              'High Var / High Lev': '#3498db', 'Low Var / Low Lev': '#ecf0f1'}

    # Panel A: Brain quadrant plot
    ax1 = fig.add_subplot(gs[0, 0])

    for quadrant, color in colors.items():
        subset = brain_df[brain_df['quadrant'] == quadrant]
        ax1.scatter(subset['log_variance'], subset['log_leverage'],
                   c=color, alpha=0.4, s=8, label=quadrant, rasterized=True)

    # Highlight brain markers
    for _, row in brain_markers_df.iterrows():
        if row['quadrant'] != 'NOT FOUND':
            if 'symbol' in brain_df.columns:
                gene_data = brain_df[brain_df['symbol'].str.upper() == row['marker'].upper()]
            else:
                gene_data = brain_df[brain_df['gene'].str.upper() == row['marker'].upper()]
            if len(gene_data) > 0:
                ax1.scatter(gene_data['log_variance'], gene_data['log_leverage'],
                           c='#c0392b', s=120, marker='*', zorder=5, edgecolors='white', linewidths=0.5)
                ax1.annotate(row['marker'],
                            (gene_data['log_variance'].values[0] + 0.15,
                             gene_data['log_leverage'].values[0]),
                            fontsize=9, fontweight='bold', color='#c0392b')

    ax1.axvline(brain_df['log_variance'].median(), color='gray', linestyle='--', alpha=0.5, linewidth=1)
    ax1.axhline(brain_df['log_leverage'].median(), color='gray', linestyle='--', alpha=0.5, linewidth=1)
    ax1.set_xlabel('log₁₀(Variance)', fontsize=11)
    ax1.set_ylabel('log₁₀(Leverage)', fontsize=11)
    ax1.set_title('Mouse Brain', fontsize=12, fontweight='bold')
    ax1.legend(loc='lower left', fontsize=8, framealpha=0.9)
    ax1.text(-0.12, 1.05, 'a', transform=ax1.transAxes, fontsize=16, fontweight='bold', va='top')

    # Panel B: Kidney quadrant plot
    ax2 = fig.add_subplot(gs[0, 1])

    for quadrant, color in colors.items():
        subset = kidney_df[kidney_df['quadrant'] == quadrant]
        ax2.scatter(subset['log_variance'], subset['log_leverage'],
                   c=color, alpha=0.4, s=8, rasterized=True)

    # Highlight kidney markers - separate podocyte (in GOLD) from tubular (not in GOLD)
    for _, row in kidney_markers_df.iterrows():
        if row['quadrant'] != 'NOT FOUND':
            gene_data = kidney_df[kidney_df['gene'].str.upper() == row['marker'].upper()]
            if len(gene_data) > 0:
                marker_color = '#c0392b' if row['in_gold'] else '#2c3e50'
                ax2.scatter(gene_data['log_variance'], gene_data['log_leverage'],
                           c=marker_color, s=120, marker='*', zorder=5, edgecolors='white', linewidths=0.5)
                ax2.annotate(row['marker'],
                            (gene_data['log_variance'].values[0] + 0.03,
                             gene_data['log_leverage'].values[0]),
                            fontsize=9, fontweight='bold', color=marker_color)

    ax2.axvline(kidney_df['log_variance'].median(), color='gray', linestyle='--', alpha=0.5, linewidth=1)
    ax2.axhline(kidney_df['log_leverage'].median(), color='gray', linestyle='--', alpha=0.5, linewidth=1)
    ax2.set_xlabel('log₁₀(Variance)', fontsize=11)
    ax2.set_ylabel('log₁₀(Leverage)', fontsize=11)
    ax2.set_title('Mouse Kidney', fontsize=12, fontweight='bold')
    ax2.text(-0.12, 1.05, 'b', transform=ax2.transAxes, fontsize=16, fontweight='bold', va='top')

    # Panel C: Brain GO enrichment
    ax3 = fig.add_subplot(gs[1, 0])
    if len(brain_go) > 0:
        top_terms = brain_go.head(10)
        y_pos = range(len(top_terms))
        bars = ax3.barh(y_pos, -np.log10(top_terms['p_value']), color='#f39c12', edgecolor='#d68910', linewidth=0.5)
        ax3.set_yticks(y_pos)
        labels = [t[:45] + '...' if len(t) > 45 else t for t in top_terms['name']]
        ax3.set_yticklabels(labels, fontsize=9)
        ax3.set_xlabel('-log₁₀(p-value)', fontsize=11)
        ax3.invert_yaxis()
    ax3.set_title('Brain GOLD genes: GO enrichment', fontsize=12, fontweight='bold')
    ax3.text(-0.12, 1.05, 'c', transform=ax3.transAxes, fontsize=16, fontweight='bold', va='top')

    # Panel D: Kidney GO enrichment
    ax4 = fig.add_subplot(gs[1, 1])
    if len(kidney_go) > 0:
        top_terms = kidney_go.head(10)
        y_pos = range(len(top_terms))
        bars = ax4.barh(y_pos, -np.log10(top_terms['p_value']), color='#f39c12', edgecolor='#d68910', linewidth=0.5)
        ax4.set_yticks(y_pos)
        labels = [t[:45] + '...' if len(t) > 45 else t for t in top_terms['name']]
        ax4.set_yticklabels(labels, fontsize=9)
        ax4.set_xlabel('-log₁₀(p-value)', fontsize=11)
        ax4.invert_yaxis()
    ax4.set_title('Kidney GOLD genes: GO enrichment', fontsize=12, fontweight='bold')
    ax4.text(-0.12, 1.05, 'd', transform=ax4.transAxes, fontsize=16, fontweight='bold', va='top')

    plt.savefig(OUTPUT_DIR / "cross_tissue_validation.png", dpi=300,
                bbox_inches='tight', facecolor='white')
    plt.savefig(OUTPUT_DIR / "cross_tissue_validation.pdf",
                bbox_inches='tight', facecolor='white')
    plt.close()

    print(f"  Saved: {OUTPUT_DIR / 'cross_tissue_validation.png'}")


def main():
    print("=" * 60)
    print("CROSS-TISSUE VALIDATION EXPERIMENT")
    print("Proving GOLD/NOISE distinction generalizes across tissues")
    print("=" * 60)

    # Define known markers for each tissue
    brain_markers = {
        'Cldn5': 'Endothelial tight junction',
        'Rgs5': 'Pericyte marker',
        'Kdr': 'VEGF receptor (endothelial)',
        'Cdh5': 'VE-cadherin (endothelial)',
        'Ly6a': 'Endothelial marker',
    }

    kidney_markers = {
        'Nphs1': 'Podocyte (nephrin)',
        'Nphs2': 'Podocyte (podocin)',
        'Podxl': 'Podocyte (podocalyxin)',
        'Synpo': 'Podocyte (synaptopodin)',
        'Wt1': 'Podocyte transcription factor',
        'Pecam1': 'Endothelial (CD31)',
        'Slc12a1': 'Loop of Henle (NKCC2)',
        'Aqp2': 'Collecting duct',
    }

    # Load brain data
    brain_adata = load_brain_reference()
    brain_df, _, _ = compute_quadrants(brain_adata, "Brain")
    brain_markers_df = check_known_markers(brain_df, "Brain", brain_markers, symbol_col='symbol')
    brain_go = run_go_enrichment(brain_df, "Brain")

    # Load kidney data
    kidney_adata = load_kidney_reference()
    kidney_df, _, _ = compute_quadrants(kidney_adata, "Kidney")
    kidney_markers_df = check_known_markers(kidney_df, "Kidney", kidney_markers, symbol_col=None)
    kidney_go = run_go_enrichment(kidney_df, "Kidney")

    # Create comparison figure
    create_comparison_figure(brain_df, kidney_df, brain_markers_df, kidney_markers_df,
                            brain_go, kidney_go)

    # Save results
    brain_df.to_csv(OUTPUT_DIR / "brain_quadrants.csv", index=False)
    kidney_df.to_csv(OUTPUT_DIR / "kidney_quadrants.csv", index=False)
    brain_markers_df.to_csv(OUTPUT_DIR / "brain_markers.csv", index=False)
    kidney_markers_df.to_csv(OUTPUT_DIR / "kidney_markers.csv", index=False)

    if len(brain_go) > 0:
        brain_go.to_csv(OUTPUT_DIR / "brain_go_enrichment.csv", index=False)
    if len(kidney_go) > 0:
        kidney_go.to_csv(OUTPUT_DIR / "kidney_go_enrichment.csv", index=False)

    # Print summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)

    brain_gold_markers = brain_markers_df[brain_markers_df['in_gold']]['marker'].tolist()
    kidney_gold_markers = kidney_markers_df[kidney_markers_df['in_gold']]['marker'].tolist()

    print(f"\nBrain GOLD markers: {brain_gold_markers}")
    print(f"Kidney GOLD markers: {kidney_gold_markers}")

    if len(kidney_gold_markers) >= 2:
        print("\n✓ SUCCESS: Leverage scores identify rare cell markers")
        print("  across multiple tissues, validating cross-tissue generalization!")
    else:
        print("\n○ PARTIAL: Some markers found, further investigation needed")

    print(f"\nAll results saved to: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
