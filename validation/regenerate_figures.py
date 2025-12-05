"""
Regenerate all Deep Dive figures with UNIFIED aspect ratios for clean composition.

Key principle: All panels should have similar width:height ratios
Target: 10:8 (1.25:1) for all panels
"""

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

from flashdeconv.utils.genes import compute_leverage_scores

OUTPUT_DIR = Path("validation/results/leverage_deep_dive")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Unified figure size and DPI
FIGSIZE = (10, 8)  # 1.25:1 aspect ratio for all panels
DPI = 300


def load_data():
    """Load data with gene symbols."""
    print("Loading data...")
    adata = sc.read_h5ad("validation/mouse_brain/scrna_reference.h5ad")
    annotations = pd.read_csv("validation/mouse_brain/cell_annotation.csv", index_col=0)
    common_cells = adata.obs_names.intersection(annotations.index)
    adata = adata[common_cells].copy()
    adata.obs['cell_type'] = annotations.loc[common_cells, 'annotation_1'].values
    ensembl_to_symbol = dict(zip(adata.var_names, adata.var['SYMBOL']))
    return adata, ensembl_to_symbol


def regenerate_panel_a(adata):
    """Panel A: Abundance Invariance - unified aspect ratio."""
    print("\nRegenerating Panel A...")

    df = pd.read_csv(OUTPUT_DIR / "experiment1_abundance_invariance.csv")
    ct_counts = adata.obs['cell_type'].value_counts()
    dominant_ct = ct_counts.index[0]
    for ct in ct_counts.index[::-1]:
        if ct_counts[ct] >= 100:
            rare_ct = ct
            break

    # Single unified figure with 2 subplots side by side
    fig, axes = plt.subplots(1, 2, figsize=FIGSIZE)

    # Left plot
    ax = axes[0]
    ax.plot(df['abundance_dominant_pct'], df['avg_var_rank_dominant'], 'o-',
            label='Variance Rank', color='#3498db', linewidth=2.5, markersize=10)
    ax.plot(df['abundance_dominant_pct'], df['avg_lev_rank_dominant'], 's-',
            label='Leverage Rank', color='#e74c3c', linewidth=2.5, markersize=10)
    ax.set_xlabel(f'{dominant_ct} Abundance (%)', fontsize=12)
    ax.set_ylabel('Average Marker Gene Rank', fontsize=12)
    ax.set_title(f'Dominant Type ({dominant_ct})\nMarker Stability', fontsize=11, fontweight='bold')
    ax.legend(fontsize=10)
    ax.invert_yaxis()
    ax.grid(True, alpha=0.3)

    # Right plot
    ax = axes[1]
    ax.plot(df['abundance_dominant_pct'], df['avg_var_rank_rare'], 'o-',
            label='Variance Rank', color='#3498db', linewidth=2.5, markersize=10)
    ax.plot(df['abundance_dominant_pct'], df['avg_lev_rank_rare'], 's-',
            label='Leverage Rank', color='#e74c3c', linewidth=2.5, markersize=10)
    ax.set_xlabel(f'{dominant_ct} Abundance (%)', fontsize=12)
    ax.set_ylabel('Average Marker Gene Rank', fontsize=12)
    ax.set_title(f'Rare Type ({rare_ct})\nMarker Stability', fontsize=11, fontweight='bold')
    ax.legend(fontsize=10)
    ax.invert_yaxis()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "panel_a_abundance.png", dpi=DPI, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.savefig(OUTPUT_DIR / "panel_a_abundance.pdf", bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.close()
    print(f"  Saved: panel_a_abundance.png ({FIGSIZE})")


def regenerate_panel_b(adata, ensembl_to_symbol):
    """Panel B: Gene Quadrant - unified aspect ratio."""
    print("\nRegenerating Panel B...")

    # Build reference
    cell_types = adata.obs['cell_type'].unique()
    X = np.zeros((len(cell_types), adata.n_vars))
    for i, ct in enumerate(cell_types):
        mask = adata.obs['cell_type'] == ct
        X[i, :] = np.asarray(adata[mask].X.mean(axis=0)).flatten()

    X_dense = adata.X.toarray() if hasattr(adata.X, 'toarray') else np.asarray(adata.X)
    gene_var = np.var(X_dense, axis=0)
    leverage = compute_leverage_scores(X)
    gene_symbols = [ensembl_to_symbol.get(g, g) for g in adata.var_names]

    log_var = np.log10(gene_var + 1e-10)
    log_lev = np.log10(leverage + 1e-10)
    var_thresh = np.median(log_var)
    lev_thresh = np.median(log_lev)

    quadrants = []
    for i in range(len(gene_symbols)):
        if log_var[i] >= var_thresh and log_lev[i] >= lev_thresh:
            quadrants.append('High Var / High Lev')
        elif log_var[i] >= var_thresh and log_lev[i] < lev_thresh:
            quadrants.append('High Var / Low Lev')
        elif log_var[i] < var_thresh and log_lev[i] >= lev_thresh:
            quadrants.append('Low Var / High Lev (GOLD)')
        else:
            quadrants.append('Low Var / Low Lev')

    df = pd.DataFrame({
        'gene_symbol': gene_symbols,
        'log_variance': log_var,
        'log_leverage': log_lev,
        'leverage': leverage,
        'quadrant': quadrants,
    })

    gold_genes = df[df['quadrant'] == 'Low Var / High Lev (GOLD)'].nlargest(5, 'leverage')

    # Unified figure size
    fig, ax = plt.subplots(figsize=FIGSIZE)

    colors = {
        'High Var / High Lev': '#2ecc71',
        'High Var / Low Lev': '#e74c3c',
        'Low Var / High Lev (GOLD)': '#f39c12',
        'Low Var / Low Lev': '#bdc3c7',
    }

    for q in colors:
        mask = df['quadrant'] == q
        ax.scatter(df.loc[mask, 'log_variance'], df.loc[mask, 'log_leverage'],
                   alpha=0.5, s=12, c=colors[q], label=f'{q} ({mask.sum()})')

    ax.axvline(var_thresh, color='black', linestyle='--', alpha=0.5)
    ax.axhline(lev_thresh, color='black', linestyle='--', alpha=0.5)

    # Labels with manual offsets
    offsets = [(15, 10), (15, -15), (-50, 5), (15, 5), (-55, -10)]
    for i, (_, row) in enumerate(gold_genes.iterrows()):
        ax.annotate(row['gene_symbol'],
                    (row['log_variance'], row['log_leverage']),
                    fontsize=10, fontweight='bold',
                    xytext=offsets[i], textcoords='offset points',
                    arrowprops=dict(arrowstyle='-', color='gray', alpha=0.5),
                    bbox=dict(boxstyle='round,pad=0.3', facecolor='#f39c12', alpha=0.4, edgecolor='none'))

    ax.set_xlabel('Log10(Variance)', fontsize=13)
    ax.set_ylabel('Log10(Leverage Score)', fontsize=13)
    ax.set_title('Gene Quadrant Analysis: Variance vs Leverage', fontsize=13, fontweight='bold')
    ax.legend(loc='upper left', fontsize=9, framealpha=0.9)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "panel_b_quadrant.png", dpi=DPI, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.savefig(OUTPUT_DIR / "panel_b_quadrant.pdf", bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.close()
    print(f"  Saved: panel_b_quadrant.png ({FIGSIZE})")


def regenerate_panel_c():
    """Panel C: GO Enrichment - unified aspect ratio."""
    print("\nRegenerating Panel C...")

    go_gold = pd.read_csv(OUTPUT_DIR / "go_enrichment_gold.csv")
    top_gold = go_gold.head(12)  # Fewer terms to fit better

    terms = []
    for t in top_gold['Term']:
        clean = t.split('(GO:')[0].strip()
        if len(clean) > 35:
            clean = clean[:32] + '...'
        terms.append(clean)

    # Unified figure size
    fig, ax = plt.subplots(figsize=FIGSIZE)

    y_pos = np.arange(len(terms))
    pvals = -np.log10(top_gold['Adjusted P-value'])

    ax.barh(y_pos, pvals, color='#f39c12', alpha=0.85, edgecolor='#d68910', linewidth=1.5)
    ax.axvline(x=-np.log10(0.05), color='red', linestyle='--', alpha=0.7, linewidth=2, label='p=0.05')
    ax.axvline(x=-np.log10(0.01), color='darkred', linestyle=':', alpha=0.7, linewidth=2, label='p=0.01')

    ax.set_yticks(y_pos)
    ax.set_yticklabels(terms, fontsize=10)
    ax.set_xlabel('-log10(Adjusted P-value)', fontsize=13)
    ax.set_title('GOLD Genes: GO Biological Process Enrichment', fontsize=13, fontweight='bold')
    ax.legend(loc='lower right', fontsize=10)
    ax.invert_yaxis()

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "panel_c_go.png", dpi=DPI, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.savefig(OUTPUT_DIR / "panel_c_go.pdf", bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.close()
    print(f"  Saved: panel_c_go.png ({FIGSIZE})")


def regenerate_panel_d():
    """Panel D: Spatial Visualization - unified aspect ratio."""
    print("\nRegenerating Panel D...")

    visium = sc.read_h5ad("validation/mouse_brain/visium_mouse_brain.h5ad")
    visium.var_names_make_unique()

    if 'log1p' not in visium.uns:
        sc.pp.normalize_total(visium, target_sum=1e4)
        sc.pp.log1p(visium)

    gold_genes = ['Cldn5', 'Ly6a', 'Rgs5']
    noise_genes = ['Gm42418', 'Gcg', 'Gm26917']

    # Unified figure size - 2x3 grid
    fig, axes = plt.subplots(2, 3, figsize=FIGSIZE)

    # GOLD genes: warm orange colormap (matches gold theme in Panel B/C)
    for i, gene in enumerate(gold_genes):
        if gene in visium.var_names:
            ax = axes[0, i]
            sc.pl.spatial(visium, color=gene, ax=ax, show=False,
                         title=f'GOLD: {gene}', cmap='YlOrRd',
                         colorbar_loc=None, frameon=False)

    # NOISE genes: grey colormap (emphasizes "noise" concept)
    for i, gene in enumerate(noise_genes):
        if gene in visium.var_names:
            ax = axes[1, i]
            sc.pl.spatial(visium, color=gene, ax=ax, show=False,
                         title=f'NOISE: {gene}', cmap='Greys',
                         colorbar_loc=None, frameon=False)

    plt.suptitle('Spatial Expression: GOLD vs NOISE Genes', fontsize=12, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "panel_d_spatial.png", dpi=DPI, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.savefig(OUTPUT_DIR / "panel_d_spatial.pdf", bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.close()
    print(f"  Saved: panel_d_spatial.png ({FIGSIZE})")


def main():
    print("="*60)
    print("REGENERATING ALL PANELS WITH UNIFIED ASPECT RATIO")
    print(f"Target size: {FIGSIZE} (ratio 1.25:1)")
    print("="*60)

    adata, ensembl_to_symbol = load_data()

    regenerate_panel_a(adata)
    regenerate_panel_b(adata, ensembl_to_symbol)
    regenerate_panel_c()
    regenerate_panel_d()

    print("\n" + "="*60)
    print("ALL PANELS REGENERATED WITH UNIFIED ASPECT RATIO")
    print("="*60)


if __name__ == "__main__":
    main()
