"""
Benchmark: Leverage-weighted Sketching vs Uniform Sketching
==========================================================

This script tests the TRUE advantage of FlashDeconv:
Using leverage scores as importance weights in the sketching step.

Both runs use the SAME gene set (HVG + markers), but differ in how
the sketching matrix weights the genes.

Experiments:
- Run A: HVG + markers genes, UNIFORM weights in sketching
- Run B: HVG + markers genes, LEVERAGE weights in sketching (FlashDeconv default)

This tests whether leverage-weighted sketching improves rare cell type accuracy.
"""

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from scipy import sparse
from scipy.stats import pearsonr, spearmanr
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

from flashdeconv.utils.genes import select_hvg, compute_leverage_scores, select_markers
from flashdeconv.core.sketching import sketch_data
from flashdeconv.core.solver import bcd_solve, normalize_proportions
from flashdeconv.utils.graph import coords_to_adjacency

# Output directory
OUTPUT_DIR = Path("validation/results/benchmark_sketching_weights")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def load_and_prepare_data():
    """Load data and prepare common gene matrices."""
    print("Loading data...")

    # Load scRNA reference
    scrna = sc.read_h5ad("validation/mouse_brain/scrna_reference.h5ad")
    annotations = pd.read_csv("validation/mouse_brain/cell_annotation.csv", index_col=0)

    common_cells = scrna.obs_names.intersection(annotations.index)
    scrna = scrna[common_cells].copy()
    scrna.obs['cell_type'] = annotations.loc[common_cells, 'annotation_1'].values

    # Load Visium
    visium = sc.read_h5ad("validation/mouse_brain/visium_mouse_brain.h5ad")
    visium.var_names_make_unique()

    # Gene symbol mapping
    if 'SYMBOL' in scrna.var.columns:
        scrna_symbols = scrna.var['SYMBOL'].fillna(scrna.var_names.to_series()).tolist()
    else:
        scrna_symbols = scrna.var_names.tolist()

    # Find common genes
    visium_symbols = set(visium.var_names)
    common_mask = [s in visium_symbols for s in scrna_symbols]
    common_idx = np.where(common_mask)[0]
    common_symbols = [scrna_symbols[i] for i in common_idx]

    # Map to visium indices
    visium_idx_map = {s: i for i, s in enumerate(visium.var_names)}
    visium_idx = [visium_idx_map[s] for s in common_symbols]

    # Subset data
    scrna_common = scrna[:, common_idx].copy()
    visium_common = visium[:, visium_idx].copy()

    print(f"  Common genes: {len(common_symbols)}")

    # Build reference matrix
    cell_types = sorted(scrna_common.obs['cell_type'].unique())
    X = np.zeros((len(cell_types), len(common_symbols)))

    for i, ct in enumerate(cell_types):
        mask = scrna_common.obs['cell_type'] == ct
        if hasattr(scrna_common.X, 'toarray'):
            X[i, :] = scrna_common[mask].X.toarray().mean(axis=0).flatten()
        else:
            X[i, :] = scrna_common[mask].X.mean(axis=0).flatten()

    # Spatial matrix
    if hasattr(visium_common.X, 'toarray'):
        Y = visium_common.X.toarray()
    else:
        Y = np.asarray(visium_common.X)

    coords = visium_common.obsm['spatial']

    return Y, X, coords, cell_types, common_symbols, visium_common


def run_deconv_with_weights(Y, X, coords, gene_idx, leverage_scores,
                            use_leverage_weights=True, sketch_dim=512,
                            lambda_spatial=5000, verbose=False):
    """
    Run deconvolution with specified weights.

    Parameters
    ----------
    use_leverage_weights : bool
        If True, use leverage scores for sketching weights.
        If False, use uniform weights.
    """
    Y_sub = Y[:, gene_idx]
    X_sub = X[:, gene_idx]

    # Preprocessing
    Y_cpm = Y_sub / (Y_sub.sum(axis=1, keepdims=True) + 1e-10) * 1e4
    X_cpm = X_sub / (X_sub.sum(axis=1, keepdims=True) + 1e-10) * 1e4
    Y_tilde = np.log1p(Y_cpm)
    X_tilde = np.log1p(X_cpm)

    # Choose weights
    if use_leverage_weights:
        weights = leverage_scores
    else:
        weights = np.ones(len(gene_idx)) / len(gene_idx)

    # Sketching
    Y_sketch, X_sketch, Omega = sketch_data(
        Y_tilde, X_tilde,
        sketch_dim=sketch_dim,
        leverage_scores=weights,
        random_state=42
    )

    # Spatial graph
    A = coords_to_adjacency(coords, method='knn', k=6)

    # Solve
    beta, info = bcd_solve(
        Y_sketch, X_sketch, A,
        lambda_=lambda_spatial,
        rho=0.01,
        max_iter=100,
        tol=1e-4,
        verbose=verbose
    )

    proportions = normalize_proportions(beta)
    return proportions, info


def main():
    print("="*70)
    print("BENCHMARK: Leverage-weighted vs Uniform Sketching")
    print("="*70)

    # Load data
    Y, X, coords, cell_types, gene_symbols, visium = load_and_prepare_data()

    print(f"\nData shapes:")
    print(f"  Y (spatial): {Y.shape}")
    print(f"  X (reference): {X.shape}")
    print(f"  Cell types: {len(cell_types)}")

    # Gene selection: HVG + markers (same for both runs)
    print("\nSelecting genes (HVG + markers)...")
    hvg_idx = select_hvg(Y, n_top=2000)
    marker_idx, _ = select_markers(X, n_markers=50)
    gene_idx = np.union1d(hvg_idx, marker_idx)
    print(f"  Selected {len(gene_idx)} genes")

    # Compute leverage scores for selected genes
    leverage_scores = compute_leverage_scores(X[:, gene_idx])

    # Check key markers
    key_markers = ['Cldn5', 'Rgs5', 'Kdr', 'Cdh5', 'Sox17']
    print("\nKey markers in selected genes:")
    for m in key_markers:
        if m in gene_symbols:
            idx = gene_symbols.index(m)
            in_selected = idx in gene_idx
            if in_selected:
                local_idx = np.where(gene_idx == idx)[0][0]
                lev = leverage_scores[local_idx]
                print(f"  {m}: IN (leverage={lev:.2e})")
            else:
                print(f"  {m}: NOT IN")

    # =========================================================================
    # Run Deconvolution
    # =========================================================================
    print("\n" + "="*70)
    print("RUNNING DECONVOLUTION")
    print("="*70)

    # Run A: Uniform weights
    print("\nRun A: Uniform weights (no leverage)...")
    prop_uniform, info_uniform = run_deconv_with_weights(
        Y, X, coords, gene_idx, leverage_scores,
        use_leverage_weights=False, verbose=True
    )

    # Run B: Leverage weights
    print("\nRun B: Leverage weights (FlashDeconv default)...")
    prop_leverage, info_leverage = run_deconv_with_weights(
        Y, X, coords, gene_idx, leverage_scores,
        use_leverage_weights=True, verbose=True
    )

    # =========================================================================
    # Evaluation
    # =========================================================================
    print("\n" + "="*70)
    print("EVALUATION")
    print("="*70)

    # Find endothelial cell type
    endo_idx = None
    for i, ct in enumerate(cell_types):
        if 'endo' in ct.lower():
            endo_idx = i
            endo_name = ct
            break

    if endo_idx is None:
        print("Endothelial cell type not found!")
        return

    print(f"\nEvaluating: {endo_name}")

    # Get Cldn5 expression
    if 'Cldn5' in visium.var_names:
        cldn5_expr = visium[:, 'Cldn5'].X
        if hasattr(cldn5_expr, 'toarray'):
            cldn5_expr = cldn5_expr.toarray().flatten()
        else:
            cldn5_expr = np.asarray(cldn5_expr).flatten()
        cldn5_expr = cldn5_expr / (cldn5_expr.max() + 1e-10)

        # Correlations
        r_uniform, p_uniform = pearsonr(cldn5_expr, prop_uniform[:, endo_idx])
        r_leverage, p_leverage = pearsonr(cldn5_expr, prop_leverage[:, endo_idx])

        print(f"\nCorrelation with Cldn5 marker:")
        print(f"  Uniform weights:  r = {r_uniform:.4f} (p = {p_uniform:.2e})")
        print(f"  Leverage weights: r = {r_leverage:.4f} (p = {p_leverage:.2e})")
        print(f"  Improvement: {(r_leverage - r_uniform) / abs(r_uniform) * 100:.1f}%")

    # =========================================================================
    # Visualization
    # =========================================================================
    print("\nGenerating visualizations...")

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # Cldn5 expression
    sc = axes[0].scatter(coords[:, 0], coords[:, 1], c=cldn5_expr,
                        cmap='plasma', s=10, alpha=0.8)
    axes[0].set_title(f'Cldn5 Expression\n(Endothelial Marker)', fontsize=12)
    axes[0].axis('off')
    plt.colorbar(sc, ax=axes[0], shrink=0.5)

    # Uniform weights prediction
    sc = axes[1].scatter(coords[:, 0], coords[:, 1], c=prop_uniform[:, endo_idx],
                        cmap='plasma', s=10, alpha=0.8)
    axes[1].set_title(f'Uniform Weights\nr = {r_uniform:.3f}', fontsize=12)
    axes[1].axis('off')
    plt.colorbar(sc, ax=axes[1], shrink=0.5)

    # Leverage weights prediction
    sc = axes[2].scatter(coords[:, 0], coords[:, 1], c=prop_leverage[:, endo_idx],
                        cmap='plasma', s=10, alpha=0.8)
    axes[2].set_title(f'Leverage Weights\nr = {r_leverage:.3f}', fontsize=12)
    axes[2].axis('off')
    plt.colorbar(sc, ax=axes[2], shrink=0.5)

    plt.suptitle(f'Endothelial Cell Prediction: Sketching Weight Comparison\n{endo_name}',
                fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "sketching_comparison_endo.png", dpi=200, bbox_inches='tight')
    plt.savefig(OUTPUT_DIR / "sketching_comparison_endo.pdf", bbox_inches='tight')

    # =========================================================================
    # Test multiple rare cell types
    # =========================================================================
    print("\n" + "="*70)
    print("MULTI-CELL TYPE EVALUATION")
    print("="*70)

    results = []

    # Define markers for different cell types
    ct_markers = {
        'Endo': 'Cldn5',
        'Astro': 'Gfap',
        'Oligo': 'Mbp',
        'Micro': 'Cx3cr1',
    }

    for ct_prefix, marker in ct_markers.items():
        # Find cell type
        ct_idx = None
        ct_full_name = None
        for i, ct in enumerate(cell_types):
            if ct_prefix.lower() in ct.lower():
                ct_idx = i
                ct_full_name = ct
                break

        if ct_idx is None:
            continue

        # Get marker expression
        if marker not in visium.var_names:
            continue

        marker_expr = visium[:, marker].X
        if hasattr(marker_expr, 'toarray'):
            marker_expr = marker_expr.toarray().flatten()
        else:
            marker_expr = np.asarray(marker_expr).flatten()

        if marker_expr.max() == 0:
            continue

        marker_expr = marker_expr / marker_expr.max()

        # Compute correlations
        r_uni, _ = pearsonr(marker_expr, prop_uniform[:, ct_idx])
        r_lev, _ = pearsonr(marker_expr, prop_leverage[:, ct_idx])

        # Count cells (rarity)
        cell_count = 0
        for ct in cell_types:
            if ct_prefix.lower() in ct.lower():
                # This is approximate - just count matching cell types
                cell_count += 1

        results.append({
            'cell_type': ct_full_name,
            'marker': marker,
            'r_uniform': r_uni,
            'r_leverage': r_lev,
            'improvement_pct': (r_lev - r_uni) / abs(r_uni) * 100 if r_uni != 0 else 0,
        })

        print(f"\n{ct_full_name}:")
        print(f"  Uniform:  r = {r_uni:.4f}")
        print(f"  Leverage: r = {r_lev:.4f}")

    # Save results
    results_df = pd.DataFrame(results)
    results_df.to_csv(OUTPUT_DIR / "multi_celltype_results.csv", index=False)

    # Bar chart comparison
    if len(results_df) > 0:
        fig, ax = plt.subplots(figsize=(10, 6))

        x = np.arange(len(results_df))
        width = 0.35

        bars1 = ax.bar(x - width/2, results_df['r_uniform'], width,
                      label='Uniform Weights', color='#3498db', alpha=0.8)
        bars2 = ax.bar(x + width/2, results_df['r_leverage'], width,
                      label='Leverage Weights', color='#e74c3c', alpha=0.8)

        ax.set_xlabel('Cell Type', fontsize=12)
        ax.set_ylabel('Pearson Correlation with Marker', fontsize=12)
        ax.set_title('Sketching Weight Comparison: Impact on Cell Type Prediction\n'
                    'Same genes, different importance weights', fontsize=13)
        ax.set_xticks(x)
        ax.set_xticklabels([ct[:20] + '...' if len(ct) > 20 else ct
                          for ct in results_df['cell_type']], rotation=45, ha='right')
        ax.legend()
        ax.grid(True, alpha=0.3, axis='y')

        plt.tight_layout()
        plt.savefig(OUTPUT_DIR / "multi_celltype_comparison.png", dpi=200, bbox_inches='tight')

    print(f"\nResults saved to: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
