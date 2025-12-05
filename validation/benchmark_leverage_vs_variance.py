"""
Benchmark: Leverage-based vs Variance-based Gene Selection for Deconvolution
============================================================================

This script demonstrates that Leverage-based gene selection improves
deconvolution accuracy, especially for rare cell types.

Experiments:
1. Run deconvolution with HVG (Variance-based) gene selection
2. Run deconvolution with Leverage-based gene selection
3. Compare predictions on rare cell types (Endothelial, Pericytes)

Key metrics:
- Correlation with marker gene expression (Cldn5 for Endothelial)
- Spatial structure of predictions
- Co-localization of related cell types
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

from flashdeconv import FlashDeconv
from flashdeconv.utils.genes import select_hvg, compute_leverage_scores, select_markers
from flashdeconv.core.sketching import sketch_data
from flashdeconv.core.solver import bcd_solve, normalize_proportions
from flashdeconv.utils.graph import coords_to_adjacency

# Output directory
OUTPUT_DIR = Path("validation/results/benchmark_leverage_vs_variance")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# =============================================================================
# Data Loading
# =============================================================================
def load_data():
    """Load scRNA reference and Visium spatial data."""
    print("Loading data...")

    # Load scRNA reference
    scrna = sc.read_h5ad("validation/mouse_brain/scrna_reference.h5ad")
    annotations = pd.read_csv("validation/mouse_brain/cell_annotation.csv", index_col=0)

    # Align cells
    common_cells = scrna.obs_names.intersection(annotations.index)
    scrna = scrna[common_cells].copy()
    scrna.obs['cell_type'] = annotations.loc[common_cells, 'annotation_1'].values

    # Load Visium
    visium = sc.read_h5ad("validation/mouse_brain/visium_mouse_brain.h5ad")
    visium.var_names_make_unique()

    # Add gene symbol mapping to scrna
    if 'SYMBOL' in scrna.var.columns:
        scrna.var['gene_symbol'] = scrna.var['SYMBOL'].fillna(scrna.var_names.to_series())
    else:
        scrna.var['gene_symbol'] = scrna.var_names

    print(f"  scRNA: {scrna.shape[0]} cells, {scrna.shape[1]} genes")
    print(f"  Visium: {visium.shape[0]} spots, {visium.shape[1]} genes")
    print(f"  Cell types: {scrna.obs['cell_type'].nunique()}")

    return scrna, visium


def build_reference_matrix(adata, cell_type_col='cell_type'):
    """Build cell type signature matrix (K x G)."""
    cell_types = sorted(adata.obs[cell_type_col].unique())
    n_genes = adata.n_vars

    X = np.zeros((len(cell_types), n_genes))

    for i, ct in enumerate(cell_types):
        mask = adata.obs[cell_type_col] == ct
        if hasattr(adata.X, 'toarray'):
            X[i, :] = np.asarray(adata[mask].X.toarray().mean(axis=0)).flatten()
        else:
            X[i, :] = np.asarray(adata[mask].X.mean(axis=0)).flatten()

    gene_names = adata.var_names.tolist()
    gene_symbols = adata.var['gene_symbol'].tolist()

    return X, cell_types, gene_names, gene_symbols


def find_common_genes(scrna, visium):
    """Find genes common to both datasets."""
    # Get gene symbols from scrna
    if 'SYMBOL' in scrna.var.columns:
        scrna_symbols = scrna.var['SYMBOL'].fillna(scrna.var_names.to_series()).tolist()
    else:
        scrna_symbols = scrna.var_names.tolist()

    # Visium uses symbols as var_names
    visium_symbols = set(visium.var_names)

    # Find common genes (by symbol)
    common_symbols = [s for s in scrna_symbols if s in visium_symbols]

    # Get indices
    scrna_idx = [i for i, s in enumerate(scrna_symbols) if s in visium_symbols]

    # Map symbol to visium index
    visium_symbol_to_idx = {s: i for i, s in enumerate(visium.var_names)}
    visium_idx = [visium_symbol_to_idx[s] for s in common_symbols]

    print(f"  Common genes: {len(common_symbols)}")

    return scrna_idx, visium_idx, common_symbols


# =============================================================================
# Gene Selection Methods
# =============================================================================
def select_genes_variance(Y, n_genes=2000):
    """Select top genes by variance (HVG method)."""
    hvg_idx = select_hvg(Y, n_top=n_genes)
    return hvg_idx


def select_genes_leverage(X, n_genes=2000):
    """Select top genes by leverage score."""
    leverage = compute_leverage_scores(X)
    top_idx = np.argsort(-leverage)[:n_genes]
    return np.sort(top_idx), leverage


def select_genes_hybrid(Y, X, n_genes=2000, leverage_weight=0.5):
    """
    Hybrid selection: combine variance and leverage scores.
    This is closer to what FlashDeconv actually does.
    """
    # Compute variance ranks
    if hasattr(Y, 'toarray'):
        Y_dense = Y.toarray()
    else:
        Y_dense = np.asarray(Y)

    gene_var = np.var(Y_dense, axis=0)
    var_rank = np.argsort(np.argsort(-gene_var))

    # Compute leverage ranks
    leverage = compute_leverage_scores(X)
    lev_rank = np.argsort(np.argsort(-leverage))

    # Combined score (lower rank = better)
    combined_rank = (1 - leverage_weight) * var_rank + leverage_weight * lev_rank

    top_idx = np.argsort(combined_rank)[:n_genes]
    return np.sort(top_idx), leverage[top_idx]


# =============================================================================
# Custom Deconvolution Runner
# =============================================================================
def run_deconvolution(Y, X, coords, gene_idx, leverage_scores=None,
                      sketch_dim=512, lambda_spatial=5000, verbose=False):
    """
    Run deconvolution with specific gene selection.

    Parameters
    ----------
    Y : array (n_spots, n_genes)
        Spatial data
    X : array (n_cell_types, n_genes)
        Reference signatures
    coords : array (n_spots, 2)
        Spatial coordinates
    gene_idx : array
        Indices of genes to use
    leverage_scores : array, optional
        Pre-computed leverage scores for selected genes
    """
    # Subset to selected genes
    if sparse.issparse(Y):
        Y_sub = Y[:, gene_idx].toarray()
    else:
        Y_sub = Y[:, gene_idx]
    X_sub = X[:, gene_idx]

    # Preprocessing (log-CPM)
    Y_cpm = Y_sub / (Y_sub.sum(axis=1, keepdims=True) + 1e-10) * 1e4
    X_cpm = X_sub / (X_sub.sum(axis=1, keepdims=True) + 1e-10) * 1e4
    Y_tilde = np.log1p(Y_cpm)
    X_tilde = np.log1p(X_cpm)

    # Compute leverage scores if not provided
    if leverage_scores is None:
        leverage_scores = compute_leverage_scores(X_sub)

    # Sketching
    Y_sketch, X_sketch, Omega = sketch_data(
        Y_tilde, X_tilde,
        sketch_dim=sketch_dim,
        leverage_scores=leverage_scores,
        random_state=42
    )

    # Build spatial graph
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


# =============================================================================
# Evaluation Metrics
# =============================================================================
def evaluate_marker_correlation(proportions, visium, cell_type_idx, marker_gene,
                                 cell_type_name):
    """
    Evaluate correlation between predicted cell type proportion
    and marker gene expression.
    """
    # Get marker expression
    if marker_gene in visium.var_names:
        if hasattr(visium.X, 'toarray'):
            marker_expr = visium[:, marker_gene].X.toarray().flatten()
        else:
            marker_expr = visium[:, marker_gene].X.flatten()

        # Normalize marker expression
        marker_expr = marker_expr / (marker_expr.max() + 1e-10)

        # Get predicted proportions
        pred_prop = proportions[:, cell_type_idx]

        # Compute correlations
        pearson_r, pearson_p = pearsonr(marker_expr, pred_prop)
        spearman_r, spearman_p = spearmanr(marker_expr, pred_prop)

        return {
            'cell_type': cell_type_name,
            'marker': marker_gene,
            'pearson_r': pearson_r,
            'pearson_p': pearson_p,
            'spearman_r': spearman_r,
            'spearman_p': spearman_p,
        }
    else:
        print(f"  Warning: Marker {marker_gene} not found in Visium data")
        return None


def compute_spatial_smoothness(proportions, coords, k=6):
    """
    Compute spatial smoothness of predictions.
    Smoother predictions (higher autocorrelation) are generally better.
    """
    from scipy.spatial.distance import cdist

    dist_matrix = cdist(coords, coords)
    n_spots = len(coords)

    smoothness_scores = []

    for ct_idx in range(proportions.shape[1]):
        prop = proportions[:, ct_idx]

        # For each spot, compute correlation with neighbors
        local_corrs = []
        for i in range(n_spots):
            neighbors = np.argsort(dist_matrix[i])[1:k+1]
            if np.std(prop[neighbors]) > 1e-10 and np.std([prop[i]]) > 1e-10:
                corr = np.corrcoef(prop[i], prop[neighbors].mean())[0, 1]
                if not np.isnan(corr):
                    local_corrs.append(corr)

        if len(local_corrs) > 0:
            smoothness_scores.append(np.mean(local_corrs))
        else:
            smoothness_scores.append(0)

    return np.array(smoothness_scores)


# =============================================================================
# Visualization
# =============================================================================
def plot_spatial_comparison(visium, proportions_var, proportions_lev,
                            cell_type_idx, cell_type_name, output_path):
    """Create side-by-side spatial comparison plot."""

    coords = visium.obsm['spatial']

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # Get marker gene expression (Cldn5 for endothelial)
    marker = 'Cldn5' if 'Cldn5' in visium.var_names else None

    if marker:
        if hasattr(visium.X, 'toarray'):
            marker_expr = visium[:, marker].X.toarray().flatten()
        else:
            marker_expr = visium[:, marker].X.flatten()

        # Plot 1: Marker expression (ground truth proxy)
        sc = axes[0].scatter(coords[:, 0], coords[:, 1], c=marker_expr,
                            cmap='plasma', s=10, alpha=0.8)
        axes[0].set_title(f'Marker: {marker}\n(Ground Truth Proxy)', fontsize=12)
        axes[0].axis('off')
        plt.colorbar(sc, ax=axes[0], shrink=0.5)
    else:
        axes[0].text(0.5, 0.5, 'Marker not found', ha='center', va='center')
        axes[0].axis('off')

    # Plot 2: Variance-based prediction
    pred_var = proportions_var[:, cell_type_idx]
    sc = axes[1].scatter(coords[:, 0], coords[:, 1], c=pred_var,
                        cmap='plasma', s=10, alpha=0.8)
    axes[1].set_title(f'HVG (Variance)\n{cell_type_name}', fontsize=12)
    axes[1].axis('off')
    plt.colorbar(sc, ax=axes[1], shrink=0.5)

    # Plot 3: Leverage-based prediction
    pred_lev = proportions_lev[:, cell_type_idx]
    sc = axes[2].scatter(coords[:, 0], coords[:, 1], c=pred_lev,
                        cmap='plasma', s=10, alpha=0.8)
    axes[2].set_title(f'Leverage-based\n{cell_type_name}', fontsize=12)
    axes[2].axis('off')
    plt.colorbar(sc, ax=axes[2], shrink=0.5)

    plt.suptitle(f'Spatial Prediction Comparison: {cell_type_name}',
                fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(output_path, dpi=200, bbox_inches='tight')
    plt.close()


def plot_correlation_comparison(results_df, output_path):
    """Plot bar chart comparing correlations."""

    fig, ax = plt.subplots(figsize=(10, 6))

    cell_types = results_df['cell_type'].unique()
    x = np.arange(len(cell_types))
    width = 0.35

    var_corrs = results_df[results_df['method'] == 'HVG']['pearson_r'].values
    lev_corrs = results_df[results_df['method'] == 'Leverage']['pearson_r'].values

    bars1 = ax.bar(x - width/2, var_corrs, width, label='HVG (Variance)',
                   color='#3498db', alpha=0.8)
    bars2 = ax.bar(x + width/2, lev_corrs, width, label='Leverage-based',
                   color='#e74c3c', alpha=0.8)

    ax.set_xlabel('Cell Type', fontsize=12)
    ax.set_ylabel('Pearson Correlation with Marker', fontsize=12)
    ax.set_title('Deconvolution Accuracy: Leverage vs Variance Gene Selection\n'
                 'Correlation between Predicted Proportion and Marker Expression',
                fontsize=13, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(cell_types, rotation=45, ha='right')
    ax.legend()
    ax.grid(True, alpha=0.3, axis='y')

    # Add value labels
    for bar in bars1:
        height = bar.get_height()
        ax.annotate(f'{height:.3f}',
                   xy=(bar.get_x() + bar.get_width()/2, height),
                   xytext=(0, 3), textcoords="offset points",
                   ha='center', va='bottom', fontsize=9)
    for bar in bars2:
        height = bar.get_height()
        ax.annotate(f'{height:.3f}',
                   xy=(bar.get_x() + bar.get_width()/2, height),
                   xytext=(0, 3), textcoords="offset points",
                   ha='center', va='bottom', fontsize=9)

    plt.tight_layout()
    plt.savefig(output_path, dpi=200, bbox_inches='tight')
    plt.close()


# =============================================================================
# Main Benchmark
# =============================================================================
def main():
    print("="*70)
    print("BENCHMARK: Leverage vs Variance Gene Selection for Deconvolution")
    print("="*70)

    # Load data
    scrna, visium = load_data()

    # Find common genes
    scrna_idx, visium_idx, common_symbols = find_common_genes(scrna, visium)

    # Subset to common genes
    scrna_common = scrna[:, scrna_idx].copy()
    visium_common = visium[:, visium_idx].copy()

    # Build reference matrix
    print("\nBuilding reference matrix...")
    X, cell_types, _, _ = build_reference_matrix(scrna_common)
    print(f"  Reference: {X.shape[0]} cell types x {X.shape[1]} genes")

    # Get spatial data matrix and coordinates
    if hasattr(visium_common.X, 'toarray'):
        Y = visium_common.X.toarray()
    else:
        Y = np.asarray(visium_common.X)
    coords = visium_common.obsm['spatial']

    print(f"  Spatial: {Y.shape[0]} spots x {Y.shape[1]} genes")

    # =========================================================================
    # Gene Selection
    # =========================================================================
    print("\n" + "="*70)
    print("GENE SELECTION")
    print("="*70)

    n_genes = 2000

    # Method A: Variance-based (HVG)
    print(f"\nMethod A: Top {n_genes} HVG (Variance-based)...")
    hvg_idx = select_genes_variance(Y, n_genes=n_genes)
    print(f"  Selected {len(hvg_idx)} genes")

    # Check if key markers are included
    cldn5_idx = common_symbols.index('Cldn5') if 'Cldn5' in common_symbols else -1
    if cldn5_idx >= 0:
        print(f"  Cldn5 (endothelial marker) in HVG: {cldn5_idx in hvg_idx}")

    # Method B: Leverage-based
    print(f"\nMethod B: Top {n_genes} Leverage genes...")
    lev_idx, leverage_scores_all = select_genes_leverage(X, n_genes=n_genes)
    print(f"  Selected {len(lev_idx)} genes")

    if cldn5_idx >= 0:
        print(f"  Cldn5 (endothelial marker) in Leverage: {cldn5_idx in lev_idx}")

    # Save gene lists
    hvg_genes = [common_symbols[i] for i in hvg_idx]
    lev_genes = [common_symbols[i] for i in lev_idx]

    pd.DataFrame({'gene': hvg_genes}).to_csv(OUTPUT_DIR / "hvg_genes.csv", index=False)
    pd.DataFrame({'gene': lev_genes}).to_csv(OUTPUT_DIR / "leverage_genes.csv", index=False)

    # =========================================================================
    # Run Deconvolution
    # =========================================================================
    print("\n" + "="*70)
    print("RUNNING DECONVOLUTION")
    print("="*70)

    # Run A: HVG (Variance) - use uniform leverage for fair comparison
    print("\nRun A: HVG (Variance-based) selection...")
    uniform_leverage = np.ones(len(hvg_idx)) / len(hvg_idx)
    proportions_var, info_var = run_deconvolution(
        Y, X, coords, hvg_idx,
        leverage_scores=uniform_leverage,  # Uniform weights
        verbose=True
    )
    print(f"  Converged: {info_var['converged']}, Iterations: {info_var['n_iterations']}")

    # Run B: Leverage-based
    print("\nRun B: Leverage-based selection...")
    leverage_scores_selected = compute_leverage_scores(X[:, lev_idx])
    proportions_lev, info_lev = run_deconvolution(
        Y, X, coords, lev_idx,
        leverage_scores=leverage_scores_selected,  # Use leverage weights
        verbose=True
    )
    print(f"  Converged: {info_lev['converged']}, Iterations: {info_lev['n_iterations']}")

    # =========================================================================
    # Evaluation
    # =========================================================================
    print("\n" + "="*70)
    print("EVALUATION: Rare Cell Types")
    print("="*70)

    # Define rare cell types and their markers
    rare_cells = {
        'Endo': 'Cldn5',      # Endothelial cells
        'Peri': 'Rgs5',       # Pericytes
        'VLMC': 'Slc6a13',    # Vascular leptomeningeal cells
    }

    results = []

    for ct_name, marker in rare_cells.items():
        # Find cell type index
        ct_matches = [i for i, ct in enumerate(cell_types) if ct_name.lower() in ct.lower()]

        if not ct_matches:
            print(f"\n  Warning: Cell type '{ct_name}' not found")
            continue

        ct_idx = ct_matches[0]
        full_ct_name = cell_types[ct_idx]

        print(f"\n  Evaluating: {full_ct_name} (marker: {marker})")

        # Evaluate HVG method
        eval_var = evaluate_marker_correlation(
            proportions_var, visium_common, ct_idx, marker, full_ct_name
        )
        if eval_var:
            eval_var['method'] = 'HVG'
            results.append(eval_var)
            print(f"    HVG:      Pearson r = {eval_var['pearson_r']:.4f}")

        # Evaluate Leverage method
        eval_lev = evaluate_marker_correlation(
            proportions_lev, visium_common, ct_idx, marker, full_ct_name
        )
        if eval_lev:
            eval_lev['method'] = 'Leverage'
            results.append(eval_lev)
            print(f"    Leverage: Pearson r = {eval_lev['pearson_r']:.4f}")

        # Plot spatial comparison
        if ct_matches:
            plot_spatial_comparison(
                visium_common, proportions_var, proportions_lev,
                ct_idx, full_ct_name,
                OUTPUT_DIR / f"spatial_comparison_{ct_name}.png"
            )

    # Save results
    results_df = pd.DataFrame(results)
    results_df.to_csv(OUTPUT_DIR / "correlation_results.csv", index=False)

    # Plot correlation comparison
    if len(results_df) > 0:
        plot_correlation_comparison(results_df, OUTPUT_DIR / "correlation_comparison.png")

    # =========================================================================
    # Summary Statistics
    # =========================================================================
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)

    if len(results_df) > 0:
        hvg_mean = results_df[results_df['method'] == 'HVG']['pearson_r'].mean()
        lev_mean = results_df[results_df['method'] == 'Leverage']['pearson_r'].mean()

        print(f"\nMean Pearson correlation with markers:")
        print(f"  HVG (Variance):  {hvg_mean:.4f}")
        print(f"  Leverage-based:  {lev_mean:.4f}")
        print(f"  Improvement:     {(lev_mean - hvg_mean) / hvg_mean * 100:.1f}%")

    # Save proportions
    np.save(OUTPUT_DIR / "proportions_hvg.npy", proportions_var)
    np.save(OUTPUT_DIR / "proportions_leverage.npy", proportions_lev)

    # Save cell type names
    pd.DataFrame({'cell_type': cell_types}).to_csv(
        OUTPUT_DIR / "cell_types.csv", index=False
    )

    print(f"\nResults saved to: {OUTPUT_DIR}")
    print("\nKey outputs:")
    print("  - correlation_results.csv: Marker correlation metrics")
    print("  - correlation_comparison.png: Bar chart comparison")
    print("  - spatial_comparison_*.png: Spatial prediction plots")


if __name__ == "__main__":
    main()
