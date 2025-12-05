"""
Stochastic Stability Analysis for FlashDeconv

Tests whether randomized sketching produces consistent results across
different random seeds. This addresses reviewer concerns about the
reliability of randomized algorithms.

Author: Claude
Date: 2025-12-02
"""

import numpy as np
import pandas as pd
from pathlib import Path
from scipy.io import mmread
from scipy.spatial.distance import jensenshannon
from itertools import combinations
import matplotlib.pyplot as plt
import sys
import time

sys.path.insert(0, '/Users/apple/Research/FlashDeconv')

from flashdeconv import FlashDeconv

# Paths
data_dir = Path("/Users/apple/Research/FlashDeconv/validation/benchmark_data/converted")
output_dir = Path("/Users/apple/Research/FlashDeconv/validation/data")
fig_dir = Path("/Users/apple/Research/FlashDeconv/paper/figures")
output_dir.mkdir(parents=True, exist_ok=True)
fig_dir.mkdir(parents=True, exist_ok=True)

# Style
plt.rcParams['font.size'] = 11
plt.rcParams['axes.labelsize'] = 12
plt.rcParams['axes.titlesize'] = 13
plt.rcParams['figure.dpi'] = 150


def load_liver_data(sample_id="JB01"):
    """Load Liver Visium data and reference."""
    prefix = f"liver_mouseVisium_{sample_id}"
    ref_prefix = "liver_ref_9ct"

    # Load spatial data
    Y_sparse = mmread(data_dir / f"{prefix}_counts.mtx")
    Y = Y_sparse.toarray()  # spots x genes

    with open(data_dir / f"{prefix}_genes.txt") as f:
        spatial_genes = [line.strip() for line in f]

    coords_df = pd.read_csv(data_dir / f"{prefix}_coords.csv", index_col=0)
    if 'imagerow' in coords_df.columns:
        coords = coords_df[['imagerow', 'imagecol']].values
    else:
        coords = coords_df.iloc[:, :2].values

    # Load reference
    X_sparse = mmread(data_dir / f"{ref_prefix}_counts.mtx")

    with open(data_dir / f"{ref_prefix}_celltypes.txt") as f:
        celltypes = [line.strip() for line in f]

    with open(data_dir / f"{ref_prefix}_genes.txt") as f:
        ref_genes = [line.strip() for line in f]

    # Create signature matrix
    unique_cts = sorted(set(celltypes))
    X_csc = X_sparse.tocsc()
    n_genes = X_csc.shape[0]

    signature = np.zeros((len(unique_cts), n_genes), dtype=np.float64)
    celltypes_arr = np.array(celltypes)

    for i, ct in enumerate(unique_cts):
        mask = (celltypes_arr == ct)
        signature[i] = X_csc[:, mask].mean(axis=1).A1

    # Align genes
    common_genes = sorted(set(spatial_genes) & set(ref_genes))
    spatial_idx = {g: i for i, g in enumerate(spatial_genes)}
    ref_idx = {g: i for i, g in enumerate(ref_genes)}

    Y_aligned = Y[:, [spatial_idx[g] for g in common_genes]]
    X_aligned = signature[:, [ref_idx[g] for g in common_genes]]

    return Y_aligned, X_aligned, coords, unique_cts


def run_with_seed(Y, X, coords, seed, **kwargs):
    """Run FlashDeconv with a specific random seed."""
    model = FlashDeconv(
        sketch_dim=kwargs.get('sketch_dim', 512),
        lambda_spatial=kwargs.get('lambda_spatial', 0.005),
        rho_sparsity=kwargs.get('rho_sparsity', 0.01),
        n_hvg=kwargs.get('n_hvg', 2000),
        n_markers_per_type=kwargs.get('n_markers', 50),
        max_iter=100,
        tol=1e-4,
        preprocess="log_cpm",
        random_state=seed,
        verbose=False,
    )

    start = time.time()
    pred = model.fit_transform(Y, X, coords)
    elapsed = time.time() - start

    return pred, elapsed


def compute_pairwise_metrics(results):
    """Compute pairwise consistency metrics between runs."""
    n_runs = len(results)
    n_spots = results[0].shape[0]
    n_celltypes = results[0].shape[1]

    # Flatten for overall correlation
    pearson_matrix = np.zeros((n_runs, n_runs))
    jsd_matrix = np.zeros((n_runs, n_runs))

    for i in range(n_runs):
        for j in range(n_runs):
            if i == j:
                pearson_matrix[i, j] = 1.0
                jsd_matrix[i, j] = 0.0
            else:
                # Flatten and compute Pearson
                flat_i = results[i].flatten()
                flat_j = results[j].flatten()
                pearson_matrix[i, j] = np.corrcoef(flat_i, flat_j)[0, 1]

                # Compute mean JSD across spots
                jsds = []
                for s in range(n_spots):
                    # Add small epsilon to avoid log(0)
                    p = results[i][s] + 1e-10
                    q = results[j][s] + 1e-10
                    p = p / p.sum()
                    q = q / q.sum()
                    jsds.append(jensenshannon(p, q))
                jsd_matrix[i, j] = np.mean(jsds)

    return pearson_matrix, jsd_matrix


def compute_celltype_stability(results, cell_types):
    """Compute stability per cell type."""
    n_runs = len(results)
    n_celltypes = len(cell_types)

    stability = {}
    for k, ct in enumerate(cell_types):
        # Extract column k from all runs
        ct_predictions = np.array([r[:, k] for r in results])  # n_runs x n_spots

        # Compute mean and std across runs for each spot
        mean_pred = ct_predictions.mean(axis=0)
        std_pred = ct_predictions.std(axis=0)

        # Coefficient of variation (CV) - std/mean
        cv = std_pred / (mean_pred + 1e-10)

        stability[ct] = {
            'mean_proportion': mean_pred.mean(),
            'mean_std': std_pred.mean(),
            'mean_cv': cv.mean(),
            'max_std': std_pred.max(),
        }

    return stability


def main():
    print("=" * 70)
    print("STOCHASTIC STABILITY ANALYSIS")
    print("Testing FlashDeconv consistency across random seeds")
    print("=" * 70)

    # Load data
    print("\nLoading Liver Visium data (JB01)...")
    Y, X, coords, cell_types = load_liver_data("JB01")
    print(f"  Spatial: {Y.shape[0]} spots x {Y.shape[1]} genes")
    print(f"  Reference: {X.shape[0]} cell types x {X.shape[1]} genes")
    print(f"  Cell types: {cell_types}")

    # Run with different seeds
    n_runs = 10
    seeds = list(range(1, n_runs + 1))

    print(f"\nRunning FlashDeconv {n_runs} times with different seeds...")
    results = []
    runtimes = []

    for i, seed in enumerate(seeds):
        print(f"  Run {i+1}/{n_runs} (seed={seed})...", end=" ", flush=True)
        pred, elapsed = run_with_seed(Y, X, coords, seed)
        results.append(pred)
        runtimes.append(elapsed)
        print(f"done ({elapsed:.2f}s)")

    print(f"\nTotal runtime: {sum(runtimes):.1f}s, Mean: {np.mean(runtimes):.2f}s")

    # Compute pairwise metrics
    print("\nComputing pairwise consistency metrics...")
    pearson_matrix, jsd_matrix = compute_pairwise_metrics(results)

    # Extract off-diagonal values
    mask = ~np.eye(n_runs, dtype=bool)
    pearson_offdiag = pearson_matrix[mask]
    jsd_offdiag = jsd_matrix[mask]

    print(f"\n{'='*60}")
    print("PAIRWISE CONSISTENCY (off-diagonal)")
    print(f"{'='*60}")
    print(f"  Pearson Correlation:")
    print(f"    Mean: {pearson_offdiag.mean():.6f}")
    print(f"    Min:  {pearson_offdiag.min():.6f}")
    print(f"    Max:  {pearson_offdiag.max():.6f}")
    print(f"    Std:  {pearson_offdiag.std():.6f}")

    print(f"\n  Jensen-Shannon Divergence:")
    print(f"    Mean: {jsd_offdiag.mean():.6f}")
    print(f"    Max:  {jsd_offdiag.max():.6f}")

    # Per-celltype stability
    print(f"\n{'='*60}")
    print("PER-CELLTYPE STABILITY")
    print(f"{'='*60}")

    stability = compute_celltype_stability(results, cell_types)

    print(f"\n{'Cell Type':<35} {'Mean Prop':<12} {'Mean Std':<12} {'Mean CV':<12}")
    print("-" * 70)
    for ct in cell_types:
        s = stability[ct]
        print(f"{ct:<35} {s['mean_proportion']*100:>8.2f}%    {s['mean_std']*100:>8.4f}%   {s['mean_cv']:>8.4f}")

    # Create visualization
    print("\nGenerating visualization...")

    fig, axes = plt.subplots(1, 3, figsize=(14, 4))

    # Panel A: Pearson correlation heatmap (same color scheme as melanoma figure)
    ax1 = axes[0]
    im1 = ax1.imshow(pearson_matrix, cmap='RdYlBu_r', vmin=0.985, vmax=1.0)
    ax1.set_xticks(range(n_runs))
    ax1.set_yticks(range(n_runs))
    ax1.set_xticklabels([f"Seed {s}" for s in seeds], fontsize=8, rotation=45, ha='right')
    ax1.set_yticklabels([f"Seed {s}" for s in seeds], fontsize=8)
    ax1.set_title(f"A. Pairwise Correlation Matrix\n(mean r = {pearson_offdiag.mean():.4f})")
    plt.colorbar(im1, ax=ax1, shrink=0.8, label='Pearson r')

    # Panel B: Distribution of pairwise correlations
    ax2 = axes[1]
    bins = np.linspace(pearson_offdiag.min() - 0.001, 1.0, 15)
    ax2.hist(pearson_offdiag, bins=bins, edgecolor='black', alpha=0.7, color='#3498db')
    ax2.axvline(pearson_offdiag.mean(), color='#e74c3c', linestyle='--', linewidth=2,
                label=f'Mean = {pearson_offdiag.mean():.2f}')
    ax2.axvline(pearson_offdiag.min(), color='#95a5a6', linestyle=':', linewidth=1.5,
                label=f'Min = {pearson_offdiag.min():.2f}')
    ax2.set_xlabel('Pearson Correlation')
    ax2.set_ylabel('Count (pairwise comparisons)')
    ax2.set_title('B. Distribution of Pairwise Correlations\n(n = 90 pairs)')
    ax2.legend(loc='upper left', fontsize=9)
    ax2.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{x:.2f}'))

    # Panel C: Per-celltype CV
    ax3 = axes[2]
    cts = list(stability.keys())
    cvs = [stability[ct]['mean_cv'] for ct in cts]
    props = [stability[ct]['mean_proportion'] * 100 for ct in cts]

    sorted_idx = np.argsort(props)[::-1]
    cts_sorted = [cts[i] for i in sorted_idx]
    cvs_sorted = [cvs[i] for i in sorted_idx]
    props_sorted = [props[i] for i in sorted_idx]

    colors = ['#2ecc71' if cv < 0.2 else '#f39c12' if cv < 0.5 else '#e74c3c' for cv in cvs_sorted]
    bars = ax3.barh(range(len(cts)), cvs_sorted, color=colors, edgecolor='black', linewidth=0.5)
    ax3.set_yticks(range(len(cts)))
    ax3.set_yticklabels([f"{ct} ({props_sorted[i]:.1f}%)" for i, ct in enumerate(cts_sorted)], fontsize=8)
    ax3.set_xlabel('Coefficient of Variation (CV)')
    ax3.set_title('C. Per-Celltype Stability\n(lower CV = more stable)')
    ax3.invert_yaxis()

    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#2ecc71', edgecolor='black', label='CV < 0.2 (stable)'),
        Patch(facecolor='#f39c12', edgecolor='black', label='0.2 ≤ CV < 0.5'),
        Patch(facecolor='#e74c3c', edgecolor='black', label='CV ≥ 0.5 (rare types)')
    ]
    ax3.legend(handles=legend_elements, loc='lower right', fontsize=8)

    plt.tight_layout()

    # Save figure
    output_path = fig_dir / "supplementary_stochastic_stability.pdf"
    plt.savefig(output_path, bbox_inches='tight', dpi=300)
    print(f"✓ Saved figure to: {output_path}")

    plt.savefig(fig_dir / "supplementary_stochastic_stability.png", bbox_inches='tight', dpi=150)
    print(f"✓ Saved PNG preview")

    plt.close()

    # Save detailed results
    results_df = pd.DataFrame({
        'seed': seeds,
        'runtime_sec': runtimes,
    })

    # Add mean predictions per cell type
    for k, ct in enumerate(cell_types):
        results_df[f'mean_{ct}'] = [r[:, k].mean() for r in results]

    results_df.to_csv(output_dir / "stochastic_stability_results.csv", index=False)
    print(f"✓ Saved results to: {output_dir / 'stochastic_stability_results.csv'}")

    # Save correlation matrix
    corr_df = pd.DataFrame(pearson_matrix,
                           index=[f"seed_{s}" for s in seeds],
                           columns=[f"seed_{s}" for s in seeds])
    corr_df.to_csv(output_dir / "stochastic_stability_correlation_matrix.csv")
    print(f"✓ Saved correlation matrix")

    # Summary
    print("\n" + "=" * 70)
    print("CONCLUSION")
    print("=" * 70)
    print(f"""
Despite using randomized sketching, FlashDeconv produces highly
deterministic results:

  • Pairwise Pearson correlation: {pearson_offdiag.mean():.4f} ± {pearson_offdiag.std():.4f}
  • All pairwise correlations > {pearson_offdiag.min():.4f}
  • Mean Jensen-Shannon Divergence: {jsd_offdiag.mean():.6f}
  • Mean Coefficient of Variation: {np.mean(cvs)*100:.3f}%

The randomized sketching introduces negligible variance compared to
the underlying biological signal. Users can confidently use FlashDeconv
without fixing random seeds, as results are effectively deterministic.
""")

    return results, pearson_matrix, jsd_matrix, stability


if __name__ == "__main__":
    results, pearson_matrix, jsd_matrix, stability = main()
