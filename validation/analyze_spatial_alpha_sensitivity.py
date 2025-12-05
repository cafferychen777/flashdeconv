"""
Spatial Regularization Parameter (α) Sensitivity Analysis

Investigates how the spatial regularization strength affects:
1. Deconvolution accuracy (RMSE)
2. Spatial coherence (boundary preservation)

Uses Brain Cortex Silver Standard with SYNTHETIC coordinates based on Region labels.
This creates a controlled experiment where we KNOW the ground truth boundaries.

Author: Claude
Date: 2025-12-02
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.io import mmread
import sys

sys.path.insert(0, '/Users/apple/Research/FlashDeconv')

from flashdeconv import FlashDeconv
from flashdeconv.utils.metrics import compute_rmse

# Paths
data_dir = Path("/Users/apple/Research/FlashDeconv/validation/benchmark_data")
converted_dir = data_dir / "converted"
output_dir = Path("/Users/apple/Research/FlashDeconv/validation/data")
fig_dir = Path("/Users/apple/Research/FlashDeconv/paper/figures")
output_dir.mkdir(parents=True, exist_ok=True)
fig_dir.mkdir(parents=True, exist_ok=True)

# Style
plt.rcParams['font.size'] = 11
plt.rcParams['axes.labelsize'] = 12
plt.rcParams['axes.titlesize'] = 13
plt.rcParams['figure.dpi'] = 150


def assign_synthetic_coordinates(regions, seed=42):
    """
    Given region labels, assign (x, y) coordinates to create spatial structure.

    Structure: Horizontal stripes (Layers) - mimicking cortical layers.
    Each region occupies a horizontal band.

    Parameters
    ----------
    regions : array-like
        Region labels for each spot (e.g., 'priorregion1', 'priorregion2', ...)
    seed : int
        Random seed for reproducibility

    Returns
    -------
    coords : ndarray
        (n_spots, 2) array of (x, y) coordinates
    """
    np.random.seed(seed)

    unique_regions = sorted(set(regions))
    n_regions = len(unique_regions)
    n_spots = len(regions)

    # Tissue is 100x100 square
    # Each region occupies a horizontal stripe
    height_per_region = 100.0 / n_regions

    coords = np.zeros((n_spots, 2))

    for i, region in enumerate(unique_regions):
        # Find spots belonging to this region
        mask = np.array([r == region for r in regions])
        n_spots_in_region = np.sum(mask)

        # Assign coordinates within the stripe
        # y: [i*h, (i+1)*h]
        # x: [0, 100]
        y_min = i * height_per_region
        y_max = (i + 1) * height_per_region

        # Random scatter within stripe (+1/-1 buffer to avoid exact boundaries)
        x = np.random.uniform(2, 98, n_spots_in_region)
        y = np.random.uniform(y_min + 2, y_max - 2, n_spots_in_region)

        coords[mask, 0] = x
        coords[mask, 1] = y

    return coords


def load_silver_standard_from_rds(rds_path):
    """Load Silver Standard data from RDS file using rpy2."""
    try:
        import rpy2.robjects as ro
        from rpy2.robjects import pandas2ri
        pandas2ri.activate()

        ro.r(f'''
        rds <- readRDS("{rds_path}")
        counts <- as.matrix(rds$counts)
        spot_comp <- rds$spot_composition
        rel_comp <- rds$relative_spot_composition
        ''')

        counts = np.array(ro.r('counts'))
        spot_comp = ro.r('spot_comp')
        rel_comp = ro.r('rel_comp')

        return counts, spot_comp, rel_comp
    except ImportError:
        return None, None, None


def load_silver_standard_converted(dataset_id="silver_1_1"):
    """
    Load Silver Standard from converted directory.

    Returns
    -------
    Y : ndarray (n_spots, n_genes)
    X : ndarray (n_celltypes, n_genes)
    ground_truth : ndarray (n_spots, n_celltypes)
    cell_types : list
    regions : list (for coordinate generation)
    """
    # Load spatial counts
    Y_sparse = mmread(converted_dir / f"{dataset_id}_counts.mtx")
    Y = Y_sparse.toarray().T  # Convert to spots x genes

    with open(converted_dir / f"{dataset_id}_genes.txt") as f:
        spatial_genes = [line.strip() for line in f]

    # Load ground truth proportions
    props_df = pd.read_csv(converted_dir / f"{dataset_id}_proportions.csv", index_col=0)
    cell_types = list(props_df.columns)
    ground_truth = props_df.values

    # Load reference
    X_sparse = mmread(converted_dir / "reference_brain_cortex_counts.mtx")

    with open(converted_dir / "reference_brain_cortex_celltypes.txt") as f:
        ref_celltypes = [line.strip() for line in f]

    with open(converted_dir / "reference_brain_cortex_genes.txt") as f:
        ref_genes = [line.strip() for line in f]

    # Build signature matrix (mean expression per cell type)
    unique_cts = sorted(set(ref_celltypes))
    X_csc = X_sparse.tocsc()
    n_genes = X_csc.shape[0]

    signature = np.zeros((len(unique_cts), n_genes), dtype=np.float64)
    ref_celltypes_arr = np.array(ref_celltypes)

    for i, ct in enumerate(unique_cts):
        mask = (ref_celltypes_arr == ct)
        signature[i] = X_csc[:, mask].mean(axis=1).A1

    # Align genes between spatial and reference
    common_genes = sorted(set(spatial_genes) & set(ref_genes))
    spatial_idx = {g: i for i, g in enumerate(spatial_genes)}
    ref_idx = {g: i for i, g in enumerate(ref_genes)}

    Y_aligned = Y[:, [spatial_idx[g] for g in common_genes]]
    X_aligned = signature[:, [ref_idx[g] for g in common_genes]]

    # Align cell types between ground truth and reference
    ct_mapping = {ct: i for i, ct in enumerate(unique_cts)}

    # Reorder ground truth to match reference cell type order
    gt_reordered = np.zeros((ground_truth.shape[0], len(unique_cts)))
    for j, ct in enumerate(cell_types):
        if ct in ct_mapping:
            gt_reordered[:, ct_mapping[ct]] = ground_truth[:, j]

    return Y_aligned, X_aligned, gt_reordered, unique_cts


def load_silver_standard_with_regions(rds_path):
    """
    Load Silver Standard from RDS and extract region information.
    """
    import subprocess
    import tempfile
    import os

    # Use R to extract data
    with tempfile.TemporaryDirectory() as tmpdir:
        r_script = f'''
        rds <- readRDS("{rds_path}")

        # Save counts
        counts <- as.matrix(rds$counts)
        write.csv(counts, "{tmpdir}/counts.csv")

        # Save spot composition (includes region info)
        write.csv(rds$spot_composition, "{tmpdir}/spot_composition.csv")

        # Save relative composition (ground truth proportions)
        write.csv(rds$relative_spot_composition, "{tmpdir}/relative_composition.csv")
        '''

        result = subprocess.run(['Rscript', '-e', r_script], capture_output=True, text=True)

        if result.returncode != 0:
            print(f"R error: {result.stderr}")
            return None, None, None, None, None

        # Load the CSVs
        counts_df = pd.read_csv(f"{tmpdir}/counts.csv", index_col=0)
        spot_comp_df = pd.read_csv(f"{tmpdir}/spot_composition.csv", index_col=0)
        rel_comp_df = pd.read_csv(f"{tmpdir}/relative_composition.csv", index_col=0)

        Y = counts_df.values.T  # genes x spots -> spots x genes
        genes = list(counts_df.index)

        # Extract regions
        regions = spot_comp_df['region'].tolist()

        # Ground truth proportions
        # Remove 'name' and 'region' columns
        ct_cols = [c for c in rel_comp_df.columns if c not in ['name', 'region']]
        ground_truth = rel_comp_df[ct_cols].values
        cell_types = ct_cols

        return Y, genes, ground_truth, cell_types, regions


def run_with_alpha(Y, X, coords, alpha, **kwargs):
    """Run FlashDeconv with specific alpha value."""
    model = FlashDeconv(
        sketch_dim=kwargs.get('sketch_dim', 512),
        lambda_spatial=alpha,  # This is the alpha parameter
        rho_sparsity=kwargs.get('rho_sparsity', 0.01),
        n_hvg=kwargs.get('n_hvg', 2000),
        n_markers_per_type=kwargs.get('n_markers', 50),
        max_iter=100,
        tol=1e-4,
        preprocess="log_cpm",
        random_state=42,
        verbose=False,
    )

    pred = model.fit_transform(Y, X, coords)
    return pred


def compute_spatial_autocorrelation(values, coords, k=6):
    """Compute Moran's I for spatial autocorrelation."""
    from scipy.spatial import cKDTree

    n = len(values)
    if n < k + 1:
        return np.nan

    tree = cKDTree(coords)

    # Build k-NN weight matrix
    _, indices = tree.query(coords, k=k+1)
    indices = indices[:, 1:]  # Exclude self

    # Compute Moran's I
    mean_val = np.mean(values)
    centered = values - mean_val

    numerator = 0
    denominator = np.sum(centered ** 2)

    if denominator == 0:
        return np.nan

    w_sum = 0

    for i in range(n):
        for j in indices[i]:
            numerator += centered[i] * centered[j]
            w_sum += 1

    morans_i = (n / w_sum) * (numerator / denominator)
    return morans_i


def main():
    print("=" * 70)
    print("SPATIAL REGULARIZATION (α) SENSITIVITY ANALYSIS")
    print("=" * 70)
    print("\nUsing Silver Standard with SYNTHETIC coordinates based on Regions")
    print("This creates a controlled experiment with known boundary ground truth.\n")

    # Load Silver Standard from RDS
    rds_path = data_dir / "silver_standard_1-1" / "brain_cortex_artificial_uniform_distinct_rep1.rds"

    print(f"Loading: {rds_path.name}")
    Y, genes, ground_truth, cell_types, regions = load_silver_standard_with_regions(rds_path)

    if Y is None:
        print("Failed to load RDS. Trying converted format...")
        # Fallback: use converted data (but won't have region info)
        Y, X, ground_truth, cell_types = load_silver_standard_converted("silver_1_1")
        # Generate fake regions based on spot index
        n_spots = Y.shape[0]
        regions = [f"region_{i // (n_spots // 5)}" for i in range(n_spots)]
    else:
        # Load reference and align
        X_sparse = mmread(converted_dir / "reference_brain_cortex_counts.mtx")

        with open(converted_dir / "reference_brain_cortex_celltypes.txt") as f:
            ref_celltypes = [line.strip() for line in f]

        with open(converted_dir / "reference_brain_cortex_genes.txt") as f:
            ref_genes = [line.strip() for line in f]

        # Build signature matrix
        unique_cts = sorted(set(ref_celltypes))
        X_csc = X_sparse.tocsc()
        n_ref_genes = X_csc.shape[0]

        signature = np.zeros((len(unique_cts), n_ref_genes), dtype=np.float64)
        ref_celltypes_arr = np.array(ref_celltypes)

        for i, ct in enumerate(unique_cts):
            mask = (ref_celltypes_arr == ct)
            signature[i] = X_csc[:, mask].mean(axis=1).A1

        # Align genes
        common_genes = sorted(set(genes) & set(ref_genes))
        print(f"  Common genes: {len(common_genes)}")

        spatial_idx = {g: i for i, g in enumerate(genes)}
        ref_idx = {g: i for i, g in enumerate(ref_genes)}

        Y = Y[:, [spatial_idx[g] for g in common_genes]]
        X = signature[:, [ref_idx[g] for g in common_genes]]

        # Align cell types
        ct_mapping = {ct: i for i, ct in enumerate(unique_cts)}
        gt_reordered = np.zeros((ground_truth.shape[0], len(unique_cts)))
        for j, ct in enumerate(cell_types):
            if ct in ct_mapping:
                gt_reordered[:, ct_mapping[ct]] = ground_truth[:, j]
        ground_truth = gt_reordered
        cell_types = unique_cts

    # Generate synthetic coordinates based on regions
    print(f"\nGenerating synthetic coordinates from {len(set(regions))} regions...")
    coords = assign_synthetic_coordinates(regions)

    print(f"\n  Spatial: {Y.shape[0]} spots x {Y.shape[1]} genes")
    print(f"  Reference: {X.shape[0]} cell types x {X.shape[1]} genes")
    print(f"  Regions: {sorted(set(regions))}")
    print(f"  Cell types: {cell_types[:5]}... ({len(cell_types)} total)")

    # Alpha values to test (log scale + 0)
    alphas = [0, 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1.0]

    print(f"\nTesting α values: {alphas}")

    results = []
    predictions = {}

    for alpha in alphas:
        print(f"\n  α = {alpha}...", end=" ", flush=True)

        pred = run_with_alpha(Y, X, coords, alpha)

        # Compute metrics
        rmse = compute_rmse(pred, ground_truth)

        # Compute per-celltype correlation
        correlations = []
        for k in range(pred.shape[1]):
            if np.std(pred[:, k]) > 0 and np.std(ground_truth[:, k]) > 0:
                r = np.corrcoef(pred[:, k], ground_truth[:, k])[0, 1]
                if not np.isnan(r):
                    correlations.append(r)
        mean_corr = np.mean(correlations) if correlations else 0

        # Compute spatial autocorrelation of predictions
        morans_i_list = []
        for k in range(min(5, pred.shape[1])):
            mi = compute_spatial_autocorrelation(pred[:, k], coords)
            if not np.isnan(mi):
                morans_i_list.append(mi)
        mean_morans = np.mean(morans_i_list) if morans_i_list else 0

        results.append({
            'alpha': alpha,
            'rmse': rmse,
            'mean_correlation': mean_corr,
            'mean_morans_i': mean_morans,
        })

        # Store predictions for visualization
        if alpha in [0, 0.001, 0.005, 0.01, 0.1, 0.5]:
            predictions[alpha] = pred.copy()

        print(f"RMSE={rmse:.4f}, Corr={mean_corr:.4f}, Moran's I={mean_morans:.3f}")

    # Convert to DataFrame
    results_df = pd.DataFrame(results)
    print("\n" + "=" * 70)
    print("RESULTS SUMMARY")
    print("=" * 70)
    print(results_df.to_string(index=False))

    # Find optimal alpha
    best_idx = results_df['rmse'].idxmin()
    best_alpha = results_df.loc[best_idx, 'alpha']
    best_rmse = results_df.loc[best_idx, 'rmse']

    print(f"\nBest α = {best_alpha} (RMSE = {best_rmse:.4f})")

    # Save results
    results_df.to_csv(output_dir / "alpha_sensitivity_results.csv", index=False)
    print(f"\n✓ Saved results to: {output_dir / 'alpha_sensitivity_results.csv'}")

    # Create visualization
    print("\nGenerating visualization...")

    fig = plt.figure(figsize=(16, 10))

    # Panel A: RMSE vs Alpha (U-curve)
    ax1 = fig.add_subplot(2, 3, 1)
    ax1.plot(results_df['alpha'], results_df['rmse'], 'o-', color='#e74c3c', linewidth=2, markersize=8)
    ax1.axvline(x=0.005, color='green', linestyle='--', alpha=0.7, label='Default α=0.005')
    if best_alpha != 0.005:
        ax1.axvline(x=best_alpha, color='blue', linestyle=':', alpha=0.7, label=f'Best α={best_alpha}')
    ax1.set_xscale('symlog', linthresh=0.0001)
    ax1.set_xlabel('Spatial Regularization Strength (α)')
    ax1.set_ylabel('RMSE')
    ax1.set_title('A. RMSE vs α\n(lower is better)')
    ax1.legend(loc='upper right', fontsize=9)
    ax1.grid(True, alpha=0.3)

    # Panel B: Correlation vs Alpha
    ax2 = fig.add_subplot(2, 3, 2)
    ax2.plot(results_df['alpha'], results_df['mean_correlation'], 'o-', color='#3498db', linewidth=2, markersize=8)
    ax2.axvline(x=0.005, color='green', linestyle='--', alpha=0.7, label='Default α=0.005')
    ax2.set_xscale('symlog', linthresh=0.0001)
    ax2.set_xlabel('Spatial Regularization Strength (α)')
    ax2.set_ylabel('Mean Pearson Correlation')
    ax2.set_title('B. Correlation vs α\n(higher is better)')
    ax2.legend(loc='lower right', fontsize=9)
    ax2.grid(True, alpha=0.3)

    # Panel C: Moran's I vs Alpha (Spatial Smoothness)
    ax3 = fig.add_subplot(2, 3, 3)
    ax3.plot(results_df['alpha'], results_df['mean_morans_i'], 'o-', color='#9b59b6', linewidth=2, markersize=8)
    ax3.axvline(x=0.005, color='green', linestyle='--', alpha=0.7, label='Default α=0.005')
    ax3.set_xscale('symlog', linthresh=0.0001)
    ax3.set_xlabel('Spatial Regularization Strength (α)')
    ax3.set_ylabel("Mean Moran's I")
    ax3.set_title("C. Spatial Autocorrelation vs α\n(higher = more smooth)")
    ax3.legend(loc='lower right', fontsize=9)
    ax3.grid(True, alpha=0.3)

    # Panels D-F: Visual comparison of spatial maps
    # Find a cell type with clear spatial pattern (highest variance in ground truth)
    ct_variances = np.var(ground_truth, axis=0)
    ct_idx = np.argmax(ct_variances)
    ct_name = cell_types[ct_idx] if ct_idx < len(cell_types) else f"CellType_{ct_idx}"

    visual_alphas = [0, 0.005, 0.5]
    titles = ['D. No regularization (α=0)', 'E. Default (α=0.005)', 'F. Over-regularized (α=0.5)']

    for i, (alpha, title) in enumerate(zip(visual_alphas, titles)):
        ax = fig.add_subplot(2, 3, 4 + i)

        if alpha in predictions:
            pred = predictions[alpha]
            scatter = ax.scatter(coords[:, 0], coords[:, 1],
                               c=pred[:, ct_idx], cmap='viridis',
                               s=15, alpha=0.8)
            plt.colorbar(scatter, ax=ax, shrink=0.8, label='Proportion')

            # Draw region boundaries
            unique_regions = sorted(set(regions))
            n_regions = len(unique_regions)
            for j in range(1, n_regions):
                y_boundary = j * (100.0 / n_regions)
                ax.axhline(y=y_boundary, color='red', linestyle='--', alpha=0.5, linewidth=1)

        ax.set_title(f'{title}\n{ct_name}')
        ax.set_xlabel('X coordinate')
        ax.set_ylabel('Y coordinate')
        ax.set_xlim(-5, 105)
        ax.set_ylim(-5, 105)
        ax.set_aspect('equal')

    plt.tight_layout()

    # Save figure
    output_path = fig_dir / "supplementary_alpha_sensitivity.pdf"
    plt.savefig(output_path, bbox_inches='tight', dpi=300)
    print(f"✓ Saved figure to: {output_path}")

    plt.savefig(fig_dir / "supplementary_alpha_sensitivity.png", bbox_inches='tight', dpi=150)
    print(f"✓ Saved PNG preview")

    plt.close()

    # Detailed analysis
    print("\n" + "=" * 70)
    print("DETAILED ANALYSIS")
    print("=" * 70)

    default_rmse = results_df[results_df['alpha'] == 0.005]['rmse'].values[0]

    # Find the range where RMSE is within 5% of best
    tolerance = 0.05
    good_region = results_df[results_df['rmse'] <= best_rmse * (1 + tolerance)]

    print(f"""
1. OPTIMAL α:
   Best α = {best_alpha} with RMSE = {best_rmse:.4f}
   Default α = 0.005 with RMSE = {default_rmse:.4f}
   Difference: {((default_rmse - best_rmse) / best_rmse * 100):.2f}%

2. ROBUSTNESS (α values within 5% of best RMSE):
   Good region: α ∈ [{good_region['alpha'].min()}, {good_region['alpha'].max()}]
   Number of good values: {len(good_region)}/{len(results_df)}

3. TRADE-OFF ANALYSIS:
   - α = 0 (no regularization): RMSE = {results_df[results_df['alpha'] == 0]['rmse'].values[0]:.4f}
   - α = 0.005 (default): RMSE = {default_rmse:.4f}
   - α = 0.5 (high regularization): RMSE = {results_df[results_df['alpha'] == 0.5]['rmse'].values[0]:.4f}

4. SPATIAL SMOOTHNESS (Moran's I):
   - α = 0: Moran's I = {results_df[results_df['alpha'] == 0]['mean_morans_i'].values[0]:.3f}
   - α = 0.005: Moran's I = {results_df[results_df['alpha'] == 0.005]['mean_morans_i'].values[0]:.3f}
   - α = 0.5: Moran's I = {results_df[results_df['alpha'] == 0.5]['mean_morans_i'].values[0]:.3f}

5. INTERPRETATION:
   The synthetic coordinates create KNOWN region boundaries.
   - If α is too low: predictions are noisy within regions
   - If α is too high: predictions bleed across region boundaries
   - Optimal α maintains within-region smoothness while preserving boundaries
""")

    return results_df, predictions, coords, regions


if __name__ == "__main__":
    results_df, predictions, coords, regions = main()
