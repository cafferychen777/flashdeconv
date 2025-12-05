"""
Spotless-style Benchmark Simulation
-----------------------------------
Replicating 'Rare' and 'Dominant' abundance patterns to benchmark
FlashDeconv against NNLS (Baseline).

This script implements the core logic of synthspot without requiring R,
focusing on two challenging scenarios:
1. Rare cell type detection (AUPR metric)
2. Dominant cell type estimation (RMSE metric)

Reference: Sang-aram et al., eLife 2024 - Spotless benchmark
"""
import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import mean_squared_error, average_precision_score
from scipy.optimize import nnls
from scipy.spatial.distance import jensenshannon
from scipy.stats import pearsonr
import time
import warnings
warnings.filterwarnings('ignore')

# Add FlashDeconv to path
sys.path.insert(0, '/Users/apple/Research/FlashDeconv')
from flashdeconv import FlashDeconv


# ============================================
# Data Generation (Synthspot-style)
# ============================================

def get_simulated_reference(n_types=10, n_genes=2000, seed=42):
    """
    Generate simulated scRNA-seq reference signatures.

    Mimics realistic gene expression patterns:
    - Log-normal distribution (most genes low, few high)
    - Sparse expression (many zeros)
    - Cell-type specific markers
    """
    np.random.seed(seed)

    # Base expression: log-normal with sparsity
    ref = np.random.lognormal(mean=1.0, sigma=1.5, size=(n_types, n_genes))

    # Sparsify: set low values to 0
    ref[ref < 1.5] = 0

    # Add cell-type specific markers (20 markers per type)
    n_markers = 20
    for k in range(n_types):
        marker_idx = np.random.choice(n_genes, n_markers, replace=False)
        ref[k, marker_idx] *= np.random.uniform(5, 15, n_markers)

    # Normalize rows to have similar total expression
    row_sums = ref.sum(axis=1, keepdims=True)
    ref = ref / row_sums * 10000  # Normalize to ~10k total counts

    return ref


def generate_spatial_coords(n_spots, grid_size=None):
    """Generate spatial coordinates on a grid."""
    if grid_size is None:
        grid_size = int(np.ceil(np.sqrt(n_spots)))

    coords = []
    for i in range(n_spots):
        x = i % grid_size
        y = i // grid_size
        coords.append([x, y])

    return np.array(coords, dtype=float)


def generate_spatial_data(ref, n_spots, mode='uniform', noise_level=0.1,
                          spatial_pattern=False, seed=None):
    """
    Core generator: Simulate Spotless-style synthetic spatial data.

    Parameters
    ----------
    ref : ndarray (n_types, n_genes)
        Reference cell type signatures
    n_spots : int
        Number of spatial spots to generate
    mode : str
        Abundance pattern:
        - 'uniform': Dirichlet(1,1,...,1) - all cell types equally likely
        - 'dominant': One cell type dominates (60-90%)
        - 'rare': One cell type is very rare (0-3%)
        - 'regional_rare': Rare cell type localized to a region
    noise_level : float
        Poisson noise scaling factor
    spatial_pattern : bool
        If True, add spatial autocorrelation to proportions
    seed : int
        Random seed

    Returns
    -------
    Y : ndarray (n_spots, n_genes)
        Synthetic spatial count matrix
    beta_true : ndarray (n_spots, n_types)
        Ground truth cell type proportions
    """
    if seed is not None:
        np.random.seed(seed)

    n_types, n_genes = ref.shape

    # 1. Generate true proportions based on mode
    if mode == 'uniform':
        # Uniform Dirichlet - all cell types equally likely
        beta = np.random.dirichlet(np.ones(n_types), size=n_spots)

    elif mode == 'dominant':
        # Dominant pattern: Type 0 takes 60-90%
        # Use skewed Dirichlet: alpha[0] >> alpha[others]
        alpha = np.ones(n_types) * 0.5
        alpha[0] = 20.0  # Dominant cell type
        beta = np.random.dirichlet(alpha, size=n_spots)

    elif mode == 'rare':
        # Rare pattern: Type 0 is very rare (1-3% on average)
        # Most spots have 0, few spots have small amounts
        beta = np.random.dirichlet(np.ones(n_types) * 2, size=n_spots)

        # Generate sparse rare cell proportions
        # Use gamma distribution for sparsity
        rare_prop = np.random.gamma(shape=0.3, scale=0.02, size=n_spots)
        rare_prop = np.clip(rare_prop, 0, 0.1)

        # Zero out most spots (make it sparse)
        mask = np.random.random(n_spots) > 0.3  # 70% spots have zero
        rare_prop[mask] = 0

        # Reassign proportions
        beta[:, 0] = rare_prop
        remaining = 1 - rare_prop
        other_props = beta[:, 1:]
        other_sums = other_props.sum(axis=1, keepdims=True)
        other_sums[other_sums == 0] = 1
        beta[:, 1:] = other_props / other_sums * remaining[:, None]

    elif mode == 'regional_rare':
        # Rare cell type localized to specific spatial region
        beta = np.random.dirichlet(np.ones(n_types) * 2, size=n_spots)

        # Get coordinates for spatial pattern
        coords = generate_spatial_coords(n_spots)
        grid_size = int(np.ceil(np.sqrt(n_spots)))

        # Define a "hot region" in the corner
        center_x, center_y = grid_size * 0.8, grid_size * 0.8

        rare_prop = np.zeros(n_spots)
        for i in range(n_spots):
            dist = np.sqrt((coords[i, 0] - center_x)**2 + (coords[i, 1] - center_y)**2)
            if dist < grid_size * 0.3:
                # In hot region: rare cell present
                rare_prop[i] = np.random.beta(2, 10) * 0.15  # 0-15%

        # Reassign proportions
        beta[:, 0] = rare_prop
        remaining = 1 - rare_prop
        other_props = beta[:, 1:]
        other_sums = other_props.sum(axis=1, keepdims=True)
        other_sums[other_sums == 0] = 1
        beta[:, 1:] = other_props / other_sums * remaining[:, None]

    else:
        raise ValueError(f"Unknown mode: {mode}")

    # 2. Add spatial autocorrelation (optional)
    if spatial_pattern:
        coords = generate_spatial_coords(n_spots)
        from scipy.ndimage import gaussian_filter
        grid_size = int(np.ceil(np.sqrt(n_spots)))

        for k in range(n_types):
            # Reshape to grid, smooth, reshape back
            prop_grid = beta[:, k].reshape(grid_size, -1)[:grid_size, :grid_size]
            smoothed = gaussian_filter(prop_grid, sigma=1.5)
            beta[:grid_size*grid_size, k] = smoothed.flatten()[:n_spots]

        # Renormalize
        beta = beta / beta.sum(axis=1, keepdims=True)

    # 3. Generate synthetic counts: Y = Beta @ Ref + Noise
    expected_counts = beta @ ref

    # Simulate variable sequencing depth
    depths = np.random.normal(20000, 5000, size=n_spots)
    depths = np.maximum(depths, 5000)  # Minimum depth

    # Scale and add Poisson noise
    Y = np.zeros((n_spots, n_genes))
    for i in range(n_spots):
        profile = expected_counts[i, :]
        if profile.sum() > 0:
            profile = profile / profile.sum() * depths[i]
        Y[i, :] = np.random.poisson(np.maximum(profile, 0))

    return Y.astype(float), beta


# ============================================
# Deconvolution Methods
# ============================================

def run_nnls(Y, X):
    """
    Baseline NNLS deconvolution.

    Parameters
    ----------
    Y : ndarray (n_spots, n_genes)
        Spatial expression matrix
    X : ndarray (n_types, n_genes)
        Reference signatures

    Returns
    -------
    beta : ndarray (n_spots, n_types)
        Estimated proportions
    """
    n_spots = Y.shape[0]
    n_types = X.shape[0]
    beta_pred = np.zeros((n_spots, n_types))

    for i in range(n_spots):
        beta_pred[i], _ = nnls(X.T, Y[i])

    # Normalize to proportions
    sums = beta_pred.sum(axis=1, keepdims=True)
    sums[sums == 0] = 1
    return beta_pred / sums


def run_flashdeconv(Y, X, coords, lambda_spatial=0.1, preprocess=False):
    """
    Run FlashDeconv deconvolution.

    Parameters
    ----------
    Y : ndarray (n_spots, n_genes)
        Spatial expression matrix
    X : ndarray (n_types, n_genes)
        Reference signatures
    coords : ndarray (n_spots, 2)
        Spatial coordinates
    lambda_spatial : float
        Spatial regularization strength
    preprocess : bool
        Whether to apply VST preprocessing

    Returns
    -------
    beta : ndarray (n_spots, n_types)
        Estimated proportions
    """
    model = FlashDeconv(
        sketch_dim=256,
        lambda_spatial=lambda_spatial,
        rho_sparsity=0.01,
        preprocess=preprocess,
        n_hvg=min(1000, Y.shape[1]),
        max_iter=100,
        verbose=False
    )

    proportions = model.fit_transform(Y, X, coords)
    return proportions


# ============================================
# Evaluation Metrics
# ============================================

def compute_all_metrics(pred, true, target_idx=0):
    """
    Compute comprehensive metrics for deconvolution evaluation.

    Parameters
    ----------
    pred : ndarray (n_spots, n_types)
        Predicted proportions
    true : ndarray (n_spots, n_types)
        Ground truth proportions
    target_idx : int
        Index of target cell type for AUPR

    Returns
    -------
    metrics : dict
        Dictionary of computed metrics
    """
    n_spots, n_types = true.shape

    # 1. Pearson correlation (overall)
    corr, _ = pearsonr(pred.flatten(), true.flatten())

    # 2. Per-cell-type correlation (mean)
    per_type_corr = []
    for k in range(n_types):
        if true[:, k].std() > 1e-6 and pred[:, k].std() > 1e-6:
            r, _ = pearsonr(pred[:, k], true[:, k])
            per_type_corr.append(r)
    mean_type_corr = np.mean(per_type_corr) if per_type_corr else 0

    # 3. RMSE (per spot, then average)
    rmse_per_spot = np.sqrt(np.mean((pred - true)**2, axis=1))
    rmse = np.mean(rmse_per_spot)

    # 4. RMSE for target cell type
    rmse_target = np.sqrt(mean_squared_error(true[:, target_idx], pred[:, target_idx]))

    # 5. JSD (Jensen-Shannon Divergence)
    eps = 1e-10
    pred_norm = pred + eps
    pred_norm = pred_norm / pred_norm.sum(axis=1, keepdims=True)
    true_norm = true + eps
    true_norm = true_norm / true_norm.sum(axis=1, keepdims=True)

    jsd_per_spot = []
    for i in range(n_spots):
        jsd = jensenshannon(pred_norm[i], true_norm[i])**2
        jsd_per_spot.append(jsd)
    jsd = np.mean(jsd_per_spot)

    # 6. AUPR for target cell type (binary: present vs absent)
    y_true_binary = (true[:, target_idx] > 0.005).astype(int)
    if len(np.unique(y_true_binary)) == 2:
        aupr = average_precision_score(y_true_binary, pred[:, target_idx])
    else:
        aupr = np.nan

    # 7. False Positive Rate for target cell type
    # Predicted present but actually absent
    pred_present = pred[:, target_idx] > 0.01
    true_absent = true[:, target_idx] < 0.005
    if true_absent.sum() > 0:
        fpr = (pred_present & true_absent).sum() / true_absent.sum()
    else:
        fpr = np.nan

    return {
        'corr': corr,
        'mean_type_corr': mean_type_corr,
        'rmse': rmse,
        'rmse_target': rmse_target,
        'jsd': jsd,
        'aupr': aupr,
        'fpr': fpr
    }


# ============================================
# Main Benchmark
# ============================================

def run_benchmark(n_replicates=5, n_spots=500, n_types=10, n_genes=2000):
    """
    Run full Spotless-style benchmark.

    Parameters
    ----------
    n_replicates : int
        Number of replicates per condition
    n_spots : int
        Number of spatial spots
    n_types : int
        Number of cell types
    n_genes : int
        Number of genes

    Returns
    -------
    results : DataFrame
        Benchmark results
    """
    print("=" * 70)
    print("Spotless-style Benchmark: FlashDeconv vs NNLS")
    print("=" * 70)
    print(f"Settings: {n_spots} spots, {n_types} cell types, {n_genes} genes")
    print(f"Replicates: {n_replicates}")
    print()

    # Generate reference (same for all replicates)
    ref = get_simulated_reference(n_types=n_types, n_genes=n_genes, seed=42)
    print(f"Reference shape: {ref.shape}")
    print(f"Reference sparsity: {(ref == 0).mean():.1%}")
    print()

    modes = ['uniform', 'rare', 'dominant', 'regional_rare']
    all_results = []

    for mode in modes:
        print(f"\n{'='*50}")
        print(f"Mode: {mode.upper()}")
        print('='*50)

        for rep in range(n_replicates):
            seed = 1000 + rep

            # Generate synthetic data
            Y, beta_true = generate_spatial_data(
                ref, n_spots=n_spots, mode=mode, seed=seed
            )
            coords = generate_spatial_coords(n_spots)

            # Print data summary
            if rep == 0:
                target_mean = beta_true[:, 0].mean()
                target_max = beta_true[:, 0].max()
                target_nonzero = (beta_true[:, 0] > 0.005).mean()
                print(f"Target cell type (idx=0):")
                print(f"  Mean proportion: {target_mean:.3f}")
                print(f"  Max proportion: {target_max:.3f}")
                print(f"  Present in: {target_nonzero:.1%} of spots")

            # Run NNLS (baseline)
            t0 = time.time()
            beta_nnls = run_nnls(Y, ref)
            time_nnls = time.time() - t0

            # Run FlashDeconv with low lambda (for synthetic data)
            t0 = time.time()
            beta_flash = run_flashdeconv(Y, ref, coords, lambda_spatial=0.1, preprocess=False)
            time_flash = time.time() - t0

            # Compute metrics
            metrics_nnls = compute_all_metrics(beta_nnls, beta_true, target_idx=0)
            metrics_flash = compute_all_metrics(beta_flash, beta_true, target_idx=0)

            # Store results
            for method, metrics, runtime in [
                ('NNLS', metrics_nnls, time_nnls),
                ('FlashDeconv', metrics_flash, time_flash)
            ]:
                result = {
                    'mode': mode,
                    'replicate': rep,
                    'method': method,
                    'time': runtime,
                    **metrics
                }
                all_results.append(result)

            if rep == 0:
                print(f"\nRep {rep} Results:")
                print(f"  NNLS:       Corr={metrics_nnls['corr']:.3f}, RMSE={metrics_nnls['rmse']:.4f}, "
                      f"AUPR={metrics_nnls['aupr']:.3f}, FPR={metrics_nnls['fpr']:.3f}")
                print(f"  FlashDeconv: Corr={metrics_flash['corr']:.3f}, RMSE={metrics_flash['rmse']:.4f}, "
                      f"AUPR={metrics_flash['aupr']:.3f}, FPR={metrics_flash['fpr']:.3f}")

    # Create results DataFrame
    df = pd.DataFrame(all_results)

    return df


def summarize_results(df):
    """Summarize benchmark results by mode and method."""
    print("\n" + "=" * 70)
    print("SUMMARY RESULTS (Mean ± Std)")
    print("=" * 70)

    summary = df.groupby(['mode', 'method']).agg({
        'corr': ['mean', 'std'],
        'rmse': ['mean', 'std'],
        'aupr': ['mean', 'std'],
        'fpr': ['mean', 'std'],
        'jsd': ['mean', 'std'],
        'time': ['mean']
    }).round(4)

    print(summary)
    return summary


def plot_results(df, save_path='validation/spotless_sim_results.png'):
    """Create visualization of benchmark results."""
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))

    # Colors
    palette = {'NNLS': '#1f77b4', 'FlashDeconv': '#ff7f0e'}

    # Plot 1: Overall Correlation by Mode
    ax = axes[0, 0]
    sns.barplot(data=df, x='mode', y='corr', hue='method', ax=ax, palette=palette)
    ax.set_title('Overall Correlation (Higher is Better)')
    ax.set_xlabel('Abundance Pattern')
    ax.set_ylabel('Pearson r')
    ax.legend(title='Method')

    # Plot 2: RMSE by Mode
    ax = axes[0, 1]
    sns.barplot(data=df, x='mode', y='rmse', hue='method', ax=ax, palette=palette)
    ax.set_title('RMSE (Lower is Better)')
    ax.set_xlabel('Abundance Pattern')
    ax.set_ylabel('RMSE')
    ax.legend(title='Method')

    # Plot 3: JSD by Mode
    ax = axes[0, 2]
    sns.barplot(data=df, x='mode', y='jsd', hue='method', ax=ax, palette=palette)
    ax.set_title('Jensen-Shannon Divergence (Lower is Better)')
    ax.set_xlabel('Abundance Pattern')
    ax.set_ylabel('JSD')
    ax.legend(title='Method')

    # Plot 4: AUPR for Rare Cell Types
    ax = axes[1, 0]
    df_rare = df[df['mode'].isin(['rare', 'regional_rare'])]
    sns.barplot(data=df_rare, x='mode', y='aupr', hue='method', ax=ax, palette=palette)
    ax.set_title('AUPR for Rare Cell Detection\n(Higher is Better)')
    ax.set_xlabel('Abundance Pattern')
    ax.set_ylabel('AUPR')
    ax.legend(title='Method')

    # Plot 5: False Positive Rate
    ax = axes[1, 1]
    sns.barplot(data=df_rare, x='mode', y='fpr', hue='method', ax=ax, palette=palette)
    ax.set_title('False Positive Rate\n(Lower is Better)')
    ax.set_xlabel('Abundance Pattern')
    ax.set_ylabel('FPR')
    ax.legend(title='Method')

    # Plot 6: Runtime
    ax = axes[1, 2]
    sns.barplot(data=df, x='mode', y='time', hue='method', ax=ax, palette=palette)
    ax.set_title('Runtime (Lower is Better)')
    ax.set_xlabel('Abundance Pattern')
    ax.set_ylabel('Time (seconds)')
    ax.legend(title='Method')

    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    print(f"\nPlot saved to: {save_path}")
    plt.close()


def main():
    """Run the full benchmark."""
    # Run benchmark
    df = run_benchmark(
        n_replicates=5,
        n_spots=500,
        n_types=10,
        n_genes=2000
    )

    # Save raw results
    df.to_csv('validation/spotless_sim_raw_results.csv', index=False)
    print(f"\nRaw results saved to: validation/spotless_sim_raw_results.csv")

    # Summarize
    summary = summarize_results(df)

    # Plot
    plot_results(df)

    # Print key findings
    print("\n" + "=" * 70)
    print("KEY FINDINGS")
    print("=" * 70)

    # Compare methods on rare cell detection
    rare_df = df[df['mode'] == 'rare']
    nnls_aupr = rare_df[rare_df['method'] == 'NNLS']['aupr'].mean()
    flash_aupr = rare_df[rare_df['method'] == 'FlashDeconv']['aupr'].mean()
    print(f"\n1. Rare Cell Detection (AUPR):")
    print(f"   NNLS: {nnls_aupr:.3f}")
    print(f"   FlashDeconv: {flash_aupr:.3f}")
    if flash_aupr > nnls_aupr:
        print(f"   -> FlashDeconv is BETTER by {(flash_aupr - nnls_aupr):.3f}")

    # Compare on dominant cell type
    dom_df = df[df['mode'] == 'dominant']
    nnls_rmse = dom_df[dom_df['method'] == 'NNLS']['rmse_target'].mean()
    flash_rmse = dom_df[dom_df['method'] == 'FlashDeconv']['rmse_target'].mean()
    print(f"\n2. Dominant Cell Type (RMSE):")
    print(f"   NNLS: {nnls_rmse:.4f}")
    print(f"   FlashDeconv: {flash_rmse:.4f}")
    if flash_rmse < nnls_rmse:
        print(f"   -> FlashDeconv is BETTER by {(nnls_rmse - flash_rmse):.4f}")

    # Compare FPR
    nnls_fpr = rare_df[rare_df['method'] == 'NNLS']['fpr'].mean()
    flash_fpr = rare_df[rare_df['method'] == 'FlashDeconv']['fpr'].mean()
    print(f"\n3. False Positive Rate (Rare Mode):")
    print(f"   NNLS: {nnls_fpr:.3f}")
    print(f"   FlashDeconv: {flash_fpr:.3f}")
    if flash_fpr < nnls_fpr:
        print(f"   -> FlashDeconv has {(nnls_fpr - flash_fpr)*100:.1f}% lower FPR")

    # Speed comparison
    nnls_time = df[df['method'] == 'NNLS']['time'].mean()
    flash_time = df[df['method'] == 'FlashDeconv']['time'].mean()
    print(f"\n4. Speed:")
    print(f"   NNLS: {nnls_time:.3f}s")
    print(f"   FlashDeconv: {flash_time:.3f}s")
    print(f"   Speedup: {nnls_time/flash_time:.1f}x" if flash_time < nnls_time else
          f"   Overhead: {flash_time/nnls_time:.1f}x")


if __name__ == '__main__':
    main()
