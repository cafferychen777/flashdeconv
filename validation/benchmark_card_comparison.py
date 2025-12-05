"""
Benchmark: FlashDeconv vs CARD comparison
-----------------------------------------
Re-test FlashDeconv on CARD's pancreatic cancer example data
with the improved preprocessing method (log_cpm).
"""
import sys
import time
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import pearsonr

sys.path.insert(0, '/Users/apple/Research/FlashDeconv')
from flashdeconv import FlashDeconv


def log(msg):
    """Print with flush for real-time output."""
    print(msg, flush=True)


def main():
    log("=" * 70)
    log("FlashDeconv vs CARD Comparison on Pancreatic Cancer Data")
    log("=" * 70)

    # Load data
    data_dir = "/Users/apple/Research/SpatialTrans_MCP/chatspatial/data/examples/card/"

    log("\n[1] Loading spatial data...")
    spatial = sc.read_h5ad(f"{data_dir}/card_spatial.h5ad")
    log(f"    Spatial: {spatial.shape[0]} spots x {spatial.shape[1]} genes")

    log("\n[2] Loading reference data...")
    reference = sc.read_h5ad(f"{data_dir}/card_reference.h5ad")
    log(f"    Reference: {reference.shape[0]} cells x {reference.shape[1]} genes")

    # Get cell type info
    if 'cellType' in reference.obs.columns:
        ct_col = 'cellType'
    elif 'cell_type' in reference.obs.columns:
        ct_col = 'cell_type'
    else:
        ct_col = reference.obs.columns[0]

    cell_types = reference.obs[ct_col].unique()
    log(f"    Cell types ({len(cell_types)}): {list(cell_types)[:5]}...")

    # Find common genes between spatial and reference
    log("\n[3] Aligning genes...")
    spatial_genes = set(spatial.var_names)
    ref_genes = set(reference.var_names)
    common_genes = list(spatial_genes & ref_genes)
    common_genes.sort()  # Ensure consistent order
    log(f"    Spatial genes: {len(spatial_genes)}")
    log(f"    Reference genes: {len(ref_genes)}")
    log(f"    Common genes: {len(common_genes)}")

    # Subset to common genes
    spatial_common = spatial[:, common_genes]
    reference_common = reference[:, common_genes]

    # Build reference matrix (mean expression per cell type)
    log("\n[4] Building reference matrix...")
    X_ref = []
    ct_names = []
    for ct in cell_types:
        mask = (reference_common.obs[ct_col] == ct).values
        if isinstance(reference_common.X, np.ndarray):
            mean_expr = reference_common.X[mask].mean(axis=0)
        else:
            mean_expr = np.asarray(reference_common.X[mask].mean(axis=0)).flatten()
        X_ref.append(mean_expr)
        ct_names.append(ct)

    X_ref = np.array(X_ref)
    log(f"    Reference matrix: {X_ref.shape[0]} types x {X_ref.shape[1]} genes")

    # Prepare spatial data
    if isinstance(spatial_common.X, np.ndarray):
        Y = spatial_common.X
    else:
        Y = np.asarray(spatial_common.X.todense())

    # Get spatial coordinates
    if 'spatial' in spatial_common.obsm:
        coords = spatial_common.obsm['spatial']
    else:
        # Use x, y from obs if available
        coords = spatial_common.obs[['x', 'y']].values if 'x' in spatial_common.obs else None

    log(f"    Y shape: {Y.shape}")
    log(f"    X shape: {X_ref.shape}")
    if coords is not None:
        log(f"    Coords shape: {coords.shape}")

    # Load CARD results for comparison
    log("\n[5] Loading CARD results...")
    card_props = pd.read_csv("/Users/apple/Research/FlashDeconv/validation/_deprecated/card_proportions.csv", index_col=0)
    log(f"    CARD proportions: {card_props.shape}")

    # Run FlashDeconv with improved preprocessing
    log("\n[6] Running FlashDeconv (preprocess='log_cpm')...")
    model = FlashDeconv(
        sketch_dim=512,
        lambda_spatial="auto",  # Auto-tune
        n_hvg=2000,
        max_iter=100,
        preprocess="log_cpm",  # The improved preprocessing
        verbose=True,
        random_state=42,
    )

    t0 = time.time()
    beta_flash = model.fit_transform(Y, X_ref, coords)
    flash_time = time.time() - t0
    log(f"\n    FlashDeconv time: {flash_time:.2f}s")

    # Create FlashDeconv results DataFrame
    flash_props = pd.DataFrame(beta_flash, columns=ct_names, index=spatial.obs_names)

    # Compare with CARD
    log("\n" + "=" * 70)
    log("COMPARISON: FlashDeconv vs CARD")
    log("=" * 70)

    # Align cell types
    common_cts = set(flash_props.columns) & set(card_props.columns)
    log(f"\nCommon cell types: {len(common_cts)}")

    results = []
    for ct in sorted(common_cts):
        flash_vals = flash_props[ct].values
        card_vals = card_props[ct].values

        # Ensure same length
        min_len = min(len(flash_vals), len(card_vals))
        flash_vals = flash_vals[:min_len]
        card_vals = card_vals[:min_len]

        r, _ = pearsonr(flash_vals, card_vals)
        results.append({
            'cell_type': ct,
            'card_mean': card_vals.mean(),
            'flash_mean': flash_vals.mean(),
            'pearson_r': r
        })

    results_df = pd.DataFrame(results)
    results_df = results_df.sort_values('pearson_r', ascending=False)

    log("\n" + "-" * 70)
    log(f"{'Cell Type':<45} {'CARD%':>8} {'Flash%':>8} {'r':>8}")
    log("-" * 70)
    for _, row in results_df.iterrows():
        log(f"{row['cell_type']:<45} {row['card_mean']*100:>7.2f}% {row['flash_mean']*100:>7.2f}% {row['pearson_r']:>8.3f}")

    # Summary statistics
    mean_r = results_df['pearson_r'].mean()
    median_r = results_df['pearson_r'].median()
    high_corr = (results_df['pearson_r'] > 0.5).sum()

    log("\n" + "=" * 70)
    log("SUMMARY")
    log("=" * 70)
    log(f"Mean Pearson r:    {mean_r:.4f}")
    log(f"Median Pearson r:  {median_r:.4f}")
    log(f"Cell types with r > 0.5: {high_corr}/{len(results_df)}")
    log(f"FlashDeconv time:  {flash_time:.2f}s")

    # Compare with old results
    old_results = pd.read_csv("/Users/apple/Research/FlashDeconv/validation/_deprecated/card_vs_flashdeconv_comparison.csv")
    old_mean_r = old_results['pearson_r'].mean()

    log(f"\n--- Improvement over old preprocessing ---")
    log(f"Old mean r:        {old_mean_r:.4f}")
    log(f"New mean r:        {mean_r:.4f}")
    log(f"Improvement:       {mean_r - old_mean_r:+.4f}")

    # Save results
    output_path = "/Users/apple/Research/FlashDeconv/validation/card_comparison_new.csv"
    results_df.to_csv(output_path, index=False)
    log(f"\nResults saved to: {output_path}")

    return results_df


if __name__ == "__main__":
    main()
