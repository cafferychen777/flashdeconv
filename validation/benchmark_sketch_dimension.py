"""
Benchmark FlashDeconv sketch dimension sensitivity analysis.

Tests the effect of sketch dimension on accuracy and runtime
using Silver Standard datasets (Brain Cortex and Kidney).

Author: Claude
Date: 2025-12-02
"""

import numpy as np
import pandas as pd
import time
from pathlib import Path
from scipy.io import mmread
from scipy import sparse
import sys
import gc
import argparse

# Add FlashDeconv to path
sys.path.insert(0, '/Users/apple/Research/FlashDeconv')

# Import FlashDeconv API
from flashdeconv import FlashDeconv
from flashdeconv.utils.metrics import compute_rmse, compute_correlation


def load_silver_standard(data_dir, tissue_id, replicate_id):
    """
    Load a Silver Standard dataset.

    Parameters:
    -----------
    data_dir : Path
        Directory containing the data files
    tissue_id : int
        Tissue ID (1=Brain Cortex, 5=Kidney)
    replicate_id : int
        Replicate ID (1-11 for Brain Cortex, 1-9 for Kidney)

    Returns:
    --------
    Y : ndarray
        Spatial count matrix (n_spots x n_genes)
    X : ndarray
        Reference signature matrix (n_cell_types x n_genes)
    gt_proportions : ndarray
        Ground truth proportions (n_spots x n_cell_types)
    coords : ndarray
        Synthetic coordinates (n_spots x 2)
    cell_type_names : list
        Cell type names
    """
    prefix = f"silver_{tissue_id}_{replicate_id}"
    ref_prefix = f"reference_{tissue_id}"

    # Load spatial data
    Y_sparse = mmread(data_dir / f"{prefix}_counts.mtx")
    Y = Y_sparse.toarray().T  # Convert to (n_spots x n_genes)

    with open(data_dir / f"{prefix}_genes.txt") as f:
        genes = [line.strip() for line in f]

    with open(data_dir / f"{prefix}_spots.txt") as f:
        spots = [line.strip() for line in f]

    # Load ground truth proportions
    gt_df = pd.read_csv(data_dir / f"{prefix}_proportions.csv", index_col=0)
    gt_proportions = gt_df.values
    cell_type_names = list(gt_df.columns)

    # Load reference
    X_sparse = mmread(data_dir / f"{ref_prefix}_counts.mtx")

    with open(data_dir / f"{ref_prefix}_celltypes.txt") as f:
        ref_celltypes = [line.strip() for line in f]

    with open(data_dir / f"{ref_prefix}_genes.txt") as f:
        ref_genes = [line.strip() for line in f]

    # Create signature matrix by averaging cells per type
    unique_cts = sorted(set(ref_celltypes))
    ct_to_idx = {ct: i for i, ct in enumerate(unique_cts)}

    X_csc = X_sparse.tocsc()
    X = np.zeros((len(unique_cts), len(ref_genes)), dtype=np.float32)

    ref_celltypes_arr = np.array(ref_celltypes)
    for i, ct in enumerate(unique_cts):
        mask = (ref_celltypes_arr == ct)
        X[i] = X_csc[:, mask].mean(axis=1).A1

    # Generate synthetic coordinates (grid layout)
    n_spots = Y.shape[0]
    side = int(np.ceil(np.sqrt(n_spots)))
    coords = np.array([[i % side, i // side] for i in range(n_spots)], dtype=float)

    # Align genes between Y and X
    gene_set = set(genes) & set(ref_genes)
    gene_idx_Y = [genes.index(g) for g in gene_set if g in genes]
    gene_idx_X = [ref_genes.index(g) for g in gene_set if g in ref_genes]

    Y_aligned = Y[:, gene_idx_Y]
    X_aligned = X[:, gene_idx_X]

    # Align cell types between ground truth and reference
    ct_order = [cell_type_names.index(ct) if ct in cell_type_names else -1
                for ct in unique_cts]
    valid_cts = [i for i, idx in enumerate(ct_order) if idx >= 0]

    X_final = X_aligned[valid_cts]
    gt_final = gt_proportions[:, [ct_order[i] for i in valid_cts]]
    final_ct_names = [unique_cts[i] for i in valid_cts]

    return Y_aligned, X_final, gt_final, coords, final_ct_names


def run_benchmark(tissue_name, tissue_id, replicates, sketch_dims, data_dir, output_dir):
    """
    Run benchmark for a tissue across multiple sketch dimensions and replicates.
    """
    results = []

    print(f"\n{'='*60}")
    print(f"Benchmarking {tissue_name} (tissue_id={tissue_id})")
    print(f"Sketch dimensions: {sketch_dims}")
    print(f"Replicates: {replicates}")
    print(f"{'='*60}")

    for rep_id in replicates:
        print(f"\n[Replicate {rep_id}] Loading data...")

        try:
            Y, X, gt, coords, ct_names = load_silver_standard(data_dir, tissue_id, rep_id)
        except Exception as e:
            print(f"  Error loading replicate {rep_id}: {e}")
            continue

        print(f"  Y: {Y.shape}, X: {X.shape}, GT: {gt.shape}")
        print(f"  Cell types: {len(ct_names)}")

        for d in sketch_dims:
            print(f"\n  [d={d}] Running FlashDeconv...")

            # Run FlashDeconv
            start_time = time.time()

            model = FlashDeconv(
                sketch_dim=d,
                lambda_spatial=0.0,  # No spatial smoothing for synthetic data
                rho_sparsity=0.01,
                n_hvg=2000,
                n_markers_per_type=50,
                max_iter=100,
                tol=1e-4,
                preprocess="log_cpm",
                random_state=42,
                verbose=False,
            )

            try:
                pred = model.fit_transform(Y, X, coords)
            except Exception as e:
                print(f"    Error: {e}")
                continue

            elapsed = time.time() - start_time

            # Compute metrics
            rmse_val = compute_rmse(pred, gt)
            corr_val = compute_correlation(pred, gt, method="pearson")

            results.append({
                'tissue': tissue_name,
                'tissue_id': tissue_id,
                'replicate': rep_id,
                'sketch_dim': d,
                'n_spots': Y.shape[0],
                'n_genes': Y.shape[1],
                'n_cell_types': X.shape[0],
                'rmse': rmse_val,
                'pearson': corr_val,
                'runtime_sec': elapsed,
            })

            print(f"    RMSE: {rmse_val:.4f}, Pearson: {corr_val:.4f}, Time: {elapsed:.2f}s")

            # Clean up
            del model, pred
            gc.collect()

    return pd.DataFrame(results)


def main():
    parser = argparse.ArgumentParser(description='Benchmark sketch dimension sensitivity')
    parser.add_argument('--tissue', type=str, choices=['brain', 'kidney', 'both'],
                        default='both', help='Which tissue to benchmark')
    parser.add_argument('--quick', action='store_true',
                        help='Quick mode: fewer replicates and dimensions')
    args = parser.parse_args()

    data_dir = Path("/Users/apple/Research/FlashDeconv/validation/benchmark_data/converted")
    output_dir = Path("/Users/apple/Research/FlashDeconv/validation/data")
    output_dir.mkdir(parents=True, exist_ok=True)

    # Configuration
    if args.quick:
        sketch_dims = [64, 256, 512, 1024]
        brain_replicates = [1, 2, 3]
        kidney_replicates = [1, 2, 3]
    else:
        sketch_dims = [64, 128, 256, 512, 1024, 2048]
        brain_replicates = [1, 2, 3, 4, 5]
        kidney_replicates = [1, 2, 3, 4, 5]

    all_results = []

    # Brain Cortex (silver_1)
    if args.tissue in ['brain', 'both']:
        brain_results = run_benchmark(
            tissue_name="Brain Cortex",
            tissue_id=1,
            replicates=brain_replicates,
            sketch_dims=sketch_dims,
            data_dir=data_dir,
            output_dir=output_dir,
        )
        all_results.append(brain_results)

    # Kidney (silver_5)
    if args.tissue in ['kidney', 'both']:
        kidney_results = run_benchmark(
            tissue_name="Kidney",
            tissue_id=5,
            replicates=kidney_replicates,
            sketch_dims=sketch_dims,
            data_dir=data_dir,
            output_dir=output_dir,
        )
        all_results.append(kidney_results)

    # Combine and save
    if all_results:
        combined = pd.concat(all_results, ignore_index=True)
        output_path = output_dir / "sketch_dimension_sensitivity.csv"
        combined.to_csv(output_path, index=False)
        print(f"\n✓ Saved results to: {output_path}")

        # Print summary
        print("\n" + "="*60)
        print("SUMMARY: Mean metrics by sketch dimension")
        print("="*60)

        for tissue in combined['tissue'].unique():
            tissue_data = combined[combined['tissue'] == tissue]
            print(f"\n{tissue}:")
            summary = tissue_data.groupby('sketch_dim').agg({
                'rmse': ['mean', 'std'],
                'pearson': ['mean', 'std'],
                'runtime_sec': ['mean', 'std'],
            }).round(4)
            print(summary)


if __name__ == "__main__":
    main()
