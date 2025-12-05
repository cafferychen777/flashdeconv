"""
Fair benchmark of FlashDeconv on ALL Spotless Silver Standards.

Tests on all 6 datasets:
1. brain_cortex (11 patterns)
2. cerebellum_cell (9 patterns)
3. cerebellum_nucleus (9 patterns)
4. hippocampus (9 patterns)
5. kidney (9 patterns)
6. scc_p5 (9 patterns)

Total: 56 dataset-pattern combinations
"""

import os
import sys
import time
import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.stats import pearsonr
from scipy.spatial.distance import jensenshannon
from sklearn.metrics import precision_recall_curve, auc
import warnings
warnings.filterwarnings('ignore')

sys.path.insert(0, '/Users/apple/Research/FlashDeconv')
from flashdeconv import FlashDeconv


def load_reference(prefix):
    """Load converted reference data."""
    counts = mmread(f"{prefix}_counts.mtx").T.tocsr()  # cells x genes
    genes = pd.read_csv(f"{prefix}_genes.txt", header=None)[0].values
    cells = pd.read_csv(f"{prefix}_cells.txt", header=None)[0].values
    celltypes = pd.read_csv(f"{prefix}_celltypes.txt", header=None)[0].values

    return {
        'counts': counts,
        'genes': genes,
        'cells': cells,
        'celltypes': celltypes
    }


def load_test_data(prefix):
    """Load converted test data."""
    counts = mmread(f"{prefix}_counts.mtx").T.tocsr()  # spots x genes
    genes = pd.read_csv(f"{prefix}_genes.txt", header=None)[0].values
    props = pd.read_csv(f"{prefix}_proportions.csv", index_col=0)

    # Filter numeric columns only
    props = props.select_dtypes(include=[np.number])

    # Generate grid coordinates
    n_spots = counts.shape[0]
    grid_size = int(np.ceil(np.sqrt(n_spots)))
    coords = np.array([[i % grid_size, i // grid_size] for i in range(n_spots)], dtype=float)

    return {
        'counts': counts.toarray(),
        'genes': genes,
        'proportions': props.values,
        'cell_types': list(props.columns),
        'coords': coords
    }


def build_signature_matrix(ref, target_genes, target_celltypes):
    """Build cell type signature matrix."""
    ref_genes = list(ref['genes'])
    common_genes = [g for g in target_genes if g in ref_genes]

    if len(common_genes) < 50:
        return None, None

    target_idx = [list(target_genes).index(g) for g in common_genes]
    ref_idx = [ref_genes.index(g) for g in common_genes]

    ref_counts = ref['counts'][:, ref_idx]
    if hasattr(ref_counts, 'toarray'):
        ref_counts = ref_counts.toarray()

    ref_cts = ref['celltypes']
    unique_ref_cts = list(set(ref_cts))

    X_ref = np.zeros((len(target_celltypes), len(common_genes)))
    for i, ct in enumerate(target_celltypes):
        if ct in unique_ref_cts:
            mask = np.array(ref_cts) == ct
        else:
            # Partial match
            matches = [rct for rct in unique_ref_cts if ct in rct or rct in ct]
            if matches:
                mask = np.array(ref_cts) == matches[0]
            else:
                mask = np.ones(len(ref_cts), dtype=bool)

        if mask.sum() > 0:
            X_ref[i] = ref_counts[mask].mean(axis=0)
        else:
            X_ref[i] = ref_counts.mean(axis=0)

    return X_ref, target_idx


def calculate_metrics(pred, true):
    """Calculate RMSE, JSD, AUPR, Pearson."""
    # RMSE
    rmse = np.sqrt(np.mean((pred - true) ** 2))

    # JSD
    jsd_values = []
    for i in range(pred.shape[0]):
        p = pred[i] + 1e-10
        q = true[i] + 1e-10
        p = p / p.sum()
        q = q / q.sum()
        jsd_values.append(jensenshannon(p, q) ** 2)
    jsd = np.mean(jsd_values)

    # AUPR
    true_binary = (true.flatten() > 0.01).astype(int)
    pred_flat = pred.flatten()
    if true_binary.sum() > 0 and true_binary.sum() < len(true_binary):
        precision, recall, _ = precision_recall_curve(true_binary, pred_flat)
        aupr = auc(recall, precision)
    else:
        aupr = np.nan

    # Pearson
    pearson = pearsonr(pred.flatten(), true.flatten())[0]

    return rmse, jsd, aupr, pearson


def run_benchmark():
    """Run full benchmark on all datasets."""
    base_dir = "/Users/apple/Research/FlashDeconv/validation/benchmark_data/converted"

    # Dataset info
    datasets = {
        "1": {"name": "brain_cortex", "patterns": 11},
        "2": {"name": "cerebellum_cell", "patterns": 9},
        "3": {"name": "cerebellum_nucleus", "patterns": 9},
        "4": {"name": "hippocampus", "patterns": 9},
        "5": {"name": "kidney", "patterns": 9},
        "6": {"name": "scc_p5", "patterns": 9},
    }

    all_results = []

    for ds_id, ds_info in datasets.items():
        ds_name = ds_info["name"]
        print(f"\n{'='*60}")
        print(f"Dataset {ds_id}: {ds_name}")
        print('='*60)

        # Load reference
        ref_path = f"{base_dir}/reference_{ds_id}"
        if not os.path.exists(f"{ref_path}_counts.mtx"):
            print(f"  Reference not found, skipping")
            continue

        ref = load_reference(ref_path)
        print(f"  Reference: {ref['counts'].shape[0]} cells, {len(set(ref['celltypes']))} cell types")

        # Test each pattern
        for pattern_id in range(1, ds_info["patterns"] + 1):
            test_path = f"{base_dir}/silver_{ds_id}_{pattern_id}"
            if not os.path.exists(f"{test_path}_counts.mtx"):
                continue

            try:
                test = load_test_data(test_path)
                print(f"\n  Pattern {pattern_id}: {test['counts'].shape[0]} spots, {len(test['cell_types'])} cell types")

                # Build signature matrix
                X_ref, gene_idx = build_signature_matrix(ref, test['genes'], test['cell_types'])
                if X_ref is None:
                    print(f"    Insufficient common genes, skipping")
                    continue

                Y = test['counts'][:, gene_idx]
                ground_truth = test['proportions']
                coords = test['coords']

                print(f"    Common genes: {len(gene_idx)}")

                # Run FlashDeconv with log_cpm
                t0 = time.time()
                model = FlashDeconv(
                    sketch_dim=min(512, Y.shape[1]),
                    lambda_spatial="auto",
                    rho_sparsity=0.01,
                    preprocess="log_cpm",
                    n_hvg=min(2000, Y.shape[1]),
                    max_iter=100,
                    verbose=False,
                    random_state=42,
                )
                pred = model.fit_transform(Y, X_ref, coords)
                elapsed = time.time() - t0

                rmse, jsd, aupr, pearson = calculate_metrics(pred, ground_truth)
                print(f"    log_cpm: Pearson={pearson:.4f}, RMSE={rmse:.4f}, JSD={jsd:.4f}, Time={elapsed:.1f}s")

                all_results.append({
                    'dataset_id': ds_id,
                    'dataset_name': ds_name,
                    'pattern': pattern_id,
                    'preprocess': 'log_cpm',
                    'pearson': pearson,
                    'rmse': rmse,
                    'jsd': jsd,
                    'aupr': aupr,
                    'time': elapsed,
                    'n_spots': Y.shape[0],
                    'n_genes': len(gene_idx),
                    'n_celltypes': len(test['cell_types']),
                })

            except Exception as e:
                print(f"    Error: {e}")

    return pd.DataFrame(all_results)


def main():
    print("="*60)
    print("Fair Benchmark: FlashDeconv on ALL Spotless Silver Standards")
    print("="*60)

    df = run_benchmark()

    if len(df) == 0:
        print("No results!")
        return

    # Save results
    df.to_csv('/Users/apple/Research/FlashDeconv/validation/fair_silver_benchmark_results.csv', index=False)

    # Summary
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)

    # Per-dataset results
    print("\nPer-dataset results:")
    for ds_id in sorted(df['dataset_id'].unique()):
        ds_df = df[df['dataset_id'] == ds_id]
        ds_name = ds_df['dataset_name'].iloc[0]
        print(f"  {ds_id}. {ds_name:20s}: Pearson={ds_df['pearson'].mean():.4f}±{ds_df['pearson'].std():.4f}, "
              f"RMSE={ds_df['rmse'].mean():.4f}, n={len(ds_df)}")

    # Overall average
    print(f"\n{'='*60}")
    print("OVERALL AVERAGE (comparable to Spotless paper)")
    print('='*60)
    print(f"  Pearson: {df['pearson'].mean():.4f} ± {df['pearson'].std():.4f}")
    print(f"  RMSE:    {df['rmse'].mean():.4f} ± {df['rmse'].std():.4f}")
    print(f"  JSD:     {df['jsd'].mean():.4f} ± {df['jsd'].std():.4f}")
    print(f"  AUPR:    {df['aupr'].mean():.4f} ± {df['aupr'].std():.4f}")
    print(f"  Total:   {len(df)} dataset-pattern combinations")

    # Comparison with published results
    print(f"\n{'='*60}")
    print("COMPARISON WITH SPOTLESS PAPER")
    print('='*60)
    print("Spotless published results (63 silver standards average):")
    print("  Method          Pearson   RMSE")
    print("  rctd            0.9046    0.0613")
    print("  cell2location   0.8953    0.0603")
    print("  music           0.8901    0.0775")
    print("  spatialdwls     0.8751    0.0654")
    print("  nnls            0.8133    0.0868")
    print(f"\n  FlashDeconv     {df['pearson'].mean():.4f}    {df['rmse'].mean():.4f}  ← Our result")

    # Rank
    our_pearson = df['pearson'].mean()
    published = [
        ('rctd', 0.9046),
        ('cell2location', 0.8953),
        ('music', 0.8901),
        ('spatialdwls', 0.8751),
        ('nnls', 0.8133),
    ]
    all_methods = published + [('FlashDeconv', our_pearson)]
    all_methods.sort(key=lambda x: x[1], reverse=True)
    rank = [i+1 for i, (m, _) in enumerate(all_methods) if m == 'FlashDeconv'][0]
    print(f"\n  FlashDeconv rank: #{rank} out of {len(all_methods)} methods")


if __name__ == '__main__':
    main()
