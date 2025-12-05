"""
Run FlashDeconv on Spotless Silver Standards with proper evaluation metrics.
Metrics: RMSE, JSD, AUPR (as used in the Spotless paper)
"""

import os
import sys
import time
import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.sparse import csr_matrix
from scipy.stats import pearsonr
from scipy.spatial.distance import jensenshannon
from scipy.optimize import nnls
from sklearn.metrics import precision_recall_curve, auc
import warnings
warnings.filterwarnings('ignore')

sys.path.insert(0, '/Users/apple/Research/FlashDeconv')
from flashdeconv import FlashDeconv


def load_dataset(prefix):
    """Load a converted Spotless dataset."""
    counts = mmread(f"{prefix}_counts.mtx").T.tocsr()  # spots x genes
    genes = pd.read_csv(f"{prefix}_genes.txt", header=None)[0].values
    spots = pd.read_csv(f"{prefix}_spots.txt", header=None)[0].values
    proportions = pd.read_csv(f"{prefix}_proportions.csv", index_col=0)

    # Filter to numeric columns only
    numeric_cols = proportions.select_dtypes(include=[np.number]).columns
    proportions = proportions[numeric_cols]

    # Generate grid coordinates
    n_spots = counts.shape[0]
    grid_size = int(np.ceil(np.sqrt(n_spots)))
    coords = np.array([[i % grid_size, i // grid_size] for i in range(n_spots)], dtype=float)

    return {
        'counts': counts.toarray(),
        'genes': genes,
        'spots': spots,
        'proportions': proportions.values,
        'cell_types': list(proportions.columns),
        'coords': coords
    }


def load_reference(prefix):
    """Load a converted reference dataset."""
    counts = mmread(f"{prefix}_counts.mtx").T.tocsr()  # cells x genes
    genes = pd.read_csv(f"{prefix}_genes.txt", header=None)[0].values
    cells = pd.read_csv(f"{prefix}_cells.txt", header=None)[0].values
    cell_types = pd.read_csv(f"{prefix}_celltypes.txt", header=None)[0].values

    return {
        'counts': counts,
        'genes': genes,
        'cells': cells,
        'cell_types': cell_types
    }


def build_reference_matrix(ref_data, target_genes, target_cell_types):
    """Build reference signature matrix from scRNA-seq data."""
    # Find common genes
    ref_genes = list(ref_data['genes'])
    common_genes = [g for g in target_genes if g in ref_genes]

    if len(common_genes) < 100:
        raise ValueError(f"Too few common genes: {len(common_genes)}")

    # Get indices
    target_idx = [list(target_genes).index(g) for g in common_genes]
    ref_idx = [ref_genes.index(g) for g in common_genes]

    # Get unique cell types in reference
    ref_cts = ref_data['cell_types']
    unique_ref_cts = list(set(ref_cts))

    # Build signature matrix
    ref_counts = ref_data['counts'][:, ref_idx]
    if hasattr(ref_counts, 'toarray'):
        ref_counts = ref_counts.toarray()

    X_ref = np.zeros((len(target_cell_types), len(common_genes)))
    for i, ct in enumerate(target_cell_types):
        # Find matching cells
        if ct in unique_ref_cts:
            mask = np.array(ref_cts) == ct
        else:
            # Try partial match
            matches = [rct for rct in unique_ref_cts if ct in rct or rct in ct]
            if matches:
                mask = np.array(ref_cts) == matches[0]
            else:
                # Use mean of all cells
                mask = np.ones(len(ref_cts), dtype=bool)

        if mask.sum() > 0:
            X_ref[i] = ref_counts[mask].mean(axis=0)
        else:
            X_ref[i] = ref_counts.mean(axis=0)

    return X_ref, target_idx


def calculate_rmse(pred, true):
    """Root Mean Squared Error."""
    return np.sqrt(np.mean((pred - true) ** 2))


def calculate_jsd(pred, true):
    """Jensen-Shannon Divergence (averaged over spots)."""
    jsd_values = []
    for i in range(pred.shape[0]):
        p = pred[i] + 1e-10
        q = true[i] + 1e-10
        p = p / p.sum()
        q = q / q.sum()
        jsd_values.append(jensenshannon(p, q) ** 2)
    return np.mean(jsd_values)


def calculate_aupr(pred, true, threshold=0.01):
    """Area Under Precision-Recall Curve."""
    true_binary = (true.flatten() > threshold).astype(int)
    pred_flat = pred.flatten()

    if true_binary.sum() == 0 or true_binary.sum() == len(true_binary):
        return np.nan

    precision, recall, _ = precision_recall_curve(true_binary, pred_flat)
    return auc(recall, precision)


def run_nnls(Y, X_ref):
    """Simple NNLS deconvolution."""
    Y_norm = Y / (Y.sum(axis=1, keepdims=True) + 1e-10) * 1e4
    X_norm = X_ref / (X_ref.sum(axis=1, keepdims=True) + 1e-10) * 1e4

    n_spots = Y.shape[0]
    n_types = X_ref.shape[0]
    beta = np.zeros((n_spots, n_types))

    for i in range(n_spots):
        beta[i], _ = nnls(X_norm.T, Y_norm[i])

    # Normalize to proportions
    row_sums = beta.sum(axis=1, keepdims=True)
    beta = beta / (row_sums + 1e-10)
    return beta


def benchmark_dataset(data, ref_data, name):
    """Benchmark FlashDeconv and NNLS on a single dataset."""
    print(f"\n{'='*60}")
    print(f"Dataset: {name}")
    print(f"{'='*60}")

    Y = data['counts']
    ground_truth = data['proportions']
    cell_types = data['cell_types']
    coords = data['coords']

    print(f"Spots: {Y.shape[0]}, Genes: {Y.shape[1]}, Cell types: {len(cell_types)}")

    # Build reference matrix
    X_ref, gene_idx = build_reference_matrix(ref_data, data['genes'], cell_types)
    Y_common = Y[:, gene_idx]
    print(f"Common genes: {len(gene_idx)}, Reference shape: {X_ref.shape}")

    results = {}

    # Run FlashDeconv with preprocess=True (Log-CPM)
    print("\n--- FlashDeconv (preprocess=True, Log-CPM) ---")
    t0 = time.time()
    try:
        model = FlashDeconv(
            sketch_dim=512,
            lambda_spatial="auto",
            rho_sparsity=0.01,
            preprocess=True,  # Log-CPM normalization
            n_hvg=min(2000, Y_common.shape[1]),
            max_iter=100,
            verbose=False,
            random_state=42,
        )
        pred_flash = model.fit_transform(Y_common, X_ref, coords)
        flash_time = time.time() - t0

        rmse = calculate_rmse(pred_flash, ground_truth)
        jsd = calculate_jsd(pred_flash, ground_truth)
        aupr = calculate_aupr(pred_flash, ground_truth)
        pearson = pearsonr(pred_flash.flatten(), ground_truth.flatten())[0]

        print(f"Time: {flash_time:.2f}s")
        print(f"RMSE: {rmse:.4f}, JSD: {jsd:.4f}, AUPR: {aupr:.4f}, Pearson: {pearson:.4f}")

        results['flashdeconv_logcpm'] = {
            'time': flash_time, 'rmse': rmse, 'jsd': jsd, 'aupr': aupr, 'pearson': pearson
        }
    except Exception as e:
        print(f"Error: {e}")
        results['flashdeconv_logcpm'] = None

    # Run FlashDeconv with preprocess=False (raw)
    print("\n--- FlashDeconv (preprocess=False, raw) ---")
    t0 = time.time()
    try:
        model = FlashDeconv(
            sketch_dim=512,
            lambda_spatial="auto",
            rho_sparsity=0.01,
            preprocess=False,
            n_hvg=min(2000, Y_common.shape[1]),
            max_iter=100,
            verbose=False,
            random_state=42,
        )
        pred_flash = model.fit_transform(Y_common, X_ref, coords)
        flash_time = time.time() - t0

        rmse = calculate_rmse(pred_flash, ground_truth)
        jsd = calculate_jsd(pred_flash, ground_truth)
        aupr = calculate_aupr(pred_flash, ground_truth)
        pearson = pearsonr(pred_flash.flatten(), ground_truth.flatten())[0]

        print(f"Time: {flash_time:.2f}s")
        print(f"RMSE: {rmse:.4f}, JSD: {jsd:.4f}, AUPR: {aupr:.4f}, Pearson: {pearson:.4f}")

        results['flashdeconv_raw'] = {
            'time': flash_time, 'rmse': rmse, 'jsd': jsd, 'aupr': aupr, 'pearson': pearson
        }
    except Exception as e:
        print(f"Error: {e}")
        results['flashdeconv_raw'] = None

    # Run NNLS
    print("\n--- NNLS ---")
    t0 = time.time()
    try:
        pred_nnls = run_nnls(Y_common, X_ref)
        nnls_time = time.time() - t0

        rmse = calculate_rmse(pred_nnls, ground_truth)
        jsd = calculate_jsd(pred_nnls, ground_truth)
        aupr = calculate_aupr(pred_nnls, ground_truth)
        pearson = pearsonr(pred_nnls.flatten(), ground_truth.flatten())[0]

        print(f"Time: {nnls_time:.2f}s")
        print(f"RMSE: {rmse:.4f}, JSD: {jsd:.4f}, AUPR: {aupr:.4f}, Pearson: {pearson:.4f}")

        results['nnls'] = {
            'time': nnls_time, 'rmse': rmse, 'jsd': jsd, 'aupr': aupr, 'pearson': pearson
        }
    except Exception as e:
        print(f"Error: {e}")
        results['nnls'] = None

    return results


def main():
    print("="*60)
    print("Spotless Silver Standard Benchmark")
    print("="*60)

    base_dir = "/Users/apple/Research/FlashDeconv/validation/benchmark_data/converted"

    # Load reference
    print("\nLoading reference...")
    ref_data = load_reference(f"{base_dir}/reference_brain_cortex")
    print(f"Reference: {ref_data['counts'].shape[0]} cells, {ref_data['counts'].shape[1]} genes")
    print(f"Cell types: {len(set(ref_data['cell_types']))}")

    # Find all silver standard datasets
    datasets = []
    for i in range(1, 12):
        prefix = f"{base_dir}/silver_standard_1-{i}"
        if os.path.exists(f"{prefix}_counts.mtx"):
            datasets.append((f"silver_1-{i}", prefix))

    print(f"\nFound {len(datasets)} datasets")

    all_results = []

    for name, prefix in datasets:
        try:
            data = load_dataset(prefix)
            results = benchmark_dataset(data, ref_data, name)

            row = {'dataset': name}
            if results.get('flashdeconv_logcpm'):
                row.update({
                    'logcpm_rmse': results['flashdeconv_logcpm']['rmse'],
                    'logcpm_jsd': results['flashdeconv_logcpm']['jsd'],
                    'logcpm_aupr': results['flashdeconv_logcpm']['aupr'],
                    'logcpm_pearson': results['flashdeconv_logcpm']['pearson'],
                    'logcpm_time': results['flashdeconv_logcpm']['time'],
                })
            if results.get('flashdeconv_raw'):
                row.update({
                    'raw_rmse': results['flashdeconv_raw']['rmse'],
                    'raw_jsd': results['flashdeconv_raw']['jsd'],
                    'raw_aupr': results['flashdeconv_raw']['aupr'],
                    'raw_pearson': results['flashdeconv_raw']['pearson'],
                    'raw_time': results['flashdeconv_raw']['time'],
                })
            if results.get('nnls'):
                row.update({
                    'nnls_rmse': results['nnls']['rmse'],
                    'nnls_jsd': results['nnls']['jsd'],
                    'nnls_aupr': results['nnls']['aupr'],
                    'nnls_pearson': results['nnls']['pearson'],
                    'nnls_time': results['nnls']['time'],
                })
            all_results.append(row)
        except Exception as e:
            print(f"Error on {name}: {e}")

    # Summary
    if all_results:
        print("\n" + "="*60)
        print("SUMMARY")
        print("="*60)

        df = pd.DataFrame(all_results)

        if 'logcpm_rmse' in df.columns:
            print("\nFlashDeconv (Log-CPM, preprocess=True):")
            print(f"  RMSE:    {df['logcpm_rmse'].mean():.4f} ± {df['logcpm_rmse'].std():.4f}")
            print(f"  JSD:     {df['logcpm_jsd'].mean():.4f} ± {df['logcpm_jsd'].std():.4f}")
            print(f"  AUPR:    {df['logcpm_aupr'].mean():.4f} ± {df['logcpm_aupr'].std():.4f}")
            print(f"  Pearson: {df['logcpm_pearson'].mean():.4f} ± {df['logcpm_pearson'].std():.4f}")
            print(f"  Time:    {df['logcpm_time'].mean():.2f}s")

        if 'raw_rmse' in df.columns:
            print("\nFlashDeconv (Raw, preprocess=False):")
            print(f"  RMSE:    {df['raw_rmse'].mean():.4f} ± {df['raw_rmse'].std():.4f}")
            print(f"  JSD:     {df['raw_jsd'].mean():.4f} ± {df['raw_jsd'].std():.4f}")
            print(f"  AUPR:    {df['raw_aupr'].mean():.4f} ± {df['raw_aupr'].std():.4f}")
            print(f"  Pearson: {df['raw_pearson'].mean():.4f} ± {df['raw_pearson'].std():.4f}")
            print(f"  Time:    {df['raw_time'].mean():.2f}s")

        if 'nnls_rmse' in df.columns:
            print("\nNNLS (baseline):")
            print(f"  RMSE:    {df['nnls_rmse'].mean():.4f} ± {df['nnls_rmse'].std():.4f}")
            print(f"  JSD:     {df['nnls_jsd'].mean():.4f} ± {df['nnls_jsd'].std():.4f}")
            print(f"  AUPR:    {df['nnls_aupr'].mean():.4f} ± {df['nnls_aupr'].std():.4f}")
            print(f"  Pearson: {df['nnls_pearson'].mean():.4f} ± {df['nnls_pearson'].std():.4f}")
            print(f"  Time:    {df['nnls_time'].mean():.2f}s")

        if 'logcpm_rmse' in df.columns and 'nnls_rmse' in df.columns:
            print("\nComparison (Log-CPM vs NNLS):")
            print(f"  RMSE:    {'Better ↓' if df['logcpm_rmse'].mean() < df['nnls_rmse'].mean() else 'Worse ↑'}")
            print(f"  Pearson: {'Better ↑' if df['logcpm_pearson'].mean() > df['nnls_pearson'].mean() else 'Worse ↓'}")

        # Save results
        df.to_csv('/Users/apple/Research/FlashDeconv/validation/spotless_benchmark_results.csv', index=False)
        print(f"\nResults saved to spotless_benchmark_results.csv")

        # Reference
        print("\n" + "="*60)
        print("Reference: Spotless Paper (Sang-Aram et al., eLife 2024)")
        print("="*60)
        print("Top methods RMSE (lower is better):")
        print("  Cell2location: ~0.05-0.08")
        print("  RCTD:          ~0.06-0.09")
        print("  DestVI:        ~0.07-0.10")
        print("  Spotlight:     ~0.08-0.12")


if __name__ == "__main__":
    main()
