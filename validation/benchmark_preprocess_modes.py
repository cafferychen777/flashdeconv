"""
Comprehensive benchmark of preprocess modes on real data.

Tests:
1. Mouse Brain Visium + scRNA-seq reference (full reference)
2. Mouse Brain Visium + Allen Cortex reference (cross-platform)
"""

import numpy as np
import pandas as pd
import scanpy as sc
import sys
import time
from scipy.stats import pearsonr
sys.path.insert(0, '/Users/apple/Research/FlashDeconv')

from flashdeconv import FlashDeconv


def load_mouse_brain_with_scrna_ref():
    """Load Mouse Brain Visium with paired scRNA-seq reference."""
    print("Loading Mouse Brain + scRNA-seq reference...")

    sp_path = '/Users/apple/Research/SpatialTrans_MCP/chatspatial/data/processed_datasets/mouse_brain/mouse_brain_sparse_slice1.h5ad'
    sc_path = '/Users/apple/Research/SpatialTrans_MCP/chatspatial/data/processed_datasets/mouse_brain/mousebrain_reference_clean.h5ad'

    adata_sp = sc.read_h5ad(sp_path)
    adata_sc = sc.read_h5ad(sc_path)
    adata_sp.var_names_make_unique()
    adata_sc.var_names_make_unique()

    # Common genes
    common_genes = list(set(adata_sp.var_names) & set(adata_sc.var_names))
    adata_sp = adata_sp[:, common_genes]
    adata_sc = adata_sc[:, common_genes]

    Y = adata_sp.X.toarray() if hasattr(adata_sp.X, 'toarray') else np.array(adata_sp.X)
    coords = adata_sp.obsm['spatial']

    # Build reference
    cell_types = list(adata_sc.obs['cell_type_clean'].unique())
    X = np.zeros((len(cell_types), len(common_genes)))
    for i, ct in enumerate(cell_types):
        mask = (adata_sc.obs['cell_type_clean'] == ct).values
        ct_data = adata_sc.X[mask]
        if hasattr(ct_data, 'toarray'):
            ct_data = ct_data.toarray()
        X[i] = np.mean(ct_data, axis=0)

    print(f"  Spatial: {Y.shape}, Reference: {X.shape}, Cell types: {len(cell_types)}")
    return Y, X, coords, cell_types, "MouseBrain_scRNA"


def load_mouse_brain_with_allen_ref():
    """Load Mouse Brain Visium with Allen Cortex reference."""
    print("Loading Mouse Brain + Allen Cortex reference...")

    sp_path = '/Users/apple/Research/SpatialTrans_MCP/chatspatial/data/processed_datasets/mouse_brain/mouse_brain_sparse_slice1.h5ad'
    sc_path = '/Users/apple/Research/FlashDeconv/validation/benchmark_data/allen_cortex.h5ad'

    adata_sp = sc.read_h5ad(sp_path)
    adata_sc = sc.read_h5ad(sc_path)
    adata_sp.var_names_make_unique()
    adata_sc.var_names_make_unique()

    # Common genes
    common_genes = list(set(adata_sp.var_names) & set(adata_sc.var_names))
    adata_sp = adata_sp[:, common_genes]
    adata_sc = adata_sc[:, common_genes]

    Y = adata_sp.X.toarray() if hasattr(adata_sp.X, 'toarray') else np.array(adata_sp.X)
    coords = adata_sp.obsm['spatial']

    # Build reference from Allen Cortex (use subclass)
    cell_types = list(adata_sc.obs['subclass'].unique())
    X = np.zeros((len(cell_types), len(common_genes)))
    for i, ct in enumerate(cell_types):
        mask = (adata_sc.obs['subclass'] == ct).values
        ct_data = adata_sc.X[mask]
        if hasattr(ct_data, 'toarray'):
            ct_data = ct_data.toarray()
        X[i] = np.mean(ct_data, axis=0)

    print(f"  Spatial: {Y.shape}, Reference: {X.shape}, Cell types: {len(cell_types)}")
    return Y, X, coords, cell_types, "MouseBrain_Allen"


def evaluate_spatial_coherence(proportions, coords):
    """Compute spatial coherence (neighbor similarity)."""
    from flashdeconv.utils.graph import coords_to_adjacency
    A = coords_to_adjacency(coords, method='knn', k=6)

    neighbor_corrs = []
    for i in range(A.shape[0]):
        neighbors = A[i].nonzero()[1]
        for j in neighbors:
            r = np.corrcoef(proportions[i], proportions[j])[0, 1]
            if not np.isnan(r):
                neighbor_corrs.append(r)

    return np.mean(neighbor_corrs)


def run_benchmark(Y, X, coords, cell_types, dataset_name):
    """Run benchmark for all preprocess modes."""
    print(f"\n{'='*60}")
    print(f"Dataset: {dataset_name}")
    print(f"{'='*60}")

    results = []

    for preprocess in ["log_cpm", "pearson", "raw"]:
        print(f"\n  Testing preprocess='{preprocess}'...")
        t0 = time.time()

        model = FlashDeconv(
            sketch_dim=512,
            lambda_spatial="auto",
            rho_sparsity=0.01,
            preprocess=preprocess,
            n_hvg=2000,
            max_iter=100,
            verbose=False,
            random_state=42,
        )
        proportions = model.fit_transform(Y, X, coords)
        elapsed = time.time() - t0

        # Metrics
        spatial_coherence = evaluate_spatial_coherence(proportions, coords)
        nonzero_per_spot = (proportions > 0.01).sum(axis=1).mean()

        # Top cell types
        est_comp = proportions.mean(axis=0)
        sorted_idx = np.argsort(est_comp)[::-1]
        top_types = [(cell_types[i], est_comp[i]) for i in sorted_idx[:3]]

        results.append({
            'dataset': dataset_name,
            'preprocess': preprocess,
            'time': elapsed,
            'spatial_coherence': spatial_coherence,
            'types_per_spot': nonzero_per_spot,
            'top_type': top_types[0][0],
            'top_prop': top_types[0][1],
        })

        print(f"    Time: {elapsed:.2f}s, Spatial Coherence: {spatial_coherence:.4f}")
        print(f"    Types/spot: {nonzero_per_spot:.1f}, Top: {top_types[0][0]} ({top_types[0][1]:.3f})")

    return results


def main():
    print("="*60)
    print("FlashDeconv Preprocess Mode Benchmark")
    print("="*60)

    all_results = []

    # Test 1: Mouse Brain + scRNA-seq reference
    try:
        Y, X, coords, cell_types, name = load_mouse_brain_with_scrna_ref()
        results = run_benchmark(Y, X, coords, cell_types, name)
        all_results.extend(results)
    except Exception as e:
        print(f"Error: {e}")

    # Test 2: Mouse Brain + Allen Cortex reference
    try:
        Y, X, coords, cell_types, name = load_mouse_brain_with_allen_ref()
        results = run_benchmark(Y, X, coords, cell_types, name)
        all_results.extend(results)
    except Exception as e:
        print(f"Error: {e}")

    # Summary
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)

    df = pd.DataFrame(all_results)
    print("\n" + df.to_string(index=False))

    # Best by dataset
    print("\n" + "-"*60)
    print("Best preprocess by dataset (highest spatial coherence):")
    for dataset in df['dataset'].unique():
        sub = df[df['dataset'] == dataset]
        best = sub.loc[sub['spatial_coherence'].idxmax()]
        print(f"  {dataset}: {best['preprocess']} (coherence={best['spatial_coherence']:.4f})")

    # Save
    df.to_csv('/Users/apple/Research/FlashDeconv/validation/preprocess_benchmark_results.csv', index=False)
    print("\nResults saved to preprocess_benchmark_results.csv")

    return df


if __name__ == "__main__":
    main()
