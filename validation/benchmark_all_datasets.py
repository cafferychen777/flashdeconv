"""
Comprehensive benchmark on all available datasets.

Datasets:
1. Spotless (synthetic, has ground truth)
2. CARD data (real, from R export)
3. DestVI data (real, has spatial)
4. Mouse Brain + scRNA-seq reference
"""

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.io import mmread
from scipy.stats import pearsonr
import sys
import time
sys.path.insert(0, '/Users/apple/Research/FlashDeconv')

from flashdeconv import FlashDeconv
from flashdeconv.utils.graph import coords_to_adjacency


def evaluate_spatial_coherence(proportions, coords):
    """Compute spatial coherence (neighbor similarity)."""
    A = coords_to_adjacency(coords, method='knn', k=6)
    neighbor_corrs = []
    for i in range(A.shape[0]):
        neighbors = A[i].nonzero()[1]
        for j in neighbors:
            r = np.corrcoef(proportions[i], proportions[j])[0, 1]
            if not np.isnan(r):
                neighbor_corrs.append(r)
    return np.mean(neighbor_corrs) if neighbor_corrs else 0.0


def test_preprocess_modes(Y, X, coords, ground_truth=None, dataset_name=""):
    """Test all preprocess modes on a dataset."""
    results = []

    for preprocess in ["log_cpm", "pearson", "raw"]:
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
        pred = model.fit_transform(Y, X, coords)
        elapsed = time.time() - t0

        # Metrics
        spatial_coh = evaluate_spatial_coherence(pred, coords)
        types_per_spot = (pred > 0.01).sum(axis=1).mean()

        res = {
            'dataset': dataset_name,
            'preprocess': preprocess,
            'time': elapsed,
            'spatial_coherence': spatial_coh,
            'types_per_spot': types_per_spot,
        }

        # If ground truth available
        if ground_truth is not None:
            pearson_r = pearsonr(pred.flatten(), ground_truth.flatten())[0]
            rmse = np.sqrt(np.mean((pred - ground_truth)**2))
            res['pearson'] = pearson_r
            res['rmse'] = rmse

        results.append(res)

        if ground_truth is not None:
            print(f"  {preprocess:<10} Pearson={res['pearson']:.4f}, RMSE={res['rmse']:.4f}, Time={elapsed:.2f}s")
        else:
            print(f"  {preprocess:<10} SpatialCoh={spatial_coh:.4f}, Types/spot={types_per_spot:.1f}, Time={elapsed:.2f}s")

    return results


def load_spotless():
    """Load Spotless dataset with ground truth."""
    print("\n" + "="*60)
    print("Dataset: Spotless (synthetic, with ground truth)")
    print("="*60)

    base = '/Users/apple/Research/FlashDeconv/validation/benchmark_data/converted'
    counts = mmread(f'{base}/silver_standard_1-1_counts.mtx').T.tocsr().toarray()
    genes = pd.read_csv(f'{base}/silver_standard_1-1_genes.txt', header=None)[0].values
    props = pd.read_csv(f'{base}/silver_standard_1-1_proportions.csv', index_col=0)
    props = props.select_dtypes(include=[np.number])

    ref_counts = mmread(f'{base}/reference_brain_cortex_counts.mtx').T.tocsr()
    ref_genes = pd.read_csv(f'{base}/reference_brain_cortex_genes.txt', header=None)[0].values
    ref_cts = pd.read_csv(f'{base}/reference_brain_cortex_celltypes.txt', header=None)[0].values

    common = [g for g in genes if g in ref_genes]
    Y = counts[:, [list(genes).index(g) for g in common]]
    ref_sub = ref_counts[:, [list(ref_genes).index(g) for g in common]].toarray()

    X = np.zeros((len(props.columns), len(common)))
    for i, ct in enumerate(props.columns):
        if ct in set(ref_cts):
            X[i] = ref_sub[np.array(ref_cts) == ct].mean(axis=0)
        else:
            X[i] = ref_sub.mean(axis=0)

    gt = props.values
    n = Y.shape[0]
    coords = np.array([[i % 29, i // 29] for i in range(n)], dtype=float)

    print(f"  Y: {Y.shape}, X: {X.shape}, GT: {gt.shape}")
    return Y, X, coords, gt, "Spotless"


def load_card():
    """Load CARD data."""
    print("\n" + "="*60)
    print("Dataset: CARD (real data)")
    print("="*60)

    try:
        spatial_count = pd.read_csv("/Users/apple/Research/FlashDeconv/validation/spatial_count.csv", index_col=0)
        sc_count = pd.read_csv("/Users/apple/Research/FlashDeconv/validation/sc_count.csv", index_col=0)
        sc_meta = pd.read_csv("/Users/apple/Research/FlashDeconv/validation/sc_meta.csv", index_col=0)
        spatial_location = pd.read_csv("/Users/apple/Research/FlashDeconv/validation/card_spatial_location.csv", index_col=0)
    except FileNotFoundError as e:
        print(f"  CARD data not found: {e}")
        return None

    # Common genes
    common = list(set(spatial_count.index) & set(sc_count.index))
    Y = spatial_count.loc[common].T.values
    ref = sc_count.loc[common].T.values

    ct_col = "cellType" if "cellType" in sc_meta.columns else "cell_type"
    cell_types = sc_meta[ct_col].values
    unique_types = np.unique(cell_types)

    X = np.zeros((len(unique_types), len(common)))
    for i, ct in enumerate(unique_types):
        X[i] = ref[cell_types == ct].mean(axis=0)

    coords = spatial_location[['x', 'y']].values

    print(f"  Y: {Y.shape}, X: {X.shape}, coords: {coords.shape}")
    return Y, X, coords, None, "CARD"


def load_destvi():
    """Load DestVI data."""
    print("\n" + "="*60)
    print("Dataset: DestVI (real data)")
    print("="*60)

    try:
        adata_sp = sc.read_h5ad('/Users/apple/Research/FlashDeconv/validation/benchmark_data/destvi_spatial.h5ad')
        adata_sc = sc.read_h5ad('/Users/apple/Research/FlashDeconv/validation/benchmark_data/destvi_scrna.h5ad')
    except FileNotFoundError as e:
        print(f"  DestVI data not found: {e}")
        return None

    adata_sp.var_names_make_unique()
    adata_sc.var_names_make_unique()

    common = list(set(adata_sp.var_names) & set(adata_sc.var_names))
    adata_sp = adata_sp[:, common]
    adata_sc = adata_sc[:, common]

    Y = adata_sp.X.toarray() if hasattr(adata_sp.X, 'toarray') else np.array(adata_sp.X)
    coords = adata_sp.obsm['spatial']

    cell_types = list(adata_sc.obs['cell_types'].unique())
    X = np.zeros((len(cell_types), len(common)))
    for i, ct in enumerate(cell_types):
        mask = (adata_sc.obs['cell_types'] == ct).values
        ct_data = adata_sc.X[mask]
        if hasattr(ct_data, 'toarray'):
            ct_data = ct_data.toarray()
        X[i] = ct_data.mean(axis=0)

    print(f"  Y: {Y.shape}, X: {X.shape}, coords: {coords.shape}, Cell types: {len(cell_types)}")
    return Y, X, coords, None, "DestVI"


def load_mouse_brain():
    """Load Mouse Brain + scRNA-seq reference."""
    print("\n" + "="*60)
    print("Dataset: Mouse Brain Visium + scRNA-seq")
    print("="*60)

    try:
        sp_path = '/Users/apple/Research/SpatialTrans_MCP/chatspatial/data/processed_datasets/mouse_brain/mouse_brain_sparse_slice1.h5ad'
        sc_path = '/Users/apple/Research/SpatialTrans_MCP/chatspatial/data/processed_datasets/mouse_brain/mousebrain_reference_clean.h5ad'

        adata_sp = sc.read_h5ad(sp_path)
        adata_sc = sc.read_h5ad(sc_path)
    except FileNotFoundError as e:
        print(f"  Mouse Brain data not found: {e}")
        return None

    adata_sp.var_names_make_unique()
    adata_sc.var_names_make_unique()

    common = list(set(adata_sp.var_names) & set(adata_sc.var_names))
    adata_sp = adata_sp[:, common]
    adata_sc = adata_sc[:, common]

    Y = adata_sp.X.toarray() if hasattr(adata_sp.X, 'toarray') else np.array(adata_sp.X)
    coords = adata_sp.obsm['spatial']

    cell_types = list(adata_sc.obs['cell_type_clean'].unique())
    X = np.zeros((len(cell_types), len(common)))
    for i, ct in enumerate(cell_types):
        mask = (adata_sc.obs['cell_type_clean'] == ct).values
        ct_data = adata_sc.X[mask]
        if hasattr(ct_data, 'toarray'):
            ct_data = ct_data.toarray()
        X[i] = ct_data.mean(axis=0)

    print(f"  Y: {Y.shape}, X: {X.shape}, coords: {coords.shape}, Cell types: {len(cell_types)}")
    return Y, X, coords, None, "MouseBrain"


def main():
    print("="*60)
    print("FlashDeconv Comprehensive Benchmark")
    print("="*60)

    all_results = []

    # 1. Spotless (with ground truth)
    data = load_spotless()
    if data:
        Y, X, coords, gt, name = data
        results = test_preprocess_modes(Y, X, coords, gt, name)
        all_results.extend(results)

    # 2. CARD
    data = load_card()
    if data:
        Y, X, coords, gt, name = data
        results = test_preprocess_modes(Y, X, coords, gt, name)
        all_results.extend(results)

    # 3. DestVI
    data = load_destvi()
    if data:
        Y, X, coords, gt, name = data
        results = test_preprocess_modes(Y, X, coords, gt, name)
        all_results.extend(results)

    # 4. Mouse Brain
    data = load_mouse_brain()
    if data:
        Y, X, coords, gt, name = data
        results = test_preprocess_modes(Y, X, coords, gt, name)
        all_results.extend(results)

    # Summary
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)

    df = pd.DataFrame(all_results)

    # Display with ground truth metrics if available
    print("\n" + df.to_string(index=False))

    # Best by dataset
    print("\n" + "-"*60)
    print("Best preprocess by dataset:")
    for dataset in df['dataset'].unique():
        sub = df[df['dataset'] == dataset]
        if 'pearson' in sub.columns and sub['pearson'].notna().any():
            best = sub.loc[sub['pearson'].idxmax()]
            print(f"  {dataset}: {best['preprocess']} (Pearson={best['pearson']:.4f})")
        else:
            best = sub.loc[sub['spatial_coherence'].idxmax()]
            print(f"  {dataset}: {best['preprocess']} (SpatialCoh={best['spatial_coherence']:.4f})")

    # Save
    df.to_csv('/Users/apple/Research/FlashDeconv/validation/all_datasets_benchmark.csv', index=False)
    print("\nResults saved to all_datasets_benchmark.csv")

    return df


if __name__ == "__main__":
    main()
