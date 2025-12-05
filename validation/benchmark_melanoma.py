"""
Benchmark FlashDeconv on Melanoma case study.

Ground truth: Molecular Cartography proportions (from Spotless paper)
Evaluation metric: Jensen-Shannon Divergence

Cell type mapping:
- 7 malignant states are aggregated into "Melanocytic"
- Other cell types: Bcell, CAF, EC, Mono/Mac, Pericyte, Tcell
"""

import os
import sys
import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.spatial.distance import jensenshannon
import warnings
warnings.filterwarnings('ignore')

sys.path.insert(0, '/Users/apple/Research/FlashDeconv')
from flashdeconv import FlashDeconv

DATA_DIR = "/Users/apple/Research/FlashDeconv/validation/benchmark_data/converted"

# Cell type mapping: raw -> aggregated
CELLTYPE_MAP = {
    # Non-melanoma cells
    'B cell': 'Bcell',
    'CAF': 'CAF',
    'DC': 'DC',  # Will be removed
    'EC': 'EC',
    'Monocyte/macrophage': 'Mono/Mac',
    'Pericyte': 'Pericyte',
    'T/NK cell': 'Tcell',
    'pDC': 'pDC',  # Will be removed
    # Melanoma cells (aggregated)
    'RNA-processing': 'Melanocytic',
    'immune-like': 'Melanocytic',
    'melanocytic/oxphos': 'Melanocytic',
    'mesenchymal': 'Melanocytic',
    'neural-like': 'Melanocytic',
    'stem-like': 'Melanocytic',
    'stress-like (hypoxia/UPR)': 'Melanocytic',
}

# 7 cell types used in benchmark (matching Spotless ground truth)
BENCHMARK_CELLTYPES = ['Bcell', 'CAF', 'EC', 'Melanocytic', 'Mono/Mac', 'Pericyte', 'Tcell']

# All 15 cell types in reference (including DC/pDC to be removed later)
# Following Spotless: deconvolve with all types, then remove DC/pDC and renormalize
ALL_CELLTYPES = [
    'Bcell', 'CAF', 'DC', 'EC', 'Melanocytic', 'Mono/Mac', 'pDC', 'Pericyte', 'Tcell'
]

# Molecular Cartography ground truth (from Spotless paper Figure 6-figure supplement 6)
MC_GROUND_TRUTH = {
    'Bcell': 0.005,
    'CAF': 0.012,
    'EC': 0.032,
    'Melanocytic': 0.848,
    'Mono/Mac': 0.039,
    'Pericyte': 0.017,
    'Tcell': 0.047,
}


def load_reference():
    """Load and process reference scRNA-seq data."""
    prefix = f"{DATA_DIR}/melanoma_ref"

    counts = mmread(f"{prefix}_counts.mtx").T.tocsr()
    with open(f"{prefix}_genes.txt") as f:
        genes = [line.strip() for line in f]
    with open(f"{prefix}_celltypes.txt") as f:
        celltypes = [line.strip() for line in f]

    return {'counts': counts, 'genes': genes, 'celltypes': celltypes}


def load_visium_sample(sample_name):
    """Load Visium spatial data."""
    prefix = f"{DATA_DIR}/{sample_name}"

    # MTX is stored as spots × genes (R used t() before saving)
    # So we read directly without transposing
    counts = mmread(f"{prefix}_counts.mtx").toarray()
    with open(f"{prefix}_genes.txt") as f:
        genes = [line.strip() for line in f]
    coords = pd.read_csv(f"{prefix}_coords.csv", index_col=0)
    coords_arr = coords[['x', 'y']].values if 'x' in coords.columns else coords.iloc[:, :2].values

    return {'counts': counts, 'genes': genes, 'coords': coords_arr}


def build_signature_matrix(ref, target_genes, target_celltypes, normalize=True):
    """Build signature matrix with cell type aggregation.

    Parameters
    ----------
    normalize : bool
        If True, normalize each cell to CPM before averaging.
        This is important for cross-platform comparisons.
    """
    ref_genes = list(ref['genes'])
    common_genes = sorted(set(ref_genes) & set(target_genes))

    print(f"  Common genes: {len(common_genes)}")

    if len(common_genes) < 100:
        return None, None

    ref_gene_idx = {g: i for i, g in enumerate(ref_genes)}
    target_gene_idx = {g: i for i, g in enumerate(target_genes)}

    common_ref_idx = [ref_gene_idx[g] for g in common_genes]
    common_target_idx = [target_gene_idx[g] for g in common_genes]

    ref_counts = ref['counts'][:, common_ref_idx]
    if hasattr(ref_counts, 'toarray'):
        ref_counts = ref_counts.toarray()

    # Don't normalize - use raw counts, let FlashDeconv handle it
    # FlashDeconv will normalize according to preprocess parameter

    # Build aggregated signature matrix from raw counts
    ref_celltypes = np.array(ref['celltypes'])
    X_ref = np.zeros((len(target_celltypes), len(common_genes)))

    for i, tct in enumerate(target_celltypes):
        # Find all raw cell types that map to this aggregated type
        matching_raw = [raw for raw, agg in CELLTYPE_MAP.items() if agg == tct]

        mask = np.zeros(len(ref_celltypes), dtype=bool)
        for raw_ct in matching_raw:
            mask |= (ref_celltypes == raw_ct)

        if mask.sum() > 0:
            X_ref[i] = ref_counts[mask].mean(axis=0)
        else:
            print(f"  Warning: No cells found for {tct}")

    return X_ref, common_target_idx


def calculate_jsd(props_true, props_pred):
    """Calculate Jensen-Shannon Divergence."""
    props_true = np.array(props_true) + 1e-10
    props_pred = np.array(props_pred) + 1e-10
    props_true = props_true / props_true.sum()
    props_pred = props_pred / props_pred.sum()
    return jensenshannon(props_true, props_pred) ** 2


def run_benchmark():
    """Run FlashDeconv benchmark on melanoma data."""
    print("="*60)
    print("Melanoma Case Study Benchmark")
    print("="*60)

    # Load reference
    print("\nLoading reference data...")
    ref = load_reference()
    print(f"  Cells: {ref['counts'].shape[0]}")
    print(f"  Genes: {len(ref['genes'])}")
    print(f"  Cell types: {len(set(ref['celltypes']))}")

    # Visium samples
    samples = ['melanoma_visium_sample02', 'melanoma_visium_sample03', 'melanoma_visium_sample04']

    all_results = []

    for sample_name in samples:
        print(f"\n{'='*40}")
        print(f"Processing: {sample_name}")
        print('='*40)

        # Load spatial data
        sp = load_visium_sample(sample_name)
        print(f"  Spots: {sp['counts'].shape[0]}")
        print(f"  Genes: {len(sp['genes'])}")

        # Build signature matrix with ALL cell types (including DC/pDC)
        # Following Spotless methodology: deconvolve with all types, then remove DC/pDC
        X_ref, gene_idx = build_signature_matrix(ref, sp['genes'], ALL_CELLTYPES)

        if X_ref is None:
            print("  Skipping: insufficient common genes")
            continue

        # Prepare spatial data - use raw counts
        Y = sp['counts'][:, gene_idx]
        coords = sp['coords']

        print(f"  Data shape: Y={Y.shape}, X_ref={X_ref.shape} (9 cell types including DC/pDC)")

        # Run FlashDeconv - parameter tuning
        configs = [
            {"name": "default", "lambda_spatial": "auto", "sketch_dim": 512, "n_hvg": 2000},
            {"name": "no_spatial", "lambda_spatial": 0, "sketch_dim": 512, "n_hvg": 2000},
            {"name": "more_hvg", "lambda_spatial": 0, "sketch_dim": 1024, "n_hvg": 5000},
            {"name": "large", "lambda_spatial": 0, "sketch_dim": 2048, "n_hvg": 10000},
        ]

        for config in configs:
            print(f"\n  Config: {config['name']}")

            model = FlashDeconv(
                sketch_dim=min(config['sketch_dim'], Y.shape[1]),
                lambda_spatial=config['lambda_spatial'],
                preprocess="log_cpm",  # Use Seurat-style normalization
                n_hvg=config.get('n_hvg', 2000),
                rho_sparsity=0,
                max_iter=200,
                verbose=False,
                random_state=42,
            )

            try:
                pred = model.fit_transform(Y, X_ref, coords)

                # Aggregate predictions to get mean proportions (9 cell types)
                mean_props_all = pred.mean(axis=0)

                # Following Spotless: remove DC and pDC, renormalize
                # ALL_CELLTYPES order: ['Bcell', 'CAF', 'DC', 'EC', 'Melanocytic', 'Mono/Mac', 'pDC', 'Pericyte', 'Tcell']
                # Remove indices 2 (DC) and 6 (pDC)
                keep_idx = [i for i, ct in enumerate(ALL_CELLTYPES) if ct not in ['DC', 'pDC']]
                mean_props_filtered = mean_props_all[keep_idx]

                # Renormalize to sum to 1
                mean_props = mean_props_filtered / (mean_props_filtered.sum() + 1e-10)

                # Calculate JSD vs ground truth
                # BENCHMARK_CELLTYPES order matches filtered ALL_CELLTYPES (without DC/pDC)
                gt_vec = np.array([MC_GROUND_TRUTH[ct] for ct in BENCHMARK_CELLTYPES])
                jsd = calculate_jsd(gt_vec, mean_props)

                all_results.append({
                    'sample': sample_name,
                    'config': config['name'],
                    'jsd': jsd,
                    **{ct: mean_props[i] for i, ct in enumerate(BENCHMARK_CELLTYPES)}
                })

                print(f"    JSD: {jsd:.4f}")
                print(f"    Melanocytic prop: {mean_props[BENCHMARK_CELLTYPES.index('Melanocytic')]:.3f}")

            except Exception as e:
                print(f"    Error: {e}")

    # Summary
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)

    if all_results:
        df = pd.DataFrame(all_results)

        # Average across samples
        summary = df.groupby('config').agg({
            'jsd': 'mean',
            'Melanocytic': 'mean',
        }).round(4)

        print("\nFlashDeconv Results (avg across 3 Visium slides):")
        print(summary)

        # Compare with other methods
        print("\n" + "-"*40)
        print("Comparison with other methods (JSD, lower is better):")
        print("-"*40)

        other_methods = {
            'cell2location': 0.0002,
            'rctd': 0.0033,
            'spotlight': 0.0050,
            'spatialdwls': 0.0063,
            'stereoscope': 0.0081,
            'nnls': 0.0232,
        }

        flash_jsd = df[df['config'] == 'default']['jsd'].mean()
        all_jsd = {**other_methods, 'FlashDeconv': flash_jsd}

        for method, jsd in sorted(all_jsd.items(), key=lambda x: x[1]):
            marker = " <-- FlashDeconv" if method == 'FlashDeconv' else ""
            print(f"  {method:20s}: {jsd:.4f}{marker}")

        # Save results
        output_file = '/Users/apple/Research/FlashDeconv/validation/melanoma_benchmark_results.csv'
        df.to_csv(output_file, index=False)
        print(f"\nResults saved to: {output_file}")

    return all_results


if __name__ == '__main__':
    results = run_benchmark()
