"""
Benchmark FlashDeconv on Liver case study.

Evaluation metrics (following Spotless paper):
1. AUPR: Area under precision-recall curve for portal/central vein EC detection
2. JSD: Jensen-Shannon Divergence vs snRNA-seq cell type proportions

Cell types (9 common types):
- Hepatocytes (dominant ~60%)
- Kupffer cells
- LSECs (liver sinusoidal endothelial cells)
- T cells
- B cells
- Cholangiocytes
- Portal Vein Endothelial cells
- Central Vein Endothelial cells
- Mesothelial cells
"""

import os
import sys
import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.spatial.distance import jensenshannon
from sklearn.metrics import precision_recall_curve, auc
import warnings
warnings.filterwarnings('ignore')

sys.path.insert(0, '/Users/apple/Research/FlashDeconv')
from flashdeconv import FlashDeconv

DATA_DIR = "/Users/apple/Research/FlashDeconv/validation/benchmark_data/converted"

# 9 cell types used in benchmark
CELLTYPES = [
    'B cells',
    'Central Vein Endothelial cells',
    'Cholangiocytes',
    'Hepatocytes',
    'Kupffer cells',
    'LSECs',
    'Mesothelial cells',
    'Portal Vein Endothelial cells',
    'T cells',
]

# snRNA-seq ground truth proportions (from Spotless paper)
# These are approximate from Figure 6-figure supplement 1
SNRNASEQ_PROPORTIONS = {
    'Hepatocytes': 0.65,
    'Kupffer cells': 0.10,
    'LSECs': 0.08,
    'T cells': 0.04,
    'B cells': 0.03,
    'Cholangiocytes': 0.04,
    'Portal Vein Endothelial cells': 0.02,
    'Central Vein Endothelial cells': 0.02,
    'Mesothelial cells': 0.02,
}


def load_reference(ref_name='liver_ref_9ct'):
    """Load reference scRNA-seq data."""
    prefix = f"{DATA_DIR}/{ref_name}"

    counts = mmread(f"{prefix}_counts.mtx").T.tocsr()
    with open(f"{prefix}_genes.txt") as f:
        genes = [line.strip() for line in f]
    with open(f"{prefix}_celltypes.txt") as f:
        celltypes = [line.strip() for line in f]

    return {'counts': counts, 'genes': genes, 'celltypes': celltypes}


def load_visium_sample(sample_name):
    """Load Visium spatial data with zonation annotations."""
    prefix = f"{DATA_DIR}/{sample_name}"

    counts = mmread(f"{prefix}_counts.mtx").toarray()
    with open(f"{prefix}_genes.txt") as f:
        genes = [line.strip() for line in f]
    coords = pd.read_csv(f"{prefix}_coords.csv", index_col=0)
    coords_arr = coords[['x', 'y']].values if 'x' in coords.columns else coords.iloc[:, :2].values

    # Load metadata with zonation
    metadata = pd.read_csv(f"{prefix}_metadata.csv", index_col=0)

    return {
        'counts': counts,
        'genes': genes,
        'coords': coords_arr,
        'metadata': metadata,
    }


def build_signature_matrix(ref, target_genes, target_celltypes):
    """Build signature matrix from raw counts (average per cell type)."""
    ref_genes = list(ref['genes'])
    common_genes = sorted(set(ref_genes) & set(target_genes))

    print(f"  Common genes: {len(common_genes)}")

    ref_gene_idx = {g: i for i, g in enumerate(ref_genes)}
    target_gene_idx = {g: i for i, g in enumerate(target_genes)}

    common_ref_idx = [ref_gene_idx[g] for g in common_genes]
    common_target_idx = [target_gene_idx[g] for g in common_genes]

    ref_counts = ref['counts'][:, common_ref_idx]
    if hasattr(ref_counts, 'toarray'):
        ref_counts = ref_counts.toarray()

    # Use raw counts - FlashDeconv will normalize internally
    ref_celltypes = np.array(ref['celltypes'])
    X_ref = np.zeros((len(target_celltypes), len(common_genes)))

    for i, ct in enumerate(target_celltypes):
        mask = ref_celltypes == ct
        if mask.sum() > 0:
            X_ref[i] = ref_counts[mask].mean(axis=0)
        else:
            print(f"  Warning: No cells found for {ct}")

    return X_ref, common_target_idx


def calculate_jsd(props_true, props_pred):
    """Calculate Jensen-Shannon Divergence."""
    props_true = np.array(props_true) + 1e-10
    props_pred = np.array(props_pred) + 1e-10
    props_true = props_true / props_true.sum()
    props_pred = props_pred / props_pred.sum()
    return jensenshannon(props_true, props_pred) ** 2


def calculate_aupr(y_true, y_score):
    """Calculate Area Under Precision-Recall Curve."""
    precision, recall, _ = precision_recall_curve(y_true, y_score)
    return auc(recall, precision)


def run_benchmark():
    """Run FlashDeconv benchmark on liver data."""
    print("=" * 60)
    print("Liver Case Study Benchmark")
    print("=" * 60)

    # Load reference (combined 3 protocols)
    print("\nLoading reference data (9 cell types, all protocols)...")
    ref = load_reference('liver_ref_9ct')
    print(f"  Cells: {ref['counts'].shape[0]}")
    print(f"  Genes: {len(ref['genes'])}")
    print(f"  Cell types: {len(set(ref['celltypes']))}")

    # Sort celltypes to match expected order
    unique_cts = sorted(set(ref['celltypes']))
    print(f"  Types: {unique_cts}")

    # Visium samples
    samples = [
        'liver_mouseVisium_JB01',
        'liver_mouseVisium_JB02',
        'liver_mouseVisium_JB03',
        'liver_mouseVisium_JB04',
    ]

    all_results = []

    for sample_name in samples:
        print(f"\n{'=' * 40}")
        print(f"Processing: {sample_name}")
        print('=' * 40)

        # Load spatial data
        sp = load_visium_sample(sample_name)
        print(f"  Spots: {sp['counts'].shape[0]}")
        print(f"  Genes: {len(sp['genes'])}")

        # Check zonation annotations
        zonation = sp['metadata']['zonationGroup'].values
        print(f"  Zonation groups: {np.unique(zonation)}")
        print(f"    Periportal: {(zonation == 'Periportal').sum()}")
        print(f"    Mid: {(zonation == 'Mid').sum()}")
        print(f"    Central (vein): {(zonation == 'Central').sum()}")
        print(f"    Portal (vein): {(zonation == 'Portal').sum()}")

        # Build signature matrix
        X_ref, gene_idx = build_signature_matrix(ref, sp['genes'], unique_cts)

        # Prepare spatial data - use raw counts
        Y = sp['counts'][:, gene_idx]
        coords = sp['coords']

        print(f"  Data shape: Y={Y.shape}, X_ref={X_ref.shape}")

        # FlashDeconv - parameter tuning
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

                # Mean proportions
                mean_props = pred.mean(axis=0)

                # Calculate JSD vs snRNA-seq ground truth
                gt_vec = np.array([SNRNASEQ_PROPORTIONS.get(ct, 0) for ct in unique_cts])
                jsd = calculate_jsd(gt_vec, mean_props)

                # Calculate AUPR for Portal/Central vein EC zonation
                # Portal Vein EC should be high in Periportal regions
                # Central Vein EC should be high in Central regions
                portal_idx = unique_cts.index('Portal Vein Endothelial cells')
                central_idx = unique_cts.index('Central Vein Endothelial cells')

                portal_pred = pred[:, portal_idx]
                central_pred = pred[:, central_idx]

                # AUPR calculation following Spotless methodology:
                # Only use spots in Portal or Central zones (exclude Mid and Periportal)
                aupr_mask = (zonation == 'Portal') | (zonation == 'Central')

                # Ground truth: Portal vein EC=1 in Portal spots, Central vein EC=1 in Central spots
                portal_label = (zonation[aupr_mask] == 'Portal').astype(int)
                central_label = (zonation[aupr_mask] == 'Central').astype(int)

                portal_pred_filtered = portal_pred[aupr_mask]
                central_pred_filtered = central_pred[aupr_mask]

                # Calculate AUPR (higher is better)
                aupr_portal = calculate_aupr(portal_label, portal_pred_filtered)
                aupr_central = calculate_aupr(central_label, central_pred_filtered)
                aupr_mean = (aupr_portal + aupr_central) / 2

                all_results.append({
                    'sample': sample_name,
                    'config': config['name'],
                    'jsd': jsd,
                    'aupr_portal': aupr_portal,
                    'aupr_central': aupr_central,
                    'aupr_mean': aupr_mean,
                    'Hepatocytes_prop': mean_props[unique_cts.index('Hepatocytes')],
                })

                print(f"    JSD: {jsd:.4f}")
                print(f"    AUPR (portal): {aupr_portal:.4f}")
                print(f"    AUPR (central): {aupr_central:.4f}")
                print(f"    AUPR (mean): {aupr_mean:.4f}")
                print(f"    Hepatocytes: {mean_props[unique_cts.index('Hepatocytes')]:.3f}")

            except Exception as e:
                print(f"    Error: {e}")
                import traceback
                traceback.print_exc()

    # Summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)

    if all_results:
        df = pd.DataFrame(all_results)

        # Average across samples
        summary = df.groupby('config').agg({
            'jsd': 'mean',
            'aupr_mean': 'mean',
            'aupr_portal': 'mean',
            'aupr_central': 'mean',
            'Hepatocytes_prop': 'mean',
        }).round(4)

        print("\nFlashDeconv Results (avg across 4 Visium slides):")
        print(summary)

        # Compare with other methods (from Spotless paper)
        print("\n" + "-" * 40)
        print("Comparison with other methods:")
        print("-" * 40)

        other_methods_jsd = {
            'rctd': 0.0334,
            'cell2location': 0.0352,
            'spatialdwls': 0.0646,
            'music': 0.0659,
            'destvi': 0.0947,
            'nnls': 0.1056,
            'seurat': 0.1065,
            'stride': 0.2138,
            'spotlight': 0.2397,
            'dstg': 0.2463,
            'stereoscope': 0.3643,
            'tangram': 0.4084,
        }

        other_methods_aupr = {
            'spotlight': 0.9732,
            'cell2location': 0.9422,
            'rctd': 0.8865,
            'music': 0.8058,
            'spatialdwls': 0.7978,
            'tangram': 0.7513,
            'stride': 0.6088,
            'destvi': 0.6039,
            'dstg': 0.5637,
            'stereoscope': 0.5151,
            'nnls': 0.5000,
            'seurat': 0.5000,
        }

        # Best FlashDeconv config
        best_config = df.groupby('config').agg({'jsd': 'mean', 'aupr_mean': 'mean'})
        best_jsd_config = best_config['jsd'].idxmin()
        flash_jsd = best_config.loc[best_jsd_config, 'jsd']
        flash_aupr = best_config.loc[best_jsd_config, 'aupr_mean']

        print(f"\nFlashDeconv ({best_jsd_config}):")
        print(f"  JSD:  {flash_jsd:.4f}")
        print(f"  AUPR: {flash_aupr:.4f}")

        # JSD ranking
        all_jsd = {**other_methods_jsd, 'FlashDeconv': flash_jsd}
        ranked_jsd = sorted(all_jsd.items(), key=lambda x: x[1])
        print("\nJSD Ranking (lower is better):")
        for i, (method, val) in enumerate(ranked_jsd, 1):
            marker = " *** FlashDeconv" if method == 'FlashDeconv' else ""
            print(f"  {i:2d}. {method:15s}: {val:.4f}{marker}")

        # AUPR ranking
        all_aupr = {**other_methods_aupr, 'FlashDeconv': flash_aupr}
        ranked_aupr = sorted(all_aupr.items(), key=lambda x: -x[1])
        print("\nAUPR Ranking (higher is better):")
        for i, (method, val) in enumerate(ranked_aupr, 1):
            marker = " *** FlashDeconv" if method == 'FlashDeconv' else ""
            print(f"  {i:2d}. {method:15s}: {val:.4f}{marker}")

        # Save results
        output_file = '/Users/apple/Research/FlashDeconv/validation/liver_benchmark_results.csv'
        df.to_csv(output_file, index=False)
        print(f"\nResults saved to: {output_file}")

    return all_results


if __name__ == '__main__':
    results = run_benchmark()
