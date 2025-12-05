"""
Evaluate FlashDeconv's reference dataset stability on Liver.

This script tests FlashDeconv's robustness to different reference sequencing
protocols by comparing predictions across 3 different references:
- exVivo: scRNA-seq following ex vivo digestion
- inVivo: scRNA-seq following in vivo liver perfusion
- nuclei: single-nucleus RNA-seq (snRNA-seq)

Stability is measured by calculating the JSD between predictions from different
reference datasets. Lower JSD indicates higher stability.

Based on Spotless paper lines 129-149, especially Figure 6 supplement 5.
"""

import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.spatial.distance import jensenshannon
from scipy.stats import pearsonr
import time
import sys
import os

# Add FlashDeconv to path
sys.path.insert(0, '/Users/apple/Research/FlashDeconv/src')
from flashdeconv import FlashDeconv

# Paths
DATA_DIR = '/Users/apple/Research/FlashDeconv/validation/benchmark_data/converted'

# Reference types to test
REFERENCE_TYPES = ['exVivo', 'inVivo', 'nuclei']

# 9 common cell types across all protocols
CELL_TYPES = [
    'B cells',
    'Central Vein Endothelial cells',
    'Cholangiocytes',
    'Hepatocytes',
    'Kupffer cells',
    'LSECs',
    'Mesothelial cells',
    'Portal Vein Endothelial cells',
    'T cells'
]

# Visium slides
SLIDES = ['JB01', 'JB02', 'JB03', 'JB04']

def load_reference(ref_type):
    """Load reference dataset for a specific protocol."""
    print(f"\n=== Loading {ref_type} reference ===")

    # Load counts matrix (genes × cells in MTX, need to transpose)
    counts_path = f'{DATA_DIR}/liver_ref_{ref_type}_counts.mtx'
    counts = mmread(counts_path).T.toarray()  # Transpose to (cells × genes)

    # Load genes
    genes_path = f'{DATA_DIR}/liver_ref_{ref_type}_genes.txt'
    with open(genes_path) as f:
        genes = [line.strip() for line in f]

    # Load cell type annotations
    celltypes_path = f'{DATA_DIR}/liver_ref_{ref_type}_celltypes.txt'
    with open(celltypes_path) as f:
        celltypes = [line.strip() for line in f]

    print(f"Reference shape (cells × genes): {counts.shape}")
    print(f"Cell types: {sorted(set(celltypes))}")

    # Build reference matrix (genes × cell types)
    # Average expression of each gene in each cell type
    X_ref = np.zeros((len(CELL_TYPES), len(genes)))

    for i, ct in enumerate(CELL_TYPES):
        mask = np.array([c == ct for c in celltypes])
        if mask.sum() > 0:
            X_ref[i] = counts[mask].mean(axis=0)
        else:
            print(f"WARNING: No cells found for {ct}")

    print(f"Reference matrix shape: {X_ref.shape}")
    print(f"Number of cells per type:")
    for ct in CELL_TYPES:
        n_cells = sum(1 for c in celltypes if c == ct)
        print(f"  {ct}: {n_cells}")

    return X_ref, genes

def load_visium_slide(slide_id):
    """Load Visium slide data."""
    # Load counts (already in spots × genes format)
    counts_path = f'{DATA_DIR}/liver_mouseVisium_{slide_id}_counts.mtx'
    counts = mmread(counts_path).toarray()  # (spots × genes)

    # Load genes
    genes_path = f'{DATA_DIR}/liver_mouseVisium_{slide_id}_genes.txt'
    with open(genes_path) as f:
        genes = [line.strip() for line in f]

    # Load coordinates
    coords_path = f'{DATA_DIR}/liver_mouseVisium_{slide_id}_coords.csv'
    coords_df = pd.read_csv(coords_path, index_col=0)
    coords = coords_df[['x', 'y']].values if 'x' in coords_df.columns else coords_df.iloc[:, :2].values

    return counts, genes, coords

def run_deconvolution(ref_type, slide_id):
    """Run FlashDeconv with a specific reference on a specific slide."""
    print(f"\n--- Running deconvolution: {ref_type} on {slide_id} ---")

    # Load reference
    X_ref, ref_genes = load_reference(ref_type)

    # Load spatial data
    Y, spatial_genes, coords = load_visium_slide(slide_id)

    # Find common genes
    common_genes = sorted(set(ref_genes) & set(spatial_genes))
    print(f"Common genes: {len(common_genes)}")

    # Build index dictionaries
    ref_gene_dict = {g: i for i, g in enumerate(ref_genes)}
    spatial_gene_dict = {g: i for i, g in enumerate(spatial_genes)}

    # Get indices for common genes
    ref_gene_idx = [ref_gene_dict[g] for g in common_genes]
    spatial_gene_idx = [spatial_gene_dict[g] for g in common_genes]

    # Subset to common genes
    X_ref_common = X_ref[:, ref_gene_idx]
    Y_common = Y[:, spatial_gene_idx]

    print(f"X_ref shape: {X_ref_common.shape} (cell_types × genes)")
    print(f"Y shape: {Y_common.shape} (spots × genes)")

    # Run FlashDeconv with default settings (from liver benchmark)
    # Use raw counts + log_cpm normalization
    # NOTE: fit_transform expects (Y, X_ref, coords) - NO transpose needed!
    model = FlashDeconv(
        preprocess="log_cpm",
        n_hvg=2000,
        sketch_dim=512,
        lambda_spatial="auto",
        verbose=True
    )

    start_time = time.time()
    props = model.fit_transform(Y_common, X_ref_common, coords)
    runtime = time.time() - start_time

    print(f"Runtime: {runtime:.2f}s")
    print(f"Predictions shape: {props.shape}")
    print(f"Mean proportions: {props.mean(axis=0)}")

    return props

def calculate_stability_jsd(props_dict):
    """
    Calculate pairwise JSD between predictions from different references.

    props_dict: {ref_type: {slide_id: proportions}}

    Returns: DataFrame with pairwise JSD for each slide and overall average
    """
    results = []

    ref_types = list(props_dict.keys())

    # For each slide
    for slide_id in SLIDES:
        slide_jsds = []

        # Calculate pairwise JSD between all reference pairs
        for i in range(len(ref_types)):
            for j in range(i + 1, len(ref_types)):
                ref1, ref2 = ref_types[i], ref_types[j]

                props1 = props_dict[ref1][slide_id]
                props2 = props_dict[ref2][slide_id]

                # Calculate average JSD across all spots
                jsds = []
                for spot_idx in range(props1.shape[0]):
                    jsd = jensenshannon(props1[spot_idx], props2[spot_idx]) ** 2
                    jsds.append(jsd)

                mean_jsd = np.mean(jsds)
                slide_jsds.append(mean_jsd)

                results.append({
                    'slide': slide_id,
                    'ref1': ref1,
                    'ref2': ref2,
                    'jsd': mean_jsd
                })

                print(f"{slide_id}: JSD({ref1} vs {ref2}) = {mean_jsd:.4f}")

        # Average JSD for this slide across all pairs
        avg_jsd = np.mean(slide_jsds)
        print(f"{slide_id}: Average stability JSD = {avg_jsd:.4f}")

    df = pd.DataFrame(results)

    # Overall average
    overall_avg = df['jsd'].mean()
    print(f"\n{'='*60}")
    print(f"Overall Stability JSD (FlashDeconv): {overall_avg:.4f}")
    print(f"{'='*60}")

    return df, overall_avg

def main():
    """Main evaluation pipeline."""
    print("="*80)
    print("FlashDeconv Reference Dataset Stability Evaluation on Liver")
    print("="*80)

    print("\nThis test evaluates FlashDeconv's robustness to different reference")
    print("sequencing protocols by measuring prediction consistency across:")
    print("  - exVivo: scRNA-seq (ex vivo digestion)")
    print("  - inVivo: scRNA-seq (in vivo perfusion)")
    print("  - nuclei: snRNA-seq")
    print("\nLower stability JSD indicates more robust predictions.")
    print(f"According to Spotless paper, RCTD and Seurat have the lowest JSD.\n")

    # Store all predictions
    # Structure: {ref_type: {slide_id: proportions}}
    all_predictions = {ref_type: {} for ref_type in REFERENCE_TYPES}

    # Run deconvolution for each reference × slide combination
    for ref_type in REFERENCE_TYPES:
        print(f"\n{'='*80}")
        print(f"Testing with {ref_type} reference")
        print(f"{'='*80}")

        for slide_id in SLIDES:
            props = run_deconvolution(ref_type, slide_id)
            all_predictions[ref_type][slide_id] = props

    # Calculate stability metrics
    print("\n" + "="*80)
    print("Calculating Stability Metrics")
    print("="*80)

    stability_df, overall_jsd = calculate_stability_jsd(all_predictions)

    # Save results
    output_file = '/Users/apple/Research/FlashDeconv/validation/liver_stability_results.csv'
    stability_df.to_csv(output_file, index=False)
    print(f"\nDetailed results saved to: {output_file}")

    # Summary
    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)
    print(f"\nFlashDeconv Stability JSD: {overall_jsd:.4f}")
    print("\nInterpretation:")
    print("  - This JSD measures prediction consistency across different references")
    print("  - Lower values indicate more stable/robust predictions")
    print("  - Compare with Spotless paper Figure 6 supplement 5:")
    print("    - RCTD & Seurat: lowest JSD (most stable)")
    print("    - cell2location: medium JSD")
    print("    - Other methods: higher JSD")

    # Save summary
    summary = {
        'method': 'FlashDeconv',
        'overall_stability_jsd': overall_jsd,
        'n_references': len(REFERENCE_TYPES),
        'n_slides': len(SLIDES),
        'n_cell_types': len(CELL_TYPES)
    }

    summary_file = '/Users/apple/Research/FlashDeconv/validation/liver_stability_summary.csv'
    pd.DataFrame([summary]).to_csv(summary_file, index=False)
    print(f"\nSummary saved to: {summary_file}")

if __name__ == '__main__':
    main()
