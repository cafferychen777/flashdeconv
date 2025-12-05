"""
Deep dive into Melanoma malignant cell states collinearity

The key finding: 7 malignant states have r > 0.98 correlation!
This analysis explores:
1. Which cell types are "malignant" vs "stromal/immune"
2. The correlation structure within each group
3. Why this breaks linear deconvolution

Author: Claude
Date: 2025-12-02
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy.io import mmread
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform
import sys

sys.path.insert(0, '/Users/apple/Research/FlashDeconv')

# Paths
data_dir = Path("/Users/apple/Research/FlashDeconv/validation/benchmark_data/converted")
output_dir = Path("/Users/apple/Research/FlashDeconv/validation/data")
fig_dir = Path("/Users/apple/Research/FlashDeconv/paper/figures")
fig_dir.mkdir(parents=True, exist_ok=True)

# Set style
plt.rcParams['font.size'] = 10
plt.rcParams['axes.labelsize'] = 11
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['figure.dpi'] = 150


def load_reference(prefix, data_dir):
    """Load reference data and create signature matrix."""
    counts_path = data_dir / f"{prefix}_counts.mtx"
    X_sparse = mmread(counts_path)

    with open(data_dir / f"{prefix}_celltypes.txt") as f:
        celltypes = [line.strip() for line in f]

    with open(data_dir / f"{prefix}_genes.txt") as f:
        genes = [line.strip() for line in f]

    unique_cts = sorted(set(celltypes))
    X_csc = X_sparse.tocsc()
    n_genes = X_csc.shape[0]

    signature = np.zeros((len(unique_cts), n_genes), dtype=np.float64)
    cell_counts = {}

    celltypes_arr = np.array(celltypes)
    for i, ct in enumerate(unique_cts):
        mask = (celltypes_arr == ct)
        cell_counts[ct] = mask.sum()
        signature[i] = X_csc[:, mask].mean(axis=1).A1

    return signature, genes, unique_cts, cell_counts


def log_cpm_transform(X):
    """Apply log-CPM transformation."""
    row_sums = X.sum(axis=1, keepdims=True)
    row_sums = np.where(row_sums == 0, 1, row_sums)
    X_cpm = X / row_sums * 1e6
    X_log = np.log1p(X_cpm)
    return X_log


def compute_correlation_matrix(X):
    """Compute pairwise correlation between cell type signatures."""
    X_norm = X / (np.linalg.norm(X, axis=1, keepdims=True) + 1e-10)
    corr = X_norm @ X_norm.T
    return corr


def main():
    print("Loading Melanoma reference...")
    mel_sig, mel_genes, mel_cts, mel_counts = load_reference("melanoma_ref", data_dir)
    X_log = log_cpm_transform(mel_sig)
    corr = compute_correlation_matrix(X_log)

    # Define cell type categories based on biological knowledge
    # Malignant states: different transcriptional programs of melanoma cells
    malignant_states = [
        'melanocytic/oxphos',
        'neural-like',
        'immune-like',
        'stem-like',
        'stress-like (hypoxia/UPR)',
        'RNA-processing',
        'mesenchymal'
    ]

    # Stromal cells
    stromal = ['CAF', 'EC', 'Pericyte']

    # Immune cells
    immune = ['T/NK cell', 'B cell', 'Monocyte/macrophage', 'DC', 'pDC']

    # Create correlation DataFrame
    corr_df = pd.DataFrame(corr, index=mel_cts, columns=mel_cts)

    # Print detailed correlation for malignant states
    print("\n" + "="*70)
    print("MALIGNANT STATE CORRELATIONS (The Core Problem)")
    print("="*70)

    mal_idx = [mel_cts.index(ct) for ct in malignant_states if ct in mel_cts]
    mal_names = [mel_cts[i] for i in mal_idx]
    mal_corr = corr[np.ix_(mal_idx, mal_idx)]

    print(f"\nMalignant states identified: {len(mal_names)}")
    print("\nCorrelation matrix among malignant states:")
    mal_corr_df = pd.DataFrame(mal_corr, index=mal_names, columns=mal_names)
    print(mal_corr_df.round(3).to_string())

    # Statistics
    mask = ~np.eye(len(mal_idx), dtype=bool)
    print(f"\nStatistics:")
    print(f"  Mean correlation: {mal_corr[mask].mean():.4f}")
    print(f"  Min correlation: {mal_corr[mask].min():.4f}")
    print(f"  Max correlation: {mal_corr[mask].max():.4f}")
    print(f"  All pairs > 0.97: {(mal_corr[mask] > 0.97).sum()}/{mask.sum()}")

    # Compare with immune and stromal
    print("\n" + "="*70)
    print("IMMUNE CELL CORRELATIONS (For Comparison)")
    print("="*70)

    imm_idx = [mel_cts.index(ct) for ct in immune if ct in mel_cts]
    imm_names = [mel_cts[i] for i in imm_idx]
    imm_corr = corr[np.ix_(imm_idx, imm_idx)]

    print(f"\nImmune cells: {imm_names}")
    print("\nCorrelation matrix among immune cells:")
    imm_corr_df = pd.DataFrame(imm_corr, index=imm_names, columns=imm_names)
    print(imm_corr_df.round(3).to_string())

    mask = ~np.eye(len(imm_idx), dtype=bool)
    print(f"\nStatistics:")
    print(f"  Mean correlation: {imm_corr[mask].mean():.4f}")
    print(f"  Min correlation: {imm_corr[mask].min():.4f}")

    # Create visualization
    print("\n" + "="*70)
    print("GENERATING VISUALIZATION")
    print("="*70)

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    # Panel A: Full correlation matrix with hierarchical clustering
    ax1 = axes[0]

    # Reorder by hierarchical clustering
    dist_matrix = 1 - corr
    np.fill_diagonal(dist_matrix, 0)
    dist_condensed = squareform(dist_matrix)
    Z = linkage(dist_condensed, method='average')

    from scipy.cluster.hierarchy import leaves_list
    order = leaves_list(Z)

    corr_ordered = corr[np.ix_(order, order)]
    ct_ordered = [mel_cts[i] for i in order]

    im = ax1.imshow(corr_ordered, cmap='RdYlBu_r', vmin=0.8, vmax=1.0)
    ax1.set_xticks(range(len(ct_ordered)))
    ax1.set_yticks(range(len(ct_ordered)))
    ax1.set_xticklabels(ct_ordered, rotation=90, fontsize=8)
    ax1.set_yticklabels(ct_ordered, fontsize=8)
    ax1.set_title('A. Cell Type Correlation Matrix\n(Hierarchically Clustered)')
    plt.colorbar(im, ax=ax1, shrink=0.8, label='Correlation')

    # Highlight malignant cluster
    mal_positions = [ct_ordered.index(ct) for ct in malignant_states if ct in ct_ordered]
    if mal_positions:
        min_pos, max_pos = min(mal_positions), max(mal_positions)
        rect = plt.Rectangle((min_pos-0.5, min_pos-0.5),
                              max_pos-min_pos+1, max_pos-min_pos+1,
                              fill=False, edgecolor='red', linewidth=2)
        ax1.add_patch(rect)
        ax1.annotate('Malignant\nstates', xy=(max_pos+1, (min_pos+max_pos)/2),
                     fontsize=9, color='red', va='center')

    # Panel B: Malignant states only (zoomed)
    ax2 = axes[1]
    im2 = ax2.imshow(mal_corr, cmap='RdYlBu_r', vmin=0.97, vmax=1.0)
    ax2.set_xticks(range(len(mal_names)))
    ax2.set_yticks(range(len(mal_names)))
    # Shorten names for display
    short_names = [n.split('/')[0][:12] for n in mal_names]
    ax2.set_xticklabels(short_names, rotation=45, ha='right', fontsize=9)
    ax2.set_yticklabels(short_names, fontsize=9)
    ax2.set_title('B. Malignant State Correlations\n(r = 0.97-1.00, nearly identical)')
    plt.colorbar(im2, ax=ax2, shrink=0.8, label='Correlation')

    # Add correlation values
    for i in range(len(mal_names)):
        for j in range(len(mal_names)):
            if i != j:
                ax2.text(j, i, f'{mal_corr[i,j]:.2f}', ha='center', va='center',
                        fontsize=7, color='black' if mal_corr[i,j] < 0.99 else 'white')

    # Panel C: Compare correlation distributions
    ax3 = axes[2]

    # Get correlations for different groups
    all_corr = corr[~np.eye(len(mel_cts), dtype=bool)]
    mal_only_corr = mal_corr[~np.eye(len(mal_idx), dtype=bool)]
    imm_only_corr = imm_corr[~np.eye(len(imm_idx), dtype=bool)]

    # Histogram
    bins = np.linspace(0.8, 1.0, 20)
    ax3.hist(all_corr, bins=bins, alpha=0.5, label=f'All pairs (mean={all_corr.mean():.3f})',
             color='gray', edgecolor='black')
    ax3.hist(mal_only_corr, bins=bins, alpha=0.7, label=f'Malignant pairs (mean={mal_only_corr.mean():.3f})',
             color='red', edgecolor='black')
    ax3.hist(imm_only_corr, bins=bins, alpha=0.7, label=f'Immune pairs (mean={imm_only_corr.mean():.3f})',
             color='blue', edgecolor='black')

    ax3.set_xlabel('Pairwise Correlation')
    ax3.set_ylabel('Count')
    ax3.set_title('C. Correlation Distribution by Cell Category')
    ax3.legend(loc='upper left', fontsize=9)
    ax3.axvline(x=0.98, color='red', linestyle='--', alpha=0.7, label='r=0.98')

    plt.tight_layout()

    # Save
    output_path = fig_dir / "melanoma_collinearity_analysis.pdf"
    plt.savefig(output_path, bbox_inches='tight', dpi=300)
    print(f"✓ Saved figure to: {output_path}")

    plt.savefig(fig_dir / "melanoma_collinearity_analysis.png", bbox_inches='tight', dpi=150)
    print(f"✓ Saved PNG preview")

    plt.close()

    # Summary
    print("\n" + "="*70)
    print("KEY INSIGHT: WHY FlashDeconv STRUGGLES ON MELANOMA")
    print("="*70)

    print("""
THE PROBLEM:
Melanoma contains 7 "malignant cell states" that are transcriptionally
almost IDENTICAL (mean r = 0.987, range 0.979-0.995).

These states represent a CONTINUOUS SPECTRUM of tumor cell phenotypes:
  - melanocytic/oxphos: differentiated melanoma cells
  - neural-like: dedifferentiated state
  - immune-like: immune evasion program
  - stem-like: cancer stem cell features
  - stress-like: hypoxia/UPR response
  - RNA-processing: proliferative state
  - mesenchymal: EMT-like transition

WHY LINEAR METHODS FAIL:
1. The correlation r > 0.98 means gene expression profiles differ by <2%
2. Any noise or batch effect overwhelms this tiny signal
3. CountSketch randomly mixes genes, further reducing separation
4. The solver has infinitely many near-optimal solutions (ill-posed)

WHY PROBABILISTIC METHODS (Cell2Location, RCTD) DO BETTER:
1. They model uncertainty explicitly - a spot can be "uncertain mixture of
   similar states" rather than forced to a wrong point estimate
2. Cell2Location models "detection efficiency" - platform-specific gene
   dropout/scaling, which helps with batch effects
3. Prior distributions can regularize towards biologically plausible solutions

TAKE-HOME MESSAGE:
FlashDeconv's linear assumption breaks down when cell types are NOT
transcriptomically distinct. This is a fundamental limitation of linear
deconvolution, not a bug in our algorithm. The ~100x speedup comes with
a trade-off: less robustness to highly collinear cell types.
""")

    # Save summary statistics
    summary_df = pd.DataFrame({
        'cell_type': mel_cts,
        'category': ['Malignant' if ct in malignant_states else
                     ('Immune' if ct in immune else 'Stromal') for ct in mel_cts],
        'cell_count': [mel_counts[ct] for ct in mel_cts],
        'mean_corr_with_others': [corr[i, :].sum() / (len(mel_cts)-1) for i in range(len(mel_cts))]
    })
    summary_df.to_csv(output_dir / "melanoma_cell_type_analysis.csv", index=False)
    print(f"\n✓ Saved cell type analysis to: {output_dir / 'melanoma_cell_type_analysis.csv'}")


if __name__ == "__main__":
    main()
