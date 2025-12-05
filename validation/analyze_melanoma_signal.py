"""
Corrected Signal-Leverage Analysis

Key Insight: The original analysis used min-max normalization which was misleading.
Raw leverage values are very uniform (0.009-0.011 range), meaning leverage
does NOT strongly differentiate between states.

The selection of Oxphos and Hypoxia is primarily driven by SIGNAL (presence
in spatial data), not leverage score differences.
"""

import sys
import numpy as np
import pandas as pd
from scipy.io import mmread
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

sys.path.insert(0, '/Users/apple/Research/FlashDeconv')

from flashdeconv.utils.genes import compute_leverage_scores, select_informative_genes

# Paths
DATA_DIR = Path("/Users/apple/Research/FlashDeconv/validation/benchmark_data/converted")
OUTPUT_DIR = Path("/Users/apple/Research/FlashDeconv/validation/results")
FIG_DIR = Path("/Users/apple/Research/FlashDeconv/paper/figures")

# Malignant states
MALIGNANT_STATES = [
    'melanocytic/oxphos', 'neural-like', 'immune-like', 'stem-like',
    'stress-like (hypoxia/UPR)', 'RNA-processing', 'mesenchymal'
]

# Short names for plotting
SHORT_NAMES = {
    'melanocytic/oxphos': 'Oxphos',
    'neural-like': 'Neural',
    'immune-like': 'Immune',
    'stem-like': 'Stem',
    'stress-like (hypoxia/UPR)': 'Hypoxia',
    'RNA-processing': 'RNA-proc',
    'mesenchymal': 'Mesench',
}


def load_reference():
    prefix = DATA_DIR / "melanoma_ref"
    counts = mmread(f"{prefix}_counts.mtx")
    with open(f"{prefix}_genes.txt") as f:
        genes = [line.strip() for line in f]
    with open(f"{prefix}_celltypes.txt") as f:
        celltypes = [line.strip() for line in f]
    return counts, genes, celltypes


def load_spatial(sample_id):
    prefix = DATA_DIR / f"melanoma_visium_sample{sample_id:02d}"
    counts = mmread(f"{prefix}_counts.mtx").toarray()
    with open(f"{prefix}_genes.txt") as f:
        genes = [line.strip() for line in f]
    coords = pd.read_csv(f"{prefix}_coords.csv", index_col=0)
    coords_arr = coords[['x', 'y']].values if 'x' in coords.columns else coords.iloc[:, :2].values
    return counts, genes, coords_arr


def build_signature(counts_sparse, celltypes, target_cts):
    X_csc = counts_sparse.tocsc()
    celltypes_arr = np.array(celltypes)
    n_genes = X_csc.shape[0]
    signature = np.zeros((len(target_cts), n_genes), dtype=np.float64)
    for i, ct in enumerate(target_cts):
        mask = (celltypes_arr == ct)
        if mask.sum() > 0:
            signature[i] = X_csc[:, mask].mean(axis=1).A1
    return signature


def align_genes(Y, sp_genes, X, ref_genes):
    common = sorted(set(sp_genes) & set(ref_genes))
    sp_idx = {g: i for i, g in enumerate(sp_genes)}
    ref_idx = {g: i for i, g in enumerate(ref_genes)}
    return (Y[:, [sp_idx[g] for g in common]],
            X[:, [ref_idx[g] for g in common]],
            common)


def log_cpm(M):
    row_sums = M.sum(axis=1, keepdims=True)
    row_sums = np.where(row_sums == 0, 1, row_sums)
    return np.log1p(M / row_sums * 1e4)


def compute_projection_strength(Y, signature_vector):
    """
    Compute how strongly a cell type signature projects onto spatial data.
    """
    v = signature_vector / (np.linalg.norm(signature_vector) + 1e-10)
    projection = Y @ v
    return np.linalg.norm(projection)


def main():
    print("=" * 70)
    print("CORRECTED SIGNAL-LEVERAGE ANALYSIS")
    print("=" * 70)

    # Load data
    print("\nLoading data...")
    ref_counts, ref_genes, ref_celltypes = load_reference()
    full_cts = sorted(set(ref_celltypes))
    X_full = build_signature(ref_counts, ref_celltypes, full_cts)

    # Load all samples
    samples = [2, 3, 4]
    all_results = []

    for sample_id in samples:
        print(f"\nProcessing Sample {sample_id}...")

        Y_raw, sp_genes, coords = load_spatial(sample_id)
        Y_aligned, X_aligned, common_genes = align_genes(Y_raw, sp_genes, X_full, ref_genes)

        # Preprocess
        Y_log = log_cpm(Y_aligned)
        X_log = log_cpm(X_aligned)

        # Get informative genes and leverage scores
        gene_idx, leverage_scores = select_informative_genes(
            Y_aligned, X_aligned, n_hvg=2000, n_markers_per_type=50
        )

        Y_sub = Y_log[:, gene_idx]
        X_sub = X_log[:, gene_idx]

        # For each malignant state
        for state in MALIGNANT_STATES:
            if state not in full_cts:
                continue

            state_idx = full_cts.index(state)
            sig_vector = X_sub[state_idx]

            # Projection strength (signal)
            proj_strength = compute_projection_strength(Y_sub, sig_vector)

            # Leverage: mean of leverage scores for this state's top markers
            top_marker_idx = np.argsort(sig_vector)[::-1][:50]
            mean_leverage = leverage_scores[top_marker_idx].mean()

            all_results.append({
                'sample': sample_id,
                'state': state,
                'short_name': SHORT_NAMES[state],
                'leverage': mean_leverage,
                'projection_strength': proj_strength,
            })

    # Aggregate across samples
    df = pd.DataFrame(all_results)

    agg = df.groupby(['state', 'short_name']).agg({
        'leverage': 'mean',
        'projection_strength': 'mean',
    }).reset_index()

    print("\n" + "=" * 70)
    print("RAW VALUES (No misleading normalization)")
    print("=" * 70)

    # Show raw values
    print("\nState             | Leverage (raw)    | Signal (raw)")
    print("-" * 60)
    for _, row in agg.sort_values('projection_strength', ascending=False).iterrows():
        print(f"{row['short_name']:17s} | {row['leverage']:.6f}        | {row['projection_strength']:.1f}")

    # Key statistics
    lev_range = agg['leverage'].max() - agg['leverage'].min()
    lev_cv = agg['leverage'].std() / agg['leverage'].mean()  # coefficient of variation

    sig_range = agg['projection_strength'].max() - agg['projection_strength'].min()
    sig_cv = agg['projection_strength'].std() / agg['projection_strength'].mean()

    print(f"\nLeverage range: {lev_range:.6f} (CV = {lev_cv:.2%})")
    print(f"Signal range: {sig_range:.1f} (CV = {sig_cv:.2%})")

    # Normalize properly for visualization (z-score or relative)
    agg['signal_norm'] = (agg['projection_strength'] - agg['projection_strength'].min()) / \
                         (agg['projection_strength'].max() - agg['projection_strength'].min())

    # For leverage, show relative deviation from mean
    lev_mean = agg['leverage'].mean()
    agg['leverage_deviation'] = (agg['leverage'] - lev_mean) / lev_mean * 100  # percentage deviation

    print("\n" + "=" * 70)
    print("INTERPRETATION")
    print("=" * 70)
    print(f"""
Key Finding: Leverage scores are NEARLY UNIFORM across all malignant states.
  - Range: {agg['leverage'].min():.6f} to {agg['leverage'].max():.6f}
  - Coefficient of Variation: {lev_cv:.2%} (very low!)

This means leverage does NOT strongly differentiate states.
The original min-max normalization created a misleading 0-1 spread.

Selection of Oxphos and Hypoxia is primarily driven by SIGNAL (presence):
""")

    selected = ['Oxphos', 'Hypoxia']
    not_selected = ['Mesench', 'Neural', 'Immune', 'Stem', 'RNA-proc']

    selected_signal = agg[agg['short_name'].isin(selected)]['signal_norm'].mean()
    not_selected_signal = agg[agg['short_name'].isin(not_selected)]['signal_norm'].mean()

    selected_lev = agg[agg['short_name'].isin(selected)]['leverage'].mean()
    not_selected_lev = agg[agg['short_name'].isin(not_selected)]['leverage'].mean()

    print(f"  Selected states (Oxphos, Hypoxia):")
    print(f"    - Mean signal (norm): {selected_signal:.3f}")
    print(f"    - Mean leverage: {selected_lev:.6f}")
    print(f"  Not selected states:")
    print(f"    - Mean signal (norm): {not_selected_signal:.3f}")
    print(f"    - Mean leverage: {not_selected_lev:.6f}")

    print(f"""
Conclusion: Selection correlates with SIGNAL, not leverage.
  - Selected states have higher signal (presence in tissue)
  - Leverage scores are essentially the same across all states
""")

    # =========================================================
    # Generate CORRECTED visualization
    # =========================================================
    print("\n" + "=" * 70)
    print("GENERATING CORRECTED VISUALIZATION")
    print("=" * 70)

    import matplotlib.pyplot as plt

    plt.rcParams.update({
        'font.family': 'Arial',
        'font.size': 8,
        'axes.linewidth': 0.5,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'figure.dpi': 300,
    })

    fig, axes = plt.subplots(1, 2, figsize=(7, 3.5))

    # Panel A: Correct representation - bar chart showing leverage is uniform
    ax1 = axes[0]
    states_sorted = agg.sort_values('leverage', ascending=False)
    colors_bar = ['#E64B35' if n in selected else '#7E7E7E' for n in states_sorted['short_name']]

    bars = ax1.barh(range(len(states_sorted)), states_sorted['leverage'] * 1000,
                    color=colors_bar, edgecolor='black', linewidth=0.5)
    ax1.set_yticks(range(len(states_sorted)))
    ax1.set_yticklabels(states_sorted['short_name'])
    ax1.set_xlabel('Leverage Score (×10⁻³)', fontsize=9)
    ax1.set_xlim(8.5, 11.5)  # Zoomed to show the small differences
    ax1.axvline(x=lev_mean * 1000, color='gray', linestyle='--', linewidth=0.5, alpha=0.7)
    ax1.invert_yaxis()

    # Panel B: Signal (projection strength) - this shows real variation
    ax2 = axes[1]
    states_sorted = agg.sort_values('projection_strength', ascending=False)
    colors_bar = ['#E64B35' if n in selected else '#7E7E7E' for n in states_sorted['short_name']]

    bars = ax2.barh(range(len(states_sorted)), states_sorted['projection_strength'],
                    color=colors_bar, edgecolor='black', linewidth=0.5)
    ax2.set_yticks(range(len(states_sorted)))
    ax2.set_yticklabels(states_sorted['short_name'])
    ax2.set_xlabel('Projection Strength', fontsize=9)
    ax2.invert_yaxis()

    plt.tight_layout()

    # Save
    plt.savefig(FIG_DIR / "signal_leverage_corrected.pdf", format='pdf', bbox_inches='tight')
    plt.savefig(OUTPUT_DIR / "signal_leverage_corrected.png", format='png', dpi=300, bbox_inches='tight')
    print(f"Saved to: {FIG_DIR / 'signal_leverage_corrected.pdf'}")

    plt.close()

    # Save results as CSV
    agg['selected'] = agg['short_name'].isin(['Oxphos', 'Hypoxia']).map({True: 'Yes', False: 'No'})
    output_cols = ['state', 'short_name', 'leverage', 'projection_strength', 'signal_norm', 'selected']
    agg[output_cols].to_csv(OUTPUT_DIR / "melanoma_signal_analysis.csv", index=False)
    print(f"Saved to: {OUTPUT_DIR / 'melanoma_signal_analysis.csv'}")

    return agg


if __name__ == "__main__":
    agg = main()
