"""
Diagnostic experiments for melanoma collinearity analysis.

Hypothesis 1: Averaging Bias - merging dilutes marker gene signals
Hypothesis 2: L1 Regularization - FlashDeconv auto-selects best state per spot

Author: Analysis following reviewer suggestion
"""

import sys
import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.spatial.distance import jensenshannon
from pathlib import Path
import json
import warnings
warnings.filterwarnings('ignore')

sys.path.insert(0, '/Users/apple/Research/FlashDeconv')
from flashdeconv import FlashDeconv

# Paths
DATA_DIR = Path("/Users/apple/Research/FlashDeconv/validation/benchmark_data/converted")
OUTPUT_DIR = Path("/Users/apple/Research/FlashDeconv/validation/results")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Malignant states
MALIGNANT_STATES = [
    'melanocytic/oxphos', 'neural-like', 'immune-like', 'stem-like',
    'stress-like (hypoxia/UPR)', 'RNA-processing', 'mesenchymal'
]

NON_MALIGNANT = ['B cell', 'CAF', 'DC', 'EC', 'Monocyte/macrophage', 'pDC', 'Pericyte', 'T/NK cell']

# Ground truth
MC_GROUND_TRUTH = {
    'Bcell': 0.005, 'CAF': 0.012, 'EC': 0.032, 'Melanocytic': 0.848,
    'Mono/Mac': 0.039, 'Pericyte': 0.017, 'Tcell': 0.047,
}
EVAL_MAP = {
    'B cell': 'Bcell', 'CAF': 'CAF', 'EC': 'EC',
    'Monocyte/macrophage': 'Mono/Mac', 'Pericyte': 'Pericyte', 'T/NK cell': 'Tcell',
}
EVAL_CELLTYPES = list(MC_GROUND_TRUTH.keys())


def load_reference():
    """Load reference data."""
    prefix = DATA_DIR / "melanoma_ref"
    counts = mmread(f"{prefix}_counts.mtx")
    with open(f"{prefix}_genes.txt") as f:
        genes = [line.strip() for line in f]
    with open(f"{prefix}_celltypes.txt") as f:
        celltypes = [line.strip() for line in f]
    return counts, genes, celltypes


def load_spatial(sample_id):
    """Load spatial data."""
    prefix = DATA_DIR / f"melanoma_visium_sample{sample_id:02d}"
    counts = mmread(f"{prefix}_counts.mtx").toarray()
    with open(f"{prefix}_genes.txt") as f:
        genes = [line.strip() for line in f]
    coords = pd.read_csv(f"{prefix}_coords.csv", index_col=0)
    coords_arr = coords[['x', 'y']].values if 'x' in coords.columns else coords.iloc[:, :2].values
    return counts, genes, coords_arr


def build_signature(counts_sparse, celltypes, genes, target_cts):
    """Build signature matrix."""
    X_csc = counts_sparse.tocsc()
    celltypes_arr = np.array(celltypes)
    signature = np.zeros((len(target_cts), X_csc.shape[0]), dtype=np.float64)
    for i, ct in enumerate(target_cts):
        mask = (celltypes_arr == ct)
        if mask.sum() > 0:
            signature[i] = X_csc[:, mask].mean(axis=1).A1
    return signature


def build_merged_signature(counts_sparse, celltypes, genes):
    """Build signature with malignant merged."""
    X_csc = counts_sparse.tocsc()
    celltypes_arr = np.array(celltypes)
    merged_cts = sorted(NON_MALIGNANT) + ['Malignant']
    signature = np.zeros((len(merged_cts), X_csc.shape[0]), dtype=np.float64)

    for i, ct in enumerate(merged_cts):
        if ct == 'Malignant':
            mask = np.zeros(len(celltypes_arr), dtype=bool)
            for mal_ct in MALIGNANT_STATES:
                mask |= (celltypes_arr == mal_ct)
        else:
            mask = (celltypes_arr == ct)
        if mask.sum() > 0:
            signature[i] = X_csc[:, mask].mean(axis=1).A1

    return signature, merged_cts


def align_genes(Y, sp_genes, X, ref_genes):
    """Align genes."""
    common = sorted(set(sp_genes) & set(ref_genes))
    sp_idx = {g: i for i, g in enumerate(sp_genes)}
    ref_idx = {g: i for i, g in enumerate(ref_genes)}
    return Y[:, [sp_idx[g] for g in common]], X[:, [ref_idx[g] for g in common]], common


def aggregate_to_eval(props, cell_types):
    """Aggregate to evaluation cell types."""
    result = {}
    mal_sum = 0
    for i, ct in enumerate(cell_types):
        if ct in MALIGNANT_STATES or ct == 'Malignant':
            mal_sum += props[i]
        elif ct in EVAL_MAP:
            result[EVAL_MAP[ct]] = props[i]
    result['Melanocytic'] = mal_sum
    total = sum(result.get(ct, 0) for ct in EVAL_CELLTYPES)
    if total > 0:
        for ct in EVAL_CELLTYPES:
            result[ct] = result.get(ct, 0) / total
    return np.array([result.get(ct, 0) for ct in EVAL_CELLTYPES])


def calculate_jsd(p_true, p_pred):
    p = np.array(p_true) + 1e-10
    q = np.array(p_pred) + 1e-10
    return jensenshannon(p / p.sum(), q / q.sum()) ** 2


def run_flashdeconv(Y, X, coords, rho_sparsity=0.01, verbose=False):
    """Run FlashDeconv with specified parameters."""
    model = FlashDeconv(
        sketch_dim=min(512, Y.shape[1]),
        lambda_spatial="auto",
        preprocess="log_cpm",
        n_hvg=2000,
        n_markers_per_type=50,
        rho_sparsity=rho_sparsity,
        max_iter=100,
        tol=1e-4,
        random_state=42,
        verbose=verbose,
    )
    props = model.fit_transform(Y, X, coords)
    return props, model


def main():
    print("=" * 70)
    print("DIAGNOSTIC EXPERIMENTS: Why does Full mode beat Merged?")
    print("=" * 70)

    # Load data
    print("\nLoading data...")
    ref_counts, ref_genes, ref_celltypes = load_reference()
    full_cts = sorted(set(ref_celltypes))

    # Build signatures
    X_full = build_signature(ref_counts, ref_celltypes, ref_genes, full_cts)
    X_merged, merged_cts = build_merged_signature(ref_counts, ref_celltypes, ref_genes)

    gt_vec = np.array([MC_GROUND_TRUTH[ct] for ct in EVAL_CELLTYPES])

    # Find malignant indices in full reference
    mal_idx_full = [full_cts.index(ct) for ct in MALIGNANT_STATES if ct in full_cts]
    print(f"Malignant state indices in full ref: {mal_idx_full}")
    print(f"Malignant states: {[full_cts[i] for i in mal_idx_full]}")

    # Use sample 2 for detailed analysis
    sample_id = 2
    print(f"\n{'=' * 60}")
    print(f"Analyzing Sample {sample_id}")
    print("=" * 60)

    Y_raw, sp_genes, coords = load_spatial(sample_id)
    Y_full, X_full_aligned, common_genes = align_genes(Y_raw, sp_genes, X_full, ref_genes)
    Y_merged, X_merged_aligned, _ = align_genes(Y_raw, sp_genes, X_merged, ref_genes)

    print(f"Spots: {Y_full.shape[0]}, Genes: {Y_full.shape[1]}")

    # =========================================================
    # HYPOTHESIS 2: L1 Regularization Effect
    # =========================================================
    print("\n" + "=" * 60)
    print("HYPOTHESIS 2: L1 Regularization Effect")
    print("=" * 60)

    results_h2 = []

    # Experiment 2a: Full + L1 (default)
    print("\n[2a] Full Mode + L1 (rho=0.01)...")
    props_full_l1, model_full_l1 = run_flashdeconv(Y_full, X_full_aligned, coords, rho_sparsity=0.01)
    eval_full_l1 = aggregate_to_eval(props_full_l1.mean(axis=0), full_cts)
    jsd_full_l1 = calculate_jsd(gt_vec, eval_full_l1)
    print(f"  JSD: {jsd_full_l1:.4f}, Melanocytic: {eval_full_l1[EVAL_CELLTYPES.index('Melanocytic')]*100:.1f}%")

    # Experiment 2b: Full + No L1
    print("\n[2b] Full Mode + No L1 (rho=0)...")
    props_full_no_l1, model_full_no_l1 = run_flashdeconv(Y_full, X_full_aligned, coords, rho_sparsity=0)
    eval_full_no_l1 = aggregate_to_eval(props_full_no_l1.mean(axis=0), full_cts)
    jsd_full_no_l1 = calculate_jsd(gt_vec, eval_full_no_l1)
    print(f"  JSD: {jsd_full_no_l1:.4f}, Melanocytic: {eval_full_no_l1[EVAL_CELLTYPES.index('Melanocytic')]*100:.1f}%")

    # Experiment 2c: Merged (for comparison)
    print("\n[2c] Merged Mode (rho=0.01)...")
    props_merged, model_merged = run_flashdeconv(Y_merged, X_merged_aligned, coords, rho_sparsity=0.01)
    eval_merged = aggregate_to_eval(props_merged.mean(axis=0), merged_cts)
    jsd_merged = calculate_jsd(gt_vec, eval_merged)
    print(f"  JSD: {jsd_merged:.4f}, Melanocytic: {eval_merged[EVAL_CELLTYPES.index('Melanocytic')]*100:.1f}%")

    print("\n" + "-" * 40)
    print("L1 Ablation Summary:")
    print("-" * 40)
    print(f"  Full + L1:     JSD = {jsd_full_l1:.4f}")
    print(f"  Full + No L1:  JSD = {jsd_full_no_l1:.4f}")
    print(f"  Merged:        JSD = {jsd_merged:.4f}")

    if jsd_full_l1 < jsd_full_no_l1:
        print("\n  => L1 helps! (Full+L1 better than Full+NoL1)")
    else:
        print("\n  => L1 does NOT explain the difference")

    # =========================================================
    # SPARSITY ANALYSIS: Does L1 select specific malignant states?
    # =========================================================
    print("\n" + "=" * 60)
    print("SPARSITY ANALYSIS: Malignant state selection per spot")
    print("=" * 60)

    # Extract malignant proportions per spot
    mal_props_per_spot = props_full_l1[:, mal_idx_full]  # N x 7

    # Count how many states are "active" per spot (> 1% threshold)
    threshold = 0.01
    active_per_spot = (mal_props_per_spot > threshold).sum(axis=1)

    print(f"\nNumber of active malignant states per spot (threshold > {threshold*100}%):")
    for n_active in range(8):
        count = (active_per_spot == n_active).sum()
        pct = count / len(active_per_spot) * 100
        print(f"  {n_active} states: {count:4d} spots ({pct:5.1f}%)")

    # Which state dominates?
    dominant_state = mal_props_per_spot.argmax(axis=1)
    print(f"\nDominant malignant state distribution:")
    for i, ct in enumerate([full_cts[j] for j in mal_idx_full]):
        count = (dominant_state == i).sum()
        pct = count / len(dominant_state) * 100
        print(f"  {ct:30s}: {count:4d} spots ({pct:5.1f}%)")

    # Average proportion when a state is dominant
    print(f"\nAverage proportion of dominant state:")
    dominant_props = mal_props_per_spot[np.arange(len(dominant_state)), dominant_state]
    print(f"  Mean: {dominant_props.mean()*100:.1f}%")
    print(f"  Std:  {dominant_props.std()*100:.1f}%")

    # =========================================================
    # HYPOTHESIS 1: Residual Analysis
    # =========================================================
    print("\n" + "=" * 60)
    print("HYPOTHESIS 1: Residual Analysis (Reconstruction Error)")
    print("=" * 60)

    # Get the preprocessed data for residual calculation
    # We need to compute Y_sketch - beta @ X_sketch
    # But FlashDeconv stores these internally, let's approximate with raw data

    # Compute reconstruction error in log-CPM space
    def log_cpm(M):
        row_sums = M.sum(axis=1, keepdims=True)
        row_sums = np.where(row_sums == 0, 1, row_sums)
        return np.log1p(M / row_sums * 1e4)

    Y_log = log_cpm(Y_full)
    X_full_log = log_cpm(X_full_aligned)
    X_merged_log = log_cpm(X_merged_aligned)

    # Reconstruction: Y_hat = beta @ X
    Y_hat_full = props_full_l1 @ X_full_log
    Y_hat_merged = props_merged @ X_merged_log

    # Residuals
    residual_full = np.linalg.norm(Y_log - Y_hat_full, axis=1)
    residual_merged = np.linalg.norm(Y_log - Y_hat_merged, axis=1)

    print(f"\nReconstruction error (per-spot L2 norm in log-CPM space):")
    print(f"  Full Mode:   mean = {residual_full.mean():.2f}, std = {residual_full.std():.2f}")
    print(f"  Merged Mode: mean = {residual_merged.mean():.2f}, std = {residual_merged.std():.2f}")

    if residual_merged.mean() > residual_full.mean():
        improvement = (residual_merged.mean() - residual_full.mean()) / residual_merged.mean() * 100
        print(f"\n  => Merged has {improvement:.1f}% higher reconstruction error!")
        print("     This confirms the 'Averaging Bias' hypothesis.")

    # =========================================================
    # SUMMARY
    # =========================================================
    print("\n" + "=" * 70)
    print("SUMMARY OF FINDINGS")
    print("=" * 70)

    summary = {
        'jsd_full_l1': float(jsd_full_l1),
        'jsd_full_no_l1': float(jsd_full_no_l1),
        'jsd_merged': float(jsd_merged),
        'melanocytic_full_l1': float(eval_full_l1[EVAL_CELLTYPES.index('Melanocytic')]),
        'melanocytic_full_no_l1': float(eval_full_no_l1[EVAL_CELLTYPES.index('Melanocytic')]),
        'melanocytic_merged': float(eval_merged[EVAL_CELLTYPES.index('Melanocytic')]),
        'residual_full_mean': float(residual_full.mean()),
        'residual_merged_mean': float(residual_merged.mean()),
        'active_states_distribution': {str(i): int((active_per_spot == i).sum()) for i in range(8)},
        'l1_helps': bool(jsd_full_l1 < jsd_full_no_l1),
        'averaging_bias_confirmed': bool(residual_merged.mean() > residual_full.mean()),
    }

    print(f"""
HYPOTHESIS 1 (Averaging Bias):
  - Full reconstruction error:   {residual_full.mean():.2f}
  - Merged reconstruction error: {residual_merged.mean():.2f}
  - Confirmed: {'YES' if summary['averaging_bias_confirmed'] else 'NO'}

HYPOTHESIS 2 (L1 Regularization):
  - Full + L1:     JSD = {jsd_full_l1:.4f}
  - Full + No L1:  JSD = {jsd_full_no_l1:.4f}
  - L1 helps: {'YES' if summary['l1_helps'] else 'NO'}

SPARSITY:
  - Spots with only 1 active malignant state: {(active_per_spot == 1).sum()} ({(active_per_spot == 1).sum()/len(active_per_spot)*100:.1f}%)
  - Spots with 2-3 active states: {((active_per_spot >= 2) & (active_per_spot <= 3)).sum()} ({((active_per_spot >= 2) & (active_per_spot <= 3)).sum()/len(active_per_spot)*100:.1f}%)
  - Average dominant state proportion: {dominant_props.mean()*100:.1f}%
""")

    # Save results
    with open(OUTPUT_DIR / "diagnostic_melanoma_results.json", 'w') as f:
        json.dump(summary, f, indent=2)
    print(f"Results saved to: {OUTPUT_DIR / 'diagnostic_melanoma_results.json'}")

    return summary


if __name__ == "__main__":
    summary = main()
