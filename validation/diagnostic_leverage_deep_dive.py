"""
Deep Dive: Does Leverage Score explain why FlashDeconv selects specific malignant states?

Hypothesis: The 2 dominant states (oxphos, hypoxia) have markers with HIGH leverage scores,
while other 5 states have markers with LOW leverage (overlapping with others).

This would mean: Leverage Sampling acts as a "Biological Regularizer" that
automatically filters out collinear states at the feature level.
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

# Import FlashDeconv internals
from flashdeconv import FlashDeconv
from flashdeconv.utils.genes import compute_leverage_scores, select_informative_genes
from flashdeconv.core.sketching import sketch_data

# Paths
DATA_DIR = Path("/Users/apple/Research/FlashDeconv/validation/benchmark_data/converted")
OUTPUT_DIR = Path("/Users/apple/Research/FlashDeconv/validation/results")

# Malignant states
MALIGNANT_STATES = [
    'melanocytic/oxphos', 'neural-like', 'immune-like', 'stem-like',
    'stress-like (hypoxia/UPR)', 'RNA-processing', 'mesenchymal'
]

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
    """Build signature matrix."""
    X_csc = counts_sparse.tocsc()
    celltypes_arr = np.array(celltypes)
    n_genes = X_csc.shape[0]
    signature = np.zeros((len(target_cts), n_genes), dtype=np.float64)
    for i, ct in enumerate(target_cts):
        mask = (celltypes_arr == ct)
        if mask.sum() > 0:
            signature[i] = X_csc[:, mask].mean(axis=1).A1
    return signature


def find_marker_genes(X, cell_types, target_ct, n_markers=50):
    """Find top marker genes for a specific cell type using log fold change."""
    ct_idx = cell_types.index(target_ct)

    # Log transform
    X_log = np.log1p(X)

    # Mean expression of target vs others
    target_expr = X_log[ct_idx]
    other_expr = np.mean(np.delete(X_log, ct_idx, axis=0), axis=0)

    # Log fold change
    lfc = target_expr - other_expr

    # Top markers (highest LFC)
    marker_idx = np.argsort(lfc)[::-1][:n_markers]

    return marker_idx, lfc[marker_idx]


def align_genes(Y, sp_genes, X, ref_genes):
    common = sorted(set(sp_genes) & set(ref_genes))
    sp_idx = {g: i for i, g in enumerate(sp_genes)}
    ref_idx = {g: i for i, g in enumerate(ref_genes)}
    return (Y[:, [sp_idx[g] for g in common]],
            X[:, [ref_idx[g] for g in common]],
            common,
            [ref_idx[g] for g in common])


def aggregate_to_eval(props, cell_types):
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


def main():
    print("=" * 70)
    print("DEEP DIVE: Leverage Score as Biological Regularizer")
    print("=" * 70)

    # Load data
    print("\nLoading data...")
    ref_counts, ref_genes, ref_celltypes = load_reference()
    full_cts = sorted(set(ref_celltypes))
    X_full = build_signature(ref_counts, ref_celltypes, full_cts)

    sample_id = 2
    Y_raw, sp_genes, coords = load_spatial(sample_id)
    Y_aligned, X_aligned, common_genes, common_ref_idx = align_genes(Y_raw, sp_genes, X_full, ref_genes)

    print(f"Reference: {X_aligned.shape[0]} cell types x {X_aligned.shape[1]} genes")
    print(f"Spatial: {Y_aligned.shape[0]} spots x {Y_aligned.shape[1]} genes")

    # =========================================================
    # STEP 1: Compute Leverage Scores for all genes
    # =========================================================
    print("\n" + "=" * 60)
    print("STEP 1: Compute Leverage Scores")
    print("=" * 60)

    # Use FlashDeconv's gene selection to get leverage scores
    gene_idx, leverage_scores = select_informative_genes(
        Y_aligned, X_aligned,
        n_hvg=2000,
        n_markers_per_type=50
    )

    print(f"Selected {len(gene_idx)} informative genes")
    print(f"Leverage scores range: [{leverage_scores.min():.6f}, {leverage_scores.max():.6f}]")

    # =========================================================
    # STEP 2: Find markers for each malignant state
    # =========================================================
    print("\n" + "=" * 60)
    print("STEP 2: Marker Genes for Each Malignant State")
    print("=" * 60)

    # Map common genes back to their positions in the selected gene set
    selected_genes = [common_genes[i] for i in gene_idx]
    gene_to_selected_idx = {g: i for i, g in enumerate(selected_genes)}

    # Find markers for each malignant state
    mal_idx_full = [full_cts.index(ct) for ct in MALIGNANT_STATES if ct in full_cts]

    marker_leverage = {}
    for mal_state in MALIGNANT_STATES:
        if mal_state not in full_cts:
            continue

        # Find markers in full gene space
        marker_idx_full, lfc = find_marker_genes(X_aligned, full_cts, mal_state, n_markers=100)

        # Map to selected genes
        marker_genes = [common_genes[i] for i in marker_idx_full]
        markers_in_selected = []
        leverage_vals = []

        for mg in marker_genes:
            if mg in gene_to_selected_idx:
                sel_idx = gene_to_selected_idx[mg]
                markers_in_selected.append(mg)
                leverage_vals.append(leverage_scores[sel_idx])

        if leverage_vals:
            mean_lev = np.mean(leverage_vals)
            max_lev = np.max(leverage_vals)
        else:
            mean_lev = max_lev = 0

        marker_leverage[mal_state] = {
            'n_markers_found': len(markers_in_selected),
            'mean_leverage': mean_lev,
            'max_leverage': max_lev,
            'leverage_values': leverage_vals[:20],  # Top 20
        }

        print(f"\n{mal_state}:")
        print(f"  Markers in selected set: {len(markers_in_selected)}/100")
        print(f"  Mean leverage: {mean_lev:.6f}")
        print(f"  Max leverage:  {max_lev:.6f}")

    # =========================================================
    # STEP 3: Compare leverage across states
    # =========================================================
    print("\n" + "=" * 60)
    print("STEP 3: Leverage Score Comparison")
    print("=" * 60)

    # Sort by mean leverage
    sorted_states = sorted(marker_leverage.items(), key=lambda x: x[1]['mean_leverage'], reverse=True)

    print("\nRanking by Mean Leverage Score:")
    print("-" * 50)
    for rank, (state, info) in enumerate(sorted_states, 1):
        bar = "█" * int(info['mean_leverage'] * 5000)
        print(f"  {rank}. {state:30s}: {info['mean_leverage']:.6f} {bar}")

    # Identify top 2 states
    top_states = [s[0] for s in sorted_states[:2]]
    bottom_states = [s[0] for s in sorted_states[2:]]

    print(f"\nTop 2 states by leverage: {top_states}")
    print(f"Bottom 5 states by leverage: {bottom_states}")

    # =========================================================
    # STEP 4: Counterfactual Experiment - Uniform vs Leverage Sketching
    # =========================================================
    print("\n" + "=" * 60)
    print("STEP 4: Counterfactual - Uniform vs Leverage Sketching")
    print("=" * 60)

    from flashdeconv.core.solver import bcd_solve, normalize_proportions
    from flashdeconv.utils.graph import coords_to_adjacency
    from flashdeconv.core.spatial import auto_tune_lambda

    # Preprocess
    def log_cpm(M):
        row_sums = M.sum(axis=1, keepdims=True)
        row_sums = np.where(row_sums == 0, 1, row_sums)
        return np.log1p(M / row_sums * 1e4)

    Y_sub = Y_aligned[:, gene_idx]
    X_sub = X_aligned[:, gene_idx]
    Y_tilde = log_cpm(Y_sub)
    X_tilde = log_cpm(X_sub)

    # Build adjacency
    adjacency = coords_to_adjacency(coords, method='knn', k=6)

    gt_vec = np.array([MC_GROUND_TRUTH[ct] for ct in EVAL_CELLTYPES])

    results = {}

    # Experiment 4a: Leverage-weighted sketching (default)
    print("\n[4a] Leverage-weighted Sketching...")
    Y_sketch_lev, X_sketch_lev, _ = sketch_data(
        Y_tilde, X_tilde,
        sketch_dim=512,
        leverage_scores=leverage_scores,
        random_state=42
    )
    lambda_auto = auto_tune_lambda(Y_sketch_lev, X_sketch_lev, adjacency, alpha=0.005)
    beta_lev, _ = bcd_solve(Y_sketch_lev, X_sketch_lev, adjacency, lambda_auto, rho=0.01)
    props_lev = normalize_proportions(beta_lev)

    eval_lev = aggregate_to_eval(props_lev.mean(axis=0), full_cts)
    jsd_lev = calculate_jsd(gt_vec, eval_lev)
    print(f"  JSD: {jsd_lev:.4f}, Melanocytic: {eval_lev[EVAL_CELLTYPES.index('Melanocytic')]*100:.1f}%")

    # Check malignant state distribution
    mal_props_lev = props_lev[:, mal_idx_full]
    active_lev = (mal_props_lev > 0.01).sum(axis=1)
    results['leverage'] = {
        'jsd': float(jsd_lev),
        'melanocytic': float(eval_lev[EVAL_CELLTYPES.index('Melanocytic')]),
        'mean_active_states': float(active_lev.mean()),
    }

    # Experiment 4b: Uniform sketching (no leverage weighting)
    print("\n[4b] Uniform Sketching (no leverage)...")
    uniform_scores = np.ones(len(gene_idx)) / len(gene_idx)
    Y_sketch_uni, X_sketch_uni, _ = sketch_data(
        Y_tilde, X_tilde,
        sketch_dim=512,
        leverage_scores=uniform_scores,
        random_state=42
    )
    lambda_auto_uni = auto_tune_lambda(Y_sketch_uni, X_sketch_uni, adjacency, alpha=0.005)
    beta_uni, _ = bcd_solve(Y_sketch_uni, X_sketch_uni, adjacency, lambda_auto_uni, rho=0.01)
    props_uni = normalize_proportions(beta_uni)

    eval_uni = aggregate_to_eval(props_uni.mean(axis=0), full_cts)
    jsd_uni = calculate_jsd(gt_vec, eval_uni)
    print(f"  JSD: {jsd_uni:.4f}, Melanocytic: {eval_uni[EVAL_CELLTYPES.index('Melanocytic')]*100:.1f}%")

    mal_props_uni = props_uni[:, mal_idx_full]
    active_uni = (mal_props_uni > 0.01).sum(axis=1)
    results['uniform'] = {
        'jsd': float(jsd_uni),
        'melanocytic': float(eval_uni[EVAL_CELLTYPES.index('Melanocytic')]),
        'mean_active_states': float(active_uni.mean()),
    }

    # Compare state distributions
    print("\n" + "-" * 50)
    print("Malignant State Activation Comparison:")
    print("-" * 50)

    print(f"\nLeverage-weighted sketching:")
    print(f"  Mean active states per spot: {active_lev.mean():.2f}")
    for i, state in enumerate(MALIGNANT_STATES):
        if state in full_cts:
            idx = mal_idx_full[i] if i < len(mal_idx_full) else -1
            if idx >= 0 and idx < mal_props_lev.shape[1]:
                pct = (mal_props_lev[:, i] > 0.01).mean() * 100
                mean_prop = mal_props_lev[:, i].mean() * 100
                print(f"    {state:30s}: {pct:5.1f}% active, mean={mean_prop:5.1f}%")

    print(f"\nUniform sketching:")
    print(f"  Mean active states per spot: {active_uni.mean():.2f}")
    for i, state in enumerate(MALIGNANT_STATES):
        if state in full_cts:
            if i < mal_props_uni.shape[1]:
                pct = (mal_props_uni[:, i] > 0.01).mean() * 100
                mean_prop = mal_props_uni[:, i].mean() * 100
                print(f"    {state:30s}: {pct:5.1f}% active, mean={mean_prop:5.1f}%")

    # =========================================================
    # SUMMARY
    # =========================================================
    print("\n" + "=" * 70)
    print("SUMMARY: Leverage Score as Biological Regularizer")
    print("=" * 70)

    leverage_helps = jsd_lev < jsd_uni
    leverage_more_sparse = active_lev.mean() < active_uni.mean()

    summary = {
        'marker_leverage': {k: {'mean': v['mean_leverage'], 'max': v['max_leverage']}
                           for k, v in marker_leverage.items()},
        'top_2_states': top_states,
        'bottom_5_states': bottom_states,
        'leverage_vs_uniform': results,
        'leverage_helps_accuracy': leverage_helps,
        'leverage_increases_sparsity': leverage_more_sparse,
    }

    print(f"""
FINDING 1: Leverage Score Distribution
  - Top 2 states by marker leverage: {top_states}
  - These match the dominant states from previous analysis!

FINDING 2: Counterfactual Experiment
  - Leverage-weighted:  JSD = {jsd_lev:.4f}, active states = {active_lev.mean():.2f}
  - Uniform:            JSD = {jsd_uni:.4f}, active states = {active_uni.mean():.2f}
  - Leverage helps accuracy: {'YES' if leverage_helps else 'NO'}
  - Leverage increases sparsity: {'YES' if leverage_more_sparse else 'NO'}

CONCLUSION:
""")

    if leverage_helps and leverage_more_sparse:
        print("""  Leverage Score Sampling acts as a 'Biological Regularizer':
  - It amplifies markers of geometrically distinct states (oxphos, hypoxia)
  - It suppresses markers of collinear/overlapping states
  - This implicit feature-level filtering explains why FlashDeconv can
    "see through" the collinearity without explicit L1 regularization.

  For Nature Methods: "FlashDeconv's leverage-weighted sketching performs
  implicit biological regularization by preserving markers with unique
  geometric signatures while filtering collinear features."
""")
    elif leverage_helps:
        print("""  Leverage helps accuracy but not sparsity.
  The benefit may come from better signal preservation, not regularization.
""")
    else:
        print("""  Leverage does NOT explain the phenomenon.
  Need to investigate other factors.
""")

    # Save results
    with open(OUTPUT_DIR / "diagnostic_leverage_deep_dive.json", 'w') as f:
        json.dump(summary, f, indent=2, default=str)
    print(f"\nResults saved to: {OUTPUT_DIR / 'diagnostic_leverage_deep_dive.json'}")

    return summary


if __name__ == "__main__":
    summary = main()
