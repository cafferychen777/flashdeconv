"""
Benchmark: Full (15 cell types) vs Merged (9 cell types) on Melanoma

This experiment validates the Discussion recommendation:
"When cell types are highly correlated (r > 0.98), users should merge them."

Hypothesis: Merging 7 malignant states into 1 "Melanocytic" should improve JSD.

Ground truth: Molecular Cartography proportions from Spotless paper
"""

import sys
import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.spatial.distance import jensenshannon
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

sys.path.insert(0, '/Users/apple/Research/FlashDeconv')
from flashdeconv import FlashDeconv

# Paths
DATA_DIR = Path("/Users/apple/Research/FlashDeconv/validation/benchmark_data/converted")
OUTPUT_DIR = Path("/Users/apple/Research/FlashDeconv/validation/results")
FIG_DIR = Path("/Users/apple/Research/FlashDeconv/paper/figures")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# 7 malignant states in original reference
MALIGNANT_STATES = [
    'melanocytic/oxphos', 'neural-like', 'immune-like', 'stem-like',
    'stress-like (hypoxia/UPR)', 'RNA-processing', 'mesenchymal'
]

# Non-malignant cell types
NON_MALIGNANT = ['B cell', 'CAF', 'DC', 'EC', 'Monocyte/macrophage', 'pDC', 'Pericyte', 'T/NK cell']

# Cell type mapping for evaluation (to match ground truth)
EVAL_MAP = {
    'B cell': 'Bcell',
    'CAF': 'CAF',
    'EC': 'EC',
    'Monocyte/macrophage': 'Mono/Mac',
    'Pericyte': 'Pericyte',
    'T/NK cell': 'Tcell',
}

# Ground truth from Molecular Cartography (Spotless paper)
MC_GROUND_TRUTH = {
    'Bcell': 0.005,
    'CAF': 0.012,
    'EC': 0.032,
    'Melanocytic': 0.848,
    'Mono/Mac': 0.039,
    'Pericyte': 0.017,
    'Tcell': 0.047,
}
EVAL_CELLTYPES = list(MC_GROUND_TRUTH.keys())


def load_reference():
    """Load raw reference data."""
    prefix = DATA_DIR / "melanoma_ref"
    counts = mmread(f"{prefix}_counts.mtx")  # genes x cells

    with open(f"{prefix}_genes.txt") as f:
        genes = [line.strip() for line in f]
    with open(f"{prefix}_celltypes.txt") as f:
        celltypes = [line.strip() for line in f]

    return counts, genes, celltypes


def load_spatial(sample_id):
    """Load spatial data."""
    prefix = DATA_DIR / f"melanoma_visium_sample{sample_id:02d}"

    counts = mmread(f"{prefix}_counts.mtx").toarray()  # spots x genes
    with open(f"{prefix}_genes.txt") as f:
        genes = [line.strip() for line in f]
    coords = pd.read_csv(f"{prefix}_coords.csv", index_col=0)
    coords_arr = coords[['x', 'y']].values if 'x' in coords.columns else coords.iloc[:, :2].values

    return counts, genes, coords_arr


def build_signature_matrix(counts_sparse, celltypes, genes, target_celltypes):
    """Build signature matrix for given cell types."""
    X_csc = counts_sparse.tocsc()
    celltypes_arr = np.array(celltypes)

    signature = np.zeros((len(target_celltypes), X_csc.shape[0]), dtype=np.float64)

    for i, ct in enumerate(target_celltypes):
        mask = (celltypes_arr == ct)
        if mask.sum() > 0:
            signature[i] = X_csc[:, mask].mean(axis=1).A1

    return signature, genes


def build_merged_signature(counts_sparse, celltypes, genes, method='signature_avg'):
    """Build signature matrix with malignant states merged.

    Parameters
    ----------
    method : str
        'cell_avg': Average all malignant cells (weighted by cell count)
        'signature_avg': Average the 7 malignant signatures (equal weight)
    """
    X_csc = counts_sparse.tocsc()
    celltypes_arr = np.array(celltypes)

    # Merged cell types: non-malignant + "Malignant"
    merged_cts = sorted(NON_MALIGNANT) + ['Malignant']

    signature = np.zeros((len(merged_cts), X_csc.shape[0]), dtype=np.float64)

    for i, ct in enumerate(merged_cts):
        if ct == 'Malignant':
            if method == 'cell_avg':
                # Average all malignant cells (weighted by cell count)
                mask = np.zeros(len(celltypes_arr), dtype=bool)
                for mal_ct in MALIGNANT_STATES:
                    mask |= (celltypes_arr == mal_ct)
                if mask.sum() > 0:
                    signature[i] = X_csc[:, mask].mean(axis=1).A1
            else:  # signature_avg
                # Average the 7 malignant signatures (equal weight per state)
                mal_sigs = []
                for mal_ct in MALIGNANT_STATES:
                    mask = (celltypes_arr == mal_ct)
                    if mask.sum() > 0:
                        mal_sigs.append(X_csc[:, mask].mean(axis=1).A1)
                if mal_sigs:
                    signature[i] = np.mean(mal_sigs, axis=0)
        else:
            mask = (celltypes_arr == ct)
            if mask.sum() > 0:
                signature[i] = X_csc[:, mask].mean(axis=1).A1

    return signature, genes, merged_cts


def align_genes(Y, spatial_genes, X, ref_genes):
    """Align genes between spatial and reference."""
    common = sorted(set(spatial_genes) & set(ref_genes))

    sp_idx = {g: i for i, g in enumerate(spatial_genes)}
    ref_idx = {g: i for i, g in enumerate(ref_genes)}

    sp_common = [sp_idx[g] for g in common]
    ref_common = [ref_idx[g] for g in common]

    return Y[:, sp_common], X[:, ref_common], common


def aggregate_to_eval(props, cell_types):
    """Aggregate predictions to evaluation cell types (7 types)."""
    result = {}

    # Sum malignant states
    mal_sum = 0
    for i, ct in enumerate(cell_types):
        if ct in MALIGNANT_STATES:
            mal_sum += props[i]
        elif ct == 'Malignant':
            mal_sum += props[i]
        elif ct in EVAL_MAP:
            result[EVAL_MAP[ct]] = props[i]

    result['Melanocytic'] = mal_sum

    # Remove DC and pDC (not in ground truth)
    # Renormalize
    total = sum(result.get(ct, 0) for ct in EVAL_CELLTYPES)
    if total > 0:
        for ct in EVAL_CELLTYPES:
            result[ct] = result.get(ct, 0) / total

    return np.array([result.get(ct, 0) for ct in EVAL_CELLTYPES])


def calculate_jsd(props_true, props_pred):
    """Calculate Jensen-Shannon Divergence."""
    p = np.array(props_true) + 1e-10
    q = np.array(props_pred) + 1e-10
    p, q = p / p.sum(), q / q.sum()
    return jensenshannon(p, q) ** 2


def run_deconvolution(Y, X, coords, verbose=False):
    """Run FlashDeconv."""
    model = FlashDeconv(
        sketch_dim=min(512, Y.shape[1]),
        lambda_spatial="auto",
        preprocess="log_cpm",
        n_hvg=2000,
        n_markers_per_type=50,
        rho_sparsity=0.01,
        max_iter=100,
        tol=1e-4,
        random_state=42,
        verbose=verbose,
    )
    return model.fit_transform(Y, X, coords)


def main():
    print("=" * 70)
    print("MELANOMA MERGING EXPERIMENT")
    print("Hypothesis: Merging collinear malignant states improves accuracy")
    print("=" * 70)

    # Load reference
    print("\nLoading reference data...")
    ref_counts, ref_genes, ref_celltypes = load_reference()
    print(f"  Cells: {ref_counts.shape[1]}")
    print(f"  Genes: {len(ref_genes)}")
    print(f"  Cell types: {len(set(ref_celltypes))}")

    # Build signature matrices
    print("\nBuilding signature matrices...")

    # Full: all 15 cell types
    full_cts = sorted(set(ref_celltypes))
    X_full, _ = build_signature_matrix(ref_counts, ref_celltypes, ref_genes, full_cts)
    print(f"  Full reference: {len(full_cts)} cell types")

    # Merged: 8 non-malignant + 1 Malignant = 9 cell types
    X_merged, _, merged_cts = build_merged_signature(ref_counts, ref_celltypes, ref_genes)
    print(f"  Merged reference: {len(merged_cts)} cell types")

    # Ground truth
    gt_vec = np.array([MC_GROUND_TRUTH[ct] for ct in EVAL_CELLTYPES])
    print(f"\nGround truth (Molecular Cartography):")
    for ct, val in MC_GROUND_TRUTH.items():
        print(f"  {ct}: {val*100:.1f}%")

    # Test on each sample
    samples = [2, 3, 4]
    results = []

    for sample_id in samples:
        print(f"\n{'=' * 60}")
        print(f"Sample {sample_id}")
        print("=" * 60)

        # Load spatial data
        Y_raw, sp_genes, coords = load_spatial(sample_id)
        print(f"  Spots: {Y_raw.shape[0]}, Genes: {Y_raw.shape[1]}")

        # === Mode 1: Full (15 cell types) ===
        print("\n  [1] Full mode (15 cell types)...")
        Y_full, X_full_aligned, _ = align_genes(Y_raw, sp_genes, X_full, ref_genes)
        print(f"      Common genes: {Y_full.shape[1]}")

        try:
            pred_full = run_deconvolution(Y_full, X_full_aligned, coords)
            mean_props_full = pred_full.mean(axis=0)
            eval_props_full = aggregate_to_eval(mean_props_full, full_cts)
            jsd_full = calculate_jsd(gt_vec, eval_props_full)
            print(f"      JSD: {jsd_full:.4f}")
            print(f"      Melanocytic: {eval_props_full[EVAL_CELLTYPES.index('Melanocytic')]*100:.1f}%")
        except Exception as e:
            print(f"      Error: {e}")
            jsd_full = np.nan
            eval_props_full = np.zeros(len(EVAL_CELLTYPES))

        # === Mode 2a: Merged (cell average) ===
        print("\n  [2a] Merged mode (cell average)...")
        X_merged_cell, _, merged_cts_cell = build_merged_signature(
            ref_counts, ref_celltypes, ref_genes, method='cell_avg')
        Y_merged, X_merged_aligned, _ = align_genes(Y_raw, sp_genes, X_merged_cell, ref_genes)
        print(f"      Common genes: {Y_merged.shape[1]}")

        try:
            pred_merged = run_deconvolution(Y_merged, X_merged_aligned, coords)
            mean_props_merged = pred_merged.mean(axis=0)
            eval_props_merged = aggregate_to_eval(mean_props_merged, merged_cts_cell)
            jsd_merged_cell = calculate_jsd(gt_vec, eval_props_merged)
            print(f"      JSD: {jsd_merged_cell:.4f}")
            print(f"      Melanocytic: {eval_props_merged[EVAL_CELLTYPES.index('Melanocytic')]*100:.1f}%")
        except Exception as e:
            print(f"      Error: {e}")
            jsd_merged_cell = np.nan
            eval_props_merged = np.zeros(len(EVAL_CELLTYPES))

        # === Mode 2b: Merged (signature average) ===
        print("\n  [2b] Merged mode (signature average)...")
        X_merged_sig, _, merged_cts_sig = build_merged_signature(
            ref_counts, ref_celltypes, ref_genes, method='signature_avg')
        Y_merged2, X_merged_aligned2, _ = align_genes(Y_raw, sp_genes, X_merged_sig, ref_genes)

        try:
            pred_merged2 = run_deconvolution(Y_merged2, X_merged_aligned2, coords)
            mean_props_merged2 = pred_merged2.mean(axis=0)
            eval_props_merged2 = aggregate_to_eval(mean_props_merged2, merged_cts_sig)
            jsd_merged_sig = calculate_jsd(gt_vec, eval_props_merged2)
            print(f"      JSD: {jsd_merged_sig:.4f}")
            print(f"      Melanocytic: {eval_props_merged2[EVAL_CELLTYPES.index('Melanocytic')]*100:.1f}%")
        except Exception as e:
            print(f"      Error: {e}")
            jsd_merged_sig = np.nan
            eval_props_merged2 = np.zeros(len(EVAL_CELLTYPES))

        # Store results
        results.append({
            'sample': sample_id,
            'jsd_full': jsd_full,
            'jsd_merged_cell': jsd_merged_cell,
            'jsd_merged_sig': jsd_merged_sig,
            'full_melanocytic': eval_props_full[EVAL_CELLTYPES.index('Melanocytic')],
            'merged_cell_melanocytic': eval_props_merged[EVAL_CELLTYPES.index('Melanocytic')],
            'merged_sig_melanocytic': eval_props_merged2[EVAL_CELLTYPES.index('Melanocytic')],
        })

    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)

    df = pd.DataFrame(results)

    mean_jsd_full = df['jsd_full'].mean()
    mean_jsd_merged_cell = df['jsd_merged_cell'].mean()
    mean_jsd_merged_sig = df['jsd_merged_sig'].mean()

    print(f"\nAverage JSD (lower is better):")
    print(f"  Full (15 types):            {mean_jsd_full:.4f}")
    print(f"  Merged (cell avg):          {mean_jsd_merged_cell:.4f}")
    print(f"  Merged (signature avg):     {mean_jsd_merged_sig:.4f}")

    print(f"\nAverage Melanocytic proportion (ground truth: 84.8%):")
    print(f"  Full:            {df['full_melanocytic'].mean()*100:.1f}%")
    print(f"  Merged (cell):   {df['merged_cell_melanocytic'].mean()*100:.1f}%")
    print(f"  Merged (sig):    {df['merged_sig_melanocytic'].mean()*100:.1f}%")

    print(f"\nPer-sample JSD:")
    print(df[['sample', 'jsd_full', 'jsd_merged_cell', 'jsd_merged_sig']].to_string(index=False))

    # Comparison with other methods
    print("\n" + "-" * 40)
    print("Comparison with other methods:")
    print("-" * 40)

    other_methods = {
        'cell2location': 0.0002,
        'rctd': 0.0033,
        'spotlight': 0.0050,
        'spatialdwls': 0.0063,
        'stereoscope': 0.0081,
        'nnls': 0.0232,
    }

    # Use best FlashDeconv result
    best_flash = min(mean_jsd_full, mean_jsd_merged_cell, mean_jsd_merged_sig)
    all_jsd = {
        **other_methods,
        'FlashDeconv (full)': mean_jsd_full,
        'FlashDeconv (merged-cell)': mean_jsd_merged_cell,
        'FlashDeconv (merged-sig)': mean_jsd_merged_sig,
    }

    for method, jsd in sorted(all_jsd.items(), key=lambda x: x[1]):
        marker = " <--" if 'FlashDeconv' in method else ""
        print(f"  {method:28s}: {jsd:.4f}{marker}")

    # Save results
    df.to_csv(OUTPUT_DIR / "melanoma_merged_comparison.csv", index=False)
    print(f"\nResults saved to: {OUTPUT_DIR / 'melanoma_merged_comparison.csv'}")

    # Save summary for visualization
    summary = {
        'mean_jsd_full': float(mean_jsd_full),
        'mean_jsd_merged_cell': float(mean_jsd_merged_cell),
        'mean_jsd_merged_sig': float(mean_jsd_merged_sig),
        'mean_melanocytic_full': float(df['full_melanocytic'].mean()),
        'mean_melanocytic_merged_cell': float(df['merged_cell_melanocytic'].mean()),
        'mean_melanocytic_merged_sig': float(df['merged_sig_melanocytic'].mean()),
        'other_methods': other_methods,
    }

    import json
    with open(OUTPUT_DIR / "melanoma_merged_summary.json", 'w') as f:
        json.dump(summary, f, indent=2)

    print(f"\n{'=' * 70}")
    print("ANALYSIS")
    print("=" * 70)

    if mean_jsd_full < mean_jsd_merged_cell and mean_jsd_full < mean_jsd_merged_sig:
        print("""
FINDING: Full mode (15 cell types) outperforms Merged mode!

This counter-intuitive result suggests that:
1. Although 7 malignant states are highly correlated (r > 0.98),
   FlashDeconv can still correctly estimate their TOTAL contribution.
2. Merging states changes the reference signature in ways that may
   not match the actual tumor expression profile.
3. The "collinearity problem" affects STATE-LEVEL resolution, not
   the aggregate MALIGNANT proportion.

IMPLICATION FOR USERS:
- If you need to distinguish individual malignant states: use probabilistic
  methods (Cell2Location, RCTD) that model uncertainty.
- If you only need total malignant proportion: FlashDeconv with full
  reference works well despite collinearity.
""")
    else:
        print(f"""
Merging collinear cell types improves accuracy.
Best approach: {'Merged (cell avg)' if mean_jsd_merged_cell < mean_jsd_merged_sig else 'Merged (signature avg)'}
""")

    return df


if __name__ == "__main__":
    df = main()
