"""
Verify which malignant states FlashDeconv actually selects in Full mode.
"""

import sys
import numpy as np
import pandas as pd
from scipy.io import mmread
from pathlib import Path

sys.path.insert(0, '/Users/apple/Research/FlashDeconv')
from flashdeconv import FlashDeconv

DATA_DIR = Path("/Users/apple/Research/FlashDeconv/validation/benchmark_data/converted")

MALIGNANT_STATES = [
    'melanocytic/oxphos', 'neural-like', 'immune-like', 'stem-like',
    'stress-like (hypoxia/UPR)', 'RNA-processing', 'mesenchymal'
]

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


def main():
    print("=" * 70)
    print("VERIFY: Which malignant states does FlashDeconv select?")
    print("=" * 70)

    ref_counts, ref_genes, ref_celltypes = load_reference()
    full_cts = sorted(set(ref_celltypes))
    X_full = build_signature(ref_counts, ref_celltypes, full_cts)

    # Get malignant state indices
    mal_indices = [full_cts.index(s) for s in MALIGNANT_STATES if s in full_cts]

    all_mal_props = {s: [] for s in MALIGNANT_STATES}

    for sample_id in [2, 3, 4]:
        print(f"\nSample {sample_id}:")

        Y_raw, sp_genes, coords = load_spatial(sample_id)
        Y_aligned, X_aligned, _ = align_genes(Y_raw, sp_genes, X_full, ref_genes)

        model = FlashDeconv(
            sketch_dim=min(512, Y_aligned.shape[1]),
            lambda_spatial="auto",
            preprocess="log_cpm",
            n_hvg=2000,
            n_markers_per_type=50,
            rho_sparsity=0.01,
            max_iter=100,
            random_state=42,
        )

        props = model.fit_transform(Y_aligned, X_aligned, coords)

        # Average proportions per malignant state
        print(f"  Mean proportions for malignant states:")
        for state in MALIGNANT_STATES:
            if state in full_cts:
                idx = full_cts.index(state)
                mean_prop = props[:, idx].mean()
                all_mal_props[state].append(mean_prop)
                print(f"    {SHORT_NAMES[state]:10s}: {mean_prop*100:5.2f}%")

    print("\n" + "=" * 70)
    print("AVERAGE ACROSS 3 SAMPLES")
    print("=" * 70)

    results = []
    for state in MALIGNANT_STATES:
        if all_mal_props[state]:
            avg = np.mean(all_mal_props[state])
            results.append((SHORT_NAMES[state], avg))

    results.sort(key=lambda x: x[1], reverse=True)

    print("\nMalignant state proportions (sorted by proportion):")
    total_mal = sum(r[1] for r in results)
    for name, prop in results:
        pct_of_mal = prop / total_mal * 100 if total_mal > 0 else 0
        print(f"  {name:10s}: {prop*100:5.2f}% (= {pct_of_mal:5.1f}% of total malignant)")

    print(f"\nTotal malignant: {total_mal*100:.1f}%")

    # Check if top 2 match our "selected" states
    top2 = [r[0] for r in results[:2]]
    print(f"\nTop 2 states: {top2}")
    print(f"Expected (from signal analysis): ['Oxphos', 'Hypoxia']")
    print(f"Match: {set(top2) == {'Oxphos', 'Hypoxia'}}")


if __name__ == "__main__":
    main()
