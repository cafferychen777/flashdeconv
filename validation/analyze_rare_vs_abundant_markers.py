"""
Memory-optimized version: Analyze Rare vs Abundant Cell Type Markers
Using Real FlashDeconv API with Memory Efficiency

Author: Claude
Date: 2025-12-01
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.io import mmread
from scipy.stats import mannwhitneyu
from scipy import sparse
from pathlib import Path
from collections import Counter
import sys
import gc

# Add FlashDeconv to path
sys.path.insert(0, '/Users/apple/Research/FlashDeconv')

# Import ACTUAL FlashDeconv APIs
from flashdeconv.utils.genes import compute_leverage_scores

print("=" * 80)
print("Variance vs Leverage Analysis - Memory-Optimized with Real API")
print("=" * 80)

# ========================================
# 1. Load Data Efficiently
# ========================================
print("\n" + "=" * 80)
print("Step 1: Loading Liver scRNA-seq Reference")
print("=" * 80)

data_dir = Path("/Users/apple/Research/FlashDeconv/validation/benchmark_data/converted")

# Load gene names first (small)
print("Loading gene names...")
with open(data_dir / "liver_ref_9ct_genes.txt") as f:
    genes = np.array([line.strip() for line in f])
n_genes = len(genes)
print(f"  Number of genes: {n_genes}")

# Load cell types (small)
print("Loading cell type labels...")
with open(data_dir / "liver_ref_9ct_celltypes.txt") as f:
    cell_types = np.array([line.strip() for line in f])
n_cells = len(cell_types)
print(f"  Number of cells: {n_cells}")

# Compute cell type distribution
ct_counts = Counter(cell_types)
unique_cts = sorted(set(cell_types))
n_cell_types = len(unique_cts)
ct_to_idx = {ct: i for i, ct in enumerate(unique_cts)}

print("\nCell type distribution:")
ct_abundance = {}
for ct, count in sorted(ct_counts.items(), key=lambda x: -x[1]):
    pct = 100 * count / n_cells
    ct_abundance[ct] = pct
    print(f"  {ct:40s}: {count:6d} cells ({pct:5.2f}%)")

# Categorize
rare_cts = [ct for ct, pct in ct_abundance.items() if pct < 3.0]
abundant_cts = [ct for ct, pct in ct_abundance.items() if pct > 15.0]
print(f"\n🔴 RARE: {rare_cts}")
print(f"🟢 ABUNDANT: {abundant_cts}")

# ========================================
# 2. Create Signature Matrix (Memory-Efficient)
# ========================================
print("\n" + "=" * 80)
print("Step 2: Creating Reference Signature Matrix (Memory-Efficient)")
print("=" * 80)

print("Loading counts matrix in sparse format...")
counts_coo = mmread(data_dir / "liver_ref_9ct_counts.mtx")
print(f"  Counts shape: {counts_coo.shape} (genes x cells)")
print(f"  Converting to CSC format for efficient column indexing...")
counts_csc = counts_coo.tocsc()
del counts_coo
print(f"  Format: CSC")

# Compute cell type signatures directly from sparse matrix
print(f"Computing {n_cell_types} x {n_genes} signature matrix...")
X_raw = np.zeros((n_cell_types, n_genes), dtype=np.float32)

for i, ct in enumerate(unique_cts):
    mask = (cell_types == ct)
    n_ct = mask.sum()
    # Extract columns for this cell type, average
    ct_cols = counts_csc[:, mask]
    X_raw[i] = ct_cols.mean(axis=1).A1  # .A1 converts to 1D array
    if (i + 1) % 3 == 0:
        print(f"  Processed {i+1}/{n_cell_types} cell types...")

print(f"  X_raw shape: {X_raw.shape}")

# Free memory
del counts_csc
gc.collect()

# Log-CPM normalization
print("Applying Log-CPM normalization...")
lib_size_X = X_raw.sum(axis=1, keepdims=True) + 1e-10
X_cpm = X_raw / lib_size_X * 1e4
X = np.log1p(X_cpm)
print(f"  X (normalized) shape: {X.shape}")

del X_raw, X_cpm
gc.collect()

# ========================================
# 3. Compute Leverage Scores using Real API
# ========================================
print("\n" + "=" * 80)
print("Step 3: Computing Leverage Scores using FlashDeconv API")
print("=" * 80)

print("Calling flashdeconv.utils.genes.compute_leverage_scores()...")
leverage_scores = compute_leverage_scores(
    X,
    regularization=1e-6
)
print(f"  Leverage scores shape: {leverage_scores.shape}")
print(f"  Leverage range: [{leverage_scores.min():.6f}, {leverage_scores.max():.6f}]")
print(f"  Leverage sum: {leverage_scores.sum():.4f}")

# ========================================
# 4. Compute Variance
# ========================================
print("\n" + "=" * 80)
print("Step 4: Computing Gene Variance")
print("=" * 80)

gene_variance = np.var(X, axis=0)
print(f"  Variance range: [{gene_variance.min():.4f}, {gene_variance.max():.4f}]")

# ========================================
# 5. Find Top Markers
# ========================================
print("\n" + "=" * 80)
print("Step 5: Identifying Top Markers")
print("=" * 80)

def find_top_markers(X, ct_idx, n_top=3):
    expr_in_ct = X[ct_idx]
    expr_in_others = np.delete(X, ct_idx, axis=0).mean(axis=0)
    specificity = expr_in_ct - expr_in_others
    valid_mask = expr_in_ct > 0.5
    specificity[~valid_mask] = -np.inf
    top_idx = np.argsort(specificity)[::-1][:n_top]
    return top_idx, specificity[top_idx]

all_markers = []

for ct in rare_cts + abundant_cts:
    ct_idx = ct_to_idx[ct]
    ct_pct = ct_abundance[ct]
    category = "RARE" if ct in rare_cts else "ABUNDANT"

    print(f"\n{category}: {ct} ({ct_pct:.2f}%)")

    top_idx, specificities = find_top_markers(X, ct_idx, n_top=3)

    for rank, (idx, spec) in enumerate(zip(top_idx, specificities), 1):
        gene = genes[idx]
        var = gene_variance[idx]
        lev = leverage_scores[idx]
        expr = X[ct_idx, idx]

        all_markers.append({
            'cell_type': ct,
            'abundance_pct': ct_pct,
            'category': category,
            'gene': gene,
            'rank': rank,
            'specificity': spec,
            'variance': var,
            'leverage': lev,
            'expression': expr,
        })

        print(f"  {rank}. {gene:15s}  Var={var:6.4f}  Lev={lev:.6f}")

markers_df = pd.DataFrame(all_markers)

# ========================================
# 6. Statistical Comparison
# ========================================
print("\n" + "=" * 80)
print("Step 6: Statistical Comparison")
print("=" * 80)

rare_markers = markers_df[markers_df['category'] == 'RARE']
abundant_markers = markers_df[markers_df['category'] == 'ABUNDANT']

print(f"\nRARE (n={len(rare_markers)}):")
print(f"  Variance: {rare_markers['variance'].mean():.4f}")
print(f"  Leverage: {rare_markers['leverage'].mean():.6f}")

print(f"\nABUNDANT (n={len(abundant_markers)}):")
print(f"  Variance: {abundant_markers['variance'].mean():.4f}")
print(f"  Leverage: {abundant_markers['leverage'].mean():.6f}")

var_stat, var_p = mannwhitneyu(rare_markers['variance'], abundant_markers['variance'])
lev_stat, lev_p = mannwhitneyu(rare_markers['leverage'], abundant_markers['leverage'])

print(f"\nMann-Whitney U test:")
print(f"  Variance: p={var_p:.4f}")
print(f"  Leverage: p={lev_p:.4f}")

lev_fold = rare_markers['leverage'].mean() / abundant_markers['leverage'].mean()
print(f"\nLeverage fold change: {lev_fold:.2f}x")

# ========================================
# 7. Create Figure
# ========================================
print("\n" + "=" * 80)
print("Step 7: Creating Figure")
print("=" * 80)

fig, axes = plt.subplots(1, 2, figsize=(16, 7))

colors = {'RARE': '#E64B35', 'ABUNDANT': '#00A087'}

# Panel A: Scatter
ax = axes[0]
ax.scatter(np.log10(gene_variance + 1e-6), leverage_scores, s=1, alpha=0.2, color='lightgray', rasterized=True)

for category, color in colors.items():
    subset = markers_df[markers_df['category'] == category]
    ax.scatter(np.log10(subset['variance'] + 1e-6), subset['leverage'], s=150, color=color,
               edgecolor='black', linewidth=2, label=f'{category} markers', alpha=0.8, zorder=10)

    for ct in subset['cell_type'].unique():
        top = subset[subset['cell_type'] == ct].iloc[0]
        ax.annotate(f"{top['gene']}\n({ct})", xy=(np.log10(top['variance'] + 1e-6), top['leverage']),
                   xytext=(10, 10), textcoords='offset points', fontsize=9, fontweight='bold', color=color,
                   bbox=dict(boxstyle='round,pad=0.4', facecolor='white', edgecolor=color, linewidth=1.5, alpha=0.9),
                   arrowprops=dict(arrowstyle='->', color=color, lw=1.5))

ax.set_xlabel('Log10(Gene Variance)', fontsize=13, fontweight='bold')
ax.set_ylabel('Leverage Score', fontsize=13, fontweight='bold')
ax.set_title('A. Variance vs Leverage\n(Using FlashDeconv Real API)', fontsize=14, fontweight='bold', pad=15)
ax.legend(loc='upper left', fontsize=11)
ax.grid(True, alpha=0.3, linestyle=':', linewidth=0.5)

# Panel B: Boxplot
ax = axes[1]

plot_data = []
for category in ['RARE', 'ABUNDANT']:
    subset = markers_df[markers_df['category'] == category]
    for _, row in subset.iterrows():
        plot_data.append({'Category': category, 'Metric': 'Variance', 'Value': row['variance']})
        plot_data.append({'Category': category, 'Metric': 'Leverage', 'Value': row['leverage'] * 1000})

plot_df = pd.DataFrame(plot_data)

positions = {'RARE': [0.8, 2.8], 'ABUNDANT': [1.2, 3.2]}

for category, color in colors.items():
    ax.boxplot([plot_df[(plot_df['Category']==category) & (plot_df['Metric']=='Variance')]['Value']],
               positions=[positions[category][0]], widths=0.3, patch_artist=True,
               boxprops=dict(facecolor=color, alpha=0.7), medianprops=dict(color='black', linewidth=2),
               whiskerprops=dict(color=color, linewidth=1.5), capprops=dict(color=color, linewidth=1.5))

    ax.boxplot([plot_df[(plot_df['Category']==category) & (plot_df['Metric']=='Leverage')]['Value']],
               positions=[positions[category][1]], widths=0.3, patch_artist=True,
               boxprops=dict(facecolor=color, alpha=0.7), medianprops=dict(color='black', linewidth=2),
               whiskerprops=dict(color=color, linewidth=1.5), capprops=dict(color=color, linewidth=1.5))

ax.set_xticks([1, 3])
ax.set_xticklabels(['Variance', 'Leverage (×1000)'], fontsize=12, fontweight='bold')
ax.set_ylabel('Value', fontsize=13, fontweight='bold')
ax.set_title(f'B. Comparison\n(p_var={var_p:.3f}, p_lev={lev_p:.4f})', fontsize=14, fontweight='bold', pad=15)

import matplotlib.patches as mpatches
ax.legend(handles=[mpatches.Patch(color=colors[c], alpha=0.7, label=c) for c in ['RARE', 'ABUNDANT']],
         loc='upper left', fontsize=11)
ax.grid(True, alpha=0.3, axis='y', linestyle=':', linewidth=0.5)

plt.tight_layout()

output_path = Path("/Users/apple/Research/FlashDeconv/paper/figures/supplementary_rare_vs_abundant_final.pdf")
output_path.parent.mkdir(parents=True, exist_ok=True)
plt.savefig(output_path, dpi=300, bbox_inches='tight')
print(f"✓ Saved PDF to: {output_path}")

png_path = output_path.with_suffix('.png')
plt.savefig(png_path, dpi=150, bbox_inches='tight')
print(f"✓ Saved PNG to: {png_path}")

# ========================================
# 8. Save Results
# ========================================
print("\n" + "=" * 80)
print("Step 8: Saving Complete Data for Future Use")
print("==" * 80)

data_dir_out = Path("/Users/apple/Research/FlashDeconv/validation/data")
data_dir_out.mkdir(parents=True, exist_ok=True)

# 8.1 Save marker genes with detailed metadata
markers_path = data_dir_out / "rare_vs_abundant_markers.csv"
markers_df.to_csv(markers_path, index=False)
print(f"✓ Saved marker genes to: {markers_path}")

# 8.2 Save ALL genes' variance and leverage (for plotting background)
print("Saving complete gene-level data...")
all_genes_df = pd.DataFrame({
    'gene': genes,
    'variance': gene_variance,
    'leverage': leverage_scores,
})
all_genes_path = data_dir_out / "all_genes_variance_leverage.csv"
all_genes_df.to_csv(all_genes_path, index=False)
print(f"✓ Saved all genes data to: {all_genes_path}")
print(f"  (n={len(all_genes_df)} genes)")

# 8.3 Save statistical summary
print("Saving statistical summary...")
stats_df = pd.DataFrame({
    'metric': ['variance_mean', 'leverage_mean', 'variance_p_value', 'leverage_p_value', 'leverage_fold_change'],
    'RARE': [rare_markers['variance'].mean(), rare_markers['leverage'].mean(), var_p, lev_p, lev_fold],
    'ABUNDANT': [abundant_markers['variance'].mean(), abundant_markers['leverage'].mean(), np.nan, np.nan, np.nan],
})
stats_path = data_dir_out / "rare_vs_abundant_statistics.csv"
stats_df.to_csv(stats_path, index=False)
print(f"✓ Saved statistics to: {stats_path}")

# 8.4 Save metadata for future reference
metadata = {
    'analysis_date': '2025-12-01',
    'dataset': 'Liver scRNA-seq (liver_ref_9ct)',
    'n_cells': n_cells,
    'n_genes': n_genes,
    'n_cell_types': n_cell_types,
    'normalization': 'Log-CPM (log1p(CPM/10000))',
    'leverage_regularization': 1e-6,
    'rare_threshold_pct': 3.0,
    'abundant_threshold_pct': 15.0,
    'rare_cell_types': ', '.join(rare_cts),
    'abundant_cell_types': ', '.join(abundant_cts),
    'n_rare_markers': len(rare_markers),
    'n_abundant_markers': len(abundant_markers),
    'variance_p_value': var_p,
    'leverage_p_value': lev_p,
    'leverage_fold_change': lev_fold,
    'flashdeconv_api': 'flashdeconv.utils.genes.compute_leverage_scores',
}
metadata_df = pd.DataFrame([metadata]).T
metadata_df.columns = ['value']
metadata_path = data_dir_out / "rare_vs_abundant_metadata.csv"
metadata_df.to_csv(metadata_path)
print(f"✓ Saved metadata to: {metadata_path}")

print("\n✅ All data saved successfully! Future plotting can use these files directly.")

# ========================================
# 9. Final Summary
# ========================================
print("\n" + "=" * 80)
print("FINAL SUMMARY")
print("=" * 80)

print(f"\n📊 Results (Using Real FlashDeconv API):")
print(f"  RARE markers:      Var={rare_markers['variance'].mean():.4f}, Lev={rare_markers['leverage'].mean():.6f}")
print(f"  ABUNDANT markers:  Var={abundant_markers['variance'].mean():.4f}, Lev={abundant_markers['leverage'].mean():.6f}")
print(f"\n  Statistical tests:")
print(f"    Variance:  p={var_p:.4f}  {'✓ sig' if var_p < 0.05 else '✗ ns'}")
print(f"    Leverage:  p={lev_p:.4f}  {'✓ sig' if lev_p < 0.05 else '✗ ns'}")
print(f"    Fold:      {lev_fold:.2f}x")

if lev_p < 0.05 and lev_fold > 1.5:
    print("\n🎯 ✓✓✓ NARRATIVE VALIDATED:")
    print("    → Leverage scores capture discriminative power INDEPENDENTLY of abundance")
    print("    → Variance-based methods are BLIND to rare vs abundant distinction")
    print("    → This justifies FlashDeconv's leverage-score-based sketching!")

print("\n" + "=" * 80)
print("Analysis Complete!")
print("=" * 80)
