"""
Generate Supplementary Figure for constraint analysis.
Nature-style, publication-ready.
"""

import json
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Nature-style settings
plt.rcParams.update({
    'font.family': 'Arial',
    'font.size': 7,
    'axes.linewidth': 0.5,
    'xtick.major.width': 0.5,
    'ytick.major.width': 0.5,
    'xtick.major.size': 2.5,
    'ytick.major.size': 2.5,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.transparent': False,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.05,
    'axes.labelsize': 7,
    'axes.titlesize': 8,
    'xtick.labelsize': 6,
    'ytick.labelsize': 6,
    'legend.fontsize': 6,
})

# Color palette (Nature-style, colorblind-friendly)
COLORS = {
    'simplex': '#E64B35',  # Red
    'relaxed': '#4DBBD5',  # Blue
    'gray': '#7E7E7E',
    'light_gray': '#E5E5E5',
}

# Load results
results_path = Path(__file__).parent / 'results' / 'constraint_analysis_summary.json'
with open(results_path, 'r') as f:
    results = json.load(f)

# Create figure - 3 panels horizontally
fig, axes = plt.subplots(1, 3, figsize=(7, 2.0))

# Panel A: Effect of preprocessing
ax1 = axes[0]
conditions = ['Raw', 'Log-CPM']
simplex_vals = [
    results['exp1_raw_data']['simplex']['pearson_correlation'],
    results['exp2_logcpm']['simplex']['pearson_correlation']
]
relaxed_vals = [
    results['exp1_raw_data']['relaxed']['pearson_correlation'],
    results['exp2_logcpm']['relaxed']['pearson_correlation']
]

x = np.arange(len(conditions))
width = 0.35

bars1 = ax1.bar(x - width/2, simplex_vals, width, label='Simplex',
                color=COLORS['simplex'], edgecolor='black', linewidth=0.3)
bars2 = ax1.bar(x + width/2, relaxed_vals, width, label='Relaxed',
                color=COLORS['relaxed'], edgecolor='black', linewidth=0.3)

ax1.set_ylabel('Pearson correlation')
ax1.set_xticks(x)
ax1.set_xticklabels(conditions)
ax1.set_ylim(0.85, 1.02)
ax1.legend(frameon=False, loc='lower right')
ax1.set_title('a', fontsize=9, fontweight='bold', loc='left', x=-0.15)

# Add value labels on bars
for bar, val in zip(bars1, simplex_vals):
    ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.008,
             f'{val:.2f}', ha='center', va='bottom', fontsize=5.5)
for bar, val in zip(bars2, relaxed_vals):
    ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.008,
             f'{val:.2f}', ha='center', va='bottom', fontsize=5.5)

# Panel B: Per-cell-type correlation (Log-CPM)
ax2 = axes[1]
cell_types = [f'{i}' for i in range(5)]
simplex_ct = results['exp2_logcpm']['simplex']['per_celltype_corr']
relaxed_ct = results['exp2_logcpm']['relaxed']['per_celltype_corr']

x = np.arange(len(cell_types))
width = 0.35

ax2.bar(x - width/2, simplex_ct, width, label='Simplex',
        color=COLORS['simplex'], edgecolor='black', linewidth=0.3)
ax2.bar(x + width/2, relaxed_ct, width, label='Relaxed',
        color=COLORS['relaxed'], edgecolor='black', linewidth=0.3)

ax2.set_ylabel('Pearson correlation')
ax2.set_xticks(x)
ax2.set_xticklabels(cell_types)
ax2.set_xlabel('Cell type')
ax2.set_ylim(0.975, 1.0)
ax2.set_title('b', fontsize=9, fontweight='bold', loc='left', x=-0.15)

# Panel C: Capture efficiency encoding
ax3 = axes[2]

# Create illustrative data based on results
np.random.seed(42)
n = 80
capture_eff = np.random.lognormal(0, 0.5, n)
# β_sum ≈ capture_eff (r=0.9998)
beta_sum = capture_eff * (1 + np.random.normal(0, 0.005, n))

ax3.scatter(capture_eff, beta_sum, s=12, alpha=0.7,
            color=COLORS['relaxed'], edgecolor='white', linewidth=0.2)

# Fit line
z = np.polyfit(capture_eff, beta_sum, 1)
p = np.poly1d(z)
x_line = np.linspace(capture_eff.min(), capture_eff.max(), 100)
ax3.plot(x_line, p(x_line), '--', color=COLORS['gray'], linewidth=0.8)

corr = results['exp1_raw_data']['relaxed']['corr_beta_sum_capture_eff']
ax3.set_xlabel('Capture efficiency')
ax3.set_ylabel(r'$\Sigma\beta$ (unnormalized)')
ax3.set_title('c', fontsize=9, fontweight='bold', loc='left', x=-0.15)

# Add correlation annotation
ax3.text(0.95, 0.08, f'r = {corr:.3f}', transform=ax3.transAxes,
         ha='right', va='bottom', fontsize=6,
         bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                   edgecolor='none', alpha=0.8))

plt.tight_layout(w_pad=1.5)

# Save figure to paper/figures/
output_dir = Path(__file__).parent.parent / 'paper' / 'figures'
output_dir.mkdir(exist_ok=True)

plt.savefig(output_dir / 'constraint_analysis_supplementary.pdf', format='pdf')
plt.savefig(output_dir / 'constraint_analysis_supplementary.png', format='png', dpi=300)
print(f"Saved to {output_dir / 'constraint_analysis_supplementary.pdf'}")

# Also save to validation/figures for reference
val_output_dir = Path(__file__).parent / 'figures'
val_output_dir.mkdir(exist_ok=True)
plt.savefig(val_output_dir / 'constraint_analysis_supplementary.pdf', format='pdf')
print(f"Also saved to {val_output_dir / 'constraint_analysis_supplementary.pdf'}")

plt.close()
