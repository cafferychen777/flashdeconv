"""
Generate Supplementary Figure for Brain Cortex Leverage Analysis
Shows that leverage-abundance relationship is tissue-dependent

Author: Claude
Date: 2025-12-01
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# Set style
plt.rcParams['font.size'] = 11
plt.rcParams['axes.labelsize'] = 12
plt.rcParams['axes.titlesize'] = 13
plt.rcParams['figure.dpi'] = 150

# Load data
data_dir = Path("/Users/apple/Research/FlashDeconv/validation/data")
output_dir = Path("/Users/apple/Research/FlashDeconv/paper/figures")
output_dir.mkdir(parents=True, exist_ok=True)

# Load Brain Cortex stats
brain_stats = pd.read_csv(data_dir / "brain_cortex_deep_dive_ct_stats.csv")

# Sort by abundance
brain_stats = brain_stats.sort_values('abundance_pct', ascending=False)

# Create figure
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Panel A: Bar plot of leverage by cell type
ax1 = axes[0]

colors = []
for _, row in brain_stats.iterrows():
    if row['category'] == 'ABUNDANT':
        colors.append('#e74c3c')  # Red for abundant
    elif row['category'] == 'RARE':
        colors.append('#3498db')  # Blue for rare
    else:
        colors.append('#95a5a6')  # Gray for medium

bars = ax1.barh(range(len(brain_stats)), brain_stats['leverage_mean'],
                xerr=brain_stats['leverage_std'], capsize=3,
                color=colors, edgecolor='black', linewidth=0.5)

ax1.set_yticks(range(len(brain_stats)))
ax1.set_yticklabels([f"{row['cell_type']} ({row['abundance_pct']:.1f}%)"
                     for _, row in brain_stats.iterrows()], fontsize=9)
ax1.set_xlabel('Mean Leverage Score')
ax1.set_title('A. Leverage by Cell Type (Brain Cortex)')
ax1.invert_yaxis()

# Add legend
from matplotlib.patches import Patch
legend_elements = [Patch(facecolor='#e74c3c', edgecolor='black', label='Abundant (>10%)'),
                   Patch(facecolor='#95a5a6', edgecolor='black', label='Medium (3-10%)'),
                   Patch(facecolor='#3498db', edgecolor='black', label='Rare (<3%)')]
ax1.legend(handles=legend_elements, loc='lower right', fontsize=9)

# Panel B: Scatter plot of abundance vs leverage
ax2 = axes[1]

for _, row in brain_stats.iterrows():
    if row['category'] == 'ABUNDANT':
        color = '#e74c3c'
        marker = 's'
    elif row['category'] == 'RARE':
        color = '#3498db'
        marker = 'o'
    else:
        color = '#95a5a6'
        marker = '^'

    ax2.scatter(row['abundance_pct'], row['leverage_mean'],
                c=color, marker=marker, s=100, edgecolor='black', linewidth=0.5)

# Add labels for key cell types
key_types = ['Vip', 'Sst', 'Macrophage', 'Endo', 'Oligo']
for _, row in brain_stats.iterrows():
    if row['cell_type'] in key_types:
        offset = (5, 5) if row['abundance_pct'] > 5 else (5, -10)
        ax2.annotate(row['cell_type'],
                     (row['abundance_pct'], row['leverage_mean']),
                     xytext=offset, textcoords='offset points',
                     fontsize=9, style='italic')

ax2.set_xlabel('Cell Type Abundance (%)')
ax2.set_ylabel('Mean Leverage Score')
ax2.set_title('B. Leverage vs Abundance (Brain Cortex)')

# Add trend annotation
ax2.annotate('Abundant neurons\nhave HIGH leverage\n(specialized programs)',
             xy=(12, 0.0012), fontsize=9, ha='center',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

ax2.annotate('Rare support cells\nhave LOW leverage\n(generic programs)',
             xy=(2, 0.0004), fontsize=9, ha='center',
             bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.5))

plt.tight_layout()

# Save figure
output_path = output_dir / "supplementary_brain_leverage.pdf"
plt.savefig(output_path, bbox_inches='tight', dpi=300)
print(f"✓ Saved figure to: {output_path}")

# Also save as PNG for preview
plt.savefig(output_dir / "supplementary_brain_leverage.png", bbox_inches='tight', dpi=150)
print(f"✓ Saved PNG preview")

plt.close()

# Print summary statistics
print("\n" + "=" * 60)
print("BRAIN CORTEX LEVERAGE SUMMARY")
print("=" * 60)

for cat in ['ABUNDANT', 'MEDIUM', 'RARE']:
    cat_data = brain_stats[brain_stats['category'] == cat]
    if len(cat_data) > 0:
        mean_lev = cat_data['leverage_mean'].mean()
        print(f"{cat}: Mean leverage = {mean_lev:.6f} (n={len(cat_data)} cell types)")

print("\nKey observation: In Brain Cortex, ABUNDANT neurons have HIGHER")
print("leverage than RARE support cells - opposite to Liver!")
