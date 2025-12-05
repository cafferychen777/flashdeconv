"""
Cortex Main Figure - Final Version
==================================
5x3 grid with unified colormap for negative controls
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import LinearSegmentedColormap
from scipy.ndimage import gaussian_filter1d
import matplotlib.patches as mpatches

# Load data
data = np.load("validation/level2_v3_data.npz", allow_pickle=True)
props = data['proportions']
coords = data['coordinates']
cell_types = list(data['cell_types'])

def get_idx(col):
    return cell_types.index(col) if col in cell_types else None

# Row 1: Cortical layers
layer_spatial = [
    ('L2/3', 'Ext_L23', '#2166ac'),
    ('L2-5', 'Ext_L25', '#4393c3'),
    ('L5-6', 'Ext_L56', '#f4a582'),
    ('L6', 'Ext_L6', '#d6604d'),
]

# Row 3: Context maps - distinct colormaps per cell type
context_maps = [
    ('Thalamus', 'Ext_Thal_1', 'Oranges'),
    ('Hippocampus CA1', 'Ext_Hpc_CA1', 'Greens'),
    ('Hippocampus DG', 'Ext_Hpc_DG2', 'Purples'),
    ('Oligodendrocytes', 'Oligo_2', 'YlOrBr'),
    ('Microglia', 'Micro', 'BuPu'),
]

# Layer config for profiles
layer_profile = [
    ('L2/3', 'Ext_L23', '#2166ac'),
    ('L2-5', 'Ext_L25', '#4393c3'),
    ('L5-6', 'Ext_L56', '#f4a582'),
    ('L6', 'Ext_L6', '#d6604d'),
    ('L6b', 'Ext_L6B', '#b2182b'),
]

# Dominant layer data
layer_cols = ['Ext_L23', 'Ext_L25', 'Ext_L56', 'Ext_L6', 'Ext_L6B']
layer_names = ['L2/3', 'L2-5', 'L5-6', 'L6', 'L6b']
layer_colors = ['#2166ac', '#4393c3', '#f4a582', '#d6604d', '#b2182b']
layer_indices = [get_idx(col) for col in layer_cols]
layer_data = np.column_stack([props[:, idx] for idx in layer_indices if idx is not None])
dominant_layer = np.argmax(layer_data, axis=1)
total_cortical = layer_data.sum(axis=1)

# Line profile data
cortex_mask = total_cortical > 0.05
x_center = np.median(coords[cortex_mask, 1])
near_line = np.abs(coords[:, 1] - x_center) < 150
line_coords = coords[near_line]
line_props = props[near_line]
sort_idx = np.argsort(line_coords[:, 0])
y_sorted = line_coords[sort_idx, 0]
y_normalized = (y_sorted - y_sorted.min()) / (y_sorted.max() - y_sorted.min()) * 100

# =============================================================================
# Create Figure: 5 columns, 3 rows
# =============================================================================
fig = plt.figure(figsize=(16, 12))

gs = GridSpec(3, 5, figure=fig,
              height_ratios=[1.0, 0.75, 1.0],
              width_ratios=[1, 1, 1, 1, 1],
              hspace=0.25, wspace=0.15,
              left=0.03, right=0.97, top=0.93, bottom=0.03)

# -----------------------------------------------------------------------------
# Row 1: Cortical Layers A-D + Dominant Layer E
# -----------------------------------------------------------------------------
for i, (name, col, color) in enumerate(layer_spatial):
    ax = fig.add_subplot(gs[0, i])
    idx = get_idx(col)
    values = props[:, idx]
    vmax = np.percentile(values, 98)
    cmap = LinearSegmentedColormap.from_list('custom', ['#f7f7f7', color])

    sc = ax.scatter(coords[:, 1], -coords[:, 0], c=values, cmap=cmap,
                    s=8, alpha=0.9, vmin=0, vmax=max(vmax, 0.01))
    cbar = plt.colorbar(sc, ax=ax, shrink=0.6, pad=0.01, aspect=15)
    cbar.ax.tick_params(labelsize=6)

    ax.set_title(f'{name}', fontsize=10, fontweight='bold', color=color)
    ax.axis('off')
    ax.set_aspect('equal')
    ax.text(-0.02, 1.02, chr(97+i), transform=ax.transAxes, fontsize=12,
            fontweight='bold', va='bottom', ha='right')

# Cortical Depth Index (E) - weighted average depth
ax = fig.add_subplot(gs[0, 4])
depth_weights = np.array([0, 0.25, 0.5, 0.75, 1.0])  # L2/3 -> L6b (surface to deep)
depth_index = np.sum(layer_data * depth_weights, axis=1) / (total_cortical + 1e-10)
# Mask non-cortical regions
cortex_threshold = 0.02
depth_index_masked = np.where(total_cortical > cortex_threshold, depth_index, np.nan)

sc = ax.scatter(coords[:, 1], -coords[:, 0], c=depth_index_masked,
                cmap='RdYlBu_r', s=8, alpha=0.9, vmin=0, vmax=1)
# Grey background for non-cortical
non_cortex = total_cortical <= cortex_threshold
ax.scatter(coords[non_cortex, 1], -coords[non_cortex, 0], c='#e0e0e0', s=8, alpha=0.5)
cbar = plt.colorbar(sc, ax=ax, shrink=0.6, pad=0.01, aspect=15)
cbar.ax.tick_params(labelsize=6)
cbar.set_label('Depth', fontsize=7)
ax.set_title('Cortical Depth', fontsize=10, fontweight='bold')
ax.axis('off')
ax.set_aspect('equal')
ax.text(-0.02, 1.02, 'e', transform=ax.transAxes, fontsize=12, fontweight='bold', va='bottom', ha='right')

# -----------------------------------------------------------------------------
# Row 2: Line Profile F + Stacked Area G
# -----------------------------------------------------------------------------
ax = fig.add_subplot(gs[1, 0:3])
for name, col, color in layer_profile:
    idx = get_idx(col)
    if idx is None: continue
    values = line_props[sort_idx, idx]
    values_smooth = gaussian_filter1d(values, sigma=3)
    ax.plot(y_normalized, values_smooth, label=name, color=color, linewidth=2)

ax.set_xlabel('Cortical Depth (% from surface)', fontsize=9)
ax.set_ylabel('Proportion', fontsize=9)
ax.set_title('Layer Proportions Along Cortical Depth', fontsize=10, fontweight='bold')
ax.legend(loc='upper right', fontsize=7, framealpha=0.9)
ax.set_xlim(0, 100)
ax.set_ylim(0, None)
ax.grid(True, alpha=0.3)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(labelsize=8)
ax.text(-0.05, 1.02, 'f', transform=ax.transAxes, fontsize=12, fontweight='bold', va='bottom')

ax = fig.add_subplot(gs[1, 3:5])
stack_data, labels, colors_stack = [], [], []
for name, col, color in layer_profile:
    idx = get_idx(col)
    if idx is None: continue
    values = line_props[sort_idx, idx]
    values_smooth = gaussian_filter1d(values, sigma=3)
    stack_data.append(values_smooth)
    labels.append(name)
    colors_stack.append(color)

ax.stackplot(y_normalized, np.array(stack_data), labels=labels, colors=colors_stack, alpha=0.85)
ax.set_xlabel('Cortical Depth (% from surface)', fontsize=9)
ax.set_ylabel('Cumulative Proportion', fontsize=9)
ax.set_title('Laminar Composition', fontsize=10, fontweight='bold')
ax.legend(loc='upper right', fontsize=7, framealpha=0.9)
ax.set_xlim(0, 100)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(labelsize=8)
ax.text(-0.05, 1.02, 'g', transform=ax.transAxes, fontsize=12, fontweight='bold', va='bottom')

# -----------------------------------------------------------------------------
# Row 3: Context Maps H-L - distinct colormaps per cell type
# -----------------------------------------------------------------------------
for i, (name, col, cmap_name) in enumerate(context_maps):
    ax = fig.add_subplot(gs[2, i])
    idx = get_idx(col)
    if idx is None:
        ax.axis('off')
        continue

    values = props[:, idx]
    vmax = np.percentile(values, 98)

    sc = ax.scatter(coords[:, 1], -coords[:, 0], c=values, cmap=cmap_name,
                    s=8, alpha=0.9, vmin=0, vmax=max(vmax, 0.01))
    cbar = plt.colorbar(sc, ax=ax, shrink=0.6, pad=0.01, aspect=15)
    cbar.ax.tick_params(labelsize=6)

    ax.set_title(f'{name}', fontsize=10, fontweight='bold')
    ax.axis('off')
    ax.set_aspect('equal')
    ax.text(-0.02, 1.02, chr(104+i), transform=ax.transAxes, fontsize=12,
            fontweight='bold', va='bottom', ha='right')

# Main title
fig.suptitle('FlashDeconv Accurately Reconstructs Cortical Laminar Organization\n'
             'Mouse Brain Sagittal Section (2,987 spots, 59 cell types, 1.3 seconds)',
             fontsize=13, fontweight='bold', y=0.98)

# Save
plt.savefig('validation/cortex_main_figure.png', dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig('paper/figures/cortex_main_figure.png', dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig('paper/figures/cortex_main_figure.pdf', bbox_inches='tight', facecolor='white')
plt.close()

print("✅ Saved with distinct colormaps for negative controls (H-L)")
