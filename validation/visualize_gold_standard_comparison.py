"""
Visualize Gold Standard benchmark results with original Spotless results.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['font.size'] = 10
plt.rcParams['axes.labelsize'] = 11
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['figure.dpi'] = 150


def create_gold_standard_figure():
    """Create comprehensive Gold Standard comparison figure."""

    # Original Spotless results (Pearson correlation)
    seqfish_cortex = {
        'music': 0.7219, 'rctd': 0.7202, 'cell2location': 0.6903, 'nnls': 0.6726,
        'seurat': 0.6137, 'dstg': 0.5450, 'spatialdwls': 0.4257, 'stride': 0.2984,
        'stereoscope': 0.2801, 'spotlight': 0.1483, 'tangram': 0.0141, 'destvi': -0.1200
    }

    seqfish_ob = {
        'stride': 0.8068, 'dstg': 0.7773, 'cell2location': 0.7543, 'music': 0.7481,
        'rctd': 0.7412, 'nnls': 0.7074, 'seurat': 0.6604, 'spotlight': 0.5977,
        'spatialdwls': 0.5058, 'destvi': 0.2895, 'stereoscope': 0.2649, 'tangram': 0.0999
    }

    starmap = {
        'spatialdwls': 0.7116, 'rctd': 0.6792, 'spotlight': 0.6235, 'music': 0.6033,
        'tangram': 0.5733, 'cell2location': 0.5422, 'nnls': 0.3661, 'dstg': 0.2883,
        'seurat': 0.2509, 'stereoscope': 0.1789, 'stride': -0.1033, 'destvi': -0.1012
    }

    # FlashDeconv results
    flashdeconv = {
        'cortex': 0.5809,
        'ob': 0.6185,
        'starmap': 0.7426
    }

    fig = plt.figure(figsize=(16, 10))
    gs = fig.add_gridspec(2, 3, hspace=0.35, wspace=0.3)

    # Color setup
    def get_colors(methods, flash_val):
        sorted_methods = sorted(methods.items(), key=lambda x: x[1], reverse=True)
        all_methods = sorted_methods + [('FlashDeconv', flash_val)]
        all_methods = sorted(all_methods, key=lambda x: x[1], reverse=True)
        colors = []
        for m, v in all_methods:
            if m == 'FlashDeconv':
                colors.append('#2ecc71')  # Green for FlashDeconv
            elif m in ['rctd', 'cell2location', 'music']:
                colors.append('#3498db')  # Blue for top methods
            else:
                colors.append('#95a5a6')  # Gray for others
        return all_methods, colors

    # A: seqFISH+ Cortex
    ax1 = fig.add_subplot(gs[0, 0])
    methods, colors = get_colors(seqfish_cortex, flashdeconv['cortex'])
    methods_sorted = sorted(methods, key=lambda x: x[1], reverse=False)  # Ascending for barh
    colors_sorted = [colors[methods.index(m)] for m in methods_sorted]

    y_pos = range(len(methods_sorted))
    bars = ax1.barh(y_pos, [v for _, v in methods_sorted], color=colors_sorted, edgecolor='white', linewidth=0.5)
    ax1.set_yticks(y_pos)
    ax1.set_yticklabels([m for m, _ in methods_sorted], fontsize=9)
    ax1.set_xlabel('Pearson Correlation')
    ax1.set_title('A. seqFISH+ Cortex\n(7 FOVs, 9 spots each)', fontweight='bold')
    ax1.set_xlim(-0.2, 0.85)
    ax1.axvline(x=0, color='black', linestyle='-', linewidth=0.5, alpha=0.3)

    # Add FlashDeconv rank annotation
    flash_rank = [i for i, (m, _) in enumerate(methods, 1) if m == 'FlashDeconv'][0]
    ax1.text(0.95, 0.05, f'FlashDeconv: #{flash_rank}/13', transform=ax1.transAxes,
             fontsize=10, fontweight='bold', color='#2ecc71', ha='right')

    # B: seqFISH+ Olfactory Bulb
    ax2 = fig.add_subplot(gs[0, 1])
    methods, colors = get_colors(seqfish_ob, flashdeconv['ob'])
    methods_sorted = sorted(methods, key=lambda x: x[1], reverse=False)
    colors_sorted = [colors[methods.index(m)] for m in methods_sorted]

    y_pos = range(len(methods_sorted))
    bars = ax2.barh(y_pos, [v for _, v in methods_sorted], color=colors_sorted, edgecolor='white', linewidth=0.5)
    ax2.set_yticks(y_pos)
    ax2.set_yticklabels([m for m, _ in methods_sorted], fontsize=9)
    ax2.set_xlabel('Pearson Correlation')
    ax2.set_title('B. seqFISH+ Olfactory Bulb\n(7 FOVs, 9 spots each)', fontweight='bold')
    ax2.set_xlim(-0.2, 0.9)
    ax2.axvline(x=0, color='black', linestyle='-', linewidth=0.5, alpha=0.3)

    flash_rank = [i for i, (m, _) in enumerate(methods, 1) if m == 'FlashDeconv'][0]
    ax2.text(0.95, 0.05, f'FlashDeconv: #{flash_rank}/13', transform=ax2.transAxes,
             fontsize=10, fontweight='bold', color='#2ecc71', ha='right')

    # C: STARMap
    ax3 = fig.add_subplot(gs[0, 2])
    methods, colors = get_colors(starmap, flashdeconv['starmap'])
    methods_sorted = sorted(methods, key=lambda x: x[1], reverse=False)
    colors_sorted = [colors[methods.index(m)] for m in methods_sorted]

    y_pos = range(len(methods_sorted))
    bars = ax3.barh(y_pos, [v for _, v in methods_sorted], color=colors_sorted, edgecolor='white', linewidth=0.5)
    ax3.set_yticks(y_pos)
    ax3.set_yticklabels([m for m, _ in methods_sorted], fontsize=9)
    ax3.set_xlabel('Pearson Correlation')
    ax3.set_title('C. STARMap\n(108 spots)', fontweight='bold')
    ax3.set_xlim(-0.2, 0.85)
    ax3.axvline(x=0, color='black', linestyle='-', linewidth=0.5, alpha=0.3)

    flash_rank = [i for i, (m, _) in enumerate(methods, 1) if m == 'FlashDeconv'][0]
    ax3.text(0.95, 0.05, f'FlashDeconv: #{flash_rank}/13', transform=ax3.transAxes,
             fontsize=10, fontweight='bold', color='#2ecc71', ha='right')

    # D: Summary across all 3 Gold Standards
    ax4 = fig.add_subplot(gs[1, 0])

    # Calculate average rank for each method
    all_methods_set = set(seqfish_cortex.keys()) | set(seqfish_ob.keys()) | set(starmap.keys()) | {'FlashDeconv'}

    def get_rank(methods_dict, flash_val, method):
        all_m = list(methods_dict.items()) + [('FlashDeconv', flash_val)]
        all_m = sorted(all_m, key=lambda x: x[1], reverse=True)
        for i, (m, _) in enumerate(all_m, 1):
            if m == method:
                return i
        return len(all_m)

    avg_ranks = {}
    for method in all_methods_set:
        if method == 'FlashDeconv':
            r1 = get_rank(seqfish_cortex, flashdeconv['cortex'], 'FlashDeconv')
            r2 = get_rank(seqfish_ob, flashdeconv['ob'], 'FlashDeconv')
            r3 = get_rank(starmap, flashdeconv['starmap'], 'FlashDeconv')
        else:
            r1 = get_rank(seqfish_cortex, flashdeconv['cortex'], method)
            r2 = get_rank(seqfish_ob, flashdeconv['ob'], method)
            r3 = get_rank(starmap, flashdeconv['starmap'], method)
        avg_ranks[method] = (r1 + r2 + r3) / 3

    # Sort by average rank
    sorted_ranks = sorted(avg_ranks.items(), key=lambda x: x[1])
    methods_names = [m for m, _ in sorted_ranks]
    ranks_vals = [r for _, r in sorted_ranks]

    colors = ['#2ecc71' if m == 'FlashDeconv' else '#3498db' if m in ['rctd', 'cell2location', 'music'] else '#95a5a6'
              for m in methods_names]

    bars = ax4.barh(range(len(methods_names)), ranks_vals, color=colors, edgecolor='white', linewidth=0.5)
    ax4.set_yticks(range(len(methods_names)))
    ax4.set_yticklabels(methods_names, fontsize=9)
    ax4.set_xlabel('Average Rank (lower is better)')
    ax4.set_title('D. Average Rank Across 3 Gold Standards', fontweight='bold')
    ax4.invert_xaxis()

    flash_avg_rank = avg_ranks['FlashDeconv']
    ax4.text(0.05, 0.05, f'FlashDeconv avg rank: {flash_avg_rank:.1f}', transform=ax4.transAxes,
             fontsize=10, fontweight='bold', color='#2ecc71')

    # E: Per-FOV performance for seqFISH+ (FlashDeconv)
    ax5 = fig.add_subplot(gs[1, 1])

    # Load FlashDeconv per-FOV results
    df = pd.read_csv('/Users/apple/Research/FlashDeconv/validation/gold_benchmark_correct_ref.csv')
    df['tissue'] = df['dataset'].apply(lambda x: 'cortex_svz' if 'cortex_svz' in x else ('ob' if '_ob_' in x else 'visp'))

    # Use best preprocess (raw for seqFISH+)
    df_raw = df[df['preprocess'] == 'raw']

    cortex_fovs = df_raw[df_raw['tissue'] == 'cortex_svz'].sort_values('dataset')['pearson'].values
    ob_fovs = df_raw[df_raw['tissue'] == 'ob'].sort_values('dataset')['pearson'].values

    x = range(7)
    ax5.plot(x, cortex_fovs, 'o-', color='#3498db', label='Cortex/SVZ', markersize=8, linewidth=2)
    ax5.plot(x, ob_fovs, 's-', color='#e74c3c', label='Olfactory Bulb', markersize=8, linewidth=2)
    ax5.axhline(y=0.7219, color='#3498db', linestyle='--', alpha=0.5, label='music (cortex best)')
    ax5.axhline(y=0.8068, color='#e74c3c', linestyle='--', alpha=0.5, label='stride (ob best)')
    ax5.set_xticks(x)
    ax5.set_xticklabels([f'FOV{i}' for i in range(7)])
    ax5.set_xlabel('Field of View')
    ax5.set_ylabel('Pearson Correlation')
    ax5.set_title('E. FlashDeconv Per-FOV Performance\n(seqFISH+)', fontweight='bold')
    ax5.legend(loc='lower right', fontsize=8)
    ax5.set_ylim(-0.1, 1.05)

    # F: Summary table
    ax6 = fig.add_subplot(gs[1, 2])
    ax6.axis('off')

    summary_text = """
    GOLD STANDARD BENCHMARK SUMMARY
    ════════════════════════════════════════

    Dataset              FlashDeconv    Rank    Best Method
    ─────────────────────────────────────────────────────────
    STARMap              0.743          #1      -
    seqFISH+ Cortex      0.581          #6      music (0.722)
    seqFISH+ OB          0.619          #8      stride (0.807)

    ════════════════════════════════════════

    Key Findings:
    • STARMap (108 spots): FlashDeconv ranks #1
    • seqFISH+ has only 9 spots per FOV
    • High variance across FOVs (0.08 - 0.95)
    • Performance correlates with sample size

    Legend:
    ■ FlashDeconv    ■ Top methods    ■ Other methods
    """

    ax6.text(0.05, 0.95, summary_text, transform=ax6.transAxes, fontsize=10,
             verticalalignment='top', fontfamily='monospace',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))

    plt.savefig('/Users/apple/Research/FlashDeconv/validation/figures/gold_standard_comparison.png',
                dpi=200, bbox_inches='tight', facecolor='white')
    plt.close()
    print("Saved: gold_standard_comparison.png")


def main():
    import os
    os.makedirs('/Users/apple/Research/FlashDeconv/validation/figures', exist_ok=True)
    create_gold_standard_figure()
    print("Done!")


if __name__ == '__main__':
    main()
