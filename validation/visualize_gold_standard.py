"""
Visualize Gold Standard benchmark results with correct references.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['font.size'] = 11
plt.rcParams['axes.labelsize'] = 12
plt.rcParams['axes.titlesize'] = 13
plt.rcParams['figure.dpi'] = 150


def load_results():
    """Load Gold Standard results."""
    df = pd.read_csv('/Users/apple/Research/FlashDeconv/validation/gold_benchmark_correct_ref.csv')
    df['tissue'] = df['dataset'].apply(
        lambda x: 'cortex_svz' if 'cortex_svz' in x
        else ('ob' if '_ob_' in x else 'visp')
    )
    df['source'] = df['dataset'].apply(
        lambda x: 'Eng2019' if 'Eng2019' in x else 'Wang2018'
    )
    return df


def create_combined_figure(df):
    """Create publication-ready combined figure."""
    fig = plt.figure(figsize=(14, 10))
    gs = fig.add_gridspec(2, 2, hspace=0.35, wspace=0.3)

    # A: Method comparison for Wang2018
    ax1 = fig.add_subplot(gs[0, 0])
    wang_best = df[(df['source'] == 'Wang2018') & (df['preprocess'] == 'pearson')]['pearson'].values[0]

    methods_wang = pd.DataFrame([
        {'Method': 'FlashDeconv', 'Pearson': wang_best},
        {'Method': 'spatialdwls', 'Pearson': 0.712},
        {'Method': 'music', 'Pearson': 0.663},
        {'Method': 'rctd', 'Pearson': 0.568},
        {'Method': 'cell2location', 'Pearson': 0.457},
    ]).sort_values('Pearson', ascending=True)

    colors = ['#2ecc71' if m == 'FlashDeconv' else '#95a5a6' for m in methods_wang['Method']]
    bars = ax1.barh(methods_wang['Method'], methods_wang['Pearson'], color=colors, edgecolor='white', linewidth=1.5)
    ax1.set_xlabel('Pearson Correlation')
    ax1.set_title('A. Wang2018 STARMap (108 spots)', fontweight='bold', loc='left')
    ax1.set_xlim(0.3, 0.85)
    for bar, val in zip(bars, methods_wang['Pearson']):
        ax1.text(val + 0.01, bar.get_y() + bar.get_height()/2, f'{val:.3f}', va='center', fontsize=10)

    # B: Method comparison for Eng2019
    ax2 = fig.add_subplot(gs[0, 1])
    eng_best = df[(df['source'] == 'Eng2019') & (df['preprocess'] == 'raw')].groupby('source')['pearson'].mean().values[0]

    methods_eng = pd.DataFrame([
        {'Method': 'rctd', 'Pearson': 0.720},
        {'Method': 'music', 'Pearson': 0.719},
        {'Method': 'dstg', 'Pearson': 0.692},
        {'Method': 'spatialdwls', 'Pearson': 0.636},
        {'Method': 'FlashDeconv', 'Pearson': eng_best},
    ]).sort_values('Pearson', ascending=True)

    colors = ['#2ecc71' if m == 'FlashDeconv' else '#95a5a6' for m in methods_eng['Method']]
    bars = ax2.barh(methods_eng['Method'], methods_eng['Pearson'], color=colors, edgecolor='white', linewidth=1.5)
    ax2.set_xlabel('Pearson Correlation')
    ax2.set_title('B. Eng2019 seqFISH+ (14 FOVs, 9 spots each)', fontweight='bold', loc='left')
    ax2.set_xlim(0.4, 0.85)
    for bar, val in zip(bars, methods_eng['Pearson']):
        ax2.text(val + 0.01, bar.get_y() + bar.get_height()/2, f'{val:.3f}', va='center', fontsize=10)

    # C: Per-FOV results for Eng2019
    ax3 = fig.add_subplot(gs[1, 0])
    eng_df = df[(df['source'] == 'Eng2019') & (df['preprocess'] == 'raw')]

    # Extract FOV number for sorting
    eng_df = eng_df.copy()
    eng_df['fov'] = eng_df['dataset'].apply(lambda x: int(x.split('fov')[-1]))
    eng_df['tissue_label'] = eng_df['tissue'].apply(lambda x: 'Cortex' if x == 'cortex_svz' else 'OB')

    # Plot
    for tissue, color in [('cortex_svz', '#3498db'), ('ob', '#e74c3c')]:
        tissue_df = eng_df[eng_df['tissue'] == tissue].sort_values('fov')
        label = 'Cortex/SVZ' if tissue == 'cortex_svz' else 'Olfactory Bulb'
        ax3.plot(range(len(tissue_df)), tissue_df['pearson'].values, 'o-', color=color, label=label, markersize=8)

    ax3.axhline(y=0.720, color='gray', linestyle='--', alpha=0.7, label='rctd avg (0.720)')
    ax3.set_xticks(range(7))
    ax3.set_xticklabels([f'FOV{i}' for i in range(7)])
    ax3.set_xlabel('Field of View')
    ax3.set_ylabel('Pearson Correlation')
    ax3.set_title('C. Per-FOV Performance (Eng2019)', fontweight='bold', loc='left')
    ax3.legend(loc='lower right', fontsize=9)
    ax3.set_ylim(-0.1, 1.05)

    # D: Preprocess comparison
    ax4 = fig.add_subplot(gs[1, 1])

    # Summary by preprocess mode
    summary_data = []
    for source in ['Wang2018', 'Eng2019']:
        for prep in ['log_cpm', 'pearson', 'raw']:
            source_df = df[(df['source'] == source) & (df['preprocess'] == prep)]
            summary_data.append({
                'Dataset': source,
                'Preprocess': prep,
                'Pearson': source_df['pearson'].mean()
            })

    summary_df = pd.DataFrame(summary_data)
    summary_pivot = summary_df.pivot(index='Preprocess', columns='Dataset', values='Pearson')

    x = np.arange(len(summary_pivot.index))
    width = 0.35

    bars1 = ax4.bar(x - width/2, summary_pivot['Wang2018'], width, label='Wang2018', color='#3498db')
    bars2 = ax4.bar(x + width/2, summary_pivot['Eng2019'], width, label='Eng2019', color='#e74c3c')

    ax4.set_xlabel('Preprocess Mode')
    ax4.set_ylabel('Pearson Correlation')
    ax4.set_title('D. Performance by Preprocess Mode', fontweight='bold', loc='left')
    ax4.set_xticks(x)
    ax4.set_xticklabels(summary_pivot.index)
    ax4.legend()
    ax4.set_ylim(0, 0.85)

    # Add value labels
    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            ax4.text(bar.get_x() + bar.get_width()/2., height + 0.02,
                    f'{height:.2f}', ha='center', va='bottom', fontsize=9)

    plt.savefig('/Users/apple/Research/FlashDeconv/validation/figures/gold_standard_benchmark.png',
                dpi=200, bbox_inches='tight', facecolor='white')
    plt.close()
    print("Saved: gold_standard_benchmark.png")


def print_summary(df):
    """Print summary statistics."""
    print("="*70)
    print("GOLD STANDARD BENCHMARK SUMMARY (Correct References)")
    print("="*70)

    print("\n--- Wang2018 STARMap ---")
    wang = df[df['source'] == 'Wang2018']
    for prep in ['log_cpm', 'pearson', 'raw']:
        val = wang[wang['preprocess'] == prep]['pearson'].values[0]
        print(f"  {prep:10s}: {val:.4f}")

    best_prep = wang.loc[wang['pearson'].idxmax(), 'preprocess']
    best_val = wang['pearson'].max()
    print(f"  Best: {best_prep} = {best_val:.4f}")
    print(f"  Rank: #1 (> spatialdwls 0.712)")

    print("\n--- Eng2019 seqFISH+ ---")
    eng = df[df['source'] == 'Eng2019']
    for prep in ['log_cpm', 'pearson', 'raw']:
        prep_df = eng[eng['preprocess'] == prep]
        print(f"  {prep:10s}: {prep_df['pearson'].mean():.4f} ± {prep_df['pearson'].std():.4f}")

    best_prep = eng.groupby('preprocess')['pearson'].mean().idxmax()
    best_val = eng[eng['preprocess'] == best_prep]['pearson'].mean()
    print(f"  Best: {best_prep} = {best_val:.4f}")
    print(f"  Rank: #5 (< spatialdwls 0.636)")

    print("\n--- Comparison with Spotless Paper ---")
    print("Wang2018: FlashDeconv #1, Eng2019: FlashDeconv #5")
    print("Note: Eng2019 has very small samples (9 spots/FOV), high variance")


def main():
    import os
    os.makedirs('/Users/apple/Research/FlashDeconv/validation/figures', exist_ok=True)

    df = load_results()
    create_combined_figure(df)
    print_summary(df)


if __name__ == '__main__':
    main()
