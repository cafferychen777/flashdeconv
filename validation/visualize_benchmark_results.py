"""
Visualize FlashDeconv benchmark results on Spotless Silver Standards.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Set style
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['font.size'] = 11
plt.rcParams['axes.labelsize'] = 12
plt.rcParams['axes.titlesize'] = 13
plt.rcParams['figure.dpi'] = 150


def load_results():
    """Load benchmark results."""
    df = pd.read_csv('/Users/apple/Research/FlashDeconv/validation/fair_silver_benchmark_results.csv')
    return df


def plot_method_comparison(df, save_path):
    """Compare FlashDeconv with published methods."""
    # Published results from Spotless paper (63 silver standards average)
    published = pd.DataFrame([
        {'Method': 'rctd', 'Pearson': 0.9046, 'RMSE': 0.0613},
        {'Method': 'cell2location', 'Pearson': 0.8953, 'RMSE': 0.0603},
        {'Method': 'music', 'Pearson': 0.8901, 'RMSE': 0.0775},
        {'Method': 'spatialdwls', 'Pearson': 0.8751, 'RMSE': 0.0654},
        {'Method': 'nnls', 'Pearson': 0.8133, 'RMSE': 0.0868},
        {'Method': 'FlashDeconv', 'Pearson': df['pearson'].mean(), 'RMSE': df['rmse'].mean()},
    ])

    # Sort by Pearson
    published = published.sort_values('Pearson', ascending=True)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Pearson correlation
    colors = ['#2ecc71' if m == 'FlashDeconv' else '#3498db' for m in published['Method']]
    ax1 = axes[0]
    bars1 = ax1.barh(published['Method'], published['Pearson'], color=colors, edgecolor='white', linewidth=1.5)
    ax1.set_xlabel('Pearson Correlation')
    ax1.set_title('Method Comparison: Pearson Correlation')
    ax1.set_xlim(0.75, 1.0)

    # Add value labels
    for bar, val in zip(bars1, published['Pearson']):
        ax1.text(val + 0.005, bar.get_y() + bar.get_height()/2, f'{val:.4f}',
                va='center', fontsize=10)

    # RMSE (lower is better)
    published_rmse = published.sort_values('RMSE', ascending=False)
    colors_rmse = ['#2ecc71' if m == 'FlashDeconv' else '#e74c3c' for m in published_rmse['Method']]
    ax2 = axes[1]
    bars2 = ax2.barh(published_rmse['Method'], published_rmse['RMSE'], color=colors_rmse, edgecolor='white', linewidth=1.5)
    ax2.set_xlabel('RMSE (lower is better)')
    ax2.set_title('Method Comparison: RMSE')
    ax2.set_xlim(0, 0.12)

    # Add value labels
    for bar, val in zip(bars2, published_rmse['RMSE']):
        ax2.text(val + 0.002, bar.get_y() + bar.get_height()/2, f'{val:.4f}',
                va='center', fontsize=10)

    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"Saved: {save_path}")


def plot_dataset_performance(df, save_path):
    """Box plot of performance by dataset."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # Dataset name mapping for better labels
    name_map = {
        'brain_cortex': 'Brain Cortex',
        'cerebellum_cell': 'Cerebellum (Cell)',
        'cerebellum_nucleus': 'Cerebellum (Nucleus)',
        'hippocampus': 'Hippocampus',
        'kidney': 'Kidney',
        'scc_p5': 'SCC P5'
    }
    df['Dataset'] = df['dataset_name'].map(name_map)

    # Order by mean Pearson
    order = df.groupby('Dataset')['pearson'].mean().sort_values(ascending=False).index

    # Pearson correlation
    sns.boxplot(data=df, x='Dataset', y='pearson', ax=axes[0], order=order, palette='viridis')
    axes[0].set_ylabel('Pearson Correlation')
    axes[0].set_title('Pearson by Dataset')
    axes[0].tick_params(axis='x', rotation=45)
    axes[0].axhline(y=0.9046, color='red', linestyle='--', alpha=0.7, label='rctd (0.9046)')
    axes[0].legend(loc='lower right')

    # RMSE
    sns.boxplot(data=df, x='Dataset', y='rmse', ax=axes[1], order=order, palette='viridis')
    axes[1].set_ylabel('RMSE')
    axes[1].set_title('RMSE by Dataset')
    axes[1].tick_params(axis='x', rotation=45)
    axes[1].axhline(y=0.0613, color='red', linestyle='--', alpha=0.7, label='rctd (0.0613)')
    axes[1].legend(loc='upper right')

    # JSD
    sns.boxplot(data=df, x='Dataset', y='jsd', ax=axes[2], order=order, palette='viridis')
    axes[2].set_ylabel('Jensen-Shannon Divergence')
    axes[2].set_title('JSD by Dataset')
    axes[2].tick_params(axis='x', rotation=45)

    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"Saved: {save_path}")


def plot_heatmap(df, save_path):
    """Heatmap of Pearson correlation by dataset and pattern."""
    # Pivot table
    pivot = df.pivot_table(values='pearson', index='dataset_name', columns='pattern', aggfunc='mean')

    # Rename for better labels
    name_map = {
        'brain_cortex': 'Brain Cortex',
        'cerebellum_cell': 'Cerebellum (Cell)',
        'cerebellum_nucleus': 'Cerebellum (Nucleus)',
        'hippocampus': 'Hippocampus',
        'kidney': 'Kidney',
        'scc_p5': 'SCC P5'
    }
    pivot.index = pivot.index.map(name_map)

    fig, ax = plt.subplots(figsize=(12, 6))

    sns.heatmap(pivot, annot=True, fmt='.3f', cmap='RdYlGn', vmin=0.8, vmax=1.0,
                linewidths=0.5, ax=ax, cbar_kws={'label': 'Pearson Correlation'})

    ax.set_xlabel('Abundance Pattern')
    ax.set_ylabel('Dataset')
    ax.set_title('FlashDeconv Performance: Pearson Correlation Across Datasets and Patterns')

    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"Saved: {save_path}")


def plot_summary_stats(df, save_path):
    """Summary statistics visualization."""
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # 1. Overall distribution of Pearson
    ax1 = axes[0, 0]
    ax1.hist(df['pearson'], bins=20, color='#3498db', edgecolor='white', alpha=0.8)
    ax1.axvline(df['pearson'].mean(), color='red', linestyle='--', linewidth=2, label=f'Mean: {df["pearson"].mean():.4f}')
    ax1.axvline(0.9046, color='orange', linestyle='--', linewidth=2, label='rctd: 0.9046')
    ax1.set_xlabel('Pearson Correlation')
    ax1.set_ylabel('Count')
    ax1.set_title('Distribution of Pearson Correlation (n=56)')
    ax1.legend()

    # 2. Scatter: Pearson vs RMSE
    ax2 = axes[0, 1]
    scatter = ax2.scatter(df['pearson'], df['rmse'], c=df['dataset_id'], cmap='tab10', s=60, alpha=0.7)
    ax2.set_xlabel('Pearson Correlation')
    ax2.set_ylabel('RMSE')
    ax2.set_title('Pearson vs RMSE')

    # Add colorbar legend
    cbar = plt.colorbar(scatter, ax=ax2)
    cbar.set_label('Dataset ID')

    # 3. Runtime by dataset
    ax3 = axes[1, 0]
    runtime_by_ds = df.groupby('dataset_name')['time'].mean().sort_values()
    name_map = {
        'brain_cortex': 'Brain Cortex',
        'cerebellum_cell': 'Cerebellum (Cell)',
        'cerebellum_nucleus': 'Cerebellum (Nucleus)',
        'hippocampus': 'Hippocampus',
        'kidney': 'Kidney',
        'scc_p5': 'SCC P5'
    }
    runtime_by_ds.index = runtime_by_ds.index.map(name_map)
    ax3.barh(runtime_by_ds.index, runtime_by_ds.values, color='#9b59b6', edgecolor='white')
    ax3.set_xlabel('Runtime (seconds)')
    ax3.set_title('Average Runtime by Dataset')
    for i, (name, val) in enumerate(runtime_by_ds.items()):
        ax3.text(val + 0.01, i, f'{val:.2f}s', va='center')

    # 4. Summary table as text
    ax4 = axes[1, 1]
    ax4.axis('off')

    summary_text = f"""
    BENCHMARK SUMMARY
    {'='*40}

    Total Tests: {len(df)} dataset-pattern combinations

    Overall Performance:
      Pearson:  {df['pearson'].mean():.4f} ± {df['pearson'].std():.4f}
      RMSE:     {df['rmse'].mean():.4f} ± {df['rmse'].std():.4f}
      JSD:      {df['jsd'].mean():.4f} ± {df['jsd'].std():.4f}
      AUPR:     {df['aupr'].mean():.4f} ± {df['aupr'].std():.4f}

    Best Dataset:  {df.groupby('dataset_name')['pearson'].mean().idxmax()}
                   ({df.groupby('dataset_name')['pearson'].mean().max():.4f})

    Worst Dataset: {df.groupby('dataset_name')['pearson'].mean().idxmin()}
                   ({df.groupby('dataset_name')['pearson'].mean().min():.4f})

    Comparison with Published Methods:
      FlashDeconv: {df['pearson'].mean():.4f} (Rank #1)
      rctd:        0.9046
      cell2loc:    0.8953
      music:       0.8901
    """

    ax4.text(0.1, 0.9, summary_text, transform=ax4.transAxes, fontsize=11,
             verticalalignment='top', fontfamily='monospace',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"Saved: {save_path}")


def create_combined_figure(df, save_path):
    """Create a publication-ready combined figure."""
    fig = plt.figure(figsize=(14, 10))

    # Create grid
    gs = fig.add_gridspec(2, 2, hspace=0.3, wspace=0.3)

    # A: Method comparison (Pearson)
    ax1 = fig.add_subplot(gs[0, 0])
    published = pd.DataFrame([
        {'Method': 'rctd', 'Pearson': 0.9046},
        {'Method': 'cell2location', 'Pearson': 0.8953},
        {'Method': 'music', 'Pearson': 0.8901},
        {'Method': 'spatialdwls', 'Pearson': 0.8751},
        {'Method': 'nnls', 'Pearson': 0.8133},
        {'Method': 'FlashDeconv', 'Pearson': df['pearson'].mean()},
    ]).sort_values('Pearson', ascending=True)

    colors = ['#2ecc71' if m == 'FlashDeconv' else '#95a5a6' for m in published['Method']]
    bars = ax1.barh(published['Method'], published['Pearson'], color=colors, edgecolor='white', linewidth=1.5)
    ax1.set_xlabel('Pearson Correlation')
    ax1.set_title('A. FlashDeconv vs Published Methods', fontweight='bold', loc='left')
    ax1.set_xlim(0.75, 1.0)
    for bar, val in zip(bars, published['Pearson']):
        ax1.text(val + 0.005, bar.get_y() + bar.get_height()/2, f'{val:.3f}', va='center', fontsize=9)

    # B: Box plot by dataset
    ax2 = fig.add_subplot(gs[0, 1])
    name_map = {
        'brain_cortex': 'Brain\nCortex',
        'cerebellum_cell': 'Cerebellum\n(Cell)',
        'cerebellum_nucleus': 'Cerebellum\n(Nucleus)',
        'hippocampus': 'Hippo-\ncampus',
        'kidney': 'Kidney',
        'scc_p5': 'SCC P5'
    }
    df_plot = df.copy()
    df_plot['Dataset'] = df_plot['dataset_name'].map(name_map)
    order = df_plot.groupby('Dataset')['pearson'].mean().sort_values(ascending=False).index

    sns.boxplot(data=df_plot, x='Dataset', y='pearson', ax=ax2, order=order, palette='viridis')
    ax2.axhline(y=0.9046, color='red', linestyle='--', alpha=0.7, label='rctd (0.9046)')
    ax2.set_ylabel('Pearson Correlation')
    ax2.set_xlabel('')
    ax2.set_title('B. Performance by Dataset', fontweight='bold', loc='left')
    ax2.legend(loc='lower right', fontsize=9)
    ax2.tick_params(axis='x', rotation=0)

    # C: Heatmap
    ax3 = fig.add_subplot(gs[1, 0])
    pivot = df.pivot_table(values='pearson', index='dataset_name', columns='pattern', aggfunc='mean')
    name_map_full = {
        'brain_cortex': 'Brain Cortex',
        'cerebellum_cell': 'Cerebellum (Cell)',
        'cerebellum_nucleus': 'Cerebellum (Nucleus)',
        'hippocampus': 'Hippocampus',
        'kidney': 'Kidney',
        'scc_p5': 'SCC P5'
    }
    pivot.index = pivot.index.map(name_map_full)

    sns.heatmap(pivot, annot=True, fmt='.2f', cmap='RdYlGn', vmin=0.8, vmax=1.0,
                linewidths=0.5, ax=ax3, cbar_kws={'label': 'Pearson', 'shrink': 0.8},
                annot_kws={'size': 8})
    ax3.set_xlabel('Abundance Pattern')
    ax3.set_ylabel('')
    ax3.set_title('C. Performance Heatmap', fontweight='bold', loc='left')

    # D: Distribution
    ax4 = fig.add_subplot(gs[1, 1])
    ax4.hist(df['pearson'], bins=15, color='#3498db', edgecolor='white', alpha=0.8, density=True)
    ax4.axvline(df['pearson'].mean(), color='#2ecc71', linestyle='-', linewidth=2.5,
                label=f'FlashDeconv mean: {df["pearson"].mean():.3f}')
    ax4.axvline(0.9046, color='red', linestyle='--', linewidth=2, label='rctd: 0.9046')
    ax4.set_xlabel('Pearson Correlation')
    ax4.set_ylabel('Density')
    ax4.set_title('D. Distribution of Results (n=56)', fontweight='bold', loc='left')
    ax4.legend(loc='upper left', fontsize=9)

    plt.savefig(save_path, dpi=200, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"Saved: {save_path}")


def main():
    print("="*60)
    print("Visualizing FlashDeconv Benchmark Results")
    print("="*60)

    # Load data
    df = load_results()
    print(f"Loaded {len(df)} results")

    # Output directory
    out_dir = '/Users/apple/Research/FlashDeconv/validation/figures'
    import os
    os.makedirs(out_dir, exist_ok=True)

    # Generate plots
    print("\nGenerating visualizations...")

    plot_method_comparison(df, f'{out_dir}/method_comparison.png')
    plot_dataset_performance(df, f'{out_dir}/dataset_performance.png')
    plot_heatmap(df, f'{out_dir}/performance_heatmap.png')
    plot_summary_stats(df, f'{out_dir}/summary_stats.png')
    create_combined_figure(df, f'{out_dir}/benchmark_combined.png')

    print(f"\nAll figures saved to: {out_dir}")

    # Print summary
    print("\n" + "="*60)
    print("SUMMARY STATISTICS")
    print("="*60)
    print(f"Total tests: {len(df)}")
    print(f"\nOverall metrics:")
    print(f"  Pearson: {df['pearson'].mean():.4f} ± {df['pearson'].std():.4f}")
    print(f"  RMSE:    {df['rmse'].mean():.4f} ± {df['rmse'].std():.4f}")
    print(f"  JSD:     {df['jsd'].mean():.4f} ± {df['jsd'].std():.4f}")
    print(f"  AUPR:    {df['aupr'].mean():.4f} ± {df['aupr'].std():.4f}")

    print(f"\nPer-dataset Pearson:")
    for ds in df['dataset_name'].unique():
        ds_df = df[df['dataset_name'] == ds]
        print(f"  {ds:25s}: {ds_df['pearson'].mean():.4f} ± {ds_df['pearson'].std():.4f} (n={len(ds_df)})")


if __name__ == '__main__':
    main()
