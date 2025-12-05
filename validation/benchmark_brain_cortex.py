"""
Level 2 (v3): Mouse Brain with Cell2location PAIRED Visium Data
----------------------------------------------------------------
Uses the actual paired Visium data (ST8059048) from Cell2location's paper
which was generated from adjacent tissue sections to the scRNA-seq reference.
"""
import os
import sys
import time
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from scipy.sparse import issparse

sys.path.insert(0, '/Users/apple/Research/FlashDeconv')
from flashdeconv import FlashDeconv


def load_paired_data():
    """Load Cell2location's paired Visium and scRNA-seq data."""
    print("="*70)
    print("Loading Cell2location PAIRED Mouse Brain Data")
    print("="*70)

    # Load scRNA-seq reference (already downloaded)
    print("\n1. Loading scRNA-seq reference...")
    sc_adata = sc.read_h5ad("validation/mouse_brain/scrna_reference.h5ad")
    anno = pd.read_csv("validation/mouse_brain/cell_annotation.csv", index_col=0)

    # Subset to cells with annotation
    common_cells = list(set(sc_adata.obs_names) & set(anno.index))
    sc_adata = sc_adata[common_cells].copy()
    sc_adata.obs['cell_type'] = anno.loc[sc_adata.obs_names, 'annotation_1'].values
    print(f"   scRNA: {sc_adata.shape}, {sc_adata.obs['cell_type'].nunique()} cell types")

    # Convert Ensembl to Symbol
    print("\n2. Converting gene IDs to symbols...")
    symbols = sc_adata.var['SYMBOL'].astype(str).values
    has_symbol = (symbols != '') & (symbols != 'nan') & (symbols != 'None')
    sc_adata = sc_adata[:, has_symbol].copy()

    # Make unique names
    new_names = symbols[has_symbol]
    from collections import Counter
    name_counts = Counter()
    unique_names = []
    for name in new_names:
        if name_counts[name] > 0:
            unique_names.append(f"{name}-{name_counts[name]}")
        else:
            unique_names.append(name)
        name_counts[name] += 1
    sc_adata.var_names = pd.Index(unique_names)
    print(f"   After symbol conversion: {sc_adata.shape}")

    # Load PAIRED Visium spatial data (ST8059048 - section 48)
    print("\n3. Loading PAIRED Visium spatial data (ST8059048)...")

    # Read 10X h5 file
    sp_adata = sc.read_10x_h5("validation/mouse_brain/C2L/ST/48/ST8059048_filtered_feature_bc_matrix.h5")
    sp_adata.var_names_make_unique()

    # Load spatial coordinates
    coords_df = pd.read_csv(
        "validation/mouse_brain/C2L/ST/48/spatial/tissue_positions_list.csv",
        header=None,
        index_col=0
    )
    coords_df.columns = ['in_tissue', 'array_row', 'array_col', 'pxl_row', 'pxl_col']

    # Filter to in_tissue spots
    in_tissue = coords_df[coords_df['in_tissue'] == 1].index
    common_spots = list(set(sp_adata.obs_names) & set(in_tissue))
    sp_adata = sp_adata[common_spots].copy()

    # Add spatial coordinates
    sp_adata.obsm['spatial'] = coords_df.loc[sp_adata.obs_names, ['pxl_row', 'pxl_col']].values
    print(f"   Spatial: {sp_adata.shape}, {len(common_spots)} spots in tissue")

    # Find common genes
    common_genes = list(set(sc_adata.var_names) & set(sp_adata.var_names))
    print(f"   Common genes: {len(common_genes)}")

    # Subset to common genes
    sc_adata = sc_adata[:, common_genes].copy()
    sp_adata = sp_adata[:, common_genes].copy()

    return sc_adata, sp_adata, common_genes


def create_reference_signatures(sc_adata):
    """Create cell type reference signatures from scRNA-seq data."""
    print("\n4. Creating reference signatures...")

    cell_types = sc_adata.obs['cell_type'].unique()
    n_types = len(cell_types)
    n_genes = sc_adata.n_vars

    # Average expression per cell type
    X_ref = np.zeros((n_types, n_genes))

    for i, ct in enumerate(cell_types):
        mask = (sc_adata.obs['cell_type'] == ct).values
        if issparse(sc_adata.X):
            X_ref[i] = sc_adata.X[mask].toarray().mean(axis=0)
        else:
            X_ref[i] = sc_adata.X[mask].mean(axis=0)

    print(f"   Reference shape: {X_ref.shape}")
    print(f"   Cell types: {n_types}")

    return X_ref, list(cell_types)


def run_flashdeconv(sp_adata, X_ref, cell_types):
    """Run FlashDeconv deconvolution."""
    print("\n" + "="*70)
    print("Running FlashDeconv")
    print("="*70)

    # Prepare spatial data
    if issparse(sp_adata.X):
        Y = sp_adata.X.toarray()
    else:
        Y = np.array(sp_adata.X)

    coords = sp_adata.obsm['spatial']

    # Run FlashDeconv
    t0 = time.time()
    model = FlashDeconv(
        sketch_dim=512,
        lambda_spatial=5000.0,
        rho_sparsity=0.01,
        preprocess="log_cpm",
        n_hvg=2000,
        max_iter=100,
        verbose=True,
        random_state=42,
    )

    proportions = model.fit_transform(Y, X_ref, coords)
    elapsed = time.time() - t0

    print(f"\nCompleted in {elapsed:.1f}s")

    # Create DataFrame
    prop_df = pd.DataFrame(
        proportions,
        index=sp_adata.obs_names,
        columns=cell_types
    )

    return prop_df, elapsed


def evaluate_with_markers(sp_adata, prop_df, sc_adata):
    """Evaluate predictions using marker gene correlation."""
    print("\n" + "="*70)
    print("Evaluating Predictions (Marker Gene Correlation)")
    print("="*70)

    # Normalize for marker finding (on a copy)
    sc_adata_eval = sc_adata.copy()
    sc.pp.normalize_total(sc_adata_eval, target_sum=1e4)
    sc.pp.log1p(sc_adata_eval)

    try:
        sc.tl.rank_genes_groups(sc_adata_eval, 'cell_type', method='wilcoxon', n_genes=20)
        markers_available = True
    except Exception as e:
        print(f"Marker finding failed: {e}")
        markers_available = False

    results = []

    if markers_available:
        for ct in prop_df.columns:
            try:
                # Get top markers for this cell type
                markers = sc_adata_eval.uns['rank_genes_groups']['names'][ct][:10]
                markers = [g for g in markers if g in sp_adata.var_names]

                if len(markers) < 2:
                    continue

                # Get marker expression in spatial data
                marker_idx = [list(sp_adata.var_names).index(g) for g in markers]
                if issparse(sp_adata.X):
                    marker_expr = sp_adata.X[:, marker_idx].toarray().mean(axis=1)
                else:
                    marker_expr = sp_adata.X[:, marker_idx].mean(axis=1)

                # Correlation with prediction
                pred = prop_df[ct].values
                r, p = pearsonr(marker_expr.flatten(), pred.flatten())

                results.append({
                    'cell_type': ct,
                    'correlation': r,
                    'p_value': p,
                    'n_markers': len(markers)
                })
            except Exception as e:
                continue

    if results:
        results_df = pd.DataFrame(results).sort_values('correlation', ascending=False)
        print(f"\nTop 15 cell types by marker correlation:")
        print(results_df.head(15).to_string(index=False))

        mean_r = results_df['correlation'].mean()
        print(f"\nMean correlation: {mean_r:.3f}")
        return results_df
    else:
        print("Could not compute marker correlations")
        return None


def plot_results(sp_adata, prop_df, results_df, save_prefix):
    """Create comprehensive visualization."""
    coords = sp_adata.obsm['spatial']

    # Plot 1: Top cell types by proportion
    top_types = prop_df.mean().sort_values(ascending=False).head(12).index

    fig, axes = plt.subplots(3, 4, figsize=(16, 12))
    axes = axes.flatten()

    for i, ct in enumerate(top_types):
        ax = axes[i]
        sc_plot = ax.scatter(
            coords[:, 0], coords[:, 1],
            c=prop_df[ct].values,
            cmap='Reds',
            s=8,
            alpha=0.8
        )
        plt.colorbar(sc_plot, ax=ax, shrink=0.8)

        if results_df is not None and ct in results_df['cell_type'].values:
            r = results_df[results_df['cell_type']==ct]['correlation'].values[0]
            ax.set_title(f'{ct[:20]}\\n(r={r:.2f})', fontsize=9)
        else:
            ax.set_title(f'{ct[:20]}', fontsize=9)
        ax.axis('off')

    plt.suptitle('Level 2 (v3): FlashDeconv on Cell2location PAIRED Data\\n(ST8059048 with matched scRNA-seq reference)',
                 fontsize=12, fontweight='bold')
    plt.tight_layout()
    plt.savefig(f'{save_prefix}_spatial.png', dpi=150, bbox_inches='tight')
    plt.close()
    print(f"\nSaved: {save_prefix}_spatial.png")

    # Plot 2: Detailed analysis
    if results_df is not None:
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))

        # Histogram of correlations
        ax = axes[0, 0]
        valid_r = results_df['correlation'].dropna()
        ax.hist(valid_r, bins=20, color='steelblue', edgecolor='black', alpha=0.7)
        ax.axvline(valid_r.mean(), color='orange', linestyle='--', linewidth=2,
                   label=f'Mean ({valid_r.mean():.2f})')
        ax.axvline(0.5, color='green', linestyle=':', linewidth=2, label='Good (r>0.5)')
        ax.axvline(0.3, color='red', linestyle=':', linewidth=2, label='Fair (r>0.3)')
        ax.set_xlabel('Marker Correlation (r)')
        ax.set_ylabel('Number of Cell Types')
        ax.set_title(f'Distribution of Marker Correlations\\n({len(valid_r)} cell types)')
        ax.legend()

        # Top 20 cell types bar chart
        ax = axes[0, 1]
        top20 = results_df.head(20)
        colors = ['green' if r > 0.5 else 'orange' if r > 0.3 else 'red'
                  for r in top20['correlation']]
        ax.barh(range(len(top20)), top20['correlation'].values, color=colors)
        ax.set_yticks(range(len(top20)))
        ax.set_yticklabels(top20['cell_type'].values, fontsize=8)
        ax.axvline(0.5, color='green', linestyle=':', alpha=0.5)
        ax.axvline(0.3, color='orange', linestyle=':', alpha=0.5)
        ax.set_xlabel('Marker Correlation (r)')
        ax.set_title('Top 20 Cell Types by Correlation')
        ax.invert_yaxis()

        # Cell type category analysis
        ax = axes[1, 0]
        results_df_copy = results_df.copy()

        def get_category(ct):
            ct_lower = ct.lower()
            if 'ext' in ct_lower or 'excit' in ct_lower:
                return 'Excitatory'
            elif 'inh' in ct_lower:
                return 'Inhibitory'
            elif 'astro' in ct_lower:
                return 'Astrocytes'
            elif 'oligo' in ct_lower or 'opc' in ct_lower:
                return 'Oligodendrocytes'
            else:
                return 'Other'

        results_df_copy['category'] = results_df_copy['cell_type'].apply(get_category)
        cat_stats = results_df_copy.groupby('category').agg({
            'correlation': ['mean', 'count']
        }).round(3)
        cat_stats.columns = ['Mean r', 'Count']
        cat_stats = cat_stats.sort_values('Mean r', ascending=False)

        colors_cat = {'Excitatory': 'blue', 'Inhibitory': 'red',
                      'Astrocytes': 'purple', 'Oligodendrocytes': 'green', 'Other': 'gray'}
        bar_colors = [colors_cat.get(c, 'gray') for c in cat_stats.index]

        bars = ax.bar(range(len(cat_stats)), cat_stats['Mean r'].values, color=bar_colors)
        ax.set_xticks(range(len(cat_stats)))
        ax.set_xticklabels(cat_stats.index, rotation=45, ha='right')
        ax.set_ylabel('Mean Correlation')
        ax.set_title('Mean Correlation by Cell Type Category')
        ax.axhline(0.3, color='orange', linestyle=':', alpha=0.5)

        # Add count labels
        for i, (idx, row) in enumerate(cat_stats.iterrows()):
            ax.annotate(f'n={int(row["Count"])}', xy=(i, row['Mean r']),
                        ha='center', va='bottom', fontsize=9)

        # Correlation vs Abundance
        ax = axes[1, 1]
        mean_props = prop_df.mean()
        for _, row in results_df.iterrows():
            ct = row['cell_type']
            if ct in mean_props.index:
                cat = get_category(ct)
                ax.scatter(mean_props[ct], row['correlation'],
                           c=colors_cat.get(cat, 'gray'), alpha=0.6, s=50)
        ax.set_xlabel('Mean Proportion')
        ax.set_ylabel('Marker Correlation')
        ax.set_title('Correlation vs Abundance\\n(r = nan)')

        plt.tight_layout()
        plt.savefig(f'{save_prefix}_analysis.png', dpi=150, bbox_inches='tight')
        plt.close()
        print(f"Saved: {save_prefix}_analysis.png")


def main():
    # Load PAIRED data
    sc_adata, sp_adata, common_genes = load_paired_data()

    # Create reference
    X_ref, cell_types = create_reference_signatures(sc_adata)

    # Run FlashDeconv
    prop_df, elapsed = run_flashdeconv(sp_adata, X_ref, cell_types)

    # Evaluate
    results_df = evaluate_with_markers(sp_adata, prop_df, sc_adata)

    # Plot
    plot_results(sp_adata, prop_df, results_df, 'validation/level2_v3_paired')

    # Save results
    prop_df.to_csv('validation/level2_v3_proportions.csv')
    if results_df is not None:
        results_df.to_csv('validation/level2_v3_correlations.csv', index=False)

    # Verdict
    print("\n" + "="*70)
    print("LEVEL 2 (v3) VALIDATION RESULT - PAIRED DATA")
    print("="*70)

    if results_df is not None:
        mean_r = results_df['correlation'].mean()
        top_r = results_df.head(10)['correlation'].mean()
        good_count = (results_df['correlation'] > 0.5).sum()
        fair_count = ((results_df['correlation'] > 0.3) & (results_df['correlation'] <= 0.5)).sum()

        print(f"Mean marker correlation (all types): {mean_r:.3f}")
        print(f"Mean marker correlation (top 10):    {top_r:.3f}")
        print(f"Cell types with r > 0.5: {good_count}/{len(results_df)}")
        print(f"Cell types with r > 0.3: {good_count + fair_count}/{len(results_df)}")
        print(f"Number of cell types evaluated: {len(results_df)}")
        print(f"Time: {elapsed:.1f}s")

        if top_r > 0.5:
            print("\nSTATUS: PASSED (top 10 cell types r > 0.5)")
        elif top_r > 0.4:
            print("\nSTATUS: GOOD")
        elif top_r > 0.3:
            print("\nSTATUS: MARGINAL")
        else:
            print("\nSTATUS: NEEDS IMPROVEMENT")

    return prop_df, results_df


if __name__ == "__main__":
    main()
