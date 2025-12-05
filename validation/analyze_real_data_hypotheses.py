"""
Validate preprocessing hypotheses on REAL scRNA-seq data
"""

import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
from pathlib import Path

def analyze_real_data():
    """Use real scRNA-seq reference data to validate hypotheses"""

    # Try to load real data
    data_paths = [
        "/Users/apple/Research/FlashDeconv/validation/spotless/gold_ref_1.h5ad",
        "/Users/apple/Research/FlashDeconv/validation/spotless/gold_ref_3.h5ad",
    ]

    adata = None
    for path in data_paths:
        if Path(path).exists():
            print(f"Loading real data from: {path}")
            adata = sc.read_h5ad(path)
            break

    if adata is None:
        print("No real data found, generating realistic simulation...")
        # Generate more realistic data based on real scRNA-seq characteristics
        np.random.seed(42)
        n_cells = 5000
        n_genes = 15000

        # Real scRNA-seq has ~95% zeros
        X_raw = np.zeros((n_cells, n_genes))

        # For each cell, only ~5% genes are expressed
        for i in range(n_cells):
            n_expressed = int(n_genes * 0.05)
            expressed_genes = np.random.choice(n_genes, n_expressed, replace=False)
            # Log-normal distribution for expressed genes
            X_raw[i, expressed_genes] = np.random.lognormal(mean=2, sigma=2, size=n_expressed)

        # Add housekeeping genes (expressed in all cells at high levels)
        housekeeping = np.arange(50)  # First 50 genes are housekeeping
        X_raw[:, housekeeping] = np.random.lognormal(mean=6, sigma=1, size=(n_cells, 50))

        X_raw = X_raw.astype(np.float32)
    else:
        print(f"Data shape: {adata.shape}")
        X_raw = adata.X
        if hasattr(X_raw, 'toarray'):
            X_raw = X_raw.toarray()
        X_raw = X_raw.astype(np.float32)

    print(f"Matrix shape: {X_raw.shape}")
    print(f"Sparsity: {(X_raw == 0).mean()*100:.1f}%")
    print(f"Max value: {X_raw.max():.0f}")
    print()

    # ==========================================
    # Preprocessing functions
    # ==========================================
    def preprocess_logcpm(mat):
        """Log-CPM normalization"""
        lib_size = mat.sum(axis=1, keepdims=True) + 1
        return np.log1p(1e4 * mat / lib_size)

    def preprocess_pearson(mat, theta=100):
        """Uncentered Pearson residuals"""
        # Gene-wise statistics
        mu = mat.mean(axis=0) + 1e-6
        var = mu + (mu**2) / theta
        sigma = np.sqrt(var)
        return mat / sigma

    # Apply preprocessing
    print("Applying preprocessing...")
    X_log = preprocess_logcpm(X_raw)
    X_pearson = preprocess_pearson(X_raw)

    # ==========================================
    # Hypothesis 1: Energy Concentration
    # ==========================================
    def gini(array):
        """Calculate Gini coefficient"""
        array = np.abs(array.flatten())
        array = array + 1e-10
        array = np.sort(array)
        n = len(array)
        index = np.arange(1, n + 1)
        return (np.sum((2 * index - n - 1) * array)) / (n * np.sum(array))

    # Per-gene L2 norms (what CountSketch sees)
    l2_raw = np.linalg.norm(X_raw, axis=0)
    l2_log = np.linalg.norm(X_log, axis=0)
    l2_pearson = np.linalg.norm(X_pearson, axis=0)

    print("=" * 70)
    print("HYPOTHESIS 1: Sketching Sensitivity (L2 Energy Distribution)")
    print("=" * 70)
    print(f"Gini of per-gene L2 norms:")
    print(f"  Raw:     {gini(l2_raw):.4f}")
    print(f"  Pearson: {gini(l2_pearson):.4f}")
    print(f"  Log-CPM: {gini(l2_log):.4f}")
    print()

    # Top genes' share of total L2 energy
    def top_k_share(norms, k):
        sorted_norms = np.sort(norms)[::-1]
        return np.sum(sorted_norms[:k]**2) / np.sum(sorted_norms**2)

    print("Top genes' share of total L2² energy:")
    for k in [10, 50, 100]:
        print(f"  Top {k:3d}: Raw={top_k_share(l2_raw,k)*100:5.1f}%, "
              f"Pearson={top_k_share(l2_pearson,k)*100:5.1f}%, "
              f"Log={top_k_share(l2_log,k)*100:5.1f}%")
    print()

    # ==========================================
    # Hypothesis 2: Information Saturation
    # ==========================================
    print("=" * 70)
    print("HYPOTHESIS 2: Information Saturation (High-expression blindness)")
    print("=" * 70)

    # Find genes with different expression levels
    gene_means = X_raw.mean(axis=0)

    # Group genes by expression level
    low_expr = np.where((gene_means > 1) & (gene_means < 10))[0]
    mid_expr = np.where((gene_means >= 100) & (gene_means < 500))[0]
    high_expr = np.where(gene_means >= 1000)[0]

    print(f"Gene groups: Low (1-10): {len(low_expr)}, Mid (100-500): {len(mid_expr)}, High (>1000): {len(high_expr)}")

    if len(high_expr) > 0 and len(mid_expr) > 0:
        # Compare variance in transformed space
        var_log_mid = np.var(X_log[:, mid_expr])
        var_log_high = np.var(X_log[:, high_expr])
        var_pearson_mid = np.var(X_pearson[:, mid_expr])
        var_pearson_high = np.var(X_pearson[:, high_expr])

        print(f"\nVariance in transformed space:")
        print(f"  Log-CPM:  Mid={var_log_mid:.4f}, High={var_log_high:.4f}, Ratio={var_log_high/var_log_mid:.2f}")
        print(f"  Pearson:  Mid={var_pearson_mid:.4f}, High={var_pearson_high:.4f}, Ratio={var_pearson_high/var_pearson_mid:.2f}")
        print("→ If Pearson ratio << Log ratio, Pearson loses discriminative power in high-expression region")

    # Check actual value ranges
    print(f"\nValue ranges in transformed space:")
    print(f"  Log-CPM:  [{X_log.min():.2f}, {X_log.max():.2f}], std={X_log.std():.2f}")
    print(f"  Pearson:  [{X_pearson.min():.2f}, {X_pearson.max():.2f}], std={X_pearson.std():.2f}")
    print()

    # ==========================================
    # Hypothesis 3: Condition Number (on signature matrix)
    # ==========================================
    print("=" * 70)
    print("HYPOTHESIS 3: Condition Number (use cell-type centroids)")
    print("=" * 70)

    # Create pseudo-signature matrix by clustering
    from sklearn.cluster import KMeans

    n_types = 10
    print(f"Creating {n_types} cell-type centroids via K-means...")

    # Subsample for speed
    n_sample = min(2000, X_raw.shape[0])
    idx = np.random.choice(X_raw.shape[0], n_sample, replace=False)

    kmeans = KMeans(n_clusters=n_types, random_state=42, n_init=3)
    kmeans.fit(X_log[idx])

    # Get centroids in each space
    sig_log = kmeans.cluster_centers_

    # Transform centroids back to get raw centroids (approximate)
    kmeans_raw = KMeans(n_clusters=n_types, random_state=42, n_init=3)
    kmeans_raw.fit(X_raw[idx])
    sig_raw = kmeans_raw.cluster_centers_

    kmeans_pearson = KMeans(n_clusters=n_types, random_state=42, n_init=3)
    kmeans_pearson.fit(X_pearson[idx])
    sig_pearson = kmeans_pearson.cluster_centers_

    # Compute condition numbers
    def safe_cond(mat):
        s = np.linalg.svd(mat, compute_uv=False)
        return s[0] / (s[-1] + 1e-10)

    cond_raw = safe_cond(sig_raw)
    cond_pearson = safe_cond(sig_pearson)
    cond_log = safe_cond(sig_log)

    print(f"Condition number of signature matrix:")
    print(f"  Raw:     {cond_raw:.2e}")
    print(f"  Pearson: {cond_pearson:.2e}")
    print(f"  Log-CPM: {cond_log:.2e}")
    print()

    # ==========================================
    # Key Additional Test: Deconvolution Accuracy
    # ==========================================
    print("=" * 70)
    print("DECONVOLUTION TEST: Simulated spots with known proportions")
    print("=" * 70)

    # Create synthetic spots by mixing cell-type centroids
    np.random.seed(123)
    n_spots = 100
    true_props = np.random.dirichlet(np.ones(n_types), size=n_spots)

    # Create spots in log space and raw space
    Y_log = true_props @ sig_log
    Y_raw = true_props @ sig_raw
    Y_pearson = true_props @ sig_pearson

    # Add noise
    noise_level = 0.1
    Y_log += np.random.normal(0, noise_level * Y_log.std(), Y_log.shape)
    Y_raw += np.random.normal(0, noise_level * Y_raw.std(), Y_raw.shape)
    Y_pearson += np.random.normal(0, noise_level * Y_pearson.std(), Y_pearson.shape)

    # Solve NNLS
    from scipy.optimize import nnls

    def deconvolve_nnls(Y, X):
        """Simple NNLS deconvolution"""
        n_spots = Y.shape[0]
        n_types = X.shape[0]
        props = np.zeros((n_spots, n_types))
        for i in range(n_spots):
            props[i], _ = nnls(X.T, Y[i])
        # Normalize
        props = props / (props.sum(axis=1, keepdims=True) + 1e-10)
        return props

    pred_log = deconvolve_nnls(Y_log, sig_log)
    pred_raw = deconvolve_nnls(Y_raw, sig_raw)
    pred_pearson = deconvolve_nnls(Y_pearson, sig_pearson)

    # Compute correlations
    def mean_pearson(pred, true):
        corrs = []
        for i in range(pred.shape[0]):
            corrs.append(np.corrcoef(pred[i], true[i])[0, 1])
        return np.mean(corrs)

    corr_log = mean_pearson(pred_log, true_props)
    corr_raw = mean_pearson(pred_raw, true_props)
    corr_pearson = mean_pearson(pred_pearson, true_props)

    print(f"NNLS Deconvolution Accuracy (Pearson correlation):")
    print(f"  Raw:     {corr_raw:.4f}")
    print(f"  Pearson: {corr_pearson:.4f}")
    print(f"  Log-CPM: {corr_log:.4f}")
    print()

    # ==========================================
    # Visualization
    # ==========================================
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # A: L2 norm distribution
    ax1 = axes[0, 0]
    for norms, label, color in [(l2_raw, 'Raw', 'red'),
                                 (l2_pearson, 'Pearson', 'orange'),
                                 (l2_log, 'Log-CPM', 'green')]:
        sorted_norms = np.sort(norms)[::-1]
        cumsum = np.cumsum(sorted_norms**2) / np.sum(sorted_norms**2)
        ax1.plot(cumsum[:500], label=label, linewidth=2, color=color)
    ax1.axhline(y=0.5, color='gray', linestyle='--', alpha=0.5)
    ax1.set_title("A. Cumulative L2² Energy by Gene Rank", fontweight='bold')
    ax1.set_xlabel("Number of Top Genes")
    ax1.set_ylabel("Cumulative Energy Share")
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # B: Value distribution
    ax2 = axes[0, 1]
    # Only plot non-zero values
    log_nonzero = X_log[X_log > 0].flatten()
    pearson_nonzero = X_pearson[X_pearson > 0].flatten()

    # Subsample for plotting
    n_plot = min(100000, len(log_nonzero))
    ax2.hist(np.random.choice(log_nonzero, n_plot), bins=50, alpha=0.6,
             label='Log-CPM', density=True, color='green')
    ax2.hist(np.random.choice(pearson_nonzero, n_plot), bins=50, alpha=0.6,
             label='Pearson', density=True, color='orange')
    ax2.set_title("B. Value Distribution (non-zero)", fontweight='bold')
    ax2.set_xlabel("Transformed Value")
    ax2.set_ylabel("Density")
    ax2.legend()
    ax2.set_xlim(0, 15)

    # C: Deconvolution accuracy
    ax3 = axes[1, 0]
    methods = ['Raw', 'Pearson', 'Log-CPM']
    corrs = [corr_raw, corr_pearson, corr_log]
    colors = ['red', 'orange', 'green']
    bars = ax3.bar(methods, corrs, color=colors, alpha=0.7, edgecolor='black')
    ax3.set_title("C. Deconvolution Accuracy (NNLS)", fontweight='bold')
    ax3.set_ylabel("Mean Pearson Correlation")
    ax3.set_ylim(0, 1)
    for bar, c in zip(bars, corrs):
        ax3.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02,
                f'{c:.3f}', ha='center', fontsize=11, fontweight='bold')

    # D: Condition number
    ax4 = axes[1, 1]
    conds = [cond_raw, cond_pearson, cond_log]
    bars = ax4.bar(methods, np.log10(conds), color=colors, alpha=0.7, edgecolor='black')
    ax4.set_title("D. Condition Number (log10)", fontweight='bold')
    ax4.set_ylabel("log10(Condition Number)")
    for bar, c in zip(bars, conds):
        ax4.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.1,
                f'{c:.1e}', ha='center', fontsize=10)

    plt.tight_layout()
    plt.savefig('/Users/apple/Research/FlashDeconv/validation/real_data_hypotheses.png', dpi=150, bbox_inches='tight')
    plt.show()

    print("=" * 70)
    print("CONCLUSION")
    print("=" * 70)

if __name__ == "__main__":
    analyze_real_data()
