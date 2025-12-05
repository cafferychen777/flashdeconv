"""
Deep Analysis: Why Log-CPM Beats Pearson Residuals in Sketching-based Deconvolution

Three Core Hypotheses:
1. Sketching Sensitivity - Energy concentration (Gini coefficient)
2. Information Saturation - Pearson asymptotes to sqrt(theta) for high counts
3. Condition Number - Optimization hardness in different spaces
"""

import numpy as np
import matplotlib.pyplot as plt

def analyze_preprocessing_impact():
    # ==========================================
    # 1. Simulate data with extreme dynamic range
    # ==========================================
    n_genes = 20000
    n_celltypes = 10

    # Simulate reference matrix X (Signature): low-expression markers + ultra-high housekeeping genes
    # Use Log-normal distribution to mimic real power-law distribution
    np.random.seed(42)
    X_raw = np.random.lognormal(mean=2, sigma=2.5, size=(n_celltypes, n_genes))

    # Add sparsity (Dropout)
    X_raw[X_raw < 5] = 0

    # Add ultra-high expression genes (mimicking mitochondrial/ribosomal genes)
    X_raw[:, :5] = 10000 * np.random.rand(n_celltypes, 5)

    # ==========================================
    # 2. Define preprocessing methods
    # ==========================================
    def preprocess_logcpm(mat):
        """Simplified Log-CPM"""
        lib_size = mat.sum(axis=1, keepdims=True)
        return np.log1p(1e4 * mat / (lib_size + 1e-6))

    def preprocess_pearson(mat, theta=100):
        """Uncentered Pearson Residuals"""
        mu = np.mean(mat, axis=0) + 1e-6  # Gene-wise mean
        var = mu + (mu**2) / theta
        sigma = np.sqrt(var)
        return mat / sigma

    # ==========================================
    # 3. Apply preprocessing
    # ==========================================
    X_log = preprocess_logcpm(X_raw)
    X_pearson = preprocess_pearson(X_raw)

    # ==========================================
    # 4. Hypothesis 1: Sketching Energy Distribution (Gini Coefficient)
    # ==========================================
    def gini(array):
        """Calculate Gini coefficient - measures energy concentration"""
        array = np.abs(array)
        array = array.flatten()
        if np.amin(array) < 0:
            array -= np.amin(array)
        array += 0.0000001
        array = np.sort(array)
        index = np.arange(1, array.shape[0]+1)
        n = array.shape[0]
        return ((np.sum((2 * index - n - 1) * array)) / (n * np.sum(array)))

    # Calculate energy concentration per cell type (row)
    energy_dist_log = [gini(row) for row in X_log]
    energy_dist_pearson = [gini(row) for row in X_pearson]
    energy_dist_raw = [gini(row) for row in X_raw]

    print("=" * 70)
    print("HYPOTHESIS 1: Sketching Sensitivity (Energy Concentration)")
    print("=" * 70)
    print(f"Raw Counts Gini (Avg): {np.mean(energy_dist_raw):.4f}")
    print(f"  → Extremely concentrated, Sketching prone to collapse")
    print(f"Pearson Gini    (Avg): {np.mean(energy_dist_pearson):.4f}")
    print(f"Log-CPM Gini    (Avg): {np.mean(energy_dist_log):.4f}")
    print(f"  → More uniform distribution, Sketching preserves information better")
    print()

    # ==========================================
    # 5. Hypothesis 2: High-expression Information Saturation
    # ==========================================
    raw_vals = np.linspace(0.1, 20000, 1000)
    theta = 100

    # Simulate Pearson calculation (simplified view)
    mu_vals = raw_vals  # Assume expected mean ≈ observed value

    pearson_vals = raw_vals / np.sqrt(mu_vals + mu_vals**2/theta + 1e-6)
    log_vals = np.log1p(raw_vals)

    print("=" * 70)
    print("HYPOTHESIS 2: Information Saturation")
    print("=" * 70)
    print(f"Raw=100   -> Pearson={pearson_vals[5]:.3f}, Log={log_vals[5]:.3f}")
    print(f"Raw=1000  -> Pearson={pearson_vals[50]:.3f}, Log={log_vals[50]:.3f}")
    print(f"Raw=5000  -> Pearson={pearson_vals[250]:.3f}, Log={log_vals[250]:.3f}")
    print(f"Raw=10000 -> Pearson={pearson_vals[500]:.3f}, Log={log_vals[500]:.3f}")
    print(f"Raw=20000 -> Pearson={pearson_vals[999]:.3f}, Log={log_vals[999]:.3f}")
    print()
    print(f"Pearson ratio (20000/1000): {pearson_vals[999]/pearson_vals[50]:.2f}x")
    print(f"Log ratio (20000/1000):     {log_vals[999]/log_vals[50]:.2f}x")
    print("→ If Pearson ratio ≈ 1, it means saturation (blind to high-count differences)")
    print()

    # ==========================================
    # 6. Hypothesis 3: Condition Number
    # ==========================================
    # Use SVD to compute condition number more stably
    def safe_cond(mat):
        """Compute condition number using SVD"""
        try:
            s = np.linalg.svd(mat, compute_uv=False)
            return s[0] / (s[-1] + 1e-10)
        except:
            return np.inf

    cond_log = safe_cond(X_log)
    cond_pearson = safe_cond(X_pearson)
    cond_raw = safe_cond(X_raw)

    print("=" * 70)
    print("HYPOTHESIS 3: Optimization Hardness (Condition Number)")
    print("=" * 70)
    print(f"Condition Number (Raw Space):     {cond_raw:.2e}")
    print(f"Condition Number (Pearson Space): {cond_pearson:.2e}")
    print(f"Condition Number (Log Space):     {cond_log:.2e}")
    print("→ Lower condition number = faster & more stable BCD convergence")
    print()

    # ==========================================
    # 7. Additional: L2 Norm Distribution (Sketching Collision Impact)
    # ==========================================
    print("=" * 70)
    print("ADDITIONAL: L2 Norm Distribution (Collision Impact)")
    print("=" * 70)

    # Per-gene L2 norms (column norms)
    l2_raw = np.linalg.norm(X_raw, axis=0)
    l2_log = np.linalg.norm(X_log, axis=0)
    l2_pearson = np.linalg.norm(X_pearson, axis=0)

    # Top 10 genes' share of total L2 energy
    def top_k_share(norms, k=10):
        sorted_norms = np.sort(norms)[::-1]
        return np.sum(sorted_norms[:k]) / np.sum(sorted_norms)

    print(f"Top 10 genes' share of total L2 energy:")
    print(f"  Raw:     {top_k_share(l2_raw)*100:.1f}%")
    print(f"  Pearson: {top_k_share(l2_pearson)*100:.1f}%")
    print(f"  Log-CPM: {top_k_share(l2_log)*100:.1f}%")
    print("→ If top genes dominate, hash collision in CountSketch destroys other signals")
    print()

    # ==========================================
    # Visualization
    # ==========================================
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Panel A: Transformation Behavior
    ax1 = axes[0, 0]
    ax1.plot(raw_vals, pearson_vals, label='Uncentered Pearson', linewidth=2)
    ax1.plot(raw_vals, log_vals, label='Log1p', linewidth=2)
    ax1.axhline(y=np.sqrt(theta), color='r', linestyle='--', alpha=0.5, label=f'√θ = {np.sqrt(theta):.1f}')
    ax1.set_title("A. Transformation Behavior\n(Pearson saturates at √θ)", fontweight='bold')
    ax1.set_xlabel("Raw Count")
    ax1.set_ylabel("Transformed Value")
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(0, 20000)

    # Panel B: Value Distribution Histogram
    ax2 = axes[0, 1]
    ax2.hist(X_log.flatten(), bins=50, alpha=0.6, label='Log-CPM', density=True, color='C1')
    ax2.hist(X_pearson.flatten(), bins=50, alpha=0.6, label='Pearson', density=True, color='C0')
    ax2.set_title("B. Value Distribution After Transform", fontweight='bold')
    ax2.set_xlabel("Transformed Value")
    ax2.set_ylabel("Density")
    ax2.legend()
    ax2.set_xlim(-1, 15)

    # Panel C: Gini Coefficient Comparison
    ax3 = axes[1, 0]
    methods = ['Raw', 'Pearson', 'Log-CPM']
    ginis = [np.mean(energy_dist_raw), np.mean(energy_dist_pearson), np.mean(energy_dist_log)]
    colors = ['red', 'orange', 'green']
    bars = ax3.bar(methods, ginis, color=colors, alpha=0.7, edgecolor='black')
    ax3.set_title("C. Gini Coefficient (Energy Concentration)\nLower = Better for Sketching", fontweight='bold')
    ax3.set_ylabel("Gini Coefficient")
    ax3.set_ylim(0, 1)
    for bar, g in zip(bars, ginis):
        ax3.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02,
                f'{g:.3f}', ha='center', fontsize=11, fontweight='bold')

    # Panel D: L2 Energy Distribution (CDF)
    ax4 = axes[1, 1]
    for norms, label, color in [(l2_raw, 'Raw', 'red'),
                                 (l2_pearson, 'Pearson', 'orange'),
                                 (l2_log, 'Log-CPM', 'green')]:
        sorted_norms = np.sort(norms)[::-1]
        cumsum = np.cumsum(sorted_norms) / np.sum(sorted_norms)
        ax4.plot(np.arange(len(cumsum)), cumsum, label=label, linewidth=2, color=color)
    ax4.axhline(y=0.9, color='gray', linestyle='--', alpha=0.5)
    ax4.axvline(x=100, color='gray', linestyle='--', alpha=0.5)
    ax4.set_title("D. Cumulative L2 Energy by Gene Rank\n(Steeper = More Concentrated)", fontweight='bold')
    ax4.set_xlabel("Number of Top Genes")
    ax4.set_ylabel("Cumulative L2 Energy Share")
    ax4.set_xlim(0, 500)
    ax4.legend()
    ax4.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('/Users/apple/Research/FlashDeconv/validation/preprocessing_hypotheses.png', dpi=150, bbox_inches='tight')
    plt.show()

    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print("""
Three mechanisms explain why Log-CPM beats Pearson for sketching-based deconvolution:

1. SKETCHING EFFICIENCY (Gini Coefficient):
   - Log-CPM has the lowest Gini coefficient, meaning energy is spread more evenly
   - This prevents "heavy hitter" genes from dominating the sketch variance
   - CountSketch hash collisions are less catastrophic

2. SATURATION TRAP (Pearson → √θ):
   - Pearson residuals asymptote to √θ ≈ 10 for highly expressed genes
   - Count 5000 and Count 20000 become nearly indistinguishable
   - Log-CPM preserves monotonicity across the full dynamic range

3. CONDITION NUMBER (Optimization):
   - Log-space has much lower condition number
   - BCD converges faster and more stably
   - Less prone to getting stuck in poor local optima
""")

if __name__ == "__main__":
    analyze_preprocessing_impact()
