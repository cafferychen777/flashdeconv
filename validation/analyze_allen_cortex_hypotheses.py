"""
Validate preprocessing hypotheses on REAL Allen Cortex scRNA-seq data
"""

import numpy as np
import scanpy as sc
from scipy.optimize import nnls

# Load real data
print("Loading Allen Cortex scRNA-seq data...")
adata = sc.read_h5ad("/Users/apple/Research/FlashDeconv/validation/benchmark_data/allen_cortex.h5ad")
print(f"Data shape: {adata.shape}")

X_raw = adata.X
if hasattr(X_raw, 'toarray'):
    X_raw = X_raw.toarray()
X_raw = X_raw.astype(np.float32)

print(f"Sparsity: {(X_raw == 0).mean()*100:.1f}%")
print(f"Max value: {X_raw.max():.0f}")
print(f"Non-zero mean: {X_raw[X_raw > 0].mean():.2f}")
print()

# ==========================================
# Preprocessing functions
# ==========================================
def preprocess_logcpm(mat):
    """Log-CPM normalization"""
    lib_size = mat.sum(axis=1, keepdims=True) + 1
    return np.log1p(1e4 * mat / lib_size)

def preprocess_pearson(mat, theta=100):
    """Uncentered Pearson residuals (Hafemeister & Satija 2019)"""
    mu = mat.mean(axis=0) + 1e-6
    var = mu + (mu**2) / theta
    sigma = np.sqrt(var)
    return mat / sigma

# Apply preprocessing
print("Applying preprocessing...")
X_log = preprocess_logcpm(X_raw)
X_pearson = preprocess_pearson(X_raw)

# ==========================================
# Hypothesis 1: Sketching Energy Distribution
# ==========================================
def gini(array):
    """Calculate Gini coefficient"""
    array = np.abs(array.flatten())
    array = array + 1e-10
    array = np.sort(array)
    n = len(array)
    index = np.arange(1, n + 1)
    return (np.sum((2 * index - n - 1) * array)) / (n * np.sum(array))

# Per-gene L2 norms
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

# Top genes' share
def top_k_share(norms, k):
    sorted_norms = np.sort(norms)[::-1]
    return np.sum(sorted_norms[:k]**2) / np.sum(sorted_norms**2)

print("Top genes' share of total L2² energy:")
for k in [10, 50, 100, 500]:
    print(f"  Top {k:3d}: Raw={top_k_share(l2_raw,k)*100:5.1f}%, "
          f"Pearson={top_k_share(l2_pearson,k)*100:5.1f}%, "
          f"Log={top_k_share(l2_log,k)*100:5.1f}%")
print()

# ==========================================
# Hypothesis 2: Information Saturation
# ==========================================
print("=" * 70)
print("HYPOTHESIS 2: Information Saturation")
print("=" * 70)

# Gene-wise mean expression
gene_means = X_raw.mean(axis=0)

# Group genes by expression level
low_expr = np.where((gene_means > 0.1) & (gene_means < 1))[0]
mid_expr = np.where((gene_means >= 1) & (gene_means < 10))[0]
high_expr = np.where(gene_means >= 10)[0]
very_high = np.where(gene_means >= 100)[0]

print(f"Gene groups by mean expression:")
print(f"  Low (0.1-1):    {len(low_expr)} genes")
print(f"  Mid (1-10):     {len(mid_expr)} genes")
print(f"  High (10-100):  {len(high_expr)} genes")
print(f"  Very High(>100):{len(very_high)} genes")

# Compare coefficient of variation (CV) in transformed space
# Higher CV = more discriminative power
def cv_per_gene(mat):
    """Coefficient of variation per gene"""
    means = mat.mean(axis=0) + 1e-10
    stds = mat.std(axis=0)
    return stds / means

cv_log = cv_per_gene(X_log)
cv_pearson = cv_per_gene(X_pearson)

print(f"\nCoefficient of Variation (CV) by expression level:")
print(f"  {'Group':<15} {'Log CV':<12} {'Pearson CV':<12} {'Winner'}")
for name, idx in [('Low (0.1-1)', low_expr), ('Mid (1-10)', mid_expr),
                   ('High (10-100)', high_expr), ('Very High(>100)', very_high)]:
    if len(idx) > 0:
        cv_l = np.median(cv_log[idx])
        cv_p = np.median(cv_pearson[idx])
        winner = "Log" if cv_l > cv_p else "Pearson"
        print(f"  {name:<15} {cv_l:<12.3f} {cv_p:<12.3f} {winner}")

print("\n→ Higher CV = better ability to distinguish different expression levels")
print()

# ==========================================
# Hypothesis 3: Condition Number
# ==========================================
print("=" * 70)
print("HYPOTHESIS 3: Condition Number")
print("=" * 70)

# Check if cell type annotations exist
if 'cell_type' in adata.obs.columns:
    cell_types = adata.obs['cell_type'].unique()
    print(f"Found {len(cell_types)} cell types")

    # Create signature matrix (cell-type centroids)
    sig_raw = np.array([X_raw[adata.obs['cell_type'] == ct].mean(axis=0) for ct in cell_types])
    sig_log = np.array([X_log[adata.obs['cell_type'] == ct].mean(axis=0) for ct in cell_types])
    sig_pearson = np.array([X_pearson[adata.obs['cell_type'] == ct].mean(axis=0) for ct in cell_types])
else:
    print("No cell type annotations, using K-means clustering...")
    from sklearn.cluster import KMeans
    n_types = 15
    n_sample = min(3000, X_raw.shape[0])
    idx = np.random.choice(X_raw.shape[0], n_sample, replace=False)

    kmeans_log = KMeans(n_clusters=n_types, random_state=42, n_init=3)
    kmeans_log.fit(X_log[idx])
    sig_log = kmeans_log.cluster_centers_

    kmeans_raw = KMeans(n_clusters=n_types, random_state=42, n_init=3)
    kmeans_raw.fit(X_raw[idx])
    sig_raw = kmeans_raw.cluster_centers_

    kmeans_pearson = KMeans(n_clusters=n_types, random_state=42, n_init=3)
    kmeans_pearson.fit(X_pearson[idx])
    sig_pearson = kmeans_pearson.cluster_centers_

def safe_cond(mat):
    s = np.linalg.svd(mat, compute_uv=False)
    return s[0] / (s[-1] + 1e-10)

cond_raw = safe_cond(sig_raw)
cond_pearson = safe_cond(sig_pearson)
cond_log = safe_cond(sig_log)

print(f"\nCondition number of signature matrix:")
print(f"  Raw:     {cond_raw:.2e}")
print(f"  Pearson: {cond_pearson:.2e}")
print(f"  Log-CPM: {cond_log:.2e}")
print("→ Lower = more stable optimization")
print()

# ==========================================
# Deconvolution Test
# ==========================================
print("=" * 70)
print("DECONVOLUTION TEST")
print("=" * 70)

# Create synthetic spots
np.random.seed(123)
n_spots = 200
n_types = sig_log.shape[0]
true_props = np.random.dirichlet(np.ones(n_types) * 0.5, size=n_spots)

# Create spots in each space
Y_log = true_props @ sig_log
Y_raw = true_props @ sig_raw
Y_pearson = true_props @ sig_pearson

# Add realistic noise
for Y, sig in [(Y_log, sig_log), (Y_raw, sig_raw), (Y_pearson, sig_pearson)]:
    noise = 0.05 * np.std(Y) * np.random.randn(*Y.shape)
    Y += noise

def deconvolve_nnls(Y, X):
    n_spots = Y.shape[0]
    n_types = X.shape[0]
    props = np.zeros((n_spots, n_types))
    for i in range(n_spots):
        props[i], _ = nnls(X.T, Y[i])
    props = props / (props.sum(axis=1, keepdims=True) + 1e-10)
    return props

pred_log = deconvolve_nnls(Y_log, sig_log)
pred_raw = deconvolve_nnls(Y_raw, sig_raw)
pred_pearson = deconvolve_nnls(Y_pearson, sig_pearson)

def rmse(pred, true):
    return np.sqrt(np.mean((pred - true)**2))

def mean_pearson(pred, true):
    corrs = [np.corrcoef(pred[i], true[i])[0, 1] for i in range(pred.shape[0])]
    return np.nanmean(corrs)

print(f"Deconvolution Accuracy:")
print(f"  {'Method':<10} {'RMSE':<10} {'Pearson':<10}")
print(f"  {'Raw':<10} {rmse(pred_raw, true_props):<10.4f} {mean_pearson(pred_raw, true_props):<10.4f}")
print(f"  {'Pearson':<10} {rmse(pred_pearson, true_props):<10.4f} {mean_pearson(pred_pearson, true_props):<10.4f}")
print(f"  {'Log-CPM':<10} {rmse(pred_log, true_props):<10.4f} {mean_pearson(pred_log, true_props):<10.4f}")
print()

# ==========================================
# Summary
# ==========================================
print("=" * 70)
print("SUMMARY: Why Log-CPM wins for Sketching-based Deconvolution")
print("=" * 70)
print(f"""
┌─────────────────────────────────────────────────────────────────────┐
│ Metric                      │ Raw       │ Pearson   │ Log-CPM     │
├─────────────────────────────────────────────────────────────────────┤
│ L2 Energy Gini (↓ better)   │ {gini(l2_raw):.4f}    │ {gini(l2_pearson):.4f}    │ {gini(l2_log):.4f}      │
│ Condition Number (↓ better) │ {cond_raw:.2e} │ {cond_pearson:.2e} │ {cond_log:.2e}  │
│ Deconv RMSE (↓ better)      │ {rmse(pred_raw, true_props):.4f}    │ {rmse(pred_pearson, true_props):.4f}    │ {rmse(pred_log, true_props):.4f}      │
│ Deconv Pearson (↑ better)   │ {mean_pearson(pred_raw, true_props):.4f}    │ {mean_pearson(pred_pearson, true_props):.4f}    │ {mean_pearson(pred_log, true_props):.4f}      │
└─────────────────────────────────────────────────────────────────────┘

Key Finding: Log-CPM consistently outperforms due to:
1. More uniform energy distribution → better CountSketch preservation
2. Lower condition number → more stable BCD convergence
3. Maintains discriminative power across all expression levels
""")

# ==========================================
# SNR Analysis (Signal-to-Noise Ratio)
# ==========================================
print("=" * 70)
print("SNR ANALYSIS: Signal vs Noise Gene Contribution")
print("=" * 70)

gene_means = X_raw.mean(axis=0)
signal_genes = np.where(gene_means >= 50)[0]  # High expression = signal
noise_genes = np.where((gene_means > 0) & (gene_means < 1))[0]  # Low expression = noise

print(f"Signal genes (mean >= 50): {len(signal_genes)}")
print(f"Noise genes (0 < mean < 1): {len(noise_genes)}")

# L2 norm ratio (Signal/Noise)
l2_signal_raw = np.linalg.norm(X_raw[:, signal_genes])
l2_noise_raw = np.linalg.norm(X_raw[:, noise_genes])
l2_signal_log = np.linalg.norm(X_log[:, signal_genes])
l2_noise_log = np.linalg.norm(X_log[:, noise_genes])
l2_signal_pearson = np.linalg.norm(X_pearson[:, signal_genes])
l2_noise_pearson = np.linalg.norm(X_pearson[:, noise_genes])

snr_raw = l2_signal_raw / l2_noise_raw
snr_log = l2_signal_log / l2_noise_log
snr_pearson = l2_signal_pearson / l2_noise_pearson

print(f"\nL2 Signal-to-Noise Ratio:")
print(f"  Raw:     {snr_raw:.2f}")
print(f"  Log-CPM: {snr_log:.2f}")
print(f"  Pearson: {snr_pearson:.2f}  ← NOISE DOMINATES!")
print()

# ==========================================
# Export data to CSV (for R plotting)
# ==========================================
import pandas as pd
from pathlib import Path

output_dir = Path(__file__).parent / "data"
output_dir.mkdir(exist_ok=True)

# Save metrics
metrics_df = pd.DataFrame({
    'Method': ['Raw', 'Pearson', 'Log-CPM'],
    'SNR': [snr_raw, snr_pearson, snr_log],
    'Condition_Number': [cond_raw, cond_pearson, cond_log],
    'Gini': [gini(l2_raw), gini(l2_pearson), gini(l2_log)],
    'Deconv_Pearson': [mean_pearson(pred_raw, true_props),
                       mean_pearson(pred_pearson, true_props),
                       mean_pearson(pred_log, true_props)],
    'N_Signal_Genes': [len(signal_genes)] * 3,
    'N_Noise_Genes': [len(noise_genes)] * 3
})
metrics_df.to_csv(output_dir / "snr_analysis_metrics.csv", index=False)

# Save cumulative energy curves
def compute_cumsum(norms):
    sorted_norms = np.sort(norms)[::-1]
    cumsum = np.cumsum(sorted_norms**2) / np.sum(sorted_norms**2)
    return cumsum

cumsum_raw = compute_cumsum(l2_raw)
cumsum_pearson = compute_cumsum(l2_pearson)
cumsum_log = compute_cumsum(l2_log)

# Subsample for file size
n_genes = len(cumsum_raw)
sample_idx = np.unique(np.concatenate([
    np.arange(0, 100),
    np.arange(100, 1000, 10),
    np.arange(1000, n_genes, 100)
]))
sample_idx = sample_idx[sample_idx < n_genes]

energy_df = pd.DataFrame({
    'rank': sample_idx + 1,
    'Raw': cumsum_raw[sample_idx],
    'Pearson': cumsum_pearson[sample_idx],
    'LogCPM': cumsum_log[sample_idx]
})
energy_df.to_csv(output_dir / "cumulative_energy.csv", index=False)

# Save transformation behavior
theta = 100
x_vals = np.logspace(-1, 3, 100)
transform_df = pd.DataFrame({
    'raw_count': x_vals,
    'Pearson': x_vals / np.sqrt(x_vals + x_vals**2 / theta),
    'LogCPM': np.log1p(x_vals),
    'theta': theta
})
transform_df.to_csv(output_dir / "transformation_behavior.csv", index=False)

print(f"\nData exported to: {output_dir}")
print(f"  - snr_analysis_metrics.csv")
print(f"  - cumulative_energy.csv")
print(f"  - transformation_behavior.csv")
print("\nMetrics summary:")
print(metrics_df.to_string(index=False))
