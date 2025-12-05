"""
Export SNR Analysis Data for Supplementary Figure
Allen Cortex scRNA-seq (14,249 cells × 34,617 genes)

Outputs:
- snr_analysis_metrics.csv: Summary metrics (SNR, Condition Number, etc.)
- snr_analysis_cumulative_energy.csv: Full cumulative L2² energy curves
"""

import numpy as np
import pandas as pd
import scanpy as sc
from pathlib import Path

# Output directory
output_dir = Path(__file__).parent / "data"
output_dir.mkdir(exist_ok=True)

print("Loading Allen Cortex data...")
adata = sc.read_h5ad("/Users/apple/Research/FlashDeconv/validation/benchmark_data/allen_cortex.h5ad")
print(f"Data shape: {adata.shape[0]} cells × {adata.shape[1]} genes")

# Get raw counts
X_raw = adata.X.toarray() if hasattr(adata.X, 'toarray') else adata.X

# ============================================================
# Preprocessing transformations
# ============================================================
print("Computing transformations...")

# Log-CPM (same as original: lib_size + 1 to avoid division by zero)
lib_size = X_raw.sum(axis=1, keepdims=True) + 1
X_log = np.log1p(1e4 * X_raw / lib_size)

# Pearson residuals (uncentered, fixed theta=100 as in original)
THETA = 100  # Fixed overdispersion parameter
gene_means = X_raw.mean(axis=0) + 1e-6
expected_var = gene_means + (gene_means**2) / THETA
X_pearson = X_raw / np.sqrt(expected_var)

# ============================================================
# 1. Cumulative L2² Energy (for Panel A)
# ============================================================
print("Computing cumulative L2² energy...")

def compute_cumulative_energy(X):
    """Compute cumulative L2² energy by gene rank"""
    l2_per_gene = np.linalg.norm(X, axis=0)
    l2_squared = l2_per_gene ** 2
    sorted_idx = np.argsort(l2_squared)[::-1]
    sorted_l2sq = l2_squared[sorted_idx]
    cumsum = np.cumsum(sorted_l2sq)
    cumsum_normalized = cumsum / cumsum[-1]
    return cumsum_normalized

cumsum_raw = compute_cumulative_energy(X_raw)
cumsum_pearson = compute_cumulative_energy(X_pearson)
cumsum_log = compute_cumulative_energy(X_log)

n_genes = len(cumsum_raw)

# Save full curves (subsampled for file size)
sample_indices = np.unique(np.concatenate([
    np.arange(0, 100),  # First 100 genes in detail
    np.arange(100, 1000, 10),  # Every 10th gene up to 1000
    np.arange(1000, n_genes, 100)  # Every 100th gene after that
]))
sample_indices = sample_indices[sample_indices < n_genes]

energy_df = pd.DataFrame({
    'rank': sample_indices + 1,  # 1-indexed
    'Raw': cumsum_raw[sample_indices],
    'Pearson': cumsum_pearson[sample_indices],
    'LogCPM': cumsum_log[sample_indices]
})
energy_df.to_csv(output_dir / "cumulative_energy.csv", index=False)
print(f"  Saved cumulative energy data ({len(sample_indices)} points)")

# ============================================================
# 2. SNR Analysis (for Panel C)
# ============================================================
print("Computing SNR...")

signal_genes = np.where(gene_means >= 50)[0]
noise_genes = np.where((gene_means > 0) & (gene_means < 1))[0]

def compute_snr(X, signal_idx, noise_idx):
    l2_signal = np.linalg.norm(X[:, signal_idx])
    l2_noise = np.linalg.norm(X[:, noise_idx])
    return l2_signal / l2_noise if l2_noise > 0 else np.inf

snr_raw = compute_snr(X_raw, signal_genes, noise_genes)
snr_pearson = compute_snr(X_pearson, signal_genes, noise_genes)
snr_log = compute_snr(X_log, signal_genes, noise_genes)

print(f"  Signal genes (mean >= 50): {len(signal_genes)}")
print(f"  Noise genes (0 < mean < 1): {len(noise_genes)}")
print(f"  SNR - Raw: {snr_raw:.2f}, Pearson: {snr_pearson:.2f}, Log-CPM: {snr_log:.2f}")

# ============================================================
# 3. Condition Number (for Panel E)
# ============================================================
print("Computing condition numbers...")

# Use K-means to get pseudo cell-type centroids (same as original script)
from sklearn.cluster import KMeans

n_clusters = 15  # Same as original
n_sample = min(3000, X_raw.shape[0])  # Sample cells for clustering
np.random.seed(42)
sample_idx = np.random.choice(X_raw.shape[0], n_sample, replace=False)

# Cluster in each space separately (same as original)
kmeans_log = KMeans(n_clusters=n_clusters, random_state=42, n_init=3)
kmeans_log.fit(X_log[sample_idx])
centroids_log = kmeans_log.cluster_centers_

kmeans_raw = KMeans(n_clusters=n_clusters, random_state=42, n_init=3)
kmeans_raw.fit(X_raw[sample_idx])
centroids_raw = kmeans_raw.cluster_centers_

kmeans_pearson = KMeans(n_clusters=n_clusters, random_state=42, n_init=3)
kmeans_pearson.fit(X_pearson[sample_idx])
centroids_pearson = kmeans_pearson.cluster_centers_

# Condition number of centroid matrix (using SVD like original)
def safe_cond(mat):
    s = np.linalg.svd(mat, compute_uv=False)
    return s[0] / (s[-1] + 1e-10)

cond_raw = safe_cond(centroids_raw)
cond_pearson = safe_cond(centroids_pearson)
cond_log = safe_cond(centroids_log)

print(f"  Condition Number - Raw: {cond_raw:.1f}, Pearson: {cond_pearson:.1f}, Log-CPM: {cond_log:.1f}")

# ============================================================
# 4. Gini Coefficient (for energy distribution analysis)
# ============================================================
print("Computing Gini coefficients...")

def gini(x):
    x = np.sort(x)
    n = len(x)
    cumsum = np.cumsum(x)
    return (2 * np.sum((np.arange(1, n+1) * x)) - (n + 1) * cumsum[-1]) / (n * cumsum[-1])

l2_raw = np.linalg.norm(X_raw, axis=0)
l2_pearson = np.linalg.norm(X_pearson, axis=0)
l2_log = np.linalg.norm(X_log, axis=0)

gini_raw = gini(l2_raw)
gini_pearson = gini(l2_pearson)
gini_log = gini(l2_log)

print(f"  Gini - Raw: {gini_raw:.4f}, Pearson: {gini_pearson:.4f}, Log-CPM: {gini_log:.4f}")

# ============================================================
# 5. Deconvolution Accuracy (synthetic test)
# ============================================================
print("Computing deconvolution accuracy...")

# Create synthetic mixture (same as original script)
np.random.seed(123)  # Same seed as original
n_test = 200  # Same as original
true_props = np.random.dirichlet(np.ones(n_clusters) * 0.5, n_test)  # Same alpha

# Generate synthetic bulk from centroids
bulk_raw = true_props @ centroids_raw
bulk_pearson = true_props @ centroids_pearson
bulk_log = true_props @ centroids_log

# Add realistic noise (5%, same as original)
bulk_raw += 0.05 * np.std(bulk_raw) * np.random.randn(*bulk_raw.shape)
bulk_pearson += 0.05 * np.std(bulk_pearson) * np.random.randn(*bulk_pearson.shape)
bulk_log += 0.05 * np.std(bulk_log) * np.random.randn(*bulk_log.shape)

# Simple NNLS deconvolution
from scipy.optimize import nnls

def deconvolve_nnls(bulk, centroids):
    pred_props = np.zeros((bulk.shape[0], centroids.shape[0]))
    for i in range(bulk.shape[0]):
        pred, _ = nnls(centroids.T, bulk[i])
        pred_props[i] = pred / (pred.sum() + 1e-10)
    return pred_props

pred_raw = deconvolve_nnls(bulk_raw, centroids_raw)
pred_pearson = deconvolve_nnls(bulk_pearson, centroids_pearson)
pred_log = deconvolve_nnls(bulk_log, centroids_log)

# Compute correlation per sample (same as original)
def mean_pearson_corr(true, pred):
    corrs = [np.corrcoef(pred[i], true[i])[0, 1] for i in range(pred.shape[0])]
    return np.nanmean(corrs)

deconv_raw = mean_pearson_corr(true_props, pred_raw)
deconv_pearson = mean_pearson_corr(true_props, pred_pearson)
deconv_log = mean_pearson_corr(true_props, pred_log)

print(f"  Deconv Pearson r - Raw: {deconv_raw:.4f}, Pearson: {deconv_pearson:.4f}, Log-CPM: {deconv_log:.4f}")

# ============================================================
# Save Summary Metrics
# ============================================================
print("\nSaving summary metrics...")

metrics_df = pd.DataFrame({
    'Method': ['Raw', 'Pearson', 'Log-CPM'],
    'SNR': [snr_raw, snr_pearson, snr_log],
    'Condition_Number': [cond_raw, cond_pearson, cond_log],
    'Gini': [gini_raw, gini_pearson, gini_log],
    'Deconv_Pearson': [deconv_raw, deconv_pearson, deconv_log],
    'N_Signal_Genes': [len(signal_genes)] * 3,
    'N_Noise_Genes': [len(noise_genes)] * 3
})

metrics_df.to_csv(output_dir / "snr_analysis_metrics.csv", index=False)
print(f"  Saved: {output_dir / 'snr_analysis_metrics.csv'}")

# ============================================================
# Save Transformation Behavior Data (for Panel B)
# ============================================================
print("Saving transformation behavior data...")

# Use fixed theta=100 (same as in preprocessing)
x_vals = np.logspace(-1, 3, 100)
pearson_transform = x_vals / np.sqrt(x_vals + x_vals**2 / THETA)
log_transform = np.log1p(x_vals)

transform_df = pd.DataFrame({
    'raw_count': x_vals,
    'Pearson': pearson_transform,
    'LogCPM': log_transform,
    'theta': THETA
})
transform_df.to_csv(output_dir / "transformation_behavior.csv", index=False)
print(f"  Saved: {output_dir / 'transformation_behavior.csv'}")

# ============================================================
# Print Final Summary
# ============================================================
print("\n" + "=" * 70)
print("DATA EXPORT COMPLETE")
print("=" * 70)
print(f"\nDataset: Allen Cortex scRNA-seq")
print(f"  Cells: {adata.shape[0]:,}")
print(f"  Genes: {adata.shape[1]:,}")
print(f"\nFiles saved to: {output_dir.absolute()}")
print(f"  - snr_analysis_metrics.csv")
print(f"  - cumulative_energy.csv")
print(f"  - transformation_behavior.csv")

print("\n" + "-" * 70)
print("SUMMARY TABLE (for paper)")
print("-" * 70)
print(metrics_df.to_string(index=False))
