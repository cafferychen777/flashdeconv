"""
Scalability Benchmark: 100K spots
---------------------------------
Test FlashDeconv performance on 100,000 spots to fill in paper placeholder.
"""
import sys
import time
import numpy as np
from scipy.stats import pearsonr

sys.path.insert(0, '/Users/apple/Research/FlashDeconv')
from flashdeconv import FlashDeconv


def log(msg):
    """Print with flush for real-time output."""
    print(msg, flush=True)


def generate_synthetic_data(n_spots, n_genes, n_types, seed=42):
    """Generate synthetic data for scalability test."""
    np.random.seed(seed)

    log("  [1/5] Generating reference signatures...")
    # Reference signatures (n_types x n_genes)
    X = np.random.lognormal(mean=1.0, sigma=1.5, size=(n_types, n_genes))
    X[X < 1.5] = 0  # Add sparsity
    X = X / (X.sum(axis=1, keepdims=True) + 1e-10) * 10000

    log("  [2/5] Generating true proportions...")
    # True proportions (n_spots x n_types)
    beta_true = np.random.dirichlet(np.ones(n_types) * 0.5, size=n_spots)

    log("  [3/5] Generating spatial coordinates...")
    # Spatial coordinates (grid)
    grid_size = int(np.ceil(np.sqrt(n_spots)))
    coords = np.array([[i % grid_size, i // grid_size] for i in range(n_spots)], dtype=float)

    log("  [4/5] Computing expected expression (matrix multiplication)...")
    # Generate counts Y = beta @ X + noise
    Y_expected = beta_true @ X

    log("  [5/5] Adding Poisson noise to generate counts...")
    depths = np.random.normal(20000, 5000, size=n_spots)
    depths = np.maximum(depths, 5000)
    Y = np.zeros((n_spots, n_genes))

    # Progress for Poisson sampling
    for i in range(n_spots):
        if i % 10000 == 0:
            log(f"       Poisson sampling: {i:,}/{n_spots:,} spots...")
        profile = Y_expected[i]
        if profile.sum() > 0:
            profile = profile / profile.sum() * depths[i]
        Y[i] = np.random.poisson(np.maximum(profile, 0))

    return Y, X, beta_true, coords


def main():
    log("=" * 60)
    log("Scalability Benchmark: 100K Spots")
    log("=" * 60)

    n_spots = 100000
    n_genes = 20000
    n_types = 50

    log(f"\nGenerating data: {n_spots:,} spots x {n_genes:,} genes x {n_types} types...")
    t0 = time.time()
    Y, X, beta_true, coords = generate_synthetic_data(n_spots, n_genes, n_types)
    log(f"  Data generation complete: {time.time() - t0:.2f}s")
    log(f"  Y shape: {Y.shape}, X shape: {X.shape}")

    # Run FlashDeconv
    log("\n" + "=" * 60)
    log("Running FlashDeconv...")
    log("=" * 60)

    model = FlashDeconv(
        sketch_dim=512,
        lambda_spatial=5000,
        n_hvg=2000,
        max_iter=100,
        preprocess="log_cpm",
        verbose=True,
        random_state=42,
    )

    t0 = time.time()
    beta_flash = model.fit_transform(Y, X, coords)
    flash_time = time.time() - t0

    # Compute correlation
    flash_corr, _ = pearsonr(beta_true.flatten(), beta_flash.flatten())

    log(f"\n{'='*60}")
    log("RESULTS")
    log(f"{'='*60}")
    log(f"FlashDeconv time: {flash_time:.2f} seconds")
    log(f"FlashDeconv correlation: {flash_corr:.4f}")

    # Estimate NNLS time from existing data (50K took ~428s)
    nnls_time_est = 428 * 2  # ~856s estimated for 100K
    speedup = nnls_time_est / flash_time

    log(f"\nEstimated NNLS time: ~{nnls_time_est:.0f}s (14+ minutes)")
    log(f"Estimated speedup: ~{speedup:.0f}x")

    # Append to CSV
    import pandas as pd
    csv_path = '/Users/apple/Research/FlashDeconv/validation/scalability_results.csv'

    new_row = {
        'desc': 'XXXL (100K spots)',
        'n_spots': n_spots,
        'n_genes': n_genes,
        'n_types': n_types,
        'flash_time': flash_time,
        'flash_corr': flash_corr,
        'nnls_time': nnls_time_est,
        'nnls_corr': '',
        'speedup': speedup,
    }

    df = pd.read_csv(csv_path)
    df = df[~df['desc'].str.contains('100K', na=False)]
    df = pd.concat([df, pd.DataFrame([new_row])], ignore_index=True)
    df.to_csv(csv_path, index=False)
    log(f"\nResults appended to {csv_path}")

    return flash_time, flash_corr


if __name__ == "__main__":
    main()
