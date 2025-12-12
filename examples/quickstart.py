"""
FlashDeconv Quick Start Example

This script demonstrates how to use FlashDeconv for spatial transcriptomics
deconvolution with synthetic data.
"""

import numpy as np
import time


def generate_synthetic_data(
    n_spots: int = 1000,
    n_genes: int = 2000,
    n_cell_types: int = 10,
    random_state: int = 42,
):
    """Generate synthetic spatial transcriptomics data."""
    print(f"Generating synthetic data...")
    print(f"  Spots: {n_spots}")
    print(f"  Genes: {n_genes}")
    print(f"  Cell types: {n_cell_types}")

    np.random.seed(random_state)

    # Reference signatures
    X = np.exp(np.random.randn(n_cell_types, n_genes) * 0.5 + 2)

    # Add cell-type specific markers
    for k in range(n_cell_types):
        markers = np.random.choice(n_genes, size=30, replace=False)
        X[k, markers] *= 5

    # Spatial coordinates (grid)
    side = int(np.ceil(np.sqrt(n_spots)))
    x = np.tile(np.arange(side), side)[:n_spots]
    y = np.repeat(np.arange(side), side)[:n_spots]
    coords = np.column_stack([x, y]).astype(float)
    coords += np.random.randn(n_spots, 2) * 0.1

    # True proportions with spatial patterns
    beta_true = np.zeros((n_spots, n_cell_types))
    for k in range(n_cell_types):
        center = np.random.rand(2) * side
        dist = np.sqrt(np.sum((coords - center) ** 2, axis=1))
        beta_true[:, k] = np.exp(-dist / (side / 3))

    beta_true = beta_true / beta_true.sum(axis=1, keepdims=True)

    # Generate counts
    expected = beta_true @ X
    depth = np.random.gamma(shape=5, scale=2000, size=n_spots)
    expected = expected * depth[:, np.newaxis]
    Y = np.random.poisson(expected).astype(float)

    # Cell type names
    cell_type_names = np.array([f"CellType_{i+1}" for i in range(n_cell_types)])

    return Y, X, coords, beta_true, cell_type_names


def main():
    """Run FlashDeconv on synthetic data."""
    print("=" * 60)
    print("FlashDeconv Quick Start Example")
    print("=" * 60)

    # Generate data
    Y, X, coords, beta_true, cell_type_names = generate_synthetic_data(
        n_spots=1000,
        n_genes=2000,
        n_cell_types=10,
    )

    print()

    # Import FlashDeconv
    from flashdeconv import FlashDeconv
    from flashdeconv.utils.metrics import evaluate_deconvolution

    # Create model
    model = FlashDeconv(
        sketch_dim=256,
        lambda_spatial="auto",
        rho_sparsity=0.01,
        n_hvg=1000,
        n_markers_per_type=30,
        k_neighbors=6,
        max_iter=100,
        verbose=True,
    )

    # Fit model
    print()
    print("Running FlashDeconv...")
    start_time = time.time()

    proportions = model.fit_transform(
        Y, X, coords,
        cell_type_names=cell_type_names,
    )

    elapsed = time.time() - start_time
    print(f"\nTotal time: {elapsed:.2f} seconds")
    print(f"Time per 1000 spots: {elapsed * 1000 / Y.shape[0]:.2f} seconds")

    # Evaluate results
    print()
    print("=" * 60)
    print("Evaluation Results")
    print("=" * 60)

    metrics = evaluate_deconvolution(proportions, beta_true, cell_type_names)

    print(f"\nOverall Metrics:")
    print(f"  RMSE: {metrics['overall']['rmse']:.4f}")
    print(f"  MAE: {metrics['overall']['mae']:.4f}")
    print(f"  Pearson Correlation: {metrics['overall']['pearson']:.4f}")
    print(f"  Spearman Correlation: {metrics['overall']['spearman']:.4f}")
    print(f"  Mean JSD: {metrics['overall']['mean_jsd']:.4f}")

    print(f"\nPer Cell Type Performance:")
    print(f"{'Cell Type':<15} {'RMSE':<10} {'Pearson':<10} {'True Mean':<10}")
    print("-" * 45)
    for name in cell_type_names:
        ct_metrics = metrics['per_cell_type'][name]
        print(
            f"{name:<15} "
            f"{ct_metrics['rmse']:<10.4f} "
            f"{ct_metrics['pearson']:<10.4f} "
            f"{ct_metrics['mean_proportion_true']:<10.4f}"
        )

    # Model summary
    print()
    print("Model Summary:")
    summary = model.summary()
    for key, value in summary.items():
        if key != 'fitted':
            print(f"  {key}: {value}")

    print()
    print("=" * 60)
    print("Example completed successfully!")
    print("=" * 60)


if __name__ == "__main__":
    main()
