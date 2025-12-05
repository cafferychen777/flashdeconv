"""
Create combined figures for FlashDeconv paper.
"""

import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.gridspec import GridSpec
import numpy as np
from pathlib import Path

# Paths
RESULTS_DIR = Path("../validation/results/leverage_deep_dive")
BENCHMARK_DIR = Path("../validation/results/benchmark_sketching_weights")
OUTPUT_DIR = Path("figures")

def create_figure2():
    """
    Create Figure 2: Leverage scores decouple biological identity from population abundance.

    Panel A: Abundance Invariance Test
    Panel B: Variance-Leverage Plane (Gene Quadrant)
    Panel C: GO Enrichment
    Panel D: Spatial Verification

    All panels generated with unified 10:8 aspect ratio for clean composition.
    """
    print("Creating Figure 2...")

    # Load UNIFIED aspect ratio images (all 10:8)
    img_a = mpimg.imread(RESULTS_DIR / "panel_a_abundance.png")
    img_b = mpimg.imread(RESULTS_DIR / "panel_b_quadrant.png")
    img_c = mpimg.imread(RESULTS_DIR / "panel_c_go.png")
    img_d = mpimg.imread(RESULTS_DIR / "panel_d_spatial.png")

    # Create figure - 2x2 grid with equal cells
    fig, axes = plt.subplots(2, 2, figsize=(16, 13))

    # Panel a
    axes[0, 0].imshow(img_a)
    axes[0, 0].axis('off')
    axes[0, 0].set_title('a', fontsize=22, fontweight='bold', loc='left', x=-0.02, y=1.0)

    # Panel b
    axes[0, 1].imshow(img_b)
    axes[0, 1].axis('off')
    axes[0, 1].set_title('b', fontsize=22, fontweight='bold', loc='left', x=-0.02, y=1.0)

    # Panel c
    axes[1, 0].imshow(img_c)
    axes[1, 0].axis('off')
    axes[1, 0].set_title('c', fontsize=22, fontweight='bold', loc='left', x=-0.02, y=1.0)

    # Panel d
    axes[1, 1].imshow(img_d)
    axes[1, 1].axis('off')
    axes[1, 1].set_title('d', fontsize=22, fontweight='bold', loc='left', x=-0.02, y=1.0)

    plt.tight_layout(pad=0.5)

    # Save at high resolution
    plt.savefig(OUTPUT_DIR / "figure2_mechanism.png", dpi=600, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.savefig(OUTPUT_DIR / "figure2_mechanism.pdf", dpi=600, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    print(f"  Saved: {OUTPUT_DIR / 'figure2_mechanism.pdf'}")
    plt.close()


def create_supplementary_sketching():
    """
    Create Supplementary Figure: Leverage-weighted sketching benchmark.

    Panel A: Spatial comparison (Endo)
    Panel B: Multi-cell type comparison bar chart
    """
    print("Creating Supplementary Figure (Sketching Benchmark)...")

    # Load images
    img_a = mpimg.imread(BENCHMARK_DIR / "sketching_comparison_endo.png")
    img_b = mpimg.imread(BENCHMARK_DIR / "multi_celltype_comparison.png")

    # Create figure
    fig = plt.figure(figsize=(16, 12))
    gs = GridSpec(2, 1, figure=fig, height_ratios=[1, 1], hspace=0.25)

    # Panel a - Spatial comparison
    ax_a = fig.add_subplot(gs[0])
    ax_a.imshow(img_a)
    ax_a.axis('off')
    ax_a.set_title('a', fontsize=20, fontweight='bold', loc='left', x=-0.01, y=1.02)

    # Panel b - Bar chart
    ax_b = fig.add_subplot(gs[1])
    ax_b.imshow(img_b)
    ax_b.axis('off')
    ax_b.set_title('b', fontsize=20, fontweight='bold', loc='left', x=-0.01, y=1.02)

    plt.tight_layout()

    # Save at high resolution
    plt.savefig(OUTPUT_DIR / "supplementary_sketching_benchmark.png", dpi=600,
                bbox_inches='tight', facecolor='white', edgecolor='none')
    plt.savefig(OUTPUT_DIR / "supplementary_sketching_benchmark.pdf", dpi=600,
                bbox_inches='tight', facecolor='white', edgecolor='none')
    print(f"  Saved: {OUTPUT_DIR / 'supplementary_sketching_benchmark.pdf'}")
    plt.close()


def main():
    print("="*60)
    print("Creating Paper Figures")
    print("="*60)

    # Ensure output directory exists
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Create figures
    create_figure2()
    create_supplementary_sketching()

    print("\nDone!")
    print(f"Figures saved to: {OUTPUT_DIR.absolute()}")


if __name__ == "__main__":
    main()
