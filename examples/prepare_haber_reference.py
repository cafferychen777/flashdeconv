#!/usr/bin/env python3
"""
Download Haber et al. (2017) scRNA-seq reference for Visium HD mouse intestine deconvolution.

This script downloads the pre-processed mouse intestinal epithelium single-cell RNA-seq data
from Zenodo (processed by Besca with cell type annotations).

Reference:
    Haber AL, Biton M, Rogel N, Herbst RH et al.
    A single-cell survey of the small intestinal epithelium.
    Nature 2017 Nov 16;551(7680):333-339. PMID: 29144463

Data source:
    Zenodo: https://zenodo.org/records/4447233
    Contains 10,896 cells with detailed cell type annotations

Cell types available (celltype1):
    - epithelial fate stem cell (3,309 cells)
    - transit amplifying cell (2,776 cells)
    - enterocyte progenitor (952 cells)
    - immature goblet cell (933 cells)
    - enterocyte (876 cells)
    - goblet cell (663 cells)
    - paneth cell (647 cells)
    - enteroendocrine cell (325 cells)
    - brush cell (271 cells)
    - immature enterocyte (144 cells)

Usage:
    python prepare_haber_reference.py [--output OUTPUT_PATH]

Output:
    Creates `haber_intestine_reference.h5ad` containing:
    - Raw UMI counts (cells x genes)
    - Cell type annotations in .obs['celltype1']
    - Gene symbols as .var_names
"""

import argparse
import os
import urllib.request
from pathlib import Path

import scanpy as sc


# Zenodo URLs for Haber et al. data (processed by Besca)
ZENODO_RAW_URL = "https://zenodo.org/records/4447233/files/haber_raw.h5ad?download=1"
ZENODO_PROCESSED_URL = "https://zenodo.org/records/4447233/files/haber_processed.h5ad?download=1"


def download_with_progress(url: str, output_path: Path, description: str = "file") -> None:
    """Download a file with progress indicator."""
    if output_path.exists():
        print(f"  {description} already exists at {output_path}")
        return

    print(f"  Downloading {description}...")
    print(f"  URL: {url}")
    print(f"  This may take a few minutes...")

    def progress_hook(count, block_size, total_size):
        percent = int(count * block_size * 100 / total_size)
        print(f"\r  Progress: {percent}%", end="", flush=True)

    urllib.request.urlretrieve(url, output_path, progress_hook)
    print(f"\n  Downloaded to {output_path}")


def prepare_reference(
    processed_path: Path,
    raw_path: Path,
    output_path: Path,
    cell_type_key: str = "celltype1",
) -> None:
    """
    Prepare the reference by combining raw counts with cell type annotations.

    The processed file has cell type annotations but only HVGs.
    The raw file has all genes but no cell type annotations.
    We combine them to get raw counts with cell type annotations.
    """
    print("\n[2/3] Preparing reference data...")

    # Load processed data (has cell type annotations)
    print("  Loading processed data (for cell type annotations)...")
    adata_processed = sc.read_h5ad(processed_path)
    print(f"    Shape: {adata_processed.shape}")
    print(f"    Cell types: {adata_processed.obs[cell_type_key].nunique()}")

    # Load raw data (has all genes)
    print("  Loading raw data (for full gene set)...")
    adata_raw = sc.read_h5ad(raw_path)
    print(f"    Shape: {adata_raw.shape}")

    # Match cells using CELL column (obs_names differ between files)
    # processed has numeric index, raw has barcode as index
    # Both have CELL column with matching barcodes
    processed_cells = set(adata_processed.obs["CELL"])
    raw_cells = set(adata_raw.obs["CELL"])
    common_cells = processed_cells.intersection(raw_cells)
    print(f"    Common cells (by CELL column): {len(common_cells)}")

    # Filter raw data to common cells
    raw_mask = adata_raw.obs["CELL"].isin(common_cells)
    adata_raw = adata_raw[raw_mask].copy()

    # Create mapping from CELL to processed obs for annotation transfer
    processed_cell_to_idx = {cell: idx for idx, cell in enumerate(adata_processed.obs["CELL"])}

    # Add cell type annotations from processed data
    cell_types = []
    for cell in adata_raw.obs["CELL"]:
        proc_idx = processed_cell_to_idx[cell]
        cell_types.append(adata_processed.obs[cell_type_key].iloc[proc_idx])
    adata_raw.obs[cell_type_key] = cell_types

    # Also copy other useful columns
    for col in ["CONDITION", "sample_type", "donor", "region", "sample", "percent_mito"]:
        col_name = col
        if col not in adata_processed.obs.columns:
            if col + "_x" in adata_processed.obs.columns:
                col_name = col + "_x"
            else:
                continue
        values = []
        for cell in adata_raw.obs["CELL"]:
            proc_idx = processed_cell_to_idx[cell]
            values.append(adata_processed.obs[col_name].iloc[proc_idx])
        adata_raw.obs[col] = values

    return adata_raw


def main():
    parser = argparse.ArgumentParser(
        description="Download and prepare Haber et al. scRNA-seq reference for FlashDeconv",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--output",
        "-o",
        type=str,
        default="./haber_intestine_reference.h5ad",
        help="Output path for the h5ad file (default: ./haber_intestine_reference.h5ad)",
    )
    parser.add_argument(
        "--data-dir",
        type=str,
        default="./zenodo_data",
        help="Directory to store downloaded files (default: ./zenodo_data)",
    )
    parser.add_argument(
        "--keep-downloads",
        action="store_true",
        help="Keep downloaded intermediate files",
    )

    args = parser.parse_args()

    # Create directories
    data_dir = Path(args.data_dir)
    data_dir.mkdir(parents=True, exist_ok=True)
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("Preparing Haber et al. (2017) scRNA-seq Reference")
    print("=" * 60)
    print("\nData source: https://zenodo.org/records/4447233")
    print("Original publication: Haber et al., Nature 2017")

    # Download files
    print("\n[1/3] Downloading data from Zenodo...")
    raw_path = data_dir / "haber_raw.h5ad"
    processed_path = data_dir / "haber_processed.h5ad"

    download_with_progress(ZENODO_RAW_URL, raw_path, "haber_raw.h5ad (181 MB)")
    download_with_progress(ZENODO_PROCESSED_URL, processed_path, "haber_processed.h5ad (279 MB)")

    # Prepare reference
    adata = prepare_reference(processed_path, raw_path, output_path)

    # Save
    print("\n[3/3] Saving reference...")
    adata.write_h5ad(output_path)

    # Print summary
    print("\n" + "=" * 60)
    print("Reference Summary")
    print("=" * 60)
    print(f"  Output file: {output_path}")
    print(f"  Total cells: {adata.n_obs:,}")
    print(f"  Total genes: {adata.n_vars:,}")
    print(f"\n  Cell type distribution (celltype1):")
    for ct, count in adata.obs["celltype1"].value_counts().items():
        pct = count / adata.n_obs * 100
        print(f"    {ct}: {count:,} ({pct:.1f}%)")

    # Cleanup
    if not args.keep_downloads:
        print(f"\n  Cleaning up temporary files in {data_dir}...")
        raw_path.unlink(missing_ok=True)
        processed_path.unlink(missing_ok=True)
        try:
            data_dir.rmdir()
        except OSError:
            pass

    print("\n" + "=" * 60)
    print("Done! You can now use this reference with FlashDeconv:")
    print("=" * 60)
    print(
        f"""
    import scanpy as sc
    from flashdeconv import FlashDeconv

    # Load your Visium HD data
    adata_st = sc.read_10x_h5("path/to/filtered_feature_bc_matrix.h5")

    # Load reference
    adata_ref = sc.read_h5ad("{output_path}")

    # Run deconvolution
    model = FlashDeconv()
    proportions = model.deconvolve(adata_st, adata_ref, cell_type_key='celltype1')
    """
    )


if __name__ == "__main__":
    main()
