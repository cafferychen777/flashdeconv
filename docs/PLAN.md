# FlashDeconv Implementation Plan

## Project Structure

```
FlashDeconv/
├── flashdeconv/                    # Main package
│   ├── __init__.py                # Package init with version and API exports
│   ├── core/                       # Core algorithms
│   │   ├── __init__.py
│   │   ├── deconv.py              # Main FlashDeconv class (entry point)
│   │   ├── preprocessing.py       # Variance stabilization & Pearson residuals
│   │   ├── platform_correction.py # Platform effect estimation (γ_g)
│   │   ├── sketching.py           # Structure-preserving randomized sketching
│   │   ├── spatial.py             # Spatial graph & Laplacian regularization
│   │   └── solver.py              # Block Coordinate Descent solver (Numba JIT)
│   ├── utils/                      # Utility functions
│   │   ├── __init__.py
│   │   ├── sparse.py              # Sparse matrix operations
│   │   ├── graph.py               # Spatial graph construction
│   │   ├── genes.py               # HVG & marker gene selection
│   │   └── metrics.py             # Evaluation metrics (RMSE, correlation)
│   └── io/                         # Input/Output
│       ├── __init__.py
│       └── loader.py              # AnnData/scanpy compatible loader
├── tests/                          # Unit tests
│   ├── __init__.py
│   ├── test_preprocessing.py
│   ├── test_sketching.py
│   ├── test_spatial.py
│   ├── test_solver.py
│   └── test_integration.py
├── examples/                       # Usage examples
│   └── quickstart.py
├── docs/                           # Documentation
│   ├── FlashDeconv_methods.md     # (existing)
│   └── FlashDeconv_methods.tex    # (existing)
├── setup.py
├── pyproject.toml
├── requirements.txt
└── README.md
```

## Implementation Steps

### Phase 1: Project Setup
1. Create directory structure
2. Create `pyproject.toml` with dependencies
3. Create `requirements.txt`
4. Create package `__init__.py` files

### Phase 2: Core Preprocessing Module
5. Implement `preprocessing.py`:
   - `compute_size_factors(Y)`: Calculate N_i (sequencing depth per spot)
   - `estimate_overdispersion(Y)`: Estimate global θ parameter
   - `pearson_residuals(Y, mu, theta)`: Transform raw counts to residuals

### Phase 3: Platform Effect Correction
6. Implement `platform_correction.py`:
   - `pseudo_bulk_aggregate(Y)`: Aggregate ST data to pseudo-bulk
   - `estimate_platform_effect(Y_bulk, X_ref)`: Estimate γ_g via weighted regression
   - `apply_platform_correction(X, gamma)`: Scale reference by platform factors

### Phase 4: Structure-Preserving Sketching
7. Implement `genes.py`:
   - `select_hvg(Y, n_top)`: Select highly variable genes
   - `select_markers(X, n_markers)`: Select cell-type specific markers
   - `compute_leverage_scores(X)`: Calculate gene importance weights

8. Implement `sketching.py`:
   - `build_countsketch_matrix(n_genes, d, leverage_scores)`: Sparse sketch matrix Ω
   - `project_to_sketch(Y_tilde, X_tilde, Omega)`: Compute Y_sketch, X_sketch

### Phase 5: Spatial Graph Module
9. Implement `graph.py`:
   - `build_knn_graph(coords, k)`: K-nearest neighbor graph
   - `build_radius_graph(coords, radius)`: Radius-based neighbor graph
   - `coords_to_adjacency(coords, method)`: Convert coordinates to adjacency

10. Implement `spatial.py`:
    - `compute_laplacian(A)`: L = D - A
    - `neighbor_indices(A)`: Sparse neighbor lookup structure

### Phase 6: BCD Solver (Numba-Optimized)
11. Implement `solver.py`:
    - `precompute_gram_matrix(X_sketch)`: X_sketch^T @ X_sketch (K×K)
    - `bcd_update_spot(...)`: Single spot update (Numba JIT)
    - `bcd_solve(Y_sketch, X_sketch, L, lambda_, rho, max_iter, tol)`: Full solver
    - Parallelization with `numba.prange`

### Phase 7: Main API Class
12. Implement `deconv.py`:
    - `FlashDeconv` class with:
      - `__init__(sketch_dim, lambda_spatial, rho_sparsity, ...)`
      - `fit(Y, X, coords)`: Main deconvolution method
      - `fit_transform(...)`: Fit and return β
      - `get_cell_type_proportions()`: Return normalized β

### Phase 8: I/O and Integration
13. Implement `loader.py`:
    - `load_spatial_data(adata)`: Extract Y, coords from AnnData
    - `load_reference(adata_ref)`: Extract X from scRNA-seq reference
    - `result_to_anndata(beta, adata)`: Store results in AnnData.obsm

### Phase 9: Testing
14. Create unit tests for each module
15. Create integration test with synthetic data

### Phase 10: Documentation
16. Create README.md with installation and usage
17. Move existing .tex/.md to docs/

## Key Dependencies
- numpy >= 1.21
- scipy >= 1.7
- numba >= 0.56
- scikit-learn >= 1.0
- anndata >= 0.8
- scanpy >= 1.9 (optional, for preprocessing)

## Algorithm Parameters (Defaults)
- `sketch_dim`: 512
- `lambda_spatial`: auto-tuned
- `rho_sparsity`: 0.01
- `max_iter`: 100
- `tol`: 1e-4
- `n_hvg`: 2000
- `n_markers_per_type`: 50
