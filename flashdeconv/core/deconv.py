"""
Main FlashDeconv class - the primary API for spatial transcriptomics deconvolution.

FlashDeconv combines:
1. Variance-stabilizing transformation with platform effect correction
2. Structure-preserving randomized sketching
3. Spatial graph Laplacian regularization
4. Numba-accelerated Block Coordinate Descent solver
"""

import numpy as np
from scipy import sparse
from typing import Union, Optional, Tuple, Dict, Any, Literal

# Type alias
ArrayLike = Union[np.ndarray, sparse.spmatrix]
PreprocessMethod = Literal["log_cpm", "pearson", "raw"]


class FlashDeconv:
    """
    Fast spatial transcriptomics deconvolution with spatial regularization.

    FlashDeconv efficiently estimates cell type proportions from spatial
    transcriptomics data using single-cell reference signatures.

    Parameters
    ----------
    sketch_dim : int, default=512
        Dimension of the randomized sketch space. Higher values preserve
        more information but increase computation.
    lambda_spatial : float or "auto", default=5000.0
        Spatial regularization strength. Higher values encourage smoother
        spatial patterns. Tuned default (5000) works well for 10x Visium data.
        Set to "auto" for automatic tuning (may underestimate for real data).
        Typical range: 1000-10000 for Visium, lower for synthetic/clean data.
    rho_sparsity : float, default=0.01
        L1 sparsity regularization to encourage sparse solutions.
    n_hvg : int, default=2000
        Number of highly variable genes to select.
    n_markers_per_type : int, default=50
        Number of marker genes per cell type.
    spatial_method : str, default="knn"
        Method for spatial graph construction:
        - "knn": k-nearest neighbors
        - "radius": fixed radius
        - "grid": auto-detect grid structure
    k_neighbors : int, default=6
        Number of neighbors for KNN graph.
    max_iter : int, default=100
        Maximum iterations for BCD solver.
    tol : float, default=1e-4
        Convergence tolerance.
    preprocess : {"log_cpm", "pearson", "raw"}, default="log_cpm"
        Preprocessing method for Y and X matrices:
        - "log_cpm": Log1p of counts per million (recommended for most data)
        - "pearson": Uncentered Pearson residuals Y/σ, X/σ (variance stabilizing)
        - "raw": No preprocessing (for pre-normalized or synthetic data)
    random_state : int, optional
        Random seed for reproducibility.
    verbose : bool, default=False
        Whether to print progress.

    Attributes
    ----------
    beta_ : ndarray of shape (n_spots, n_cell_types)
        Estimated cell type abundances (after fitting).
    proportions_ : ndarray of shape (n_spots, n_cell_types)
        Normalized cell type proportions (sum to 1).
    gene_idx_ : ndarray
        Indices of genes used for deconvolution.
    info_ : dict
        Optimization information.

    Example
    -------
    >>> from flashdeconv import FlashDeconv
    >>> model = FlashDeconv(sketch_dim=512, lambda_spatial=5000)
    >>> proportions = model.fit_transform(Y, X, coords)
    """

    def __init__(
        self,
        sketch_dim: int = 512,
        lambda_spatial: Union[float, str] = 5000.0,
        rho_sparsity: float = 0.01,
        n_hvg: int = 2000,
        n_markers_per_type: int = 50,
        spatial_method: str = "knn",
        k_neighbors: int = 6,
        max_iter: int = 100,
        tol: float = 1e-4,
        preprocess: PreprocessMethod = "log_cpm",
        random_state: Optional[int] = None,
        verbose: bool = False,
    ):
        self.sketch_dim = sketch_dim
        self.lambda_spatial = lambda_spatial
        self.rho_sparsity = rho_sparsity
        self.n_hvg = n_hvg
        self.n_markers_per_type = n_markers_per_type
        self.spatial_method = spatial_method
        self.k_neighbors = k_neighbors
        self.max_iter = max_iter
        self.tol = tol
        self.preprocess = preprocess
        self.random_state = random_state
        self.verbose = verbose

        # Fitted attributes
        self.beta_ = None
        self.proportions_ = None
        self.gene_idx_ = None
        self.info_ = None
        self._fitted = False

    def _preprocess_data(
        self,
        Y: ArrayLike,
        X: np.ndarray,
        method: PreprocessMethod,
    ) -> Tuple[ArrayLike, np.ndarray]:
        """
        Preprocess Y and X matrices.

        Supports sparse Y matrices for memory efficiency.
        Key insight: log1p(0) = 0, so sparsity is preserved.

        Parameters
        ----------
        Y : array-like of shape (n_spots, n_genes)
            Spatial count matrix (sparse or dense).
        X : ndarray of shape (n_cell_types, n_genes)
            Reference signature matrix.
        method : {"log_cpm", "pearson", "raw"}
            Preprocessing method.

        Returns
        -------
        Y_norm : array-like
            Preprocessed Y matrix (sparse if input was sparse, for log_cpm).
        X_norm : ndarray
            Preprocessed X matrix.
        """
        from scipy.sparse import diags, issparse

        if method == "log_cpm":
            # CPM normalization + log1p (recommended)
            # Supports sparse matrices efficiently

            if issparse(Y):
                # Sparse-friendly implementation
                lib_size = np.array(Y.sum(axis=1)).flatten()
                lib_size[lib_size == 0] = 1.0
                D = diags(1e4 / lib_size)
                Y_norm = D @ Y
                # Log1p in-place on non-zero values (preserves sparsity, no copy!)
                Y_norm.data = np.log1p(Y_norm.data)
            else:
                Y_cpm = Y / (Y.sum(axis=1, keepdims=True) + 1e-10) * 1e4
                Y_norm = np.log1p(Y_cpm)

            # X is always small, dense is fine
            X_cpm = X / (X.sum(axis=1, keepdims=True) + 1e-10) * 1e4
            X_norm = np.log1p(X_cpm)

            return Y_norm, X_norm

        elif method == "pearson":
            # Uncentered Pearson residuals: divide by σ only (no centering)
            # This keeps values non-negative for NNLS to work
            # σ² = μ + μ²/θ (NB variance model)
            theta = 100.0  # dispersion parameter

            # For Y (spots): compute per-gene variance
            Y_mean = Y.mean(axis=0, keepdims=True) + 1e-6
            Y_var = Y_mean + Y_mean**2 / theta
            Y_sigma = np.sqrt(Y_var)
            Y_norm = Y / Y_sigma

            # For X (reference): use same formula
            X_mean = X.mean(axis=0, keepdims=True) + 1e-6
            X_var = X_mean + X_mean**2 / theta
            X_sigma = np.sqrt(X_var)
            X_norm = X / X_sigma

            return Y_norm, X_norm

        elif method == "raw":
            # No preprocessing
            return Y.astype(float), X.astype(float)

        else:
            raise ValueError(
                f"Unknown preprocess method: {method}. "
                f"Choose from 'log_cpm', 'pearson', or 'raw'."
            )

    def fit(
        self,
        Y: ArrayLike,
        X: np.ndarray,
        coords: np.ndarray,
        gene_names: Optional[np.ndarray] = None,
        cell_type_names: Optional[np.ndarray] = None,
    ) -> "FlashDeconv":
        """
        Fit the deconvolution model.

        Parameters
        ----------
        Y : array-like of shape (n_spots, n_genes)
            Spatial transcriptomics count matrix.
        X : ndarray of shape (n_cell_types, n_genes)
            Reference cell type signature matrix.
        coords : ndarray of shape (n_spots, 2) or (n_spots, 3)
            Spatial coordinates of spots.
        gene_names : ndarray of shape (n_genes,), optional
            Gene names for logging.
        cell_type_names : ndarray of shape (n_cell_types,), optional
            Cell type names.

        Returns
        -------
        self : FlashDeconv
            Fitted model.
        """
        from flashdeconv.core.sketching import sketch_data
        from flashdeconv.core.spatial import auto_tune_lambda
        from flashdeconv.core.solver import bcd_solve, normalize_proportions
        from flashdeconv.utils.genes import select_informative_genes
        from flashdeconv.utils.graph import coords_to_adjacency

        if self.verbose:
            print("FlashDeconv: Starting deconvolution...")
            print(f"  Spatial data: {Y.shape[0]} spots x {Y.shape[1]} genes")
            print(f"  Reference: {X.shape[0]} cell types x {X.shape[1]} genes")

        # Store metadata
        self.n_spots_ = Y.shape[0]
        self.n_genes_ = Y.shape[1]
        self.n_cell_types_ = X.shape[0]
        self.cell_type_names_ = cell_type_names

        # Step 1: Select informative genes
        if self.verbose:
            print("Step 1: Selecting informative genes...")

        gene_idx, leverage_scores = select_informative_genes(
            Y, X,
            n_hvg=self.n_hvg,
            n_markers_per_type=self.n_markers_per_type,
        )
        self.gene_idx_ = gene_idx
        n_selected = len(gene_idx)

        if self.verbose:
            print(f"  Selected {n_selected} genes (HVG + markers)")

        # Subset to selected genes (keep sparse if input was sparse)
        Y_subset = Y[:, gene_idx]
        if sparse.issparse(Y_subset) and not sparse.isspmatrix_csr(Y_subset):
            Y_subset = Y_subset.tocsr()  # Ensure CSR for efficient row operations
        X_subset = X[:, gene_idx]

        # Step 2: Preprocessing
        if self.verbose:
            print(f"Step 2: Preprocessing with method='{self.preprocess}'...")

        Y_tilde, X_tilde = self._preprocess_data(Y_subset, X_subset, self.preprocess)

        if self.verbose:
            if self.preprocess == "log_cpm":
                print("  Y and X normalized to log-CPM space")
            elif self.preprocess == "pearson":
                print("  Y and X transformed with uncentered Pearson residuals")
            else:
                print("  No preprocessing applied (raw)")

        # Step 4: Structure-preserving sketching
        if self.verbose:
            print(f"Step 4: Sketching to {self.sketch_dim} dimensions...")

        Y_sketch, X_sketch, Omega = sketch_data(
            Y_tilde, X_tilde,
            sketch_dim=self.sketch_dim,
            leverage_scores=leverage_scores,
            random_state=self.random_state,
        )

        if self.verbose:
            print(f"  Compressed {n_selected} genes -> {self.sketch_dim} dims")

        # Step 5: Build spatial graph
        if self.verbose:
            print("Step 5: Building spatial graph...")

        A = coords_to_adjacency(
            coords,
            method=self.spatial_method,
            k=self.k_neighbors,
        )
        self.adjacency_ = A

        avg_neighbors = np.mean(np.asarray(A.sum(axis=1)).flatten())
        if self.verbose:
            print(f"  Average neighbors per spot: {avg_neighbors:.1f}")

        # Step 6: Auto-tune lambda if needed
        if self.lambda_spatial == "auto":
            lambda_ = auto_tune_lambda(Y_sketch, X_sketch, A)
            if self.verbose:
                print(f"Step 6: Auto-tuned lambda = {lambda_:.4f}")
        else:
            lambda_ = float(self.lambda_spatial)
            if self.verbose:
                print(f"Step 6: Using lambda = {lambda_:.4f}")

        self.lambda_used_ = lambda_

        # Step 7: Solve via BCD
        if self.verbose:
            print("Step 7: Solving via Block Coordinate Descent...")

        beta, info = bcd_solve(
            Y_sketch, X_sketch, A,
            lambda_=lambda_,
            rho=self.rho_sparsity,
            max_iter=self.max_iter,
            tol=self.tol,
            verbose=self.verbose,
        )

        self.beta_ = beta
        self.proportions_ = normalize_proportions(beta)
        self.info_ = info
        self._fitted = True

        if self.verbose:
            print(f"  Converged: {info['converged']}")
            print(f"  Iterations: {info['n_iterations']}")
            print("FlashDeconv: Done!")

        return self

    def fit_transform(
        self,
        Y: ArrayLike,
        X: np.ndarray,
        coords: np.ndarray,
        **kwargs,
    ) -> np.ndarray:
        """
        Fit the model and return cell type proportions.

        Parameters
        ----------
        Y : array-like of shape (n_spots, n_genes)
            Spatial count matrix.
        X : ndarray of shape (n_cell_types, n_genes)
            Reference signatures.
        coords : ndarray of shape (n_spots, 2)
            Spatial coordinates.
        **kwargs
            Additional arguments passed to fit().

        Returns
        -------
        proportions : ndarray of shape (n_spots, n_cell_types)
            Cell type proportions (sum to 1 per spot).
        """
        self.fit(Y, X, coords, **kwargs)
        return self.proportions_

    def get_cell_type_proportions(self) -> np.ndarray:
        """
        Get normalized cell type proportions.

        Returns
        -------
        proportions : ndarray of shape (n_spots, n_cell_types)
            Cell type proportions.

        Raises
        ------
        RuntimeError
            If model has not been fitted.
        """
        if not self._fitted:
            raise RuntimeError("Model has not been fitted. Call fit() first.")
        return self.proportions_

    def get_abundances(self) -> np.ndarray:
        """
        Get raw (unnormalized) cell type abundances.

        Returns
        -------
        beta : ndarray of shape (n_spots, n_cell_types)
            Raw abundances.
        """
        if not self._fitted:
            raise RuntimeError("Model has not been fitted. Call fit() first.")
        return self.beta_

    def get_dominant_cell_type(self) -> np.ndarray:
        """
        Get the dominant cell type for each spot.

        Returns
        -------
        dominant : ndarray of shape (n_spots,)
            Index of dominant cell type per spot.
        """
        if not self._fitted:
            raise RuntimeError("Model has not been fitted. Call fit() first.")
        return np.argmax(self.proportions_, axis=1)

    def summary(self) -> Dict[str, Any]:
        """
        Get summary of fitted model.

        Returns
        -------
        summary : dict
            Model summary including parameters and fit statistics.
        """
        if not self._fitted:
            return {"fitted": False}

        return {
            "fitted": True,
            "n_spots": self.n_spots_,
            "n_cell_types": self.n_cell_types_,
            "n_genes_used": len(self.gene_idx_),
            "sketch_dim": self.sketch_dim,
            "lambda_spatial": self.lambda_used_,
            "rho_sparsity": self.rho_sparsity,
            "preprocess_method": self.preprocess,
            "converged": self.info_["converged"],
            "n_iterations": self.info_["n_iterations"],
            "final_objective": self.info_["final_objective"],
        }

    def __repr__(self) -> str:
        status = "fitted" if self._fitted else "not fitted"
        return (
            f"FlashDeconv(sketch_dim={self.sketch_dim}, "
            f"lambda_spatial={self.lambda_spatial}, "
            f"status={status})"
        )
