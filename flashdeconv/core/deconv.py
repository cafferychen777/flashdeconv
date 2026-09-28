"""
Main FlashDeconv class - the primary API for spatial transcriptomics deconvolution.

FlashDeconv combines:
1. Gene selection (spatial highly variable genes + reference markers) and
   log-normalization of spatial and reference expression
2. Structure-preserving gene representation (deterministic expected
   leverage-weighted CountSketch weights by default; the randomized
   CountSketch projection is available as a legacy option)
3. Spatial graph Laplacian regularization
4. Numba-accelerated Block Coordinate Descent solver
"""

import numbers
import warnings
from typing import Any, Dict, Literal, Optional, Tuple, Union

import numpy as np
from scipy import sparse

from flashdeconv.core.preprocessing import PREPROCESS_METHODS, preprocess_expression

# Type aliases
ArrayLike = Union[np.ndarray, sparse.spmatrix]
PreprocessMethod = Literal["log_cpm", "pearson", "raw"]
GeneWeighting = Literal["expected", "countsketch"]

GENE_WEIGHTING_OPTIONS = ("expected", "countsketch")
SPATIAL_METHODS = ("knn", "radius", "grid")


def _check_finite_nonnegative(name: str, M: ArrayLike) -> None:
    """Raise if ``M`` holds NaN/inf or negative values (O(nnz), no copies)."""
    values = M.data if sparse.issparse(M) else M
    if values.size == 0:
        return
    lo, hi = values.min(), values.max()
    if not (np.isfinite(lo) and np.isfinite(hi)):
        raise ValueError(f"{name} contains NaN or infinite values.")
    if lo < 0:
        raise ValueError(
            f"{name} contains negative values; FlashDeconv expects "
            f"non-negative counts (Y) and expression profiles (X)."
        )


def _validate_fit_inputs(
    Y: Any, X: Any, coords: Any
) -> Tuple[ArrayLike, np.ndarray, np.ndarray]:
    """
    Coerce and validate the inputs of :meth:`FlashDeconv.fit`.

    Dense inputs are converted with ``np.asarray`` (no copy or dtype change
    for ndarrays); sparse ``Y`` in a format without efficient row/column
    indexing (COO, LIL, DOK, BSR, DIA) is converted to CSR; a sparse ``X`` is
    densified (it is small: n_cell_types x n_genes).
    """
    if sparse.issparse(Y):
        if Y.format not in ("csr", "csc"):
            Y = Y.tocsr()
    else:
        Y = np.asarray(Y)
    X = X.toarray() if sparse.issparse(X) else np.asarray(X)
    coords = np.asarray(coords)

    if Y.ndim != 2:
        raise ValueError(f"Y must be 2D (n_spots, n_genes), got shape {Y.shape}")
    if X.ndim != 2:
        raise ValueError(f"X must be 2D (n_cell_types, n_genes), got shape {X.shape}")
    if coords.ndim != 2 or coords.shape[1] == 0:
        raise ValueError(
            f"coords must be 2D (n_spots, n_dims), got shape {coords.shape}"
        )
    if Y.shape[1] != X.shape[1]:
        raise ValueError(
            f"Gene dimension mismatch: Y has {Y.shape[1]} genes but "
            f"X has {X.shape[1]} genes. They must share the same gene "
            f"space (align before calling fit)."
        )
    if coords.shape[0] != Y.shape[0]:
        raise ValueError(
            f"Spot count mismatch: Y has {Y.shape[0]} spots but "
            f"coords has {coords.shape[0]} rows. Each spot needs "
            f"exactly one coordinate."
        )
    if Y.shape[0] == 0:
        raise ValueError("Y must contain at least one spot.")
    if X.shape[0] == 0:
        raise ValueError(
            "Reference matrix X must contain at least one cell type "
            "(X.shape[0] > 0). Check your reference filtering and "
            "cell_type_key mapping."
        )
    _check_finite_nonnegative("Y", Y)
    _check_finite_nonnegative("X", X)
    if not np.all(np.isfinite(coords)):
        raise ValueError("coords contains NaN or infinite values.")
    return Y, X, coords


class FlashDeconv:
    """
    Fast spatial transcriptomics deconvolution with spatial regularization.

    FlashDeconv efficiently estimates cell type proportions from spatial
    transcriptomics data using single-cell reference signatures.

    Parameters
    ----------
    sketch_dim : int, default=512
        Number of CountSketch buckets d. With ``gene_weighting="expected"``
        it enters only through the expected per-gene weights (the regression
        runs on the selected genes); with ``gene_weighting="countsketch"`` it
        is the dimension of the randomized sketch space.
    lambda_spatial : float or "auto", default="auto"
        Spatial regularization strength. Higher values encourage smoother
        spatial patterns. Default "auto" automatically tunes based on data scale.
    rho_sparsity : float, default=0.01
        L1 sparsity penalty as a dimensionless fraction. Internally scaled
        by mean(diag(G)) so that soft-thresholding is commensurate with the
        partial residual magnitude. Note: the non-negativity constraint already
        provides the primary sparsity; L1 offers marginal refinement, most
        effective when K >> number of true types per spot (e.g., fine-grained
        taxonomy with K > 30).
    n_hvg : int, default=2000
        Number of highly variable genes to select.
    n_markers_per_type : int, default=50
        Number of marker genes per cell type.
    spatial_method : str, default="knn"
        Method for spatial graph construction:
        - "knn": k-nearest neighbors
        - "radius": fixed radius (requires ``radius`` parameter)
        - "grid": auto-detect grid structure
    k_neighbors : int, default=6
        Number of neighbors for KNN graph (used when spatial_method="knn").
    radius : float, optional
        Radius for spatial graph construction (required when spatial_method="radius").
    max_iter : int, default=1000
        Maximum iterations for the BCD solver.
    tol : float, default=1e-4
        Convergence tolerance: the solver stops when the largest absolute
        change of any abundance, relative to the largest abundance, falls
        below ``tol``.
    preprocess : {"log_cpm", "pearson", "raw"}, default="log_cpm"
        Preprocessing of Y and X on the selected genes:

        - "log_cpm": ``log1p`` of expression normalized to 10,000 counts per
          spot / cell type (CP10k; the option name is historical).
          Recommended for count data.
        - "pearson": uncentered Pearson residuals Y/sigma, X/sigma with
          ``sigma^2 = mu + mu^2/100``.
        - "raw": no transformation (for pre-normalized, non-negative data).
    random_state : int, RandomState or None, default=0
        Random seed for the legacy CountSketch projection
        (``gene_weighting="countsketch"``). It has no effect with the
        deterministic default ``gene_weighting="expected"``.
    verbose : bool, default=False
        Whether to print progress.
    gene_weighting : {"expected", "countsketch"}, default="expected"
        Gene representation fed to the solver:

        - "expected": keep the selected genes (no projection) and scale gene
          ``j`` by ``w_j``, the exact expectation over random bucket
          assignment of its weight in the column-normalized leverage-weighted
          CountSketch with ``sketch_dim`` buckets (see
          :func:`flashdeconv.core.sketching.expected_countsketch_weights`).
          Deterministic; ``random_state`` is not used.
        - "countsketch": legacy randomized leverage-weighted CountSketch
          projection to ``sketch_dim`` dimensions (behaviour of versions
          <= 0.1.6, bit-identical for a given ``random_state``).

        With ``lambda_spatial="auto"`` the regularization adapts to the
        representation's scale; a fixed numeric ``lambda_spatial`` acts on
        the Gram matrix of the chosen representation, whose scale differs
        between the two options (roughly by a factor ``n_genes / sketch_dim``).

    Attributes
    ----------
    beta_ : ndarray of shape (n_spots, n_cell_types)
        Non-negative regression coefficients (abundances before row
        normalization; not absolute cell counts).
    proportions_ : ndarray of shape (n_spots, n_cell_types)
        Normalized cell type proportions (sum to 1; an all-zero row of
        ``beta_`` is assigned the uniform distribution).
    gene_idx_ : ndarray
        Indices (columns of ``Y``/``X``) of the genes used for deconvolution.
    gene_weights_ : ndarray or None
        Per-gene weights of the selected genes (``gene_weighting="expected"``),
        None for the legacy CountSketch projection.
    lambda_used_ : float
        Spatial regularization used by the solver.
    adjacency_ : scipy.sparse.csr_matrix of shape (n_spots, n_spots)
        Binary spatial adjacency matrix.
    Y_repr_, X_repr_ : array-like of shape (n_spots, p), (n_cell_types, p)
        Spatial data and reference in the solver's gene representation
        (p = number of selected genes for ``"expected"``, sparse when ``Y``
        is sparse; p = ``sketch_dim`` for ``"countsketch"``). Also available
        under the pre-0.2.0 names ``Y_sketch_`` / ``X_sketch_`` (deprecated
        aliases).
    info_ : dict
        Optimization information: ``converged``, ``n_iterations``,
        ``final_objective``, ``final_change``, ``objectives`` (verbose only),
        ``XtX`` (Gram matrix) and ``n_neighbors`` (per-spot degree).

    Example
    -------
    >>> from flashdeconv import FlashDeconv
    >>> model = FlashDeconv(sketch_dim=512)
    >>> proportions = model.fit_transform(Y, X, coords)
    """

    def __init__(
        self,
        sketch_dim: int = 512,
        lambda_spatial: Union[float, str] = "auto",
        rho_sparsity: float = 0.01,
        n_hvg: int = 2000,
        n_markers_per_type: int = 50,
        spatial_method: str = "knn",
        k_neighbors: int = 6,
        radius: Optional[float] = None,
        max_iter: int = 1000,
        tol: float = 1e-4,
        preprocess: PreprocessMethod = "log_cpm",
        random_state: Optional[int] = 0,
        verbose: bool = False,
        gene_weighting: GeneWeighting = "expected",
    ):
        # Parameter validation
        if isinstance(sketch_dim, float) and sketch_dim.is_integer():
            sketch_dim = int(sketch_dim)
        if not isinstance(sketch_dim, numbers.Integral) or isinstance(sketch_dim, bool):
            raise TypeError(f"sketch_dim must be an integer, got {sketch_dim!r}")
        if sketch_dim <= 0:
            raise ValueError(f"sketch_dim must be positive, got {sketch_dim}")
        if k_neighbors < 0:
            raise ValueError(f"k_neighbors must be non-negative, got {k_neighbors}")
        if max_iter < 0:
            raise ValueError(f"max_iter must be non-negative, got {max_iter}")
        if tol <= 0:
            raise ValueError(f"tol must be positive, got {tol}")
        if isinstance(lambda_spatial, str):
            if lambda_spatial != "auto":
                raise ValueError(
                    f"lambda_spatial must be 'auto' or a non-negative number, "
                    f"got {lambda_spatial!r}"
                )
        elif not lambda_spatial >= 0:
            raise ValueError(f"lambda_spatial must be non-negative, got {lambda_spatial}")
        if rho_sparsity < 0:
            raise ValueError(f"rho_sparsity must be non-negative, got {rho_sparsity}")
        if n_hvg < 0:
            raise ValueError(f"n_hvg must be non-negative, got {n_hvg}")
        if n_markers_per_type < 0:
            raise ValueError(f"n_markers_per_type must be non-negative, got {n_markers_per_type}")
        if spatial_method not in SPATIAL_METHODS:
            raise ValueError(
                f"spatial_method must be one of {SPATIAL_METHODS}, got {spatial_method!r}"
            )
        if spatial_method == "radius" and radius is None:
            raise ValueError("radius must be specified when spatial_method='radius'")
        if radius is not None and radius <= 0:
            raise ValueError(f"radius must be positive, got {radius}")
        if preprocess not in PREPROCESS_METHODS:
            raise ValueError(
                f"preprocess must be one of {PREPROCESS_METHODS}, got {preprocess!r}"
            )
        if gene_weighting not in GENE_WEIGHTING_OPTIONS:
            raise ValueError(
                f"gene_weighting must be 'expected' or 'countsketch', "
                f"got {gene_weighting!r}"
            )

        self.sketch_dim = sketch_dim
        self.lambda_spatial = lambda_spatial
        self.rho_sparsity = rho_sparsity
        self.n_hvg = n_hvg
        self.n_markers_per_type = n_markers_per_type
        self.spatial_method = spatial_method
        self.k_neighbors = k_neighbors
        self.radius = radius
        self.max_iter = max_iter
        self.tol = tol
        self.preprocess = preprocess
        self.random_state = random_state
        self.verbose = verbose
        self.gene_weighting = gene_weighting

        # Fitted attributes
        self.beta_ = None
        self.proportions_ = None
        self.gene_idx_ = None
        self.gene_weights_ = None
        self.Y_repr_ = None
        self.X_repr_ = None
        self.info_ = None
        self._sketch_matrix = None
        self._fitted = False

    # Deprecated aliases (pre-0.2.0 names; with the default
    # gene_weighting="expected" the representation is not a sketch).
    @property
    def Y_sketch_(self) -> Optional[ArrayLike]:
        """Deprecated alias of :attr:`Y_repr_`."""
        return self.Y_repr_

    @property
    def X_sketch_(self) -> Optional[np.ndarray]:
        """Deprecated alias of :attr:`X_repr_`."""
        return self.X_repr_

    def _preprocess_data(
        self,
        Y: ArrayLike,
        X: np.ndarray,
        method: PreprocessMethod,
    ) -> Tuple[ArrayLike, np.ndarray]:
        """Preprocess Y and X (see :func:`~flashdeconv.core.preprocessing.preprocess_expression`)."""
        return preprocess_expression(Y, X, method)

    def fit(
        self,
        Y: ArrayLike,
        X: np.ndarray,
        coords: np.ndarray,
        cell_type_names: Optional[np.ndarray] = None,
    ) -> "FlashDeconv":
        """
        Fit the deconvolution model.

        Parameters
        ----------
        Y : array-like of shape (n_spots, n_genes)
            Spatial transcriptomics count matrix (dense or sparse;
            non-negative and finite).
        X : ndarray of shape (n_cell_types, n_genes)
            Reference cell type signature matrix (same genes as ``Y``, same
            order; non-negative and finite).
        coords : ndarray of shape (n_spots, 2) or (n_spots, 3)
            Spatial coordinates of spots.
        cell_type_names : ndarray of shape (n_cell_types,), optional
            Cell type names.

        Returns
        -------
        self : FlashDeconv
            Fitted model.
        """
        from flashdeconv.core.sketching import (
            apply_gene_weights,
            expected_countsketch_weights,
            sketch_data,
        )
        from flashdeconv.core.spatial import auto_tune_lambda
        from flashdeconv.core.solver import bcd_solve, normalize_proportions
        from flashdeconv.utils.genes import select_informative_genes
        from flashdeconv.utils.graph import coords_to_adjacency

        Y, X, coords = _validate_fit_inputs(Y, X, coords)
        if cell_type_names is not None and len(cell_type_names) != X.shape[0]:
            raise ValueError(
                f"cell_type_names length ({len(cell_type_names)}) does not "
                f"match number of cell types in X ({X.shape[0]})."
            )

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
        if sparse.issparse(Y_subset) and Y_subset.format != "csr":
            Y_subset = Y_subset.tocsr()  # CSR for efficient row operations
        X_subset = X[:, gene_idx]

        # Step 2: Preprocessing
        if self.verbose:
            print(f"Step 2: Preprocessing with method='{self.preprocess}'...")

        Y_tilde, X_tilde = preprocess_expression(Y_subset, X_subset, self.preprocess)

        # Step 3: Structure-preserving gene representation
        if self.gene_weighting == "expected":
            if self.verbose:
                print("Step 3: Expected leverage-weighted CountSketch gene "
                      f"weights (d={self.sketch_dim})...")
            gene_weights = expected_countsketch_weights(
                leverage_scores, n_selected, sketch_dim=self.sketch_dim,
            )
            Y_repr, X_repr = apply_gene_weights(Y_tilde, X_tilde, gene_weights)
            self.gene_weights_ = gene_weights
            self._sketch_matrix = None
        else:
            if self.verbose:
                print(f"Step 3: Sketching {n_selected} genes to "
                      f"{self.sketch_dim} dimensions...")
            Y_repr, X_repr, sketch_matrix = sketch_data(
                Y_tilde, X_tilde,
                sketch_dim=self.sketch_dim,
                leverage_scores=leverage_scores,
                random_state=self.random_state,
            )
            self.gene_weights_ = None
            # Kept so that the bootstrap reuses the exact hash functions even
            # when random_state is None or a RandomState instance.
            self._sketch_matrix = sketch_matrix

        # Step 4: Build spatial graph
        if self.verbose:
            print("Step 4: Building spatial graph...")

        A = coords_to_adjacency(
            coords,
            method=self.spatial_method,
            k=self.k_neighbors,
            radius=self.radius,
        )
        self.adjacency_ = A

        if self.verbose:
            avg_neighbors = np.mean(np.asarray(A.sum(axis=1)).flatten())
            print(f"  Average neighbors per spot: {avg_neighbors:.1f}")

        # Step 5: Spatial regularization strength
        if self.lambda_spatial == "auto":
            lambda_ = auto_tune_lambda(Y_repr, X_repr, A)
        else:
            lambda_ = float(self.lambda_spatial)
        self.lambda_used_ = lambda_
        if self.verbose:
            print(f"Step 5: lambda = {lambda_:.4f}"
                  + (" (auto)" if self.lambda_spatial == "auto" else ""))

        # Step 6: Solve via BCD
        if self.verbose:
            print("Step 6: Solving via Block Coordinate Descent...")

        beta, info = bcd_solve(
            Y_repr, X_repr, A,
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

        # Representation-space data for uncertainty quantification
        self.Y_repr_ = Y_repr
        self.X_repr_ = X_repr
        # Original inputs and leverage scores for the bootstrap and the
        # reference diagnostics (references, not copies)
        self._Y_raw = Y
        self._X_raw = X
        self._leverage_scores = leverage_scores

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
            "gene_weighting": self.gene_weighting,
            "representation_dim": int(self.X_repr_.shape[1]),
            "lambda_spatial": self.lambda_used_,
            "rho_sparsity": self.rho_sparsity,
            "preprocess_method": self.preprocess,
            "converged": self.info_["converged"],
            "n_iterations": self.info_["n_iterations"],
            "final_objective": self.info_["final_objective"],
        }

    def compute_uncertainty(
        self,
        alpha: float = 0.05,
        method: str = "sandwich",
        fdr_q: float = 0.1,
    ) -> Dict[str, Any]:
        """
        Model-based analytical uncertainty (opt-in; not computed by ``fit``).

        The intervals quantify sampling uncertainty of the estimate under the
        fitted model. They do not include error from reference-tissue
        mismatch, cell types missing from the reference, or the log-space
        mixture approximation, which in benchmarks dominate the error
        against true cell-type fractions. They are therefore not calibrated
        confidence intervals for true proportions; use them to compare the
        stability of estimates across spots and cell types.

        Parameters
        ----------
        alpha : float, default=0.05
            Interval level (0.05 = 95%).
        method : {"sandwich", "model", "laplace_diag"}, default="sandwich"
            - "sandwich": per-spot active-set sandwich (HC0) variance using
              the full inverse Hessian ``(G_AA + lambda * n_neighbors I)^-1``,
              delta method to proportions. Types estimated at zero get a
              one-sided interval ``[0, upper]`` (Wald bound from a Newton
              step on the augmented active set, whose gradient includes the
              spatial pull ``lambda * sum_j beta_jk`` of the neighbours,
              truncated at 0). See
              :func:`flashdeconv.core.uncertainty.compute_active_set_uncertainty`.
            - "model": as "sandwich" but with the homoscedastic variance
              ``s^2 H_A^-1``.
            - "laplace_diag": deprecated pre-0.2.0 diagonal-Hessian
              approximation (``1/diag(H)``); anti-conservative even when the
              model holds, zero-width intervals for types estimated at zero.
        fdr_q : float, default=0.1
            Benjamini-Hochberg level for ``detected`` (sandwich/model only),
            applied over all spot x type entries to one-sided p-values for
            abundance > 0. In benchmarks the realised FDR against true cell
            presence was close to ``q`` only when the reference matched the
            tissue (about 0.03 at q=0.1 on colorectal Xenium bins with a
            matched reference) and 0.18-0.25 under reference mismatch.

        Returns
        -------
        uq : dict with keys:
            'entropy', 'residual_ss', 'residual_norm': ndarray (n_spots,)
            'se_prop': ndarray (n_spots, n_cell_types), standard error on the
                proportion scale (NaN for spots with no fitted signal)
            'var_prop': se_prop ** 2
            'ci_lower', 'ci_upper': ndarray (n_spots, n_cell_types)
            'ci_half_width': z * se_prop (before truncation to [0, 1])
            'cv': se_prop / proportion (0 where proportion is 0)
            'z_score', 'p_value': one-sided Wald test of abundance > 0
            'detected': bool, BH rejections at ``fdr_q``
            'detection_confident': alias of 'detected' (for "laplace_diag":
                the old rule lower bound > 0.01)
            'mean_ci_width': float
            'method', 'alpha', 'fdr_q'
        """
        if not self._fitted:
            raise RuntimeError("Model has not been fitted. Call fit() first.")

        from .uncertainty import (
            compute_entropy,
            compute_reconstruction_residual,
            compute_hessian_variance,
            compute_confidence_scores,
            compute_active_set_uncertainty,
            bh_detection,
        )

        if method not in ("sandwich", "model", "laplace_diag"):
            raise ValueError(
                "method must be 'sandwich', 'model' or 'laplace_diag', "
                f"got {method!r}"
            )
        if not 0.0 < alpha < 1.0:
            raise ValueError(f"alpha must be in (0, 1), got {alpha}")

        entropy = compute_entropy(self.proportions_)
        residual_ss, residual_norm = compute_reconstruction_residual(
            self.Y_repr_, self.X_repr_, self.beta_,
        )
        XtX = self.info_['XtX']
        n_neighbors = self.info_['n_neighbors']

        uq: Dict[str, Any] = {
            'entropy': entropy,
            'residual_ss': residual_ss,
            'residual_norm': residual_norm,
            'method': method,
            'alpha': alpha,
        }

        if method == "laplace_diag":
            warnings.warn(
                "method='laplace_diag' ignores collinearity between cell "
                "types and under-covers even when the model holds; use "
                "method='sandwich'.",
                DeprecationWarning,
                stacklevel=2,
            )
            _, var_prop = compute_hessian_variance(
                self.beta_, XtX, residual_ss, self.Y_repr_.shape[1],
                self.lambda_used_, n_neighbors,
            )
            conf = compute_confidence_scores(self.proportions_, var_prop, alpha)
            uq.update(conf)
            uq['var_prop'] = var_prop
            uq['se_prop'] = np.sqrt(np.maximum(var_prop, 0))
        else:
            from scipy.stats import norm
            rho_abs = self.rho_sparsity * float(np.mean(np.diag(XtX)))
            res = compute_active_set_uncertainty(
                self.Y_repr_, self.X_repr_, self.beta_,
                lambda_=self.lambda_used_, rho_abs=rho_abs,
                n_neighbors=n_neighbors, alpha=alpha, variance=method,
                neighbor_sum=np.asarray(self.adjacency_ @ self.beta_),
            )
            se = res['se_prop']
            z = norm.ppf(1 - alpha / 2)
            p_value = norm.sf(res['z_score'])
            detected = bh_detection(p_value, q=fdr_q)
            cv = np.zeros_like(self.proportions_)
            nz = (self.proportions_ > 1e-6) & np.isfinite(se)
            cv[nz] = se[nz] / self.proportions_[nz]
            uq.update({
                'se_prop': se,
                'var_prop': se ** 2,
                'ci_lower': res['ci_lower'],
                'ci_upper': res['ci_upper'],
                'ci_half_width': z * se,
                'cv': cv,
                'z_score': res['z_score'],
                'p_value': p_value,
                'detected': detected,
                'detection_confident': detected,
                'fdr_q': fdr_q,
            })

        uq['mean_ci_width'] = float(np.mean(uq['ci_upper'] - uq['ci_lower']))
        self.uncertainty_ = uq
        return uq

    def bootstrap_uncertainty(
        self,
        n_bootstrap: int = 100,
        max_iter_boot: Optional[int] = None,
        seed: int = 42,
        verbose: bool = False,
        tol: Optional[float] = None,
        alpha: float = 0.05,
    ) -> Dict[str, Any]:
        """
        Poisson bootstrap of the observed counts (model-based; opt-in).

        Each replicate resamples ``Y* ~ Poisson(Y)`` on the selected genes,
        re-applies the same preprocessing and the fitted gene representation
        (fixed gene weights, or the fitted CountSketch matrix), and
        re-solves with the fitted lambda and rho, warm-started at the fitted
        abundances. Genes, weights, penalties and the reference are held
        fixed, so the spread reflects count sampling under the fitted model
        only; like :meth:`compute_uncertainty`, it does not include
        reference-tissue mismatch or model misspecification.

        Parameters
        ----------
        n_bootstrap : int, default=100
        max_iter_boot : int, optional
            Iteration limit per refit; default is the model's ``max_iter``.
        seed : int, default=42
        verbose : bool
        tol : float, optional
            Convergence tolerance per refit; default is the model's ``tol``.
        alpha : float, default=0.05
            Percentile interval level.

        Returns
        -------
        boot : dict with keys 'boot_mean', 'boot_std', 'boot_ci_lower',
            'boot_ci_upper', 'boot_cv', 'n_bootstrap', 'n_converged',
            'n_iterations', 'final_change', 'warm_start'
        """
        if not self._fitted:
            raise RuntimeError("Model has not been fitted. Call fit() first.")

        from .uncertainty import poisson_bootstrap

        model_params = {
            'gene_idx': self.gene_idx_,
            'sketch_dim': self.sketch_dim,
            'preprocess': self.preprocess,
            'lambda_': self.lambda_used_,
            'rho': self.rho_sparsity,
            'leverage_scores': self._leverage_scores,
            'random_state': self.random_state,
            'adjacency': self.adjacency_,
            'gene_weighting': self.gene_weighting,
            'gene_weights': self.gene_weights_,
            'sketch_matrix': self._sketch_matrix,
        }

        boot = poisson_bootstrap(
            self._Y_raw, self._X_raw,
            coords=None,  # not needed; adjacency already stored
            model_params=model_params,
            beta_init=self.beta_,
            n_bootstrap=n_bootstrap,
            max_iter_boot=self.max_iter if max_iter_boot is None else max_iter_boot,
            tol=self.tol if tol is None else tol,
            seed=seed,
            verbose=verbose,
            alpha=alpha,
        )

        self.bootstrap_ = boot
        return boot

    def __repr__(self) -> str:
        status = "fitted" if self._fitted else "not fitted"
        return (
            f"FlashDeconv(sketch_dim={self.sketch_dim}, "
            f"gene_weighting={self.gene_weighting!r}, "
            f"lambda_spatial={self.lambda_spatial}, "
            f"status={status})"
        )
