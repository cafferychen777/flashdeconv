"""Input validation and edge cases of FlashDeconv.fit and its helpers."""

import warnings

import numpy as np
import pytest
from scipy import sparse

from flashdeconv import FlashDeconv, reference_fit_scores
from flashdeconv.utils.graph import build_knn_graph


@pytest.fixture(scope="module")
def small():
    rng = np.random.RandomState(0)
    G = 60
    X = rng.gamma(1.0, 1.0, (3, G))
    Y = rng.poisson(5, (40, G)).astype(float)
    coords = rng.rand(40, 2)
    return Y, X, coords


class TestInputValidation:
    @pytest.mark.parametrize("which", ["Y", "X"])
    @pytest.mark.parametrize("bad", [np.nan, np.inf])
    def test_non_finite_rejected(self, small, which, bad):
        Y, X, coords = (a.copy() for a in small)
        (Y if which == "Y" else X)[0, 0] = bad
        with pytest.raises(ValueError, match=f"{which} contains NaN or infinite"):
            FlashDeconv().fit(Y, X, coords)

    def test_non_finite_sparse_rejected(self, small):
        Y, X, coords = small
        Ys = sparse.csr_matrix(Y)
        Ys.data[0] = np.nan
        with pytest.raises(ValueError, match="NaN or infinite"):
            FlashDeconv().fit(Ys, X, coords)

    @pytest.mark.parametrize("which", ["Y", "X"])
    def test_negative_rejected(self, small, which):
        Y, X, coords = (a.copy() for a in small)
        (Y if which == "Y" else X)[0, 0] = -1.0
        with pytest.raises(ValueError, match="negative"):
            FlashDeconv().fit(Y, X, coords)

    def test_bad_coords(self, small):
        Y, X, coords = small
        with pytest.raises(ValueError, match="coords must be 2D"):
            FlashDeconv().fit(Y, X, coords[:, 0])
        c = coords.copy()
        c[0, 0] = np.nan
        with pytest.raises(ValueError, match="coords contains NaN"):
            FlashDeconv().fit(Y, X, c)

    def test_empty_spots_rejected(self, small):
        Y, X, coords = small
        with pytest.raises(ValueError, match="at least one spot"):
            FlashDeconv().fit(Y[:0], X, coords[:0])

    def test_no_genes_selected(self, small):
        Y, X, coords = small
        with pytest.raises(ValueError, match="No genes selected"):
            FlashDeconv(n_hvg=0, n_markers_per_type=0).fit(Y, X, coords)

    def test_cell_type_names_length(self, small):
        Y, X, coords = small
        with pytest.raises(ValueError, match="cell_type_names length"):
            FlashDeconv().fit(Y, X, coords, cell_type_names=["a", "b"])

    @pytest.mark.parametrize(
        "kwargs, match",
        [
            (dict(lambda_spatial="Auto"), "lambda_spatial"),
            (dict(lambda_spatial=-1.0), "lambda_spatial"),
            (dict(lambda_spatial=float("nan")), "lambda_spatial"),
            (dict(preprocess="log"), "preprocess"),
            (dict(spatial_method="delaunay"), "spatial_method"),
            (dict(gene_weighting="random"), "gene_weighting"),
            (dict(sketch_dim=0), "sketch_dim"),
            (dict(tol=0), "tol"),
            (dict(max_iter=-1), "max_iter"),
            (dict(k_neighbors=-1), "k_neighbors"),
            (dict(rho_sparsity=-0.1), "rho_sparsity"),
            (dict(radius=-1.0), "radius"),
        ],
    )
    def test_constructor_validation(self, kwargs, match):
        with pytest.raises(ValueError, match=match):
            FlashDeconv(**kwargs)

    def test_sketch_dim_type(self):
        with pytest.raises(TypeError):
            FlashDeconv(sketch_dim=8.5)
        assert FlashDeconv(sketch_dim=8.0).sketch_dim == 8

    def test_uq_argument_validation(self, small):
        m = FlashDeconv().fit(*small)
        with pytest.raises(ValueError, match="alpha"):
            m.compute_uncertainty(alpha=1.5)
        with pytest.raises(ValueError, match="method"):
            m.compute_uncertainty(method="bayes")
        with pytest.raises(ValueError, match="n_bootstrap"):
            m.bootstrap_uncertainty(n_bootstrap=0)


class TestInputCoercion:
    def test_lists_and_matrix_and_sparse_reference(self, small):
        Y, X, coords = small
        ref = FlashDeconv().fit(Y, X, coords).beta_
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", PendingDeprecationWarning)
            variants = [
                (Y.tolist(), X.tolist(), coords.tolist()),
                (np.asmatrix(Y), np.asmatrix(X), coords),
                (Y, sparse.csr_matrix(X), coords),
                (sparse.csr_array(Y), X, coords),
            ]
            for Yv, Xv, cv in variants:
                beta = FlashDeconv().fit(Yv, Xv, cv).beta_
                np.testing.assert_allclose(beta, ref, rtol=1e-12, atol=1e-12)


class TestDegenerateShapes:
    def test_single_spot(self, small):
        Y, X, coords = small
        m = FlashDeconv().fit(Y[:1], X, coords[:1])
        assert m.proportions_.shape == (1, 3)
        np.testing.assert_allclose(m.proportions_.sum(), 1.0)
        assert m.adjacency_.nnz == 0
        uq = m.compute_uncertainty()
        assert uq["ci_lower"].shape == (1, 3)
        reference_fit_scores(m)

    def test_single_cell_type(self, small):
        Y, X, coords = small
        m = FlashDeconv().fit(Y, X[:1], coords)
        np.testing.assert_array_equal(m.proportions_, np.ones((40, 1)))

    def test_few_genes(self, small):
        Y, X, coords = small
        m = FlashDeconv().fit(Y[:, :3], X[:, :3], coords)
        assert np.all(np.isfinite(m.proportions_))
        np.testing.assert_allclose(m.proportions_.sum(1), 1.0)

    def test_all_zero_spots_get_uniform(self, small):
        Y, X, coords = small
        m = FlashDeconv(lambda_spatial=0.0).fit(np.zeros_like(Y), X, coords)
        np.testing.assert_allclose(m.proportions_, 1.0 / 3)

    def test_all_zero_reference_type_gets_zero(self, small):
        Y, X, coords = small
        X = X.copy()
        X[1] = 0.0
        m = FlashDeconv().fit(Y, X, coords)
        np.testing.assert_array_equal(m.beta_[:, 1], 0.0)

    def test_no_spatial_graph(self, small):
        m = FlashDeconv(k_neighbors=0).fit(*small)
        assert m.adjacency_.nnz == 0
        assert np.all(np.isfinite(m.proportions_))

    def test_max_iter_zero_returns_start(self, small):
        m = FlashDeconv(max_iter=0).fit(*small)
        np.testing.assert_allclose(m.proportions_, 1.0 / 3)
        assert m.info_["n_iterations"] == 0


class TestKNNDuplicates:
    def test_unique_coords_have_k_out_neighbors(self):
        coords = np.random.RandomState(1).rand(50, 2)
        A = build_knn_graph(coords, k=4)
        assert A.diagonal().sum() == 0
        assert np.asarray(A.sum(1)).min() >= 4

    def test_many_duplicates_keep_k_neighbors(self):
        coords = np.zeros((10, 2))
        A = build_knn_graph(coords, k=3)
        assert A.diagonal().sum() == 0
        # Each spot selects exactly 3 others; symmetrization can only add.
        assert np.asarray(A.sum(1)).min() >= 3
        A_self = build_knn_graph(coords, k=3, include_self=True)
        np.testing.assert_array_equal(A_self.diagonal(), 1.0)
        np.testing.assert_array_equal((A_self - sparse.eye(10)).toarray(), A.toarray())


def test_deprecated_sketch_aliases(small):
    m = FlashDeconv()
    assert m.Y_sketch_ is None and m.X_sketch_ is None
    m.fit(*small)
    assert m.Y_sketch_ is m.Y_repr_
    assert m.X_sketch_ is m.X_repr_
    assert m.X_repr_.shape == (3, len(m.gene_idx_))
