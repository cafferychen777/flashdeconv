"""Tests for preprocessing, gene selection, metrics and AnnData I/O helpers."""

import numpy as np
import pytest
from scipy import sparse

from flashdeconv.core.preprocessing import PREPROCESS_METHODS, preprocess_expression
from flashdeconv.utils import genes, metrics


@pytest.fixture(scope="module")
def counts():
    rng = np.random.RandomState(3)
    Y = rng.poisson(rng.gamma(0.5, 2.0, size=(30, 50))).astype(float)
    Y[4] = 0.0
    X = rng.gamma(1.0, 1.0, size=(3, 50))
    return Y, X


class TestPreprocessing:
    def test_log_cpm_is_log1p_cp10k(self, counts):
        Y, X = counts
        Yt, Xt = preprocess_expression(Y, X, "log_cpm")
        lib = Y.sum(1, keepdims=True)
        expected = np.log1p(np.divide(Y, lib, out=np.zeros_like(Y), where=lib > 0) * 1e4)
        np.testing.assert_allclose(Yt, expected, rtol=1e-12)
        np.testing.assert_allclose(np.expm1(Xt).sum(1), 1e4, rtol=1e-10)
        np.testing.assert_array_equal(Yt[4], 0.0)

    @pytest.mark.parametrize("method", PREPROCESS_METHODS)
    def test_sparse_matches_dense(self, counts, method):
        Y, X = counts
        Yd, Xd = preprocess_expression(Y, X, method)
        Ys, Xs = preprocess_expression(sparse.csr_matrix(Y), X, method)
        assert sparse.issparse(Ys)
        np.testing.assert_allclose(Ys.toarray(), Yd, rtol=1e-12, atol=1e-12)
        np.testing.assert_array_equal(Xs, Xd)

    def test_pearson_nonnegative(self, counts):
        Yt, Xt = preprocess_expression(*counts, "pearson")
        assert Yt.min() >= 0 and Xt.min() >= 0

    def test_unknown_method(self, counts):
        with pytest.raises(ValueError, match="Unknown preprocess method"):
            preprocess_expression(*counts, "tpm")

    def test_class_wrapper_delegates(self, counts):
        from flashdeconv import FlashDeconv

        a = FlashDeconv()._preprocess_data(*counts, "pearson")
        b = preprocess_expression(*counts, "pearson")
        np.testing.assert_array_equal(a[0], b[0])


class TestGeneSelection:
    def test_hvg_sparse_matches_dense(self, counts):
        Y, _ = counts
        np.testing.assert_array_equal(
            genes.select_hvg(Y, n_top=10), genes.select_hvg(sparse.csr_matrix(Y), n_top=10)
        )

    def test_hvg_single_spot(self, counts):
        idx = genes.select_hvg(counts[0][:1], n_top=5)
        assert len(idx) == 5 and np.all(np.diff(idx) > 0)

    @pytest.mark.parametrize("method", ["diff", "ratio", "specificity"])
    def test_markers_are_specific(self, method):
        X = np.ones((3, 30))
        for k in range(3):
            X[k, 10 * k:10 * k + 5] = 50.0
        idx, assign = genes.select_markers(X, n_markers=5, method=method)
        assert len(assign) == 15
        for k in range(3):
            assert set(range(10 * k, 10 * k + 5)) <= set(idx.tolist())

    def test_markers_bad_method(self):
        with pytest.raises(ValueError):
            genes.select_markers(np.ones((2, 5)), method="foo")

    def test_leverage_scores_sum_to_one(self, counts):
        lev = genes.compute_leverage_scores(counts[1])
        assert lev.shape == (50,) and np.all(lev >= 0)
        np.testing.assert_allclose(lev.sum(), 1.0, atol=1e-5)

    def test_informative_genes_union(self, counts):
        Y, X = counts
        idx, lev = genes.select_informative_genes(Y, X, n_hvg=10, n_markers_per_type=3)
        assert idx.dtype == np.intp and np.all(np.diff(idx) > 0)
        assert len(lev) == len(idx)


class TestMetrics:
    def test_perfect_prediction(self):
        P = np.random.RandomState(0).dirichlet(np.ones(4), size=20)
        assert metrics.compute_rmse(P, P) == 0
        assert metrics.compute_mae(P, P) == 0
        np.testing.assert_allclose(metrics.compute_correlation(P, P), 1.0)
        np.testing.assert_allclose(metrics.compute_jsd(P, P), 0.0, atol=1e-12)
        res = metrics.evaluate_deconvolution(P, P, cell_type_names=list("abcd"))
        assert set(res["per_cell_type"]) == set("abcd")
        np.testing.assert_allclose(res["overall"]["spearman"], 1.0)

    def test_per_cell_type_shapes_and_constant_input(self):
        rng = np.random.RandomState(1)
        P, Q = rng.dirichlet(np.ones(3), size=10), rng.dirichlet(np.ones(3), size=10)
        assert metrics.compute_rmse(P, Q, per_cell_type=True).shape == (3,)
        Q[:, 0] = 0.3  # constant column -> correlation defined as 0
        assert metrics.compute_correlation(P, Q, per_cell_type=True)[0] == 0.0

    def test_jsd_bounded(self):
        P = np.array([[1.0, 0.0], [0.5, 0.5]])
        Q = np.array([[0.0, 1.0], [0.5, 0.5]])
        jsd = metrics.compute_jsd(P, Q)
        assert jsd[0] <= np.log(2) + 1e-9 and abs(jsd[1]) < 1e-12

    def test_rare_cell_detection(self):
        true = np.array([[0.02, 0.98], [0.0, 1.0]])
        pred = np.array([[0.03, 0.97], [0.2, 0.8]])
        p, r, f = metrics.compute_rare_cell_detection(pred, true, threshold=0.05)
        np.testing.assert_allclose([p, r], [0.5, 1.0], atol=1e-8)
        assert np.isnan(metrics.compute_rare_cell_detection(pred, np.ones_like(true))[0])


anndata = pytest.importorskip("anndata")
pd = pytest.importorskip("pandas")


def _adata(X, obs=None, var_names=None, obsm=None):
    a = anndata.AnnData(X=X, obs=obs)
    if var_names is not None:
        a.var_names = var_names
    for k, v in (obsm or {}).items():
        a.obsm[k] = v
    return a


class TestIO:
    def test_load_spatial_coordinate_fallbacks(self):
        from flashdeconv.io import load_spatial_data

        Y = np.ones((4, 3))
        xy = np.arange(8.0).reshape(4, 2)
        a = _adata(Y, obsm={"X_spatial": xy})
        np.testing.assert_array_equal(load_spatial_data(a)[1], xy)
        a = _adata(Y, obs=pd.DataFrame({"x": xy[:, 0], "y": xy[:, 1]},
                                       index=[f"s{i}" for i in range(4)]))
        np.testing.assert_array_equal(load_spatial_data(a)[1], xy)
        a = _adata(Y)
        with pytest.raises(ValueError, match="spatial coordinates"):
            load_spatial_data(a)

    def test_load_reference_mean_sum_layer(self):
        from flashdeconv.io import load_reference

        X = sparse.csr_matrix(np.array([[1.0, 0], [3, 2], [5, 5]]))
        obs = pd.DataFrame({"ct": ["b", "a", "b"]}, index=["c0", "c1", "c2"])
        a = _adata(X, obs=obs, var_names=["g0", "g1"])
        a.layers["counts"] = X * 2
        M, names, gn = load_reference(a, "ct")
        np.testing.assert_array_equal(names, ["a", "b"])
        np.testing.assert_allclose(M, [[3, 2], [3, 2.5]])
        S, _, _ = load_reference(a, "ct", method="sum", layer="counts")
        np.testing.assert_allclose(S, [[6, 4], [12, 10]])
        with pytest.raises(ValueError, match="aggregation"):
            load_reference(a, "ct", method="median")
        with pytest.raises(ValueError, match="not found"):
            load_reference(a, "celltype")

    def test_load_reference_missing_labels(self):
        from flashdeconv.io import load_reference

        obs = pd.DataFrame({"ct": pd.Categorical(["a", None, "b"])}, index=["c0", "c1", "c2"])
        a = _adata(np.ones((3, 2)), obs=obs)
        with pytest.raises(ValueError, match="1 cells have no 'ct' label"):
            load_reference(a, "ct")

    def test_align_genes(self):
        from flashdeconv.io import align_genes

        Y = sparse.csr_matrix(np.arange(12.0).reshape(3, 4))
        X = np.arange(6.0).reshape(2, 3)
        Ya, Xa, common = align_genes(Y, X, np.array(["d", "b", "a", "c"]),
                                     np.array(["c", "a", "e"]))
        np.testing.assert_array_equal(common, ["a", "c"])
        np.testing.assert_array_equal(Ya.toarray(), Y.toarray()[:, [2, 3]])
        np.testing.assert_array_equal(Xa, X[:, [1, 0]])
        with pytest.raises(ValueError, match="No common genes"):
            align_genes(Y, X, np.array(["x"] * 4), np.array(["y"] * 3))

    def test_deconvolve_forwards_solver_settings(self):
        import flashdeconv as fd

        rng = np.random.RandomState(0)
        G = 40
        genes_ = [f"g{i}" for i in range(G)]
        ref = _adata(rng.poisson(3, (60, G)).astype(np.float32),
                     obs=pd.DataFrame({"cell_type": np.repeat(["A", "B", "C"], 20)},
                                      index=[f"c{i}" for i in range(60)]),
                     var_names=genes_)
        st = _adata(sparse.csr_matrix(rng.poisson(5, (25, G)).astype(np.float32)),
                    var_names=genes_, obsm={"spatial": rng.rand(25, 2)})
        st.obs_names = [f"s{i}" for i in range(25)]
        out = fd.tl.deconvolve(st, ref, max_iter=3, tol=1e-12, copy=True)
        params = out.uns["flashdeconv_params"]
        assert params["max_iter"] == 3 and params["n_iterations"] == 3
        assert "flashdeconv" not in st.obsm  # copy=True leaves the input untouched
        np.testing.assert_allclose(out.obsm["flashdeconv"].sum(1), 1.0)
