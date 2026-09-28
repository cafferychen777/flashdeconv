"""Tests for the expected leverage-weighted CountSketch gene representation."""

import numpy as np
import pytest
from scipy import sparse

from flashdeconv import FlashDeconv
from flashdeconv.core.deconv import FlashDeconv as _FD
from flashdeconv.core.sketching import (
    apply_gene_weights,
    build_countsketch_matrix,
    expected_countsketch_weights,
    sketch_data,
)
from flashdeconv.core.solver import bcd_solve, normalize_proportions
from flashdeconv.core.spatial import auto_tune_lambda
from flashdeconv.utils.genes import select_informative_genes
from flashdeconv.utils.graph import coords_to_adjacency

from tests.test_integration import generate_synthetic_data


def _skewed_leverage(g, seed=0):
    rng = np.random.default_rng(seed)
    lev = rng.gamma(0.5, size=g)
    lev[:5] *= 200.0  # a few high-leverage genes hit the upper clip
    return lev


@pytest.fixture(scope="module")
def data():
    return generate_synthetic_data(
        n_spots=150, n_genes=600, n_cell_types=5, noise_level=0.1, random_state=3
    )


class TestExpectedWeights:
    """Properties of expected_countsketch_weights."""

    @pytest.mark.parametrize("g,d", [(50, 512), (400, 64), (2000, 512)])
    def test_positive_finite_length(self, g, d):
        w = expected_countsketch_weights(_skewed_leverage(g), g, sketch_dim=d)
        assert w.shape == (g,)
        assert np.all(np.isfinite(w))
        assert np.all(w > 0)
        # A column-normalized weight a_j / ||bucket|| never exceeds 1.
        assert np.all(w <= 1.0 + 1e-8)

    def test_uniform_leverage_gives_equal_weights(self):
        g = 300
        w1 = expected_countsketch_weights(None, g, sketch_dim=64)
        w2 = expected_countsketch_weights(np.ones(g), g, sketch_dim=64)
        np.testing.assert_allclose(w1, w1[0], rtol=1e-12)
        np.testing.assert_allclose(w1, w2, rtol=1e-12)

    def test_matches_monte_carlo_hashing(self):
        """Formula equals the average over many random CountSketch hashings."""
        g, d, n_draws = 300, 32, 3000
        lev = _skewed_leverage(g, seed=1)
        w = expected_countsketch_weights(lev, g, sketch_dim=d)
        scale = np.sqrt(g / d)  # global factor inside build_countsketch_matrix
        mc = np.zeros(g)
        for s in range(n_draws):
            Om = build_countsketch_matrix(g, d, lev, random_state=s)
            mc += np.asarray(abs(Om).sum(axis=1)).ravel() / scale
        mc /= n_draws
        assert np.corrcoef(w, mc)[0, 1] > 0.999
        np.testing.assert_allclose(w, mc, rtol=0.03)

    def test_single_bucket_is_exact(self):
        g = 40
        lev = _skewed_leverage(g)
        w = expected_countsketch_weights(lev, g, sketch_dim=1)
        p = lev / (lev.sum() + 1e-10)
        a = np.clip(np.sqrt(p * g + 1e-10), 0.1, 10.0)
        np.testing.assert_allclose(w, a / np.linalg.norm(a), rtol=1e-12)

    def test_blocked_quadrature_matches_unblocked(self):
        g = 500
        lev = _skewed_leverage(g, seed=2)
        full = expected_countsketch_weights(lev, g, sketch_dim=128)
        blocked = expected_countsketch_weights(
            lev, g, sketch_dim=128, max_block_elements=g * 97
        )
        np.testing.assert_allclose(blocked, full, rtol=1e-12)

    def test_bad_leverage_shape_raises(self):
        with pytest.raises(ValueError):
            expected_countsketch_weights(np.ones(5), 6)


class TestGeneWeightingFit:
    """End-to-end behaviour of the gene_weighting option."""

    def test_default_is_expected(self):
        assert FlashDeconv().gene_weighting == "expected"

    def test_invalid_option_raises(self):
        with pytest.raises(ValueError):
            FlashDeconv(gene_weighting="random")

    def test_deterministic_across_random_state(self, data):
        Y, X, coords, _ = data
        p = [FlashDeconv(random_state=rs).fit(Y, X, coords).proportions_
             for rs in (0, 1, 12345, None)]
        for q in p[1:]:
            np.testing.assert_array_equal(q, p[0])

    def test_matches_manual_weighted_pipeline(self, data):
        """Default = selected genes scaled by the expected weights, then the
        unchanged auto-lambda / BCD pipeline."""
        Y, X, coords, _ = data
        m = FlashDeconv().fit(Y, X, coords)

        idx, lev = select_informative_genes(Y, X, n_hvg=2000, n_markers_per_type=50)
        Yt, Xt = _FD._preprocess_data(None, Y[:, idx], X[:, idx], "log_cpm")
        w = expected_countsketch_weights(lev, len(idx), sketch_dim=512)
        Yw, Xw = apply_gene_weights(Yt, Xt, w)
        A = coords_to_adjacency(coords, method="knn", k=6)
        lam = auto_tune_lambda(Yw, Xw, A)
        beta, _ = bcd_solve(Yw, Xw, A, lambda_=lam, rho=0.01, max_iter=m.max_iter, tol=m.tol)

        np.testing.assert_array_equal(m.gene_weights_, w)
        np.testing.assert_array_equal(m.proportions_, normalize_proportions(beta))
        assert m.X_sketch_.shape[1] == len(idx)

    def test_legacy_countsketch_matches_previous_pipeline(self, data):
        """gene_weighting='countsketch' reproduces the <=0.1.6 pipeline bit-for-bit."""
        Y, X, coords, _ = data
        m = FlashDeconv(gene_weighting="countsketch", random_state=7).fit(Y, X, coords)

        idx, lev = select_informative_genes(Y, X, n_hvg=2000, n_markers_per_type=50)
        Yt, Xt = _FD._preprocess_data(None, Y[:, idx], X[:, idx], "log_cpm")
        Ys, Xs, _ = sketch_data(Yt, Xt, sketch_dim=512, leverage_scores=lev,
                                random_state=7)
        A = coords_to_adjacency(coords, method="knn", k=6)
        lam = auto_tune_lambda(Ys, Xs, A)
        beta, _ = bcd_solve(Ys, Xs, A, lambda_=lam, rho=0.01, max_iter=m.max_iter, tol=m.tol)

        assert m.gene_weights_ is None
        np.testing.assert_array_equal(m.beta_, beta)
        assert m.X_sketch_.shape[1] == 512

    def test_legacy_depends_on_random_state(self, data):
        Y, X, coords, _ = data
        p0 = FlashDeconv(gene_weighting="countsketch", random_state=0).fit(Y, X, coords)
        p1 = FlashDeconv(gene_weighting="countsketch", random_state=1).fit(Y, X, coords)
        assert not np.array_equal(p0.proportions_, p1.proportions_)

    def test_sparse_matches_dense(self, data):
        Y, X, coords, _ = data
        pd_ = FlashDeconv().fit(Y, X, coords).proportions_
        ps = FlashDeconv().fit(sparse.csr_matrix(Y), X, coords).proportions_
        np.testing.assert_allclose(ps, pd_, atol=1e-10)

    def test_recovery_accuracy(self, data):
        Y, X, coords, beta_true = data
        p = FlashDeconv().fit(Y, X, coords).proportions_
        r = np.corrcoef(p.ravel(), beta_true.ravel())[0, 1]
        assert r > 0.8
        np.testing.assert_allclose(p.sum(axis=1), 1.0)
        assert np.all(p >= 0)

    @pytest.mark.parametrize("to_sparse", [False, True])
    def test_uncertainty(self, data, to_sparse):
        Y, X, coords, _ = data
        Yin = sparse.csr_matrix(Y) if to_sparse else Y
        m = FlashDeconv().fit(Yin, X, coords)
        uq = m.compute_uncertainty()
        for k in ("var_prop", "ci_lower", "ci_upper", "residual_ss", "entropy"):
            assert np.all(np.isfinite(uq[k]))
        assert np.all(uq["ci_lower"] <= m.proportions_ + 1e-12)
        assert np.all(uq["ci_upper"] >= m.proportions_ - 1e-12)
        boot = m.bootstrap_uncertainty(n_bootstrap=3)
        assert boot["n_converged"] == 3
        assert boot["boot_mean"].shape == m.proportions_.shape
        assert np.all(np.isfinite(boot["boot_std"]))

    def test_sparse_uncertainty_matches_dense(self, data):
        Y, X, coords, _ = data
        ud = FlashDeconv().fit(Y, X, coords).compute_uncertainty()
        us = FlashDeconv().fit(sparse.csr_matrix(Y), X, coords).compute_uncertainty()
        np.testing.assert_allclose(us["residual_ss"], ud["residual_ss"], rtol=1e-8)
        np.testing.assert_allclose(us["var_prop"], ud["var_prop"], rtol=1e-6, atol=1e-14)

    def test_summary_reports_representation(self, data):
        Y, X, coords, _ = data
        s = FlashDeconv().fit(Y, X, coords).summary()
        assert s["gene_weighting"] == "expected"
        assert s["representation_dim"] == s["n_genes_used"]


def test_deconvolve_forwards_gene_weighting():
    anndata = pytest.importorskip("anndata")
    import pandas as pd
    from flashdeconv.tl import deconvolve

    rng = np.random.default_rng(0)
    n_spots, n_genes, n_cells, n_types = 60, 200, 120, 3
    adata_st = anndata.AnnData(
        X=rng.poisson(5, (n_spots, n_genes)).astype(float),
        obsm={"spatial": rng.random((n_spots, 2)) * 10},
    )
    adata_st.var_names = [f"g{i}" for i in range(n_genes)]
    labels = np.repeat([f"t{i}" for i in range(n_types)], n_cells // n_types)
    adata_ref = anndata.AnnData(
        X=rng.poisson(5, (n_cells, n_genes)).astype(float),
        obs=pd.DataFrame({"cell_type": labels}),
    )
    adata_ref.var_names = adata_st.var_names.copy()

    deconvolve(adata_st, adata_ref)
    assert adata_st.uns["flashdeconv_params"]["gene_weighting"] == "expected"
    deconvolve(adata_st, adata_ref, gene_weighting="countsketch", sketch_dim=32,
               key_added="fd_legacy")
    assert adata_st.uns["fd_legacy_params"]["gene_weighting"] == "countsketch"
