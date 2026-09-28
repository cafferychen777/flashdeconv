"""
Regression tests pinning numerical outputs on a small deterministic fixture.

The pinned values were produced by the v0.2.0 package that generated the
manuscript results. A failure means a change altered the numerical output
of the estimator; if that is intended, update the values deliberately.
Tolerances allow for BLAS / platform floating-point differences only.
"""

import warnings

import numpy as np
import pytest
from scipy import sparse

from flashdeconv import FlashDeconv, reference_fit_scores

RTOL = 1e-6
ATOL = 1e-9


def _fixture():
    # Legacy RandomState: its streams are frozen across numpy versions.
    rng = np.random.RandomState(12345)
    K, G, side = 4, 300, 10
    X = rng.gamma(0.5, 2.0, size=(K, G))
    for k in range(K):
        X[k, rng.choice(G, 15, replace=False)] *= 20
    coords = np.array([(i % side, i // side) for i in range(side * side)], dtype=float)
    P = rng.dirichlet(np.full(K, 0.5), size=side * side)
    mu = P @ X
    mu = mu / mu.sum(1, keepdims=True) * rng.uniform(500, 2000, size=(side * side, 1))
    Y = rng.poisson(mu).astype(np.float64)
    return Y, X, coords, P


PARAMS = dict(n_hvg=150, n_markers_per_type=20, sketch_dim=64)


@pytest.fixture(scope="module")
def data():
    return _fixture()


@pytest.fixture(scope="module")
def expected_model(data):
    Y, X, coords, _ = data
    return FlashDeconv(**PARAMS).fit(sparse.csr_matrix(Y), X, coords)


def test_gene_selection_pinned(expected_model):
    assert len(expected_model.gene_idx_) == 172
    assert int(expected_model.gene_idx_.sum()) == 25640


def test_default_path_pinned(expected_model):
    m = expected_model
    assert abs(m.info_["n_iterations"] - 19) <= 1
    np.testing.assert_allclose(m.lambda_used_, 0.4165116115284268, rtol=RTOL)
    np.testing.assert_allclose(m.info_["final_objective"], 3883.15349442539, rtol=RTOL)
    np.testing.assert_allclose(
        m.proportions_.mean(0),
        [0.18483184909099534, 0.3078592968350551, 0.2502890855638255, 0.2570197685101242],
        rtol=RTOL, atol=ATOL,
    )
    np.testing.assert_allclose(
        m.proportions_[0],
        [0.0, 0.541732306682379, 0.39744369652669154, 0.06082399679092942],
        rtol=RTOL, atol=ATOL,
    )
    np.testing.assert_allclose(m.gene_weights_.sum(), 84.16469334153982, rtol=1e-10)
    np.testing.assert_allclose(m.gene_weights_.min(), 0.30551721287730166, rtol=1e-10)
    np.testing.assert_allclose(m.gene_weights_.max(), 0.975180413628658, rtol=1e-10)


def test_legacy_countsketch_pinned(data):
    Y, X, coords, _ = data
    m = FlashDeconv(gene_weighting="countsketch", random_state=0, **PARAMS).fit(
        sparse.csr_matrix(Y), X, coords
    )
    assert m.gene_weights_ is None
    assert m.X_repr_.shape == (4, 64)
    assert abs(m.info_["n_iterations"] - 22) <= 1
    np.testing.assert_allclose(m.lambda_used_, 1.3032669003233368, rtol=RTOL)
    np.testing.assert_allclose(m.info_["final_objective"], 12069.236174651196, rtol=RTOL)
    np.testing.assert_allclose(
        m.proportions_.mean(0),
        [0.17329584768405376, 0.3508521118554823, 0.2124579394404577, 0.26339410102000604],
        rtol=RTOL, atol=ATOL,
    )
    np.testing.assert_allclose(
        m.proportions_[0],
        [0.03326953215463459, 0.5455128311304942, 0.3179292242409505, 0.10328841247392077],
        rtol=RTOL, atol=ATOL,
    )


def test_uncertainty_pinned(expected_model):
    uq = expected_model.compute_uncertainty()
    np.testing.assert_allclose(np.nanmean(uq["se_prop"]), 0.04372690396696026, rtol=1e-5)
    np.testing.assert_allclose(uq["mean_ci_width"], 0.1419375688831919, rtol=1e-5)
    assert abs(int(uq["detected"].sum()) - 299) <= 1


def test_bootstrap_pinned(expected_model):
    boot = expected_model.bootstrap_uncertainty(n_bootstrap=3, seed=0)
    assert boot["n_converged"] == 3
    np.testing.assert_allclose(
        boot["boot_mean"].mean(0),
        [0.17881909, 0.311783, 0.24631207, 0.26308593],
        rtol=1e-5,
    )


def test_refcheck_pinned(expected_model):
    s = reference_fit_scores(expected_model)
    np.testing.assert_allclose(
        s["z_raw"][:3], [1.2315759066979923, -2.26143907481568, -0.6942997459353999],
        rtol=1e-6,
    )
    assert s["null"]["score"]["method"] == "left_half"
    np.testing.assert_allclose(s["null"]["score"]["center"], -0.2689992645342816, rtol=1e-6)
    np.testing.assert_allclose(s["null"]["score"]["scale"], 0.8700260184058313, rtol=1e-6)
    np.testing.assert_allclose(s["null"]["score_pooled"]["center"], -0.8606888499204057, rtol=1e-6)
    np.testing.assert_allclose(s["null"]["score_pooled"]["scale"], 0.6675619273335311, rtol=1e-6)


def test_recovers_ground_truth(data, expected_model):
    P = data[3]
    r = np.corrcoef(expected_model.proportions_.ravel(), P.ravel())[0, 1]
    assert r > 0.9


class TestDeterminism:
    """Default output is bit-identical across runs, random_state and input format."""

    def test_repeat_and_random_state(self, data, expected_model):
        Y, X, coords, _ = data
        for rs in (0, 1, None):
            m = FlashDeconv(random_state=rs, **PARAMS).fit(sparse.csr_matrix(Y), X, coords)
            np.testing.assert_array_equal(m.beta_, expected_model.beta_)

    def test_sparse_formats_identical(self, data, expected_model):
        Y, X, coords, _ = data
        for fmt in ("csc", "coo", "lil"):
            m = FlashDeconv(**PARAMS).fit(sparse.csr_matrix(Y).asformat(fmt), X, coords)
            np.testing.assert_array_equal(m.beta_, expected_model.beta_)

    def test_dense_matches_sparse(self, data, expected_model):
        Y, X, coords, _ = data
        m = FlashDeconv(**PARAMS).fit(Y, X, coords)
        np.testing.assert_array_equal(m.gene_idx_, expected_model.gene_idx_)
        np.testing.assert_allclose(m.proportions_, expected_model.proportions_, atol=1e-10)

    def test_integer_counts_match_float(self, data):
        Y, X, coords, _ = data
        a = FlashDeconv(**PARAMS).fit(Y, X, coords)
        b = FlashDeconv(**PARAMS).fit(Y.astype(np.int64), X, coords)
        np.testing.assert_array_equal(a.beta_, b.beta_)

    def test_legacy_random_state_instance_equals_int_seed(self, data):
        Y, X, coords, _ = data
        a = FlashDeconv(gene_weighting="countsketch", random_state=5, **PARAMS).fit(Y, X, coords)
        b = FlashDeconv(
            gene_weighting="countsketch",
            random_state=np.random.RandomState(5),
            **PARAMS,
        ).fit(Y, X, coords)
        np.testing.assert_array_equal(a.beta_, b.beta_)


def test_bootstrap_reuses_fitted_sketch(data, monkeypatch):
    """With random_state=None the bootstrap must use the fitted hash functions."""
    import flashdeconv.core.sketching as sk

    Y, X, coords, _ = data
    m = FlashDeconv(gene_weighting="countsketch", random_state=None, **PARAMS).fit(
        Y, X, coords
    )
    assert m._sketch_matrix is not None

    def _no_resketch(*args, **kwargs):
        raise AssertionError("bootstrap drew a new CountSketch")

    monkeypatch.setattr(sk, "sketch_data", _no_resketch)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        boot = m.bootstrap_uncertainty(n_bootstrap=3, seed=1)
    assert boot["boot_mean"].shape == m.proportions_.shape
    assert np.all(np.isfinite(boot["boot_mean"]))
