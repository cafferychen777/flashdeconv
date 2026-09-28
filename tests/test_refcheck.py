"""Tests for the opt-in incomplete-reference diagnostic (flashdeconv.core.refcheck)."""

import numpy as np
import pytest
from scipy import sparse

from flashdeconv import FlashDeconv
from flashdeconv.core.refcheck import (
    reference_fit_scores,
    suggest_missing_types,
    unexplained_genes,
)


def _simulate(seed=0, n_side=30, n_genes=300, K=4, depth=300):
    """K+1 cell types with disjoint marker blocks; the last type fills a square region."""
    rng = np.random.default_rng(seed)
    base = rng.gamma(2.0, 1.0, size=(K + 1, n_genes))
    block = n_genes // (K + 1)
    for k in range(K + 1):
        base[k, k * block:(k + 1) * block] *= 20.0
    xs, ys = np.meshgrid(np.arange(n_side), np.arange(n_side))
    coords = np.c_[xs.ravel(), ys.ravel()].astype(float)
    N = coords.shape[0]
    props = rng.dirichlet(np.ones(K), size=N)
    region = (coords[:, 0] < n_side // 3) & (coords[:, 1] < n_side // 3)
    full = np.zeros((N, K + 1))
    full[:, :K] = props
    full[region] = 0.3 * full[region]
    full[region, K] = 0.7
    prof = base / base.sum(axis=1, keepdims=True)
    mu = depth * (full @ prof)
    Y = sparse.csr_matrix(rng.poisson(mu).astype(np.float64))
    genes = np.array([f"g{j}" for j in range(n_genes)])
    return Y, base, coords, region, genes, block


def _auroc(labels, scores):
    from scipy.stats import rankdata

    r = rankdata(scores)
    n1 = labels.sum()
    n0 = len(labels) - n1
    return (r[labels].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)


def _fit(Y, X, coords):
    m = FlashDeconv(sketch_dim=128, n_hvg=200, n_markers_per_type=20, max_iter=50)
    m.fit(Y, X, coords)
    return m


def test_complete_reference_flags_near_nominal():
    Y, base, coords, region, genes, _ = _simulate()
    m = _fit(Y, base, coords)
    s = reference_fit_scores(m)
    assert s["score"].shape == (Y.shape[0],)
    assert np.all(np.isfinite(s["score"]))
    # empirical-null calibration: roughly the nominal 5% upper tail
    assert 0.0 <= s["flag"].mean() < 0.15
    assert np.allclose(s["mixture_weights"].sum(axis=1), 1.0)
    np.testing.assert_allclose(s["n_umi"], np.asarray(Y[:, m.gene_idx_].sum(axis=1)).ravel())


def test_missing_type_is_flagged_and_named():
    Y, base, coords, region, genes, block = _simulate()
    K = base.shape[0] - 1
    m = _fit(Y, base[:K], coords)
    s = reference_fit_scores(m)
    assert _auroc(region, s["score"]) > 0.95
    assert _auroc(region, s["score_pooled"]) > 0.95
    assert s["flag_pooled"][region].mean() > 0.8
    assert s["flag_pooled"][~region].mean() < 0.15

    g = unexplained_genes(m, s["flag_pooled"], gene_names=genes)
    top = set(g["gene"][:20])
    missing_markers = set(genes[K * block:(K + 1) * block])
    assert len(top & missing_markers) >= 15

    # genes below min_count are scored nan and sorted last
    g_hi = unexplained_genes(m, s["flag_pooled"], gene_names=genes, min_count=10 ** 9)
    assert np.all(np.isnan(g_hi["score"]))
    finite = np.isfinite(g["score"])
    assert not np.any(finite[np.argmin(finite):]) or finite.all()
    assert np.all(np.diff(g["score"][finite]) <= 0)

    r = suggest_missing_types(g, base, genes, [f"t{k}" for k in range(K + 1)], n_markers=20)
    assert r["type"][0] == f"t{K}"


def test_existing_fit_unchanged_and_inputs_validated():
    Y, base, coords, region, genes, _ = _simulate(seed=1)
    m = _fit(Y, base, coords)
    before = m.proportions_.copy()
    reference_fit_scores(m)
    np.testing.assert_array_equal(before, m.proportions_)
    with pytest.raises(RuntimeError):
        reference_fit_scores(FlashDeconv())
    with pytest.raises(ValueError):
        reference_fit_scores(m, eta=0.0)
    with pytest.raises(ValueError):
        unexplained_genes(m, np.zeros(Y.shape[0], dtype=bool))
    no_pool = reference_fit_scores(m, pool=False)
    assert "score_pooled" not in no_pool


def test_null_methods_rank_identically_and_report_params():
    Y, base, coords, region, genes, _ = _simulate(seed=2)
    m = _fit(Y, base[:-1], coords)
    out = {meth: reference_fit_scores(m, null=meth) for meth in ("auto", "left_half", "central")}
    for s in out.values():
        info = s["null"]["score"]
        assert info["scale"] > 0 and np.isfinite(info["center"])
        assert "score_pooled" in s["null"]
        # the null only moves the threshold: calibrated scores are an increasing map of z_raw
        o = np.argsort(s["z_raw"], kind="stable")
        assert np.all(np.diff(s["score"][o]) >= -1e-9)
    a = out["auto"]["null"]["score"]
    assert a["scale"] == min(out["left_half"]["null"]["score"]["scale"],
                             out["central"]["null"]["score"]["scale"])
    with pytest.raises(ValueError):
        reference_fit_scores(m, null="median")


def test_auto_null_robust_to_lower_shoulder():
    from flashdeconv.core.refcheck import _empirical_null

    rng = np.random.default_rng(0)
    null = rng.normal(0.0, 1.0, 70000)
    # a block of bins explained better than the multinomial expectation (heavy lower shoulder)
    shoulder = rng.normal(-8.0, 3.0, 30000)
    z = np.r_[null, shoulder]
    thr = 1.6448536269514722
    for meth in ("auto", "left_half", "central"):
        s, _ = _empirical_null(rng.normal(2.0, 3.0, 50000), meth)
        assert abs((s > thr).mean() - 0.05) < 0.01  # Gaussian bulk: all methods nominal
    s_lh, _ = _empirical_null(z, "left_half")
    s_auto, info = _empirical_null(z, "auto")
    assert (s_lh[:70000] > thr).mean() < 0.01  # left-half scale inflated: almost no flags
    assert info["method"] == "central"
    assert abs((s_auto[:70000] > thr).mean() - 0.05) < 0.015
