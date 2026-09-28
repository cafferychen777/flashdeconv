"""Tests for model-based uncertainty quantification."""

import warnings

import numpy as np
import pytest
from scipy import sparse
from scipy.stats import norm

from flashdeconv import FlashDeconv
from flashdeconv.core.solver import bcd_solve, normalize_proportions
from flashdeconv.core.uncertainty import (
    bh_detection,
    compute_active_set_uncertainty,
)
from flashdeconv.utils.graph import coords_to_adjacency


def _explicit_spot(y, Xs, b, lam, deg, rho, z, sandwich=True, nsum=None):
    """Reference implementation with explicit matrix inverses.

    ``nsum`` is the sum of the neighbours' abundances; at b_k = 0 the
    gradient of 0.5 lam sum_j ||b - b_j||^2 in k is -lam * nsum_k.
    """
    if nsum is None:
        nsum = np.zeros_like(b)
    K, p = Xs.shape
    e = y - b @ Xs
    G = Xs @ Xs.T
    Hf = G + lam * deg * np.eye(K)
    M = (Xs * e ** 2) @ Xs.T
    act = np.where(b > 0)[0]
    S = b[act].sum()
    HA_inv = np.linalg.inv(Hf[np.ix_(act, act)])
    if sandwich:
        VA = HA_inv @ M[np.ix_(act, act)] @ HA_inv
    else:
        VA = (e @ e) / max(p - len(act), 1) * HA_inv
    pA = b[act] / S
    J = (np.eye(len(act)) - pA[:, None]) / S  # d p_u / d b_v
    se = np.full(K, np.nan)
    up = np.full(K, np.nan)
    se[act] = np.sqrt(np.diag(J @ VA @ J.T))
    for k in np.where(b <= 0)[0]:
        idx = np.append(act, k)
        HB_inv = np.linalg.inv(Hf[np.ix_(idx, idx)])
        gB = np.zeros(len(idx))
        gB[-1] = Xs[k] @ e + lam * nsum[k] - rho
        bk = (HB_inv @ gB)[-1]
        if sandwich:
            vk = (HB_inv @ M[np.ix_(idx, idx)] @ HB_inv)[-1, -1]
        else:
            vk = (e @ e) / max(p - len(idx), 1) * HB_inv[-1, -1]
        sk = np.sqrt(vk)
        se[k] = sk / S
        u = max(bk + z * sk, 0.0)
        up[k] = u / (S + u)
    return se, up, VA, HA_inv


def _simulate(n=400, p=300, K=5, seed=0, zero_frac=0.0, collinear=True):
    rng = np.random.default_rng(seed)
    Xs = rng.gamma(2.0, 1.0, (K, p))
    if collinear:
        Xs[1] = 0.7 * Xs[0] + 0.3 * Xs[1]  # strongly correlated signatures
    beta = rng.dirichlet(np.full(K, 3.0), n) * 5.0
    if zero_frac > 0:
        beta[rng.random((n, K)) < zero_frac] = 0.0
        beta[beta.sum(1) == 0, 0] = 5.0
    mu = beta @ Xs
    sd = 0.3 * np.sqrt(mu)  # heteroscedastic noise
    return rng, Xs, beta, mu, sd


def _fit_nnls(Y, Xs, lam=0.0, rho=0.0, return_graph=False):
    coords = np.column_stack([np.arange(Y.shape[0]) % 20, np.arange(Y.shape[0]) // 20])
    A = coords_to_adjacency(coords.astype(float), method="knn", k=4)
    beta, info = bcd_solve(Y, Xs, A, lambda_=lam, rho=rho, max_iter=5000, tol=1e-10)
    return (beta, info, A) if return_graph else (beta, info)


class TestSandwichFormula:
    @pytest.mark.parametrize("variance", ["sandwich", "model"])
    def test_matches_explicit_inverse(self, variance):
        rng, Xs, _, mu, sd = _simulate(n=60, zero_frac=0.3, seed=1)
        Y = mu + sd * rng.standard_normal(mu.shape)
        beta, info, A = _fit_nnls(Y, Xs, lam=0.5, rho=0.0, return_graph=True)
        rho_abs = 0.0
        deg = info["n_neighbors"]
        nsum = np.asarray(A @ beta)
        res = compute_active_set_uncertainty(
            Y, Xs, beta, lambda_=0.5, rho_abs=rho_abs, n_neighbors=deg,
            variance=variance, neighbor_sum=nsum,
        )
        z = norm.ppf(0.975)
        for i in range(10):
            se, up, _, _ = _explicit_spot(Y[i], Xs, beta[i], 0.5, deg[i], rho_abs, z,
                                          sandwich=(variance == "sandwich"),
                                          nsum=nsum[i])
            np.testing.assert_allclose(res["se_prop"][i], se, rtol=1e-7)
            inact = beta[i] <= 0
            np.testing.assert_allclose(res["ci_upper"][i, inact], up[inact], rtol=1e-7)
            assert np.all(res["ci_lower"][i, inact] == 0)

    def test_full_inverse_not_diagonal(self):
        # With collinear signatures, the model-based variance must use
        # diag(H^-1), which exceeds 1/diag(H).
        rng, Xs, beta_true, mu, sd = _simulate(n=40, seed=2)
        Y = mu + sd * rng.standard_normal(mu.shape)
        beta, info = _fit_nnls(Y, Xs)
        i = int(np.where((beta > 0).all(1))[0][0])
        _, _, VA, HA_inv = _explicit_spot(Y[i], Xs, beta[i], 0.0, 0, 0.0, 1.96,
                                          sandwich=False)
        G = Xs @ Xs.T
        assert HA_inv[0, 0] > 2.0 / G[0, 0]
        e = Y[i] - beta[i] @ Xs
        s2 = e @ e / (Xs.shape[1] - Xs.shape[0])
        np.testing.assert_allclose(np.diag(VA), s2 * np.diag(np.linalg.inv(G)), rtol=1e-10)


class TestSpatialTermForZeroTypes:
    """Newton step for types estimated at zero includes the spatial gradient."""

    def test_matches_finite_difference_gradient(self):
        from flashdeconv.core.solver import compute_objective
        from flashdeconv.core.spatial import compute_laplacian

        rng, Xs, _, mu, sd = _simulate(n=80, p=120, K=4, seed=5, zero_frac=0.4)
        Y = mu + 2.0 * sd * rng.standard_normal(mu.shape)
        lam, rho = 2.0, 0.01
        beta, info, A = _fit_nnls(Y, Xs, lam=lam, rho=rho, return_graph=True)
        G = Xs @ Xs.T
        rho_abs = rho * np.mean(np.diag(G))
        nsum = np.asarray(A @ beta)
        deg = info["n_neighbors"]
        res = compute_active_set_uncertainty(
            Y, Xs, beta, lam, rho_abs, deg, variance="model", neighbor_sum=nsum,
        )
        L = compute_laplacian(A)
        H, YtY = Xs @ Y.T, float(np.sum(Y ** 2))

        def F(b):
            return compute_objective(b, H, G, YtY, L, lam, rho_abs)

        cand = [(i, k) for i, k in zip(*np.where(beta <= 0))
                if nsum[i, k] > 0 and beta[i].sum() > 0]
        assert len(cand) >= 5
        for i, k in cand[:10]:
            # Richardson extrapolation of the one-sided difference (F is
            # quadratic along b_ik >= 0, so this is exact up to rounding).
            h = 1e-3

            def D(t, i=i, k=k):
                b = beta.copy()
                b[i, k] = t
                return (F(b) - F(beta)) / t

            g_fd = 2 * D(h) - D(2 * h)
            act = np.where(beta[i] > 0)[0]
            Hf = G + lam * deg[i] * np.eye(G.shape[0])
            hA = Hf[act, k]
            s = Hf[k, k] - hA @ np.linalg.solve(Hf[np.ix_(act, act)], hA)
            e = Y[i] - beta[i] @ Xs
            sk = np.sqrt((e @ e) / max(Xs.shape[1] - len(act) - 1, 1) / s)
            np.testing.assert_allclose(res["z_score"][i, k], (-g_fd / s) / sk,
                                       rtol=1e-5, atol=1e-6)

    def test_zero_type_coverage_with_spatial_smoothing(self):
        # Low-abundance type present everywhere; strong smoothing zeroes it in
        # some spots although the neighbours carry it. The one-sided upper
        # bound must still cover the true proportion near the nominal level.
        side, p, K = 30, 200, 4
        xy = np.array([(i % side, i // side) for i in range(side * side)], float)
        A = coords_to_adjacency(xy, method="knn", k=4)
        cov_new, cov_old = [], []
        for seed in range(3):
            rng = np.random.default_rng(seed)
            Xs = rng.gamma(2.0, 1.0, (K, p))
            beta_t = np.tile([3.0, 2.0, 1.0, 0.1], (side * side, 1))
            beta_t *= np.exp(0.2 * np.sin(xy[:, :1] / 5.0))
            mu = beta_t @ Xs
            Y = mu + np.sqrt(mu) * rng.standard_normal(mu.shape)
            lam = np.mean(np.diag(Xs @ Xs.T)) / np.mean(np.asarray(A.sum(1)))
            beta, info = bcd_solve(Y, Xs, A, lambda_=lam, rho=0.0, max_iter=5000, tol=1e-10)
            p_true = normalize_proportions(beta_t)
            inact = beta <= 0
            for nsum, out in ((np.asarray(A @ beta), cov_new), (np.zeros_like(beta), cov_old)):
                r = compute_active_set_uncertainty(
                    Y, Xs, beta, lam, 0.0, info["n_neighbors"], variance="model",
                    neighbor_sum=nsum,
                )
                out.extend((p_true <= r["ci_upper"] + 1e-12)[inact])
        assert len(cov_new) > 50
        assert np.mean(cov_new) >= 0.8, np.mean(cov_new)
        assert np.mean(cov_old) < np.mean(cov_new) - 0.2

    def test_missing_neighbor_sum_warns(self):
        rng, Xs, _, mu, sd = _simulate(n=40, seed=6)
        beta, info = _fit_nnls(mu, Xs, lam=0.5)
        with pytest.warns(RuntimeWarning, match="neighbor_sum"):
            compute_active_set_uncertainty(mu, Xs, beta, 0.5, 0.0, info["n_neighbors"])
        with pytest.raises(ValueError, match="neighbor_sum"):
            compute_active_set_uncertainty(mu, Xs, beta, 0.5, 0.0, info["n_neighbors"],
                                           neighbor_sum=beta[:3])


class TestModelBasedCoverage:
    def test_coverage_near_nominal(self):
        rng, Xs, beta_true, mu, sd = _simulate(n=1500, p=300, K=5, seed=3)
        Y = mu + sd * rng.standard_normal(mu.shape)
        beta, info = _fit_nnls(Y, Xs)
        res = compute_active_set_uncertainty(
            Y, Xs, beta, 0.0, 0.0, info["n_neighbors"], alpha=0.05,
        )
        p_true = normalize_proportions(beta_true)
        cov = (p_true >= res["ci_lower"]) & (p_true <= res["ci_upper"])
        assert 0.90 <= cov.mean() <= 0.99, cov.mean()

    def test_inactive_upper_bounds_cover_zero_types(self):
        rng, Xs, beta_true, mu, sd = _simulate(n=1500, seed=4, zero_frac=0.4)
        Y = mu + sd * rng.standard_normal(mu.shape)
        beta, info = _fit_nnls(Y, Xs)
        res = compute_active_set_uncertainty(Y, Xs, beta, 0.0, 0.0, info["n_neighbors"])
        p_true = normalize_proportions(beta_true)
        inact = beta <= 0
        assert inact.any()
        # one-sided: lower bound 0, upper bound >= estimate and not degenerate
        assert np.all(res["ci_lower"][inact] == 0)
        assert np.all(res["ci_upper"][inact] >= 0)
        assert np.mean(res["ci_upper"][inact] > 0) > 0.5
        cov = (p_true >= res["ci_lower"] - 1e-12) & (p_true <= res["ci_upper"] + 1e-12)
        assert cov.mean() >= 0.90, cov.mean()


class TestDetection:
    def test_bh_shape_and_null_control(self):
        rng = np.random.default_rng(0)
        fdps = []
        for _ in range(20):
            pv = rng.uniform(size=(500, 8))
            pv[:50, 0] = rng.uniform(0, 1e-4, 50)  # true signals
            rej = bh_detection(pv, q=0.1)
            assert rej.shape == pv.shape and rej.dtype == bool
            null = np.ones_like(rej)
            null[:50, 0] = False
            fdps.append((rej & null).sum() / max(rej.sum(), 1))
        assert np.mean(fdps) <= 0.1 + 0.03
        # global null: almost never any rejection
        rej0 = bh_detection(np.random.default_rng(1).uniform(size=(200, 5)), q=0.1)
        assert rej0.sum() <= 2

    def test_bh_handles_nan_and_monotone_in_q(self):
        pv = np.array([[0.001, np.nan], [0.02, 0.5]])
        r1 = bh_detection(pv, q=0.01)
        r2 = bh_detection(pv, q=0.2)
        assert not r1[0, 1] and not r2[0, 1]
        assert r1.sum() <= r2.sum()

    def test_null_simulation_z_scores(self):
        # Types absent everywhere: one-sided p-values should give few calls.
        rng, Xs, beta_true, mu, sd = _simulate(n=800, K=5, seed=5, collinear=False)
        beta_true[:, 4] = 0.0
        mu = beta_true @ Xs
        sd = 0.3 * np.sqrt(np.maximum(mu, 1e-6))
        Y = mu + sd * rng.standard_normal(mu.shape)
        beta, info = _fit_nnls(Y, Xs)
        res = compute_active_set_uncertainty(Y, Xs, beta, 0.0, 0.0, info["n_neighbors"])
        rej = bh_detection(norm.sf(res["z_score"]), q=0.1)
        false_calls = rej[:, 4].sum()
        assert false_calls / max(rej.sum(), 1) <= 0.1


@pytest.fixture(scope="module")
def fitted():
    rng = np.random.default_rng(0)
    n_spots, n_genes, K = 150, 300, 4
    X = rng.gamma(2, 1, (K, n_genes)) * 5
    X[:, :60] *= rng.choice([0.1, 8.0], (K, 60))
    beta = rng.dirichlet(np.ones(K) * 0.7, n_spots)
    Y = rng.poisson(beta @ X * 3).astype(float)
    coords = np.column_stack([np.arange(n_spots) % 15, np.arange(n_spots) // 15]).astype(float)
    return Y, X, coords


class TestEstimatorAPI:
    @pytest.mark.parametrize("to_sparse", [False, True])
    def test_default_output(self, fitted, to_sparse):
        Y, X, coords = fitted
        Yin = sparse.csr_matrix(Y) if to_sparse else Y
        m = FlashDeconv().fit(Yin, X, coords)
        uq = m.compute_uncertainty(fdr_q=0.1)
        assert uq["method"] == "sandwich"
        P = m.proportions_
        for k in ("ci_lower", "ci_upper", "se_prop", "z_score", "p_value", "detected"):
            assert uq[k].shape == P.shape
        assert np.all(uq["ci_lower"] <= P + 1e-12)
        assert np.all(uq["ci_upper"] >= P - 1e-12)
        inact = m.beta_ <= 0
        assert np.all(uq["ci_lower"][inact] == 0)
        assert uq["detected"].dtype == bool
        assert np.array_equal(uq["detection_confident"], uq["detected"])

    def test_sparse_matches_dense(self, fitted):
        Y, X, coords = fitted
        ud = FlashDeconv().fit(Y, X, coords).compute_uncertainty()
        us = FlashDeconv().fit(sparse.csr_matrix(Y), X, coords).compute_uncertainty()
        for k in ("ci_lower", "ci_upper", "z_score"):
            np.testing.assert_allclose(us[k], ud[k], rtol=1e-6, atol=1e-10)

    def test_deterministic(self, fitted):
        Y, X, coords = fitted
        m = FlashDeconv().fit(Y, X, coords)
        a = m.compute_uncertainty()
        b = m.compute_uncertainty()
        for k in ("ci_lower", "ci_upper", "se_prop", "detected"):
            np.testing.assert_array_equal(a[k], b[k])
        ba = m.bootstrap_uncertainty(n_bootstrap=3, seed=7)
        bb = m.bootstrap_uncertainty(n_bootstrap=3, seed=7)
        np.testing.assert_array_equal(ba["boot_ci_upper"], bb["boot_ci_upper"])

    def test_laplace_diag_deprecated(self, fitted):
        Y, X, coords = fitted
        m = FlashDeconv().fit(Y, X, coords)
        with pytest.warns(DeprecationWarning):
            uq = m.compute_uncertainty(method="laplace_diag")
        assert "detection_confident" in uq
        with pytest.raises(ValueError):
            m.compute_uncertainty(method="bogus")

    def test_countsketch_representation(self, fitted):
        Y, X, coords = fitted
        m = FlashDeconv(gene_weighting="countsketch", sketch_dim=64).fit(Y, X, coords)
        uq = m.compute_uncertainty()
        assert np.all(uq["ci_upper"] >= m.proportions_ - 1e-12)


class TestBootstrap:
    def test_warm_start_and_convergence(self, fitted, monkeypatch):
        Y, X, coords = fitted
        m = FlashDeconv().fit(Y, X, coords)
        import flashdeconv.core.solver as solver
        seen = []
        orig = solver.bcd_solve

        def spy(*args, **kwargs):
            seen.append(kwargs.get("beta_init"))
            return orig(*args, **kwargs)

        monkeypatch.setattr(solver, "bcd_solve", spy)
        with warnings.catch_warnings():
            warnings.simplefilter("error", RuntimeWarning)
            boot = m.bootstrap_uncertainty(n_bootstrap=4, seed=0)
        assert len(seen) == 4
        assert all(s is m.beta_ for s in seen)
        assert boot["warm_start"] is True
        assert boot["n_converged"] == 4
        assert np.all(boot["n_iterations"] < m.max_iter)
        assert boot["boot_mean"].shape == m.proportions_.shape
        assert np.all(boot["boot_ci_lower"] <= boot["boot_ci_upper"])

    def test_warm_start_converges_faster(self, fitted):
        Y, X, coords = fitted
        m = FlashDeconv().fit(Y, X, coords)
        rho = m.rho_sparsity
        _, cold = bcd_solve(m.Y_sketch_, m.X_sketch_, m.adjacency_, m.lambda_used_,
                            rho, max_iter=m.max_iter, tol=m.tol)
        _, warm = bcd_solve(m.Y_sketch_, m.X_sketch_, m.adjacency_, m.lambda_used_,
                            rho, max_iter=m.max_iter, tol=m.tol, beta_init=m.beta_)
        assert warm["converged"]
        assert warm["n_iterations"] <= 2 < cold["n_iterations"]

    def test_unconverged_refits_warn(self, fitted):
        Y, X, coords = fitted
        m = FlashDeconv().fit(Y, X, coords)
        with pytest.warns(RuntimeWarning):
            boot = m.bootstrap_uncertainty(n_bootstrap=2, max_iter_boot=1, tol=1e-14)
        assert boot["n_converged"] == 0
