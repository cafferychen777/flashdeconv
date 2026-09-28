"""
Incomplete-reference diagnostic (opt-in, post hoc on a fitted FlashDeconv model).

A cell type that is absent from the reference leaves bins whose counts cannot be reproduced by
any mixture of the reference profiles: its marker genes appear although no reference type
expresses them. This module scores that per bin and names the unexplained genes.

Per-bin score (multinomial-mixture surprisal)
---------------------------------------------
On the genes selected by the fit, bin ``i`` with ``n_i`` UMIs is compared with the best
count-space mixture of reference profiles,
``q_i = (1 - eta) * sum_k p_ik xbar_k + eta / G`` (``xbar_k`` = reference profile normalised to
sum 1, ``eta`` = small uniform ambient floor, ``p_i`` = per-bin maximum-likelihood mixture weights
obtained by EM started from the fitted proportions). With
``l_i = (1/n_i) sum_g y_ig log q_ig`` the per-UMI log-likelihood, under
``y_i ~ Multinomial(n_i, q_i)`` one has ``E[l_i] = -H(q_i)`` and
``Var[l_i] = V(q_i) / n_i`` with ``V`` the variance of ``log q`` under ``q``. The raw score

    z_i = -(l_i + H(q_i)) * sqrt(n_i / V(q_i))

is therefore depth- and composition-aware (positive = worse than expected). Counts simulated
from the fitted mixture give ``z ~ N(0, 1)`` at every depth, but real counts deviate from the
multinomial by a per-UMI amount that ``z`` multiplies by ``sqrt(n_i)``, so ``z`` is re-centred and
re-scaled with an empirical null that ignores the upper tail created by missing cell types:
by default the smaller-scale of (i) the median and left-half spread and (ii) central matching
(mode and curvature of the density near its peak), see ``null``. The pooled score sums the
likelihood deviations of each bin and its spatial neighbours (the fitted kNN graph), which is
exact under independence and gains power when the missing type forms a region.

Gene-level residuals
--------------------
For a set of bins ``R`` (e.g. flagged bins) and its complement ``U``, observed counts are compared
with the counts expected from the fitted proportions on all genes shared by the data and the
reference; the score ``log2(O/E)_R - log2(O/E)_U`` cancels gene-level platform effects common
to the section. Top genes can be mapped to cell types of any broader atlas with
:func:`suggest_missing_types`.

Cost: one pass over the non-zeros of Y for the EM and likelihood, plus a chunked
``(n_bins x K) @ (K x n_selected_genes)`` product for ``H`` and ``V``.
"""

from typing import Dict, Optional, Sequence

import numpy as np
from numba import njit, prange
from scipy import sparse

_Z95 = 1.6448536269514722


@njit(parallel=True, cache=True)
def _mixture_loglik(indptr, indices, data, Xbar, P0, eta, n_em, em_tol):
    """Per-bin EM for the mixture weights and the per-UMI log-likelihood (sparse rows)."""
    N = len(indptr) - 1
    K, G = Xbar.shape
    ell = np.zeros(N)
    nn = np.zeros(N)
    P = P0.copy()
    fl = eta / G
    for i in prange(N):
        a, b = indptr[i], indptr[i + 1]
        n = 0.0
        for t in range(a, b):
            n += data[t]
        nn[i] = n
        p = P[i]
        s = 0.0
        for k in range(K):
            if p[k] < 0.0:
                p[k] = 0.0
            s += p[k]
        if s <= 0.0:
            for k in range(K):
                p[k] = 1.0 / K
        else:
            for k in range(K):
                p[k] /= s
        if n <= 0.0:
            continue
        acc = np.zeros(K)
        for _ in range(n_em):
            for k in range(K):
                acc[k] = 0.0
            for t in range(a, b):
                g = indices[t]
                q = 0.0
                for k in range(K):
                    q += p[k] * Xbar[k, g]
                q = (1.0 - eta) * q + fl
                f = data[t] * (1.0 - eta) / q
                for k in range(K):
                    acc[k] += f * p[k] * Xbar[k, g]
            tot = 0.0
            for k in range(K):
                tot += acc[k]
            if tot <= 0.0:
                break
            dmax = 0.0
            for k in range(K):
                nv = acc[k] / tot
                d = abs(nv - p[k])
                if d > dmax:
                    dmax = d
                p[k] = nv
            if dmax < em_tol:
                break
        ll = 0.0
        for t in range(a, b):
            g = indices[t]
            q = 0.0
            for k in range(K):
                q += p[k] * Xbar[k, g]
            q = (1.0 - eta) * q + fl
            ll += data[t] * np.log(q)
        ell[i] = ll / n
    return ell, nn, P


@njit(parallel=True, cache=True)
def _entropy_rows(Q):
    N, G = Q.shape
    H = np.empty(N)
    V = np.empty(N)
    for i in prange(N):
        h = 0.0
        e2 = 0.0
        for g in range(G):
            q = Q[i, g]
            lq = np.log(q)
            h -= q * lq
            e2 += q * lq * lq
        H[i] = h
        V[i] = max(e2 - h * h, 1e-12)
    return H, V


def _entropy_var(P, Xbar, eta, chunk=8192):
    N = P.shape[0]
    G = Xbar.shape[1]
    H = np.empty(N)
    V = np.empty(N)
    for s in range(0, N, chunk):
        e = min(s + chunk, N)
        Q = (1.0 - eta) * (P[s:e] @ Xbar) + eta / G
        H[s:e], V[s:e] = _entropy_rows(Q)
    return H, V


_NULL_METHODS = ("auto", "left_half", "central")


def _left_half_null(z: np.ndarray, lo_q: float = 0.16):
    """Median and left-half robust spread (median - q_lo)."""
    mu = float(np.median(z))
    return mu, max(mu - float(np.quantile(z, lo_q)), 1e-9)


def _central_null(z: np.ndarray, drop: float = 1.0):
    """Central matching (Efron 2004): mode and curvature of the log density near its peak.

    Histogram of the central 99.8% of ``z``, lightly Gaussian-smoothed; quadratic fit of the log
    counts over the contiguous window around the mode where the smoothed density is at least
    ``exp(-drop)`` times its peak (``drop = 1``: about +-1.4 null SD). Robust to both tails.
    """
    from scipy.ndimage import gaussian_filter1d

    lo, hi = np.quantile(z, [0.001, 0.999])
    if not hi > lo:
        return _left_half_null(z)
    nb = int(np.clip(len(z) // 100, 30, 400))
    h, e = np.histogram(z, bins=nb, range=(lo, hi))
    c = 0.5 * (e[1:] + e[:-1])
    w = e[1] - e[0]
    q25, q75 = np.quantile(z, [0.25, 0.75])
    bw = max(0.1 * (q75 - q25) / 1.349, w)
    hs = gaussian_filter1d(h.astype(np.float64), bw / w, mode="constant")
    i0 = int(np.argmax(hs))
    thr = hs[i0] * np.exp(-drop)
    a = b = i0
    while a > 0 and hs[a - 1] >= thr:
        a -= 1
    while b < nb - 1 and hs[b + 1] >= thr:
        b += 1
    if b - a < 4:
        a, b = max(0, i0 - 3), min(nb - 1, i0 + 3)
    y = hs[a:b + 1]
    d2, d1, _ = np.polyfit(c[a:b + 1], np.log(np.maximum(y, 1e-12)), 2,
                           w=np.sqrt(np.maximum(y, 1e-12)))
    if d2 >= 0:
        return _left_half_null(z)
    s2 = -1.0 / (2.0 * d2) - bw ** 2  # remove the smoothing kernel's variance
    return float(-d1 / (2.0 * d2)), float(np.sqrt(max(s2, 1e-18)))


def _empirical_null(z: np.ndarray, method: str = "auto", keep: Optional[np.ndarray] = None):
    """Calibrate ``z`` against its own bulk; returns (calibrated z, info dict).

    ``left_half``: median and left-half spread; inflated when many bins fit BETTER than the
    multinomial expectation (a heavy lower shoulder, e.g. deep, homogeneous tissue). ``central``:
    central matching; inflated when the peak is flattened by a broad mixture. ``auto``: the
    estimate with the smaller scale, i.e. the one less inflated by either departure.
    """
    zk = z if keep is None or not np.any(keep) else z[keep]
    if method == "left_half":
        mu, sd, used = *_left_half_null(zk), "left_half"
    elif method == "central":
        mu, sd, used = *_central_null(zk), "central"
    else:
        m1, s1 = _left_half_null(zk)
        m2, s2 = _central_null(zk)
        mu, sd, used = (m2, s2, "central") if s2 < s1 else (m1, s1, "left_half")
    return (z - mu) / sd, {"method": used, "center": mu, "scale": sd}


def reference_fit_scores(
    model,
    eta: float = 0.01,
    max_em_iter: int = 10,
    em_tol: float = 1e-4,
    pool: bool = True,
    alpha: float = 0.05,
    null: str = "auto",
) -> Dict[str, np.ndarray]:
    """
    Per-bin score of how poorly the reference explains each bin (opt-in diagnostic).

    Parameters
    ----------
    model : FlashDeconv
        A fitted model (any ``preprocess`` / ``gene_weighting``; the score uses the raw counts
        of the selected genes and the reference profiles, not the solver representation).
    eta : float, default=0.01
        Weight of a uniform ambient component in the expected profile; bounds the surprisal of
        a single UMI at ``log(G / eta)``.
    max_em_iter : int, default=10
        Maximum EM iterations for the per-bin mixture weights (started from
        ``model.proportions_``; 0 uses the fitted proportions as they are).
    em_tol : float, default=1e-4
        EM stops for a bin when no mixture weight changes by more than ``em_tol``.
    pool : bool, default=True
        Also compute the neighbourhood-pooled score over the fitted spatial graph.
    alpha : float, default=0.05
        Nominal upper-tail level of the flags (one-sided, on the calibrated scale).
    null : {'auto', 'left_half', 'central'}, default='auto'
        Empirical null used to calibrate the raw scores (separately for the unpooled and the
        pooled score). ``'left_half'``: median and left-half spread (median - 16th percentile;
        the previous default). ``'central'``: central matching (mode and curvature of the log
        density near its peak). ``'auto'``: the one with the smaller scale. The null only sets
        the flag threshold: all three are monotone in ``z_raw`` and rank bins identically.

    Returns
    -------
    scores : dict
        ``'z_raw'`` : analytic multinomial score (not calibrated for over-dispersion);
        ``'score'`` : calibrated score (empirical null centre 0, scale 1);
        ``'flag'`` : ``score > Phi^{-1}(1 - alpha)``;
        ``'score_pooled'``, ``'flag_pooled'`` : same for the neighbourhood-pooled score
        (when ``pool`` and a spatial graph is available);
        ``'n_umi'`` : UMIs on the selected genes; ``'mixture_weights'`` : per-bin count-space
        mixture weights used for the expectation (n_bins x n_cell_types);
        ``'null'`` : dict with the estimator used, centre and scale for ``'score'`` (and
        ``'score_pooled'``).

    Notes
    -----
    Flags are relative to the bulk of the section: the empirical null assumes that most bins are
    explained by the reference. The score ranks bins; the flag level is approximate.
    Real sections can contain many bins that the reference mixture explains *better* than the
    multinomial expectation (e.g. deep, homogeneous tumour epithelium); because the raw score
    scales with sqrt(depth), they form a heavy lower shoulder that inflates the left-half spread
    and suppresses almost all flags. ``'auto'`` then falls back to central matching.
    """
    from scipy.stats import norm

    if not getattr(model, "_fitted", False):
        raise RuntimeError("Model has not been fitted. Call fit() first.")
    if not 0.0 < eta < 1.0:
        raise ValueError(f"eta must be in (0, 1), got {eta}")
    if max_em_iter < 0:
        raise ValueError(f"max_em_iter must be non-negative, got {max_em_iter}")
    if not 0.0 < alpha < 1.0:
        raise ValueError(f"alpha must be in (0, 1), got {alpha}")
    if null not in _NULL_METHODS:
        raise ValueError(f"null must be one of {_NULL_METHODS}, got {null!r}")

    gidx = model.gene_idx_
    Y = model._Y_raw[:, gidx]
    Y = sparse.csr_matrix(Y, dtype=np.float64)
    Y.sort_indices()
    X = np.asarray(model._X_raw[:, gidx], dtype=np.float64)
    Xbar = np.ascontiguousarray(X / np.maximum(X.sum(axis=1, keepdims=True), 1e-300))
    P0 = np.ascontiguousarray(model.proportions_, dtype=np.float64)

    ell, n, P = _mixture_loglik(
        Y.indptr.astype(np.int64), Y.indices.astype(np.int64), Y.data,
        Xbar, P0, float(eta), int(max_em_iter), float(em_tol),
    )
    H, V = _entropy_var(P, Xbar, float(eta))
    dev = n * (ell + H)
    z = -dev / np.sqrt(np.maximum(n * V, 1e-12))
    z[n <= 0] = 0.0
    thr = norm.ppf(1.0 - alpha)

    out = {"z_raw": z, "n_umi": n, "mixture_weights": P}
    out["score"], info = _empirical_null(z, null, keep=n > 0)
    out["flag"] = out["score"] > thr
    out["null"] = {"score": info}

    A = getattr(model, "adjacency_", None)
    if pool and A is not None:
        A = sparse.csr_matrix(A, dtype=np.float64)
        A = (A + sparse.identity(A.shape[0], format="csr")).tocsr()
        A.data[:] = 1.0
        zp = -(A @ dev) / np.sqrt(np.maximum(A @ (n * V), 1e-12))
        out["z_raw_pooled"] = zp
        out["score_pooled"], out["null"]["score_pooled"] = _empirical_null(zp, null)
        out["flag_pooled"] = out["score_pooled"] > thr
    return out


def unexplained_genes(
    model,
    bins,
    gene_names: Optional[Sequence[str]] = None,
    min_count: int = 20,
    pseudocount: float = 0.5,
) -> Dict[str, np.ndarray]:
    """
    Genes over-represented in ``bins`` relative to what the fitted reference mixture predicts.

    Parameters
    ----------
    model : FlashDeconv
        A fitted model.
    bins : array-like of bool (n_bins,) or integer indices
        Bins to characterise (e.g. ``reference_fit_scores(model)['flag_pooled']``).
    gene_names : sequence of str, optional
        Names of the genes of ``Y``/``X`` as passed to ``fit`` (same order).
    min_count : int, default=20
        Genes with fewer observed counts inside ``bins`` get ``score = nan``.
    pseudocount : float, default=0.5

    Returns
    -------
    genes : dict of arrays sorted by decreasing ``score``
        ``'gene'`` (index or name), ``'score'`` (log2 O/E inside minus outside), ``'observed'``,
        ``'expected'`` (inside ``bins``), ``'log2_oe_in'``, ``'log2_oe_out'``.
    """
    if not getattr(model, "_fitted", False):
        raise RuntimeError("Model has not been fitted. Call fit() first.")
    Y = model._Y_raw
    Y = sparse.csr_matrix(Y, dtype=np.float64)
    X = np.asarray(model._X_raw, dtype=np.float64)
    n_bins, n_genes = Y.shape
    mask = np.zeros(n_bins, dtype=bool)
    b = np.asarray(bins)
    if b.dtype == bool:
        if b.shape[0] != n_bins:
            raise ValueError("boolean `bins` must have length n_bins")
        mask = b
    else:
        mask[b.astype(np.int64)] = True
    if mask.sum() == 0 or mask.sum() == n_bins:
        raise ValueError("`bins` must select a non-empty proper subset of bins")
    if gene_names is None:
        gene_names = np.arange(n_genes)
    gene_names = np.asarray(gene_names)
    if gene_names.shape[0] != n_genes:
        raise ValueError("gene_names length does not match the number of genes")

    Xbar = X / np.maximum(X.sum(axis=1, keepdims=True), 1e-300)
    N_i = np.asarray(Y.sum(axis=1)).ravel()
    P = model.proportions_
    O_in = np.asarray(Y[mask].sum(axis=0)).ravel()
    O_out = np.asarray(Y[~mask].sum(axis=0)).ravel()
    E_in = (N_i[mask] @ P[mask]) @ Xbar
    E_out = (N_i[~mask] @ P[~mask]) @ Xbar
    E_in *= O_in.sum() / max(E_in.sum(), 1e-300)
    E_out *= O_out.sum() / max(E_out.sum(), 1e-300)
    lr_in = np.log2((O_in + pseudocount) / (E_in + pseudocount))
    lr_out = np.log2((O_out + pseudocount) / (E_out + pseudocount))
    score = lr_in - lr_out
    score[O_in < min_count] = np.nan
    order = np.argsort(np.where(np.isnan(score), np.inf, -score), kind="stable")
    return {
        "gene": gene_names[order], "score": score[order], "observed": O_in[order],
        "expected": E_in[order], "log2_oe_in": lr_in[order], "log2_oe_out": lr_out[order],
    }


def suggest_missing_types(
    genes: Dict[str, np.ndarray],
    atlas_profiles: np.ndarray,
    atlas_genes: Sequence[str],
    atlas_names: Sequence[str],
    n_markers: int = 50,
) -> Dict[str, np.ndarray]:
    """
    Rank the cell types of a broader atlas by how strongly their markers are unexplained.

    Parameters
    ----------
    genes : dict
        Output of :func:`unexplained_genes` (``'gene'`` must hold gene names).
    atlas_profiles : ndarray of shape (n_atlas_types, n_atlas_genes)
        Mean expression (counts or normalised) per atlas cell type.
    atlas_genes, atlas_names : sequences
        Gene names (columns) and type names (rows) of ``atlas_profiles``.
    n_markers : int, default=50
        Markers per atlas type: top genes by log2 fold change of CP10k versus the mean of the
        other atlas types.

    Returns
    -------
    ranking : dict of arrays sorted by decreasing ``mean_score``
        ``'type'``, ``'mean_score'`` (mean residual score over the type's scored markers),
        ``'n_markers_scored'``.
    """
    A = np.asarray(atlas_profiles, dtype=np.float64)
    atlas_genes = np.asarray(atlas_genes)
    if A.shape != (len(atlas_names), len(atlas_genes)):
        raise ValueError("atlas_profiles must have shape (len(atlas_names), len(atlas_genes))")
    lookup = {str(g): s for g, s in zip(genes["gene"], genes["score"]) if np.isfinite(s)}
    cp = A / np.maximum(A.sum(axis=1, keepdims=True), 1e-300) * 1e4
    K = cp.shape[0]
    tot = cp.sum(axis=0)
    means, counts = np.full(K, np.nan), np.zeros(K, dtype=np.int64)
    for k in range(K):
        other = (tot - cp[k]) / max(K - 1, 1)
        lfc = np.log2((cp[k] + 0.1) / (other + 0.1))
        mk = atlas_genes[np.argsort(-lfc, kind="stable")[:n_markers]]
        v = [lookup[str(g)] for g in mk if str(g) in lookup]
        counts[k] = len(v)
        if v:
            means[k] = float(np.mean(v))
    order = np.argsort(np.where(np.isnan(means), np.inf, -means), kind="stable")
    return {"type": np.asarray(atlas_names)[order], "mean_score": means[order],
            "n_markers_scored": counts[order]}
