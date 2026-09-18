# =============================================================================
# 03b_mc_estimators.py -- fast MC-specific estimator paths (optimized v2)
# =============================================================================
"""
Same math as 01_estimators (validated to 1e-15), optimized for the MC loop:
 - balanced-panel cumsum demeaning: full-sample X-demean per horizon and the
   two SPJ half-demeanings are computed from prefix sums (O(N) vs O(TN));
 - closed-form K<=2 OLS (solve of X'X), np.bincount cluster sums;
 - SPJ cut for balanced complete panels = (T-h-1)//2 (their floor(median));
   falls back to the exact per-entity searchsorted rule if unbalanced.
Cross-checked against 01_estimators to <1e-9 on 4 panel shapes x 6 horizons.
"""
import numpy as np


class _K1OLS:
    pass


def _ols1(yv, xv):
    d = xv @ xv
    return float((xv @ yv) / d) if d > 0 else np.nan


class FEFast:
    """FE with K=1; X-demeaned per horizon under the joint complete-case mask
    (== isfinite(dep) rows for MC panels where x has no NaNs)."""
    def __init__(self, X):
        x = np.asarray(X, float)
        self.x = x
        self.xcum = np.vstack([np.zeros((1, x.shape[1])), np.cumsum(x, axis=0)])

    def fit_h(self, dep):
        T, N = dep.shape
        mask = np.isfinite(dep)                      # x finite everywhere (MC)
        nrow = mask.sum(axis=0)                      # balanced: all == T-h
        if nrow.min() < 2:
            return np.nan, np.nan
        m = int(nrow[0]) if nrow.min() == nrow.max() else None
        if m is not None:
            xm = self.xcum[m] / m                    # (N,)
        else:  # unbalanced fallback (exact joint-mask means)
            sm = np.where(mask, self.x, 0.0).sum(axis=0)
            xm = sm / nrow
        xd = self.x - xm[None, :]
        ym = np.where(mask, dep, 0.0).sum(axis=0) / nrow
        yd = np.where(mask, dep - ym[None, :], np.nan)
        ok = mask.ravel() & np.isfinite(xd).ravel()
        yv = yd.ravel()[ok]; xv = xd.ravel()[ok]
        beta = _ols1(yv, xv)
        if not np.isfinite(beta):
            return np.nan, np.nan
        res = yv - beta * xv
        ents = np.repeat(np.arange(N), T)[ok]
        U = np.bincount(ents, weights=xv * res, minlength=N)
        Q = float(xv @ xv)
        V00 = (U @ U) / Q ** 2
        return beta, float(np.sqrt(max(V00, 0)))


class SPJFast:
    """SPJ with K=1. Balanced-panel cut c=(T-h-1)//2; half-X demeanings from
    prefix sums (constant across iterations); DD-sandwich SE (Q from the
    original full design, meat from DD values) -- identical to their R."""
    def __init__(self, X):
        x = np.asarray(X, float)
        self.x = x
        self.xcum = np.vstack([np.zeros((1, x.shape[1])), np.cumsum(x, axis=0)])

    def fit_h(self, dep):
        T, N = dep.shape
        nfin = np.isfinite(dep).sum(axis=0)
        if nfin.min() < 4:
            return np.nan, np.nan
        if nfin.min() == nfin.max():
            c = (int(nfin[0]) - 1) // 2              # floor(median) of 0..n-1
        else:
            cc = np.isfinite(dep) & np.isfinite(self.x)
            cpos = (cc.sum(axis=0) - 1) // 2
            cum = np.cumsum(cc, axis=0)
            c = int(np.median([np.searchsorted(cum[:, i], cpos[i] + 1)
                               for i in range(N)]))
        # ---- full (rows 0..n-1, n = T-h)
        n = int(nfin[0])
        xf = self.xcum[n] / n
        xdf = self.x - xf[None, :]
        ysum = np.nancumsum(np.where(np.isfinite(dep), dep, 0.0), axis=0)
        yf = ysum[n - 1] / n
        ydf = np.where(np.isfinite(dep), dep - yf[None, :], np.nan)
        okF = np.isfinite(ydf) & np.isfinite(xdf)
        yv = ydf.ravel()[okF.ravel()]; xv = xdf.ravel()[okF.ravel()]
        if len(yv) < 4:
            return np.nan, np.nan
        cf = _ols1(yv, xv)
        if not np.isfinite(cf):
            return np.nan, np.nan
        # ---- halves (A: rows 0..c ; B: rows c+1..n-1)
        na = c + 1
        xa = self.xcum[na] / na
        xdA = self.x - xa[None, :]
        xdA[c + 1:] = np.nan
        xb = (self.xcum[n] - self.xcum[na]) / (n - na)
        xdB = self.x - xb[None, :]
        xdB[:na] = np.nan
        ydA = np.where(np.arange(T)[:, None] <= c,
                       dep - (ysum[c] / na)[None, :], np.nan)
        ydB = np.where(np.arange(T)[:, None] > c,
                       dep - ((ysum[n - 1] - ysum[c]) / (n - na))[None, :], np.nan)
        okA = np.isfinite(ydA) & np.isfinite(xdA)
        okB = np.isfinite(ydB) & np.isfinite(xdB)
        if okA.sum() < 2 or okB.sum() < 2:
            return np.nan, np.nan
        ya = ydA.ravel()[okA.ravel()]; xav = xdA.ravel()[okA.ravel()]
        yb = ydB.ravel()[okB.ravel()]; xbv = xdB.ravel()[okB.ravel()]
        ca = _ols1(ya, xav)
        cb = _ols1(yb, xbv)
        if not (np.isfinite(ca) and np.isfinite(cb)):
            return np.nan, np.nan
        beta = 2 * cf - 0.5 * (ca + cb)
        # ---- DD sandwich: dd = 2*XdF - Xd_sub; Q from ORIGINAL full design
        Xd_sub = np.full_like(xdf, np.nan)
        Xd_sub[:na] = xdA[:na]
        Xd_sub[na:] = xdB[na:]
        dd = 2 * xdf - Xd_sub
        res_full = ydf - xdf * beta
        okm = np.isfinite(res_full) & np.isfinite(dd)
        okf = okm.ravel()
        ddv = dd.ravel()[okf]
        rv = res_full.ravel()[okf]
        ents2 = np.repeat(np.arange(N), T)[okf]
        U = np.bincount(ents2, weights=ddv * rv, minlength=N)
        Q = float(xv @ xv)                       # original-design X'X
        V00 = (U @ U) / Q ** 2
        return float(beta), float(np.sqrt(max(V00, 0)))


class HPJFast:
    """4.5*FE(N) - 4*FE(2N/3) + 0.5*FE(N/3); SE = full-sample FE sandwich."""
    def __init__(self, X):
        X = np.asarray(X, float)
        T, N = X.shape
        n2, n3 = max(int(round(2 * N / 3)), 2), max(int(round(N / 3)), 2)
        self.n2, self.n3 = n2, n3
        self.fe1 = FEFast(X)
        self.fe2 = FEFast(X[:, :n2])
        self.fe3 = FEFast(X[:, :n3])

    def fit_h(self, dep):
        b1, se1 = self.fe1.fit_h(dep)
        b2, _ = self.fe2.fit_h(dep[:, :self.n2])
        b3, _ = self.fe3.fit_h(dep[:, :self.n3])
        return float(4.5 * b1 - 4.0 * b2 + 0.5 * b3), se1


class CDOLSFast:
    """y_{t+h} - y_{t-1} on [s_t, s_{t-1}], pooled, entity-clustered (K<=2)."""
    def __init__(self, s, include_lag=True):
        self.s = np.asarray(s, float)
        self.lag = include_lag

    def fit_h(self, y, h):
        T, N = y.shape
        tt = np.arange(1, T - h)
        if len(tt) < 4:
            return np.nan, np.nan
        L = (y[tt + h, :] - y[tt - 1, :]).T          # N x Tt
        C = np.stack([self.s[tt, :]] + ([self.s[tt - 1, :]] if self.lag else []),
                     axis=2)                          # (Tt, N, K)
        K = C.shape[2]
        Yv = L.reshape(-1)
        Xv = C.transpose(1, 0, 2).reshape(N * len(tt), K)
        ok = np.isfinite(Yv) & np.isfinite(Xv).all(axis=1)
        if ok.sum() < 3 * K:
            return np.nan, np.nan
        yv, xv = Yv[ok], Xv[ok]
        ents = np.repeat(np.arange(N), len(tt))[ok]
        Q = xv.T @ xv
        b = np.linalg.solve(Q, xv.T @ yv)
        res = yv - xv @ b
        Uk = np.zeros((N, K))
        for j in range(K):
            Uk[:, j] = np.bincount(ents, weights=xv[:, j] * res, minlength=N)
        try:
            Qinv = np.linalg.inv(Q)
        except np.linalg.LinAlgError:
            Qinv = np.linalg.pinv(Q)
        V = Qinv @ (Uk.T @ Uk) @ Qinv
        return float(b[0]), float(np.sqrt(max(V[0, 0], 0)))


# =================================================================== S9 utils
def fe_se_variants(dep, x, beta=None):
    """FE (K=1) beta + three SE treatments: entity / time / two-way clusters.
    Returns (beta, se_ent, se_time, se_tw)."""
    T, N = dep.shape
    mask = np.isfinite(dep) & np.isfinite(x)
    nrow = mask.sum(axis=0)
    if nrow.min() < 2:
        return np.nan, np.nan, np.nan, np.nan
    xm = np.where(mask, x, 0.0).sum(axis=0) / nrow
    xd = np.where(mask, x - xm[None, :], np.nan)
    ym = np.where(mask, dep, 0.0).sum(axis=0) / nrow
    yd = np.where(mask, dep - ym[None, :], np.nan)
    ok = (mask & np.isfinite(xd) & np.isfinite(yd)).ravel()
    yv = yd.ravel()[ok]; xv = xd.ravel()[ok]
    if beta is None or not np.isfinite(beta):
        beta = _ols1(yv, xv)
    res = yv - beta * xv
    ents = np.repeat(np.arange(N), T)[ok]
    times = np.tile(np.arange(T), N)[ok]
    Ue = np.bincount(ents, weights=xv * res, minlength=N)
    Ut = np.bincount(times, weights=xv * res, minlength=T)
    Q = float(xv @ xv)
    Q2 = Q ** 2
    se_e = np.sqrt((Ue @ Ue) / Q2)
    se_t = np.sqrt((Ut @ Ut) / Q2)
    # two-way (Cameron-Miller): intersection clusters are singletons ->
    # U_et = sum_g (x_g r_g)^2
    Uet = float(np.sum(xv ** 2 * res ** 2))
    se_tw = np.sqrt(max((Ue @ Ue) + (Ut @ Ut) - Uet, 0.0) / Q2)
    return float(beta), float(se_e), float(se_t), float(se_tw)


def spj_se_variants(dep, x):
    """SPJ (K=1) beta + entity/time/two-way SEs on the DD design
    (Q from the original full design; exact port semantics)."""
    T, N = dep.shape
    f = SPJFast(x)
    b, se_e = f.fit_h(dep)
    if not np.isfinite(b):
        return np.nan, np.nan, np.nan, np.nan
    # rebuild the DD design pieces exactly as fit_h does
    nfin = np.isfinite(dep).sum(axis=0)
    n = int(nfin[0]) if nfin.min() == nfin.max() else T
    if nfin.min() == nfin.max():
        c = (int(nfin[0]) - 1) // 2
    else:
        cc = np.isfinite(dep) & np.isfinite(x)
        cpos = (cc.sum(axis=0) - 1) // 2
        cum = np.cumsum(cc, axis=0)
        c = int(np.median([np.searchsorted(cum[:, i], cpos[i] + 1)
                           for i in range(N)]))
    xc = f.xcum
    na = c + 1
    xf = xc[n] / n
    xdf = x - xf[None, :]
    ysum = np.nancumsum(np.where(np.isfinite(dep), dep, 0.0), axis=0)
    yf = ysum[n - 1] / n
    ydf = np.where(np.isfinite(dep), dep - yf[None, :], np.nan)
    xa = xc[na] / na
    xdA = x - xa[None, :]; xdA[na:] = np.nan
    xb = (xc[n] - xc[na]) / (n - na)
    xdB = x - xb[None, :]; xdB[:na] = np.nan
    ydA = np.where(np.arange(T)[:, None] <= c, dep - (ysum[c] / na)[None, :], np.nan)
    ydB = np.where(np.arange(T)[:, None] > c,
                   dep - ((ysum[n - 1] - ysum[c]) / (n - na))[None, :], np.nan)
    Xd_sub = np.full_like(xdf, np.nan)
    Xd_sub[:na] = xdA[:na]
    Xd_sub[na:] = xdB[na:]
    dd = 2 * xdf - Xd_sub
    res_full = ydf - xdf * b
    okm = (np.isfinite(res_full) & np.isfinite(dd)).ravel()
    ddv = dd.ravel()[okm]; rv = res_full.ravel()[okm]
    ents2 = np.repeat(np.arange(N), T)[okm]
    times2 = np.tile(np.arange(T), N)[okm]
    Ue = np.bincount(ents2, weights=ddv * rv, minlength=N)
    Ut = np.bincount(times2, weights=ddv * rv, minlength=T)
    xv = xdf.ravel()[okm]
    Q = float(xv @ xv)
    Q2 = Q ** 2
    se_t = np.sqrt((Ut @ Ut) / Q2)
    Uet = float(np.sum(ddv ** 2 * rv ** 2))
    se_tw = np.sqrt(max((Ue @ Ue) + (Ut @ Ut) - Uet, 0.0) / Q2)
    return float(b), float(se_e), float(se_t), float(se_tw)
