# =============================================================================
# 01_estimators.py -- panel LP estimators: FE, SPJ, HPJ, DB, CD-OLS
# =============================================================================
"""
Ported to numpy from the replication package's LP_panel_all.r
(Nickell Bias in Panel Local Projection, Mei-Sheng-Sheng 2026 JIE).
Matrix layout: wide T x N with NaN padding (their missing-data padding).
Estimators:
  FE   : within (entity) demeaning [+ time demeaning if te], cluster sandwich
  SPJ  : 2*full - 0.5*(halfA + halfB), per-entity median complete-case split,
         DD (2*XX - XX_sub) sandwich  [exact port]
  HPJ  : second-order 1/N correction: 4.5*b(N) - 4*b(2N/3) + 0.5*b(N/3)
         (eliminates B/N and C/N^2 bias terms; full-sample FE sandwich SE)
  DB   : data-based analytical correction, tau=0 case of their simul_LP.R:
         IRF_DB_h = IRF_FE_h + b0 * f(rho,h) / ((T0-h) * sx^2),
         f(rho,h) = su^2 * ((1-rho^h) - h/(T0-h)) / (1-rho)^2
  CDOLS: cumulative-difference pooled OLS (Ugarte-Ruiz 2026 BBVA WP 26/09,
         eq. 20): y_{i,t+h} - y_{i,t-1} on s_{i,t}, s_{i,t-1}; pooled OLS,
         entity-clustered SE. No within transformation.
"""
import numpy as np


# ----------------------------------------------------------------- utilities
def _demean_i(y, X):
    """Entity demeaning with joint complete-case mask (their want.mean logic).
    y: T x N, X: T x N x K. Returns (yd, Xd) with NaNs preserved."""
    mask = np.isfinite(y) & np.isfinite(X).all(axis=2)
    ym = np.where(mask, y, np.nan)
    yd = y - np.nanmean(ym, axis=0, keepdims=True)
    yd = np.where(mask, yd, np.nan)
    Xd = X - np.nanmean(np.where(mask[:, :, None], X, np.nan), axis=0, keepdims=True)
    Xd = np.where(mask[:, :, None], Xd, np.nan)
    return yd, Xd


def _demean_t(y, X):
    """Time demeaning (second step when te=True); mask across entities per row."""
    mask = np.isfinite(y) & np.isfinite(X).all(axis=2)
    ym = np.where(mask, y, np.nan)
    yd = y - np.nanmean(ym, axis=1, keepdims=True)
    yd = np.where(mask, yd, np.nan)
    Xd = X - np.nanmean(np.where(mask[:, :, None], X, np.nan), axis=1, keepdims=True)
    Xd = np.where(mask[:, :, None], Xd, np.nan)
    return yd, Xd


def _ols_fit(Y, X):
    """OLS on complete rows; returns (coef, res_full, mask)."""
    mask = np.isfinite(Y).ravel() & np.isfinite(X).all(axis=1)
    if mask.sum() == 0 or mask.sum() <= X.shape[1]:
        return np.full(X.shape[1], np.nan), np.full(len(Y), np.nan), mask
    coef, *_ = np.linalg.lstsq(X[mask], Y[mask], rcond=None)
    coef = coef.ravel()
    res = Y.ravel() - (X @ coef)
    return coef, res, mask


def _cluster_sandwich(Xv, rv, ents, N, times, T_eff, mode, robust):
    """Cluster sandwich on complete rows. ents/times are group ids per row."""
    k = Xv.shape[1]
    Q = Xv.T @ Xv

    def csum(groups, G):
        S = np.zeros((k, k))
        for g in range(G):
            m = groups == g
            if m.sum():
                u = Xv[m].T @ rv[m]
                S += np.outer(u, u)
        return S

    S = csum(ents, N)
    if mode == "twoway":
        S = S + csum(times, T_eff) - Xv.T @ (rv[:, None] * Xv)
        if robust:
            smp = len(rv)
            S = S * (min(N, T_eff) / (min(N, T_eff) - 1) *
                     (smp - 1) / (smp - N - k))
    elif robust:
        S = S * (N / (N - 1) * (len(rv) - 1) / (len(rv) - k))
    try:
        Qinv = np.linalg.inv(Q)
    except np.linalg.LinAlgError:
        Qinv = np.linalg.pinv(Q)
    return Qinv @ S @ Qinv


def _design(dep, Xmats, Wmats, te):
    """Build demeaned design. dep: T x N; Xmats: list of T x N;
    Wmats: list of T x N (controls). Returns (yd T x N, Xd T x N x K, K)."""
    X3 = np.stack([np.asarray(x, dtype=float) for x in Xmats], axis=2)
    if Wmats:
        X3 = np.concatenate([X3, np.stack([np.asarray(w, dtype=float)
                                           for w in Wmats], axis=2)], axis=2)
    yd, Xd = _demean_i(dep, X3)
    if te:
        yd, Xd = _demean_t(yd, Xd)
    return yd, Xd, X3.shape[2]


def _stack(yd, Xd):
    """Entity-major stacked vectors: row j = (entity i, time t), i-major."""
    T, N, K = Xd.shape
    ok = np.isfinite(yd).ravel() & np.isfinite(Xd).all(axis=2).ravel()
    Xv = Xd.reshape(T * N, K)[ok]
    yv = yd.reshape(T * N)[ok]
    ents = np.repeat(np.arange(N), T)[ok]
    times = np.tile(np.arange(T), N)[ok]
    return yv[:, None], Xv, ents, times, ok


def fe_core(dep, Xmats, Wmats=None, te=False, cluster="entity", robust=False,
            eigen=False, return_V=False):
    """FE panel LP for one horizon. dep: T x N (y_{t+h}); Xmats[0]: s_t."""
    yd, Xd, K = _design(dep, Xmats, Wmats, te)
    yv, Xv, ents, times, ok = _stack(yd, Xd)
    if len(yv) == 0:
        kx = len(Xmats)
        out = (np.full(kx, np.nan), np.full(kx, np.nan))
        return out + (None,) if return_V else out
    coef, res, mask = _ols_fit(yv, Xv)
    T_eff = dep.shape[0]
    V = _cluster_sandwich(Xv, res[mask], ents[mask], dep.shape[1], times[mask],
                          T_eff, cluster, robust)
    if eigen and V is not None:
        vals, vecs = np.linalg.eigh(V)
        V = vecs @ np.diag(np.clip(vals, 0, None)) @ vecs.T
    se = np.sqrt(np.clip(np.diag(V), 0, None))
    out = (coef[:len(Xmats)].copy(), se[:len(Xmats)].copy())
    return out + (V,) if return_V else out


def spj_core(dep, Xmats, Wmats=None, te=False, cluster="entity", robust=False,
             eigen=False):
    """Split-panel jackknife, exact port of LP_panel_all.r SPJ branch."""
    T, N = dep.shape[0], dep.shape[1]
    X3 = np.stack([np.asarray(x, dtype=float) for x in Xmats], axis=2)
    if Wmats:
        X3 = np.concatenate([X3, np.stack([np.asarray(w, dtype=float)
                                           for w in Wmats], axis=2)], axis=2)
    # ---- full-sample estimate
    ydF, XdF, K = _design(dep, Xmats, Wmats, te)
    yv, Xv, ents, times, ok = _stack(ydF, XdF)
    if len(yv) == 0:
        kx = len(Xmats)
        return np.full(kx, np.nan), np.full(kx, np.nan)
    coef_f, res_f, mask_f = _ols_fit(yv, Xv)
    if not np.isfinite(coef_f).all():
        kx = len(Xmats)
        return np.full(kx, np.nan), np.full(kx, np.nan)
    # ---- per-entity median cut of the complete-case indicator (their `cut`)
    cc = np.isfinite(dep) & np.isfinite(X3).all(axis=2)
    cut = np.zeros(N, dtype=int)
    for i in range(N):
        idx = np.where(cc[:, i])[0]
        cut[i] = int(np.median(idx)) if len(idx) else 0
    # ---- halves: mask rows > / <= cut per entity, re-demean within each half
    rows = np.arange(T)
    depA, depB = dep.copy(), dep.copy()
    XA, XB = X3.copy(), X3.copy()
    for i in range(N):
        c = cut[i]
        depA[rows > c, i] = np.nan
        depB[rows <= c, i] = np.nan
        XA[rows > c, i, :] = np.nan
        XB[rows <= c, i, :] = np.nan

    def half_beta(dep_h, X_h):
        Ktot = X_h.shape[2]
        yd, Xd, _ = _design(dep_h, [X_h[:, :, j] for j in range(len(Xmats))],
                            [X_h[:, :, len(Xmats) + j] for j in range(Ktot - len(Xmats))]
                            if Ktot > len(Xmats) else None, te)
        yv2, Xv2, _, _, ok2 = _stack(yd, Xd)
        if len(yv2) == 0:
            return np.full(Ktot, np.nan)
        c2, _, m2 = _ols_fit(yv2, Xv2)
        return c2

    beta_a = half_beta(depA, XA)
    beta_b = half_beta(depB, XB)
    beta = 2 * coef_f - 0.5 * (beta_a + beta_b)
    # ---- DD sandwich: dd = 2*XdF - Xd_sub  (Xd_sub assembled from half-demeaned)
    ydA, XdA, _ = _design(depA, [XA[:, :, j] for j in range(len(Xmats))],
                          [XA[:, :, len(Xmats) + j] for j in range(X3.shape[2] - len(Xmats))]
                          if X3.shape[2] > len(Xmats) else None, te)
    ydB, XdB, _ = _design(depB, [XB[:, :, j] for j in range(len(Xmats))],
                          [XB[:, :, len(Xmats) + j] for j in range(X3.shape[2] - len(Xmats))]
                          if X3.shape[2] > len(Xmats) else None, te)
    Xd_sub = np.full_like(XdF, np.nan)
    for i in range(N):
        c = cut[i]
        Xd_sub[:c + 1, i, :] = XdA[:c + 1, i, :]
        Xd_sub[c + 1:, i, :] = XdB[c + 1:, i, :]
    dd = 2 * XdF - Xd_sub
    res_full = ydF - np.einsum("tnk,k->tn", XdF, beta)
    okm = np.isfinite(res_full) & np.isfinite(dd).all(axis=2)
    okf = okm.ravel()
    # Q from the ORIGINAL full design (their Q.hat = indep_var' indep_var)
    Xo = XdF.reshape(T * N, K)[okf]
    Qo = Xo.T @ Xo
    Xv2 = dd.reshape(T * N, K)[okf]
    rv = res_full.ravel()[okf]
    ents2 = np.repeat(np.arange(N), T)[okf]
    times2 = np.tile(np.arange(T), N)[okf]
    # meat from DD values
    def csum(groups, G):
        S = np.zeros((K, K))
        for g in range(G):
            m = groups == g
            if m.sum():
                u = Xv2[m].T @ rv[m]
                S += np.outer(u, u)
        return S
    S = csum(ents2, N)
    if cluster == "twoway":
        S = S + csum(times2, T) - Xv2.T @ (rv[:, None] * Xv2)
    if robust:
        smp = int(okf.sum())
        S = S * (min(N, T) / (min(N, T) - 1) * (smp - 1) / (smp - N - K))
    try:
        Qinv = np.linalg.inv(Qo)
    except np.linalg.LinAlgError:
        Qinv = np.linalg.pinv(Qo)
    V = Qinv @ S @ Qinv
    if eigen and V is not None:
        vals, vecs = np.linalg.eigh(V)
        V = vecs @ np.diag(np.clip(vals, 0, None)) @ vecs.T
    se = np.sqrt(np.clip(np.diag(V), 0, None))
    return beta.copy(), se[:len(Xmats)].copy()


def hpj_core(dep, Xmats, Wmats=None, te=False, cluster="entity", robust=False,
             eigen=False, id_order=None):
    """Higher-order jackknife: 4.5*b(N) - 4*b(2N/3) + 0.5*b(N/3).
    Nested entity prefixes; kills 1/N and 1/N^2 bias terms.
    SE = full-sample FE sandwich (Jacobian -> 1 asymptotically)."""
    N = dep.shape[1]
    ids = np.arange(N) if id_order is None else np.asarray(id_order)
    n2, n3 = max(int(round(2 * N / 3)), 2), max(int(round(N / 3)), 2)
    sels = [ids, ids[:n2], ids[:n3]]

    def sub(arr, s):
        return None if arr is None else arr[:, s]

    b1, se1, V1 = fe_core(sub(dep, sels[0]), [sub(x, sels[0]) for x in Xmats],
                          [sub(w, sels[0]) for w in Wmats] if Wmats else None,
                          te, cluster, robust, eigen, return_V=True)
    b2, _ = fe_core(sub(dep, sels[1]), [sub(x, sels[1]) for x in Xmats],
                    [sub(w, sels[1]) for w in Wmats] if Wmats else None, te,
                    cluster, robust, eigen)
    b3, _ = fe_core(sub(dep, sels[2]), [sub(x, sels[2]) for x in Xmats],
                    [sub(w, sels[2]) for w in Wmats] if Wmats else None, te,
                    cluster, robust, eigen)
    beta = 4.5 * b1 - 4.0 * b2 + 0.5 * b3
    return beta, se1


def db_terms(s, y0):
    """Precompute DB correction ingredients (their simul_LP.R tau=0 branch):
    pooled AR(1) rho on demeaned s; b0, sx2 from demeaned (y0, s).
    s, y0: T x N wide (y0 = the 'raw' level-equivalent outcome)."""
    s = np.asarray(s, dtype=float); y0 = np.asarray(y0, dtype=float)
    xd = s - np.nanmean(s, axis=0, keepdims=True)
    yd = y0 - np.nanmean(y0, axis=0, keepdims=True)
    xl, xf = xd[:-1].ravel(), xd[1:].ravel()
    ok = np.isfinite(xl) & np.isfinite(xf)
    terms = dict(rho=np.nan, b0=0.0, sx2=np.nan, s_u2=np.nan)
    if ok.sum() < 10 or np.allclose(xl[ok], 0):
        return terms
    A = np.vstack([xl[ok], np.ones(int(ok.sum()))]).T
    sol, *_ = np.linalg.lstsq(A, xf[ok], rcond=None)
    rho_hat, c0 = sol[0], sol[1]
    s_u2 = np.nanmean((xf[ok] - rho_hat * xl[ok] - c0) ** 2)
    xv, yv = xd.ravel(), yd.ravel()
    m = np.isfinite(xv) & np.isfinite(yv)
    den = np.nansum(xv[m] ** 2)
    if den <= 0:
        return terms
    terms.update(rho=float(rho_hat), b0=float(np.nansum(xv[m] * yv[m]) / den),
                 sx2=float(np.nanmean(xv[m] ** 2)), s_u2=float(s_u2))
    return terms


def db_apply(terms, T, h, coef_fe_h, se_fe_h):
    """DB correction for horizon h: FE_h + b0*f(rho,h)/((T-h)*sx^2),
    f(rho,h) = s_u2*((1-rho^h) - h/(T-h))/(1-rho)^2. SE = FE SE (their choice)."""
    if terms["rho"] is None or not np.isfinite(terms["rho"]):
        return float(coef_fe_h), float(se_fe_h)
    h = np.asarray(h)
    denom = np.maximum(T - h, 1)
    with np.errstate(divide="ignore", invalid="ignore"):
        f_r = terms["s_u2"] * ((1 - terms["rho"] ** h) - h / denom) / (1 - terms["rho"]) ** 2
    return float(coef_fe_h + terms["b0"] * f_r / (denom * terms["sx2"])), float(se_fe_h)


def db_correction(dep, s, coef_fe, se_fe):
    """DB analytical correction (their simul_LP.R, tau=0 branch).
    dep: wide y (T x N); s: wide shock (T x N)."""
    dep = np.asarray(dep, dtype=float)
    coef_fe = np.atleast_1d(np.asarray(coef_fe, dtype=float))
    se_fe = np.atleast_1d(np.asarray(se_fe, dtype=float))
    T, N = dep.shape
    H = len(coef_fe)
    hh = np.arange(H)
    denom = np.maximum(T - hh, 1)
    irf_db = coef_fe.copy()
    # their formula is single-x; caller passes the primary shock column
    # (documented in spec log for multi-shock apps)
    s_j = s
    xd = s_j - np.nanmean(s_j, axis=0, keepdims=True)
    yd = dep - np.nanmean(dep, axis=0, keepdims=True)
    xl, xf = xd[:-1].ravel(), xd[1:].ravel()
    ok = np.isfinite(xl) & np.isfinite(xf)
    if ok.sum() < 10 or np.allclose(xl[ok], 0):
        return irf_db, se_fe.copy()
    A = np.vstack([xl[ok], np.ones(int(ok.sum()))]).T
    sol, *_ = np.linalg.lstsq(A, xf[ok], rcond=None)
    rho_hat, c0 = sol[0], sol[1]
    s_u2 = np.nanmean((xf[ok] - rho_hat * xl[ok] - c0) ** 2)
    xv, yv = xd.ravel(), yd.ravel()
    m = np.isfinite(xv) & np.isfinite(yv)
    den = np.nansum(xv[m] ** 2)
    if den <= 0:
        return irf_db, se_fe.copy()
    b0 = np.nansum(xv[m] * yv[m]) / den
    sx2 = np.nanmean(xv[m] ** 2)
    with np.errstate(divide="ignore", invalid="ignore"):
        f_r = s_u2 * ((1 - rho_hat ** hh) - hh / denom) / (1 - rho_hat) ** 2
    irf_db = coef_fe + b0 * f_r / (denom * sx2)
    return irf_db, se_fe.copy()


def cdols_ols(lhs, s_t, s_lag=None):
    """CD-OLS generic: regress lhs (T x N) on [s_t, s_lag] pooled, entity-clustered.
    Returns (beta_s, se_s) for the s_t coefficient."""
    T, N = lhs.shape
    cols = [s_t] + ([s_lag] if s_lag is not None else [])
    C = np.stack(cols, axis=2)                  # (T, N, K)
    K = C.shape[2]
    Yv = lhs.T.reshape(-1)                      # entity-major
    Xv = C.transpose(1, 0, 2).reshape(N * T, K)
    ok = np.isfinite(Yv) & np.isfinite(Xv).all(axis=1)
    if ok.sum() < 3 * K:
        return np.nan, np.nan
    yv, xv = Yv[ok], Xv[ok]
    ents = np.repeat(np.arange(N), T)[ok]
    coef, *_ = np.linalg.lstsq(xv, yv, rcond=None)
    res = yv - xv @ coef
    Q = xv.T @ xv
    S = np.zeros((K, K))
    for g in range(N):
        m = ents == g
        if m.sum():
            u = xv[m].T @ res[m]
            S += np.outer(u, u)
    try:
        Qinv = np.linalg.inv(Q)
    except np.linalg.LinAlgError:
        Qinv = np.linalg.pinv(Q)
    V = Qinv @ S @ Qinv
    return float(coef[0]), float(np.sqrt(np.clip(V[0, 0], 0, None)))


def cdols_h(y, s, h, include_lag=True, cluster="entity"):
    """CD-OLS (Ugarte-Ruiz 2026 eq. 20) for horizon h, wide inputs.
    LHS: y_{i,t+h} - y_{i,t-1}; RHS: s_{i,t} [, s_{i,t-1}].
    y, s: T x N wide. Returns (beta, se) for the s_t coefficient (+lag)."""
    T, N = y.shape
    tt = np.arange(1, T - h)          # t-1 >= 1 (0-based t from 1); t+h <= T-1
    if len(tt) < 3:
        return np.nan, np.nan
    # Entity-major stacked rows: entity i, time tt (Tt rows each)
    L = y[tt + h, :] - y[tt - 1, :]          # Tt x N
    C = np.stack([s[tt, :]] + ([s[tt - 1, :]] if include_lag else []), axis=2)
    K = C.shape[2]                           # (Tt, N, K)
    Yv = L.T.reshape(-1)                     # N*Tt, entity-major
    Xv = C.transpose(1, 0, 2).reshape(N * len(tt), K)
    ok = np.isfinite(Yv) & np.isfinite(Xv).all(axis=1)
    if ok.sum() < 3 * K:
        return np.full(1, np.nan), np.full(1, np.nan)
    yv, xv = Yv[ok], Xv[ok]
    ents = np.repeat(np.arange(N), len(tt))[ok]
    coef, *_ = np.linalg.lstsq(xv, yv, rcond=None)
    res = yv - xv @ coef
    Q = xv.T @ xv
    S = np.zeros((K, K))
    for g in range(N):
        m = ents == g
        if m.sum():
            u = xv[m].T @ res[m]
            S += np.outer(u, u)
    try:
        Qinv = np.linalg.inv(Q)
    except np.linalg.LinAlgError:
        Qinv = np.linalg.pinv(Q)
    V = Qinv @ S @ Qinv
    return float(coef[0]), float(np.sqrt(np.clip(V[0, 0], 0, None)))
