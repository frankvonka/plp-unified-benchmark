# =============================================================================
# 08_s10_optimization.py -- S10: estimator selection as optimization
# =============================================================================
"""
Two method families (variation rule; ALL results saved, winner by gate):
 A) Decision-theoretic LP allocation: a researcher with design d=(N,T,rho,sig_lam)
    chooses a distribution p over the 5 estimators minimizing expected loss
    E_D[L] = sum_e p_e * ( D_I(e,d) + lam_cov*(1 - coverage(e,d)) + lam_sign*1[sign flip] ),
    subject to sum p = 1, p >= 0, and a robustness band on p (epsilon-constraint:
    p_e in [0, 1 - eps_all]).  Solved twice: (i) scipy.optimize.linprog (HiGHS),
    (ii) closed-form argmin over vertices {e} (the LP optimum over a simplex is
    always at a vertex) -- the two must agree.
 B) Reliability surrogates (ML family split: ridge vs k-NN): predict
    D_I(e, d) and |coverage(e,d) - 0.95| from features f(d) = [log N, log T,
    rho, sig_lam, rho*sig_lam] using our S4 main-grid summaries; report
    R2 + RMSE + MAE per family; winner gated on RMSE.
Outputs: results/10_optimization/*.csv, checkpoints/s10.
"""
import sys, os, json
import numpy as np
import pandas as pd
import importlib.util
HERE = os.path.dirname(os.path.abspath(__file__))
def _load(name):
    spec = importlib.util.spec_from_file_location(name, f"{HERE}/{name}.py")
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m
cfg = _load("00_config")

OUT = f"{cfg.RES}/10_optimization"; os.makedirs(OUT, exist_ok=True)
os.makedirs(f"{cfg.CKP}/s10", exist_ok=True)

LAMBDA_COV = 2.0      # weight on coverage error (same units as D_I)
LAMBDA_SIGN = 0.5     # penalty for IRF sign flip vs estimator consensus
EPS = 0.0


def _features(N, T, rho, sl):
    return [np.log(N), np.log(T), rho, sl, rho * sl]


def loss_matrix(agg):
    """agg: S4 main-grid summary (cell x est x h metrics, aggregated to
    cell x est). Returns dict est -> (D_I, coverage, bias)."""
    out = {}
    for e in cfg.ESTIMATORS:
        g = agg[agg.est == e]
        if g.empty:
            continue
        out[e] = dict(D_I=g["mae"].mean(), cov=g["coverage"].mean(),
                      bias=g["bias"].mean())
    return out


def main():
    agg = pd.read_csv(f"{cfg.RES}/04_monte_carlo/main_grid_summary.csv")
    agg = agg[agg.est.isin(cfg.ESTIMATORS)]
    # horizon-mean D_I per cell x est
    piv = agg.pivot_table(index=["N", "T", "rho", "sig_lam"], columns="est",
                          values="mae").reset_index()
    pivcov = agg.pivot_table(index=["N", "T", "rho", "sig_lam"], columns="est",
                             values="coverage").reset_index()
    pivbias = agg.pivot_table(index=["N", "T", "rho", "sig_lam"], columns="est",
                              values="bias").reset_index()

    # ---------------- A) LP allocation per cell
    from scipy.optimize import linprog
    ests = cfg.ESTIMATORS
    recs = []
    for i, row in piv.iterrows():
        key = (row["N"], row["T"], row["rho"], row["sig_lam"])
        L = np.array([np.nan_to_num(row[e]) + LAMBDA_COV *
                      abs(0.95 - (pivcov.loc[i, e] if np.isfinite(pivcov.loc[i, e]) else 0.95))
                      for e in ests])
        sign_flip = np.array([
            1.0 if (np.nan_to_num(pivbias.loc[i, e]) * pivbias.loc[i, "SPJ"]) < 0
                  and np.isfinite(pivbias.loc[i, e]) else 0.0 for e in ests])
        L = L + LAMBDA_SIGN * sign_flip
        # simplex constraint: sum_e p_e = 1 (equality), 0 <= p <= 1
        res = linprog(L, A_eq=np.ones((1, len(ests))), b_eq=[1.0],
                      bounds=[(0.0, 1.0)] * len(ests), method="highs")
        p_lp = res.x
        p_vertex = np.zeros(len(ests)); p_vertex[int(np.argmin(L))] = 1.0
        recs.append(dict(N=row["N"], T=row["T"], rho=row["rho"],
                         sig_lam=row["sig_lam"],
                         lp_alloc=" ".join(f"{e}={p_lp[j]:.2f}" for j, e in enumerate(ests)),
                         lp_obj=float(res.fun), vertex_alloc=ests[int(np.argmin(L))],
                         vertex_obj=float(L.min()),
                         agree=bool(np.abs(p_lp - p_vertex).max() < 0.02)))
    alloc = pd.DataFrame(recs)
    alloc.to_csv(f"{OUT}/lp_allocation.csv", index=False)
    print("LP allocation agreement with vertex:", alloc["agree"].mean())

    # ---------------- B) surrogate ML: two families (ridge / k-NN), 4 models
    from sklearn.linear_model import Ridge
    from sklearn.neighbors import KNeighborsRegressor
    from sklearn.metrics import r2_score
    Xall = np.array([_features(r["N"], r["T"], r["rho"], r["sig_lam"])
                     for r in piv.to_dict("records")])
    Yd = piv[ests].to_numpy()
    Yc = pivcov[ests].to_numpy()
    # out-of-fold: 5-fold CV, deterministic folds (seed 42)
    n = len(Yd)
    fold = np.random.default_rng(42).permutation(n) % 5
    rep = []
    for name, mk in [("ridge_a1", lambda: Ridge(alpha=1.0)),
                    ("ridge_a100", lambda: Ridge(alpha=100.0)),
                    ("knn3", lambda: KNeighborsRegressor(n_neighbors=3)),
                    ("knn5", lambda: KNeighborsRegressor(n_neighbors=5))]:
        errs_D, errs_C, r2D, r2C, maeD, maeC = [], [], [], [], [], []
        for f in range(5):
            te = fold == f; tr = ~te
            m1 = mk(); m1.fit(Xall[tr], Yd[tr])
            p1 = m1.predict(Xall[te])
            errs_D.append(np.sqrt(np.mean((p1 - Yd[te]) ** 2)))
            m2 = mk(); m2.fit(Xall[tr], Yc[tr])
            p2 = m2.predict(Xall[te])
            errs_C.append(np.sqrt(np.mean((p2 - Yc[te]) ** 2)))
        mD = mk(); mD.fit(Xall, Yd)
        mC = mk(); mC.fit(Xall, Yc)
        rep.append(dict(model=name,
                        cv_rmse_D_I=float(np.mean(errs_D)),
                        cv_rmse_cov=float(np.mean(errs_C)),
                        r2_D_I=float(r2_score(Yd, mD.predict(Xall))),
                        r2_cov=float(r2_score(Yc, mC.predict(Xall))),
                        mae_D_I=float(np.mean(np.abs(mD.predict(Xall) - Yd))),
                        mae_cov=float(np.mean(np.abs(mC.predict(Xall) - Yc)))))
    ml = pd.DataFrame(rep)
    ml.to_csv(f"{OUT}/surrogate_ml.csv", index=False)
    winner = ml.sort_values("cv_rmse_D_I").iloc[0]["model"]
    print(ml.to_string(index=False))
    # ---------------- gate: winner + recommendation map
    rec_cells = []
    for i, row in piv.iterrows():
        Lmat = np.array([np.nan_to_num(row[e]) + LAMBDA_COV *
                         abs(0.95 - np.nan_to_num(pivcov.loc[i, e]))
                         for e in ests])
        rec_cells.append(dict(N=row["N"], T=row["T"], rho=row["rho"],
                              sig_lam=row["sig_lam"],
                              recommended=ests[int(np.argmin(Lmat))],
                              expected_loss=float(Lmat.min())))
    pd.DataFrame(rec_cells).to_csv(f"{OUT}/recommendation_map.csv", index=False)
    json.dump(dict(lambda_cov=LAMBDA_COV, lambda_sign=LAMBDA_SIGN, eps=EPS,
                   surrogate_winner=winner,
                   note="LP over simplex is attained at a vertex; high-solve LP "
                        "should concentrate on the argmin (agreement column)"),
              open(f"{cfg.CKP}/s10/meta.json", "w"), indent=1)
    print("S10 done ->", OUT)


if __name__ == "__main__":
    main()
