# =============================================================================
# 03_s4_monte_carlo.py -- S4: unified Monte Carlo benchmark (fast paths)
# =============================================================================
"""
DGP (their toy, extended per user brief sec. 7/9/10):
  x_it = delta + rho_i*x_{i,t-1} + lam_i*F_t + mu_it
  y_it = alpha_i + beta_i*x_it + gam_i*F_t + eps_it
  alpha_i = eta*sqrt(T)*mean(x_i) + chi        (their FE construction)
Conventions (their main-text sim, spec log): te=F, robust=F, entity SE,
lagY=0, lagX=1, hmax = min(H, T//3).  Balanced panels -> cut=(T-1)//2.
Metrics: bias, RMSE, MAE, coverage(95%), rejection rate
(FPR when beta_true=0; FNR = 1-rej when beta_true!=0), D_I = horizon-mean
IRF shape distortion.  Seeds 42/7/123.
"""
import sys, os, json, time
import numpy as np
import pandas as pd
import importlib.util
HERE = os.path.dirname(os.path.abspath(__file__))
def _load(name):
    spec = importlib.util.spec_from_file_location(name, f"{HERE}/{name}.py")
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m
cfg = _load("00_config")
est = _load("01_estimators")
fast = _load("03b_mc_estimators")

OUT = f"{cfg.RES}/04_monte_carlo"; os.makedirs(OUT, exist_ok=True)
os.makedirs(f"{cfg.CKP}/s4", exist_ok=True)
H = cfg.HORIZON
COLS = ["cell", "N", "T", "rho", "sig_lam", "het_rho", "het_beta",
        "beta_true", "beta_mean", "seed", "iter", "est", "h",
        "beta_hat", "se", "bias", "sq_err", "abs_err", "cover", "rej_null"]


def gen_panel(rng, N, T, rho, beta, sig_lam=0.0, sig_gam=0.0, het_rho=False,
              het_beta=False, beta_mean=None, sigma_beta=cfg.SIGMA_BETA):
    delta, eta = cfg.DELTA, cfg.ETA
    F = rng.normal(0, 1, T) if sig_lam > 0 else np.zeros(T)
    lam = rng.normal(0, sig_lam, N)
    gam = rng.normal(0, sig_gam, N)
    rhos = np.clip(rng.uniform(0, 2 * rho, N), 0, 0.99) if het_rho else np.full(N, rho)
    betas = rng.normal(beta_mean, sigma_beta, N) if het_beta else np.full(N, beta)
    X = np.zeros((T, N)); Y = np.zeros((T, N))
    for i in range(N):
        x = np.zeros(T)
        x[0] = delta + rhos[i] * rng.normal() + rng.normal() + lam[i] * F[0]
        for t in range(1, T):
            x[t] = delta + rhos[i] * x[t - 1] + rng.normal() + lam[i] * F[t]
        alpha = eta * np.sqrt(T) * x.mean() + rng.normal()
        X[:, i] = x
        Y[:, i] = alpha + betas[i] * x + gam[i] * F + rng.normal(size=T)
    return Y, X, rhos, betas


def true_irf(rhos, beta_true_h, h):
    return np.nanmean(np.asarray(rhos) ** h * np.asarray(beta_true_h))


def run_cell(tag, N, T, rho, sig_lam, sig_gam, het_rho, het_beta,
             beta_true, beta_mean, niter, quiet=True):
    t0 = time.time()
    rows = []
    for seed in cfg.SEEDS:
        rng = np.random.default_rng(seed)
        for it in range(niter):
            Y, X, rhos, betas = gen_panel(rng, N, T, rho, beta_true, sig_lam,
                                          sig_gam, het_rho, het_beta, beta_mean)
            bth = float(np.nanmean(betas)) if het_beta else float(beta_true)
            FE = fast.FEFast(X)
            SPJ = fast.SPJFast(X)
            HPJ = fast.HPJFast(X)
            CDO = fast.CDOLSFast(X, include_lag=True)
            hmax = min(H, T // 3)
            dbt = est.db_terms(X, Y)
            for h in range(hmax + 1):
                dep = np.full((T, N), np.nan)
                tt = np.arange(T - h)
                dep[tt, :] = Y[tt + h, :]
                depS = dep[:, :N]
                b = {}
                b["FE"] = FE.fit_h(depS)
                b["SPJ"] = SPJ.fit_h(depS)
                b["HPJ"] = HPJ.fit_h(depS)
                b["CDOLS"] = CDO.fit_h(Y, h)
                bfe, sfe = b["FE"]
                if np.isfinite(bfe):
                    b["DB"] = est.db_apply(dbt, T, h, bfe, sfe)
                else:
                    b["DB"] = (np.nan, np.nan)
                for e in cfg.ESTIMATORS:
                    be, se = b[e]
                    if not np.isfinite(be):
                        continue
                    t = true_irf(rhos, bth, h)
                    rows.append([tag, N, T, rho, sig_lam, int(het_rho), int(het_beta),
                                 beta_true, beta_mean, seed, it, e, h, be, se,
                                 be - t, (be - t) ** 2, abs(be - t),
                                 int(abs(be - t) <= 1.96 * se if np.isfinite(se) else 0),
                                 int(abs(be / se) > 1.96 if np.isfinite(se) and se > 0 else 0)])
            if not quiet and (it % 50 == 0):
                print(f"  {tag} {seed} it{it} {time.time()-t0:.0f}s", flush=True)
    df = pd.DataFrame(rows, columns=COLS)
    df.to_csv(f"{OUT}/{tag}_trials.csv", index=False)
    g = df.groupby(["est", "h"])
    summ = pd.DataFrame({
        "bias": g["bias"].mean(),
        "rmse": g["sq_err"].mean().map(np.sqrt),
        "mae": g["abs_err"].mean(),
        "coverage": g["cover"].mean(),
        "rej_rate": g["rej_null"].mean(),
        "n": g.size()}).reset_index()
    meta = dict(cell=tag, N=N, T=T, rho=rho, sig_lam=sig_lam, sig_gam=sig_gam,
                het_rho=bool(het_rho), het_beta=bool(het_beta), beta_true=beta_true,
                beta_mean=beta_mean, niter=niter, secs=round(time.time() - t0, 1))
    json.dump(meta, open(f"{cfg.CKP}/s4/{tag}.json", "w"), indent=1)
    print(f"cell {tag} done in {meta['secs']}s ({len(df)} trials)", flush=True)
    return summ


def main(which="all"):
    if which in ("all", "main"):
        parts = []
        for N in cfg.N_GRID:
            for T in cfg.T_GRID:
                for rho in cfg.RHO_GRID:
                    for sl in cfg.XSD_GRID:
                        tag = f"main_N{N}_T{T}_rho{rho}_xsd{int(sl)}"
                        s = run_cell(tag, N, T, rho, sl, 0.5 * sl, False, False,
                                     cfg.BETA0, cfg.BETA0, cfg.NITER_MAIN)
                        s = s.assign(cell=tag, N=N, T=T, rho=rho, sig_lam=sl)
                        s.to_csv(f"{OUT}/{tag}_summary.csv", index=False)
                        parts.append(s)
        pd.concat(parts).to_csv(f"{OUT}/main_grid_summary.csv", index=False)
        print("S4 main grid done", flush=True)
    if which in ("all", "power"):
        parts = []
        for (N, T, rho, sl) in cfg.POWER_CELLS:
            for mode in ("null", "power"):
                bt = 0.0 if mode == "null" else cfg.BETA0
                tag = f"power_{mode}_N{N}_T{T}_rho{rho}_xsd{int(sl)}"
                s = run_cell(tag, N, T, rho, sl, 0.5 * sl, False, False,
                             bt, bt, cfg.NITER_POWER)
                s = s.assign(cell=tag, N=N, T=T, rho=rho, sig_lam=sl, mode=mode)
                s.to_csv(f"{OUT}/{tag}_summary.csv", index=False)
                parts.append(s)
        pd.concat(parts).to_csv(f"{OUT}/power_summary.csv", index=False)
        print("S4 power block done", flush=True)
    if which in ("all", "het"):
        parts = []
        for (N, T, rho, sl) in cfg.HET_RHO_CELLS:
            tag = f"hetrho_N{N}_T{T}_rho{rho}_xsd{int(sl)}"
            s = run_cell(tag, N, T, rho, sl, 0.5 * sl, True, False,
                         cfg.BETA0, cfg.BETA0, cfg.NITER_MAIN)
            s = s.assign(cell=tag, N=N, T=T, rho=rho, sig_lam=sl, het_rho=1)
            s.to_csv(f"{OUT}/{tag}_summary.csv", index=False); parts.append(s)
        for (N, T, rho, sl) in cfg.HET_BETA_CELLS:
            tag = f"hetbeta_N{N}_T{T}_rho{rho}_xsd{int(sl)}"
            s = run_cell(tag, N, T, rho, sl, 0.5 * sl, False, True,
                         cfg.BETA0, cfg.BETA0, cfg.NITER_MAIN)
            s = s.assign(cell=tag, N=N, T=T, rho=rho, sig_lam=sl, het_beta=1)
            s.to_csv(f"{OUT}/{tag}_summary.csv", index=False); parts.append(s)
        pd.concat(parts).to_csv(f"{OUT}/het_summary.csv", index=False)
        print("S4 het block done", flush=True)


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    main(which)
