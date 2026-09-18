# =============================================================================
# 07_s9_xsd_inference.py -- S9: inference under cross-sectional dependence
# =============================================================================
"""
User brief sec. 9: XSD across countries invalidates entity-clustered SEs.
Experiment (4 representative cells from the main grid; XSD=1 factors in x AND y):
  For FE, SPJ, HPJ, CDOLS compare nominal-95% coverage of:
    ent  = entity-clustered SE
    time = time-clustered (Driscoll-Kraay lag-0 kernel) SE
    tw   = two-way cluster SE (Cameron-Miller)
  on the SAME draws. Question: which SE treatment restores valid inference?
Outputs: results/09_xsd_inference/xsd_coverage.csv, checkpoint s9.
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
mc = _load("03_s4_monte_carlo")
fast = _load("03b_mc_estimators")

OUT = f"{cfg.RES}/09_xsd_inference"; os.makedirs(OUT, exist_ok=True)
os.makedirs(f"{cfg.CKP}/s9", exist_ok=True)

CELLS = [(20, 20, 0.8, 1.0), (100, 50, 0.95, 1.0), (20, 50, 0.8, 1.0),
         (200, 50, 0.2, 1.0)]
COLS = ["cell", "N", "T", "rho", "seed", "iter", "h", "truth",
        "b_FE", "se_FE_ent", "se_FE_time", "se_FE_tw",
        "b_SPJ", "se_SPJ_ent", "se_SPJ_time", "se_SPJ_tw",
        "b_HPJ", "se_HPJ_ent", "se_HPJ_time", "se_HPJ_tw",
        "b_CDO", "se_CDO_ent"]


def main():
    parts = []
    for (N, T, rho, sl) in CELLS:
        tag = f"xsd_N{N}_T{T}_rho{rho}"
        hmax = min(cfg.HORIZON, T // 3)
        rows = []
        for seed in cfg.SEEDS:
            rng = np.random.default_rng(seed + 999)
            for it in range(cfg.NITER_POWER):
                Y, X, rhos, betas = mc.gen_panel(rng, N, T, rho, cfg.BETA0,
                                                 1.0, 0.5, False, False, cfg.BETA0)
                HPJ = fast.HPJFast(X)
                CDO = fast.CDOLSFast(X)
                for h in range(hmax + 1):
                    dep = np.full((T, N), np.nan)
                    tt = np.arange(T - h)
                    dep[tt, :] = Y[tt + h, :]
                    truth = cfg.BETA0 * rho ** h
                    rec = [tag, N, T, rho, seed, it, h, truth]
                    b, se_e, se_t, se_tw = fast.fe_se_variants(dep, X)
                    rec += [b, se_e, se_t, se_tw]
                    b, se_e, se_t, se_tw = fast.spj_se_variants(dep, X)
                    rec += [b, se_e, se_t, se_tw]
                    bh, seh = HPJ.fit_h(dep)
                    # HPJ SE = FE sandwich (same 3 treatments)
                    _, se_e2, se_t2, se_tw2 = fast.fe_se_variants(dep, X, beta=bh)
                    rec += [bh, se_e2, se_t2, se_tw2]
                    bc, sec = CDO.fit_h(Y, h)
                    rec += [bc, sec]
                    rows.append(rec)
        d = pd.DataFrame(rows, columns=COLS)
        d.to_csv(f"{OUT}/{tag}_trials.csv", index=False)
        parts.append(d)
        print(f"cell {tag} done ({len(d)} rows)", flush=True)
    df = pd.concat(parts)
    df.to_csv(f"{OUT}/xsd_all_trials.csv", index=False)
    recs = []
    for e, bcol, scols in [
            ("FE", "b_FE", ["se_FE_ent", "se_FE_time", "se_FE_tw"]),
            ("SPJ", "b_SPJ", ["se_SPJ_ent", "se_SPJ_time", "se_SPJ_tw"]),
            ("HPJ", "b_HPJ", ["se_HPJ_ent", "se_HPJ_time", "se_HPJ_tw"]),
            ("CDOLS", "b_CDO", ["se_CDO_ent"])]:
        for c in scols:
            se = df[c].to_numpy()
            ok = np.isfinite(se) & (se > 0)
            cov = ((np.abs(df[bcol] - df["truth"]) <= 1.96 * se) & ok).sum() / max(ok.sum(), 1)
            recs.append(dict(est=e, se_type=c.split("_")[-1], coverage=cov, n=int(ok.sum())))
    out = pd.DataFrame(recs)
    out.to_csv(f"{OUT}/xsd_coverage.csv", index=False)
    print(out.pivot_table(index="est", columns="se_type", values="coverage").round(3).to_string())
    json.dump(dict(cells=CELLS, niter=cfg.NITER_POWER, note="tw=Cameron-Miller two-way; time=DK lag-0"),
              open(f"{cfg.CKP}/s9/meta.json", "w"), indent=1)
    print("S9 done ->", OUT)


if __name__ == "__main__":
    main()
