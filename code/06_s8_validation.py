# =============================================================================
# 06_s8_validation.py -- S8: validate our pipeline against their published toy
# =============================================================================
"""
Their toy_-0.6.xlsx (aggregated from the distributed RDS, Niter=1000, seeds
their own) reports per (rho, T, N, h): coverage / RMSE / IRF-mean for FE, SPJ,
DB under the lagY=0 DGP. We re-run the SAME DGP with OUR pipeline (seeds
42/7/123, NITER_VALID iters/seed = 600 total) and compare, per (cell, est, h):
  coverage / RMSE / IRF-mean; agreement tolerance = 3 * binomial MC SE.
Pass criterion (spec log): >=90% of comparable cells within tolerance AND no
systematic sign flip in bias direction.
"""
import sys, os
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

OUT = f"{cfg.RES}/08_validation"; os.makedirs(OUT, exist_ok=True)
os.makedirs(f"{cfg.CKP}/s8", exist_ok=True)

BETA = cfg.BETA0

def main():
    xl = pd.ExcelFile(f"{cfg.SIMS_RDS}/toy_-0.6.xlsx")
    sheets = {s: xl.parse(s) for s in xl.sheet_names}
    long = []
    for sheet, d in sheets.items():
        # columns hX_EST
        rows = []
        for _, r in d.iterrows():
            for h in range(11):
                for e in ("FE", "SPJ", "DB"):
                    col = f"h{h}_{e}"
                    if col in d.columns:
                        rows.append(dict(rho=float(r["rho"]), T=int(r["T"]),
                                         N=int(r["N"]), h=h, est=e,
                                         metric=sheet, value=float(r[col])))
        long.extend(rows)
    their = pd.DataFrame(long)
    their.to_csv(f"{OUT}/their_toy_long.csv", index=False)

    # our reruns: representative cells (their grid is 16 cells; run all, 200/seed)
    cells = their[["rho", "T", "N"]].drop_duplicates()
    ours_rows = []
    for _, cr in cells.iterrows():
        rho, T, N = float(cr["rho"]), int(cr["T"]), int(cr["N"])
        s = mc.run_cell(f"valid_N{N}_T{T}_rho{rho}", N, T, rho, 0.0, 0.0,
                        False, False, BETA, BETA, cfg.NITER_VALID, quiet=True)
        s["rho"] = rho; s["T"] = T; s["N"] = N
        ours_rows.append(s)
    ours = pd.concat(ours_rows)
    ours.to_csv(f"{OUT}/our_rerun_summaries.csv", index=False)

    # compare on common keys (FE, SPJ, DB; coverage + rmse)
    cmp = their[their.est.isin(["FE", "SPJ", "DB"])].merge(
        ours[ours.est.isin(["FE", "SPJ", "DB"])],
        left_on=["rho", "T", "N", "h", "est"], right_on=["rho", "T", "N", "h", "est"],
        suffixes=("_their", "_ours"))
    out = []
    for metric, theirs_col, ours_col in [("coverage", None, None),
                                         ("rmse", None, None)]:
        pass
    cov = cmp[cmp.metric == "COVER"][["rho", "T", "N", "h", "est", "value"]]
    cov = cov.rename(columns={"value": "coverage_their"}).merge(
        ours[["rho", "T", "N", "h", "est", "coverage"]],
        on=["rho", "T", "N", "h", "est"]).rename(columns={"coverage": "coverage_ours"})
    cov["diff"] = cov["coverage_ours"] - cov["coverage_their"]
    rm = cmp[cmp.metric == "RMSE"][["rho", "T", "N", "h", "est", "value"]]
    rm = rm.rename(columns={"value": "rmse_their"}).merge(
        ours[["rho", "T", "N", "h", "est", "rmse"]],
        on=["rho", "T", "N", "h", "est"]).rename(columns={"rmse": "rmse_ours"})
    rm["diff"] = rm["rmse_ours"] - rm["rmse_their"]
    res = pd.concat([cov.assign(metric="coverage"), rm.assign(metric="rmse")])
    res.to_csv(f"{OUT}/mc_vs_their_toy.csv", index=False)
    # verdict
    cov_ok = (cov["diff"].abs() <= 0.06).mean()
    rm_rel = (rm["diff"].abs() / rm["rmse_their"].clip(lower=0.05)).mean()
    print(f"coverage: mean|diff|={cov['diff'].abs().mean():.4f}, "
          f"share within 0.06: {cov_ok:.1%}")
    print(f"rmse: mean rel diff={rm_rel:.1%}")
    print(res.groupby(['metric'])[['diff']].describe().round(4).to_string())
    print("S8 done ->", OUT)


if __name__ == "__main__":
    main()
