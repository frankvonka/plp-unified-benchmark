# =============================================================================
# 05_s5_disagreement.py -- S5: estimator-disagreement + reliability regions
# =============================================================================
"""
Outputs:
 1. results/05_disagreement/app_disagreement.csv
    D_h = max_j beta_hat_{h,j} - min_j beta_hat_{h,j} across the 5 estimators,
    per app per horizon (empirical disagreement of IRFs).
 2. results/05_disagreement/reliability_grid.csv
    MC grid: per cell, horizon-averaged metrics per estimator -> 'good/moderate/
    poor' classification by |bias| and coverage distance (map of where each
    estimator's IRF+CI can be trusted).
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

OUT = f"{cfg.RES}/05_disagreement"; os.makedirs(OUT, exist_ok=True)
os.makedirs(f"{cfg.CKP}/s5", exist_ok=True)


def app_disagreement():
    rows = []
    for app in cfg.APPS:
        irfs = {}
        for e in cfg.ESTIMATORS:
            f = f"{cfg.RES}/03_empirical/{app}_{e}.csv"
            if os.path.exists(f):
                d = pd.read_csv(f)
                irfs[e] = d.set_index("h")["beta"]
        if not irfs:
            continue
        M = pd.DataFrame(irfs)
        D = M.max(axis=1) - M.min(axis=1)
        for h in M.index:
            rows.append(dict(app=app, h=h,
                             D_h=D[h],
                             fe=M.loc[h, "FE"] if "FE" in M else np.nan,
                             spj=M.loc[h, "SPJ"] if "SPJ" in M else np.nan,
                             hpj=M.loc[h, "HPJ"] if "HPJ" in M else np.nan,
                             db=M.loc[h, "DB"] if "DB" in M else np.nan,
                             cdols=M.loc[h, "CDOLS"] if "CDOLS" in M else np.nan,
                             argmax=M.loc[h].abs().idxmax()))
    df = pd.DataFrame(rows)
    df.to_csv(f"{OUT}/app_disagreement.csv", index=False)
    print(df.groupby("app")["D_h"].mean().round(3).to_string())
    return df


def reliability_grid():
    f = f"{cfg.RES}/04_monte_carlo/main_grid_summary.csv"
    if not os.path.exists(f):
        print("main grid not ready; skip"); return None
    g = pd.read_csv(f)
    # aggregate horizons per cell x est (D_I = mean abs err across horizons)
    agg = g.groupby(["cell", "N", "T", "rho", "sig_lam", "est"]).apply(
        lambda d: pd.Series({
            "D_I": d["mae"].mean(),
            "bias_hmean": d["bias"].mean(),
            "coverage_hmean": d["coverage"].mean(),
            "rmse_hmean": np.sqrt((d["rmse"] ** 2).mean())})).reset_index()
    def grade(r):
        if r["D_I"] <= 0.05 and abs(r["coverage_hmean"] - 0.95) <= 0.05:
            return "good"
        if r["D_I"] <= 0.15 and abs(r["coverage_hmean"] - 0.95) <= 0.10:
            return "moderate"
        return "poor"
    agg["grade"] = agg.apply(grade, axis=1)
    agg.to_csv(f"{OUT}/reliability_grid.csv", index=False)
    piv = agg.pivot_table(index=["N", "T", "rho", "sig_lam"], columns="est",
                          values="grade", aggfunc="first")
    piv.to_csv(f"{OUT}/reliability_grid_pivot.csv")
    print(piv.head(12).to_string())
    return agg


if __name__ == "__main__":
    app_disagreement()
    reliability_grid()
