# =============================================================================
# 09_s11_final.py -- S11: freeze FINAL results + ALL_RESULTS.md + spec log
# =============================================================================
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

FINAL = f"{cfg.RES}/FINAL"; os.makedirs(FINAL, exist_ok=True)


def collect_apps():
    rows = []
    for app in cfg.APPS:
        for e in cfg.ESTIMATORS:
            f = f"{cfg.RES}/03_empirical/{app}_{e}.csv"
            if os.path.exists(f):
                d = pd.read_csv(f)
                for _, r in d.iterrows():
                    rows.append(dict(app=app, est=e, h=int(r["h"]), beta=r["beta"],
                                     se=r["se"], lo=r["lo"], hi=r["hi"]))
    appdf = pd.DataFrame(rows)
    appdf.to_csv(f"{FINAL}/empirical_irfs.csv", index=False)
    return appdf


def collect_gpr():
    f = f"{cfg.RES}/07_gpr/gpr_irf.csv"
    if os.path.exists(f):
        d = pd.read_csv(f)
        d.to_csv(f"{FINAL}/gpr_irfs.csv", index=False)
        return d
    return None


def collect_mc():
    fs = [f"{cfg.RES}/04_monte_carlo/main_grid_summary.csv",
          f"{cfg.RES}/04_monte_carlo/power_summary.csv",
          f"{cfg.RES}/04_monte_carlo/het_summary.csv"]
    out = []
    for f in fs:
        if os.path.exists(f):
            d = pd.read_csv(f)
            out.append(d)
    if out:
        m = pd.concat(out)
        m.to_csv(f"{FINAL}/mc_grid.csv", index=False)
        return m
    return None


def collect_misc():
    misc = {}
    for name, path in [("disagreement", f"{cfg.RES}/05_disagreement/app_disagreement.csv"),
                       ("reliability", f"{cfg.RES}/05_disagreement/reliability_grid.csv"),
                       ("xsd_coverage", f"{cfg.RES}/09_xsd_inference/xsd_coverage.csv"),
                       ("lp_allocation", f"{cfg.RES}/10_optimization/lp_allocation.csv"),
                       ("surrogate_ml", f"{cfg.RES}/10_optimization/surrogate_ml.csv"),
                       ("validation_toy", f"{cfg.RES}/08_validation/mc_vs_their_toy.csv")]:
        if os.path.exists(path):
            misc[name] = pd.read_csv(path)
    return misc


def main():
    t0 = time.time()
    appdf = collect_apps()
    gpr = collect_gpr()
    mc = collect_mc()
    misc = collect_misc()

    # headline numbers (R2+RMSE+MAE where applicable)
    head = []
    if mc is not None:
        for e in cfg.ESTIMATORS:
            g = mc[mc.est == e]
            head.append(dict(est=e, n_cells=int(g.cell.nunique()),
                             mean_abs_bias=float(g.bias.mean()),
                             mean_RMSE=float(g.rmse.mean()),
                             mean_MAE=float(g.mae.mean()),
                             mean_coverage=float(g.coverage.mean())))
    pd.DataFrame(head).to_csv(f"{FINAL}/headline_mc.csv", index=False)

    # ALL_RESULTS.md (user's requested dump)
    L = ["# plp-unified-benchmark — ALL RESULTS",
         f"_generated {time.strftime('%Y-%m-%d %H:%M UTC')}; seeds {cfg.SEEDS}_" + "\n"]
    L.append("\n## 1. Empirical applications (5 estimators x 5 apps)\n")
    if appdf is not None:
        piv = appdf.pivot_table(index=["app", "h"], columns="est", values="beta")
        L.append(piv.round(4).to_string())
    L.append("\n## 2. GPR application\n")
    if gpr is not None:
        L.append(gpr[gpr.est.isin(["FE", "SPJ", "HPJ"])].round(4).to_string(index=False))
    L.append("\n## 3. Monte Carlo benchmark\n")
    if head:
        L.append(pd.DataFrame(head).round(4).to_string(index=False))
    if mc is not None:
        L.append("\n### coverage by estimator (main grid)\n")
        L.append(mc.groupby("est").coverage.mean().round(4).to_string())
    L.append("\n## 4. Estimator disagreement (empirical D_h, app means)\n")
    if "disagreement" in misc:
        L.append(misc["disagreement"].groupby("app")["D_h"].mean().round(3).to_string())
    L.append("\n## 5. XSD inference: coverage by SE treatment\n")
    if "xsd_coverage" in misc:
        L.append(misc["xsd_coverage"].pivot(index="est", columns="se_type", values="coverage").round(4).to_string())
    L.append("\n## 6. Optimization / surrogate ML\n")
    if "lp_allocation" in misc:
        L.append("LP-vs-vertex agreement: " + str(misc["lp_allocation"].agree.mean()))
    if "surrogate_ml" in misc:
        L.append(misc["surrogate_ml"].round(4).to_string(index=False))
    L.append("\n## 7. Validation vs their published toy\n")
    if "validation_toy" in misc:
        v = misc["validation_toy"]
        v2 = v[v.metric == "coverage"]
        L.append(f"cells={v.rho.nunique()}, mean|coverage diff|={np.abs(v2['diff']).mean():.4f}")
    L.append("\n## 8. Publication tables (S13)\n")
    tabdir = f"{cfg.P}/tables"
    if os.path.isdir(tabdir):
        tabs = sorted(f for f in os.listdir(tabdir) if f.startswith("T") and f.endswith(".csv"))
        for t in tabs:
            d = pd.read_csv(f"{tabdir}/{t}")
            L.append(f"\n### {t.replace('.csv','')}\n")
            L.append(d.head(12).to_markdown(index=False))
        L.append("\n_full tables in tables/ (CSV + TeX); all embedded in manuscript.html_")
    open(f"{cfg.RES}/ALL_RESULTS.md", "w").write("\n".join(L))
    print("S11 done ->", FINAL, f"({time.time()-t0:.0f}s)")


if __name__ == "__main__":
    main()
