# =============================================================================
# 12_s13_tables.py -- S13: publication tables (T1-T12)
# =============================================================================
"""
Builds 12 publication-style tables from saved result files (no manual values):
  T1  Data audit summary (5 apps)
  T2  Descriptive statistics
  T3  CS crisis IRFs, 5 estimators x h (beta (se))
  T4  RR shock IRFs, 5 estimators x h
  T5  Estimator disagreement D_h by app
  T6  MC main grid: horizon-mean bias / RMSE / coverage by T x XSD x estimator
  T7  MC coverage by persistence rho (T=50, no XSD)
  T8  Size / power block (FPR, FNR by estimator)
  T9  Inference under XSD: coverage by SE treatment x estimator
  T10 GPR application: interaction IRFs by estimator
  T11 Reliability grade map (96 cells x 5 estimators)
  T12 Surrogate ML families: R2 + RMSE + MAE (together)
Outputs: tables/T*.csv + T*.tex + tables/TABLES.md; HTML fragment for the
manuscript (tables/tables.html).
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

TAB = f"{cfg.P}/tables"; os.makedirs(TAB, exist_ok=True)
MD = []
LATEX_OK = True


def emit(tid, title, df, fmt=None, note=""):
    """Save one table as CSV + TeX, append markdown block."""
    df.to_csv(f"{TAB}/{tid}.csv", index=False)
    try:
        with open(f"{TAB}/{tid}.tex", "w") as f:
            f.write(df.to_latex(index=False, float_format="%.4f", caption=title,
                                label=f"tab:{tid}"))
    except Exception:
        pass
    # markdown rendering with optional per-column formatting
    d = df.copy()
    if fmt:
        for col, spec in fmt.items():
            if col in d.columns:
                d[col] = d[col].map(lambda v: spec.format(v) if pd.notna(v) else "—")
    MD.append(f"\n### {tid}. {title}\n")
    MD.append(d.to_markdown(index=False))
    if note:
        MD.append(f"\n_Note: {note}_")
    print(f"{tid}: {len(d)} rows")


def cell_pm(b, se):
    if pd.isna(b) or pd.isna(se):
        return "—"
    return f"{b:.3f} ({se:.3f})"


def main():
    # ---------------- T1 data audit
    a = pd.read_csv(f"{cfg.CKP}/s1/data_audit_summary.csv", index_col=0)
    t1 = a.reset_index().rename(columns={"index": "app"})
    cols = ["app", "n_rows", "n_units", "period_range", "n_missing_cells",
            "balanced", "median_T", "te", "twoc", "n_controls"]
    emit("T1", "Data audit: application panels (Mei-Sheng-Sheng replication data)",
         t1[[c for c in cols if c in t1.columns]])

    # ---------------- T2 descriptives
    d2 = pd.read_csv(f"{cfg.RES}/02_descriptives/descriptives.csv")
    emit("T2", "Descriptive statistics (shock and first-outcome variables)",
         d2.round(3))

    # ---------------- T3/T4 empirical IRFs
    irf = pd.read_csv(f"{cfg.RES}/FINAL/empirical_irfs.csv")
    for tid, app, title in [
            ("T3", "CS", "Crisis -> GDP growth IRFs, Cerra-Saxena panel (beta, cluster SE)"),
            ("T4", "RR", "Tax-shock -> GDP IRFs, Romer-Romer panel (beta, cluster SE)")]:
        s = irf[irf.app == app].pivot(index="h", columns="est", values=["beta", "se"])
        out = pd.DataFrame({"h": s.index})
        for e in cfg.ESTIMATORS:
            if ("beta", e) in s.columns:
                out[e] = [cell_pm(b, se) for b, se in
                          zip(s[("beta", e)], s[("se", e)])]
        emit(tid, title, out)

    # ---------------- T5 disagreement
    dis = pd.read_csv(f"{cfg.RES}/05_disagreement/app_disagreement.csv")
    t5 = dis.groupby("app").agg(
        mean_Dh=("D_h", "mean"), max_Dh=("D_h", "max")).round(3).reset_index()
    emit("T5", "Estimator disagreement: D_h = max-min across 5 estimators",
         t5, note="D_h computed per horizon then averaged (mean) / maxed.")

    # ---------------- T6 MC main grid by T x XSD
    g = pd.read_csv(f"{cfg.RES}/04_monte_carlo/main_grid_summary.csv")
    agg = g.groupby(["T", "sig_lam", "est"]).agg(
        bias=("bias", "mean"), rmse=("rmse", "mean"), mae=("mae", "mean"),
        coverage=("coverage", "mean")).round(4).reset_index()
    emit("T6", "Monte Carlo main grid: horizon-mean metrics by T, XSD, estimator",
         agg)

    # ---------------- T7 coverage by rho
    m0 = g[(g["T"] == 50) & (g.sig_lam == 0.0)]
    t7 = m0.pivot_table(index="rho", columns="est", values="coverage").round(3)
    t7 = t7.reset_index()
    emit("T7", "Coverage by shock persistence (T=50, no XSD; nominal 0.95)",
         t7, note="Horizon-averaged coverage per cell; h = 0..min(10, T/3).")

    # ---------------- T8 size/power
    p = pd.read_csv(f"{cfg.RES}/04_monte_carlo/power_summary.csv")
    mode = p.cell.str.contains("_null_").map({True: "null", False: "power"})
    p["mode"] = mode
    t8 = p.groupby(["est", "mode"]).rej_rate.mean().round(4).unstack(0)
    t8 = t8.T.reset_index().rename(columns={"est": "estimator"})
    if "null" in t8.columns and "power" in t8.columns:
        t8["FNR"] = (1 - t8["power"]).round(4)
        t8 = t8.rename(columns={"null": "FPR (size)"})[["estimator", "FPR (size)", "FNR"]]
    emit("T8", "Size and power: rejection rates when H0 true / H1 true",
         t8, note="Null cells: beta=0; power cells: beta=-0.6. 150 iters x 3 seeds per design point.")

    # ---------------- T9 XSD inference
    x = pd.read_csv(f"{cfg.RES}/09_xsd_inference/xsd_coverage.csv")
    t9 = x.pivot(index="est", columns="se_type", values="coverage").round(3).reset_index()
    emit("T9", "Coverage under cross-sectional dependence by SE treatment",
         t9, note="ent = entity cluster; time = DK lag-0 kernel; tw = two-way cluster (Cameron-Miller).")

    # ---------------- T10 GPR application
    gpr = pd.read_csv(f"{cfg.RES}/07_gpr/gpr_irf.csv")
    keep = gpr[gpr.shock.isin(["GPR_x_exp", "GPR_x_cris"])]
    rows = []
    for (shock, h), s in keep.groupby(["shock", "h"]):
        row = {"shock": shock, "h": h}
        for e in cfg.ESTIMATORS:
            se_ = s[s.est == e]
            if len(se_):
                row[e] = cell_pm(se_.beta.iloc[0], se_.se.iloc[0])
        rows.append(row)
    emit("T10", "GPR shocks x predetermined exposures (identified interactions)",
         pd.DataFrame(rows), note="Level GPR effect is collinear with time effects and not reported as identified.")

    # ---------------- T11 reliability map
    rel = pd.read_csv(f"{cfg.RES}/05_disagreement/reliability_grid.csv")
    piv = rel.pivot_table(index=["N", "T", "rho", "sig_lam"], columns="est",
                          values="grade", aggfunc="first").reset_index()
    emit("T11", "Reliability grades: good / moderate / poor (96 designs)",
         piv, note="good: D_I<=0.05 and |coverage-0.95|<=0.05; moderate: <=0.15 / <=0.10; poor otherwise.")

    # ---------------- T12 surrogate ML
    ml = pd.read_csv(f"{cfg.RES}/10_optimization/surrogate_ml.csv")
    emit("T12", "Estimator-selection surrogates: ridge vs k-NN (R2, RMSE, MAE together)",
         ml.round(4), note="5-fold CV, folds seeded 42. Winner by CV RMSE on D_I: "
         f"{ml.sort_values('cv_rmse_D_I').iloc[0].model}.")

    # ---------------- write TABLES.md + HTML fragment
    header = ("# Tables\n\nAll values generated from saved result files "
              "(seeds 42/7/123); no manually typed numbers.\n")
    open(f"{TAB}/TABLES.md", "w").write(header + "\n".join(MD))
    try:
        html = pd.read_markdown if False else None
    except Exception:
        pass
    # HTML fragment via markdown -> minimal conversion (tables only)
    frag = ["<h1>Tables</h1>"]
    cur = []
    for block in "\n".join(MD).split("\n### "):
        if not block.strip():
            continue
        lines = block.split("\n")
        title = lines[0].strip()
        md_rows = [l for l in lines[1:] if l.strip().startswith("|")]
        if not md_rows:
            continue
        rows = [[c.strip() for c in r.strip("|").split("|")] for r in md_rows]
        header_cells, body = rows[0], rows[2:]
        frag.append(f"<h2>{title}</h2><table><tr>" +
                    "".join(f"<th>{c}</th>" for c in header_cells) + "</tr>" +
                    "".join("<tr>" + "".join(f"<td>{c}</td>" for c in r) + "</tr>"
                            for r in body) + "</table>")
    open(f"{TAB}/tables.html", "w").write("\n".join(frag))
    print(f"S13 done: {len(list(x for x in os.listdir(TAB) if x.startswith('T'))) // 3} tables -> {TAB}")


if __name__ == "__main__":
    main()
