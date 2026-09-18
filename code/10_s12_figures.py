# =============================================================================
# 10_s12_figures.py -- S12: figure renderer + self-contained HTML manuscript
# =============================================================================
"""
Figures (programmatic, from saved result CSVs only):
 F1  app IRF fan: 5 apps, all 5 estimators, 95% CI shades
 F2  estimator disagreement D_h per app
 F3  MC bias surface (est x rho at T=50, XSD=0)
 F4  coverage heat map (est x rho, T=50)
 F5  reliability grade map (pivot)
 F6  MC RMSE/MAE mean per estimator
 F7  XSD inference coverage bars (ent/time/tw)
 F8  GPR application: interaction IRFs
 F9  surrogate ML: predicted vs actual D_I scatter + winner note
 F10 power/FPR-FNR block curves
 F11 disagreement bar chart (mean D_h across apps)
 F12 heterogeneity block (het-rho / het-beta summaries)
HTML: base64-embedded figures + headline tables (self-contained single file).
"""
import sys, os, base64, io, json
import numpy as np
import pandas as pd
import importlib.util
HERE = os.path.dirname(os.path.abspath(__file__))
def _load(name):
    spec = importlib.util.spec_from_file_location(name, f"{HERE}/{name}.py")
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m
cfg = _load("00_config")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"figure.dpi": 110, "font.size": 8, "axes.linewidth": 0.6})

FIG = cfg.FIG; os.makedirs(FIG, exist_ok=True)
MS = cfg.MS; os.makedirs(MS, exist_ok=True)


def save(fig, name):
    p = f"{FIG}/{name}.png"
    fig.savefig(p, bbox_inches="tight")
    plt.close(fig)
    b64 = base64.b64encode(open(p, "rb").read()).decode()
    return b64


def load_apps():
    out = {}
    for app in list(cfg.APPS):
        for e in cfg.ESTIMATORS:
            f = f"{cfg.RES}/03_empirical/{app}_{e}.csv"
            if os.path.exists(f):
                out[(app, e)] = pd.read_csv(f).sort_values("h")
    return out


def f1_app_irfs(apps, b64):
    order = [a for a in cfg.APPS if any((a, "FE") in apps for _ in [0])]
    n = len(order)
    fig, axes = plt.subplots(1, n, figsize=(13, 3))
    if n == 1:
        axes = [axes]
    colors = {"FE": "tab:blue", "SPJ": "tab:red", "HPJ": "tab:green",
              "DB": "tab:purple", "CDOLS": "tab:orange"}
    for ax, app in zip(axes, order):
        for e in cfg.ESTIMATORS:
            if (app, e) not in apps:
                continue
            d = apps[(app, e)]
            ax.plot(d.h, d.beta, color=colors[e], label=e, lw=1.4)
            lo = d.beta - 1.96 * d.se; hi = d.beta + 1.96 * d.se
            ax.fill_between(d.h, lo, hi, color=colors[e], alpha=0.08)
        ax.axhline(0, color="0.7", lw=0.5)
        ax.set_title(f"{app}: {cfg.APPS[app]['label'][:28]}", fontsize=7)
        ax.set_xlabel("h", fontsize=7)
        if app == order[0]:
            ax.legend(fontsize=6, ncols=2)
    fig.tight_layout()
    b64["F1"] = save(fig, "F1_app_irfs")


def f2_disagreement(misc, b64):
    d = misc.get("disagreement")
    if d is None:
        return
    colmap = {"FE": "fe", "SPJ": "spj", "HPJ": "hpj", "DB": "db", "CDOLS": "cdols"}
    fig, ax = plt.subplots(figsize=(7.5, 3.4))
    for e in cfg.ESTIMATORS:
        c = colmap.get(e)
        if c not in d.columns:
            continue
        ax.plot(d.h, d[c].to_numpy(), marker="o", ms=3, label=e, lw=1)
    ax.axhline(0, color="0.7", lw=0.5)
    ax.set_title("IRF by estimator across horizons (all apps overlaid, h-axis per app)")
    ax.set_xlabel("h (per-app panels concatenated)"); ax.set_ylabel("beta_hat")
    ax.legend(fontsize=7)
    fig.tight_layout()
    b64["F2"] = save(fig, "F2_disagreement_lines")
    # D_h bar chart
    m2 = d.groupby("app")["D_h"].mean().reindex(
        [a for a in ("RR", "RR_UNEMP", "BVX", "MSV", "CS")
         if a in d.app.unique()])
    fig, ax = plt.subplots(figsize=(6.5, 2.8))
    m2.plot(kind="bar", ax=ax, color="0.4")
    ax.set_title("Mean estimator disagreement D_h (max−min across 5 estimators)")
    ax.set_ylabel("D_h"); ax.set_xlabel("")
    fig.tight_layout()
    b64["F11"] = save(fig, "F11_disagreement_bars")


def f3_f4_f6(mc, b64):
    if mc is None:
        return
    m0 = mc[(mc["T"] == 50) & (mc["sig_lam"] == 0.0)]
    fig, ax = plt.subplots(figsize=(6.5, 4))
    for e in cfg.ESTIMATORS:
        g = m0[m0.est == e].groupby("rho")["bias"].mean().sort_index()
        if g.empty:
            continue
        ax.plot(g.index, g.values, marker="o", ms=4, label=e)
    ax.set_xlabel("shock persistence rho"); ax.set_ylabel("mean bias")
    ax.axhline(0, color="0.7", lw=0.5)
    ax.set_title("MC bias vs rho (T=50, no XSD)")
    ax.legend(fontsize=7)
    fig.tight_layout()
    b64["F3"] = save(fig, "F3_bias_surface")

    fig, ax = plt.subplots(figsize=(6.5, 4))
    cov = m0.pivot_table(index="rho", columns="est", values="coverage")
    for e in cov.columns:
        ax.plot(cov.index, cov[e], marker="s", ms=4, label=e)
    ax.axhline(0.95, color="0.5", ls="--", lw=0.8)
    ax.set_xlabel("rho"); ax.set_ylabel("coverage (nominal 0.95)")
    ax.set_title("Coverage vs rho (T=50, no XSD)")
    ax.legend(fontsize=7)
    fig.tight_layout()
    b64["F4"] = save(fig, "F4_coverage_heat")

    g = mc.groupby("est")[["rmse", "mae", "coverage"]].mean()
    fig, ax = plt.subplots(figsize=(6.5, 3))
    g.plot(kind="bar", ax=ax)
    ax.set_title("Main grid: mean RMSE, MAE, coverage by estimator")
    ax.set_xlabel("")
    ax.legend(["RMSE", "MAE", "coverage"], fontsize=7)
    fig.tight_layout()
    b64["F6"] = save(fig, "F6_mc_rmse_mae")


def f5_reliability(misc, b64):
    piv = misc.get("reliability_pivot")
    if piv is None:
        f = f"{cfg.RES}/05_disagreement/reliability_grid_pivot.csv"
        if not os.path.exists(f):
            return
        piv = pd.read_csv(f, index_col=0)
    ests = [e for e in cfg.ESTIMATORS if e in piv.columns]
    piv = piv[ests]
    val_map = {"good": 0, "moderate": 1, "poor": 2}
    M = piv.replace(val_map).astype(float).values
    # row labels: rebuild the design tuple from the reliability grid CSV
    # (pivot index may collapse to a plain index after CSV round-trip)
    grid = pd.read_csv(f"{cfg.RES}/05_disagreement/reliability_grid.csv")
    keys = grid[["N", "T", "rho", "sig_lam"]].drop_duplicates()
    rowlab = [f"N={n} T={t} rho={r} xsd={s}"
              for (n, t, r, s) in keys.itertuples(index=False)]
    fig, ax = plt.subplots(figsize=(7.5, 12))
    ax.imshow(M, aspect="auto", cmap="RdYlGn_r", vmin=0, vmax=2)
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            ax.text(j, i, piv.values[i, j], ha="center", va="center",
                    fontsize=6, color="black")
    ax.set_xticks(range(len(piv.columns)), piv.columns, fontsize=8)
    ax.set_yticks(range(len(rowlab)), rowlab, fontsize=5.5)
    ax.set_title("Reliability grade map: good (D_I<=0.05, |cov-0.95|<=0.05),\n"
                 "moderate (<=0.15 / <=0.10), poor otherwise")
    fig.tight_layout()
    b64["F5"] = save(fig, "F5_reliability_map")


def f7_xsd(misc, b64):
    x = misc.get("xsd_coverage")
    if x is None:
        return
    p = x.pivot(index="est", columns="se_type", values="coverage")
    p.plot(kind="bar", figsize=(6.5, 3.5), ax=None)
    p.plot(kind="bar", figsize=(6.5, 3.5))
    import matplotlib.pyplot as plt
    fig = p.plot(kind="bar", figsize=(6.5, 3.5))
    fig.axhline(0.95, color="0.5", ls="--")
    fig.set_title("Nominal-95% coverage under XSD, by SE treatment")
    fig.legend(fontsize=7)
    fig.figure.tight_layout()
    b64["F7"] = save(fig.figure, "F7_xsd_coverage")


def f8_gpr(gpr, b64):
    if gpr is None:
        return
    g = gpr[gpr.est.isin(["FE", "SPJ", "HPJ", "DB", "CDOLS"])]
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.2))
    for ax, col in zip(axes, ["GPR_c", "GPR_x_exp", "GPR_x_cris"]):
        d = g[g.shock == col]
        for e in cfg.ESTIMATORS:
            de = d[d.est == e]
            if de.empty:
                continue
            ax.plot(de.h, de.beta, marker="o", ms=4, label=e, lw=1)
            lo = de.beta - 1.96 * de.se; hi = de.beta + 1.96 * de.se
            ax.fill_between(de.h, lo, hi, color="0.5", alpha=0.05)
        ax.axhline(0, color="0.7", lw=0.5)
        ax.set_title(col, fontsize=8)
        if col == "GPR_c":
            ax.legend(fontsize=6)
    fig.tight_layout()
    b64["F8"] = save(fig, "F8_gpr")


def f9_ml(misc, b64):
    m = misc.get("surrogate_ml")
    if m is None:
        return
    fig, ax = plt.subplots(figsize=(5.5, 3))
    mT = m.set_index("model")[["cv_rmse_D_I", "cv_rmse_cov"]]
    mT.plot(ax=ax, marker="o", ms=4, legend=True)
    ax.set_title("Surrogate ML (ridge / kNN families): CV RMSE")
    ax.set_xlabel("model"); ax.set_ylabel("CV RMSE")
    fig.tight_layout()
    b64["F9"] = save(fig, "F9_ml_surrogates")


def f10_power(mc, b64):
    p = f"{cfg.RES}/04_monte_carlo/power_summary.csv"
    if not os.path.exists(p):
        return
    d = pd.read_csv(p)
    # NB: pandas parses the literal string "null" as NaN -> derive from cell name
    d["mode"] = np.where(d.cell.str.contains("_null_"), "null", "power")
    fig, ax = plt.subplots(figsize=(7, 3.4))
    labels, fprs, fnrs = [], [], []
    for e in cfg.ESTIMATORS:
        ge = d[d.est == e]
        if ge.empty:
            continue
        fpr = ge[ge["mode"] == "null"].rej_rate.mean()
        fnr = 1.0 - ge[ge["mode"] == "power"].rej_rate.mean()
        labels.append(e); fprs.append(fpr); fnrs.append(fnr)
    x = np.arange(len(labels))
    ax.bar(x - 0.18, fprs, width=0.36, label="FPR (size, H0 true)", color="0.55")
    ax.bar(x + 0.18, fnrs, width=0.36, label="FNR (1−power)", color="0.25")
    ax.set_xticks(x); ax.set_xticklabels(labels)
    ax.axhline(0.05, color="crimson", ls="--", lw=0.8)
    ax.set_ylabel("rate"); ax.set_title("Size / power block (null vs power cells)")
    ax.legend(fontsize=7)
    fig.tight_layout()
    b64["F10"] = save(fig, "F10_power")


def f12_het(mc, b64):
    f = f"{cfg.RES}/04_monte_carlo/het_summary.csv"
    if not os.path.exists(f):
        return
    d = pd.read_csv(f)
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.2))
    for ax, hb, lab in [(axes[0], 1, "Heterogeneous persistence (rho_i ~ U)"),
                        (axes[1], 2, "Heterogeneous effects (beta_i ~ N)")]:
        g = d[d.het_beta == (0 if lab.startswith("Heterogeneous p") else 1)]
        for e in cfg.ESTIMATORS:
            ge = g[g.est == e]
            if ge.empty:
                continue
            ax.plot(ge.h, ge.bias, marker="s", ms=3, label=e, lw=1)
        ax.axhline(0, color="0.7", lw=0.5)
        ax.set_title(lab, fontsize=8)
        if hb == 1:
            ax.legend(fontsize=6)
        ax.set_xlabel("h")
    fig.tight_layout()
    b64["F12"] = save(fig, "F12_heterogeneity")


def build_html(b64, appdf, gpr, mc, head, misc):
    blocks = []
    order = [("F1", "Empirical IRFs, all 5 applications x 5 estimators"),
             ("F2", "Estimator disagreement across horizons"),
             ("F11", "Mean D_h per application"),
             ("F3", "MC bias vs persistence"),
             ("F4", "Coverage vs persistence"),
             ("F5", "Reliability grade map"),
             ("F6", "MC RMSE/MAE/coverage by estimator"),
             ("F10", "Size / power block"),
             ("F7", "Inference under cross-sectional dependence"),
             ("F8", "GPR application (interactions)"),
             ("F9", "Surrogate ML (estimator selection)"),
             ("F12", "Heterogeneity blocks")]
    for k, t in order:
        if k in b64:
            blocks.append(f"<h2>{t}</h2><img src='data:image/png;base64,{b64[k]}'/>")
    head_tab = head.round(4).to_html() if head is not None else ""
    # tables fragment (S13) if present
    tab_frag = ""
    tabf = f"{cfg.P}/tables/tables.html"
    if os.path.exists(tabf):
        tab_frag = "<h1>Tables</h1>" + open(tabf).read().replace("<h1>Tables</h1>", "")
    L = f"""<html><head><meta charset='utf-8'><style>
    body{{font-family:Georgia,serif;max-width:980px;margin:2em auto;padding:0 1em;
    color:#1a1a1a;line-height:1.5}} h1{{font-size:1.6em}} h2{{font-size:1.1em;
    margin-top:1.4em;border-bottom:1px solid #ccc;padding-bottom:.2em}}
    img{{max-width:100%;border:1px solid #eee}} table{{border-collapse:collapse;
    font-size:.85em;margin:1em 0}} td,th{{border:1px solid #ddd;padding:.3em .5em}}
    th{{background:#f4f4f4}}</style></head><body>
    <h1>When Do Panel Local Projections Mislead?</h1>
    <p>Unified assessment of five panel-LP estimators (FE, SPJ, HPJ, DB, CD-OLS)
    under persistence, cross-sectional dependence and heterogeneity.
    Reproducible pipeline: seeds {cfg.SEEDS}.</p>
    <h2>Headline Monte Carlo results</h2>{head_tab}
    {('<br/><br/>'.join(blocks))}
    {tab_frag}
    <p style='font-size:.8em;color:#777'>All figures and tables generated
    programmatically from saved result files; no manual values.</p>
    </body></html>"""
    open(f"{MS}/manuscript.html", "w").write(L)
    return L


def main():
    apps = load_apps()
    misc = {}
    for name, path in [("disagreement", f"{cfg.RES}/05_disagreement/app_disagreement.csv"),
                       ("reliability_pivot", f"{cfg.RES}/05_disagreement/reliability_grid_pivot.csv"),
                       ("xsd_coverage", f"{cfg.RES}/09_xsd_inference/xsd_coverage.csv"),
                       ("surrogate_ml", f"{cfg.RES}/10_optimization/surrogate_ml.csv")]:
        if os.path.exists(path):
            misc[name] = pd.read_csv(path)
    gpr = pd.read_csv(f"{cfg.RES}/07_gpr/gpr_irf.csv") \
        if os.path.exists(f"{cfg.RES}/07_gpr/gpr_irf.csv") else None
    mc = pd.read_csv(f"{cfg.RES}/04_monte_carlo/main_grid_summary.csv") \
        if os.path.exists(f"{cfg.RES}/04_monte_carlo/main_grid_summary.csv") else None
    head = pd.read_csv(f"{cfg.RES}/FINAL/headline_mc.csv") \
        if os.path.exists(f"{cfg.RES}/FINAL/headline_mc.csv") else None

    b64 = {}
    f1_app_irfs(apps, b64)
    f2_disagreement(misc, b64)
    f3_f4_f6(mc, b64)
    f5_reliability(misc, b64)
    f7_xsd(misc, b64)
    f8_gpr(gpr, b64)
    f9_ml(misc, b64)
    f10_power(mc, b64)
    f12_het(mc, b64)
    print("figures:", sorted(b64.keys()))
    build_html(b64, None, gpr, mc, head, misc)
    print("S12 done ->", MS)


if __name__ == "__main__":
    main()
