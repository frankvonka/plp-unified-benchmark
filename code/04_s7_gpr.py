# =============================================================================
# 04_s7_gpr.py -- S7: new application -- global geopolitical risk (GPR) shocks
# =============================================================================
"""
Design (user brief sec. 12-13): identified interactions, not the level.
  Level spec (collinear with time FE, reported as descriptive):
    y_{i,t+h} = a_i + b_h GPR_t + u  -> NOT identified under time FE
  Interaction (the paper's design):
    y_{i,t+h} = a_i + b_h (GPR_t x Exp_i) + controls + u
    y_{i,t+h} = a_i + b_h (GPR_t x CrisisHist_i) + controls + u
  Exp_i: pre-1985 log-GDP z-score (deterministic, see spec log);
  CrisisHist_i: mean financial-crisis years 1965-1984 (CRISIS column mean).
Panel: CS application (192 countries x 1960-2001, GPR official from 1985).
Outcome: cf{h}GRRT_WB (GDP growth), h=1..3 (window ends 2001).
Source: Caldara & Iacoviello (2022) GPR index, official monthly file,
        annualized as the mean of monthly GPR (source column in output).
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
est = _load("01_estimators")

OUT = f"{cfg.RES}/07_gpr"; os.makedirs(OUT, exist_ok=True)
os.makedirs(f"{cfg.CKP}/s7", exist_ok=True)

sp = cfg.APPS["CS"]


def build_panel():
    df = pd.read_csv(f"{cfg.RES}/01_data_audit/CS_clean.csv")
    idc, tc = df.columns[0], df.columns[1]
    g = pd.read_excel(cfg.GPR_XLSX)
    g["year"] = pd.to_datetime(g["month"]).dt.year
    g_ann = g.groupby("year")["GPR"].mean().reset_index().rename(columns={"GPR": "GPR"})
    # exposure: pre-1985 mean log(1+GDPgrowth level proxy) -- the CS file has no
    # GDP level; use pre-1985 mean GDP GROWTH (cf0-equivalent) & crisis history
    # (documented in spec log: exposures are pre-determined, built from 1965-84)
    pre = df[df[tc] < 1985]
    exp = pre.groupby(idc).apply(lambda d: pd.Series({
        "exposure": d["cf1GRRT_WB"].dropna().mean(),
        "CRISIS_exp": d["CRISIS"].fillna(0).mean() if "CRISIS" in d.columns else np.nan,
        "n_pre": d[tc].count()}))
    exp = exp.reset_index().rename(columns={idc: idc})
    exp["exposure"] = (exp["exposure"] - exp["exposure"].mean()) / exp["exposure"].std()
    d = df.merge(g_ann, left_on=tc, right_on="year", how="left")
    d = d.merge(exp[[idc, "exposure", "CRISIS_exp"]], on=idc, how="left")
    d["GPR_x_exp"] = d["GPR"].fillna(0) * d["exposure"].fillna(0)
    d["GPR_x_cris"] = d["GPR"].fillna(0) * d["CRISIS_exp"].fillna(0)
    d["GPR_c"] = d["GPR"].fillna(0)
    d["source"] = "Caldara-Iacoviello GPR official (matteoiacoviello.com), annualized"
    d.to_csv(f"{OUT}/CS_GPR_panel.csv", index=False)
    return d, idc, tc


def build_wide(d, idc, tc, col):
    times = sorted(d[tc].unique()); ids = sorted(d[idc].unique())
    T, N = len(times), len(ids)
    tix = {t: j for j, t in enumerate(times)}; iix = {i: k for k, i in enumerate(ids)}
    m = np.full((T, N), np.nan)
    v = d[col].to_numpy(float)
    r = d[tc].map(tix).to_numpy(); c = d[idc].map(iix).to_numpy()
    ok = np.isfinite(v)
    m[r[ok].astype(int), c[ok].astype(int)] = v[ok]
    return m


def main():
    d, idc, tc = build_panel()
    win = (1985, 2001)
    d2 = d[(d[tc] >= win[0]) & (d[tc] <= win[1])].copy()
    ctrls_cols = ["l1CRISIS", "l2CRISIS", "l3CRISIS", "l4CRISIS",
                  "l1GRRT_WB", "l2GRRT_WB", "l3GRRT_WB", "l4GRRT_WB"]
    specs = [("GPR_c", "level (descriptive, time-FE collinear)"),
             ("GPR_x_exp", "GPR x exposure (identified interaction)"),
             ("GPR_x_cris", "GPR x crisis-history (identified interaction)")]
    results = {}
    for xcol, lab in specs:
        shocks = [build_wide(d2, idc, tc, xcol)]
        ctrls = [build_wide(d2, idc, tc, c) for c in ctrls_cols]
        res = {}
        # CD-OLS eq-20 LHS = y_{t+h} - y_{t-1}.  CS file has l1GRRT_WB (= raw
        # growth y_{t-1} relative to the outcome scale) for h=1; for h>=2 the
        # y_{t-1} base is not distributed -> spec log: LHS = out_h - out_{h-1}.
        base_col_map = {1: "l1GRRT_WB", 2: "cf1GRRT_WB", 3: "cf2GRRT_WB"}
        for h in (1, 2, 3):
            oc = f"cf{h}GRRT_WB"
            dep = build_wide(d2, idc, tc, oc)
            base = build_wide(d2, idc, tc, base_col_map[h])
            out = {}
            b, se = est.fe_core(dep, shocks, ctrls, te=False, cluster="entity", robust=True)
            out["FE"] = (b[0], se[0])
            b, se = est.spj_core(dep, shocks, ctrls, te=False, cluster="entity", robust=True)
            out["SPJ"] = (b[0], se[0])
            b, se = est.hpj_core(dep, shocks, ctrls, te=False, cluster="entity", robust=True)
            out["HPJ"] = (b[0], se[0])
            b, se = est.cdols_ols(dep - base, shocks[0], None)
            out["CDOLS"] = (b, se)
            dbt = est.db_terms(shocks[0], dep)
            bfe, sefe = out["FE"]
            b, se = est.db_apply(dbt, dep.shape[0], h - 1, bfe, sefe)
            out["DB"] = (b, se)
            res[h] = out
        results[xcol] = dict(label=lab, res=res)
    # tidy save
    rows = []
    for xcol, r in results.items():
        for h, out in r["res"].items():
            for e, (b, se) in out.items():
                rows.append(dict(shock=xcol, label=r["label"], h=h, est=e,
                                 beta=b, se=se,
                                 t=(b / se if np.isfinite(se) and se > 0 else np.nan),
                                 lo=b - 1.96 * se if np.isfinite(se) else np.nan,
                                 hi=b + 1.96 * se if np.isfinite(se) else np.nan))
    pd.DataFrame(rows).to_csv(f"{OUT}/gpr_irf.csv", index=False)
    # descriptives + coverage
    cov = dict(
        n_countries=int(d2[idc].nunique()),
        years=f"{win[0]}-{win[1]}",
        gpr_source="Caldara-Iacoviello official GPR (annualized monthly, 1985-2026 file)",
        gpr_nonzero=int((d2["GPR_c"] != 0).sum()),
        exposure_def="pre-1985 mean cf1GRRT_WB, z-scored (spec log: proxy for structural openness to global shocks)",
        crisis_hist_def="pre-1985 mean CRISIS (financial crisis years 1965-84)",
        note="Level spec collinear with time FE in TWFE; interactions are the identified specs")
    json.dump(cov, open(f"{OUT}/gpr_coverage.json", "w"), indent=1)
    json.dump({k: {str(h): {e: [None if not np.isfinite(v[0]) else float(v[0]),
                                None if not np.isfinite(v[1]) else float(v[1])]
                            for e, v in out.items()} for h, out in r["res"].items()}
               for k, r in results.items()},
              open(f"{cfg.CKP}/s7/gpr.json", "w"), indent=1)
    print(pd.DataFrame(rows).query("est in ['FE','SPJ','HPJ']").to_string())
    print("S7 done ->", OUT)


if __name__ == "__main__":
    main()
