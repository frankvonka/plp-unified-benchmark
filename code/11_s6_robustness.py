# =============================================================================
# 11_s6_robustness.py -- S6: robustness of empirical IRFs
# =============================================================================
"""
Two pre-registered robustness checks on the benchmark apps (spec log):
 R1. Leave-one-country-out (LOCO): re-estimate FE and SPJ dropping each
     country in turn; report the max |delta beta_h| across dropped countries
     (influence diagnostic). Apps: RR, CS (largest + smallest panels).
 R2. Alternative SE treatments: FE with entity vs two-way clustering
     (where feasible) -- does significance ordering FE<SPJ survive?
Outputs: results/06_robustness/{app}_loco.csv, {app}_se_alt.csv
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

OUT = f"{cfg.RES}/06_robustness"; os.makedirs(OUT, exist_ok=True)
os.makedirs(f"{cfg.CKP}/s6", exist_ok=True)


def wide(df, col):
    idc, tc = df.columns[0], df.columns[1]
    times = sorted(df[tc].unique()); ids = sorted(df[idc].unique())
    T, N = len(times), len(ids)
    tix = {t: j for j, t in enumerate(times)}; iix = {i: k for k, i in enumerate(ids)}
    m = np.full((T, N), np.nan)
    v = df[col].to_numpy(float)
    r = df[tc].map(tix).to_numpy(); c = df[idc].map(iix).to_numpy()
    ok = np.isfinite(v)
    m[r[ok].astype(int), c[ok].astype(int)] = v[ok]
    return m, ids


def loco(app, max_h=6):
    sp = cfg.APPS[app]
    df = pd.read_csv(f"{cfg.RES}/01_data_audit/{app}_clean.csv")
    idc, tc = df.columns[0], df.columns[1]
    if sp["y_scale"] != 1.0:
        for oc in sp["outcome_cols"]:
            if oc in df.columns:
                df[oc] = df[oc] * sp["y_scale"]
    dep0, ids = wide(df, sp["outcome_cols"][min(max_h - (0 if sp["h_list"][0] == 0 else 1), len(sp["outcome_cols"]) - 1)])
    shocks = [wide(df, c)[0] for c in sp["shock_cols"]]
    ctrls = [wide(df, c)[0] for c in sp["control_cols"]]
    opts = dict(te=sp["two_way_fe"], cluster=("twoway" if sp["two_way_cluster"] else "entity"),
                robust=sp["robust_scale"], eigen=sp["eigen_fix"])
    rows = []
    N = dep0.shape[1]
    keep = list(range(N))
    for drop in range(N):
        sel = [i for i in keep if i != drop]
        d = dep0[:, sel]; sh = [x[:, sel] for x in shocks]; ct = [x[:, sel] for x in ctrls]
        b_fe, _ = est.fe_core(d, sh, ct, **opts)
        b_spj, _ = est.spj_core(d, sh, ct, **opts)
        rows.append(dict(app=app, dropped=ids[drop],
                         beta_FE=float(b_fe[0]), beta_SPJ=float(b_spj[0])))
    base_fe, _ = est.fe_core(dep0, shocks, ctrls, **opts)
    base_spj, _ = est.spj_core(dep0, shocks, ctrls, **opts)
    d = pd.DataFrame(rows)
    d["dFE"] = (d.beta_FE - base_fe[0]).abs()
    d["dSPJ"] = (d.beta_SPJ - base_spj[0]).abs()
    d["base_FE"] = base_fe[0]; d["base_SPJ"] = base_spj[0]
    d.to_csv(f"{OUT}/{app}_loco.csv", index=False)
    return d


def se_alt(app):
    sp = cfg.APPS[app]
    df = pd.read_csv(f"{cfg.RES}/01_data_audit/{app}_clean.csv")
    if sp["y_scale"] != 1.0:
        for oc in sp["outcome_cols"]:
            if oc in df.columns:
                df[oc] = df[oc] * sp["y_scale"]
    rows = []
    for j, oc in enumerate(sp["outcome_cols"][:4]):
        dep, _ = wide(df, oc)
        shocks = [wide(df, c)[0] for c in sp["shock_cols"]]
        ctrls = [wide(df, c)[0] for c in sp["control_cols"]]
        for cl in ("entity", "twoway"):
            try:
                b, se = est.fe_core(dep, shocks, ctrls, te=sp["two_way_fe"],
                                    cluster=cl, robust=sp["robust_scale"],
                                    eigen=sp["eigen_fix"])
                rows.append(dict(app=app, h=sp["h_list"][j], outcome=oc, cluster=cl,
                                 beta=float(b[0]), se=float(se[0])))
            except Exception as e:
                rows.append(dict(app=app, h=sp["h_list"][j], outcome=oc, cluster=cl,
                                 beta=np.nan, se=np.nan))
    d = pd.DataFrame(rows)
    d.to_csv(f"{OUT}/{app}_se_alt.csv", index=False)
    return d


if __name__ == "__main__":
    for app in ("RR", "CS"):
        d = loco(app)
        print(app, "LOCO max |dFE|=%.4f max |dSPJ|=%.4f" %
              (d.dFE.max(), d.dSPJ.max()))
    for app in ("RR", "RR_UNEMP", "BVX", "MSV", "CS"):
        se_alt(app)
        print(app, "se_alt saved")
    print("S6 done ->", OUT)
