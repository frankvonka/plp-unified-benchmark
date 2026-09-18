# =============================================================================
# 02_empirical.py -- S2 descriptives + S3: 5 apps x 5 estimators
# =============================================================================
"""
Runs the exact notebook specs (ported 1:1 from their replication.ipynb) for
FE / SPJ (validated vs linearmodels and their R formulas), plus our new HPJ,
DB and CD-OLS estimators.
CD-OLS (Ugarte-Ruiz 2026, eq. 20) LHS is y_{t+h} - y_{t-1}; with the
distributed CSVs only shifted outcome columns exist, so per spec log:
  - apps WITH an h0 column (RR, RR_UNEMP): LHS_h = f{h} - f0
  - apps WITHOUT h0 (BVX: Fd{h}y; MSV: F{h}y; CS: cf{h}): LHS_h = out{h} - out{h-1}
    for h>=2; LHS_1 = out_1 - base, base = D1y (BVX control), NaN otherwise.
Outputs: results/03_empirical/{app}_{est}.csv, checkpoints/s3/{app}.json
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

OUT = f"{cfg.RES}/03_empirical"; os.makedirs(OUT, exist_ok=True)
os.makedirs(f"{cfg.CKP}/s3", exist_ok=True)
os.makedirs(f"{cfg.RES}/02_descriptives", exist_ok=True)


def descriptives():
    rows = []
    for tag, sp in cfg.APPS.items():
        df = pd.read_csv(f"{cfg.RES}/01_data_audit/{tag}_clean.csv")
        for v in sp["shock_cols"] + [sp["outcome_cols"][0]]:
            if v in df.columns:
                s = df[v].dropna()
                rows.append(dict(app=tag, var=v, n=len(s), mean=s.mean(),
                                 sd=s.std(), min=s.min(), p50=s.median(), max=s.max()))
    pd.DataFrame(rows).to_csv(f"{cfg.RES}/02_descriptives/descriptives.csv", index=False)
    print("S2 descriptives saved")


def build_wide(df, sp, col):
    """Long -> wide T x N for one column (their NaN padding)."""
    idc, tc = df.columns[0], df.columns[1]
    times = sorted(df[tc].unique()); ids = sorted(df[idc].unique())
    T, N = len(times), len(ids)
    tix = {t: j for j, t in enumerate(times)}; iix = {i: k for k, i in enumerate(ids)}
    m = np.full((T, N), np.nan)
    v = df[col].to_numpy(float)
    r = df[tc].map(tix).to_numpy(); c = df[idc].map(iix).to_numpy()
    ok = np.isfinite(v)
    m[r[ok].astype(int), c[ok].astype(int)] = v[ok]
    return m


def run_app(tag, sp, df):
    idc, tc = df.columns[0], df.columns[1]
    opts = dict(te=sp["two_way_fe"], cluster=("twoway" if sp["two_way_cluster"] else "entity"),
                robust=sp["robust_scale"], eigen=sp["eigen_fix"])
    h_list = sp["h_list"]
    has_h0 = h_list[0] == 0
    base_col = sp["outcome_cols"][0] if has_h0 else ("D1y" if tag == "BVX" else None)
    W = {oc: build_wide(df, sp, oc) for oc in sp["outcome_cols"]}
    W.update({c: build_wide(df, sp, c) for c in sp["shock_cols"] + sp["control_cols"]
              if c not in W})
    if base_col and base_col not in W:
        W[base_col] = build_wide(df, sp, base_col)
    shocks = [W[c] for c in sp["shock_cols"]]
    ctrls = [W[c] for c in sp["control_cols"]]
    res_all = {}
    for j, oc in enumerate(sp["outcome_cols"]):
        h = h_list[j]
        dep = W[oc]
        # CD-OLS LHS base per spec log
        if has_h0:
            base = W[base_col]                       # f0
        elif j >= 1:
            base = W[sp["outcome_cols"][j - 1]]      # out{h-1}
        elif base_col:
            base = W[base_col]                        # D1y for BVX
        else:
            base = np.full_like(dep, np.nan)
        row = {}
        try:
            b, se = est.fe_core(dep, shocks, ctrls, **opts)
            row["FE"] = (b[0], se[0])
            b, se = est.spj_core(dep, shocks, ctrls, **opts)
            row["SPJ"] = (b[0], se[0])
            b, se = est.hpj_core(dep, shocks, ctrls, **opts)
            row["HPJ"] = (b[0], se[0])
            bfe, sefe = row["FE"]
            b, se = est.db_correction(dep, shocks[0], bfe, sefe)
            row["DB"] = (b[0], se[0])
            b, se = est.cdols_ols(dep - base, shocks[0],
                                  shocks[1] if len(shocks) > 1 else None)
            row["CDOLS"] = (b, se)
        except Exception as e:
            print(f"  {tag}/{oc} FAILED: {type(e).__name__}: {e}")
            row = {k: (np.nan, np.nan) for k in cfg.ESTIMATORS}
        res_all[oc] = dict(h=h, res=row)
    return res_all


def main():
    descriptives()
    for tag, sp in cfg.APPS.items():
        print(f"=== {tag}: {sp['label']}")
        df = pd.read_csv(f"{cfg.RES}/01_data_audit/{tag}_clean.csv")
        if sp["y_scale"] != 1.0:
            for oc in sp["outcome_cols"]:
                if oc in df.columns:
                    df[oc] = df[oc] * sp["y_scale"]
        res_all = run_app(tag, sp, df)
        for estn in cfg.ESTIMATORS:
            rows = []
            for oc, d in res_all.items():
                b, se = d["res"].get(estn, (np.nan, np.nan))
                rows.append(dict(h=d["h"], outcome=oc, beta=b, se=se,
                                 lo=(b - 1.96 * se if np.isfinite(se) else np.nan),
                                 hi=(b + 1.96 * se if np.isfinite(se) else np.nan)))
            pd.DataFrame(rows).sort_values("h").to_csv(f"{OUT}/{tag}_{estn}.csv", index=False)
        ck = {oc: {k: [None if not np.isfinite(v[0]) else float(v[0]),
                       None if not np.isfinite(v[1]) else float(v[1])]
                   for k, v in d["res"].items()} for oc, d in res_all.items()}
        json.dump(ck, open(f"{cfg.CKP}/s3/{tag}.json", "w"), indent=1)
        print(f"  saved {tag}: {len(res_all)} horizons")
    print("S3 core done")


if __name__ == "__main__":
    main()
