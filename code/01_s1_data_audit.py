# =============================================================================
# 01_s1_data_audit.py -- S1: data audit & cleaning of the 5 application CSVs
# =============================================================================
import sys, os, json
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import importlib.util
spec = importlib.util.spec_from_file_location("cfg", os.path.join(os.path.dirname(os.path.abspath(__file__)), "00_config.py"))
cfg = importlib.util.module_from_spec(spec); spec.loader.exec_module(cfg)

OUT = f"{cfg.RES}/01_data_audit"
os.makedirs(OUT, exist_ok=True)
audit = {}
for tag, spec_ in cfg.APPS.items():
    p = f"{cfg.APPS_CSV}/{spec_['csv']}"
    df = pd.read_csv(p)
    idc, tc = df.columns[0], df.columns[1]
    tmin, tmax = int(df[tc].min()), int(df[tc].max())
    per = df.groupby(idc)[tc].count()
    rec = {
        "file": spec_["csv"], "label": spec_["label"],
        "n_rows": len(df), "n_units": int(df[idc].nunique()),
        "n_periods_raw": int(df[tc].nunique()), "period_range": [tmin, tmax],
        "n_missing_cells": int(df.isna().sum().sum()),
        "balanced": bool(per.nunique() == 1),
        "median_T": int(per.median()), "min_T": int(per.min()), "max_T": int(per.max()),
        "shock_cols": spec_["shock_cols"], "n_controls": len(spec_["control_cols"]),
        "n_outcomes": len(spec_["outcome_cols"]), "te": spec_["two_way_fe"],
        "twoc": spec_["two_way_cluster"], "robust": spec_["robust_scale"],
        "eigen": spec_["eigen_fix"], "y_scale": spec_["y_scale"],
    }
    # verify required columns exist
    needed = spec_["shock_cols"] + spec_["control_cols"] + spec_["outcome_cols"]
    missing = [c for c in needed if c not in df.columns]
    rec["missing_columns"] = missing
    assert not missing, f"{tag}: missing {missing}"
    # duplicate key check
    rec["dup_keys"] = int(df.duplicated([idc, tc]).sum())
    # plausibility of shock share
    s0 = spec_["shock_cols"][0]
    rec["shock_share"] = round(float(df[s0].fillna(0).mean()), 4)
    # cleaned CSV: sorted, no key duplicates, numeric coerced
    d2 = df.drop_duplicates([idc, tc]).sort_values([idc, tc]).reset_index(drop=True)
    num = d2.select_dtypes("number")
    num = num.apply(pd.to_numeric, errors="coerce")
    d2 = pd.concat([d2.drop(columns=num.columns), num], axis=1)
    d2.to_csv(f"{OUT}/{tag}_clean.csv", index=False)
    audit[tag] = rec
    print(tag, {k: rec[k] for k in ["n_rows", "n_units", "period_range",
          "n_missing_cells", "balanced", "dup_keys", "shock_share"]})
json.dump(audit, open(f"{OUT}/audit.json", "w"), indent=1)
# checkpoint
os.makedirs(f"{cfg.CKP}/s1", exist_ok=True)
pd.DataFrame(audit, index=audit.keys()).T.to_csv(f"{cfg.CKP}/s1/data_audit_summary.csv")
print("S1 done ->", OUT)
