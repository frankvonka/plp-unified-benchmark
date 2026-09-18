#!/usr/bin/env python
"""
Fetch external data so the pipeline (run_all.py) runs from a clean checkout.
Downloads the two official sources used by the project:
  1) Mendeley p8sxs7sxjr (Mei-Sheng-Sheng Nickell-bias replication pkg) -> raw/pkg/
  2) Caldara-Iacoviello GPR official monthly file -> raw/GPRMonthly.xls
All >100MB-safe (both are <10MB). Idempotent: skips if present.
"""
import os, sys, urllib.request, zipfile

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # project root
RAW = f"{HERE}/raw"


def fetch_mendeley():
    out = f"{RAW}/mendeley_p8sxs7sxjr.zip"
    if not os.path.exists(out):
        url = "https://data.mendeley.com/public-api/zip/p8sxs7sxjr/download/1"
        print("downloading Mendeley p8sxs7sxjr ...", flush=True)
        urllib.request.urlretrieve(url, out)
    pkg = f"{RAW}/pkg"
    if not os.path.isdir(pkg):
        with zipfile.ZipFile(out) as z:
            z.extractall(f"{RAW}/pkg_tmp")
        os.rename(f"{RAW}/pkg_tmp", pkg)
    print("Mendeley pkg -> raw/pkg/", flush=True)


def fetch_gpr():
    out = f"{RAW}/GPRMonthly.xls"
    if not os.path.exists(out):
        url = "https://www.matteoiacoviello.com/gpr_files/data_gpr_export.xls"
        print("downloading official GPR monthly file ...", flush=True)
        urllib.request.urlretrieve(url, out)
    print("GPR -> raw/GPRMonthly.xls", flush=True)


if __name__ == "__main__":
    os.makedirs(RAW, exist_ok=True)
    fetch_mendeley()
    fetch_gpr()
    print("done")
