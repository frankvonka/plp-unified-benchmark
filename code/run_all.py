# =============================================================================
# run_all.py -- S0 pipeline runner (stage order, per-stage push+tag ready)
# =============================================================================
"""
Usage:
  python run_all.py            # run every stage in order
  python run_all.py s4 power    # single stage + its sub-argument
"""
import subprocess, sys, time
PY = "/root/.venvs/shared/bin/python"
HERE = sys.path[0] if sys.path[0].startswith("/") else \
    "/root/projects/plp-unified-benchmark/code"

STAGES = [
    ("s1",  ["code/01_s1_data_audit.py"]),
    ("s2s3", ["code/02_empirical.py"]),
    ("s4main", ["code/03_s4_monte_carlo.py", "main"]),
    ("s4power", ["code/03_s4_monte_carlo.py", "power"]),
    ("s4het", ["code/03_s4_monte_carlo.py", "het"]),
    ("s5",  ["code/05_s5_disagreement.py"]),
    ("s6",  ["code/11_s6_robustness.py"]),
    ("s7",  ["code/04_s7_gpr.py"]),
    ("s8",  ["code/06_s8_validation.py"]),
    ("s9",  ["code/07_s9_xsd_inference.py"]),
    ("s10", ["code/08_s10_optimization.py"]),
    ("s11", ["code/09_s11_final.py"]),
    ("s12", ["code/10_s12_figures.py"]),
]


def run(cmd):
    print("$ " + " ".join(cmd), flush=True)
    r = subprocess.run(cmd, cwd="/root/projects/plp-unified-benchmark")
    if r.returncode != 0:
        print(f"STAGE FAILED rc={r.returncode}: {' '.join(cmd)}", flush=True)
        sys.exit(r.returncode)


def main():
    only = sys.argv[1] if len(sys.argv) > 1 else None
    for tag, cmd in STAGES:
        if only and only not in tag:
            continue
        t0 = time.time()
        run([PY] + cmd)
        print(f"[{tag}] ok in {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
