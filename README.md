# plp-unified-benchmark

**When Do Panel Local Projections Mislead? A Unified Assessment of Bias-Corrected
Estimators under Persistence, Cross-Sectional Dependence, and Heterogeneous Dynamics.**

Reproducible research pipeline (Python, no R) benchmarking five panel local-projection
estimators on the Mei–Sheng–Sheng (2026, *JIE*) Nickell-bias replication package:

| Estimator | Description | Source |
|---|---|---|
| **FE** | panel LP with entity FE, cluster-robust SE | standard |
| **SPJ** | split-panel jackknife `2b_full − ½(b_A+b_B)`, DD-sandwich SE | exact port of their `LP_panel_all.r` |
| **HPJ** | second-order jackknife `4.5·b(N) − 4·b(2N/3) + 0.5·b(N/3)` | **new here** |
| **DB** | analytical correction `FE_h + b0·f(ρ,h)/((T−h)·s_x²)` | their `simul_LP.R` (τ=0 branch) |
| **CD-OLS** | cumulative-difference pooled OLS `y_{t+h} − y_{t−1} ~ s_t, s_{t−1}` | Ugarte-Ruiz (2026, BBVA WP 26/09, eq. 20) |

## Contents
- **S1** data audit of the 5 application panels (RR, RR-unemp, BVX, MSV, CS)
- **S2–S3** exact replication + extension: all 5 apps × 5 estimators
- **S4** Monte Carlo benchmark: 96-cell grid (N×T×ρ×XSD) + power/null + heterogeneity blocks
- **S5** estimator-disagreement (D_h) + reliability grade map
- **S6** LOCO + alternative-clustering robustness
- **S7** new application: geopolitical-risk (GPR) shocks with exposure interactions
- **S8** validation of our MC pipeline against their published toy (toy_-0.6.xlsx)
- **S9** inference under cross-sectional dependence (entity vs time vs two-way SE)
- **S10** estimator selection as optimization (LP + ridge/kNN surrogates)
- **S11** frozen FINAL results + `ALL_RESULTS.md`
- **S12** 12-figure renderer + self-contained HTML manuscript

## Reproduce
```bash
/root/.venvs/shared/bin/python code/run_all.py     # runs S1..S12 in order
```
Seeds: 42, 7, 123 (all stages deterministic; CV folds seeded 42).

Data: Mendeley `p8sxs7sxjr/1` (Mei–Sheng–Sheng replication package) in `raw/pkg/`;
official GPR index (Caldara–Iacoviello, matteoiacoviello.com) in `raw/GPRMonthly.xls`.
