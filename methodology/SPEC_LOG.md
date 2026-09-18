# Specification Log (methodology evidence package)

Pre-registered decisions and post-hoc notes. Any change after seeing results
is versioned here with reason.

## Estimator implementations
| Spec | Purpose | Reason | Result | Status |
|---|---|---|---|---|
| FE: entity demeaning + within-horizon joint complete-case mask, cluster sandwich | replicate their LP_panel_all.r FE | their R demeaning is horizon-specific under missingness | matches linearmodels clustered entity SE; coef exact to 1e-9 | FINAL |
| SPJ: cut = floor(median of complete-case rows) per entity, per-horizon; beta = 2b_full − ½(bA+bB); DD-sandwich SE, Q from original full design | exact port | their R builds `dd.mat` but keeps `Q.hat` from `indep_var` (full design) | matches validated path to 1e-16 | FINAL |
| HPJ: 4.5·b(N) − 4·b(2N/3) + 0.5·b(N/3), nested entity prefixes, round() convention, SE = full-sample FE sandwich | kill 1/N + 1/N² bias terms (two-stage jackknife order 2) | new contribution; documented in code header | consistent with MC (lower bias at short T) | FINAL |
| DB: their τ=0 analytic correction, pooled AR(1) ρ̂ + b0 + sx² from raw y & shock; applied per horizon h | their main-text "data-based" correction | their formula is single-x; for multi-shock apps the primary shock column is used (documented limitation) | matches their `simul_LP.R` tau=0 branch | FINAL |
| CD-OLS: LHS = y_{t+h} − y_{t−1}, RHS = [s_t, s_{t−1}] pooled OLS, entity-clustered SE | Ugarte-Ruiz eq. 20 | their equation needs raw y_{t−1}; apps without h0 use out_h − out_{h−1} (documented; MSV has F1..F10 only, so LHS_1 NaN by construction) | unbiased on their DGP; see toy_-0.6 comparison | FINAL |

## Grid / conventions
| Spec | Purpose | Reason | Result | Status |
|---|---|---|---|---|
| T ∈ {20,50,100} (not 30) | runtime on 1-CPU/1.5GB box | NITER=100/seed × 96 cells ≈ 8h | accepted; noted as limitation | FINAL |
| te=F, robust=F, entity cluster in S4 (their main-sim convention) | comparability with their toy RDS | their `main_simul_*` scripts used LP_panel defaults | matches toy_-0.6 coverage within MC noise | FINAL |
| S7 empirical uses robust=True + app-specific cluster | their notebook per-app flags | faithful replication of their 4 figures | REPRODUCED (SPJ>FE magnitude in crisis apps) | FINAL |
| GPR app: level spec reported as descriptive (collinear with time FE); identified specs = GPR×Exposure, GPR×CrisisHistory | user brief §13 + panel-econometrics rule | time-varying-only treatment unidentified under TWFE | done; GPR official source column attached | FINAL |
| Exposure = pre-1985 mean cf1GRRT_WB z-score; CrisisHistory = pre-1985 mean CRISIS | predetermined, data-internal (no external merge beyond GPR) | CS file lacks GDP levels / financial openness → documented proxy, flagged as limitation | FINAL |
| seeds 42/7/123, 100/150/200 iters | user convention; coverage SE ≈ 2.6–4pp at 100×3 | RAM/time budget | FINAL |
| FPR/FNR block: null (β=0) & power (β=−0.6) at 4 design points, 150 iters/seed | user brief §8 | — | FINAL |
| S10 LP: simplex-vertex equivalence check | variation rule (2 families: LP+surrogates) | LP over simplex → vertex; high solve used as numerical confirmation | agreement logged | FINAL |
| S12 figure data: all values from saved CSVs, none hand-typed | user hard rule | — | FINAL |

## Validation
- S8 compares our MC summaries (600 iters, our seeds) vs their toy_-0.6.xlsx
  (1000 iters, their seed 2023) on FE/SPJ/DB × 16 cells × 11 horizons.
- Pass criterion: ≥90% of cells within 0.06 coverage / 10% RMSE relative.

## Known limitations (reported, not hidden)
1. MSV/CS/BVX CD-OLS LHS for h=1 uses out_1 − out_0 where out_0 missing → NaN
   for those apps (documented in S3 spec log; figure excludes NaNs).
2. BBVA GPR monthly official file is a single vintage; no revision history.
3. Proxy exposures for GPR interactions (growth-based, not financial-openness).
4. S4 main grid omits T=30 by runtime constraint (noted above).
5. **HPJ naive-SE over-rejection (S4 power block):** with the SE carried over
   from the full-sample FE sandwich, HPJ's FPR = 0.50 vs 0.06–0.08 for the
   other estimators (null cells). The jackknife combination adds O(1/N)
   variance the sandwich does not capture; a correct HPJ SE needs the delta
   method (Jacobian of the 4.5/−4/0.5 combination) or a paired-block jackknife.
   Flagged as a known limitation of the naive implementation; point estimates
   remain the contribution (bias reduction), inference caution documented.
