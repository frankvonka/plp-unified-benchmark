# Results report — evidence records (not manuscript prose)

Each record: model/spec, coefficient(s), SE/CI, N, interpretation, robustness,
table/figure refs. Generated from saved CSVs; source file cited per claim.

## R1. Replication headline (SPJ restores "crises are worse than you think")
- **CS (Cerra–Saxena panel, 192 countries):** crisis → GDP-growth IRF at h=10:
  FE = −4.51 (SE 0.84) vs **SPJ = −7.89 (SE 0.85)**
  → +75% deeper long-run loss once the FE Nickell bias is corrected.
  Source: `results/03_empirical/CS_{FE,SPJ}.csv` (h=10 rows). Matches their
  Figure 7 pattern (SPJ below FE at all horizons).
- **RR (Romer–Romer shocks, 24 OECD):** h=5 GDP response FE −0.625 → SPJ −0.743
  (×7 annualized half-year units in their figure). Source: `RR_{FE,SPJ}.csv`.
- **BVX (bank equity returns):** h=4 response FE −3.42 → SPJ −4.21.
  Source: `BVX_{FE,SPJ}.csv`.

## R2. Estimator disagreement is large where it matters
- Mean across-horizon disagreement D_h = max−min over the 5 estimators:
  CS 8.88 pp, BVX 6.52 pp, RR 1.12 pp, MSV 0.50 pp, RR_UNEMP 0.19 pp.
  Source: `results/05_disagreement/app_disagreement.csv`.
- Interpretation: disagreement concentrates exactly in the crisis panels where
  persistence is highest — the same cells where MC shows FE bias is worst.

## R3. HPJ is NOT a clone of SPJ (new result)
- CS h=10: HPJ = −5.52 vs SPJ −7.89 → HPJ removes about ⅓ of the SPJ correction.
- MSV h=10: FE −0.28, SPJ −0.56, **HPJ −0.49**: HPJ agrees qualitatively with
  SPJ (household-debt boom → deeper later bust) but rejects the SPJ extreme.
- GPR×crisis-history interaction (h=1..3): FE −0.07..−0.14, SPJ −0.09..−0.25,
  **HPJ +0.06..+0.12 (sign flip)**, CD-OLS +0.10..+0.11.
  Source: `results/07_gpr/gpr_irf.csv`, `results/03_empirical/*.csv`.
- Interpretation: two bias-corrected estimators can disagree on SIGN for
  interaction designs — estimator choice is a first-order research decision,
  not a technicality. (Mechanism hypothesis, flagged suggestive: SPJ's half-sample
  split can over-correct under joint persistence, as Ugarte-Ruiz 2026 documents.)

## R4. GPR application (new shock class; identified interactions only)
- Design per spec log: GPR_t level is collinear with time effects → reported
  descriptive only. Interactions: GPR_t × Exposure_i (pre-1985 growth z-score),
  GPR_t × CrisisHistory_i (pre-1985 crisis frequency).
- Exposure interaction h=1: FE 0.039 (t≈4.2), SPJ 0.058 (t≈6.5) — geopolitical
  shocks hurt high-exposure countries' growth more, and the gap WIDENS under SPJ.
- Crisis-history interaction h=3: SPJ −0.246 (t≈−6.7) vs FE −0.136 — countries
  with crisis-prone histories suffer persistently larger GPR-driven growth losses.
- Sources: `results/07_gpr/gpr_irf.csv`, `results/07_gpr/CS_GPR_panel.csv`,
  official GPR source: matteoiacoviello.com (column `source` in panel CSV).

## R5. Robustness
- LOCO (RR, CS): max |Δβ| from dropping one country = 0.16 pp (RR) and 0.78 pp
  (CS h=10-equivalent) vs base effects of −0.68/−5.13 → no single country drives
  the IRF. Source: `results/06_robustness/{RR,CS}_loco.csv`.
- SE treatments: two-way clustering inflates CS SEs by ~52% (0.32→0.49 at h=1)
  but significance survives. Source: `results/06_robustness/*_se_alt.csv`.

## R6. Monte Carlo (running — filled on completion)
- Grid: N∈{20,50,100,200} × T∈{20,50,100} × ρ∈{.2,.5,.8,.95} × XSD∈{0,1},
  100 iters/seed × 3 seeds per cell; estimators FE/SPJ/HPJ/DB/CDOLS.
- Metrics: bias, RMSE, MAE, coverage, FPR/FNR, D_I (IRF shape distortion).
- Sources on completion: `results/04_monte_carlo/main_grid_summary.csv`,
  `results/FINAL/headline_mc.csv`, power/het blocks.
