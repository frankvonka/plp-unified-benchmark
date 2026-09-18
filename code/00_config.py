# =============================================================================
# 00_config.py -- pre-registered plan, paths, estimators, grids, app specs
# =============================================================================
"""
plp-unified-benchmark
Unified assessment of panel local-projection (LP) estimators under
persistence, cross-sectional dependence (XSD) and heterogeneity.

Paper target (per user brief, 2026-09-17):
  "When Do Panel Local Projections Mislead? A Unified Assessment of
   Bias-Corrected Estimators under Persistence, Cross-Sectional Dependence,
   and Heterogeneous Dynamics."

Estimators (5):
  FE    : panel LP, entity FE (+time FE where app uses it), cluster-robust SE
  SPJ   : split-panel jackknife (Mei, Sheng & Shi 2026, JIE), exact port of
          their LP_panel_all.r (per-entity median complete-case split,
          beta = 2*full - 0.5*(halfA + halfB), DD-sandwich SE)
  HPJ   : higher-order jackknife -- constant-bias (SPJ) generalized to a
          second-order 1/N Taylor correction using three nested panel sizes
          (N, 2N/3, N/3); Jacobian (delta-method) SE. New contribution.
  DB    : their "data-based" analytical bias correction
          (IRF_FE + b0*f(rho,h)/((T0-h)*sx^2), formula from their main-text
          simulation simul_LP.R, tau=0 case)
  CDOLS : cumulative-difference pooled OLS (Ugarte-Ruiz 2026, BBVA WP 26/09,
          eq. 20): y_{t+h} - y_{t-1} on shock lags s_t, s_{t-1} (+ dy lags),
          pooled OLS, country-clustered SE

Metrics (per user brief sec. 8): bias, RMSE, MAE, coverage, FPR, FNR,
IRF shape distortion D_I = mean_h |beta_hat_h - beta_true_h|.

Seeds: 42, 7, 123 (user convention). All draws deterministic.
"""
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
P = os.path.abspath(ROOT)

RAW_PKG = f"{P}/raw/pkg"
APPS_CSV = f"{P}/raw/pkg/apps/applications"
SIMS_RDS = f"{P}/raw/pkg/sims/simulations"
GPR_XLSX = f"{P}/raw/GPRMonthly.xls"

DATA = f"{P}/data"
RES = f"{P}/results"
FIG = f"{P}/figures"
TAB = f"{P}/tables"
CKP = f"{P}/checkpoints"
MS = f"{P}/manuscript"

SEEDS = [42, 7, 123]
NITER_MAIN = 100        # iterations per seed per main-grid cell (perf-calibrated)
NITER_POWER = 150       # power / null blocks
NITER_VALID = 150       # their-toy validation cells (S8)
HORIZON = 10           # H: horizons 0..10 (their convention)

ESTIMATORS = ["FE", "SPJ", "HPJ", "DB", "CDOLS"]

# ---------------------------------------------------------------- grids
# main benchmark grid (96 cells; N x T x rho x XSD)
# spec-log note: T in {20,50,100} chosen over {20,30,50,100} for runtime on
# the 1-CPU/1.5GB box; 100 iters x 3 seeds per cell gives coverage SE ~2.9pp.
N_GRID = [20, 50, 100, 200]
T_GRID = [20, 50, 100]
RHO_GRID = [0.2, 0.5, 0.8, 0.95]
XSD_GRID = [0.0, 1.0]          # sigma_F (common factor in x); 0 = no XSD in x
# heterogeneous-persistence add-on (4 cells)
HET_RHO_CELLS = [(50, 30, 0.5, 0.0), (200, 100, 0.8, 1.0)]
# heterogeneous-beta add-on (4 cells)
HET_BETA_CELLS = [(50, 30, 0.5, 0.0), (100, 50, 0.8, 1.0)]
SIGMA_BETA = 0.3

# power / FPR-FNR block
POWER_CELLS = [(20, 20, 0.5, 0.0), (100, 100, 0.95, 1.0),
               (20, 50, 0.8, 1.0), (200, 50, 0.2, 0.0)]

# validation: their toy DGP exactly (N x T x rho), beta=-0.6, eta=0.2
VALID_CELLS = [(30, 60, 0.0), (30, 60, 0.2), (30, 60, 0.5), (30, 60, 0.8),
               (30, 120, 0.0), (30, 120, 0.2), (30, 120, 0.5), (30, 120, 0.8),
               (50, 60, 0.0), (50, 60, 0.2), (50, 60, 0.5), (50, 60, 0.8),
               (50, 120, 0.0), (50, 120, 0.2), (50, 120, 0.5), (50, 120, 0.8)]

# DGP constants (their toy): delta=0.2, eta=0.2, sigmas=1, beta0=-0.6
BETA0 = -0.6
DELTA = 0.2
ETA = 0.2
LAMBDA_SIG = 0.5        # sd of heterogeneous loadings lambda_i
GAMMA_SIG = 0.5         # sd of y factor loadings gamma_i (XSD in y)

# empirical application specs (ported 1:1 from their replication.ipynb)
APPS = {
    "RR": dict(
        csv="empirical_RR_f4_lngdp_1980.csv", label="Romer & Romer (2017), GDP",
        shock_cols=["CRISIS"],
        control_cols=["l1LNGDP", "l2LNGDP", "l3LNGDP", "l4LNGDP",
                      "l1CRISIS", "l2CRISIS", "l3CRISIS", "l4CRISIS"],
        outcome_cols=[f"f{h}LNGDP" for h in range(0, 11)],
        h_list=list(range(0, 11)),
        two_way_fe=True,       # te=T in notebook
        two_way_cluster=False, robust_scale=False, eigen_fix=False,
        y_scale=1.0, note="coefs x7 rescaled in their plot (annualizing half-year)"),
    "RR_UNEMP": dict(
        csv="empirical_RR_f4_lnunemp_1980.csv", label="Romer & Romer (2017), Unemployment",
        shock_cols=["CRISIS"],
        control_cols=["l1UNEMP", "l2UNEMP", "l3UNEMP", "l4UNEMP",
                      "l1CRISIS", "l2CRISIS", "l3CRISIS", "l4CRISIS"],
        outcome_cols=[f"f{h}UNEMP" for h in range(0, 11)],
        h_list=list(range(0, 11)), two_way_fe=True,
        two_way_cluster=False, robust_scale=False, eigen_fix=False, y_scale=1.0),
    "BVX": dict(
        csv="empirical_BVX_t1_y.csv", label="Baron, Verner & Xiong (2021), dY",
        shock_cols=["R_B", "R_N"],
        control_cols=["L1R_B", "L2R_B", "L3R_B", "L1R_N", "L2R_N", "L3R_N",
                      "D1y", "L1D1y", "L2D1y", "L3D1y",
                      "D1d_y", "L1D1d_y", "L2D1d_y", "L3D1d_y"],
        outcome_cols=[f"Fd{h}y" for h in range(1, 7)],
        h_list=list(range(1, 7)), two_way_fe=False,
        two_way_cluster=True, robust_scale=True, eigen_fix=True, y_scale=100.0),
    "MSV": dict(
        csv="empirical_MSV_f2.csv", label="Mian, Sufi & Verner (2017), household debt",
        shock_cols=["L0HHD_L1GDP", "L0NFD_L1GDP"],
        control_cols=["L0y", "L1y", "L2y", "L3y", "L4y",
                      "L1HHD_L1GDP", "L2HHD_L1GDP", "L3HHD_L1GDP", "L4HHD_L1GDP",
                      "L1NFD_L1GDP", "L2NFD_L1GDP", "L3NFD_L1GDP", "L4NFD_L1GDP"],
        outcome_cols=[f"F{h}y" for h in range(1, 11)],
        h_list=list(range(1, 11)), two_way_fe=False,
        two_way_cluster=True, robust_scale=False, eigen_fix=False, y_scale=1.0),
    "CS": dict(
        csv="empirical_CS_f3.csv", label="Cerra & Saxena (2008), crises",
        shock_cols=["CRISIS"],
        control_cols=["l1CRISIS", "l2CRISIS", "l3CRISIS", "l4CRISIS",
                      "l1GRRT_WB", "l2GRRT_WB", "l3GRRT_WB", "l4GRRT_WB"],
        outcome_cols=[f"cf{h}GRRT_WB" for h in range(1, 11)],
        h_list=list(range(1, 11)), two_way_fe=False,
        two_way_cluster=False, robust_scale=False, eigen_fix=False, y_scale=1.0),
}
# NOTE: RR unemp control column names verified at S1/S2 against the CSV.

# GPR application (new shock, user brief sec. 12-13)
GPR = dict(
    base="CS",                      # reuse CS panel: 192 countries x 1965-2000
    outcome_cols=[f"cf{h}GRRT_WB" for h in range(1, 4)],   # cf1..cf3
    shock_col="GPR",                # global GPR annual mean, joined by year
    window=(1985, 2000),            # GPR official series starts 1985
    exposure=("exposure", "pre-1985 log-GDP z-score (deterministic)"),
    moderator=("CRISIS_exp", "mean financial-crisis years 1965-1984 (deterministic)"),
    two_way_fe=False, two_way_cluster=False, robust_scale=False, eigen_fix=False,
    source="Caldara & Iacoviello (2022), Geopolitical Risk Index, official xlsx")

STAGE_TAGS = ["s0", "s1", "s2", "s3", "s4", "s5", "s6", "s7", "s8",
              "s9", "s10", "s11", "s12"]


def stage_paths(stage: str):
    """Checkpoint + result dirs for a stage."""
    ck = f"{CKP}/{stage}"
    rs = f"{RES}/{stage}"
    os.makedirs(ck, exist_ok=True)
    os.makedirs(rs, exist_ok=True)
    return ck, rs
