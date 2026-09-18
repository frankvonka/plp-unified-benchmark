# Tables

All values generated from saved result files (seeds 42/7/123); no manually typed numbers.

### T1. Data audit: application panels (Mei-Sheng-Sheng replication data)

| app      |
|:---------|
| RR       |
| RR_UNEMP |
| BVX      |
| MSV      |
| CS       |

### T2. Descriptive statistics (shock and first-outcome variables)

| app      | var         |    n |     mean |      sd |      min |      p50 |      max |
|:---------|:------------|-----:|---------:|--------:|---------:|---------:|---------:|
| RR       | CRISIS      | 1584 |    0.715 |   1.954 |    0     |    0     |   14     |
| RR       | f0LNGDP     | 1728 | 1296.62  | 156.74  |  861.284 | 1280.02  | 1662.92  |
| RR_UNEMP | CRISIS      | 1584 |    0.715 |   1.954 |    0     |    0     |   14     |
| RR_UNEMP | f0UNEMP     | 1518 |    7.43  |   4.049 |    0.207 |    6.967 |   27.7   |
| BVX      | R_B         | 3317 |    0.067 |   0.25  |    0     |    0     |    1     |
| BVX      | R_N         | 3317 |    0.071 |   0.257 |    0     |    0     |    1     |
| BVX      | Fd1y        | 3317 |    0.033 |   0.043 |   -0.265 |    0.034 |    0.256 |
| MSV      | L0HHD_L1GDP |  910 |   50.168 |  28.69  |    0.201 |   48.717 |  147.933 |
| MSV      | L0NFD_L1GDP |  900 |   82.349 |  33.449 |    8.113 |   79.275 |  234.623 |
| MSV      | F1y         |  881 | 2824.2   | 238.954 | 2453.44  | 2781.33  | 3550.2   |
| CS       | CRISIS      | 5398 |    0.258 |   0.438 |    0     |    0     |    1     |
| CS       | cf1GRRT_WB  | 5812 |    7.49  |  11.387 |  -81.3   |    7.756 |  151.876 |

### T3. Crisis -> GDP growth IRFs, Cerra-Saxena panel (beta, cluster SE)

|   h | FE             | SPJ            | HPJ            | DB             | CDOLS         |
|----:|:---------------|:---------------|:---------------|:---------------|:--------------|
|   1 | -2.589 (0.324) | -2.619 (0.338) | -3.737 (0.324) | -2.589 (0.324) | —             |
|   2 | -3.385 (0.451) | -3.488 (0.474) | -3.423 (0.451) | -3.385 (0.451) | 3.073 (0.252) |
|   3 | -4.202 (0.546) | -4.788 (0.606) | -4.094 (0.546) | -4.202 (0.546) | 3.153 (0.243) |
|   4 | -4.543 (0.637) | -5.485 (0.714) | -5.923 (0.637) | -4.543 (0.637) | 3.316 (0.240) |
|   5 | -4.581 (0.695) | -6.088 (0.778) | -4.968 (0.695) | -4.581 (0.695) | 3.372 (0.220) |
|   6 | -5.134 (0.725) | -7.092 (0.826) | -4.353 (0.725) | -5.134 (0.725) | 3.353 (0.234) |
|   7 | -4.720 (0.734) | -7.396 (0.834) | -4.459 (0.734) | -4.720 (0.734) | 3.490 (0.231) |
|   8 | -4.613 (0.755) | -7.747 (0.876) | -3.877 (0.755) | -4.613 (0.755) | 3.365 (0.242) |
|   9 | -4.311 (0.781) | -7.613 (0.940) | -4.774 (0.781) | -4.311 (0.781) | 3.301 (0.205) |
|  10 | -4.507 (0.840) | -7.893 (1.034) | -5.519 (0.840) | -4.507 (0.840) | 3.215 (0.233) |

### T4. Tax-shock -> GDP IRFs, Romer-Romer panel (beta, cluster SE)

|   h | FE             | SPJ            | HPJ            | DB             | CDOLS          |
|----:|:---------------|:---------------|:---------------|:---------------|:---------------|
|   0 | -0.293 (0.079) | -0.310 (0.084) | -0.429 (0.079) | -0.293 (0.079) | -0.000 (0.000) |
|   1 | -0.392 (0.127) | -0.420 (0.137) | -0.567 (0.127) | -0.392 (0.127) | -0.021 (0.041) |
|   2 | -0.439 (0.143) | -0.470 (0.167) | -0.248 (0.143) | -0.439 (0.143) | 0.038 (0.086)  |
|   3 | -0.568 (0.169) | -0.647 (0.213) | -0.073 (0.169) | -0.568 (0.169) | 0.134 (0.128)  |
|   4 | -0.620 (0.216) | -0.738 (0.263) | 0.076 (0.216)  | -0.620 (0.216) | 0.251 (0.166)  |
|   5 | -0.625 (0.223) | -0.743 (0.280) | 0.239 (0.223)  | -0.625 (0.223) | 0.367 (0.200)  |
|   6 | -0.678 (0.226) | -0.790 (0.291) | 0.440 (0.226)  | -0.678 (0.226) | 0.479 (0.231)  |
|   7 | -0.776 (0.238) | -0.898 (0.298) | 0.553 (0.238)  | -0.776 (0.238) | 0.588 (0.259)  |
|   8 | -0.708 (0.258) | -0.864 (0.323) | 0.872 (0.258)  | -0.708 (0.258) | 0.740 (0.288)  |
|   9 | -0.601 (0.275) | -0.722 (0.347) | 0.836 (0.275)  | -0.601 (0.275) | 0.903 (0.311)  |
|  10 | -0.559 (0.286) | -0.715 (0.360) | 1.094 (0.286)  | -0.559 (0.286) | 1.086 (0.323)  |

### T5. Estimator disagreement: D_h = max-min across 5 estimators

| app      |   mean_Dh |   max_Dh |
|:---------|----------:|---------:|
| BVX      |     6.52  |    7.893 |
| CS       |     8.881 |   11.112 |
| MSV      |     0.495 |    0.707 |
| RR       |     1.117 |    1.809 |
| RR_UNEMP |     0.193 |    0.319 |

_Note: D_h computed per horizon then averaged (mean) / maxed._

### T6. Monte Carlo main grid: horizon-mean metrics by T, XSD, estimator

|   T |   sig_lam | est   |    bias |   rmse |    mae |   coverage |
|----:|----------:|:------|--------:|-------:|-------:|-----------:|
|  20 |         0 | CDOLS | -0.0368 | 0.0681 | 0.0558 |     0.8352 |
|  20 |         0 | DB    |  0.0569 | 0.0767 | 0.0692 |     0.6154 |
|  20 |         0 | FE    |  0.115  | 0.1265 | 0.1199 |     0.3737 |
|  20 |         0 | HPJ   |  0.1134 | 0.1654 | 0.1437 |     0.2931 |
|  20 |         0 | SPJ   |  0.0543 | 0.0824 | 0.072  |     0.7164 |
|  20 |         1 | CDOLS | -0.0156 | 0.0922 | 0.0724 |     0.6452 |
|  20 |         1 | DB    |  0.0713 | 0.124  | 0.0984 |     0.5314 |
|  20 |         1 | FE    |  0.1261 | 0.1541 | 0.1384 |     0.4037 |
|  20 |         1 | HPJ   |  0.1252 | 0.1845 | 0.1577 |     0.3174 |
|  20 |         1 | SPJ   |  0.0688 | 0.135  | 0.1114 |     0.5877 |
|  50 |         0 | CDOLS | -0.0423 | 0.0562 | 0.0488 |     0.6471 |
|  50 |         0 | DB    |  0.0217 | 0.036  | 0.0311 |     0.702  |
|  50 |         0 | FE    |  0.0614 | 0.0688 | 0.0647 |     0.4217 |
|  50 |         0 | HPJ   |  0.0621 | 0.0958 | 0.0822 |     0.311  |
|  50 |         0 | SPJ   |  0.0121 | 0.0334 | 0.0273 |     0.8369 |
|  50 |         1 | CDOLS | -0.0227 | 0.0678 | 0.0532 |     0.5781 |
|  50 |         1 | DB    |  0.0281 | 0.068  | 0.0562 |     0.6395 |
|  50 |         1 | FE    |  0.0682 | 0.0931 | 0.0813 |     0.5153 |
|  50 |         1 | HPJ   |  0.0688 | 0.1107 | 0.0932 |     0.4375 |
|  50 |         1 | SPJ   |  0.0199 | 0.0789 | 0.0633 |     0.6779 |
| 100 |         0 | CDOLS | -0.0395 | 0.0475 | 0.0428 |     0.5034 |
| 100 |         0 | DB    |  0.0069 | 0.0185 | 0.0151 |     0.7955 |
| 100 |         0 | FE    |  0.0296 | 0.0354 | 0.0324 |     0.5146 |
| 100 |         0 | HPJ   |  0.0293 | 0.0555 | 0.0463 |     0.3539 |
| 100 |         0 | SPJ   |  0.0007 | 0.0182 | 0.0146 |     0.8799 |
| 100 |         1 | CDOLS | -0.0221 | 0.0498 | 0.0396 |     0.5502 |
| 100 |         1 | DB    |  0.0102 | 0.0428 | 0.0346 |     0.7466 |
| 100 |         1 | FE    |  0.0331 | 0.0543 | 0.0456 |     0.6457 |
| 100 |         1 | HPJ   |  0.0313 | 0.0691 | 0.0558 |     0.5558 |
| 100 |         1 | SPJ   |  0.0031 | 0.049  | 0.0387 |     0.7594 |

### T7. Coverage by shock persistence (T=50, no XSD; nominal 0.95)

|   rho |   CDOLS |    DB |    FE |   HPJ |   SPJ |
|------:|--------:|------:|------:|------:|------:|
|  0.2  |   0.807 | 0.935 | 0.821 | 0.478 | 0.931 |
|  0.5  |   0.7   | 0.913 | 0.546 | 0.414 | 0.895 |
|  0.8  |   0.543 | 0.668 | 0.205 | 0.221 | 0.782 |
|  0.95 |   0.539 | 0.292 | 0.115 | 0.131 | 0.739 |

_Note: Horizon-averaged coverage per cell; h = 0..min(10, T/3)._

### T8. Size and power: rejection rates when H0 true / H1 true

| estimator   |   FPR (size) |    FNR |
|:------------|-------------:|-------:|
| CDOLS       |       0.0639 | 0.2489 |
| DB          |       0.0808 | 0.3864 |
| FE          |       0.0781 | 0.3726 |
| HPJ         |       0.5001 | 0.2586 |
| SPJ         |       0.0573 | 0.3872 |

_Note: Null cells: beta=0; power cells: beta=-0.6. 150 iters x 3 seeds per design point._

### T9. Coverage under cross-sectional dependence by SE treatment

| est   |   ent |    time |      tw |
|:------|------:|--------:|--------:|
| CDOLS | 0.591 | nan     | nan     |
| FE    | 0.457 |   0.28  |   0.468 |
| HPJ   | 0.42  |   0.267 |   0.441 |
| SPJ   | 0.655 |   0.462 |   0.685 |

_Note: ent = entity cluster; time = DK lag-0 kernel; tw = two-way cluster (Cameron-Miller)._

### T10. GPR shocks x predetermined exposures (identified interactions)

| shock      |   h | FE             | SPJ            | HPJ           | DB             | CDOLS         |
|:-----------|----:|:---------------|:---------------|:--------------|:---------------|:--------------|
| GPR_x_cris |   1 | -0.074 (0.024) | -0.087 (0.018) | 0.059 (0.024) | -0.074 (0.024) | 0.106 (0.009) |
| GPR_x_cris |   2 | -0.102 (0.035) | -0.164 (0.030) | 0.115 (0.035) | -0.114 (0.035) | 0.102 (0.009) |
| GPR_x_cris |   3 | -0.136 (0.042) | -0.246 (0.037) | 0.094 (0.042) | -0.155 (0.042) | 0.100 (0.009) |
| GPR_x_exp  |   1 | 0.030 (0.009)  | 0.030 (0.006)  | 0.021 (0.009) | 0.030 (0.009)  | 0.004 (0.004) |
| GPR_x_exp  |   2 | 0.044 (0.012)  | 0.058 (0.009)  | 0.017 (0.012) | 0.049 (0.012)  | 0.005 (0.004) |
| GPR_x_exp  |   3 | 0.055 (0.014)  | 0.089 (0.012)  | 0.022 (0.014) | 0.062 (0.014)  | 0.005 (0.004) |

_Note: Level GPR effect is collinear with time effects and not reported as identified._

### T11. Reliability grades: good / moderate / poor (96 designs)

|   N |   T |   rho |   sig_lam | CDOLS    | DB       | FE       | HPJ   | SPJ      |
|----:|----:|------:|----------:|:---------|:---------|:---------|:------|:---------|
|  20 |  20 |  0.2  |         0 | moderate | moderate | moderate | poor  | moderate |
|  20 |  20 |  0.2  |         1 | poor     | moderate | poor     | poor  | moderate |
|  20 |  20 |  0.5  |         0 | moderate | moderate | poor     | poor  | moderate |
|  20 |  20 |  0.5  |         1 | poor     | poor     | poor     | poor  | poor     |
|  20 |  20 |  0.8  |         0 | moderate | poor     | poor     | poor  | poor     |
|  20 |  20 |  0.8  |         1 | poor     | poor     | poor     | poor  | poor     |
|  20 |  20 |  0.95 |         0 | moderate | poor     | poor     | poor  | poor     |
|  20 |  20 |  0.95 |         1 | poor     | poor     | poor     | poor  | poor     |
|  20 |  50 |  0.2  |         0 | moderate | good     | moderate | poor  | good     |
|  20 |  50 |  0.2  |         1 | poor     | moderate | moderate | poor  | moderate |
|  20 |  50 |  0.5  |         0 | moderate | good     | poor     | poor  | good     |
|  20 |  50 |  0.5  |         1 | poor     | moderate | poor     | poor  | poor     |
|  20 |  50 |  0.8  |         0 | poor     | poor     | poor     | poor  | moderate |
|  20 |  50 |  0.8  |         1 | poor     | poor     | poor     | poor  | poor     |
|  20 |  50 |  0.95 |         0 | poor     | poor     | poor     | poor  | poor     |
|  20 |  50 |  0.95 |         1 | poor     | poor     | poor     | poor  | poor     |
|  20 | 100 |  0.2  |         0 | moderate | good     | good     | poor  | good     |
|  20 | 100 |  0.2  |         1 | poor     | moderate | moderate | poor  | moderate |
|  20 | 100 |  0.5  |         0 | poor     | good     | moderate | poor  | good     |
|  20 | 100 |  0.5  |         1 | poor     | moderate | poor     | poor  | moderate |
|  20 | 100 |  0.8  |         0 | poor     | moderate | poor     | poor  | moderate |
|  20 | 100 |  0.8  |         1 | poor     | poor     | poor     | poor  | poor     |
|  20 | 100 |  0.95 |         0 | poor     | poor     | poor     | poor  | moderate |
|  20 | 100 |  0.95 |         1 | poor     | poor     | poor     | poor  | poor     |
|  50 |  20 |  0.2  |         0 | good     | good     | poor     | poor  | good     |
|  50 |  20 |  0.2  |         1 | poor     | poor     | poor     | poor  | poor     |
|  50 |  20 |  0.5  |         0 | moderate | moderate | poor     | poor  | poor     |
|  50 |  20 |  0.5  |         1 | poor     | poor     | poor     | poor  | poor     |
|  50 |  20 |  0.8  |         0 | moderate | poor     | poor     | poor  | poor     |
|  50 |  20 |  0.8  |         1 | poor     | poor     | poor     | poor  | poor     |
|  50 |  20 |  0.95 |         0 | poor     | poor     | poor     | poor  | poor     |
|  50 |  20 |  0.95 |         1 | poor     | poor     | poor     | poor  | poor     |
|  50 |  50 |  0.2  |         0 | moderate | good     | moderate | poor  | good     |
|  50 |  50 |  0.2  |         1 | poor     | moderate | moderate | poor  | moderate |
|  50 |  50 |  0.5  |         0 | poor     | good     | poor     | poor  | moderate |
|  50 |  50 |  0.5  |         1 | poor     | poor     | poor     | poor  | poor     |
|  50 |  50 |  0.8  |         0 | poor     | poor     | poor     | poor  | poor     |
|  50 |  50 |  0.8  |         1 | poor     | poor     | poor     | poor  | poor     |
|  50 |  50 |  0.95 |         0 | poor     | poor     | poor     | poor  | poor     |
|  50 |  50 |  0.95 |         1 | poor     | poor     | poor     | poor  | poor     |
|  50 | 100 |  0.2  |         0 | poor     | good     | good     | poor  | good     |
|  50 | 100 |  0.2  |         1 | poor     | moderate | moderate | poor  | moderate |
|  50 | 100 |  0.5  |         0 | poor     | good     | poor     | poor  | good     |
|  50 | 100 |  0.5  |         1 | poor     | moderate | poor     | poor  | moderate |
|  50 | 100 |  0.8  |         0 | poor     | poor     | poor     | poor  | moderate |
|  50 | 100 |  0.8  |         1 | poor     | poor     | poor     | poor  | poor     |
|  50 | 100 |  0.95 |         0 | poor     | poor     | poor     | poor  | poor     |
|  50 | 100 |  0.95 |         1 | poor     | poor     | poor     | poor  | poor     |
| 100 |  20 |  0.2  |         0 | moderate | good     | poor     | poor  | moderate |
| 100 |  20 |  0.2  |         1 | poor     | poor     | poor     | poor  | poor     |
| 100 |  20 |  0.5  |         0 | moderate | poor     | poor     | poor  | poor     |
| 100 |  20 |  0.5  |         1 | poor     | poor     | poor     | poor  | poor     |
| 100 |  20 |  0.8  |         0 | poor     | poor     | poor     | poor  | poor     |
| 100 |  20 |  0.8  |         1 | poor     | poor     | poor     | poor  | poor     |
| 100 |  20 |  0.95 |         0 | poor     | poor     | poor     | poor  | poor     |
| 100 |  20 |  0.95 |         1 | poor     | poor     | poor     | poor  | poor     |
| 100 |  50 |  0.2  |         0 | poor     | good     | poor     | poor  | good     |
| 100 |  50 |  0.2  |         1 | poor     | poor     | poor     | poor  | poor     |
| 100 |  50 |  0.5  |         0 | poor     | good     | poor     | poor  | moderate |
| 100 |  50 |  0.5  |         1 | poor     | poor     | poor     | poor  | poor     |
| 100 |  50 |  0.8  |         0 | poor     | poor     | poor     | poor  | poor     |
| 100 |  50 |  0.8  |         1 | poor     | poor     | poor     | poor  | poor     |
| 100 |  50 |  0.95 |         0 | poor     | poor     | poor     | poor  | poor     |
| 100 |  50 |  0.95 |         1 | poor     | poor     | poor     | poor  | poor     |
| 100 | 100 |  0.2  |         0 | poor     | good     | moderate | poor  | good     |
| 100 | 100 |  0.2  |         1 | poor     | moderate | moderate | poor  | moderate |
| 100 | 100 |  0.5  |         0 | poor     | good     | poor     | poor  | good     |
| 100 | 100 |  0.5  |         1 | poor     | poor     | poor     | poor  | poor     |
| 100 | 100 |  0.8  |         0 | poor     | poor     | poor     | poor  | poor     |
| 100 | 100 |  0.8  |         1 | poor     | poor     | poor     | poor  | poor     |
| 100 | 100 |  0.95 |         0 | poor     | poor     | poor     | poor  | poor     |
| 100 | 100 |  0.95 |         1 | poor     | poor     | poor     | poor  | poor     |
| 200 |  20 |  0.2  |         0 | moderate | good     | poor     | poor  | poor     |
| 200 |  20 |  0.2  |         1 | poor     | poor     | poor     | poor  | poor     |
| 200 |  20 |  0.5  |         0 | poor     | poor     | poor     | poor  | poor     |
| 200 |  20 |  0.5  |         1 | poor     | poor     | poor     | poor  | poor     |
| 200 |  20 |  0.8  |         0 | poor     | poor     | poor     | poor  | poor     |
| 200 |  20 |  0.8  |         1 | poor     | poor     | poor     | poor  | poor     |
| 200 |  20 |  0.95 |         0 | poor     | poor     | poor     | poor  | poor     |
| 200 |  20 |  0.95 |         1 | poor     | poor     | poor     | poor  | poor     |
| 200 |  50 |  0.2  |         0 | poor     | good     | poor     | poor  | good     |
| 200 |  50 |  0.2  |         1 | poor     | poor     | poor     | poor  | poor     |
| 200 |  50 |  0.5  |         0 | poor     | good     | poor     | poor  | moderate |
| 200 |  50 |  0.5  |         1 | poor     | poor     | poor     | poor  | poor     |
| 200 |  50 |  0.8  |         0 | poor     | poor     | poor     | poor  | poor     |
| 200 |  50 |  0.8  |         1 | poor     | poor     | poor     | poor  | poor     |
| 200 |  50 |  0.95 |         0 | poor     | poor     | poor     | poor  | poor     |
| 200 |  50 |  0.95 |         1 | poor     | poor     | poor     | poor  | poor     |
| 200 | 100 |  0.2  |         0 | poor     | good     | poor     | poor  | good     |
| 200 | 100 |  0.2  |         1 | poor     | poor     | poor     | poor  | poor     |
| 200 | 100 |  0.5  |         0 | poor     | good     | poor     | poor  | good     |
| 200 | 100 |  0.5  |         1 | poor     | poor     | poor     | poor  | poor     |
| 200 | 100 |  0.8  |         0 | poor     | poor     | poor     | poor  | poor     |
| 200 | 100 |  0.8  |         1 | poor     | poor     | poor     | poor  | poor     |
| 200 | 100 |  0.95 |         0 | poor     | poor     | poor     | poor  | poor     |
| 200 | 100 |  0.95 |         1 | poor     | poor     | poor     | poor  | poor     |

_Note: good: D_I<=0.05 and |coverage-0.95|<=0.05; moderate: <=0.15 / <=0.10; poor otherwise._

### T12. Estimator-selection surrogates: ridge vs k-NN (R2, RMSE, MAE together)

| model      |   cv_rmse_D_I |   cv_rmse_cov |   r2_D_I |   r2_cov |   mae_D_I |   mae_cov |
|:-----------|--------------:|--------------:|---------:|---------:|----------:|----------:|
| ridge_a1   |        0.0204 |        0.0915 |   0.8239 |   0.8418 |    0.014  |    0.0671 |
| ridge_a100 |        0.0384 |        0.1885 |   0.3336 |   0.3011 |    0.0272 |    0.1511 |
| knn3       |        0.0284 |        0.1343 |   0.8325 |   0.797  |    0.0116 |    0.0703 |
| knn5       |        0.0266 |        0.1159 |   0.8154 |   0.8156 |    0.0138 |    0.0745 |

_Note: 5-fold CV, folds seeded 42. Winner by CV RMSE on D_I: ridge_a1._