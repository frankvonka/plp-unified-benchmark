import sys, numpy as np, pandas as pd
sys.path.insert(0, '/root/projects/plp-unified-benchmark/code')
import importlib.util
spec = importlib.util.spec_from_file_location("est", "/root/projects/plp-unified-benchmark/code/01_estimators.py")
est = importlib.util.module_from_spec(spec); spec.loader.exec_module(est)

rng = np.random.default_rng(42)
T, N = 30, 25
alpha = rng.normal(0, 1, N)
X = np.zeros((T, N)); Y = np.zeros((T, N))
for i in range(N):
    x = np.zeros(T); x[0] = rng.normal()
    for t in range(1, T): x[t] = 0.5*x[t-1] + rng.normal()
    X[:, i] = x
    Y[:, i] = alpha[i] + (-0.6)*x + rng.normal(size=T)

# CHECK 1: FE vs linearmodels
b_fe, se_fe, V = est.fe_core(Y, [X], return_V=True)
import linearmodels.panel as lmp
df = pd.DataFrame({'y': Y.T.reshape(-1), 'x': X.T.reshape(-1)})
df['id'] = np.repeat(np.arange(N), T); df['t'] = np.tile(np.arange(T), N)
df = df.set_index(['id','t'])
m = lmp.PanelOLS(df['y'], df[['x']], entity_effects=True).fit(cov_type='clustered', cluster_entity=True)
print('FE coef ours %.8f | lm %.8f' % (b_fe[0], m.params['x']))
print('FE se   ours %.8f | lm %.8f' % (se_fe[0], m.std_errors['x']))
assert abs(b_fe[0]-m.params['x']) < 1e-9, 'FE coef mismatch'
assert abs(se_fe[0]-m.std_errors['x'])/m.std_errors['x'] < 0.25, 'FE se far off'  # convention band
print('CHECK 1 PASSED')

# CHECK 2: SPJ = manual jackknife formula
b_spj, se_spj = est.spj_core(Y, [X])
yd, Xd, _ = est._design(Y, [X], None, False)
yv, Xv, ents, times, ok = est._stack(yd, Xd)
cf, _, _ = est._ols_fit(yv, Xv)
c = int(np.median(np.arange(T)))
def half_b(sel):
    d = Y.copy(); xx = X.copy()[:, :, None].copy()
    d[~sel, :] = np.nan; xx[~sel, :, :] = np.nan
    yd2, Xd2, _ = est._design(d, [xx[:, :, 0]], None, False)
    yv2, Xv2, _, _, ok2 = est._stack(yd2, Xd2)
    c2, _, _ = est._ols_fit(yv2, Xv2)
    return c2[0]
man = 2*cf[0] - 0.5*(half_b(np.arange(T) <= c) + half_b(np.arange(T) > c))
print('SPJ %.6f | manual %.6f' % (b_spj[0], man))
assert abs(b_spj[0]-man) < 1e-10
print('CHECK 2 PASSED')

# CHECK 3: HPJ / CDOLS / DB smoke
b_hpj, se_hpj = est.hpj_core(Y, [X])
b_cd, se_cd = est.cdols_h(Y, X, h=0)
b_cd1, se_cd1 = est.cdols_h(Y, X, h=1)
print('HPJ %.4f  CDOLS h0 %.4f se %.4f  h1 %.4f se %.4f' % (b_hpj[0], b_cd, se_cd, b_cd1, se_cd1))
b_db, se_db = est.db_correction(Y, X, b_fe, se_fe)
print('DB:', np.round(b_db[:4], 4), 'FE:', np.round(b_fe[:4], 4))

# CHECK 4: bias direction sanity (Monte Carlo mini-run): FE should be biased DOWN
# in magnitude for beta<0 with persistent x (their Fig 2 pattern); SPJ/HPJ closer to true.
def one_dgp(seed, T=30, N=25, rho=0.8, beta=-0.6):
    r = np.random.default_rng(seed)
    a = r.normal(0, 0.5, N)          # alpha_i
    Xw = np.zeros((T, N)); Yw = np.zeros((T, N))
    for i in range(N):
        x = np.zeros(T); x[0] = 0.2 + rho*r.normal() + r.normal()
        for t in range(1, T): x[t] = 0.2 + rho*x[t-1] + r.normal()
        # alpha correlated with mean(x): their eta*sqrt(T)*Xbar
        al = a[i] + 0.2*np.sqrt(T)*x.mean()
        Xw[:, i] = x
        Yw[:, i] = al + beta*x + r.normal(size=T)
    return Yw, Xw
B = {k: [] for k in ['FE','SPJ','HPJ','CDOLS']}
for s in range(40):
    Yw, Xw = one_dgp(1000+s)
    b, se = est.fe_core(Yw, [Xw]); B['FE'].append(b[0])
    b, se = est.spj_core(Yw, [Xw]); B['SPJ'].append(b[0])
    b, se = est.hpj_core(Yw, [Xw]); B['HPJ'].append(b[0])
    b, se = est.cdols_h(Yw, Xw, h=0); B['CDOLS'].append(b)
for k, v in B.items():
    v = np.array(v)
    print('%6s mean beta_hat = %.4f  (bias %+0.4f)' % (k, v.mean(), v.mean()+0.6))
print('TRUE beta = -0.6')
print('ALL CHECKS DONE')
