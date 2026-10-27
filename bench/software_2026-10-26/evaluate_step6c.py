"""Judge of step 6c: the pre-registered rules of step6c_ablation.py, applied to its rows, to the step-6 rows (oracle and
Kalman of B) and to the step-6b rows (identity check). Written with step6c_ablation.py, before any run."""
import hashlib, json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
C = json.load(open(os.path.join(HERE, 'results', 'step6c_rows.json')))
R6 = json.load(open(os.path.join(HERE, 'results', 'step6_rows.json')))
B6 = json.load(open(os.path.join(HERE, 'results', 'step6b_rows.json')))
D2R = np.pi / 180; Q_G = (0.01 * D2R) ** 2; S_SENS = (7e-5) ** 2
S0 = lambda V: np.asarray(V) * Q_G + S_SENS
P = C['P_grid']
out = []; P_ = lambda *a: out.append(' '.join(str(x) for x in a))
P_('evaluate_step6c.py sha256', hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest())
P_('step6c script sha256', C['sha256'], '| step6 script sha256', R6['sha256'], '| step6b script sha256', B6['sha256'])
assert not C['smoke']


def rows(scn):
    """per scenario: pow rows, oracle (tau, rms per seg x seed), orc grid (seg x window x seed), Kalman rows."""
    c = C[scn]; ns = len(c['seeds'])
    src = R6[scn] if scn == 'B_throttle' else c
    G = np.array(src['grids']['tau'])
    x = dict(V=c['V'], ns=ns, G=G,
             orc=np.array(src['orc']['steady']).reshape(3, len(G), ns),
             tor=np.array(src['oracle']['tau']), ror=np.array(src['oracle']['rms']),
             kfix=np.array(src['kfix']['steady']), kcl=np.array(src['kcl']['steady']))
    kopt_all = np.array(src['kopt']['steady']).reshape(3, -1, ns)
    best = np.argmin((kopt_all ** 2).sum(0), axis=0); x['kopt'] = kopt_all[:, best, np.arange(ns)]
    for p in P:
        x[p] = np.array(c[f'pow_{p}']['steady']); x['t', p] = np.array(c[f'pow_{p}']['tau_med'])
        x['tr', p] = np.array(c[f'pow_{p}']['trans'])
    return x


J = lambda rms: float(np.sum(np.mean(rms ** 2, axis=1)))          # total steady mean-square error, deg^2


def bestfix(x):
    mse = (x['orc'] ** 2).mean(2)                                   # seg x window
    i = int(np.argmin(mse.sum(0)))
    return i, x['G'][i], x['orc'][:, i, :]


X = {s: rows(s) for s in ('B_throttle', 'W_wide')}

# H6c.1 identity
idmax = float(np.max(np.abs(X['B_throttle'][0.5] - np.array(B6['B_throttle']['6b']['steady']))))
P_(f'\nH6c.1 identity, B: max |RMS(POW-0.5) - RMS(SYN-6B)| = {idmax:.2e} deg ->', 'OK' if idmax < 1e-9 else 'VOID RUN')
V = {'H6c.1': idmax < 1e-9}

for scn, x in X.items():
    i, tbf, rbf = bestfix(x); x['bf'] = rbf; x['tbf'] = tbf
    P_(f'\n== {scn}, V {x["V"]}, sqrt(S0/Q_G) {np.round(np.sqrt(S0(x["V"]) / Q_G), 3).tolist()} s;'
       f' RMS in millidegrees (mean over seeds), windows in s (seed medians)')
    P_(f'  oracle (best fixed per segment): {np.round(1e3 * x["ror"].mean(1), 1).tolist()}  windows {np.round(np.median(x["tor"], 1), 3).tolist()}')
    P_(f'  BEST-FIX (one window for the stream, hindsight): {tbf:.3f} s -> {np.round(1e3 * rbf.mean(1), 1).tolist()}  J {J(rbf):.4e}')
    for p in P:
        P_(f'  POW-{p:<5}: {np.round(1e3 * x[p].mean(1), 1).tolist()}  windows {np.round(np.median(x["t", p], 1), 3).tolist()}'
           f'  J {J(x[p]):.4e}  J / J(POW-0.5) {J(x[p]) / J(x[0.5]):.3f}')
    P_(f'  KF-OPT {np.round(1e3 * x["kopt"].mean(1), 1).tolist()}  KF-FIX {np.round(1e3 * x["kfix"].mean(1), 1).tolist()}'
       f'  KF-CLAIR {np.round(1e3 * x["kcl"].mean(1), 1).tolist()}')
    P_(f'  transients, mean RMS: POW-0.5 {np.round(1e3 * x["tr", 0.5].mean(1), 1).tolist()}  POW-0 {np.round(1e3 * x["tr", 0.0].mean(1), 1).tolist()}')

P_('\n== Hypotheses')
ok = []
for scn, x in X.items():
    a, b, c = J(x[0.5]), J(x['bf']), J(x[0.0])
    P_(f'H6c.2 {scn}: J(POW-0.5) {a:.4e} vs J(BEST-FIX) {b:.4e} (ratio {a / b:.3f}) and J(POW-0) {c:.4e} (ratio {a / c:.3f})')
    ok.append(a < b and a < c)
V['H6c.2'] = all(ok)

x = X['W_wide']; Js = np.array([J(x[p]) for p in P]); pp = np.array(P)
i = int(np.argmin(Js))
if 0 < i < len(P) - 1:
    cf = np.polyfit(pp[i - 1:i + 2], Js[i - 1:i + 2], 2); pstar = -cf[1] / (2 * cf[0]) if cf[0] > 0 else pp[i]
else:
    pstar = pp[i]
r25, r75 = Js[P.index(0.25)] / Js[P.index(0.5)], Js[P.index(0.75)] / Js[P.index(0.5)]
P_(f'H6c.3 W: J over p {dict(zip(P, np.round(Js / Js[P.index(0.5)], 3).tolist()))} (relative to p = 0.5);'
   f' grid minimum p = {pp[i]}, p* = {pstar:.3f} (in [0.4, 0.6]); J(0.25)/J(0.5) = {r25:.3f}, J(0.75)/J(0.5) = {r75:.3f} (>= 1.05)')
V['H6c.3'] = (0.4 <= pstar <= 0.6) and r25 >= 1.05 and r75 >= 1.05
xb = X['B_throttle']; Jb = np.array([J(xb[p]) for p in P])
P_(f'   reported, B: J over p relative to p = 0.5 {dict(zip(P, np.round(Jb / Jb[P.index(0.5)], 3).tolist()))}')
tor = x['tor']; Vw = x['V']
P_(f'   reported, W: oracle exponent d log tau_oracle / d log S0, seg 1 vs 0 '
   f'{np.median(np.log(tor[1] / tor[0])) / np.log(S0(Vw[1]) / S0(Vw[0])):.3f}, seg 2 vs 0 '
   f'{np.median(np.log(tor[2] / tor[0])) / np.log(S0(Vw[2]) / S0(Vw[0])):.3f}; POW-0.5 the same '
   f'{np.median(np.log(x["t", 0.5][1] / x["t", 0.5][0])) / np.log(S0(Vw[1]) / S0(Vw[0])):.3f}, '
   f'{np.median(np.log(x["t", 0.5][2] / x["t", 0.5][0])) / np.log(S0(Vw[2]) / S0(Vw[0])):.3f}')

rt = np.median(x['t', 0.5] / tor, axis=1); rg = np.median(x[0.5] / x['ror'] - 1, axis=1)
okc = (rt >= 1 / 1.5) & (rt <= 1.5) & (rg <= 0.10)
for g in range(3):
    P_(f'H6c.4 W seg {g}: tau/oracle {rt[g]:.3f}, regret {100 * rg[g]:.1f} % ->', 'pass' if okc[g] else 'FAIL')
V['H6c.4'] = bool(okc.all())

cl = [x[0.5][g].mean() / x['kcl'][g].mean() for g in range(3)]; fx = [x[0.5][g].mean() / x['kfix'][g].mean() for g in range(3)]
op = [x[0.5][g].mean() / x['kopt'][g].mean() for g in range(3)]
P_(f'H6c.5 W: POW-0.5 / KF-CLAIR {np.round(cl, 3).tolist()} (<= 1.10); POW-0.5 / KF-FIX {np.round(fx, 3).tolist()}'
   f' (< 1 in seg 1, 2); reported POW-0.5 / KF-OPT {np.round(op, 3).tolist()}')
V['H6c.5'] = all(v <= 1.10 for v in cl) and fx[1] < 1 and fx[2] < 1

P_('\n== Verdict')
for k, v in V.items():
    P_(f'{k}: {"HOLDS" if v else "FAILS"}')
txt = '\n'.join(out); print(txt)
open(os.path.join(HERE, 'results', 'step6c_results.txt'), 'w').write(txt + '\n')
