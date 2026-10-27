"""Judge of step 6e: the pre-registered rules of step6e_uprime_powered.py applied to its rows. Written with the script,
before any run."""
import hashlib, json, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
D = json.load(open(os.path.join(HERE, 'results', 'step6e_rows.json')))
assert not D['smoke']
D2R = np.pi / 180; Q_G = (0.01 * D2R) ** 2; S_SENS = (7e-5) ** 2
S0 = lambda V: np.asarray(V, float) * Q_G + S_SENS
V = D['V']; LG = np.log(np.array(D['grid'])); ARMS = ['U50', 'U60', 'U75', 'U90']; B = 2000
out = []; P_ = lambda *a: out.append(' '.join(str(x) for x in a))
P_('evaluate_step6e.py sha256', hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest())
P_('step6e script sha256', D['sha256'])


def popt(J):
    """minimiser and minimum of a mean-square curve over the log grid (parabola through the minimum and neighbours)."""
    i = int(np.argmin(J)); edge = i in (0, len(J) - 1)
    if edge:
        return np.exp(LG[i]), J[i], True
    x0, x1, x2 = LG[i - 1:i + 2]; y0, y1, y2 = J[i - 1:i + 2]
    c = np.polyfit([x0, x1, x2], [y0, y1, y2], 2)
    if c[0] <= 0:
        return np.exp(x1), y1, False
    xm = -c[1] / (2 * c[0]); return np.exp(xm), np.polyval(c, xm), False


def stats(idx):
    """everything judged, for one draw of seed indices (the same draw for all arms)."""
    s = {}
    for arm in ARMS:
        r = A[arm]; orc = r['orc'][:, :, idx] ** 2
        t, j, e = zip(*[popt(orc[g].mean(1)) for g in range(3)])
        t = np.array(t); j = np.array(j)
        J = lambda nm: (r[nm][:, idx] ** 2).mean(1)
        s[arm] = dict(t=t, j=j, edge=e, eu=np.log(t[1] / t[0]) / np.log(S0(V[1]) / S0(V[0])),
                      ed=np.log(t[2] / t[0]) / np.log(S0(V[2]) / S0(V[0])),
                      reg={nm: np.sqrt(J(nm) / j) - 1 for nm in ('POW_U', 'POW_S', 'ULAW', 'ULAW_K')},
                      ratio=J('POW_U')[1] / J('POW_S')[1],
                      wul=np.median(r['tau_ULAW'][:, idx], axis=1) / t)
    p = np.array([A[a]['p_H'] for a in ARMS])
    s['b_up'] = np.polyfit(p, [s[a]['eu'] for a in ARMS], 1)[0]
    s['b_dn'] = np.polyfit(p, [s[a]['ed'] for a in ARMS], 1)[0]
    return s


A = {}
for arm in ARMS:
    r = D[arm]
    A[arm] = dict(p_H=r['p_H'], H=r['H'], orc=np.array(r['orc']),
                  **{nm: np.array(r[nm]['steady']) for nm in ('POW_U', 'POW_S', 'ULAW', 'ULAW_K')},
                  tau_ULAW=np.array(r['ULAW']['tau_med']), H_rest=np.array(r['H_rest']), lam_rest=np.array(r['lamH_rest_ratio']))
ns = A['U50']['orc'].shape[2]
S = stats(np.arange(ns))
rng = np.random.default_rng(2026)
BS = [stats(rng.integers(0, ns, ns)) for _ in range(B)]
ci = lambda f: np.quantile([f(b) for b in BS], [0.025, 0.975])

P_(f'\n{ns} seeds per arm, paired; bootstrap {B} draws over seeds (95 % percentile intervals)')
for arm in ARMS:
    s = S[arm]; a = A[arm]
    P_(f'\n== {arm}: H = {a["H"]}, p_H = {a["p_H"]:.3f}; H at rest mean {a["H_rest"].mean():.3f} (range {a["H_rest"].min():.3f}..{a["H_rest"].max():.3f}),'
       f' lambda_H at rest / true median {np.median(a["lam_rest"]):.3f}')
    P_(f'  population optimum tau_opt (s): {np.round(s["t"], 3).tolist()}  at a grid edge: {list(s["edge"])}'
       f'  RMS at the optimum (mdeg): {np.round(1e3 * np.sqrt(s["j"]), 1).tolist()}')
    P_(f'  exponent up {s["eu"]:.3f} [{ci(lambda b: b[arm]["eu"])[0]:.3f}, {ci(lambda b: b[arm]["eu"])[1]:.3f}]'
       f'   down {s["ed"]:.3f} [{ci(lambda b: b[arm]["ed"])[0]:.3f}, {ci(lambda b: b[arm]["ed"])[1]:.3f}]   predicted {a["p_H"]:.3f}')
    for nm in ('POW_U', 'POW_S', 'ULAW', 'ULAW_K'):
        P_(f'  {nm:6s} regret % per segment {np.round(100 * s["reg"][nm], 1).tolist()}')
    P_(f'  ULAW window / tau_opt (seed median) {np.round(s["wul"], 3).tolist()}')
    P_(f'  J(POW-U) / J(POW-S) at V = 256: {s["ratio"]:.3f} [{ci(lambda b: b[arm]["ratio"])[0]:.3f}, {ci(lambda b: b[arm]["ratio"])[1]:.3f}]')

P_('\n== Hypotheses')
V_ = {}
cb = ci(lambda b: b['b_up'])
P_(f'H6e.1 slope of the upward exponents on p_H: b = {S["b_up"]:.3f} [{cb[0]:.3f}, {cb[1]:.3f}]; in [0.7, 1.3]: {0.7 <= S["b_up"] <= 1.3};'
   f' interval excludes 0: {cb[0] > 0 or cb[1] < 0}')
V_['H6e.1'] = (0.7 <= S['b_up'] <= 1.3) and (cb[0] > 0 or cb[1] < 0)
cd = ci(lambda b: b['b_dn'])
P_(f'   reported: slope of the downward exponents b = {S["b_dn"]:.3f} [{cd[0]:.3f}, {cd[1]:.3f}]')
ok = []
for arm in ARMS:
    d = abs(S[arm]['eu'] - A[arm]['p_H']); ok.append(d <= 0.05)
    P_(f'H6e.2 {arm}: |e - p_H| = {d:.3f} (<= 0.05: {d <= 0.05})')
V_['H6e.2'] = all(ok)
ok = []
for arm in ('U75', 'U90'):
    c = ci(lambda b: b[arm]['eu']); ok.append(c[1] < 0.5 or c[0] > 0.5)
    P_(f'H6e.3 {arm}: upward exponent interval [{c[0]:.3f}, {c[1]:.3f}] excludes 0.5: {ok[-1]}')
V_['H6e.3'] = all(ok)
ok = []
for arm in ('U75', 'U90'):
    c = ci(lambda b: b[arm]['ratio']); ok.append(S[arm]['ratio'] < 1 and c[1] < 1)
    P_(f'H6e.4 {arm}: J(POW-U)/J(POW-S) at V = 256 = {S[arm]["ratio"]:.3f} [{c[0]:.3f}, {c[1]:.3f}], entirely below 1: {ok[-1]}')
V_['H6e.4'] = all(ok)
ok = []
for arm in ARMS:
    rg = S[arm]['reg']['ULAW']; wr = S[arm]['wul']
    good = bool(np.all(rg <= 0.10) and np.all((wr >= 1 / 1.5) & (wr <= 1.5))); ok.append(good)
    P_(f'H6e.5 {arm}: ULAW regret % {np.round(100 * rg, 1).tolist()}, window ratio {np.round(wr, 3).tolist()} -> {"pass" if good else "FAIL"}')
V_['H6e.5'] = all(ok)
P_('\n== Verdict')
for k, v in V_.items():
    P_(f'{k}: {"HOLDS" if v else "FAILS"}')
txt = '\n'.join(out); print(txt)
open(os.path.join(HERE, 'results', 'step6e_results.txt'), 'w').write(txt + '\n')
