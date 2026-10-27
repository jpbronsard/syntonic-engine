"""Judge of step 6d: the pre-registered rules of step6d_uprime_mechanism.py applied to its rows, and to step 6c's rows
for the identity check. Written with step6d_uprime_mechanism.py, before any run."""
import hashlib, json, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
D = json.load(open(os.path.join(HERE, 'results', 'step6d_rows.json')))
C = json.load(open(os.path.join(HERE, 'results', 'step6c_rows.json')))
D2R = np.pi / 180; Q_G = (0.01 * D2R) ** 2; S_SENS = (7e-5) ** 2
S0 = lambda V: np.asarray(V, float) * Q_G + S_SENS
out = []; P_ = lambda *a: out.append(' '.join(str(x) for x in a))
P_('evaluate_step6d.py sha256', hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest())
P_('step6d script sha256', D['sha256'], '| step6c script sha256', C['sha256'])
assert not D['smoke']


def crit(tau, tor, rms, ror):
    rt = np.median(np.asarray(tau) / tor, axis=1); rg = np.median(np.asarray(rms) / ror - 1, axis=1)
    return rt, rg, (rt >= 1 / 1.5) & (rt <= 1.5) & (rg <= 0.10)


def expo(t, V, g):
    t = np.asarray(t, float); return float(np.median(np.log(t[g] / t[0])) / np.log(S0(V[g]) / S0(V[0])))


V = {}
idmax = float(np.max(np.abs(np.array(D['REPLAY']['LAW6B']['steady']) - np.array(C['W_wide']['pow_0.5']['steady']))))
P_(f'\nH6d.0 identity: max |RMS(REPLAY) - RMS(step 6c POW-0.5)| = {idmax:.2e} deg ->', 'OK' if idmax < 1e-9 else 'VOID RUN')
V['H6d.0'] = idmax < 1e-9

R = {}
for arm in ('REPLAY', 'M0', 'M1', 'U60', 'U75'):
    r = D[arm]; Vw = r['V']; tor = np.array(r['oracle']['tau']); ror = np.array(r['oracle']['rms'])
    P_(f'\n== {arm}  (V {Vw}; RMS in millidegrees, mean over seeds; windows in s, seed medians)')
    if 'H_true' in r:
        P_(f'  at rest: H {np.round(r["H_rest"], 3).tolist()} (true {r["H_true"]}); lambda_H / true '
           f'{np.round(np.array(r["lamH_rest"]) / r["lamH_true"], 3).tolist()}; U\' predicted windows {np.round(r["tau_pred"], 3).tolist()} s,'
           f' exponent 1/(2H+1) = {1 / (2 * r["H_true"] + 1):.3f}')
    P_(f'  oracle: {np.round(1e3 * ror.mean(1), 1).tolist()}  windows {np.round(np.median(tor, 1), 3).tolist()}'
       f'  exponent 1->256 {expo(tor, Vw, 1):.3f}, 1->1/16 {expo(tor, Vw, 2):.3f}')
    R[arm] = {'oexp': expo(tor, Vw, 1)}
    for nm in ('LAW6B', 'ULAW_K', 'ULAW', 'SQRT'):
        if nm not in r: continue
        rms = np.array(r[nm]['steady']); tau = np.array(r[nm]['tau_med'])
        rt, rg, ok = crit(tau, tor, rms, ror)
        R[arm][nm] = ok
        P_(f'  {nm:7s}: {np.round(1e3 * rms.mean(1), 1).tolist()}  windows {np.round(np.median(tau, 1), 3).tolist()}'
           f'  tau/oracle {np.round(rt, 3).tolist()}  regret % {np.round(100 * rg, 1).tolist()}  pass {ok.tolist()}'
           f'  exponent 1->256 {expo(tau, Vw, 1):.3f}')

P_('\n== Hypotheses')
for h, arm in (('H6d.1', 'M0'), ('H6d.2', 'M1')):
    a, e = R[arm]['LAW6B'], R[arm]['oexp']
    P_(f'{h} {arm}: 6b law passes in all three segments: {bool(a.all())}; oracle exponent 1->256 {e:.3f} in [0.45, 0.55]: {0.45 <= e <= 0.55}')
    V[h] = bool(a.all()) and 0.45 <= e <= 0.55
V['H6d.3'] = all(bool(R[a]['ULAW_K'].all()) for a in ('U60', 'U75'))
P_(f'H6d.3 ULAW-K passes in all segments: U60 {bool(R["U60"]["ULAW_K"].all())}, U75 {bool(R["U75"]["ULAW_K"].all())}')
ok4 = []
for arm in ('U60', 'U75'):
    pred = 1 / (2 * D[arm]['H_true'] + 1); e = R[arm]['oexp']
    P_(f'H6d.4 {arm}: oracle exponent 1->256 {e:.3f} against 1/(2H+1) = {pred:.3f} (within 0.05: {abs(e - pred) <= 0.05})')
    ok4.append(abs(e - pred) <= 0.05)
V['H6d.4'] = all(ok4)
V['H6d.5'] = all(bool(R[a]['ULAW'].all()) for a in ('U60', 'U75'))
P_(f'H6d.5 ULAW passes in all segments: U60 {bool(R["U60"]["ULAW"].all())}, U75 {bool(R["U75"]["ULAW"].all())}')
V['H6d.6'] = not bool(R['U75']['SQRT'][1])
P_(f'H6d.6 U75: SQRT fails the criterion at V = 256: {not bool(R["U75"]["SQRT"][1])}')

P_('\n== Verdict')
for k, v in V.items():
    P_(f'{k}: {"HOLDS" if v else "FAILS"}')
txt = '\n'.join(out); print(txt)
open(os.path.join(HERE, 'results', 'step6d_results.txt'), 'w').write(txt + '\n')
