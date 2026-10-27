"""Judge of step 6: applies the pre-registered decision rules of step6_quad_attitude_estimator.py to its rows.
Nothing here is tuned; every threshold is the one written in the step-6 docstring."""
import json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
fn = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, 'results', 'step6_rows.json')
R = json.load(open(fn))
SC = ['A_steady', 'B_throttle', 'C_gyrobias']
D2R = np.pi / 180
Q_G = (0.01 * D2R) ** 2; S_SENS = (7e-5) ** 2
S0 = lambda V: V * Q_G + S_SENS
out = []
P = lambda *a: out.append(' '.join(str(x) for x in a))


def arr(res, key, what='steady'):
    return np.array(res[key][what])                                   # (3, ncol)


def crit(tau, tau_or, rms, rms_or):
    rt = np.median(tau / tau_or, axis=1); rg = np.median(rms / rms_or - 1, axis=1)
    ok = (rt >= 1 / 1.5) & (rt <= 1.5) & (rg <= 0.10)
    return rt, rg, ok


P('step6 rows', os.path.basename(fn), '| script sha256', R['sha256'], 'smoke' if R['smoke'] else 'full')
verdict = {}
rows = {}
for s in SC:
    r = R[s]; ns = len(r['seeds'])
    tor = np.array(r['oracle']['tau']); ror = np.array(r['oracle']['rms'])
    s1 = arr(r, 's1'); t1 = np.array(r['s1']['tau_med']); H1 = np.array(r['s1']['H_med'])
    s2r = arr(r, 's2r'); t2 = np.array(r['s2r']['tau_med'])
    sk = arr(r, 'sk'); tk = np.array(r['tau_star_true'])[:, None] * np.ones((1, ns))
    snoi = arr(r, 'snoi'); tn = np.array(r['snoi']['tau_med'])
    kopt_all = arr(r, 'kopt').reshape(3, -1, ns)
    best = np.argmin((kopt_all ** 2).sum(0), axis=0)                 # one R per stream
    kopt = kopt_all[:, best, np.arange(ns)]
    kfix = arr(r, 'kfix'); kcl = arr(r, 'kcl'); cf05 = arr(r, 'cf05')
    mad_all = arr(r, 'mad').reshape(3, -1, ns)
    bm = np.argmin((mad_all[:, :-1] ** 2).sum(0), axis=0); mad = mad_all[:, bm, np.arange(ns)]; mad004 = mad_all[:, -1]
    rows[s] = dict(tor=tor, ror=ror, s1=s1, t1=t1, H1=H1, s2r=s2r, t2=t2, sk=sk, tk=tk, snoi=snoi, tn=tn, kopt=kopt,
                   kfix=kfix, kcl=kcl, cf05=cf05, mad=mad, mad004=mad004, best=best, bm=bm,
                   R_best=np.array(r['grids']['R_mult'])[best], b2_best=np.array(r['grids']['beta2'])[bm],
                   tr_s1=arr(r, 's1', 'trans'), tr_kfix=arr(r, 'kfix', 'trans'), tr_kcl=arr(r, 'kcl', 'trans'),
                   tr_or=None, V=r['V'])

P('\n== Table: RMS angle error on steady parts (deg), mean over seeds; windows (s), seed medians ==')
for s in SC:
    x = rows[s]; P(f'\n{s}  V = {x["V"]}  sqrt(S0/Q_G) = {np.round(R[s]["tau_star_true"], 3).tolist()}')
    for g in range(3):
        P(f'  seg {g}: oracle {x["ror"][g].mean():.4f} (tau {np.median(x["tor"][g]):.3f})'
          f' | SYN-1S {x["s1"][g].mean():.4f} (tau {np.median(x["t1"][g]):.3f}, H {np.median(x["H1"][g]):.2f})'
          f' | SYN-K {x["sk"][g].mean():.4f} | SYN-2R {x["s2r"][g].mean():.4f} (tau {np.median(x["t2"][g]):.4f})'
          f' | SYN-noI {x["snoi"][g].mean():.4f} (tau {np.median(x["tn"][g]):.3f})')
        P(f'          KF-OPT {x["kopt"][g].mean():.4f} | KF-FIX {x["kfix"][g].mean():.4f} | KF-CLAIR {x["kcl"][g].mean():.4f}'
          f' | MADG-best {x["mad"][g].mean():.4f} | MADG-0.04 {x["mad004"][g].mean():.4f} | CF-0.5 {x["cf05"][g].mean():.4f}')
    P(f'  KF-OPT R/R_nom chosen (per seed): {np.round(x["R_best"], 3).tolist()}')
    P(f'  MADG best 2beta (per seed): {np.round(x["b2_best"], 5).tolist()}')
    P(f'  transients (first 40 s of each segment), mean RMS: SYN-1S {np.round(x["tr_s1"].mean(1), 4).tolist()}'
      f' KF-FIX {np.round(x["tr_kfix"].mean(1), 4).tolist()} KF-CLAIR {np.round(x["tr_kcl"].mean(1), 4).tolist()}')

P('\n== Hypotheses ==')
# H6.1
cells = []
for s in SC:
    x = rows[s]; rt, rg, ok = crit(x['tk'], x['tor'], x['sk'], x['ror'])
    for g in range(3):
        cells.append(ok[g]); P(f'H6.1 {s} seg {g}: SYN-K tau/oracle {rt[g]:.3f}, regret {100 * rg[g]:.1f} % -> {"pass" if ok[g] else "FAIL"}')
verdict['H6.1'] = all(cells)
# H6.2
cells = []
for s in SC[:2]:
    x = rows[s]; rt, rg, ok = crit(x['t1'], x['tor'], x['s1'], x['ror'])
    for g in range(3):
        cells.append(ok[g]); P(f'H6.2 {s} seg {g}: SYN-1S tau/oracle {rt[g]:.3f}, regret {100 * rg[g]:.1f} % -> {"pass" if ok[g] else "FAIL"}')
verdict['H6.2'] = all(cells)
# H6.3
x = rows['B_throttle']; V = x['V']
sl = [np.median((np.log(x['t1'][g]) - np.log(x['t1'][0])) / np.log(S0(V[g]) / S0(V[0]))) for g in (1, 2)]
Hok = all(0.4 <= np.median(rows[s]['H1'][g]) <= 0.6 for s in SC[:2] for g in range(3))
P(f'H6.3 exponent d log tau / d log S0: seg1 {sl[0]:.3f}, seg2 {sl[1]:.3f} (0.5 +/- 0.1); H in [0.4, 0.6] all A, B cells: {Hok}')
verdict['H6.3'] = all(0.4 <= v <= 0.6 for v in sl) and Hok
# H6.4
x = rows['A_steady']; q = [x['s1'][g].mean() / x['kopt'][g].mean() for g in range(3)]
P(f'H6.4 A: RMS SYN-1S / KF-OPT per segment {np.round(q, 3).tolist()} (<= 1.10)')
verdict['H6.4'] = all(v <= 1.10 for v in q)
# H6.5
x = rows['B_throttle']
fx = [x['s1'][g].mean() < x['kfix'][g].mean() for g in (1, 2)]
cl = [x['s1'][g].mean() / x['kcl'][g].mean() for g in range(3)]
P(f'H6.5 B: SYN-1S < KF-FIX in seg 1, 2: {fx}  (ratios {np.round([x["s1"][g].mean() / x["kfix"][g].mean() for g in range(3)], 3).tolist()});'
  f' SYN-1S / KF-CLAIR {np.round(cl, 3).tolist()} (<= 1.10)')
verdict['H6.5'] = all(fx) and all(v <= 1.10 for v in cl)
# H6.6
rt, rg, ok = crit(x['t2'], x['tor'], x['s2r'], x['ror'])
P(f'H6.6 B: SYN-2R tau/oracle {np.round(rt, 4).tolist()}, regret % {np.round(100 * rg, 1).tolist()}; fails seg 0 and 1: {not ok[0] and not ok[1]}')
verdict['H6.6'] = (not ok[0]) and (not ok[1])
# H6.7
x = rows['C_gyrobias']; rt, rg, ok = crit(x['t1'], x['tor'], x['s1'], x['ror'])
rn = [x['snoi'][g].mean() / x['s1'][g].mean() for g in range(3)]
for g in range(3):
    P(f'H6.7 C seg {g}: SYN-1S tau/oracle {rt[g]:.3f}, regret {100 * rg[g]:.1f} % -> {"pass" if ok[g] else "FAIL"}; noI/1S {rn[g]:.2f}')
verdict['H6.7'] = all(ok) and all(v >= 2 for v in rn)
P('\n== Verdict ==')
for k, v in verdict.items():
    P(f'{k}: {"HOLDS" if v else "FAILS"}')
txt = '\n'.join(out); print(txt)
open(os.path.join(HERE, 'results', 'step6_results.txt' if not R['smoke'] else 'step6_smoke_results.txt'), 'w').write(txt + '\n')
