"""Judge of step 6b: the pre-registered rules of step6b_gyro_fact.py, applied to its rows and to the step-6 rows."""
import json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
B = json.load(open(os.path.join(HERE, 'results', 'step6b_rows.json')))
R = json.load(open(os.path.join(HERE, 'results', 'step6_rows.json')))
SC = ['A_steady', 'B_throttle', 'C_gyrobias']
D2R = np.pi / 180; Q_G = (0.01 * D2R) ** 2; S_SENS = (7e-5) ** 2
S0 = lambda V: V * Q_G + S_SENS
out = []; P = lambda *a: out.append(' '.join(str(x) for x in a))
P('step6b script sha256', B['sha256'], '| step6 script sha256', R['sha256'])
P(f'Allan at rest: lambda_hat/Q_G {np.round(np.array(B["lambda_hat"]) / Q_G, 3).tolist()}; slope median {np.median(B["allan_slope"]):.3f}')

# identity check: CF-0.5 recomputed must equal the stored step-6 rows
idmax = max(np.max(np.abs(np.array(B[s]['cf05_steady']) - np.array(R[s]['cf05']['steady']))) for s in SC)
P(f'identity check CF-0.5: max |diff| = {idmax:.2e} deg ->', 'OK' if idmax < 1e-9 else 'VOID RUN')
assert idmax < 1e-9


def crit(tau, tor, rms, ror):
    rt = np.median(tau / tor, axis=1); rg = np.median(rms / ror - 1, axis=1)
    return rt, rg, (rt >= 1 / 1.5) & (rt <= 1.5) & (rg <= 0.10)


rows = {}
for s in SC:
    r = R[s]; ns = len(r['seeds'])
    kopt_all = np.array(r['kopt']['steady']).reshape(3, -1, ns)
    best = np.argmin((kopt_all ** 2).sum(0), axis=0)
    rows[s] = dict(tor=np.array(r['oracle']['tau']), ror=np.array(r['oracle']['rms']),
                   kopt=kopt_all[:, best, np.arange(ns)], kfix=np.array(r['kfix']['steady']),
                   kcl=np.array(r['kcl']['steady']), sk=np.array(r['sk']['steady']),
                   trk=np.array(r['kfix']['trans']), trc=np.array(r['kcl']['trans']), V=r['V'])
    for v in ('6b', '6b_x2', '6b_d2'):
        rows[s][v] = np.array(B[s][v]['steady']); rows[s]['t' + v] = np.array(B[s][v]['tau_med'])
        rows[s]['tr' + v] = np.array(B[s][v]['trans'])

P('\n== RMS (millideg), steady parts, mean over seeds; windows (s), seed medians ==')
for s in SC:
    x = rows[s]
    for g in range(3):
        P(f'{s} seg {g} (V {x["V"][g]}): oracle {1e3*x["ror"][g].mean():.1f} ({np.median(x["tor"][g]):.3f} s)'
          f' | SYN-6B {1e3*x["6b"][g].mean():.1f} ({np.median(x["t6b"][g]):.3f} s)'
          f' | x2 {1e3*x["6b_x2"][g].mean():.1f} ({np.median(x["t6b_x2"][g]):.3f} s)'
          f' | /2 {1e3*x["6b_d2"][g].mean():.1f} ({np.median(x["t6b_d2"][g]):.3f} s)'
          f' | SYN-K {1e3*x["sk"][g].mean():.1f} | KF-OPT {1e3*x["kopt"][g].mean():.1f}'
          f' | KF-FIX {1e3*x["kfix"][g].mean():.1f} | KF-CLAIR {1e3*x["kcl"][g].mean():.1f}')
    P(f'   transients (first 40 s), mean RMS millideg: SYN-6B {np.round(1e3*x["tr6b"].mean(1), 1).tolist()}'
      f' KF-FIX {np.round(1e3*x["trk"].mean(1), 1).tolist()} KF-CLAIR {np.round(1e3*x["trc"].mean(1), 1).tolist()}')

V = {}
P('\n== Hypotheses ==')
cells = []
for s in SC[:2]:
    x = rows[s]; rt, rg, ok = crit(x['t6b'], x['tor'], x['6b'], x['ror'])
    for g in range(3):
        cells.append(ok[g]); P(f'H6b.1 {s} seg {g}: tau/oracle {rt[g]:.3f}, regret {100*rg[g]:.1f} % ->', 'pass' if ok[g] else 'FAIL')
V['H6b.1'] = all(cells)
x = rows['B_throttle']; Vb = x['V']
sl = [np.median((np.log(x['t6b'][g]) - np.log(x['t6b'][0])) / np.log(S0(Vb[g]) / S0(Vb[0]))) for g in (1, 2)]
P(f'H6b.2 exponent seg1 {sl[0]:.3f}, seg2 {sl[1]:.3f} (0.5 +/- 0.1)'); V['H6b.2'] = all(0.4 <= v <= 0.6 for v in sl)
x = rows['A_steady']; q = [x['6b'][g].mean() / x['kopt'][g].mean() for g in range(3)]
P(f'H6b.3 A: SYN-6B / KF-OPT {np.round(q, 3).tolist()} (<= 1.10)'); V['H6b.3'] = all(v <= 1.10 for v in q)
x = rows['B_throttle']
fx = [x['6b'][g].mean() / x['kfix'][g].mean() for g in range(3)]; cl = [x['6b'][g].mean() / x['kcl'][g].mean() for g in range(3)]
P(f'H6b.4 B: SYN-6B / KF-FIX {np.round(fx, 3).tolist()} (< 1 in seg 1, 2); SYN-6B / KF-CLAIR {np.round(cl, 3).tolist()} (<= 1.10)')
V['H6b.4'] = fx[1] < 1 and fx[2] < 1 and all(v <= 1.10 for v in cl)
x = rows['C_gyrobias']; rt, rg, ok = crit(x['t6b'], x['tor'], x['6b'], x['ror'])
for g in range(3):
    P(f'H6b.5 C seg {g}: tau/oracle {rt[g]:.3f}, regret {100*rg[g]:.1f} % ->', 'pass' if ok[g] else 'FAIL')
V['H6b.5'] = all(ok)
cells = []
for v in ('6b_x2', '6b_d2'):
    for s in SC[:2]:
        x = rows[s]; rt, rg, ok = crit(x['t' + v], x['tor'], x[v], x['ror'])
        for g in range(3):
            cells.append(ok[g]); P(f'H6b.6 {v} {s} seg {g}: tau/oracle {rt[g]:.3f}, regret {100*rg[g]:.1f} % ->', 'pass' if ok[g] else 'FAIL')
V['H6b.6'] = all(cells)
P('\n== Verdict ==')
for k, v in V.items():
    P(f'{k}: {"HOLDS" if v else "FAILS"}')
txt = '\n'.join(out); print(txt); open(os.path.join(HERE, 'results', 'step6b_results.txt'), 'w').write(txt + '\n')
