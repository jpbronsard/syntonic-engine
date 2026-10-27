"""EXPLORATORY, written after the step-6d verdict (not pre-registered). Two readings of the rows:
(a) the cost of the window's movement: in each arm M segment, the law's error against a fixed window of the same
    length (the oracle sweep interpolated, per seed, at that seed's median law window). If the law loses to a fixed
    window of its own length, the loss comes from the window moving, not from its value.
(b) the change of the optimal exponent with H, taking M1 (angle random walk drift, bias known, no bias integral: a
    first-order window on an H = 1/2 drift, the class of arm U) as the H = 1/2 reference, against the change U'
    predicts, 1/(2H+1) - 1/2.
Reads results/step6d_rows.json; writes results/step6d_exploratory.txt."""
import json, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
D = json.load(open(os.path.join(HERE, 'results', 'step6d_rows.json')))
D2R = np.pi / 180; Q_G = (0.01 * D2R) ** 2; S_SENS = (7e-5) ** 2
S0 = lambda V: np.asarray(V, float) * Q_G + S_SENS
out = []; P_ = lambda *a: out.append(' '.join(str(x) for x in a))
G = np.log(np.array(D['M0']['grids']['tau']))

P_('(a) law against a fixed window of its own length (per seed, oracle sweep interpolated in log tau); RMS in millidegrees')
for arm in ('REPLAY', 'M0', 'M1'):
    r = D[arm]; ns = len(r['seeds'])
    orc = np.array(r['orc']['steady']).reshape(3, len(G), ns)
    law = np.array(r['LAW6B']['steady']); tl = np.array(r['LAW6B']['tau_med'])
    for g in range(3):
        fx = np.array([np.interp(np.log(tl[g][j]), G, orc[g, :, j]) for j in range(ns)])
        P_(f'  {arm:6s} seg {g}: law {1e3 * law[g].mean():.1f}, fixed at its own window {1e3 * fx.mean():.1f}'
           f' (law / fixed {law[g].mean() / fx.mean():.3f}; law worse in {int((law[g] > fx).sum())} seeds of {ns})')

P_('\n(b) change of the oracle exponent with H, M1 as the H = 1/2 reference')
def ex(arm, g):
    t = np.array(D[arm]['oracle']['tau']); V = D[arm]['V']
    return float(np.median(np.log(t[g] / t[0])) / np.log(S0(V[g]) / S0(V[0])))
for g, lab in ((1, '1 -> 256'), (2, '1 -> 1/16')):
    e0 = ex('M1', g)
    for arm, H in (('U60', 0.6), ('U75', 0.75)):
        P_(f'  {lab}: {arm} {ex(arm, g):.3f} - M1 {e0:.3f} = {ex(arm, g) - e0:+.3f}; U\' predicts {1 / (2 * H + 1) - 0.5:+.3f}')
txt = '\n'.join(out); print(txt)
open(os.path.join(HERE, 'results', 'step6d_exploratory.txt'), 'w').write(txt + '\n')
