# Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
# Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
"""Exploratory, 2 Nov 2026, after the verdicts (not pre-registered): the law's best factor c regime by regime (low throttle: 40 s record; half: the 8 s looped record).
If the drift fact at half throttle is inflated, c* should sit near 1 at low throttle and well above 1 at half."""
import os, sys, json, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
sys.argv = ['x', os.path.abspath(os.path.join(HERE, '..', '..')), 'import']
import flight_bench as F
G = F.G
def rec():
    if not F._REC:
        for name, seg in (('idle', G.SEGS[0]), ('low', G.SEGS[1]), ('half', G.SEGS[2])):
            tilt, _, g1k = G.load(*seg[:4]); F._REC[name] = (tilt - tilt.mean(), g1k)
    return F._REC
F.records = rec
fx = json.load(open(os.path.join(HERE, 'results_a2', 'facts.json')))
CAL = np.array([0.5, 0.7, 1.0, 1.4, 2.0, 2.8, 4.0, 5.6, 8.0]); n = len(CAL)
def run(seed):
    gn, tn, reg = F.build_noise(seed); N = len(gn); kpre = int(round(F.T_PRE / F.DT))
    est = F.Est(n, np.full(n, 2), cfac=CAL); lam_of = {r: fx['lam'][r] for r in ('idle', 'low', 'half')}
    gs = 0.0; gc = 0; bc = 0.0; ss = {'low': np.zeros(n), 'half': np.zeros(n)}; m = {'low': 0, 'half': 0}
    for k in range(N):
        g = np.full(n, gn[k] - bc); a = np.full(n, tn[k // F.NSUB]) if k % F.NSUB == 0 else None
        th = est.step(g, a, lam_of[reg[k]], k * F.DT)
        if k < kpre:
            gs += gn[k]; gc += 1
            if k == kpre - 1: bc = gs / gc
            continue
        ss[reg[k]] += (th / F.D2R) ** 2; m[reg[k]] += 1
    return {r: ss[r] / m[r] for r in ss}
out = {}
for name, seeds in (('train', F.TRAIN), ('test', F.TEST)):
    rows = [run(s) for s in seeds]
    for r in ('low', 'half'):
        J = np.sqrt(np.mean([x[r] for x in rows], 0)); i = int(np.argmin(J))
        out[f'{name}_{r}'] = dict(c_star=float(CAL[i]), curve={str(c): round(float(v), 4) for c, v in zip(CAL, J)})
        print(name, r, 'c* =', CAL[i], {str(c): round(float(v), 3) for c, v in zip(CAL, J)}, flush=True)
print('lambda facts (rad^2/s):', fx['lam'])
json.dump(out, open(os.path.join(HERE, 'results_a2', 'explore_c_by_regime.json'), 'w'), indent=1)
