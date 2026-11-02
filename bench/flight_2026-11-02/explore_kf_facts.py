"""Exploratory, after the verdicts (not pre-registered): a Kalman set from the same facts the law uses, without any
tuning against truth (R = var(tilt) per sample, Q_G = lambda at idle / at half throttle), estimation error on the test
seeds, next to the law (c = 1 and c = 2.8) and the truth-tuned Kalman. Amendment 2's noise model."""
import os, sys, json, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
R = os.path.abspath(os.path.join(HERE, '..', '..')); sys.argv = ['x', R, 'import']
import flight_bench as F
G = F.G
def records_a2():
    if not F._REC:
        for name, seg in (('idle', G.SEGS[0]), ('low', G.SEGS[1]), ('half', G.SEGS[2])):
            tilt, _, g1k = G.load(*seg[:4]); F._REC[name] = (tilt - tilt.mean(), g1k)
    return F._REC
F.records = records_a2
OUT = os.path.join(HERE, 'results_a2'); fx = json.load(open(os.path.join(OUT, 'facts.json'))); tu = json.load(open(os.path.join(OUT, 'tuned_estimators.json')))
vt = fx['var_tilt']; li = fx['lam']['idle']; lh = fx['lam']['half']
cols = dict(KF_facts_idle=(1, 1.0, vt, li), KF_facts_half=(1, 1.0, vt, lh), KF_tuned=(1, 1.0, tu['kf_R'], tu['kf_QG']),
            LAW_c1=(2, 1.0, 1.0, 1.0), LAW_c2p8=(2, 2.8, 1.0, 1.0))
names = list(cols); kind = np.array([cols[n][0] for n in names]); cf = np.array([cols[n][1] for n in names])
kR = np.array([cols[n][2] for n in names]); kQ = np.array([cols[n][3] for n in names])
J = np.array([F.est_only(s, fx, kind, cfac=cf, kR=kR, kQ=kQ) for s in F.TEST])
P = np.sqrt(np.mean(J ** 2, 0)); out = {n: float(p) for n, p in zip(names, P)}
json.dump(dict(pooled_est_test=out, R_mult_facts=1.0, note='exploratory'), open(os.path.join(OUT, 'explore_kf_facts.json'), 'w'), indent=1)
print({n: round(v, 3) for n, v in out.items()})
