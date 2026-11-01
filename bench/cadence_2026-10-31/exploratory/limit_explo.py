# Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
# Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
"""Exploratory, after loop_dynamics_explo.py: at tempo 1 the correction of RJ4 and RJ6 sits at its 5 degree limit
(0.6 to 0.7 % of the time in excerpt A, 3.7 % in B; never at 0.7 or 0.5). Check: the noise-free bench of
loop_dynamics_explo.py, frozen window, with the engine's correction limit at 5 degrees (published) and at 30 degrees;
k, the tip error, and the share of time at the limit. Usage: python limit_explo.py <bench dir>"""
import sys, os, json
import numpy as np
sys.path.insert(0, sys.argv[1]); sys.path.insert(0, os.path.join(sys.argv[1], 'exploratory'))
import plan_intervention as P, loop_dynamics_explo as L

ARMS = ('B', 'U70', 'U50')


def run(name, arm, tf, lim_deg):
    t0, tt, ph, q, load_all, w = L.plan(name, arm)
    tauS, w0 = P.TAU_S[:6], P.W0[:6]; k_inv = (1.0 * P.DEG) / tauS; hb = 0.5 * 0.3 * P.DEG
    th_g = q[0].copy(); allowed = np.zeros((6, 5), bool); allowed[:, 2] = True
    eng = P.CadenceEngine(q[0], P.DT, w0, allowed, np.zeros(6, bool), limit=np.radians(lim_deg))
    eng.floor = tf.copy(); eng.tau_max = tf.copy(); eng.tau = tf.copy(); eng.tau_target = tf.copy()
    keep, ks, hit, n = [], [], np.zeros(6), 0
    for k in range(len(tt)):
        load = load_all[k]; th_out = th_g - load * k_inv - hb * np.tanh(load / (0.02 * tauS))
        a = np.round(th_out / P.ENC_LSB) * P.ENC_LSB; cmd = eng.step(q[k], a, a)
        w_eff = w0 * np.maximum(0.05, 1.0 - np.abs(load) / tauS)
        th_g = th_g + P.DT * np.clip((cmd - th_g) / P.TS, -w_eff, w_eff)
        if tt[k] >= t0 - 1e-9: hit += np.abs(eng.u) >= eng.limit - 1e-12; n += 1
        if k % 10 == 0 and tt[k] >= t0 - 1e-9: keep.append((th_out.copy(), q[k].copy())); ks.append(k)
    g = L.grid(name, arm, keep, ks, tt, ph, t0, w)
    return g, hit / n


if __name__ == '__main__':
    HERE = sys.argv[1]; tfs = json.load(open(os.path.join(HERE, 'results', 'plan_intervention.json')))['tau_fix']; out = {}
    for name in ('A', 'B'):
        tf = np.array(tfs[name])
        for lim in (5.0, 30.0):
            R = {a: run(name, a, tf, lim) for a in ARMS}
            for a in ARMS:
                (e, wI), hit = R[a]
                out[f'{name} limit={lim:g} {a}'] = dict(k=L.kcoef(e, R['B'][0][0], wI), tip_rms_mm=float(1e3 * np.sqrt(np.mean(np.sum(e ** 2, 1)))),
                                                        share_at_limit_by_joint=hit.tolist())
                print(name, 'limit', lim, a, 'k %.3f tip %.3f mm' % (out[f'{name} limit={lim:g} {a}']['k'], out[f'{name} limit={lim:g} {a}']['tip_rms_mm']),
                      'at limit', np.round(hit, 3), flush=True)
    json.dump(out, open(os.path.join(HERE, 'exploratory', 'limit_explo.json'), 'w'), indent=1)
