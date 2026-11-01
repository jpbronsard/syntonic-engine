# Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
# Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
"""Exploratory, after plan_intervention.py's verdict (1674429; go from JP): is the extra tip error locked to the plan
the servo's slew limit? Seeds 1 to 5 (all used before; a later test takes fresh seeds), excerpts A and B, the
frozen-window engine with the windows of results/plan_intervention.json, arms B, U70, U50. Two plants: the published
one, and the same with the slew limit removed. Reads: the share of time each joint spends at its slew limit; the share
of the tip mean square in the samples that follow a joint at its limit within 50 ms; and, with the limit removed, the
scaling coefficient k of plan_intervention.py. Not a test: a reading for the next pre-registration. The published
plant here is plan_intervention.run with a recorder (checked identical on one run).
Usage: python slew_explo.py <bench dir>      (about 15 minutes)"""
import sys, os, json, time
import numpy as np
from scipy.interpolate import CubicSpline
sys.path.insert(0, sys.argv[1])
import plan_intervention as P

SEEDS = [1, 2, 3, 4, 5]; ARMS = ('B', 'U70', 'U50'); J = 6


def run_sat(name, seed, arm, tau_fix, slew=True):
    t0 = float(P.Z[f'{name}_t0']); T, phi, om, omd = P.tempo(arm)
    tt = np.arange(t0 - P.PRE, t0 + T, P.DT)
    ph = tt if arm == 'B' else t0 + phi(tt - t0)
    q = CubicSpline(P.Z[f'{name}_t30'], P.Z[f'{name}_q30'][:, :J], axis=0)(ph)
    tau = CubicSpline(P.Z[f'{name}_t30'], P.Z[f'{name}_tau30'][:, :J], axis=0)(ph)
    sd = np.full(len(tt), P.ENC_SD); rng = np.random.default_rng(seed)
    tauS, w0 = P.TAU_S[:J], P.W0[:J]; k_inv = (1.0 * P.DEG) / tauS; hb = 0.5 * 0.3 * P.DEG
    th_g = q[0].copy()
    allowed = np.zeros((J, 5), bool); allowed[:, 2] = True
    eng = P.CadenceEngine(q[0], P.DT, w0, allowed, np.zeros(J, bool))
    tf = np.asarray(tau_fix, float); eng.floor = tf.copy(); eng.tau_max = tf.copy(); eng.tau = tf.copy(); eng.tau_target = tf.copy()
    clip = np.zeros((len(tt), J), bool); keep, ks = [], []
    for k in range(len(tt)):
        load = tau[k]
        th_out = th_g - load * k_inv - hb * np.tanh(load / (0.02 * tauS))
        a = np.round((th_out + rng.normal(0.0, 1.0, J) * sd[k]) / P.ENC_LSB) * P.ENC_LSB
        b = np.round((th_out + rng.normal(0.0, 1.0, J) * sd[k]) / P.ENC_LSB) * P.ENC_LSB
        cmd = eng.step(q[k], a, b)
        w_eff = w0 * np.maximum(0.05, 1.0 - np.abs(load) / tauS)
        rate = (cmd - th_g) / P.TS; clip[k] = np.abs(rate) > w_eff
        th_g = th_g + P.DT * (np.clip(rate, -w_eff, w_eff) if slew else rate)
        if k % 10 == 0 and tt[k] >= t0 - 1e-9:
            keep.append((th_out[:6].copy(), q[k, :6].copy())); ks.append(k)
    ks = np.array(ks); e = np.array([P.tip(x) - P.tip(y) for x, y in keep])
    trel = tt[ks] - t0; w = np.ones_like(trel) if arm == 'B' else om(trel)
    sc = tt >= t0 - 1e-9; c = np.cumsum(clip.any(1)); flag = (c[ks] - c[np.maximum(ks - 50, 0)]) > 0
    e2 = np.sum(e ** 2, 1)
    return dict(etip=e, phi=ph[ks] - t0, omega=w, tip_mm=float(1e3 * np.sqrt(e2.mean())),
                clip_share=clip[sc].mean(0), flag_time=float(flag.mean()), flag_ms_share=float(e2[flag].sum() / e2.sum()))


if __name__ == '__main__':
    HERE = sys.argv[1]; t_start = time.time()
    tf = json.load(open(os.path.join(HERE, 'results', 'plan_intervention.json')))['tau_fix']
    r0 = run_sat('A', 1, 'B', tf['A']); r1 = P.run('A', 1, 'frozen', 'B', tau_fix=tf['A'])
    print('identity with plan_intervention.run: tip diff', abs(r0['tip_mm'] - r1['tip_mm']), flush=True)
    out = {}
    for name in ('A', 'B'):
        for slew in (True, False):
            G = {}
            for arm in ARMS:
                rr = [run_sat(name, s, arm, tf[name], slew) for s in SEEDS]
                G[arm] = rr
                out[f'{name} slew={slew} {arm}'] = dict(
                    tip_mm_median=float(np.median([r['tip_mm'] for r in rr])),
                    clip_share_by_joint=np.median([r['clip_share'] for r in rr], 0).tolist(),
                    flag_time=float(np.median([r['flag_time'] for r in rr])), flag_ms_share=float(np.median([r['flag_ms_share'] for r in rr])))
                print(name, 'slew' if slew else 'no-slew', arm, {k: (np.round(v, 3).tolist() if isinstance(v, list) else round(v, 3))
                                                                 for k, v in out[f'{name} slew={slew} {arm}'].items()}, round(time.time() - t_start), 's', flush=True)
            ebB, NB, PB = P.cell_stats([P.on_grid(r) for r in G['B']])
            for arm in ARMS:
                eb, N, Pp = P.cell_stats([P.on_grid(r) for r in G[arm]]); w2 = P.on_grid(G[arm][0], 'omega') ** 2
                m = w2[:, None] * ebB
                k = float(np.sum(eb * m) / (np.sum(m ** 2) - (NB / 5.0) * np.sum(w2 ** 2)))
                out[f'{name} slew={slew} {arm}'].update(k=k, P_share=Pp / (Pp + N), N_over_NB=N / NB)
                print(name, 'slew' if slew else 'no-slew', arm, 'k %.3f  P/(P+N) %.2f  N/NB %.3f' % (k, Pp / (Pp + N), N / NB), flush=True)
    json.dump(out, open(os.path.join(HERE, 'exploratory', 'slew_explo.json'), 'w'), indent=1)
    print('done', round((time.time() - t_start) / 60, 1), 'min')
