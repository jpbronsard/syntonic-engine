# Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
# Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
"""Exploratory, after slew_explo.py found no joint at its slew limit (go from JP). Reading: the engine's correction
(u' = e / tau) and the servo (time constant TS) form a second-order loop, TS tau e'' + tau e' + e = TS tau qddot (plus
the load's deflection), with static gain TS tau, natural frequency 1 / sqrt(TS tau) and damping 0.5 sqrt(tau / TS),
about 0.17 for a 3 ms window: lightly damped. The plan's fast content is then not seen through a static acceleration
gain, and its share changes with the tempo. Checks, frozen window (tau of results/plan_intervention.json), arms B,
U70, U50, excerpts A and B: (1) the published plant with the encoder noise set to zero gives the plan-locked error
itself (quantization, frame hold and load kept); its k (no noise correction needed); (2) the linear loop above, driven
by the same plan and load, no hold, no quantization; its k and the share of the noise-free error it reproduces at the
tip; (3) the static model e = TS tau qddot (k = 1 by construction). Usage: python loop_dynamics_explo.py <bench dir>"""
import sys, os, json
import numpy as np
from scipy.interpolate import CubicSpline
sys.path.insert(0, sys.argv[1])
import plan_intervention as P

J = 6; ARMS = ('B', 'U70', 'U50')


def plan(name, arm):
    t0 = float(P.Z[f'{name}_t0']); T, phi, om, omd = P.tempo(arm)
    tt = np.arange(t0 - P.PRE, t0 + T, P.DT); ph = tt if arm == 'B' else t0 + phi(tt - t0)
    q = CubicSpline(P.Z[f'{name}_t30'], P.Z[f'{name}_q30'][:, :J], axis=0)(ph)
    load = CubicSpline(P.Z[f'{name}_t30'], P.Z[f'{name}_tau30'][:, :J], axis=0)(ph)
    w = np.ones_like(tt) if arm == 'B' else om(tt - t0)
    return t0, tt, ph, q, load, w


def bench_noise_free(name, arm, tf):
    t0, tt, ph, q, load_all, w = plan(name, arm)
    tauS, w0 = P.TAU_S[:J], P.W0[:J]; k_inv = (1.0 * P.DEG) / tauS; hb = 0.5 * 0.3 * P.DEG
    th_g = q[0].copy(); allowed = np.zeros((J, 5), bool); allowed[:, 2] = True
    eng = P.CadenceEngine(q[0], P.DT, w0, allowed, np.zeros(J, bool))
    tf = np.asarray(tf, float); eng.floor = tf.copy(); eng.tau_max = tf.copy(); eng.tau = tf.copy(); eng.tau_target = tf.copy()
    keep, ks = [], []
    for k in range(len(tt)):
        load = load_all[k]
        th_out = th_g - load * k_inv - hb * np.tanh(load / (0.02 * tauS))
        a = np.round(th_out / P.ENC_LSB) * P.ENC_LSB
        cmd = eng.step(q[k], a, a)
        w_eff = w0 * np.maximum(0.05, 1.0 - np.abs(load) / tauS)
        th_g = th_g + P.DT * np.clip((cmd - th_g) / P.TS, -w_eff, w_eff)
        if k % 10 == 0 and tt[k] >= t0 - 1e-9: keep.append((th_out.copy(), q[k].copy())); ks.append(k)
    return grid(name, arm, keep, ks, tt, ph, t0, w)


def linear_loop(name, arm, tf):
    t0, tt, ph, q, load_all, w = plan(name, arm)
    tauS = P.TAU_S[:J]; k_inv = (1.0 * P.DEG) / tauS; hb = 0.5 * 0.3 * P.DEG
    tf = np.asarray(tf, float); g = 1.0 - np.exp(-P.DT / tf)
    th_g = q[0].copy(); u = np.zeros(J); keep, ks = [], []
    for k in range(len(tt)):
        load = load_all[k]
        th_out = th_g - load * k_inv - hb * np.tanh(load / (0.02 * tauS))
        u = u + g * (q[k] - th_out)
        th_g = th_g + P.DT * (q[k] + u - th_g) / P.TS
        if k % 10 == 0 and tt[k] >= t0 - 1e-9: keep.append((th_out.copy(), q[k].copy())); ks.append(k)
    return grid(name, arm, keep, ks, tt, ph, t0, w)


def grid(name, arm, keep, ks, tt, ph, t0, w):
    e = np.array([P.tip(x) - P.tip(y) for x, y in keep]); ks = np.array(ks)
    r = dict(etip=e, phi=ph[ks] - t0, omega=w[ks])
    return P.on_grid(r), P.on_grid(r, 'omega')


def kcoef(eI, eB, wI):
    m = (wI ** 2)[:, None] * eB
    return float(np.sum(eI * m) / np.sum(m ** 2))


if __name__ == '__main__':
    HERE = sys.argv[1]; tfs = json.load(open(os.path.join(HERE, 'results', 'plan_intervention.json')))['tau_fix']
    out = {}
    for name in ('A', 'B'):
        tf = np.array(tfs[name]); zeta = 0.5 * np.sqrt(tf / P.TS); fn = 1.0 / (2 * np.pi * np.sqrt(P.TS * tf))
        print(name, 'window ms', np.round(tf * 1e3, 2), 'damping', np.round(zeta, 2), 'natural Hz', np.round(fn, 1), flush=True)
        nf = {a: bench_noise_free(name, a, tf) for a in ARMS}; lin = {a: linear_loop(name, a, tf) for a in ARMS}
        for a in ARMS:
            eT, wI = nf[a]; eL, _ = lin[a]
            share = 1.0 - np.sum((eT - eL) ** 2) / np.sum(eT ** 2)
            out[f'{name} {a}'] = dict(k_noise_free=kcoef(eT, nf['B'][0], wI), k_linear=kcoef(eL, lin['B'][0], wI),
                                      linear_share_of_noise_free=float(share), tip_rms_noise_free_mm=float(1e3 * np.sqrt(np.mean(np.sum(eT ** 2, 1)))),
                                      tip_rms_linear_mm=float(1e3 * np.sqrt(np.mean(np.sum(eL ** 2, 1)))))
            print(name, a, {k: round(v, 3) for k, v in out[f'{name} {a}'].items()}, flush=True)
        out[f'{name} window_ms'] = (tf * 1e3).tolist(); out[f'{name} damping'] = zeta.tolist(); out[f'{name} natural_Hz'] = fn.tolist()
    json.dump(out, open(os.path.join(HERE, 'exploratory', 'loop_dynamics_explo.json'), 'w'), indent=1)
    print('done')
