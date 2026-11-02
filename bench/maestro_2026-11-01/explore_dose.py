# Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
# Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
"""Exploratory (not pre-registered), 1 Nov 2026: does the effect grow with how much the conditions change?
F15D setting (window_transport.TransportEngine: engine 0.2.8, 2 kHz loop, 1.5 ms frames, lead (servo_tau + D)).
Arms: BASE (published transport), E64 / N64 (1/tau / tau relaxing, 64 ms), FIX5 (window fixed at 5 ms on every joint,
the best single fixed value found on seed 1 under steady and x10-step noise).
Noise profiles (encoder read noise sd multiplier over time, scored 30 s):
  steady   x1 throughout
  step10   variance x10 from 10 to 20 s (the bench's step)
  step100  variance x100 from 10 to 20 s
  bursts   variance x30 for 0.5 s every 2.5 s from 5 to 25 s
Usage: python explore_dose.py <repo_root> <profile> <excerpt> <seed...>"""
import os, sys, json, time
import numpy as np
from scipy.interpolate import CubicSpline
ROOT = os.path.abspath(sys.argv[1]); PROF = sys.argv[2]; NAME = sys.argv[3]; SEEDS = [int(s) for s in sys.argv[4:]]
HERE_W = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE_W)
sys.argv = [sys.argv[0], ROOT, 'import']
import window_transport as WT
H = WT.H
OUT = os.path.join(HERE_W, 'results_dose'); os.makedirs(OUT, exist_ok=True)
ARMS = {'BASE': (0, 0.0), 'E64': (1, 0.064), 'N64': (2, 0.064), 'FIX5': (0, 0.0)}


def profile(trel):
    m = np.ones_like(trel)
    if PROF == 'step10': m[(trel >= 10) & (trel < 20)] = np.sqrt(10.0)
    if PROF == 'step100': m[(trel >= 10) & (trel < 20)] = np.sqrt(100.0)
    if PROF == 'bursts':
        ph = (trel - 5.0) % 2.5
        m[(trel >= 5) & (trel < 25) & (ph < 0.5)] = np.sqrt(30.0)
    return m


def simulate(seed, arms, t_end=30.0):
    dt = 5e-4; joints = list(range(6)); J = 6; A = len(arms)
    t0 = float(H.Z[f'{NAME}_t0']); tt = np.arange(t0 - 3.0, t0 + t_end, dt)
    S = CubicSpline(H.Z[f'{NAME}_t30'], H.Z[f'{NAME}_q30'][:, joints], axis=0); qJ = S(tt); qdJ = S(tt, 1)
    tauJ = CubicSpline(H.Z[f'{NAME}_t30'], H.Z[f'{NAME}_tau30'][:, joints], axis=0)(tt)
    sd = H.ENC_SD * profile(tt - t0)
    rng = np.random.default_rng(seed)
    tauS = np.tile(H.TAU_S[joints], A); w0 = np.tile(H.W0[joints], A)
    k_inv = (1.0 * H.DEG) / tauS; hb = 0.5 * 0.3 * H.DEG
    allowed = np.zeros((A * J, len(H.FR)), bool); allowed[:, H.FR.index(1.5)] = True
    tr = np.repeat([ARMS[a][0] for a in arms], J); T = np.repeat([ARMS[a][1] for a in arms], J)
    eng = WT.TransportEngine(np.tile(qJ[0], A), dt, w0, allowed, np.zeros(A * J, int), np.full(A * J, 2), transport=tr, T=T)
    eng.tau_max = np.full(A * J, float(eng.tau_max))
    for i, lab in enumerate(arms):
        if lab == 'FIX5':
            s = slice(i * J, (i + 1) * J)
            eng.floor = np.asarray(eng.floor, float).copy(); eng.floor[s] = 0.005; eng.tau_max[s] = 0.005
            eng.tau[s] = 0.005; eng.tau_target[s] = 0.005
    th_g = np.tile(qJ[0], A); every = int(round(0.01 / dt)); keep = []
    for k in range(len(tt)):
        load = np.tile(tauJ[k], A); q = np.tile(qJ[k], A); qd = np.tile(qdJ[k], A)
        th_out = th_g - load * k_inv - hb * np.tanh(load / (0.02 * tauS))
        na = np.tile(rng.normal(0.0, 1.0, J), A); nb = np.tile(rng.normal(0.0, 1.0, J), A)
        ra = np.round((th_out + na * sd[k]) / H.ENC_LSB) * H.ENC_LSB
        rb = np.round((th_out + nb * sd[k]) / H.ENC_LSB) * H.ENC_LSB
        cmd = eng.step(q, qd, ra, rb)
        w_eff = w0 * np.maximum(0.05, 1.0 - np.abs(load) / tauS)
        th_g = th_g + dt * np.clip((cmd - th_g) / H.TS, -w_eff, w_eff)
        if tt[k] >= t0 - 1e-9 and k % every == 0: keep.append((th_out.copy(), q.copy()))
    out = {}
    for i, lab in enumerate(arms):
        s = slice(i * J, (i + 1) * J)
        e = np.array([H.tip(a[s]) - H.tip(b[s]) for a, b in keep])
        out[lab] = float(1e3 * np.sqrt(np.mean(np.sum(e ** 2, 1))))
    return out


for seed in SEEDS:
    t = time.time(); r = simulate(seed, list(ARMS))
    with open(os.path.join(OUT, 'runs.jsonl'), 'a') as f: f.write(json.dumps(dict(profile=PROF, excerpt=NAME, seed=seed, **r)) + '\n')
    print(PROF, NAME, seed, {a: round(v, 3) for a, v in r.items()}, round(time.time() - t), 's', flush=True)
