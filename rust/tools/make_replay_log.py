# Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
# Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
"""Run the Python engine in the servo model of examples/servo_bench.py (excerpt A, steady noise, 33 s, 12 joints)
and record every tick's inputs and output, so that another implementation can replay them and compare.

Writes replay_log.bin (the raw identification of 0.1.x) and replay_log_allan.bin (the block-mean identification,
default since 0.2.0), little-endian f64, one row per tick: q[12], read_a[12], read_b[12], cmd[12]; and
replay_w0.bin (the 12 no-load speeds). Run from the repository root:  python rust/tools/make_replay_log.py
"""
import os, sys
import numpy as np
from scipy.interpolate import CubicSpline
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, 'python'))
from syntonic_engine import SyntonicEngine

Z = np.load(os.path.join(ROOT, 'examples', 'data', 'maestro_excerpts.npz'))
DT, DEG = 1e-3, np.pi / 180
TAU_S, W0 = Z['stall_torque'], Z['no_load_speed']; ENC_SD, ENC_LSB = float(Z['enc_sd_rad']), float(Z['enc_lsb_rad'])
t0 = float(Z['A_t0']); tt = np.arange(t0 - 3.0, t0 + 30.0, DT)
q = CubicSpline(Z['A_t30'], Z['A_q30'], axis=0)(tt); tau = CubicSpline(Z['A_t30'], Z['A_tau30'], axis=0)(tt)
J = 12
for ident, fname in (('raw', 'replay_log.bin'), ('allan', 'replay_log_allan.bin')):
    rng = np.random.default_rng(5)
    eng = SyntonicEngine(q[0], DT, 0.003, W0, identification=ident)
    th_g = q[0].copy(); rows = []
    for k in range(len(tt)):
        load = tau[k]
        th_out = th_g - load * (1.0 * DEG) / TAU_S - 0.15 * DEG * np.tanh(load / (0.02 * TAU_S))
        a = np.round((th_out + rng.normal(0.0, 1.0, J) * ENC_SD) / ENC_LSB) * ENC_LSB
        b = np.round((th_out + rng.normal(0.0, 1.0, J) * ENC_SD) / ENC_LSB) * ENC_LSB
        cmd = eng.step(q[k], a, b)
        rows.append(np.concatenate([q[k], a, b, cmd]))
        w_eff = W0 * np.maximum(0.05, 1.0 - np.abs(load) / TAU_S)
        th_g = th_g + DT * np.clip((cmd - th_g) / 0.025, -w_eff, w_eff)
    np.array(rows, dtype='<f8').tofile(fname)
    print(f'wrote {fname} ({len(rows)} ticks, {ident}); final windows (ms):', np.round(eng.tau * 1e3, 2).tolist())
np.asarray(W0, dtype='<f8').tofile('replay_w0.bin')
