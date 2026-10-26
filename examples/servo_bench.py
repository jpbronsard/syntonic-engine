# Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
# Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
"""Self-contained version of the MAESTRO 1 joint bench: the syntonic engine (nothing tuned) against a PID tuned on
the rest of the same movement, on the twelve arm joints of MAESTRO 1 conducting Beethoven's Fifth (two 30 s
excerpts, A = 0:06-0:36 and B = 3:37-4:07).

Servo model (one per joint, as in the bench): command refreshed every 3 ms (333 Hz PWM) and held; internal position
loop first-order, 25 ms, speed limited to the no-load speed derated by load; 1 deg of compliance at stall torque;
0.3 deg of gear play; encoder on the joint output read twice per 1 ms tick (0.06 deg rms noise, 0.0219 deg steps).
Planned angles and load torques (inverse dynamics of the planned motion) come from examples/data.

Conditions: steady encoder noise, and a noise step (variance x10 from 10 to 20 s into the excerpt).
Scores, over the 30 s, baton arm:
- baton-tip error (the bench's primary score, fixed before any run): rms distance between the actual and the planned
  tip (mm), from the arm's kinematics (arm_fk.py, exported from the MuJoCo model), sampled at 100 Hz; and its jitter
  (rms above 5 Hz);
- joint error: rms of each joint's angle error, averaged over RJ1-RJ6 (deg).
The PID gains were chosen joint by joint to minimise the joint error (on the rest of the movement, clean noise), so
the joint score is the PID's own target; the tip score is what a conductor needs.

Run:  python examples/servo_bench.py          (about one minute)
"""
import os, sys, time
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.signal import butter, sosfiltfilt
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..', 'python')); sys.path.insert(0, HERE)
from syntonic_engine import SyntonicEngine
from arm_fk import tip

DT, FRAME = 1e-3, 3e-3
DEG = np.pi / 180
Z = np.load(os.path.join(HERE, 'data', 'maestro_excerpts.npz'))
JN = list(Z['joints']); TAU_S = Z['stall_torque']; W0 = Z['no_load_speed']
ENC_SD, ENC_LSB = float(Z['enc_sd_rad']), float(Z['enc_lsb_rad'])
TS, K_INV, HB, LIM = 0.025, (1.0 * DEG) / TAU_S, 0.5 * 0.3 * DEG, 5 * DEG


def run(engine, name, cond, seed):
    t0 = float(Z[f'{name}_t0']); tt = np.arange(t0 - 3.0, t0 + 30.0, DT)
    q = CubicSpline(Z[f'{name}_t30'], Z[f'{name}_q30'], axis=0)(tt)
    tau = CubicSpline(Z[f'{name}_t30'], Z[f'{name}_tau30'], axis=0)(tt)
    sd = np.full(len(tt), ENC_SD)
    if cond == 'noise step': sd[(tt >= t0 + 10) & (tt < t0 + 20)] *= np.sqrt(10.0)
    rng = np.random.default_rng(seed); J = q.shape[1]
    th_g = q[0].copy(); cmd = q[0].copy()
    if engine in ('syntonic', 'syntonic-allan'):
        ident = 'raw' if engine == 'syntonic' else 'allan'       # 'raw': the identification of the Oct 2026 bench
        eng = SyntonicEngine(q[0], DT, FRAME, W0, identification=ident)   # nothing tuned; floor = loop delay / TH41-A
    else:
        kp, ki, kd = (Z[f'{name}_pid'][:, i] for i in range(3))   # gains tuned on the rest of the movement
        I = np.zeros(J); ef = np.zeros(J); efp = np.zeros(J); corr = np.zeros(J)
    err2 = np.zeros(J); m = 0; keep = []
    for k in range(len(tt)):
        load = tau[k]
        th_out = th_g - load * K_INV - HB * np.tanh(load / (0.02 * TAU_S))
        a = np.round((th_out + rng.normal(0.0, 1.0, J) * sd[k]) / ENC_LSB) * ENC_LSB
        b = np.round((th_out + rng.normal(0.0, 1.0, J) * sd[k]) / ENC_LSB) * ENC_LSB
        if engine in ('syntonic', 'syntonic-allan'):
            cmd = eng.step(q[k], a, b)
        else:
            e = q[k] - 0.5 * (a + b)
            ef += (DT / 0.005) * (e - ef)
            raw = kp * e + ki * I + kd * (ef - efp) / DT
            I += np.where(np.abs(raw) >= LIM, 0.0, e * DT); efp = ef.copy()
            corr = np.clip(raw, -LIM, LIM)
            if k % 3 == 0: cmd = q[k] + corr
        w_eff = W0 * np.maximum(0.05, 1.0 - np.abs(load) / TAU_S)
        th_g = th_g + DT * np.clip((cmd - th_g) / TS, -w_eff, w_eff)
        if tt[k] >= t0: err2 += (q[k] - th_out) ** 2; m += 1
        if k % 10 == 0 and tt[k] >= t0 - 1e-9: keep.append((th_out[:6].copy(), q[k, :6].copy()))
    e = np.array([tip(a) - tip(b) for a, b in keep])                         # tip error, 100 Hz
    hp = sosfiltfilt(butter(4, 5.0, 'highpass', fs=100.0, output='sos'), e, axis=0)
    tip_rms = 1e3 * np.sqrt(np.mean(np.sum(e ** 2, 1))); jit = 1e3 * np.sqrt(np.mean(np.sum(hp ** 2, 1)))
    return tip_rms, jit, np.sqrt(err2 / m)[:6] / DEG


if __name__ == '__main__':
    t_start = time.time()
    print('baton arm; ratio = syntonic / PID (below 1: syntonic better)')
    print(f'{"":26s}{"tip error (mm)":>30s}{"tip jitter (mm)":>30s}{"joint error (deg)":>32s}')
    for name in ('A', 'B'):
        for cond, seed in (('steady noise', 5), ('noise step', 7)):
            pt, pj, pq = run('pid', name, cond, seed); st, sj, sq = run('syntonic', name, cond, seed)
            print(f'excerpt {name}, {cond:13s}   PID {pt:.3f}  syn {st:.3f}  x{st / pt:.3f}   PID {pj:.3f}  syn {sj:.3f}  x{sj / pj:.3f}'
                  f'   PID {pq.mean():.4f}  syn {sq.mean():.4f}  x{sq.mean() / pq.mean():.3f}', flush=True)
            if '--allan' in sys.argv:                                   # the 0.2.0 default identification, for comparison
                at, aj, aq = run('syntonic-allan', name, cond, seed)
                print(f'{"":26s}allan: syn {at:.3f}  x{at / pt:.3f}{"":16s}syn {aj:.3f}  x{aj / pj:.3f}{"":17s}syn {aq.mean():.4f}  x{aq.mean() / pq.mean():.3f}', flush=True)
    print(f'done in {time.time() - t_start:.0f} s')
