# Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
# Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
"""Cadence bench, exploratory diagnostic: where does the gain of a faster servo frame come from? Written 31 Oct 2026
after the cadence bench's verdict (P3 failed), before any run of this file. Exploratory: no prediction is judged here,
and nothing in it changes the cadence bench's verdict.

A servo frame does two things at once. It bounds the window (the floor, T / 2 + dt / 2 + 0.1 ms) and it holds the
command between two sends, and the command carries both the plan and the correction. The cadence rule covers only the
first. This diagnostic separates them, in simulation only, by giving each column its own floor, its own send period
for the command, and its own refresh period for the correction inside the command:
  A   send 3 ms, floor 1.85 ms                                     (= the cadence bench's 3 ms arm)
  B   send 1 ms, floor 1.11 ms                                     (= the cadence bench's 1 ms arm)
  C   send 3 ms, floor 1.11 ms   (window freed, command still held 3 ms; not realizable, a probe)
  D   send 1 ms, floor 1.85 ms   (command every 1 ms, window bound as at 3 ms)
  E   send 1 ms, floor 1.85 ms, correction refreshed every 3 ms   (only the plan is sent faster)
Same plant, excerpts, noise, seeds, scores and engine as cadence_bench.py (2 kHz, baton arm, 20 paired runs).
Gate: A and B reproduce the 3 ms and 1 ms arms of results/summary.json run by run to 1e-9 mm, or nothing is read.
Readings (descriptive, pooled tip error): the window's share of the 3 ms -> 1 ms gain, A - C (command held 3 ms) and
D - B (command every 1 ms); the hold's share, A - D (floor 1.85) and C - B (floor 1.11); within the hold at floor
1.85, the plan's part A - E and the correction's part E - D. If the hold carries most of the gain, a cadence rule must
account for the plan's motion within a frame, not only for the floor.

Run (from the repository root, about 10 minutes):   python bench/cadence_2026-10-31/diag_hold_window.py
"""
import os, sys, json, time, hashlib
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.signal import butter, sosfiltfilt
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, 'python')); sys.path.insert(0, os.path.join(ROOT, 'examples'))
from syntonic_engine import SyntonicEngine
from syntonic_engine.engine import kappa_H, window_floor
from arm_fk import tip

DEG = np.pi / 180; DT = 5e-4
Z = np.load(os.path.join(ROOT, 'examples', 'data', 'maestro_excerpts.npz'))
TAU_S = Z['stall_torque']; W0 = Z['no_load_speed']
ENC_SD, ENC_LSB = float(Z['enc_sd_rad']), float(Z['enc_lsb_rad'])
TS = 0.025
F1, F3 = window_floor(DT, 1e-3), window_floor(DT, 3e-3)
COLS = {'A': (3, F3, 3), 'B': (1, F1, 1), 'C': (3, F1, 3), 'D': (1, F3, 1), 'E': (1, F3, 3)}   # send ms, floor s, correction ms
SEEDS = [1, 2, 3, 4, 5]


class DiagEngine(SyntonicEngine):
    """The bench engine with, per column, its own floor, command send period and correction refresh period."""

    def __init__(self, q0, dt, w0, send_ms, floor_s, corr_ms):
        super().__init__(q0, dt, 3e-3, w0, identification='raw')
        self.floor = np.asarray(floor_s, float)
        self.send = np.array([int(round(x * 1e-3 / dt)) for x in send_ms])
        self.corr = np.array([int(round(x * 1e-3 / dt)) for x in corr_ms])
        self.u_lat = np.zeros(self.n); self.law = np.full(self.n, np.nan)

    def step(self, q_plan, read_a, read_b):
        dt = self.dt
        q = np.asarray(q_plan, dtype=float)
        a = np.asarray(read_a, dtype=float); b = np.asarray(read_b, dtype=float)
        meas = 0.5 * (a + b); e = q - meas
        d = (a - b) / np.sqrt(2.0)
        self.n2 += self.a_w * (d * d - self.n2)
        r = (q - self.th_ff) + (meas - self.th_m)
        fit = self._identify_raw(r)
        if fit is not None:
            Hn, lam, valid = fit
            s2 = 0.5 * self.n2 * dt
            ts = kappa_H(Hn) * (np.maximum(s2, 1e-30) / np.maximum(lam, 1e-30)) ** (1.0 / (2 * Hn + 1))
            self.law = np.where(valid, np.minimum(ts, self.tau_max), self.law)
            ts = np.clip(ts, self.floor, self.tau_max)
            self.tau_target = np.where(valid, ts, self.tau_target); self.H = np.where(valid, Hn, self.H)
            rr = 1.0 / self.tau; dr_max = self.tau_rate * dt / self.tau ** 2
            rr = rr + np.clip(1.0 / self.tau_target - rr, -dr_max, dr_max); self.tau = 1.0 / rr
        g = 1.0 - np.exp(-dt / self.tau)
        self.u = np.clip(self.u + g * e, -self.limit, self.limit)
        self.u_lat = np.where(self.k % self.corr == 0, self.u, self.u_lat)
        self.cmd = np.where(self.k % self.send == 0, q + self.u_lat, self.cmd)
        self.th_m = self.th_m + dt * np.clip((self.cmd - self.th_m) / self.servo_tau, -self.w0, self.w0)
        self.th_ff = self.th_ff + dt * np.clip((q - self.th_ff) / self.servo_tau, -self.w0, self.w0)
        self.k += 1
        return self.cmd


def simulate(name, cond, seed, cols=tuple(COLS), t_end=30.0):
    t0 = float(Z[f'{name}_t0']); tt = np.arange(t0 - 3.0, t0 + t_end, DT)
    joints = list(range(6)); J = 6; A = len(cols)
    qJ = CubicSpline(Z[f'{name}_t30'], Z[f'{name}_q30'][:, joints], axis=0)(tt)
    tauJ = CubicSpline(Z[f'{name}_t30'], Z[f'{name}_tau30'][:, joints], axis=0)(tt)
    sd = np.full(len(tt), ENC_SD)
    if cond == 'step': sd[(tt >= t0 + 10) & (tt < t0 + 20)] *= np.sqrt(10.0)
    rng = np.random.default_rng(seed)
    tauS = np.tile(TAU_S[joints], A); w0 = np.tile(W0[joints], A)
    k_inv = (1.0 * DEG) / tauS; hb = 0.5 * 0.3 * DEG
    send = np.repeat([COLS[c][0] for c in cols], J); floor = np.repeat([COLS[c][1] for c in cols], J)
    corr = np.repeat([COLS[c][2] for c in cols], J)
    eng = DiagEngine(np.tile(qJ[0], A), DT, w0, send, floor, corr)
    th_g = np.tile(qJ[0], A); every = int(round(0.01 / DT)); keep = []; below = np.zeros(A * J); m = 0
    for k in range(len(tt)):
        load = np.tile(tauJ[k], A); q = np.tile(qJ[k], A)
        th_out = th_g - load * k_inv - hb * np.tanh(load / (0.02 * tauS))
        na = np.tile(rng.normal(0.0, 1.0, J), A); nb = np.tile(rng.normal(0.0, 1.0, J), A)
        ra = np.round((th_out + na * sd[k]) / ENC_LSB) * ENC_LSB
        rb = np.round((th_out + nb * sd[k]) / ENC_LSB) * ENC_LSB
        cmd = eng.step(q, ra, rb)
        w_eff = w0 * np.maximum(0.05, 1.0 - np.abs(load) / tauS)
        th_g = th_g + DT * np.clip((cmd - th_g) / TS, -w_eff, w_eff)
        if tt[k] >= t0 - 1e-9:
            below += eng.law < eng.floor; m += 1
            if k % every == 0: keep.append((th_out.copy(), q.copy()))
    out = {}
    for i, c in enumerate(cols):
        s = slice(i * J, (i + 1) * J)
        e = np.array([tip(a[s]) - tip(b[s]) for a, b in keep])
        hp = sosfiltfilt(butter(4, 5.0, 'highpass', fs=100.0, output='sos'), e, axis=0)
        out[c] = dict(tip_mm=float(1e3 * np.sqrt(np.mean(np.sum(e ** 2, 1)))),
                      jitter_mm=float(1e3 * np.sqrt(np.mean(np.sum(hp ** 2, 1)))), below_floor=float(np.mean(below[s] / m)))
    return out


if __name__ == '__main__':
    sha = hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest()
    ref = json.load(open(os.path.join(HERE, 'results', 'summary.json')))['runs']
    t_start = time.time(); runs = {}; gate = 0.0
    for name in ('A', 'B'):
        for cond in ('steady', 'step'):
            for seed in SEEDS:
                key = f'{name} {cond} {seed}'
                runs[key] = r = simulate(name, cond, seed)
                gate = max(gate, abs(r['A']['tip_mm'] - ref[key]['F3']['tip_mm']), abs(r['B']['tip_mm'] - ref[key]['F1']['tip_mm']))
                print(key, {c: round(v['tip_mm'], 3) for c, v in r.items()}, 'gate', gate, round(time.time() - t_start), 's', flush=True)
    P = {c: float(np.sqrt(np.mean([runs[k][c]['tip_mm'] ** 2 for k in runs]))) for c in COLS}
    total = P['A'] - P['B']
    parts = {'window, command held 3 ms (A - C)': P['A'] - P['C'], 'window, command every 1 ms (D - B)': P['D'] - P['B'],
             'hold, floor 1.85 ms (A - D)': P['A'] - P['D'], 'hold, floor 1.11 ms (C - B)': P['C'] - P['B'],
             'hold: plan part (A - E)': P['A'] - P['E'], 'hold: correction part (E - D)': P['E'] - P['D']}
    out = dict(script_sha256=sha, gate_max_diff_mm=gate, gate_ok=gate <= 1e-9, pooled_tip_mm=P, total_gain_mm=total,
               parts_mm=parts, parts_share_of_total={k: v / total for k, v in parts.items()}, runs=runs,
               minutes=round((time.time() - t_start) / 60, 1))
    json.dump(out, open(os.path.join(HERE, 'results', 'diag_hold_window.json'), 'w'), indent=1)
    print('gate', 'ok' if gate <= 1e-9 else 'FAILED', gate)
    print('pooled tip (mm)', {c: round(v, 3) for c, v in P.items()}, 'total gain', round(total, 3))
    for k, v in parts.items(): print(f'  {k}: {v:.3f} mm ({v / total:.0%})')
