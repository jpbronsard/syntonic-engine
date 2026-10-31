# Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
# Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
"""Cadence bench: can the machine set its own rate from the law? Written 31 Oct 2026, before any scored run.

Origin (JP, 31 Oct 2026). In the MAESTRO 1 videos the shoulders' window sits on its floor during continuous motion, as
the window of the logged MAIDEN5 hover sat on its floor. On the wall the law decides nothing: the engine integrates over
a fixed window set by the hardware. JP's proposal: tune the machine live so that it always works inside the useful
range. The law gives the instruction itself: when tau_law falls below the floor D, the loop is too slow for that
motion (delay-limited); when tau_law is above D, noise limits the machine and the law works. Rule tested here: at each
servo frame, every joint takes the slowest frame period T whose floor stays at or below its law window,
    T_j = max{T in S : floor(T) <= tau_law_j},   the fastest T in S when none qualifies,
with floor(T) = max(T / 2 + dt / 2 + 0.1 ms, dt / 0.45), the engine's own floor (window_floor, unchanged). Nothing
is added: no margin, no hysteresis. The rule reads the law's window before it is clipped to the floor.

Plant and data: the self-contained MAESTRO 1 joint bench of examples/servo_bench.py (same servo model, same excerpts
A = 0:06-0:36 and B = 3:37-4:07 of Beethoven's Fifth, same encoder model, baton-tip kinematics from arm_fk.py), baton
arm only (RJ1-RJ6; the joints are independent in this model, so the six baton-arm joints run alone). Loop and world at
2 kHz (dt = 0.5 ms), so that the frame, not the tick, dominates the floor: floor(1, 2, 3, 5, 8 ms) = 1.11, 1.35, 1.85,
2.85, 4.35 ms. Engine: the bench engine, identification='raw' (the identification of the published joint bench),
everything else at its defaults. Frames faster than 3 ms are a simulation probe: the bench's servos accept about
330 Hz; a faster actuator link is a hardware change.

Arms (every arm sees the same encoder noise draws, tick for tick: paired):
  F3   fixed 3 ms frame (today's hardware)          F2   fixed 2 ms          F1   fixed 1 ms
  AD   adaptive, S = {1, 2, 3} ms (may only speed up from the hardware frame; starts at 3 ms)
  ADS  adaptive, S = {1, 2, 3, 5, 8} ms (may also slow down; starts at 3 ms)
Conditions: excerpts A and B; steady encoder noise, and a step (variance x10 from 10 to 20 s); seeds 1 to 5 in both
conditions: 20 runs per arm.
Scores over the 30 s: baton-tip error (rms distance actual to planned tip, sampled at 100 Hz; primary), its jitter
(above 5 Hz), joint error (RJ1-RJ6 mean rms, deg). Pooled tip score of an arm: root mean square of the 20 tip errors.
Cost: mean frame rate (frames sent per joint per second). Described: share of joint-time with tau_law below the floor
of the current frame; share of time at each frame; frame switches per joint per second.

Predictions:
P0  Gate. At 1 kHz with 12 joints and a fixed 3 ms frame this engine reproduces the syntonic tip error of
    examples/servo_bench.py in its four conditions to 1e-9 mm. If P0 fails, nothing else is scored.
P1  The wall. In F3, tau_law is below the floor (1.85 ms) for 15 to 70 % of the baton-arm joint-time under steady
    noise, in both excerpts, and less often under the noise step than under steady noise, in both excerpts.
P2  A faster frame helps: pooled tip score of F1 <= 0.95 x that of F3.
P3  Main. AD's pooled tip score is within 2 % of the best fixed arm's, with a mean frame rate at most 75 % of F1's
    (750 Hz).
P4  Slowing down costs: ADS's pooled tip score exceeds AD's by more than 5 % (the frame also holds the plan, and the
    floor rule covers the window only, not that hold).
P5  Speeding up only is safe: AD's pooled tip score <= F3's, and AD is within +1 % of F3 in each of the four cells
    (seed-pooled).
A failed line is reported as such. Before this file was committed it was run only as a smoke test (--smoke: 1.5 s of
simulation, no score printed).

Run (from the repository root; about 15 to 30 minutes):   python bench/cadence_2026-10-31/cadence_bench.py
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

DEG = np.pi / 180
Z = np.load(os.path.join(ROOT, 'examples', 'data', 'maestro_excerpts.npz'))
TAU_S = Z['stall_torque']; W0 = Z['no_load_speed']
ENC_SD, ENC_LSB = float(Z['enc_sd_rad']), float(Z['enc_lsb_rad'])
TS = 0.025
FRAMES_MS = [1, 2, 3, 5, 8]
ARMS = {'F3': ([3], False), 'F2': ([2], False), 'F1': ([1], False), 'AD': ([1, 2, 3], True), 'ADS': ([1, 2, 3, 5, 8], True)}
SEEDS = [1, 2, 3, 4, 5]


class CadenceEngine(SyntonicEngine):
    """The bench engine with a frame period per column, chosen at each frame send. allowed: (n, len(FRAMES_MS)) bool;
    adaptive: (n,) bool. Arithmetic of the window, the correction and the nominal model unchanged."""

    def __init__(self, q0, dt, w0, allowed, adaptive, start_ms=3, **kw):
        super().__init__(q0, dt, start_ms * 1e-3, w0, identification='raw', **kw)
        n = self.n
        self.F_ticks = np.array([int(round(f * 1e-3 / dt)) for f in FRAMES_MS])
        self.F_floor = np.array([window_floor(dt, f * 1e-3) for f in FRAMES_MS])
        self.allowed = np.asarray(allowed, bool); self.adaptive = np.asarray(adaptive, bool)
        i0 = FRAMES_MS.index(start_ms)
        self.fi = np.array([i0 if (self.adaptive[c] or not self.allowed[c].any() or self.allowed[c, i0])
                            else int(np.flatnonzero(self.allowed[c])[0]) for c in range(n)])
        self.floor = self.F_floor[self.fi].copy()
        self.next_send = np.zeros(n, dtype=int)
        self.law = np.full(n, np.nan)
        self.sent = np.zeros(n); self.switches = np.zeros(n)

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
        send = self.k >= self.next_send
        if send.any():
            self.cmd = np.where(send, q + self.u, self.cmd)
            self.sent += send
            pick = send & self.adaptive & ~np.isnan(self.law)
            if pick.any():
                ok = self.allowed & (self.F_floor[None, :] <= self.law[:, None])
                slowest = np.where(ok.any(1), len(FRAMES_MS) - 1 - np.argmax(ok[:, ::-1], axis=1), np.argmax(self.allowed, axis=1))
                new = np.where(pick, slowest, self.fi)
                self.switches += new != self.fi
                self.fi = new; self.floor = self.F_floor[self.fi]
            self.next_send = np.where(send, self.k + self.F_ticks[self.fi], self.next_send)
        self.th_m = self.th_m + dt * np.clip((self.cmd - self.th_m) / self.servo_tau, -self.w0, self.w0)
        self.th_ff = self.th_ff + dt * np.clip((q - self.th_ff) / self.servo_tau, -self.w0, self.w0)
        self.k += 1
        return self.cmd


def simulate(name, cond, seed, dt, joints, arms, t_end=30.0):
    """All arms in one pass, side by side as columns, sharing every encoder read. Returns per-arm scores."""
    t0 = float(Z[f'{name}_t0']); tt = np.arange(t0 - 3.0, t0 + t_end, dt)
    J = len(joints); A = len(arms)
    qJ = CubicSpline(Z[f'{name}_t30'], Z[f'{name}_q30'][:, joints], axis=0)(tt)
    tauJ = CubicSpline(Z[f'{name}_t30'], Z[f'{name}_tau30'][:, joints], axis=0)(tt)
    sd = np.full(len(tt), ENC_SD)
    if cond == 'step': sd[(tt >= t0 + 10) & (tt < t0 + 20)] *= np.sqrt(10.0)
    rng = np.random.default_rng(seed)
    tauS = np.tile(TAU_S[joints], A); w0 = np.tile(W0[joints], A)
    k_inv = (1.0 * DEG) / tauS; hb = 0.5 * 0.3 * DEG
    allowed = np.zeros((A * J, len(FRAMES_MS)), bool); adaptive = np.zeros(A * J, bool)
    for i, arm in enumerate(arms):
        fr, ad = ARMS[arm] if isinstance(arm, str) else arm
        for f in fr: allowed[i * J:(i + 1) * J, FRAMES_MS.index(f)] = True
        adaptive[i * J:(i + 1) * J] = ad
    eng = CadenceEngine(np.tile(qJ[0], A), dt, w0, allowed, adaptive)
    th_g = np.tile(qJ[0], A)
    every = int(round(0.01 / dt)); keep = []
    err2 = np.zeros(A * J); m = 0; below = np.zeros(A * J); at = np.zeros((A * J, len(FRAMES_MS)))
    sent0 = sw0 = None
    for k in range(len(tt)):
        load = np.tile(tauJ[k], A); q = np.tile(qJ[k], A)
        th_out = th_g - load * k_inv - hb * np.tanh(load / (0.02 * tauS))
        na = np.tile(rng.normal(0.0, 1.0, J), A); nb = np.tile(rng.normal(0.0, 1.0, J), A)
        ra = np.round((th_out + na * sd[k]) / ENC_LSB) * ENC_LSB
        rb = np.round((th_out + nb * sd[k]) / ENC_LSB) * ENC_LSB
        cmd = eng.step(q, ra, rb)
        w_eff = w0 * np.maximum(0.05, 1.0 - np.abs(load) / tauS)
        th_g = th_g + dt * np.clip((cmd - th_g) / TS, -w_eff, w_eff)
        if tt[k] >= t0 - 1e-9:
            if sent0 is None: sent0, sw0 = eng.sent.copy(), eng.switches.copy()
            err2 += (q - th_out) ** 2; m += 1
            below += eng.law < eng.floor
            at[np.arange(A * J), eng.fi] += 1
            if k % every == 0: keep.append((th_out.copy(), q.copy()))
    T = len(tt) * dt - 3.0
    out = {}
    for i, arm in enumerate(arms):
        s = slice(i * J, (i + 1) * J); lab = arm if isinstance(arm, str) else 'custom'
        res = dict(joint_deg=float(np.mean(np.sqrt(err2[s] / m)) / DEG),
                   rate_hz=float(np.mean(eng.sent[s] - sent0[s]) / T),
                   switches_per_s=float(np.mean(eng.switches[s] - sw0[s]) / T),
                   below_floor=float(np.mean(below[s] / m)),
                   share_at={f'{f} ms': float(np.mean(at[s, j] / m)) for j, f in enumerate(FRAMES_MS)})
        if J == 12 or joints[:6] == list(range(6)):
            e = np.array([tip(a[s][:6]) - tip(b[s][:6]) for a, b in keep])
            hp = sosfiltfilt(butter(4, 5.0, 'highpass', fs=100.0, output='sos'), e, axis=0)
            res['tip_mm'] = float(1e3 * np.sqrt(np.mean(np.sum(e ** 2, 1))))
            res['jitter_mm'] = float(1e3 * np.sqrt(np.mean(np.sum(hp ** 2, 1))))
        out[lab] = res
    return out


def gate_p0():
    """P0: 1 kHz, 12 joints, fixed 3 ms frame, against examples/servo_bench.py (the published bench)."""
    import servo_bench as SB
    rows = []
    for name in ('A', 'B'):
        for cond, sbcond, seed in (('steady', 'steady noise', 5), ('step', 'noise step', 7)):
            ref = SB.run('syntonic', name, sbcond, seed)[0]
            mine = simulate(name, cond, seed, 1e-3, list(range(12)), ['F3'])['F3']['tip_mm']
            rows.append(dict(cell=f'{name} {cond}', servo_bench_mm=ref, this_mm=mine, diff_mm=abs(ref - mine)))
            print('P0', rows[-1], flush=True)
    return rows


if __name__ == '__main__':
    sha = hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest()
    if '--smoke' in sys.argv:
        r = simulate('A', 'steady', 1, 5e-4, list(range(6)), list(ARMS), t_end=1.5)
        print('smoke ok:', sorted(r), sorted(r['F3'])); sys.exit(0)
    t_start = time.time(); OUT = os.path.join(HERE, 'results'); os.makedirs(OUT, exist_ok=True)
    summary = dict(script_sha256=sha, dt_s=5e-4, arms=list(ARMS), seeds=SEEDS)
    summary['P0'] = p0 = gate_p0()
    if max(r['diff_mm'] for r in p0) > 1e-9:
        print('P0 FAILED: nothing else is scored'); json.dump(summary, open(os.path.join(OUT, 'summary.json'), 'w'), indent=1); sys.exit(1)
    runs = {}
    for name in ('A', 'B'):
        for cond in ('steady', 'step'):
            for seed in SEEDS:
                runs[f'{name} {cond} {seed}'] = r = simulate(name, cond, seed, 5e-4, list(range(6)), list(ARMS))
                print(name, cond, seed, {a: round(v['tip_mm'], 3) for a, v in r.items()},
                      {a: round(v['rate_hz']) for a, v in r.items()}, round(time.time() - t_start), 's', flush=True)
    summary['runs'] = runs
    pooled = {a: float(np.sqrt(np.mean([runs[k][a]['tip_mm'] ** 2 for k in runs]))) for a in ARMS}
    cells = {a: {c: float(np.sqrt(np.mean([runs[k][a]['tip_mm'] ** 2 for k in runs if k.startswith(c)])))
                 for c in ('A steady', 'A step', 'B steady', 'B step')} for a in ARMS}
    rate = {a: float(np.mean([runs[k][a]['rate_hz'] for k in runs])) for a in ARMS}
    below = {a: {c: float(np.mean([runs[k][a]['below_floor'] for k in runs if k.startswith(c)]))
                 for c in ('A steady', 'A step', 'B steady', 'B step')} for a in ARMS}
    best_fixed = min(('F1', 'F2', 'F3'), key=lambda a: pooled[a])
    b3 = below['F3']
    verdict = {
        'P1': bool(all(0.15 <= b3[f'{x} steady'] <= 0.70 for x in 'AB') and all(b3[f'{x} step'] < b3[f'{x} steady'] for x in 'AB')),
        'P2': bool(pooled['F1'] <= 0.95 * pooled['F3']),
        'P3': bool(pooled['AD'] <= 1.02 * pooled[best_fixed] and rate['AD'] <= 0.75 * rate['F1']),
        'P4': bool(pooled['ADS'] > 1.05 * pooled['AD']),
        'P5': bool(pooled['AD'] <= pooled['F3'] and all(cells['AD'][c] <= 1.01 * cells['F3'][c] for c in cells['F3'])),
    }
    summary.update(pooled_tip_mm=pooled, cell_tip_mm=cells, mean_rate_hz=rate, below_floor=below, best_fixed=best_fixed,
                   verdict=verdict, minutes=round((time.time() - t_start) / 60, 1))
    json.dump(summary, open(os.path.join(OUT, 'summary.json'), 'w'), indent=1)
    print('pooled tip (mm)', {a: round(v, 3) for a, v in pooled.items()})
    print('mean frame rate (Hz)', {a: round(v) for a, v in rate.items()})
    print('below floor', {a: {c: round(v, 2) for c, v in d.items()} for a, d in below.items()})
    print('verdict', verdict, f'({summary["minutes"]} min)')
