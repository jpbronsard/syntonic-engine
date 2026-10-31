# Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
# Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
"""Cadence along the e-geodesic. Written 31 Oct 2026, before any run of this file; idea and go from JP: "Suivre
l'e-geodesique remplacerait l'escalier par un chemin continu, lisse et sans saut, porte par la meme geometrie que la
fenetre."

Geometry (Dual Geodesy, Theorem 1): the e-geodesic between two canonical windows is the rate-linear path, the
forgetting rate 1/tau affine along it. The engine's window already moves on that path at a bounded speed (tau_rate,
0.4 s per second). A frame cadence 1/T is a rate as well. In the cadence bench (AD) the frame was chosen as a staircase,
{1, 2, 3} ms, from the law's raw window. Here the cadence follows a copy of the law's window carried on the e-geodesic:
  tau_e  moves toward the law's raw window (bounded above by tau_max, not by the floor) along the rate-linear path,
         1/tau_e stepped toward 1/tau_law by at most tau_rate dt / tau_e^2 per tick: the engine's own window dynamics,
         without the floor (the floor depends on the frame, the frame on tau_e);
  T      at each send, the slowest frame on the loop's finest grid (0.5 ms at 2 kHz: 1, 1.5, 2, 2.5, 3 ms) whose
         floor stays at or below tau_e: the inverse of the engine's floor rule, rounded down to the grid.
No new constant: the speed is the engine's, the map is its floor rule, the bounds are the actuator's.

Arms, one pass, same 20 paired runs as cadence_bench.py (2 kHz, baton arm, excerpts A and B, steady and x10 step,
seeds 1 to 5): F3, F1 fixed (anchors); AD, the staircase over {1, 2, 3} ms from the raw law (cadence_bench.py's rule);
GR, the 0.5 ms grid from the raw law (separates the finer grid from the path); EG, the 0.5 ms grid from tau_e.
Second pass (the timing control of timing_control.py applied to EG): EG's own schedule replayed with the scored 30 s
shifted by 15 s (EGsh), and unshifted (EGr, a check).
Gates: F3, F1 and AD reproduce results/summary.json run by run to 1e-6 mm; EGr reproduces EG to 1e-6 mm.
Efficiency against the fixed frames: frontier(rate) = the pooled tip error of the fixed 3, 2 and 1 ms arms of
results/summary.json, linear in the frame rate between them (333, 500, 1000 Hz); gain(arm) = 1 - pooled / frontier
at the arm's mean rate. Known before this run, from today's results: AD's gain is 7.2 % (computed after its verdict).

Predictions (judged):
E1  Timing holds on the path: EG's pooled tip error <= 0.95 x EGsh's, and EG below EGsh in at least 18 of 20 runs.
E2  Efficiency holds: gain(EG) >= 5 %.
E3  The path neither adds nor costs: |gain(EG) - gain(AD)| <= 1 percentage point. If EG exceeds AD by more, the
    geometry adds something to the cadence; if it falls below by more, the staircase was better.
Described, not judged: GR against EG and AD (grid against path); mean rates; switches per joint per second and the
mean size of a switch (ms); time on the floor.
Before this file was committed it was run only as a smoke test (1.5 s of excerpt A; no arm compared).

Run (from the repository root, about 15 to 20 minutes):   python bench/cadence_2026-10-31/egeodesic_cadence.py
"""
import os, sys, json, time, hashlib
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import timing_control as TCM
from cadence_bench import SEEDS, Z, TAU_S, W0, ENC_SD, ENC_LSB, TS, DEG, tip
from syntonic_engine import SyntonicEngine
from syntonic_engine.engine import kappa_H, window_floor
from scipy.interpolate import CubicSpline
from scipy.signal import butter, sosfiltfilt

DT = 5e-4; SHIFT_S = 15.0
FR = [1.0, 1.5, 2.0, 2.5, 3.0]
ARMS = {'F3': ([3.0], 0), 'F1': ([1.0], 0), 'AD': ([1.0, 2.0, 3.0], 1), 'GR': (FR, 1), 'EG': (FR, 2)}   # frames, source


class EGEngine(SyntonicEngine):
    """The bench engine; per column a set of allowed frames and a source for the choice: 0 fixed, 1 the raw law,
    2 the law carried on the e-geodesic (tau_e). Window, correction and nominal model unchanged."""

    def __init__(self, q0, dt, w0, allowed, source):
        super().__init__(q0, dt, 3e-3, w0, identification='raw')
        n = self.n
        self.Fms = np.array(FR)
        self.F_ticks = np.array([int(round(f * 1e-3 / dt)) for f in FR])
        self.F_floor = np.array([window_floor(dt, f * 1e-3) for f in FR])
        self.allowed = np.asarray(allowed, bool); self.source = np.asarray(source, int)
        i3 = FR.index(3.0)
        self.fi = np.array([i3 if (self.source[c] > 0 or self.allowed[c, i3]) else int(np.flatnonzero(self.allowed[c])[0])
                            for c in range(n)])
        self.floor = self.F_floor[self.fi].copy()
        self.next_send = np.zeros(n, dtype=int)
        self.law = np.full(n, np.nan); self.tau_e = np.full(n, self.tau_max)
        self.sent = np.zeros(n); self.switches = np.zeros(n); self.jump_ms = np.zeros(n)

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
            raw = np.minimum(ts, self.tau_max)
            self.law = np.where(valid, raw, self.law)
            re = 1.0 / self.tau_e; dre = self.tau_rate * dt / self.tau_e ** 2          # e-geodesic: rate-linear step
            re = re + np.clip(1.0 / raw - re, -dre, dre)
            self.tau_e = np.where(valid, 1.0 / re, self.tau_e)
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
            src = np.where(self.source == 2, self.tau_e, self.law)
            pick = send & (self.source > 0) & ~np.isnan(src)
            if pick.any():
                ok = self.allowed & (self.F_floor[None, :] <= src[:, None])
                slowest = np.where(ok.any(1), len(FR) - 1 - np.argmax(ok[:, ::-1], axis=1), np.argmax(self.allowed, axis=1))
                new = np.where(pick, slowest, self.fi)
                self.switches += new != self.fi; self.jump_ms += np.abs(self.Fms[new] - self.Fms[self.fi])
                self.fi = new; self.floor = self.F_floor[self.fi]
            self.next_send = np.where(send, self.k + self.F_ticks[self.fi], self.next_send)
        self.th_m = self.th_m + dt * np.clip((self.cmd - self.th_m) / self.servo_tau, -self.w0, self.w0)
        self.th_ff = self.th_ff + dt * np.clip((q - self.th_ff) / self.servo_tau, -self.w0, self.w0)
        self.k += 1
        return self.cmd


def live_run(name, cond, seed, t_end=30.0, record='EG'):
    t0 = float(Z[f'{name}_t0']); tt = np.arange(t0 - 3.0, t0 + t_end, DT)
    joints = list(range(6)); J = 6; labels = list(ARMS); A = len(labels)
    qJ = CubicSpline(Z[f'{name}_t30'], Z[f'{name}_q30'][:, joints], axis=0)(tt)
    tauJ = CubicSpline(Z[f'{name}_t30'], Z[f'{name}_tau30'][:, joints], axis=0)(tt)
    sd = np.full(len(tt), ENC_SD)
    if cond == 'step': sd[(tt >= t0 + 10) & (tt < t0 + 20)] *= np.sqrt(10.0)
    rng = np.random.default_rng(seed)
    tauS = np.tile(TAU_S[joints], A); w0 = np.tile(W0[joints], A)
    k_inv = (1.0 * DEG) / tauS; hb = 0.5 * 0.3 * DEG
    allowed = np.zeros((A * J, len(FR)), bool); source = np.zeros(A * J, int)
    for i, lab in enumerate(labels):
        fr, src = ARMS[lab]
        for f in fr: allowed[i * J:(i + 1) * J, FR.index(f)] = True
        source[i * J:(i + 1) * J] = src
    eng = EGEngine(np.tile(qJ[0], A), DT, w0, allowed, source)
    th_g = np.tile(qJ[0], A); every = int(round(0.01 / DT)); keep = []
    below = np.zeros(A * J); m = 0; sent0 = sw0 = jm0 = None
    ir = labels.index(record); cs = slice(ir * J, (ir + 1) * J); rec_f, rec_s = [], []
    for k in range(len(tt)):
        load = np.tile(tauJ[k], A); q = np.tile(qJ[k], A)
        th_out = th_g - load * k_inv - hb * np.tanh(load / (0.02 * tauS))
        na = np.tile(rng.normal(0.0, 1.0, J), A); nb = np.tile(rng.normal(0.0, 1.0, J), A)
        ra = np.round((th_out + na * sd[k]) / ENC_LSB) * ENC_LSB
        rb = np.round((th_out + nb * sd[k]) / ENC_LSB) * ENC_LSB
        fl = eng.floor[cs].copy(); s0 = eng.sent[cs].copy()
        cmd = eng.step(q, ra, rb)
        rec_f.append(fl); rec_s.append(eng.sent[cs] != s0)
        w_eff = w0 * np.maximum(0.05, 1.0 - np.abs(load) / tauS)
        th_g = th_g + DT * np.clip((cmd - th_g) / TS, -w_eff, w_eff)
        if tt[k] >= t0 - 1e-9:
            if sent0 is None: sent0, sw0, jm0 = eng.sent.copy(), eng.switches.copy(), eng.jump_ms.copy()
            below += eng.law < eng.floor; m += 1
            if k % every == 0: keep.append((th_out.copy(), q.copy()))
    T = len(tt) * DT - 3.0; out = {}
    for i, lab in enumerate(labels):
        s = slice(i * J, (i + 1) * J)
        e = np.array([tip(a[s]) - tip(b[s]) for a, b in keep])
        hp = sosfiltfilt(butter(4, 5.0, 'highpass', fs=100.0, output='sos'), e, axis=0)
        nsw = float(np.sum(eng.switches[s] - sw0[s]))
        out[lab] = dict(tip_mm=float(1e3 * np.sqrt(np.mean(np.sum(e ** 2, 1)))),
                        jitter_mm=float(1e3 * np.sqrt(np.mean(np.sum(hp ** 2, 1)))),
                        rate_hz=float(np.mean(eng.sent[s] - sent0[s]) / T),
                        switches_per_s=float(np.mean(eng.switches[s] - sw0[s]) / T),
                        mean_switch_ms=float(np.sum(eng.jump_ms[s] - jm0[s]) / nsw) if nsw > 0 else 0.0,
                        below_floor=float(np.mean(below[s] / m)))
    return out, np.array(rec_f), np.array(rec_s)


def frontier(rate, fixed):
    """Pooled tip error of the fixed frames, linear in the frame rate between 333, 500 and 1000 Hz."""
    return float(np.interp(rate, [1000 / 3, 500.0, 1000.0], [fixed['F3'], fixed['F2'], fixed['F1']]))


if __name__ == '__main__':
    sha = hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest()
    S0 = json.load(open(os.path.join(HERE, 'results', 'summary.json'))); ref = S0['runs']; fixed = S0['pooled_tip_mm']
    t_start = time.time(); runs = {}; g1 = g2 = 0.0
    for name in ('A', 'B'):
        for cond in ('steady', 'step'):
            for seed in SEEDS:
                key = f'{name} {cond} {seed}'
                live, fl, sn = live_run(name, cond, seed)
                g1 = max([g1] + [abs(live[a]['tip_mm'] - ref[key][a]['tip_mm']) for a in ('F3', 'F1', 'AD')])
                n0 = int(round(3.0 / DT)); N = fl.shape[0] - n0; S = int(round(SHIFT_S / DT))
                idx = np.arange(fl.shape[0]); idx[n0:] = n0 + (np.arange(N) + S) % N
                imp = TCM.imposed_run(name, cond, seed, np.concatenate([fl, fl[idx]], 1), np.concatenate([sn, sn[idx]], 1),
                                      ['EGr', 'EGsh'])
                g2 = max(g2, abs(imp['EGr']['tip_mm'] - live['EG']['tip_mm']))
                runs[key] = dict(**live, **imp)
                print(key, {a: round(runs[key][a]['tip_mm'], 3) for a in ('F3', 'F1', 'AD', 'GR', 'EG', 'EGsh')},
                      {a: round(runs[key][a]['rate_hz']) for a in ('AD', 'GR', 'EG')}, 'g1', g1, 'g2', g2,
                      round(time.time() - t_start), 's', flush=True)
    arms = ('F3', 'F1', 'AD', 'GR', 'EG', 'EGr', 'EGsh')
    P = {a: float(np.sqrt(np.mean([runs[k][a]['tip_mm'] ** 2 for k in runs]))) for a in arms}
    R = {a: float(np.mean([runs[k][a]['rate_hz'] for k in runs])) for a in arms}
    gain = {a: 1.0 - P[a] / frontier(R[a], fixed) for a in ('AD', 'GR', 'EG')}
    wins = int(sum(runs[k]['EG']['tip_mm'] < runs[k]['EGsh']['tip_mm'] for k in runs))
    gates = dict(G1_max_diff_mm=g1, G1_ok=g1 <= 1e-6, G2_max_diff_mm=g2, G2_ok=g2 <= 1e-6)
    ok = gates['G1_ok'] and gates['G2_ok']
    verdict = dict(E1=bool(ok and P['EG'] <= 0.95 * P['EGsh'] and wins >= 18), E2=bool(ok and gain['EG'] >= 0.05),
                   E3=bool(ok and abs(gain['EG'] - gain['AD']) <= 0.01))
    desc = {a: dict(switches_per_s=float(np.mean([runs[k][a]['switches_per_s'] for k in runs])),
                    mean_switch_ms=float(np.mean([runs[k][a]['mean_switch_ms'] for k in runs])),
                    below_floor=float(np.mean([runs[k][a]['below_floor'] for k in runs]))) for a in ('AD', 'GR', 'EG')}
    out = dict(script_sha256=sha, gates=gates, pooled_tip_mm=P, mean_rate_hz=R, gain_vs_fixed_frontier=gain,
               EG_below_EGsh_runs=wins, ratio_EG_EGsh=P['EG'] / P['EGsh'], verdict=verdict, described=desc, runs=runs,
               minutes=round((time.time() - t_start) / 60, 1))
    json.dump(out, open(os.path.join(HERE, 'results', 'egeodesic_cadence.json'), 'w'), indent=1)
    print('gates', gates)
    print('pooled tip (mm)', {a: round(v, 3) for a, v in P.items()}); print('rate (Hz)', {a: round(v) for a, v in R.items()})
    print('gain vs fixed frontier', {a: f'{v:.1%}' for a, v in gain.items()}, '; EG below EGsh in', wins, 'of 20')
    print('described', {a: {k: round(v, 3) for k, v in d.items()} for a, d in desc.items()})
    print('verdict', verdict)
