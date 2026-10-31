# Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
# Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
"""Cadence bench, timing control: does the law know WHEN to speed up? Written 31 Oct 2026 after the cadence bench's
verdict, before any run of this file; design and threshold agreed with JP ("Oui").

Reading after the verdict that motivates it (not pre-registered): at about the same frame budget the law's frames did
better than fixed ones (AD 0.898 mm at 541 Hz; fixed 2 ms 0.979 mm at 500 Hz; the fixed frontier interpolated to
541 Hz about 0.968 mm). That could come from the law's timing, or only from having more fast frames on average and
spreading them over the joints. This test keeps AD's frames and takes away their timing.

Arms, on the same 20 paired runs as cadence_bench.py (2 kHz, baton arm, excerpts A and B, steady noise and the x10
step, seeds 1 to 5):
  AD   the law's live choice over {1, 2, 3} ms (cadence_bench.py, unchanged).
  SH   AD's own schedule from the same run, each joint's sends and frames (hence its floors) replayed with the scored
       30 s shifted circularly by 15 s: the same frames, the same shares of 1, 2 and 3 ms, the same number of
       changes, no longer lined up with the motion or the noise. The 3 s before the score are AD's, unshifted.
  JF   each joint at one fixed frame for the whole run, the period on the 0.5 ms grid closest to its mean period in
       AD in that run (1 to 3 ms), its floor the engine's for that period: the joint allocation of AD, no timing.
  ADr  AD's schedule replayed unshifted (a check of the replay).
Gates: G1, AD and the four other cadence arms of pass 1 reproduce results/summary.json run by run to 1e-9 mm;
G2, ADr reproduces AD run by run to 1e-6 mm. If a gate fails, nothing is judged.

Prediction (judged):
T1  AD's pooled tip error is at least 5 % below SH's (AD <= 0.95 x SH), and AD is below SH in at least 18 of the 20
    runs.
Described, not judged: JF against AD and SH; the actual mean frame rate of each arm; time on the floor.
If T1 holds, the law places its fast frames at the right moments. If AD and SH tie, the gain came from the budget.
Before this file was committed it was run only as a smoke test (1.5 s of excerpt A: AD replayed unshifted equals AD
live exactly; no arm compared).

Run (from the repository root, about 15 minutes):   python bench/cadence_2026-10-31/timing_control.py
"""
import os, sys, json, time, hashlib
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import cadence_bench as CB
from cadence_bench import CadenceEngine, ARMS, FRAMES_MS, SEEDS, Z, TAU_S, W0, ENC_SD, ENC_LSB, TS, DEG, tip
from syntonic_engine.engine import kappa_H, window_floor
from scipy.interpolate import CubicSpline
from scipy.signal import butter, sosfiltfilt

DT = 5e-4; SHIFT_S = 15.0
LAST = {}


class Recorder(CadenceEngine):
    """CadenceEngine, unchanged, that records per tick the floor in force for the window update and the sends."""

    def __init__(self, *a, **kw):
        super().__init__(*a, **kw); self.rec_floor = []; self.rec_send = []; LAST['eng'] = self

    def step(self, q_plan, read_a, read_b):
        fl = self.floor.copy(); s0 = self.sent.copy()
        cmd = super().step(q_plan, read_a, read_b)
        self.rec_floor.append(fl); self.rec_send.append(self.sent != s0)
        return cmd


class Imposed(CadenceEngine):
    """The same arithmetic, with the floor and the sends imposed per tick by a schedule."""

    def __init__(self, q0, dt, w0, sched_floor, sched_send):
        n = len(q0)
        super().__init__(q0, dt, w0, np.ones((n, len(FRAMES_MS)), bool), np.zeros(n, bool))
        self.sf, self.ss = sched_floor, sched_send

    def step(self, q_plan, read_a, read_b):
        dt = self.dt; k = self.k
        self.floor = self.sf[k]
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
        send = self.ss[k]
        self.cmd = np.where(send, q + self.u, self.cmd); self.sent += send
        self.th_m = self.th_m + dt * np.clip((self.cmd - self.th_m) / self.servo_tau, -self.w0, self.w0)
        self.th_ff = self.th_ff + dt * np.clip((q - self.th_ff) / self.servo_tau, -self.w0, self.w0)
        self.k += 1
        return self.cmd


def imposed_run(name, cond, seed, sched_floor, sched_send, labels, t_end=30.0):
    """Plant loop of cadence_bench.simulate, baton arm, columns given by the schedules (n_ticks, A*6)."""
    t0 = float(Z[f'{name}_t0']); tt = np.arange(t0 - 3.0, t0 + t_end, DT)
    joints = list(range(6)); J = 6; A = len(labels)
    qJ = CubicSpline(Z[f'{name}_t30'], Z[f'{name}_q30'][:, joints], axis=0)(tt)
    tauJ = CubicSpline(Z[f'{name}_t30'], Z[f'{name}_tau30'][:, joints], axis=0)(tt)
    sd = np.full(len(tt), ENC_SD)
    if cond == 'step': sd[(tt >= t0 + 10) & (tt < t0 + 20)] *= np.sqrt(10.0)
    rng = np.random.default_rng(seed)
    tauS = np.tile(TAU_S[joints], A); w0 = np.tile(W0[joints], A)
    k_inv = (1.0 * DEG) / tauS; hb = 0.5 * 0.3 * DEG
    eng = Imposed(np.tile(qJ[0], A), DT, w0, sched_floor, sched_send)
    th_g = np.tile(qJ[0], A); every = int(round(0.01 / DT)); keep = []; below = np.zeros(A * J); m = 0; sent0 = None
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
            if sent0 is None: sent0 = eng.sent.copy()
            below += eng.law < eng.floor; m += 1
            if k % every == 0: keep.append((th_out.copy(), q.copy()))
    T = len(tt) * DT - 3.0; out = {}
    for i, lab in enumerate(labels):
        s = slice(i * J, (i + 1) * J)
        e = np.array([tip(a[s]) - tip(b[s]) for a, b in keep])
        hp = sosfiltfilt(butter(4, 5.0, 'highpass', fs=100.0, output='sos'), e, axis=0)
        out[lab] = dict(tip_mm=float(1e3 * np.sqrt(np.mean(np.sum(e ** 2, 1)))),
                        jitter_mm=float(1e3 * np.sqrt(np.mean(np.sum(hp ** 2, 1)))),
                        rate_hz=float(np.mean(eng.sent[s] - sent0[s]) / T), below_floor=float(np.mean(below[s] / m)))
    return out


if __name__ == '__main__':
    sha = hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest()
    ref = json.load(open(os.path.join(HERE, 'results', 'summary.json')))['runs']
    CB.CadenceEngine = Recorder                                  # pass 1 runs cadence_bench.simulate itself
    ia = list(ARMS).index('AD'); cols = slice(ia * 6, (ia + 1) * 6)
    t_start = time.time(); runs = {}; g1 = g2 = 0.0
    for name in ('A', 'B'):
        for cond in ('steady', 'step'):
            for seed in SEEDS:
                key = f'{name} {cond} {seed}'
                live = CB.simulate(name, cond, seed, DT, list(range(6)), list(ARMS))
                g1 = max([g1] + [abs(live[a]['tip_mm'] - ref[key][a]['tip_mm']) for a in ARMS])
                eng = LAST['eng']
                fl = np.array(eng.rec_floor)[:, cols]; sn = np.array(eng.rec_send)[:, cols]
                n0 = int(round(3.0 / DT)); N = fl.shape[0] - n0; S = int(round(SHIFT_S / DT))
                idx = np.arange(fl.shape[0]); idx[n0:] = n0 + (np.arange(N) + S) % N
                fl_sh, sn_sh = fl[idx], sn[idx]
                scored = sn[n0:].sum(0) / (N * DT)                       # AD's mean frame rate per joint (Hz)
                P = np.clip(np.round(1.0 / scored / DT).astype(int), 2, 6)  # JF period in ticks (1 to 3 ms)
                kk = np.arange(fl.shape[0])[:, None]
                sn_jf = (kk % P[None, :]) == 0
                fl_jf = np.broadcast_to(np.array([window_floor(DT, p * DT) for p in P])[None, :], fl.shape)
                sf = np.concatenate([fl, fl_sh, fl_jf], axis=1); ss = np.concatenate([sn, sn_sh, sn_jf], axis=1)
                imp = imposed_run(name, cond, seed, sf, ss, ['ADr', 'SH', 'JF'])
                g2 = max(g2, abs(imp['ADr']['tip_mm'] - live['AD']['tip_mm']))
                runs[key] = dict(AD=live['AD'], **imp, JF_period_ms=(P * DT * 1e3).tolist())
                print(key, {a: round(runs[key][a]['tip_mm'], 3) for a in ('AD', 'ADr', 'SH', 'JF')},
                      {a: round(runs[key][a]['rate_hz']) for a in ('AD', 'SH', 'JF')}, 'g1', g1, 'g2', g2,
                      round(time.time() - t_start), 's', flush=True)
    pooled = {a: float(np.sqrt(np.mean([runs[k][a]['tip_mm'] ** 2 for k in runs]))) for a in ('AD', 'ADr', 'SH', 'JF')}
    rate = {a: float(np.mean([runs[k][a]['rate_hz'] for k in runs])) for a in ('AD', 'SH', 'JF')}
    wins = int(sum(runs[k]['AD']['tip_mm'] < runs[k]['SH']['tip_mm'] for k in runs))
    gates = dict(G1_max_diff_mm=g1, G1_ok=g1 <= 1e-9, G2_max_diff_mm=g2, G2_ok=g2 <= 1e-6)
    T1 = bool(gates['G1_ok'] and gates['G2_ok'] and pooled['AD'] <= 0.95 * pooled['SH'] and wins >= 18)
    out = dict(script_sha256=sha, gates=gates, pooled_tip_mm=pooled, mean_rate_hz=rate, AD_below_SH_runs=wins,
               ratio_AD_SH=pooled['AD'] / pooled['SH'], ratio_AD_JF=pooled['AD'] / pooled['JF'], T1=T1, runs=runs,
               minutes=round((time.time() - t_start) / 60, 1))
    json.dump(out, open(os.path.join(HERE, 'results', 'timing_control.json'), 'w'), indent=1)
    print('gates', gates)
    print('pooled tip (mm)', {a: round(v, 3) for a, v in pooled.items()}, 'rate (Hz)', {a: round(v) for a, v in rate.items()})
    print('AD below SH in', wins, 'of 20 runs; AD / SH', round(pooled['AD'] / pooled['SH'], 3), '; T1', T1)
