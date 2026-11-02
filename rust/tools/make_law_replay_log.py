# Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
# Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
"""Replay logs of the two laws of 0.3.0 (Python), for the C++ and Rust ports of 0.3.1 to replay and compare.

Runs python/syntonic_engine's SyntonicLoopLaw and RiccatiAttitudeLaw on deterministic inputs and records, every tick,
what each law was given and what it computed. Written to the current folder, little-endian f64:

  replay_loop_bench.bin         the loop law on the realistic MAESTRO 1 bench of bench/goal_2026-11-02/servo (baton arm,
                                6 joints, excerpt A, development seed 91, encoder noise x10 from 10 to 20 s, servo delay
                                2.5 ms, i.e. a fractional delay line), 33 s
  replay_loop_bench_spread.bin  the same bench with the C++/Rust option spread = 400 (reference: SpreadLoopLaw below)
  replay_loop_two.bin           two joints (as on the ARM-1 bench), a closed loop of a 25 ms servo with 2 ms of delay
                                (whole ticks), a fixed lead of 20 ms, a disturbance step and a read-noise step, 12 s
  replay_loop_one.bin           one joint, the same closed loop, 8 s (one joint sums the fit's terms in another order)
  replay_att_synth.bin          the attitude law, three levels, a synthetic flight: 5 s on the ground, arming, then
                                low, half, low (tilt noise x3 for 3 s), idle, half, low, with the attitude moving, 35 s
  replay_att_bench.bin          (with --flight-bench) the attitude law on the realistic flight bench of
                                bench/goal_2026-11-02/flight_realistic, circumstance C5 (half and low every 5 s),
                                development seed 300, its facts from the fact files; needs the paths of replay_setup.sh

Format. A header of 16 values, then the facts, then one row per tick.
  loop: header [1, N, ticks, dt, frame, servo_tau, dead, lead (-1 = None), spread, every, 0 ...]; w0[N]; q0[N];
        rows q[N] qd[N] a[N] b[N] send | cmd[N] Kp[N] Ki[N] I[N] R[N] A[N] gam[N] Oc[N] th_m[N] th_ff[N] gain_index[N]
        (gain_index: the candidate c = iKp x 60 + iKi in use, -1 before the first choice)
  attitude: header [2, L, ticks, dt, dt_accel, q_b, T_pre, ground, settle, ema, scale0, scale1, scale2, 0 ...];
        s2[L] lam[L] T_rec[L]; rows g a level arm | x[1+L] P[(1+L)^2] noise_reading[L] thr
        (a = NaN between accelerometer samples; arm = 1 when arm() is called right after this step; the outputs are
        those after it)

It also reports, at every gain recomputation of the loop law, the relative gap between the smallest and the next
smallest predicted loop error J among the robust candidates: the margin a port has before rounding could change the
chosen gains.
Run from the repository root (needs scipy for the servo bench):  python rust/tools/make_law_replay_log.py [--flight-bench]
"""
import os, sys
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, 'python'))
from syntonic_engine import SyntonicLoopLaw, RiccatiAttitudeLaw  # noqa: E402
from syntonic_engine.loop_law import K_gamma  # noqa: E402

SPREAD = 400


class SpreadLoopLaw(SyntonicLoopLaw):
    """Python reference of the C++/Rust option spread = m (0.3.1): the gains computed from the statistics of tick k (a
    gain-update tick) take effect in tick k + m instead of tick k. spread = 0 is SyntonicLoopLaw. step() is that of
    SyntonicLoopLaw with only the gain assignment changed."""

    def __init__(self, *args, spread=0, **kw):
        super().__init__(*args, **kw)
        self.spread = int(spread); self.pending = None

    def step(self, q, qd, a, b, send):
        dt = self.dt
        y = 0.5 * (a + b); dn = (a - b) / np.sqrt(2.0)
        self.n2 += self.aw * (dn * dn - self.n2); self.R = np.maximum(0.5 * self.n2, 1e-14)
        self._identify((q - self.th_ff) + (self.th_m - y))
        if self.started and (self.k % self.every == 0):
            g = self.LL.gains(self.A, self.gam, self.Oc, self.R * dt)
            if self.spread == 0:
                self.Kp, self.Ki = g
            else:
                self.pending = (self.k + self.spread, g)
        elif self.pending is not None and self.k == self.pending[0]:
            self.Kp, self.Ki = self.pending[1]; self.pending = None
        e = q - y; self.esum += e; self.ecnt += 1
        self.I += np.where(np.abs(self.Kp * e + self.I) >= self.lim, 0.0, self.Ki * e * dt)
        if send:
            eb = self.esum / self.ecnt; self.esum[:] = 0.0; self.ecnt = 0
            self.cmd = q + self.L * qd + np.clip(self.Kp * eb + self.I, -self.lim, self.lim); self.cff = q + self.L * qd
        cm = self._delayed(self.pc, self.cmd); cf = self._delayed(self.pf, self.cff)
        self.th_m = self.th_m + dt * np.clip((cm - self.th_m) / self.ts, -self.w0, self.w0)
        self.th_ff = self.th_ff + dt * np.clip((cf - self.th_ff) / self.ts, -self.w0, self.w0)
        self.k += 1
        return self.cmd


def j_gap(law):
    """relative gap between the best and the second best J per joint, for the statistics the law just used"""
    LL = law.LL; C = law.A / K_gamma(law.gam)
    Phi = C[:, None] * (LL.W[None, :] ** 2 + law.Oc[:, None] ** 2) ** (-(1 + law.gam[:, None]) / 2)
    J = Phi @ LL.MS.T + (law.R * law.dt)[:, None] * LL.MTr[None, :]
    J = np.sort(np.where(LL.ok[None, :], J, np.inf), axis=1)
    return (J[:, 1] - J[:, 0]) / np.abs(J[:, 0])


class LoopRecorder:
    """a controller for the benches: the law, with every tick's inputs and outputs recorded"""

    def __init__(self, law):
        self.law = law; self.rows = []; self.gaps = []
        self.index = {(kp, ki): c for c, (kp, ki) in enumerate(zip(law.LL.KP, law.LL.KI))}

    def step(self, q, qd, a, b, send):
        law = self.law; k = law.k
        cmd = law.step(q, qd, a, b, send)
        if law.started and k % law.every == 0:
            self.gaps.append(j_gap(law))
        gi = [float(self.index.get((kp, ki), -1)) for kp, ki in zip(law.Kp, law.Ki)]
        self.rows.append(np.concatenate([q, qd, a, b, [1.0 if send else 0.0], cmd, law.Kp, law.Ki, law.I, law.R, law.A,
                                         law.gam, law.Oc, law.th_m, law.th_ff, gi]))
        return cmd


def write_loop(fname, rec, q0, w0, dt, frame, servo_tau, dead, lead, spread, every=500):
    n = len(q0)
    head = np.zeros(16); head[:10] = [1, n, len(rec.rows), dt, frame, servo_tau, dead, -1.0 if lead is None else lead,
                                     spread, every]
    np.concatenate([head, w0, q0, np.array(rec.rows).ravel()]).astype('<f8').tofile(fname)
    g = np.concatenate(rec.gaps) if rec.gaps else np.array([np.inf])
    print(f'wrote {fname}: {n} joint(s), {len(rec.rows)} ticks, {len(rec.gaps)} gain recomputations, '
          f'smallest relative gap between the best and the next J: {g.min():.2e}')


def loop_bench(spread):
    os.environ.setdefault('SYNTONIC_REPO', ROOT)
    sys.path.insert(0, os.path.join(ROOT, 'bench', 'goal_2026-11-02', 'servo', 'bench'))
    import realistic_bench as RB
    dead = 0.0025; made = {}

    def factory(q0, w0):
        kw = dict(dt=RB.DT, frame=RB.FRAME, servo_tau=RB.TS_NOM, dead=dead)
        law = SyntonicLoopLaw(q0, w0, **kw) if spread == 0 else SpreadLoopLaw(q0, w0, spread=spread, **kw)
        made['rec'] = LoopRecorder(law); made['q0'] = q0.copy(); made['w0'] = w0.copy()
        return made['rec']

    RB.simulate_columns('A', dict(dead_time=dead, noise='step'), 91, [('LAW', factory, None)], n_joints=6, tip_score=False)
    name = 'replay_loop_bench.bin' if spread == 0 else 'replay_loop_bench_spread.bin'
    write_loop(name, made['rec'], made['q0'], made['w0'], RB.DT, RB.FRAME, RB.TS_NOM, dead, None, spread)


def loop_closed(n, seconds, fname, seed):
    """a closed loop as in python/tests/test_loop_law.py: first-order 25 ms servo with a speed limit, 2 ms of delay,
    an output disturbance (random walk + a 0.01 rad step at 1.5 s), two noisy reads per tick, read noise x10 from
    seconds / 2 for 2 s; a fixed lead of 20 ms"""
    rng = np.random.default_rng(seed); dt, frame, ts, dead_ticks, lead = 1e-3, 3e-3, 0.025, 2, 0.02
    w0 = np.full(n, 6.0) + 0.5 * np.arange(n); t = np.arange(0.0, seconds, dt); f = np.array([0.4, 0.7, 1.1])[:n]
    q = 0.3 * np.sin(2 * np.pi * f * t[:, None]); qd = 0.6 * np.pi * f * np.cos(2 * np.pi * f * t[:, None])
    dist = np.cumsum(rng.normal(0, 2e-4, (len(t), n)), 0) + np.where(t[:, None] >= 1.5, 0.01, 0.0)
    sd = np.where((t >= seconds / 2) & (t < seconds / 2 + 2.0), 2e-4 * np.sqrt(10.0), 2e-4)
    law = SyntonicLoopLaw(q[0], w0, dt=dt, frame=frame, servo_tau=ts, dead=dead_ticks * dt, lead=lead)
    rec = LoopRecorder(law)
    th = q[0].copy(); latched = q[0].copy(); pipe = [q[0].copy()] * (dead_ticks + 1)
    for k in range(len(t)):
        y = th + dist[k]
        a = y + rng.normal(0, sd[k], n); b = y + rng.normal(0, sd[k], n)
        send = k % 3 == 0
        cmd = rec.step(q[k], qd[k], a, b, send)
        if send:
            latched = cmd.copy()
        pipe.append(latched.copy()); u = pipe.pop(0)
        th = th + dt * np.clip((u - th) / ts, -w0, w0)
    write_loop(fname, rec, q[0], w0, dt, frame, ts, dead_ticks * dt, lead, 0)


LEVEL_NAMES = ('idle', 'low', 'half')


class AttRecorder:
    """an estimator for the benches: the attitude law, with every sample's inputs and outputs recorded"""

    def __init__(self, law):
        self.law = law; self.n = 1; self.rows = []

    def out(self):
        law = self.law
        return np.concatenate([law.x, law.P.ravel(), [law.noise_reading(r) for r in law.levels], [law.thr]])

    def step(self, g, a, reg_lam=None, t=None, level=None, k=None):
        th = self.law.step(g, a, level=level)
        self.rows.append(np.concatenate([[g, np.nan if a is None else a, float(self.law.index[level]), 0.0], self.out()]))
        return th

    def arm(self, bias_cal=None):
        self.law.arm(bias_cal)
        self.rows[-1][3] = 1.0; self.rows[-1][4:] = self.out()


def write_att(fname, rec, levels, q_b, T_pre, dt, dta):
    law = rec.law; L = len(law.levels)
    head = np.zeros(16); head[:13] = [2, L, len(rec.rows), dt, dta, q_b, T_pre, law.index[law.ground], law.settle, 1.0,
                                      *law.H]
    facts = [levels[r]['s2'] for r in law.levels] + [levels[r]['lam'] for r in law.levels] + \
            [levels[r].get('T_rec', 0.0) for r in law.levels]
    np.concatenate([head, facts, np.array(rec.rows).ravel()]).astype('<f8').tofile(fname)
    print(f'wrote {fname}: {L} levels, {len(rec.rows)} gyro samples, '
          f'{int(np.sum(~np.isnan(np.array(rec.rows)[:, 1])))} accelerometer samples')


def att_synth(fname='replay_att_synth.bin', seed=3):
    DT, DTA, NSUB, T_PRE, QB = 1e-3, 4e-3, 4, 5.0, 9.2e-11
    LAM = {'idle': 9e-9, 'low': 3.3e-8, 'half': 1.2e-6}; S2 = {'idle': 1.9e-5, 'low': 2.2e-5, 'half': 3.5e-3}
    TREC = {'low': 40.0, 'half': 8.0}; BIAS = {'idle': 6e-4, 'low': 6.27e-4, 'half': 8e-4}
    levels = {r: dict(s2=S2[r], lam=LAM[r], **({'T_rec': TREC[r]} if r in TREC else {})) for r in LEVEL_NAMES}
    rec = AttRecorder(RiccatiAttitudeLaw(levels, q_b=QB, T_pre=T_PRE, dt=DT, dt_accel=DTA))
    rng = np.random.default_rng(seed)
    plan = [('idle', T_PRE), ('low', 10.0), ('half', 6.0), ('low', 5.0), ('idle', 2.0), ('half', 4.0), ('low', 3.0)]
    kpre = int(round(T_PRE / DT)); gsum = cal = 0.0; k = 0
    for level, length in plan:
        for _ in range(int(round(length / DT))):
            tf = k * DT - T_PRE
            th = 0.01 + (0.02 * np.sin(0.3 * tf) if tf > 0 else 0.0)
            rate = 0.006 * np.cos(0.3 * tf) if tf > 0 else 0.0
            g_raw = rate + BIAS[level] + rng.normal(0.0, np.sqrt(LAM[level] / DT)); a = None
            if k % NSUB == 0:
                sd = np.sqrt(S2[level] / DTA) * (3.0 if 3.0 <= tf < 6.0 else 1.0)
                a = th + rng.normal(0.0, sd)
            rec.step(g_raw - cal, a, level=level)
            if k < kpre:
                gsum += g_raw
                if k == kpre - 1:
                    cal = gsum / kpre; rec.arm(cal)
            k += 1
    write_att(fname, rec, levels, QB, T_PRE, DT, DTA)


def att_bench(fname='replay_att_bench.bin', seed=300, circ='C5'):
    import json
    snap = '/home/claude/work/goal/confirm/flight_snap'
    if not os.path.isdir(snap):
        sys.exit(f'{snap} not found: run bench/goal_2026-11-02/replay_setup.sh first (see its NOTE.md)')
    here = os.path.join(ROOT, 'bench', 'goal_2026-11-02')
    sys.path.insert(0, os.path.join(here, 'flight_realistic'))
    import bench_bb as B
    B.install()
    F2 = json.load(open(os.path.join(here, 'flight_looped', 'flight_snap', 'results', 'facts2.json')))
    FX = json.load(open(os.path.join(ROOT, 'bench', 'flight_2026-11-02', 'results_a2', 'facts.json')))
    levels = {r: dict(s2=F2['curves']['raw'][r]['s2_rec'], lam=FX['lam'][r], T_rec=F2['curves']['raw'][r]['T_rec'])
              for r in LEVEL_NAMES}
    rec = AttRecorder(RiccatiAttitudeLaw(levels, q_b=F2['q_b_half'], T_pre=B.C.T_PRE, dt=B.C.DT, dt_accel=B.C.DTA))
    B.run4(seed, circ, [rec])
    write_att(fname, rec, levels, F2['q_b_half'], B.C.T_PRE, B.C.DT, B.C.DTA)


if __name__ == '__main__':
    loop_bench(0)
    loop_bench(SPREAD)
    loop_closed(2, 12.0, 'replay_loop_two.bin', seed=7)
    loop_closed(1, 8.0, 'replay_loop_one.bin', seed=8)
    att_synth()
    if '--flight-bench' in sys.argv:
        att_bench()
