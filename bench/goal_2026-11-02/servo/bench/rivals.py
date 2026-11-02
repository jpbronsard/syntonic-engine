"""Rival controllers for realistic_bench.py: the published PID (optionally with the servo lead, optionally with the
frame-mean proportional term), the published syntonic engine 0.2.9, and the PID tuning of the published bench.

Every controller has step(q_plan, qd_plan, read_a, read_b, send) -> command (rad), the interface of
realistic_bench; a factory(q0, w0) builds it (q0: first planned positions, w0: data-sheet no-load speeds).

PID (examples/servo_bench.py, arithmetic unchanged): error e = plan - mean of the two reads; derivative on the error
low-passed at 5 ms; integral frozen while the raw command exceeds the clip; correction clipped at
servo_tau x w0 (the authority limit of 0.2.8, servo_tau = the nominal 25 ms); command at a frame = plan + clipped
correction. Options:
  lead    + (servo_tau + D) x planned velocity at frames, D = loop_delay(dt, frame) = frame/2 + dt/2 + 0.1 ms
          = 2.1 ms at 1 kHz / 3 ms (27.1 ms in all), as in bench/maestro_2026-11-01/pid_lead_fair.py.
  p_mean  frame-mean P (coordinator, 2 Nov 2026): at a frame, the proportional term of the command uses the mean of e
          over the ticks since the previous frame, the current tick included (ticks k-2, k-1, k at 1 kHz / 3 ms;
          tick 0 alone at the first frame). Unchanged: the integral accumulates e dt every tick, the derivative, and
          the anti-windup test, made every tick on the published raw command (P on that tick's e). The frame
          command's I and D terms are those of the published raw command of the same tick.
  Primary PID rival: PID(lead=True, p_mean=True), tuned.

Syntonic: SyntonicEngine 0.2.9 as published, identification='raw', nothing tuned; lead='servo+delay' when lead=True
(the same 27.1 ms coefficient, computed by the engine).

Tuning (tune_pid): joint by joint, per excerpt, on the published 280-point grid (Kp in {0, 0.5, 1, 2, 3, 4, 6, 8},
Ki in {0, 80, 160, 320, 640, 1280, 2560} /s, Kd in {0, 0.5, 1, 2, 5} ms), criterion the joint's rms angle error over
the 30 s (mean square over the tuning seeds 61, 62, 63), steady noise, circumstance 'nominal'. Edge rule: if any
joint of any tuned excerpt picks the top value of an axis, that axis is extended once (Kp + {12, 16}; Ki + {5120};
Kd + {10 ms}) and every excerpt is re-tuned on the extended grid (the rule of pid_lead_fair.py). All grid points run
as columns of one simulation per excerpt and seed.
"""
import json
import time

import numpy as np

import realistic_bench as RB
from syntonic_engine import SyntonicEngine
from syntonic_engine.engine import loop_delay

__all__ = ['PID', 'Syntonic', 'pid_factory', 'syntonic_factory', 'published_gains', 'lead_coefficient',
           'tune_pid', 'grid', 'KP', 'KI', 'KD', 'EXT', 'TUNE_SEEDS', 'load_tuned']

KP = [0, 0.5, 1, 2, 3, 4, 6, 8]
KI = [0, 80, 160, 320, 640, 1280, 2560]
KD = [0, 0.0005, 0.001, 0.002, 0.005]
EXT = dict(kp=[12, 16], ki=[5120], kd=[0.01])
TUNE_SEEDS = (61, 62, 63)


def lead_coefficient(dt=RB.DT, frame=RB.FRAME, servo_tau=RB.TS_NOM):
    """(servo_tau + D) with D = loop_delay(dt, frame): 0.025 + 0.0021 = 0.0271 s at 1 kHz / 3 ms."""
    return servo_tau + loop_delay(dt, frame)


class PID:
    """The PID of examples/servo_bench.py, vectorised over columns.

    q0        initial command (rad), shape (n,)
    gains     (n, 3): Kp, Ki (1/s), Kd (s) per column
    w0        data-sheet no-load speed per column (rad/s); clip = servo_tau x w0
    lead      lead coefficient (s) on the planned velocity, added to the command at frames; 0 = off
    p_mean    frame-mean proportional term (see the module docstring)
    frozen    per column, the number of ticks the integral was frozen (anti-windup active), for diagnostics
    """

    def __init__(self, q0, gains, w0, dt=RB.DT, frame=RB.FRAME, servo_tau=RB.TS_NOM, lead=0.0, p_mean=False):
        g = np.asarray(gains, dtype=float)
        n = len(np.asarray(q0))
        if g.shape != (n, 3):
            raise ValueError(f'gains must have shape ({n}, 3), got {g.shape}')
        self.kp, self.ki, self.kd = g[:, 0].copy(), g[:, 1].copy(), g[:, 2].copy()
        self.dt = float(dt)
        self.lim = servo_tau * np.asarray(w0, dtype=float)
        self.lead = lead
        self.p_mean = bool(p_mean)
        self.I = np.zeros(n); self.ef = np.zeros(n); self.efp = np.zeros(n)
        self.acc = np.zeros(n); self.cnt = 0
        self.frozen = np.zeros(n, dtype=int)
        self.cmd = np.asarray(q0, dtype=float).copy()

    def step(self, q_plan, qd_plan, read_a, read_b, send):
        dt = self.dt
        e = q_plan - 0.5 * (read_a + read_b)
        self.ef += (dt / 0.005) * (e - self.ef)
        raw = self.kp * e + self.ki * self.I + self.kd * (self.ef - self.efp) / dt
        if self.p_mean:
            self.acc += e; self.cnt += 1
            if send:                                   # P on the mean error since the previous frame
                raw_f = self.kp * (self.acc / self.cnt) + self.ki * self.I + self.kd * (self.ef - self.efp) / dt
                self.acc = np.zeros_like(self.acc); self.cnt = 0
        fz = np.abs(raw) >= self.lim
        self.I += np.where(fz, 0.0, e * dt); self.efp = self.ef.copy()
        self.frozen += fz
        if send:
            corr = np.clip(raw_f if self.p_mean else raw, -self.lim, self.lim)
            if np.any(self.lead):
                self.cmd = q_plan + corr + self.lead * qd_plan
            else:
                self.cmd = q_plan + corr
        return self.cmd


class Syntonic:
    """SyntonicEngine 0.2.9 as published (identification 'raw', nothing tuned), optionally with lead='servo+delay'."""

    def __init__(self, q0, w0, dt=RB.DT, frame=RB.FRAME, lead=False):
        self.lead = bool(lead)
        opts = dict(lead='servo+delay') if self.lead else {}
        self.engine = SyntonicEngine(q0, dt, frame, w0, identification='raw', **opts)

    def step(self, q_plan, qd_plan, read_a, read_b, send):
        return self.engine.step(q_plan, read_a, read_b, frame=send, qd_plan=qd_plan if self.lead else None)


def published_gains(excerpt=None):
    """The published PID gains (Z['<excerpt>_pid'], 12 joints x (Kp, Ki, Kd)); a dict A/B when excerpt is None."""
    if excerpt is None:
        return {n: np.array(RB.Z[f'{n}_pid']) for n in ('A', 'B')}
    return np.array(RB.Z[f'{excerpt}_pid'])


def pid_factory(gains, lead=False, p_mean=False):
    """Factory of a PID: factory(q0, w0) -> PID. gains: (J, 3) rows RJ1, RJ2, ... (the first len(q0) rows are used,
    so the published 12-row tables can be passed), or a dict excerpt -> (J, 3), which returns a dict of factories.
    lead: False, True (lead_coefficient()) or a coefficient in seconds. p_mean: frame-mean P."""
    if isinstance(gains, dict):
        return {n: pid_factory(g, lead, p_mean) for n, g in gains.items()}
    G = np.asarray(gains, dtype=float)
    ld = lead_coefficient() if lead is True else (float(lead) if lead else 0.0)

    def make(q0, w0):
        n = len(q0)
        if len(G) < n:
            raise ValueError(f'{len(G)} gain rows for {n} joints')
        return PID(q0, G[:n], w0, RB.DT, RB.FRAME, RB.TS_NOM, ld, p_mean)
    return make


def syntonic_factory(lead=False):
    """Factory of the published engine (identification 'raw'), with or without lead='servo+delay'."""
    def make(q0, w0):
        return Syntonic(q0, w0, RB.DT, RB.FRAME, lead=lead)
    return make


# ---------------------------------------------------------------------------------------------- tuning
def grid(kp=KP, ki=KI, kd=KD):
    """(Kp, Ki, Kd) triples, Kp slowest, Kd fastest (the order of pid_lead_fair.py)."""
    return [(a, b, c) for a in kp for b in ki for c in kd]


def _tune_on_grid(G, lead, p_mean, excerpts, seeds, circumstance, log):
    G = np.array(G, dtype=float); n = len(G)
    cols = np.repeat(np.arange(6), n)                                # column = joint * n + grid point
    gains = np.column_stack([np.tile(G[:, 0], 6), np.tile(G[:, 1], 6), np.tile(G[:, 2], 6)])
    make = pid_factory(gains, lead=lead, p_mean=p_mean)
    out = {}
    for name in excerpts:
        J2, frozen = [], 0
        for s in seeds:
            t = time.time()
            r = RB.simulate_columns(name, circumstance, s, [('T', make, cols)], n_joints=6, tip_score=False)
            J2.append(r['T']['joint_rms_deg'] ** 2)
            if log:
                log(f'  tune lead={bool(lead)} p_mean={bool(p_mean)} {name} seed {s}: {len(cols)} columns, '
                    f'{time.time() - t:.1f} s')
        J2 = np.mean(J2, 0).reshape(6, n)
        best = np.argmin(J2, 1)
        out[name] = dict(gains=G[best].tolist(), joint_deg=np.sqrt(J2[np.arange(6), best]).tolist(),
                         grid_joint_deg=np.sqrt(J2).tolist())
    return out


def _edges(tuned, axes_top):
    hit = []
    for ax, (i, top) in axes_top.items():
        if any(g[i] >= top for t in tuned.values() for g in t['gains']):
            hit.append(ax)
    return hit


def tune_pid(lead, p_mean=False, excerpts=('A', 'B'), seeds=TUNE_SEEDS, circumstance='nominal', edge_rule=True,
             log=print):
    """Tune the PID (lead and p_mean as in PID) joint by joint, per excerpt, on the published grid; returns
    {'gains': {excerpt: [[Kp, Ki, Kd]] x 6}, 'joint_deg': {excerpt: [rms deg] x 6}, 'extended': {axis: values},
     'grid_size': int, 'first_pass': {...} (the 280-grid choice when the edge rule fired), 'wall_s': float,
     '_grids': per excerpt the (6, grid) rms joint errors}."""
    t = time.time()
    first = _tune_on_grid(grid(), lead, p_mean, excerpts, seeds, circumstance, log)
    tops = dict(kp=(0, max(KP)), ki=(1, max(KI)), kd=(2, max(KD)))
    axes = _edges(first, tops) if edge_rule else []
    ext = {ax: EXT[ax] for ax in axes}
    tuned, size = first, len(grid())
    if axes:                                                     # the declared edge rule, applied once
        G = grid(KP + ext.get('kp', []), KI + ext.get('ki', []), KD + ext.get('kd', []))
        tuned, size = _tune_on_grid(G, lead, p_mean, excerpts, seeds, circumstance, log), len(G)
    res = dict(gains={n: tuned[n]['gains'] for n in excerpts}, joint_deg={n: tuned[n]['joint_deg'] for n in excerpts},
               extended=ext, grid_size=size, seeds=list(seeds), circumstance=circumstance, lead=bool(lead),
               p_mean=bool(p_mean), lead_coefficient_s=lead_coefficient() if lead else 0.0, wall_s=time.time() - t)
    if axes:
        res['first_pass'] = {n: dict(gains=first[n]['gains'], joint_deg=first[n]['joint_deg']) for n in excerpts}
    res['_grids'] = {n: tuned[n]['grid_joint_deg'] for n in excerpts}
    return res


def load_tuned(path):
    """Tuned gains saved by run_baselines.py tune: {arm: {excerpt: (6, 3) array}} for the arms
    'leadmean' (PID + lead + frame-mean P), 'lead' (PID + lead) and 'nolead' (PID)."""
    T = json.load(open(path))
    return {arm: {n: np.array(g) for n, g in v['gains'].items()} for arm, v in T.items() if isinstance(v, dict)
            and 'gains' in v}
