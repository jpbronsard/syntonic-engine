# Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
# Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
"""Bench amendment 3 (2 Nov 2026, declared after the failed flight confirmation of confirm_flight.py, whose result stands):
block-bootstrap replay of the logged records instead of looping them.

Why. Amendment 2 replays each record in a loop from a random offset. The half-throttle record is 7.96 s long, so a 60 s
flight replays it 7.5 times: the tilt noise (mean-removed) then sums to exactly zero over every 7.96 s, and the gyro
integral (bias removed) never random-walks beyond 7.96 s. No real sensor does that. The sub-agent flagged it before the
confirmation ("flatters every long-window estimator at half throttle, the rival included"); the truth-tuned filters
exploit it (Q_G = 0, Q_B = 0; with the bias table the half-throttle tilt is ignored altogether).

What. Each throttle segment is filled with consecutive blocks drawn at uniform random offsets of its record (gyro and
tilt kept synchronous inside a block), block length L = T_lim of the record = T_rec / 4, the longest scale at which the
records' own noise facts are measured (idle 8.70 s, low 9.95 s, half 1.99 s). Below L the replayed noise is the record;
beyond L successive blocks are independent, i.e. the noise continues as white rate noise at the level the record shows at
L (the minimal assumption: nothing is known beyond the record). This is the moving-block bootstrap of dependent series.
Everything else is amendment 2 (records, pre-roll, arming calibration, circumstances and their modifiers).
Seeds 500-599 are reserved for the confirmation on this bench.
Run: python3 bench_bb.py validate"""
import os, sys
import numpy as np
SNAP = '/home/claude/work/goal/confirm/flight_snap'
if SNAP not in sys.path: sys.path.insert(0, SNAP)
import common as C
F = C.F; DT, DTA, NSUB, T_PRE, D2R = C.DT, C.DTA, C.NSUB, C.T_PRE, C.D2R
import variants as V
TLIM = {r: V.F2['curves']['raw'][r]['T_lim'] for r in V.REGS}
RESERVED = set(range(500, 600))
_build_loop = C.build_noise_c


def blocks(rng, n_rec, n_out, L):
    """indices into a record of n_rec samples: consecutive blocks of L samples at uniform offsets in [0, n_rec - L]"""
    out = np.empty(n_out, dtype=np.int64); k = 0
    while k < n_out:
        o = int(rng.integers(0, n_rec - L + 1)); m = min(L, n_out - k); out[k:k + m] = o + np.arange(m); k += m
    return out


def build_noise_bb(seed, circ):
    assert seed not in C.FORBIDDEN or os.environ.get('GOAL_ALLOW_CONFIRM') == '1', 'seeds 400-499: first confirmation'
    assert seed not in RESERVED or os.environ.get('GOAL_ALLOW_CONFIRM3') == '1', 'seeds 500-599 are reserved'
    spec = C.CIRC[circ]; R = F.records(); rng = np.random.default_rng(seed)
    segs = [('idle', 0.0, T_PRE)] + [(n, T_PRE + a, T_PRE + b) for n, a, b in spec['segs']]
    gn, tn, reg = [], [], []
    for name, a, b in segs:
        tilt, g = R[name]; ng = int(round((b - a) / DT)); na = ng // NSUB; assert na * NSUB == ng
        L = int(round(TLIM[name] / DTA)); idx_a = blocks(rng, len(tilt), na, L)
        idx_g = (idx_a[:, None] * NSUB + np.arange(NSUB)[None, :]).reshape(-1)
        gn.append(g[idx_g]); tn.append(tilt[idx_a]); reg += [name] * ng
    gn = np.concatenate(gn); tn = np.concatenate(tn).copy(); reg = np.array(reg)
    t_g = np.arange(len(gn)) * DT - T_PRE; t_a = np.arange(len(tn)) * DTA - T_PRE
    if 'tilt_amp' in spec:
        a, b, f = spec['tilt_amp']; tn[(t_a >= a) & (t_a < b)] *= f
    if 'gyro_step' in spec:
        a, db = spec['gyro_step']; gn = gn + np.where(t_g >= a, db * D2R, 0.0)
    return gn, tn, reg


def install():
    C.build_noise_c = build_noise_bb


def validate(nrep=200, T=60.0):
    """the looping artifact and its removal, half and low throttle, against the records' own facts: rms of the T-mean of
    the tilt and rms of the gyro integral over T (record-mean bias removed), white / random-walk theory from s2_rec, lam"""
    R = F.records(); out = {}
    for name in ('half', 'low'):
        tilt, g = R[name]; bias = g.mean(); s2 = V.F2['curves']['raw'][name]['s2_rec']; lam = V.LAM[name]
        na = int(T / DTA); ng = na * NSUB; L = int(round(TLIM[name] / DTA)); res = {}
        for mode in ('loop', 'bb'):
            mt, ig = [], []
            for k in range(nrep):
                rng = np.random.default_rng(10_000 + k)
                if mode == 'loop':
                    o = int(rng.integers(0, len(tilt))); ia = (o + np.arange(na)) % len(tilt); ig_ = (o * NSUB + np.arange(ng)) % len(g)
                else:
                    ia = blocks(rng, len(tilt), na, L); ig_ = (ia[:, None] * NSUB + np.arange(NSUB)[None, :]).reshape(-1)
                mt.append(tilt[ia].mean()); ig.append(((g[ig_] - bias) * DT).sum())
            res[mode] = dict(tilt_mean_rms=float(np.sqrt(np.mean(np.square(mt)))), gyro_int_rms=float(np.sqrt(np.mean(np.square(ig)))))
        res['theory'] = dict(tilt_mean_rms=float(np.sqrt(s2 / T)), gyro_int_rms=float(np.sqrt(lam * T)))
        out[name] = res
    return out


if __name__ == '__main__' and len(sys.argv) > 1 and sys.argv[1] == 'validate':
    import json
    for T in (15.0, 60.0):
        r = validate(T=T); print(f'T = {T:.0f} s (rad):', json.dumps(r, indent=None))


# ----------------------------------------------------------------------------------------------------------------- amendment 4
"""Amendment 4 (declared with amendment 3, after its first dev pass): the initial attitude is unknown. In amendments 2-3 the
plant rests at exactly 0 and every estimator starts at 0, so ignoring the tilt altogether is rewarded: the fixed-R Kalman
tuned on truth ran to R = 4.5e15 (gyro only) and won the low-throttle circumstances by never moving. A real flight
controller never knows the true attitude. Here the plant rests at theta0 ~ N(0, (0.5 deg)^2), drawn per flight from
default_rng([seed, 4242]); every estimator still starts at 0 with its own prior, sees tilt = theta0 + noise, and is scored
on (estimate - theta0). 0.5 deg is ten times the low-throttle error floor, small next to the vibration."""
SIG0 = 0.5 * D2R


def theta0(seed):
    return float(np.random.default_rng([seed, 4242]).normal(0.0, SIG0))


class Shift:
    def __init__(self, e, th0):
        self.e = e; self.n = e.n; self.th0 = th0
    def step(self, g, a, rl, t, reg, k):
        return np.asarray(self.e.step(g, None if a is None else a + self.th0, rl, t, reg, k), float) - self.th0
    def arm(self, b):
        if hasattr(self.e, 'arm'): self.e.arm(b)


def run4(seed, circ, ests, **kw):
    th0 = theta0(seed); return C.run(seed, circ, [Shift(e, th0) for e in ests], **kw)
