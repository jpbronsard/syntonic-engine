# Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
# Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
"""The Riccati attitude law: a one-axis attitude estimator whose gains are the Riccati gains of the model that the
machine's measured facts describe. Nothing in it is tuned against truth.

In the package since 0.3.0 (Python only). It is the column LAW_R of the pre-registered confirmation 3 of 2 Nov 2026,
bench/goal_2026-11-02/flight_realistic/confirm_flight3.py, i.e. ricv2.RicV with the configuration
dict(k=3, qt='level', prior='vb', qb=q_b at T_rec/2, sa=True, rm='max'), rewritten as a self-contained single column whose
facts are constructor arguments. Inside the frozen harness it reproduces the confirmation's LAW_R values exactly
(bench/goal_2026-11-02/package_gate.py). Results and limits: bench/goal_2026-11-02/NOTE.md.

FACTS (measured once; constructor arguments)
  per level r of a known input that changes the noise (on a multicopter: the throttle level), from a record at rest:
    s2      the accelerometer tilt's white-noise density at that level (rad^2 s), the record's fact s2_rec
    lam     the gyro's angle random walk at that level (rad^2/s)
    T_rec   the length of that level's record (s); not used for the ground level
  q_b       the gyro bias random walk (rad^2/s^3), a bound from the ground record: 3 sigma_A^2(T) / T at T = T_rec / 2
  T_pre     the length of the arming calibration on the ground (s)
  dt, dt_accel  the gyro and accelerometer sample periods (s)
  The first level (or `ground`) is the one the vehicle is armed at, on the ground, at rest.

MODEL. State x = (theta, b_1, ..., b_L): the angle and one gyro bias per level. The active level's bias drives the
integration and random-walks at q_b; the other biases are held.
  predict, every gyro sample g:   theta += (g - b_active) dt;   P <- F P F^T + diag(lam(level) dt at theta, q_b dt at b_active)
  update, every accelerometer sample a (the tilt):   z = a = theta + noise,   R = max(s2(level), s2_live) / dt_accel
  s2_live is the published noise reading (bench/software_2026-10-26/step6g_real_noise.py): second differences of the means
  of adjacent blocks of (a - integral of g) at 8, 16 and 32 accelerometer samples, EMA'd over `ema` = 1 s
  (weight min(1, dt_accel s / max(ema, 8 dt_accel s)) for a block of s samples), each value divided back to a density
  (D_H x s dt_accel / 6) and the minimum over the three scales taken. The reading keeps one memory per level, initialised
  at the fact (D_H = 6 s2 / (s dt_accel)); after every change of level (and at the start) a second difference is used
  only once `settle` = 4 blocks of its scale have completed since the change, so that its three blocks all lie after it.
  The fact is the floor of R and the live reading raises R when the noise rises.

PROTOCOL
  before arming (on the ground, at rest, while the flight controller calibrates the gyro): stationary alignment, the
  classical ground alignment of inertial navigation. theta is a constant estimated from the tilt alone, with prior
  variance s2(ground) / dt_accel; the gyro is not integrated into theta and the bias states are not used.
  arm(): biases 0, no cross terms with theta, bias covariance (lam(ground) / T_pre) x ones (the calibration's own
  error, common to every level) + diag(lam(r) / T_rec(r)) on every level other than the ground one (the record-mean
  bias variance).
  after arming: g is the gyro rate with the arming calibration (the pre-roll gyro mean) subtracted, as the flight
  controller does. Before arming, g is the raw rate; it enters only the noise reading.

Usage:
    law = RiccatiAttitudeLaw({'idle': dict(s2=..., lam=...), 'low': dict(s2=..., lam=..., T_rec=...),
                              'half': dict(s2=..., lam=..., T_rec=...)}, q_b=..., T_pre=5.0)
    every gyro sample:   theta = law.step(g, a, level='low')       # a = the tilt (rad) on accelerometer samples, else None
    at the end of the ground calibration:   law.arm()
  `step(g, a, reg_lam, t, level, k)` is also the call of the confirmation harness (reg_lam, t and k are not used).
"""
import numpy as np


class RiccatiAttitudeLaw:
    """One-axis Riccati attitude law with one gyro bias per level, fed by facts only (see the module docstring)."""

    n = 1   # one estimator column (the confirmation harness reads it)

    def __init__(self, levels, q_b, T_pre, dt=1e-3, dt_accel=4e-3, ground=None, block_scales=(8, 16, 32), ema=1.0,
                 settle=4):
        names = list(levels)
        if not names:
            raise ValueError('levels: at least one level is needed')
        self.ground = names[0] if ground is None else ground
        if self.ground not in levels:
            raise ValueError(f'ground level {self.ground!r} is not among the levels')
        self.levels = names
        self.index = {r: j for j, r in enumerate(names)}           # the bias of level r is x[1 + index[r]]
        self.s2, self.lam, self.T_rec = {}, {}, {}
        for r in names:
            f = levels[r]
            for key in ('s2', 'lam'):
                if key not in f:
                    raise ValueError(f'level {r!r}: the fact {key!r} is missing')
            if r != self.ground and 'T_rec' not in f:
                raise ValueError(f'level {r!r}: the fact \'T_rec\' is missing')
            self.s2[r] = float(f['s2']); self.lam[r] = float(f['lam'])
            self.T_rec[r] = float(f['T_rec']) if 'T_rec' in f else None
        self.q_b = float(q_b); self.T_pre = float(T_pre); self.dt = float(dt); self.dta = float(dt_accel)
        # live noise reading: one shared block pipeline, one memory of D_H per level
        self.H = np.asarray(block_scales, dtype=np.int64); m = len(self.H)
        self.hd = self.H * self.dta                                  # block lengths (s)
        self.a_blk = np.array([min(1.0, self.dta * s / max(ema, 8 * self.dta * s)) for s in self.H])
        self.settle = int(settle)
        self.acc = np.zeros(m); self.cnt = np.zeros(m, dtype=np.int64)
        self.p1 = np.full(m, np.nan); self.p2 = np.full(m, np.nan)
        self.since = np.zeros(m, dtype=np.int64); self.cur = None
        self.DH = {r: 6.0 * self.s2[r] / (self.H * self.dta) for r in names}
        self.thr = 0.0                                               # integral of the gyro, for the noise reading
        # estimator
        self.x = np.zeros(1 + len(names)); self.P = np.zeros((1 + len(names), 1 + len(names)))
        self.P[0, 0] = self.s2[self.ground] / self.dta
        self.armed = False

    @property
    def theta(self):
        """the angle estimate (rad)"""
        return float(self.x[0])

    def bias(self, level):
        """the gyro-bias estimate of a level (rad/s), relative to the arming calibration"""
        return float(self.x[1 + self.index[level]])

    def noise_reading(self, level):
        """the live noise reading of a level (rad^2 s), as the next update would floor it with the fact"""
        return float(np.min(self.DH[level] * self.hd / 6.0))

    def _reading(self, z, level):
        if level != self.cur:
            self.cur = level; self.since[:] = 0
        DH = self.DH[level]
        self.acc += z; self.cnt += 1
        for i in np.flatnonzero(self.cnt == self.H):
            mean = self.acc[i] / self.H[i]; self.since[i] += 1
            if not np.isnan(self.p2[i]) and self.since[i] >= self.settle:
                dd = mean - 2 * self.p1[i] + self.p2[i]
                DH[i] += self.a_blk[i] * (dd * dd - DH[i])
            self.p2[i] = self.p1[i]; self.p1[i] = mean; self.acc[i] = 0.0; self.cnt[i] = 0
        return float(np.min(DH * self.hd / 6.0))

    def arm(self, bias_cal=None):
        """end of the ground calibration. bias_cal (the calibration itself) is not used: the caller subtracts it from g."""
        self.armed = True
        self.x[1:] = 0.0; self.P[0, 1:] = 0.0; self.P[1:, 0] = 0.0
        Pb = np.full((len(self.levels), len(self.levels)), self.lam[self.ground] / self.T_pre)
        for r in self.levels:
            if r != self.ground:
                j = self.index[r]; Pb[j, j] += self.lam[r] / self.T_rec[r]
        self.P[1:, 1:] = Pb

    def step(self, g, a, reg_lam=None, t=None, level=None, k=None):
        """one gyro sample g (rad/s); a = the tilt (rad) when an accelerometer sample arrives, else None; level = the
        current level. Returns the angle estimate (rad)."""
        if level is None:
            raise ValueError('step needs the current level')
        i = 1 + self.index[level]; dt = self.dt
        self.thr += g * dt
        if self.armed:
            x, P = self.x, self.P
            x[0] += (g - x[i]) * dt
            c = P[:, i].copy(); Pii = P[i, i]
            P[0, :] -= dt * P[i, :]
            P[:, 0] -= dt * c
            P[0, 0] += dt * dt * Pii
            P[0, 0] += self.lam[level] * dt
            P[i, i] += self.q_b * dt
        if a is not None:
            live = self._reading(a - self.thr, level)
            R = max(self.s2[level], live) / self.dta
            P = self.P
            S = P[0, 0] + R; Kg = P[:, 0] / S; e = a - self.x[0]
            self.x += Kg * e
            self.P = P - Kg[:, None] * P[0, :][None, :]
        return float(self.x[0])
