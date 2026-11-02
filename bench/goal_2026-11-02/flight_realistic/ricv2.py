# Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
# Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
"""Vectorised Riccati form of the syntonic law, version 2 (development, 2 Nov 2026, amendments 3-4), written after the
first dev pass on amendment 4. Two options per column, both from the protocol or the facts, no new constant:
  sa  stationary alignment: the pre-roll is on the ground at rest by protocol (the gyro is being calibrated, its bias is
      unknown until arming), so before arming theta is estimated as a constant from the tilt alone (no gyro integration,
      the bias states untouched), the classical ground alignment of inertial navigation; at arming theta and P_theta are kept.
  rm  'live' R = live reading (version 1); 'fact' R = the record's s2_rec; 'max' R = max(fact, live): the commissioning
      fact is the floor and the live reading raises R when the noise rises.
Version 1 docstring follows.
The law's gains are the Riccati gains of the model the facts describe; nothing is tuned against truth. Columns differ by
structure only. State (theta, b_1..b_k): k = 3 holds one gyro bias per throttle level (idle, low, half; the active one
drives the integration and random-walks at q_b, the others are held), k = 1 a single bias.
  theta' = g - b_active (gyro after the arming calibration); P <- F P F^T + diag(q_theta(level) DT at theta, q_b DT at b_active)
  z = tilt = theta + noise, R = s2 / DTA with s2 = the live noise reading of the level (published Hadamard pipeline,
      per-level memory initialised at the record's fact s2_rec, straddling second differences skipped)
  q_theta(level) = lam(level) (the record's angle random walk at that throttle) if qt == 'level', else lam(idle)
  pre-roll: (theta, b_idle) with the bench's priors (p11 = R, sigma_b = 0.5 deg/s)
  arming: biases 0, cross terms with theta dropped, bias covariance = (lam_idle / T_pre) ones + diag(prior(level)):
      'cal'  0 (the calibration alone);  'vb' lam(level) / T_rec(level) (the record-mean bias variance);
      'mag'  vb + (b_res(level) - b_res(idle))^2 (the size of the logged change, sign unknown);  'diff' (0.5 deg/s)^2
  k = 1 with sw = True: at every throttle change the bias variance grows by the 'mag' amount of the change.
"""
import numpy as np
import bench_bb as B
C, V = B.C, B.V
import variants2 as V2
DT, DTA, D2R, T_PRE = C.DT, C.DTA, C.D2R, C.T_PRE
REGS = V.REGS; LV = {'idle': 0, 'low': 1, 'half': 2}
S2 = {r: V.F2['curves']['raw'][r]['s2_rec'] for r in REGS}
TREC = {r: V.F2['curves']['raw'][r]['T_rec'] for r in REGS}
LAM = {r: V.LAM[r] for r in REGS}; LI = LAM['idle']
VB = {r: LAM[r] / TREC[r] for r in REGS}
MAG = {r: VB[r] + (V.B_RES[r] - V.B_RES['idle']) ** 2 for r in REGS}
SB0 = (0.5 * D2R) ** 2


class RicV:
    def __init__(self, cols):
        """cols: list of dicts with k (1|3), qt ('level'|'idle'), prior ('cal'|'vb'|'mag'|'diff'), qb (float), sw (bool)"""
        self.cols = cols; n = self.n = len(cols); K = 4
        self.k = np.array([c['k'] for c in cols]); self.qb = np.array([c['qb'] for c in cols], float)
        self.qtl = np.array([c['qt'] == 'level' for c in cols]); self.sw = np.array([c.get('sw', False) for c in cols])
        self.sa = np.array([c.get('sa', False) for c in cols]); self.rm = [c.get('rm', 'live') for c in cols]
        self.x = np.zeros((n, K)); self.P = np.zeros((n, K, K))
        self.P[:, 0, 0] = S2['idle'] / DTA; self.P[:, 1, 1] = SB0; self.P[:, 2, 2] = SB0; self.P[:, 3, 3] = SB0
        self.armed = False; self.reg = None; self.thr = 0.0; self.ar = np.arange(n)
        self.had = V2.HadNoiseRegime(1); H = np.asarray(self.had.H, float)
        for r in REGS:
            self.had.DH[r] = np.zeros((len(H), 1)) + (6.0 * S2[r] / (H * C.G.DT))[:, None]; self.had.nb[r] = np.full(len(H), 2, dtype=int)
            self.had.seen.add(r)
    def _prior(self, c, r):
        p = c['prior']
        return 0.0 if p == 'cal' else VB[r] if p == 'vb' else MAG[r] if p == 'mag' else SB0
    def arm(self, bias_cal):
        self.armed = True; self.x[:, 1:] = 0.0; self.P[:, 0, 1:] = 0.0; self.P[:, 1:, 0] = 0.0
        for j, c in enumerate(self.cols):
            Pb = np.full((3, 3), LI / T_PRE)
            if c['k'] == 3:
                for r in ('low', 'half'): Pb[LV[r], LV[r]] += self._prior(c, r)
            else:
                Pb = np.zeros((3, 3)); Pb[0, 0] = LI / T_PRE
            self.P[j, 1:, 1:] = Pb
    def step(self, g, a, reg_lam, t, reg, k):
        if self.armed and self.reg is not None and reg != self.reg and self.sw.any():
            d = VB[reg] + VB[self.reg] + (V.B_RES[reg] - V.B_RES[self.reg]) ** 2
            self.P[self.sw & (self.k == 1), 1, 1] += d
        self.reg = reg
        i = np.where(self.k == 3, 1 + LV[reg], 1); ar = self.ar
        self.thr += g * DT
        hold = self.sa & (not self.armed)
        self.x[:, 0] += np.where(hold, 0.0, (g - self.x[ar, i]) * DT)
        P = self.P; P0 = P[hold].copy() if hold.any() else None; c = P[ar, :, i].copy(); Pii = P[ar, i, i].copy()
        P[:, 0, :] -= DT * P[ar, i, :]; P[:, :, 0] -= DT * c; P[:, 0, 0] += DT * DT * Pii
        P[:, 0, 0] += np.where(self.qtl, LAM[reg], LI) * DT; P[ar, i, i] += self.qb * DT
        if P0 is not None: P[hold] = P0
        if a is not None:
            s2 = self.had.update(np.array([a - self.thr]), reg)
            live = S2[reg] if s2 is None else float(s2[0])
            R = np.array([S2[reg] if m == 'fact' else max(S2[reg], live) if m == 'max' else live for m in self.rm]) / DTA
            S = P[:, 0, 0] + R; Kg = P[:, :, 0] / S[:, None]; e = a - self.x[:, 0]
            self.x += Kg * e[:, None]; self.P = P - Kg[:, :, None] * P[:, 0, :][:, None, :]
        return self.x[:, 0].copy()
