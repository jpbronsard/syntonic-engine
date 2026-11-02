# Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
# Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
"""Estimators and commissioning for the flight confirmation (2 Nov 2026). Imports the snapshot of the goal/flight harness
(common.py, variants*.py, results/facts2.json, results/rival.json) from this directory only.

KFsched    the bench's two-state Kalman with its measurement variance scheduled on the throttle level (known to the flight
           controller): R(level) = m * rho(level) * s2_rec(level) / DTA, rho = 1 at idle and low, rho_h at half; Q_G, Q_B;
           bias state reset at arming (the rival's arming). Optional bias table applied to the gyro after arming.
           Its m, rho_h, Q_G, Q_B are tuned against truth (commissioning, function tune below).
LevelLaw   the syntonic law with throttle-level memory: the Riccati recursion of the model (theta, b_idle, b_low, b_half),
           fed only by facts, nothing tuned. The bias of the active level drives the gyro integration and random-walks at
           q_b (idle rate-random-walk bound), the others are held. z = theta + noise, R = the live noise reading of the
           level (published Hadamard pipeline: block means 32/64/128 ms, 1 s EMA, min over scales of D_H l / 6), one memory
           per level initialised at the commissioning fact s2_rec(level), second differences straddling a throttle change
           skipped; Q_theta = q_theta (idle angle random walk). Pre-roll: (theta, b_idle) with the bench's priors
           (p11 = R, sigma_b = 0.5 deg/s). Arming: biases 0, their cross terms with theta dropped, bias covariance =
           (q_theta / T_pre) * ones (the shared error of the 5 s calibration) + the level prior on the diagonal:
             'btab' the logged per-throttle bias table applied to the gyro, prior = lam(level) / T_rec(level) (the variance of
                    the table entry: the record-mean of a gyro with angle random walk lam over T_rec);
             'mag'  no table applied, prior = lam/T_rec + (b_res(level) - b_res(idle))^2;
             'diff' no bias fact, prior = (0.5 deg/s)^2, the bench's own diffuse prior.
"""
import os, sys, json, itertools, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
if sys.path[0] != HERE: sys.path.insert(0, HERE)
import common as C
import variants as V
import variants2 as V2
G = C.G; DT, DTA, D2R, T_PRE = C.DT, C.DTA, C.D2R, C.T_PRE
S2 = {r: V.F2['curves']['raw'][r]['s2_rec'] for r in V.REGS}
TREC = {r: V.F2['curves']['raw'][r]['T_rec'] for r in V.REGS}
VB = {r: V.LAM[r] / TREC[r] for r in V.REGS}
LI = C.FX['lam_hat']; QB0 = C.F.Q_B; SB0 = (0.5 * D2R) ** 2
LV = {'idle': 1, 'low': 2, 'half': 3}
RV = json.load(open(os.path.join(C.OUT, 'rival.json')))['rival_params']


class KFsched(V.KF):
    def __init__(self, m, rho_h, kQ, qB, arming='reset', btab=False):
        m = np.asarray(m, float); rho_h = np.asarray(rho_h, float)
        self.Rtab = {'idle': m * S2['idle'] / DTA, 'low': m * S2['low'] / DTA, 'half': m * rho_h * S2['half'] / DTA}
        super().__init__(self.Rtab['idle'].copy(), np.asarray(kQ, float), np.asarray(qB, float), arming)
        self.btab = btab; self.armed = False
    def arm(self, bias_cal):
        super().arm(bias_cal); self.armed = True
    def step(self, g, a, reg_lam, t, reg, k):
        self.e.kR = self.Rtab[reg]
        if self.btab and self.armed: g = g - V.B_RES[reg]
        return super().step(g, a, reg_lam, t, reg, k)


class KFbtab(V.KF):
    """the fixed-R rival with the logged bias table applied to its gyro after arming"""
    def __init__(self, *a, **k):
        super().__init__(*a, **k); self.armed = False
    def arm(self, bias_cal):
        super().arm(bias_cal); self.armed = True
    def step(self, g, a, reg_lam, t, reg, k):
        return super().step(g - (V.B_RES[reg] if self.armed else 0.0), a, reg_lam, t, reg, k)


class LevelLaw:
    n = 1
    def __init__(self, prior='btab', live=True, qt=None, qb=None):
        self.qt = LI if qt is None else qt; self.qb = V.QB4 if qb is None else qb; self.prior = prior; self.live = live
        self.btab = prior == 'btab'; self.armed = False
        self.x = np.zeros(4); self.P = np.diag([S2['idle'] / DTA, SB0, SB0, SB0]); self.thr = 0.0
        self.had = V2.HadNoiseRegime(1); H = np.asarray(self.had.H, float)
        for r in V.REGS:
            self.had.DH[r] = np.zeros((len(H), 1)) + (6.0 * S2[r] / (H * G.DT))[:, None]; self.had.nb[r] = np.full(len(H), 2, dtype=int)
            self.had.seen.add(r)
    def arm(self, bias_cal):
        self.armed = True; self.x[1:] = 0.0; self.P[0, 1:] = 0.0; self.P[1:, 0] = 0.0
        def d(r):
            if self.prior == 'btab': return VB[r]
            if self.prior == 'mag': return VB[r] + (V.B_RES[r] - V.B_RES['idle']) ** 2
            return SB0
        Pb = np.full((3, 3), LI / T_PRE); Pb[1, 1] += d('low'); Pb[2, 2] += d('half'); self.P[1:, 1:] = Pb
    def step(self, g, a, reg_lam, t, reg, k):
        i = LV[reg]
        if self.btab and self.armed: g = g - V.B_RES[reg]
        self.thr += g * DT
        self.x[0] += (g - self.x[i]) * DT
        P = self.P; c = P[:, i].copy()                       # P <- F P F^T with F = I - DT e0 e_i^T
        P[0, :] -= DT * P[i, :]; P[:, 0] -= DT * c; P[0, 0] += DT * DT * P[i, i]
        P[0, 0] += self.qt * DT; P[i, i] += self.qb * DT
        if a is not None:
            s2 = self.had.update(np.array([a - self.thr]), reg) if self.live else None
            R = (S2[reg] if s2 is None else float(s2[0])) / DTA
            S = P[0, 0] + R; K = P[:, 0] / S; e = a - self.x[0]
            self.x += K * e; self.P = P - np.outer(K, P[0, :])
        return np.array([self.x[0]])


# ------------------------------------------------------------------------------------------------ commissioning
def tune_sched(btab, seeds=(100, 101, 102)):
    """KFsched tuned against truth exactly like the rival: schedule C3, training seeds, pooled J; coarse grid then
    half-step refinements recentred on the best point until no parameter sits on a grid edge (0, the physical limit of
    Q_G and Q_B, is always on the Q_B axis; Q_G at 0 stays at 0). Returns the parameter dict."""
    def grid(Mm, Rh, Qm, Bm):
        P = np.array(list(itertools.product(Mm, Rh, Qm, Bm))); Js = []
        for s in seeds:
            r = C.run(s, 'C3', [KFsched(P[:, 0], P[:, 1], LI * P[:, 2], QB0 * P[:, 3], btab=btab)]); Js.append(r['J'][0])
        J = C.pooled(np.array(Js)); i = int(np.argmin(J)); return P[i], float(J[i]), P, J
    Mm = 2.0 ** np.arange(0, 31, 2); Rh = 2.0 ** np.arange(-8, 5, 2); Qm = np.array([0.0, 2.0 ** -12, 2.0 ** -4])
    Bm = np.concatenate([[0.0], 2.0 ** np.arange(-8, 33, 4)])
    p, j, _, _ = grid(Mm, Rh, Qm, Bm); coarse = dict(m=p[0], rho_h=p[1], QG_mult=p[2], QB_mult=p[3], J=j)
    def edge(v, g):
        g = np.asarray(g); nz = g[g > 0]
        return 'at 0' if v == 0 else ('low edge' if v == nz.min() else 'high edge' if v == nz.max() else 'interior')
    for it in range(8):                                   # refine, recentring until no parameter sits on a grid edge
        lm, lr = np.log2(p[0]), np.log2(p[1])
        Mm2 = 2.0 ** np.arange(lm - 3, lm + 3.01, 0.5); Rh2 = 2.0 ** np.arange(lr - 2, lr + 2.01, 0.5)
        Qm2 = np.array([0.0]) if p[2] == 0 else p[2] * 2.0 ** np.arange(-4, 4.1, 2)
        Bm2 = Bm if p[3] == 0 else np.concatenate([[0.0], p[3] * 2.0 ** np.arange(-4, 4.01, 1)])
        p2, j2, _, _ = grid(Mm2, Rh2, Qm2, Bm2)
        edges = dict(m=edge(p2[0], Mm2), rho_h=edge(p2[1], Rh2), QG=edge(p2[2], Qm2) if p2[2] > 0 else 'at 0', QB=edge(p2[3], Bm2))
        p = p2
        if all(e in ('interior', 'at 0') for e in edges.values()): break
    return dict(m=float(p2[0]), rho_h=float(p2[1]), QG=float(LI * p2[2]), QB=float(QB0 * p2[3]), J_train=j2,
                coarse={k: float(v) for k, v in coarse.items()}, refinements=it + 1, edges=edges)


def kf_sched_from(d, btab):
    return KFsched([d['m']], [d['rho_h']], [d['QG']], [d['QB']], btab=btab)


if __name__ == '__main__' and len(sys.argv) > 1 and sys.argv[1] == 'commission':
    out = {}
    for name, bt in (('KF_sched', False), ('KF_sched_btab', True)):
        t0 = time.time(); out[name] = tune_sched(bt); print(name, out[name], f'{time.time() - t0:.0f} s', flush=True)
    json.dump(out, open(os.path.join(C.OUT, 'sched_tuned.json'), 'w'), indent=1)
