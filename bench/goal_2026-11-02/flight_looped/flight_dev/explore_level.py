"""Development round 4, 2 Nov 2026 (not pre-registered, written after round 2's dev results; round 3 not run).
Round 2: one bias state that is loosened at every throttle change fixes C2 but re-learns the bias at every switch and
loses C3/C5 twice over. The physics: the gyro bias depends on the throttle level, a known input, and is the same at
each return to a level. Law with that model, facts only, nothing tuned: state (theta, b_idle, b_low, b_half), the bias
of the active level drives the gyro integration and random-walks at q_b (idle RRW bound), the others are held;
measurement z = theta + noise with R = the live per-level noise reading (published Hadamard pipeline, per-level memory
initialised at the commissioning fact s2_rec, straddling second differences skipped), Q_theta = q_theta (idle ARW).
Pre-roll: (theta, b_idle) with the bench's priors (p11 = R, sigma_b = 0.5 deg/s). Arming: all biases 0, cross terms
with theta dropped, bias covariance = (q_theta / T_pre) * ones (the shared calibration error) + diag of the level prior:
  'diff' no bias facts: (0.5 deg/s)^2, the bench's own diffuse bias prior (the level's bias is learned at first visit);
  'btab' the logged bias table applied to the gyro, prior = lam(level) / T_rec(level) (the variance of the table entry);
  'mag'  no table applied, prior = lam/T_rec + (b_res(level) - b_res(idle))^2 (the size of the change, sign unknown).
Run: python3 explore_level.py <circ> <seeds a-b>"""
import sys, os, json, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, '/home/claude/work/goal/flight'); sys.path.insert(0, HERE)
import common as C
import variants as V
import variants2 as V2
from tune_sched import KFsched, S2
from explore_riccati import SC, LI, QB0, RV
from explore_riccati2 import VB
G = C.G; DT, DTA, D2R = C.DT, C.DTA, C.D2R
LV = {'idle': 1, 'low': 2, 'half': 3}; SB0 = (0.5 * D2R) ** 2


class LevelLaw:
    n = 1
    def __init__(self, qt=None, qb=None, prior='diff', live=True):
        self.qt = LI if qt is None else qt; self.qb = V.QB4 if qb is None else qb; self.prior = prior; self.live = live
        self.btab = prior == 'btab'; self.armed = False
        self.x = np.zeros(4); self.P = np.diag([S2['idle'] / DTA, SB0, SB0, SB0]); self.thr = 0.0
        self.had = V2.HadNoiseRegime(1); H = np.asarray(self.had.H, float)
        for r in V.REGS:
            self.had.DH[r] = np.full((len(H), 1), 0.0) + (6.0 * S2[r] / (H * G.DT))[:, None]; self.had.nb[r] = np.full(len(H), 2, dtype=int)
            self.had.seen.add(r)
    def arm(self, bias_cal):
        self.armed = True; self.x[1:] = 0.0; self.P[0, 1:] = 0.0; self.P[1:, 0] = 0.0
        d = np.array([VB[r] if self.btab else (VB[r] + (V.B_RES[r] - V.B_RES['idle']) ** 2 if self.prior == 'mag' else SB0) for r in ('low', 'half')])
        Pb = np.full((3, 3), LI / C.T_PRE); Pb[1, 1] += d[0]; Pb[2, 2] += d[1]; self.P[1:, 1:] = Pb
    def step(self, g, a, reg_lam, t, reg, k):
        i = LV[reg]
        if self.btab and self.armed: g = g - V.B_RES[reg]
        self.thr += g * DT
        self.x[0] += (g - self.x[i]) * DT
        P = self.P; c = P[:, i].copy()                        # F = I - DT e0 e_i^T
        P[0, :] -= DT * P[i, :]; P[:, 0] -= DT * c; P[0, 0] += DT * DT * P[i, i]
        P[0, 0] += self.qt * DT; P[i, i] += self.qb * DT
        if a is not None:
            s2 = self.had.update(np.array([a - self.thr]), reg) if self.live else None
            R = (S2[reg] if s2 is None else float(s2[0])) / DTA
            S = P[0, 0] + R; K = P[:, 0] / S; e = a - self.x[0]
            self.x += K * e; self.P = P - np.outer(K, P[0, :])
        return np.array([self.x[0]])


def build():
    cols, ests = [], []
    ests.append(V.KF(np.array([RV['R']]), np.array([RV['QG']]), np.array([RV['QB']]), 'reset')); cols.append('KF_rival')
    ests.append(KFsched([SC['m']], [SC['rho_h']], [LI * SC['QG_mult']], [QB0 * SC['QB_mult']])); cols.append('KF_sched')
    ests.append(KFsched([SC['m']], [SC['rho_h']], [LI * SC['QG_mult']], [QB0 * SC['QB_mult']], btab=True)); cols.append('KF_sched_btab')
    for name, kw in (('lvl_mag', dict(prior='mag')), ('lvl_btab', dict(prior='btab')), ('lvl_btab_fact', dict(prior='btab', live=False)),
                     ('lvl_btab_qb2', dict(prior='btab', qb=V.QB2))):
        ests.append(LevelLaw(**kw)); cols.append(name)
    return ests, cols


if __name__ == '__main__':
    circ = sys.argv[1]; lo, hi = map(int, sys.argv[2].split('-')); out = os.path.join(HERE, 'lvl_runs.jsonl')
    for s in range(lo, hi + 1):
        t0 = time.time(); ests, cols = build(); r = C.run(s, circ, ests)
        J = np.concatenate(r['J']).tolist()
        with open(out, 'a') as f: f.write(json.dumps(dict(circ=circ, seed=s, cols=cols, J=J)) + '\n')
        print(circ, s, f'{time.time() - t0:.0f} s', ' '.join(f'{c} {j:.3f}' for c, j in zip(cols, J)), flush=True)
