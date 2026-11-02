"""Development round 3, 2 Nov 2026 (not pre-registered, written after round 2 was launched, before its results).
The fact-based Riccati law cannot see an unannounced change of the noise (C4: tilt noise x sqrt(10) at 20-40 s), so it
trusts noisy data and loses to the fixed-R rival there. Fix: R = the live noise reading per throttle level (the published
block-mean Hadamard pipeline with per-throttle memory, variants2.HadNoiseRegime), each level's memory initialised at the
commissioning fact s2_rec(level) and every throttle change treated as a revisit (straddling second differences skipped).
Everything else as explore_riccati2.RicSw. Run: python3 explore_riccati3.py <circ> <seeds a-b>"""
import sys, os, json, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, '/home/claude/work/goal/flight'); sys.path.insert(0, HERE)
import common as C
import variants as V
import variants2 as V2
from tune_sched import KFsched, S2
from explore_riccati import RicLaw, SC, LI, QB0, RV
from explore_riccati2 import RicSw, VB
G = C.G


class RicLive(RicSw):
    def __init__(self, qG, qB, btab, mode='both'):
        super().__init__(qG, qB, btab, mode); n = self.e.n
        self.had = V2.HadNoiseRegime(n); H = np.asarray(self.had.H, float)
        for r in V.REGS:
            self.had.DH[r] = np.tile((6.0 * S2[r] / (H * G.DT))[:, None], (1, n)); self.had.nb[r] = np.full(len(H), 2, dtype=int)
            self.had.seen.add(r)
        self.thr = np.zeros(n)
    def step(self, g, a, reg_lam, t, reg, k):
        if self.armed and self.reg is not None and reg != self.reg:
            v = VB[reg] + (VB[self.reg] if self.mode == 'both' else 0.0)
            if not self.btab: v = v + (V.B_RES[reg] - V.B_RES[self.reg]) ** 2
            self.e.p22 = self.e.p22 + v
        self.reg = reg
        if self.btab and self.armed: g = g - V.B_RES[reg]
        self.thr = self.thr + g * C.DT
        if a is not None:
            s2 = self.had.update(np.full(self.e.n, a) - self.thr, reg)
            self.e.kR = (np.full(self.e.n, S2[reg]) if s2 is None else np.asarray(s2, float)) / C.DTA
        return C.BenchEst.step(self, g, a, reg_lam, t, reg, k)


def build():
    cols, ests = [], []
    ests.append(V.KF(np.array([RV['R']]), np.array([RV['QG']]), np.array([RV['QB']]), 'reset')); cols.append('KF_rival')
    ests.append(KFsched([SC['m']], [SC['rho_h']], [LI * SC['QG_mult']], [QB0 * SC['QB_mult']])); cols.append('KF_sched')
    for name, bt in (('live', False), ('live_btab', True)):
        ests.append(RicLive([LI], [V.QB4], bt)); cols.append(name)
    return ests, cols


if __name__ == '__main__':
    circ = sys.argv[1]; lo, hi = map(int, sys.argv[2].split('-')); out = os.path.join(HERE, 'ric3_runs.jsonl')
    for s in range(lo, hi + 1):
        t0 = time.time(); ests, cols = build(); r = C.run(s, circ, ests)
        J = np.concatenate(r['J']).tolist()
        with open(out, 'a') as f: f.write(json.dumps(dict(circ=circ, seed=s, cols=cols, J=J)) + '\n')
        print(circ, s, f'{time.time() - t0:.0f} s', ' '.join(f'{c} {j:.3f}' for c, j in zip(cols, J)), flush=True)
