"""Development round 2, 2 Nov 2026 (not pre-registered, written after explore_riccati.py's dev results).
Diagnosis: the fact-based Riccati law loses C2/C3 to the scheduled Kalman because its bias prior is the idle calibration
alone (p22 = q_theta / T_pre), while the bias at another throttle is a different quantity, known only to the precision of
its own commissioning record. Fix, facts only: at every throttle change after arming, the bias variance grows by the
uncertainty of the bias change, var = lam(new) / T_rec(new) + lam(old) / T_rec(old) (the variance of each record-mean bias),
plus (b_res(new) - b_res(old))^2 when the bias table is not applied (the change is then unknown in sign).
Run: python3 explore_riccati2.py <circ> <seeds a-b>"""
import sys, os, json, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, '/home/claude/work/goal/flight'); sys.path.insert(0, HERE)
import common as C
import variants as V
from tune_sched import KFsched
from explore_riccati import RicLaw, SC, LI, QB0
TREC = {r: V.F2['curves']['raw'][r]['T_rec'] for r in V.REGS}
VB = {r: V.LAM[r] / TREC[r] for r in V.REGS}          # variance of each record-mean bias (rad/s)^2


class RicSw(RicLaw):
    def __init__(self, qG, qB, btab, mode):
        super().__init__(1.0, qG, qB, btab=btab); self.mode = mode; self.reg = None
    def step(self, g, a, reg_lam, t, reg, k):
        if self.armed and self.reg is not None and reg != self.reg:
            v = VB[reg] + (VB[self.reg] if self.mode == 'both' else 0.0)
            if not self.btab: v = v + (V.B_RES[reg] - V.B_RES[self.reg]) ** 2
            self.e.p22 = self.e.p22 + v
        self.reg = reg
        return super().step(g, a, reg_lam, t, reg, k)


def build():
    cols, ests = [], []
    ests.append(KFsched([SC['m']], [SC['rho_h']], [LI * SC['QG_mult']], [QB0 * SC['QB_mult']])); cols.append('KF_sched')
    ests.append(KFsched([SC['m']], [SC['rho_h']], [LI * SC['QG_mult']], [QB0 * SC['QB_mult']], btab=True)); cols.append('KF_sched_btab')
    for name, qg, qb, bt, mode in [('sw_qt_qb4', LI, V.QB4, False, 'both'), ('sw_qt_qb4_new', LI, V.QB4, False, 'new'),
                                   ('sw_0_qb4', 0.0, V.QB4, False, 'both'), ('sw_qt_qb4_btab', LI, V.QB4, True, 'both'),
                                   ('sw_qt_qb4_btab_new', LI, V.QB4, True, 'new'), ('sw_0_qb4_btab', 0.0, V.QB4, True, 'both')]:
        ests.append(RicSw([qg], [qb], bt, mode)); cols.append(name)
    return ests, cols


if __name__ == '__main__':
    circ = sys.argv[1]; lo, hi = map(int, sys.argv[2].split('-')); out = os.path.join(HERE, 'ric2_runs.jsonl')
    for s in range(lo, hi + 1):
        t0 = time.time(); ests, cols = build(); r = C.run(s, circ, ests)
        J = np.concatenate(r['J']).tolist()
        with open(out, 'a') as f: f.write(json.dumps(dict(circ=circ, seed=s, cols=cols, J=J)) + '\n')
        print(circ, s, f'{time.time() - t0:.0f} s', ' '.join(f'{c} {j:.3f}' for c, j in zip(cols, J)), flush=True)
