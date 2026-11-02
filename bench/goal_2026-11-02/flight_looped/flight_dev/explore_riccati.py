"""Development, 2 Nov 2026 (not pre-registered). The throttle-scheduled Kalman tuned on truth (tune_sched.py) reaches
J_train 0.090 deg on C3 against 0.239 for the fixed-R rival: it carries the precision gathered at low throttle through
the half-throttle phases. Candidate law with that memory: the Riccati recursion of the two-state model fed only by facts,
R(reg) = s2_rec(reg) / DTA (the noise measured per throttle), Q_G = q_theta (idle ARW), Q_B = q_b (idle RRW bound);
at arming b = 0, p12 = 0, p22 = q_theta / T_pre (the variance of the 5 s calibration). Nothing tuned.
Run: python3 explore_riccati.py <circ> <seeds a-b>   (training/dev seeds only)"""
import sys, os, json, time
import numpy as np
sys.path.insert(0, '/home/claude/work/goal/flight'); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C
import variants as V
from tune_sched import KFsched, S2, LI, QB0
SC = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'sched_coarse.json')))['best']
RV = json.load(open('/home/claude/work/goal/flight/results/rival.json'))['rival_params']


class RicLaw(KFsched):
    """facts only; arming sets b = 0, p12 = 0, p22 = q_theta / T_pre (+ dB2: extra bias variance per column)"""
    def __init__(self, m, qG, qB, btab=False, p22arm=None):
        n = len(np.atleast_1d(qB)); super().__init__(np.broadcast_to(m, (n,)), np.ones(n), qG, qB, 'reset', btab)
        self.p22arm = np.broadcast_to(LI / C.T_PRE if p22arm is None else p22arm, (n,)).astype(float)
    def arm(self, bias_cal):
        super().arm(bias_cal); self.e.p12 = np.zeros(self.e.n); self.e.p22 = self.p22arm.copy()


def build():
    QB4, QB2 = V.QB4, V.QB2; dB = abs(V.B_RES['half'] - V.B_RES['idle'])
    cols, ests = [], []
    ests.append(V.KF(np.array([RV['R']]), np.array([RV['QG']]), np.array([RV['QB']]), 'reset')); cols.append('KF_rival')
    ests.append(KFsched([SC['m']], [SC['rho_h']], [LI * SC['QG_mult']], [QB0 * SC['QB_mult']])); cols.append('KF_sched')
    ests.append(KFsched([SC['m']], [SC['rho_h']], [LI * SC['QG_mult']], [QB0 * SC['QB_mult']], btab=True)); cols.append('KF_sched_btab')
    spec = [('ric_qt_qb4', LI, QB4, False, None), ('ric_qt_qb2', LI, QB2, False, None), ('ric_0_qb4', 0.0, QB4, False, None),
            ('ric_qt_qb4_db', LI, QB4, False, LI / C.T_PRE + dB ** 2), ('ric_0_qb4_db', 0.0, QB4, False, LI / C.T_PRE + dB ** 2),
            ('ric_qt_qb4_btab', LI, QB4, True, None), ('ric_0_qb4_btab', 0.0, QB4, True, None), ('ric_qt_qb2_btab', LI, QB2, True, None)]
    for name, qg, qb, bt, p in spec:
        ests.append(RicLaw(1.0, [qg], [qb], btab=bt, p22arm=p)); cols.append(name)
    return ests, cols


if __name__ == '__main__':
    circ = sys.argv[1]; lo, hi = map(int, sys.argv[2].split('-')); out = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'ric_runs.jsonl')
    for s in range(lo, hi + 1):
        t0 = time.time(); ests, cols = build(); r = C.run(s, circ, ests)
        J = np.concatenate(r['J']).tolist()
        with open(out, 'a') as f: f.write(json.dumps(dict(circ=circ, seed=s, cols=cols, J=J)) + '\n')
        print(circ, s, f'{time.time() - t0:.0f} s', ' '.join(f'{c} {j:.3f}' for c, j in zip(cols, J)), flush=True)
