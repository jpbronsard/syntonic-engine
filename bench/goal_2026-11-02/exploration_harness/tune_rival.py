"""The rival: the bench's two-state Kalman (F.Est kind 1) tuned against truth on the training seeds 100-102, training
schedule (C3: half/low every 15 s), criterion = pooled estimation error. Grid over R, Q_G and Q_B, each including the
physical limit 0 for Q_G and Q_B, wide enough to be interior; then a refinement around the best point.
Three arming variants, all tuned the same way (the strongest is kept as the rival):
  'bench'  as F.Est: the bias state learned on the uncalibrated pre-roll gyro is kept after arming (double counting);
  'shift'  at arming b <- b - bias_cal (the state follows the calibrated gyro: consistent frame);
  'reset'  at arming b <- 0 (the arming calibration replaces the filter's own bias estimate).
Run: python tune_rival.py coarse | fine   -> results/rival_coarse.json, results/rival.json"""
import sys, time, itertools, numpy as np
import common as C
F = C.F
VT, LI, QB0 = C.FX['var_tilt'], C.FX['lam_hat'], F.Q_B


class KF(C.BenchEst):
    def __init__(self, kR, kQ, qB, arming):
        super().__init__(np.ones(len(kR), int), kR=kR, kQ=kQ, qB=qB); self.arming = arming
    def arm(self, bias_cal):
        if self.arming == 'shift': self.e.b = self.e.b - bias_cal
        elif self.arming == 'reset': self.e.b = np.zeros(self.n)


def grid_eval(Rm, Qm, Bm, seeds=C.TRAIN, circ='C3', armings=('bench', 'shift', 'reset')):
    """Rm, Qm, Bm: multiplier arrays (Qm, Bm may contain 0). Returns dict arming -> pooled J array (len grid)."""
    P = np.array(list(itertools.product(Rm, Qm, Bm)))
    kR, kQ, qB = VT * P[:, 0], LI * P[:, 1], QB0 * P[:, 2]
    Js = {a: [] for a in armings}
    for s in seeds:
        t0 = time.time(); ests = [KF(kR, kQ, qB, a) for a in armings]
        r = C.run(s, circ, ests)
        for a, J in zip(armings, r['J']): Js[a].append(J)
        print(f'  seed {s}: {len(P)} x {len(armings)} columns, {time.time() - t0:.0f} s', flush=True)
    return P, {a: C.pooled(v) for a, v in Js.items()}


def best_of(P, J, Rm, Qm, Bm):
    i = int(np.argmin(J)); p = P[i]
    def edge(v, g):
        g = np.asarray(g); nz = g[g > 0]
        if v == 0: return 'at 0 (physical limit)'
        return 'low edge' if v == nz.min() and not (g == 0).any() else ('high edge' if v == nz.max() else ('lowest nonzero (next to 0)' if v == nz.min() else 'interior'))
    return dict(i=i, R_mult=float(p[0]), QG_mult=float(p[1]), QB_mult=float(p[2]), R=float(VT * p[0]), QG=float(LI * p[1]),
                QB=float(QB0 * p[2]), J_train=float(J[i]),
                edges=dict(R=edge(p[0], Rm), QG=edge(p[1], Qm), QB=edge(p[2], Bm)))


if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'coarse':
        Rm = 2.0 ** np.arange(-2, 31, 2); Qm = np.concatenate([[0.0], 2.0 ** np.arange(-24, 7, 2)]); Bm = np.concatenate([[0.0], 2.0 ** np.arange(-16, 17, 2)])
        P, J = grid_eval(Rm, Qm, Bm)
        out = {a: best_of(P, J[a], Rm, Qm, Bm) for a in J}
        for a in J:
            o = out[a]; print(a, {k: v for k, v in o.items() if k != 'i'})
        C.dump('rival_coarse.json', dict(grid=dict(R_mult=Rm, QG_mult=Qm, QB_mult=Bm), best=out,
                                          J={a: J[a] for a in J}, P=P))
    if mode == 'coarse2':                     # widened after the coarse run put Q_B on its top edge
        Rm = 2.0 ** np.arange(4, 37, 2); Qm = np.concatenate([[0.0], 2.0 ** np.arange(-24, 7, 2)]); Bm = np.concatenate([[0.0], 2.0 ** np.arange(-8, 33, 2)])
        P, J = grid_eval(Rm, Qm, Bm, armings=('bench', 'reset'))
        out = {a: best_of(P, J[a], Rm, Qm, Bm) for a in J}
        for a in J:
            o = out[a]; print(a, {k: v for k, v in o.items() if k != 'i'})
        import json, os
        co = json.load(open(os.path.join(C.OUT, 'rival_coarse.json')))
        co['best_widened'] = out; co['grid_widened'] = dict(R_mult=Rm, QG_mult=Qm, QB_mult=Bm)
        co['best']['bench'] = out['bench']; co['best']['reset'] = out['reset']
        C.dump('rival_coarse.json', co)
    if mode == 'ridge':                       # 2-D map of the (R, Q_B) ridge at Q_G = 0, plus small Q_G > 0; final rival
        import json, os
        Rm = 2.0 ** np.arange(10, 22.01, 0.5); Qm = np.concatenate([[0.0], 2.0 ** np.arange(-20, -3, 4.0)]); Bm = np.concatenate([[0.0], 2.0 ** np.arange(10, 24.01, 0.5)])
        P, J = grid_eval(Rm, Qm, Bm, armings=('bench', 'reset'))
        out = {a: best_of(P, J[a], Rm, Qm, Bm) for a in J}
        for a in J:
            print(a, {k: v for k, v in out[a].items() if k != 'i'})
            j = J[a].reshape(len(Rm), len(Qm), len(Bm))[:, 0, :]
            print('   J at Q_G = 0, rows log2 R_mult, cols log2 Q_B_mult (every 2nd):')
            print('        ' + ' '.join(f'{x:6.1f}' for x in np.log2(Bm[1::2])))
            for i in range(0, len(Rm), 2): print(f'   {np.log2(Rm[i]):5.1f} ' + ' '.join(f'{x:6.4f}' for x in j[i, 1::2]) + f'   (Q_B=0: {j[i, 0]:.4f})')
        win = min(out, key=lambda a: out[a]['J_train'])
        C.dump('rival.json', dict(variants=out, rival=win, rival_params=out[win], grid=dict(R_mult=Rm, QG_mult=Qm, QB_mult=Bm),
                                  J={a: J[a] for a in J}, note='tuned against truth on training seeds 100-102, schedule C3; ridge map'))
        print('rival =', win, {k: v for k, v in out[win].items() if k != 'i'})
    if mode == 'fine':
        import json, os
        co = json.load(open(os.path.join(C.OUT, 'rival_coarse.json')))
        out = {}
        for a, b in co['best'].items():
            lr, lq, lb = np.log2(b['R_mult']), b['QG_mult'], b['QB_mult']
            Rm = 2.0 ** (lr + np.arange(-2, 2.01, 0.5))
            Qm = np.concatenate([[0.0], LQ := (2.0 ** (np.log2(lq) + np.arange(-2, 2.01, 0.5)) if lq > 0 else 2.0 ** np.arange(-30, -23, 1.0))])
            Bm = np.concatenate([[0.0], (2.0 ** (np.log2(lb) + np.arange(-2, 2.01, 0.5)) if lb > 0 else 2.0 ** np.arange(-22, -15, 1.0))])
            P, J = grid_eval(Rm, Qm, Bm, armings=(a,))
            out[a] = best_of(P, J[a], Rm, Qm, Bm); out[a]['grid'] = dict(R_mult=Rm, QG_mult=Qm, QB_mult=Bm)
            print(a, {k: v for k, v in out[a].items() if k not in ('i', 'grid')}, flush=True)
        win = min(out, key=lambda a: out[a]['J_train'])
        C.dump('rival.json', dict(variants=out, rival=win, rival_params=out[win], note='tuned against truth on training seeds 100-102, schedule C3'))
        print('rival =', win, out[win])
