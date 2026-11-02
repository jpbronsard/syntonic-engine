"""Commissioning of the Kalman rivals on the block-bootstrap bench with unknown initial attitude (amendments 3-4), against truth, training seeds 100-102,
schedule C3, pooled J, bias state reset at arming (the strongest arming on amendment 2). Grids widened and recentred until
every parameter is interior or at its physical limit 0.
  KF_rival3  fixed R, Q_G, Q_B                                         (the corpus's tuned Kalman)
  KF_sched3  R = m rho s2_rec(level) / DTA (rho = 1 idle/low), Q_G, Q_B   (R scheduled on the throttle level)
  KF_full3   as KF_sched3 with Q_G = mq lam(level) also scheduled        (R and Q scheduled; the strongest)
Run: python3 rivals_bb.py  -> results/rivals_bb.json"""
import os, sys, json, itertools, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import bench_bb as B
B.install()
C, V = B.C, B.V
import law_flight as L
OUT = os.path.join(HERE, 'results'); os.makedirs(OUT, exist_ok=True)
VT, LI, QB0 = C.FX['var_tilt'], C.FX['lam_hat'], C.F.Q_B
LAMR = {r: V.LAM[r] for r in V.REGS}


class KFfull(L.KFsched):
    def __init__(self, m, rho_h, mq, qB):
        super().__init__(m, rho_h, np.zeros(len(np.atleast_1d(m))), qB); self.mq = np.asarray(mq, float)
    def step(self, g, a, reg_lam, t, reg, k):
        self.e.kQ = self.mq * LAMR[reg]
        return super().step(g, a, reg_lam, t, reg, k)


def make(kind, P):
    if kind == 'rival': return V.KF(VT * P[:, 0], LI * P[:, 1], QB0 * P[:, 2], 'reset')
    if kind == 'sched': return L.KFsched(P[:, 0], P[:, 1], LI * P[:, 2], QB0 * P[:, 3])
    return KFfull(P[:, 0], P[:, 1], P[:, 2], QB0 * P[:, 3])


def evaluate(kind, axes, seeds=(100, 101, 102)):
    P = np.array(list(itertools.product(*axes))); Js = []
    for s in seeds:
        Js.append(B.run4(s, 'C3', [make(kind, P)])['J'][0])
    J = C.pooled(np.array(Js)); i = int(np.argmin(J)); return P[i], float(J[i])


def edge(v, g):
    g = np.asarray(g); nz = g[g > 0]
    if v == 0: return 'at 0'
    return 'low edge' if v == nz.min() else ('high edge' if v == nz.max() else 'interior')


def around(v, half, step, zero):
    if v == 0: return None
    ax = v * 2.0 ** np.arange(-half, half + 1e-9, step)
    return np.concatenate([[0.0], ax]) if zero else ax


def tune(kind):
    Z = lambda lo, hi, st: np.concatenate([[0.0], 2.0 ** np.arange(lo, hi + 1e-9, st)])
    if kind == 'rival': axes = [2.0 ** np.arange(-2, 37, 2), Z(-24, 6, 2), Z(-16, 32, 2)]; zero = [False, True, True]
    elif kind == 'sched': axes = [2.0 ** np.arange(0, 31, 2), 2.0 ** np.arange(-8, 9, 2), Z(-24, 6, 3), Z(-8, 32, 4)]; zero = [False, False, True, True]
    else: axes = [2.0 ** np.arange(0, 31, 2), 2.0 ** np.arange(-8, 9, 2), Z(-10, 6, 2), Z(-8, 32, 4)]; zero = [False, False, True, True]
    p, j = evaluate(kind, axes); hist = [dict(p=p.tolist(), J=j)]
    for it in range(10):
        ax2 = []
        for d, v in enumerate(p):
            a = around(v, 2.0, 0.5, zero[d])
            ax2.append(axes[d] if a is None else a)
        p2, j2 = evaluate(kind, ax2); hist.append(dict(p=p2.tolist(), J=j2))
        e = [edge(v, a) for v, a in zip(p2, ax2)]; p = p2
        if all(x in ('interior', 'at 0') for x in e): break
    return dict(params=p.tolist(), J_train=j2, edges=e, refinements=it + 1, history=hist)


if __name__ == '__main__':
    res = {}
    for kind in (sys.argv[1:] or ['rival', 'sched', 'full']):
        t0 = time.time(); res[kind] = tune(kind); print(kind, {k: v for k, v in res[kind].items() if k != 'history'}, f'{time.time() - t0:.0f} s', flush=True)
        json.dump(res, open(os.path.join(OUT, 'rivals_r4_' + '_'.join(res) + '.json'), 'w'), indent=1)
