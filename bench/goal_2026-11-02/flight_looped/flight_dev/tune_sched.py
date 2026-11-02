"""Development, 2 Nov 2026 (not pre-registered): a stronger rival for the flight confirmation.
KF_sched = the bench's two-state Kalman with its measurement variance scheduled on the throttle level, which the flight
controller knows: R(reg) = m * rho(reg) * s2_rec(reg) / DTA, s2_rec = the record's white-equivalent tilt noise (a fact),
rho = 1 for idle and low, rho_h for half (a free ratio so the truth decides the half/low weighting). Tuned exactly like
the rival (results/rival.json): against truth, training seeds 100-102, schedule C3, pooled J, reset arming, Q_G in {0, ...}.
Run: python3 tune_sched.py coarse | fine"""
import sys, os, json, itertools, time
import numpy as np
sys.path.insert(0, '/home/claude/work/goal/flight')
import common as C
import variants as V
OUT = os.path.dirname(os.path.abspath(__file__))
S2 = {r: V.F2['curves']['raw'][r]['s2_rec'] for r in V.REGS}
LI = C.FX['lam_hat']; QB0 = C.F.Q_B


class KFsched(V.KF):
    """R per throttle level, switched at every step (the bench's Kalman uses kR only in its update)."""
    def __init__(self, m, rho_h, kQ, qB, arming='reset', btab=False):
        m = np.asarray(m, float); rho_h = np.asarray(rho_h, float)
        self.Rtab = {'idle': m * S2['idle'] / C.DTA, 'low': m * S2['low'] / C.DTA, 'half': m * rho_h * S2['half'] / C.DTA}
        super().__init__(self.Rtab['idle'].copy(), np.asarray(kQ, float), np.asarray(qB, float), arming)
        self.btab = btab; self.armed = False
    def arm(self, bias_cal):
        super().arm(bias_cal); self.armed = True
    def step(self, g, a, reg_lam, t, reg, k):
        self.e.kR = self.Rtab[reg]
        if self.btab and self.armed: g = g - V.B_RES[reg]
        return super().step(g, a, reg_lam, t, reg, k)


def grid(Mm, Rh, Qm, Bm, seeds=C.TRAIN):
    P = np.array(list(itertools.product(Mm, Rh, Qm, Bm)))
    Js = []
    for s in seeds:
        t0 = time.time()
        r = C.run(s, 'C3', [KFsched(P[:, 0], P[:, 1], LI * P[:, 2], QB0 * P[:, 3])])
        Js.append(r['J'][0]); print(f'  seed {s}: {len(P)} columns {time.time() - t0:.0f} s', flush=True)
    return P, C.pooled(np.array(Js))


def report(P, J, axes):
    i = int(np.argmin(J)); p = P[i]
    out = dict(m=float(p[0]), rho_h=float(p[1]), QG_mult=float(p[2]), QB_mult=float(p[3]), J_train=float(J[i]), edges={})
    for k, (name, g) in enumerate(zip(('m', 'rho_h', 'QG_mult', 'QB_mult'), axes)):
        g = np.asarray(g); v = p[k]; nz = g[g > 0]
        out['edges'][name] = ('at 0' if v == 0 else 'low edge' if v == nz.min() and len(nz) > 1 else 'high edge' if v == nz.max() and len(nz) > 1 else 'interior')
    return out


if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'coarse':
        Mm = 2.0 ** np.arange(0, 31, 2); Rh = 2.0 ** np.arange(-8, 5, 2); Qm = np.array([0.0, 2.0 ** -12, 2.0 ** -4]); Bm = np.concatenate([[0.0], 2.0 ** np.arange(-8, 33, 4)])
        P, J = grid(Mm, Rh, Qm, Bm); o = report(P, J, (Mm, Rh, Qm, Bm)); print(o)
        json.dump(dict(best=o, axes=[x.tolist() for x in (Mm, Rh, Qm, Bm)]), open(os.path.join(OUT, 'sched_coarse.json'), 'w'), indent=1)
    if mode == 'fine':
        c = json.load(open(os.path.join(OUT, 'sched_coarse.json')))['best']
        lm, lr = np.log2(c['m']), np.log2(c['rho_h'])
        Mm = 2.0 ** np.arange(lm - 3, lm + 3.01, 0.5); Rh = 2.0 ** np.arange(lr - 2, lr + 2.01, 0.5)
        Qm = np.array([0.0]) if c['QG_mult'] == 0 else c['QG_mult'] * 2.0 ** np.arange(-4, 4.1, 2)
        Bm = (np.concatenate([[0.0], 2.0 ** np.arange(-8, 33, 4)]) if c['QB_mult'] == 0 else c['QB_mult'] * 2.0 ** np.arange(-4, 4.01, 1))
        P, J = grid(Mm, Rh, Qm, Bm); o = report(P, J, (Mm, Rh, Qm, Bm)); print(o)
        rv = json.load(open('/home/claude/work/goal/flight/results/rival.json'))['rival_params']
        o['R'] = {r: float(o['m'] * (o['rho_h'] if r == 'half' else 1.0) * S2[r] / C.DTA) for r in V.REGS}
        o['QG'] = float(LI * o['QG_mult']); o['QB'] = float(QB0 * o['QB_mult']); o['rival_J_train'] = rv['J_train']
        json.dump(dict(best=o, axes=[x.tolist() for x in (Mm, Rh, Qm, Bm)]), open(os.path.join(OUT, 'sched.json'), 'w'), indent=1)
