"""Round 3 (written after round 2's dev results, labelled post hoc). IC2: a robust form of the information-weighted
start-up cap. Round 2's IC summed DTA/s2 over every instantaneous live reading, first (single-block) readings included;
1/s2 is convex and the reading is noisy, so t_eff was inflated and in many seeds the cap vanished at the start.
IC2 counts each closed throttle visit with the noise level at its end (the EMA'd reading when the throttle changed) and
the current visit by its elapsed time:
    t_eff = s2_now * sum_closed(T_visit / s2_end_of_visit) + (t - t_start_of_current_visit),   cap = max(0.5 s, t_eff).
With no throttle change it equals the published cap max(0.5 s, t). No new constant.
Variants: NM + IC2 on the published law and on both two-state laws, with and without the logged bias table."""
import math, numpy as np
import common as C
import variants as V
import variants2 as V2


class _Cap2:
    def __init__(self, n):
        self.n = n; self.Ic = np.zeros(n); self.t0 = 0.0; self.reg = None; self.s2 = None
    def set_reg(self, t, reg):
        if reg != self.reg:
            if self.reg is not None:
                if self.s2 is not None: self.Ic = self.Ic + (t - self.t0) / self.s2
                else: self.Ic = None                      # a visit closed before any reading: fall back to t
            self.reg = reg; self.t0 = t
    def on_reading(self, s2, t): self.s2 = np.array(s2, float)
    def cap(self, t):
        if self.s2 is None or self.Ic is None: return max(0.5, t)
        return np.maximum(0.5, self.s2 * self.Ic + (t - self.t0))


class Law1R3(V2.Law1R):
    def __init__(self, rules, nm=True, btab=None):
        super().__init__(rules, nm=nm, ic=False, btab=btab); self.capr = _Cap2(self.n)
    def step(self, g, a, reg_lam, t, reg, k):
        self.capr.set_reg(t, reg); return super().step(g, a, reg_lam, t, reg, k)


class Law2SSR3(V2.Law2SSR):
    def __init__(self, rules, nm=True, bm=False, btab=None):
        super().__init__(rules, nm=nm, bm=bm, ic=False, btab=btab); self.capr = _Cap2(self.n)
    def step(self, g, a, reg_lam, t, reg, k):
        self.capr.set_reg(t, reg); return super().step(g, a, reg_lam, t, reg, k)


def build(names=None):
    pub = lambda s2, reg: math.sqrt(max(s2, 1e-300) / V.LAM[reg])
    b = lambda qb: (lambda s2, reg: V.gains_closed(max(s2, 1e-300), V.QT, qb))
    spec = [('L_nm_ic2', 1, pub, dict(), 'R3 published law + NM + IC2'),
            ('b4_nm_ic2', 2, b(V.QB4), dict(), 'R3 two-state law (q_b at T_rec/4) + NM + IC2'),
            ('b2_nm_ic2', 2, b(V.QB2), dict(), 'R3 two-state law (q_b at T_rec/2) + NM + IC2'),
            ('b4_nm_ic2_btab', 2, b(V.QB4), dict(btab=[True]), 'R3 two-state law (T_rec/4) + NM + IC2 + logged bias table'),
            ('b2_nm_ic2_btab', 2, b(V.QB2), dict(btab=[True]), 'R3 two-state law (T_rec/2) + NM + IC2 + logged bias table'),
            ('L_nm_ic2_btab', 1, pub, dict(btab=[True]), 'R3 published law + NM + IC2 + logged bias table')]
    spec = [s for s in spec if names is None or s[0] in names]
    ests, cols, D = [], [], {}
    for name, kind, rule, opts, desc in spec:
        ests.append(Law1R3([rule], **opts) if kind == 1 else Law2SSR3([rule], **opts)); cols.append(name); D[name] = desc
    if names is None or 'KF_rival' in names:
        e1, c1, d1, _ = V.build(['KF_rival']); ests += e1; cols += c1; D.update(d1)
    return ests, cols, D, None


class KFbtab(V.KF):
    """reference only: the rival with the logged per-throttle bias table applied to its gyro after arming"""
    def __init__(self, *a, **k):
        super().__init__(*a, **k); self.armed = False
    def arm(self, bias_cal):
        super().arm(bias_cal); self.armed = True
    def step(self, g, a, reg_lam, t, reg, k):
        return super().step(g - (V.B_RES[reg] if self.armed else 0.0), a, reg_lam, t, reg, k)


def build_ref(names=None):
    import json, os
    rv = json.load(open(os.path.join(C.OUT, 'rival.json')))['rival_params']
    ests = [KFbtab(np.array([rv['R']]), np.array([rv['QG']]), np.array([rv['QB']]), 'reset')]; cols = ['KF_rival_btab']
    D = {'KF_rival_btab': 'reference: the rival + the logged per-throttle bias table (same table as the *_btab laws)'}
    e1, c1, d1, _ = V.build(['KF_rival']); return ests + e1, cols + c1, {**D, **d1}, None
