"""Round 2 (written after round 1's dev results, labelled post hoc). Three structural changes, no new constant:
NM  per-throttle memory of the live noise reading: the block-mean pipeline is shared, the EMA'd Hadamard values D_H are
    kept per throttle level; a second difference is used only if its three blocks all lie after the last throttle change
    (4 complete blocks since the change); first visit of a level: continue from the current values (no prior); later
    visits: resume that level's values (the throttle is a known input, like the per-throttle lambda).
BM  per-throttle memory of the bias state (two-state law only, after arming): on a throttle change the current bias
    estimate is stored for the level being left, and the stored estimate of the level entered is restored (first visit:
    continuity).
IC  information-weighted start-up cap: the published cap tau <= max(0.5 s, t) is the growing-memory window of an
    estimator that started with no knowledge at constant noise; with a varying noise level the same principle gives
    tau <= max(0.5 s, t_eff), t_eff = r_now * integral dt / r(t) (elapsed time counted in units of the current noise:
    the 1/n gain of a Kalman filter with Q_theta ~ 0 and varying R). r = the live reading s2.
"""
import math, numpy as np
import common as C
import variants as V
G = C.G
DT, DTA, T_TR, FLOOR_A, TAU_MAX, WARMUP = C.DT, C.DTA, V.T_TR, V.FLOOR_A, V.TAU_MAX, V.WARMUP


class HadNoiseRegime:
    """G.HadNoise with per-throttle EMA state (NM)."""
    def __init__(self, n):
        H = G.HSCALES; m = len(H); self.n = n; self.H = H
        self.acc = np.zeros((m, n)); self.cnt = np.zeros(m, dtype=int)
        self.p1 = np.full((m, n), np.nan); self.p2 = np.full((m, n), np.nan)
        self.a_blk = np.array([min(1.0, G.DT * s / max(G.EMA, 8 * G.DT * s)) for s in H])
        self.DH = {}; self.nb = {}; self.cur = None; self.since = np.zeros(m, dtype=int); self.seen = set()
    def update(self, z, reg):
        if reg != self.cur:
            if reg not in self.DH:
                if self.cur is None: self.DH[reg] = np.zeros((len(self.H), self.n)); self.nb[reg] = np.zeros(len(self.H), dtype=int)
                else: self.DH[reg] = self.DH[self.cur].copy(); self.nb[reg] = self.nb[self.cur].copy()
            revisit = reg in self.seen; self.seen.add(reg); self.cur = reg
            self.since[:] = 0 if revisit else 10 ** 6       # skip contaminated second differences on revisits only
        DH, nb = self.DH[reg], self.nb[reg]
        self.acc += z[None, :]; self.cnt += 1
        for i in np.flatnonzero(self.cnt == self.H):
            mean = self.acc[i] / self.H[i]; self.since[i] += 1
            if not np.isnan(self.p2[i, 0]) and self.since[i] >= 4:
                dd = mean - 2 * self.p1[i] + self.p2[i]; aw = self.a_blk[i] if nb[i] > 0 else 1.0
                DH[i] += aw * (dd * dd - DH[i]); nb[i] += 1
            self.p2[i] = self.p1[i]; self.p1[i] = mean; self.acc[i] = 0.0; self.cnt[i] = 0
        if (nb >= 2).all(): return np.min(DH * (self.H * G.DT)[:, None] / 6.0, axis=0)
        return None


class _Cap:
    """start-up cap: published max(0.5, t) or information-weighted max(0.5, s2_now * sum(DTA / s2))"""
    def __init__(self, n, ic):
        self.ic = ic; self.I = None; self.t_eff = None; self.n = n
    def on_reading(self, s2, t):
        if not self.ic: return
        if self.I is None: self.I = t / s2
        else: self.I = self.I + DTA / s2
        self.t_eff = s2 * self.I
    def cap(self, t):
        if self.ic and self.t_eff is not None: return np.maximum(0.5, self.t_eff)
        return max(0.5, t)


class Law1R(V.Law1):
    def __init__(self, rules, nm=False, ic=False, ceiling=None, btab=None):
        super().__init__(rules, ceiling=ceiling, btab=btab); self.nm = nm
        if nm: self.had = HadNoiseRegime(self.n)
        self.capr = _Cap(self.n, ic)
    def step(self, g, a, reg_lam, t, reg, k):
        gg = np.full(self.n, g)
        if self.armed: gg = gg - np.where(self.btab, V.B_RES[reg], 0.0)
        self.th = self.th + gg * DT; self.thr = self.thr + gg * DT
        tgt = np.clip(np.minimum(self.taul, self.capr.cap(t)), FLOOR_A, self.ceil)
        rr = 1.0 / self.tau; rr = rr + (DT / T_TR) * (1.0 / tgt - rr); self.tau = 1.0 / rr
        if a is not None:
            e = a - self.th
            s2 = self.had.update(a - self.thr, reg) if self.nm else self.had.update(a - self.thr)
            if s2 is not None and t > WARMUP:
                self.taul = np.array([self.rules[i](s2[i], reg) for i in range(self.n)]); self.capr.on_reading(s2, t)
            self.th = self.th + np.minimum(DTA / self.tau, 1.0) * e
        return self.th


class Law2SSR(V.Law2SS):
    def __init__(self, rules, nm=False, bm=False, ic=False, ceiling=None, btab=None):
        super().__init__(rules, ceiling=ceiling, btab=btab); self.nm = nm; self.bm = bm
        if nm: self.had = HadNoiseRegime(self.n)
        self.capr = _Cap(self.n, ic); self.bstore = {}; self.reg = None
    def step(self, g, a, reg_lam, t, reg, k):
        if self.bm and self.armed and reg != self.reg:
            if self.reg is not None: self.bstore[self.reg] = self.b.copy()
            if reg in self.bstore: self.b = self.bstore[reg].copy()
        self.reg = reg
        gg = np.full(self.n, g)
        if self.armed: gg = gg - np.where(self.btab, V.B_RES[reg], 0.0)
        self.th = self.th + (gg - self.b) * DT; self.thr = self.thr + gg * DT
        k1t = np.clip(np.maximum(self.k1l, 1.0 / self.capr.cap(t)), self.kmin, 1.0 / FLOOR_A)
        self.k1 = self.k1 + (DT / T_TR) * (k1t - self.k1)
        w2t = self.w2l if self.armed else np.zeros(self.n)
        self.w2 = self.w2 + (DT / T_TR) * (w2t - self.w2)
        if a is not None:
            e = a - self.th
            s2 = self.had.update(a - self.thr, reg) if self.nm else self.had.update(a - self.thr)
            if s2 is not None and t > WARMUP:
                kk = np.array([self.rules[i](s2[i], reg) for i in range(self.n)])
                self.k1l = kk[:, 0]; self.w2l = np.sqrt(kk[:, 1]); self.capr.on_reading(s2, t)
            k1, k2 = self.k1, self.w2 ** 2
            sq = np.sqrt((k1 * k1 - 4 * k2).astype(complex)); z1 = np.exp((-k1 + sq) / 2 * DTA); z2 = np.exp((-k1 - sq) / 2 * DTA)
            K1 = np.real(1 - z1 * z2); K2 = np.real((1 - z1) * (1 - z2)) / DTA
            self.th = self.th + K1 * e
            if self.armed: self.b = self.b - K2 * e
        return self.th
    def arm(self, bias_cal):
        super().arm(bias_cal); self.bstore = {}


def build(names=None):
    """round-2 registry: (objects, columns, descriptions). Each option set is its own object (own noise pipeline)."""
    pub = lambda s2, reg: math.sqrt(max(s2, 1e-300) / V.LAM[reg])
    b = lambda qb: (lambda s2, reg: V.gains_closed(max(s2, 1e-300), V.QT, qb))
    spec = []   # (name, kind, rule, opts, desc)
    spec += [('L_nm', 1, pub, dict(nm=True), 'R2 published law + NM'),
             ('L_ic', 1, pub, dict(ic=True), 'R2 published law + IC'),
             ('L_nm_ic', 1, pub, dict(nm=True, ic=True), 'R2 published law + NM + IC')]
    for tag, qb in (('b4', V.QB4), ('b2', V.QB2)):
        q = 'q_b at T_rec/4' if tag == 'b4' else 'q_b at T_rec/2'
        spec += [(f'{tag}_nm', 2, b(qb), dict(nm=True), f'R2 two-state law ({q}) + NM'),
                 (f'{tag}_nm_bm', 2, b(qb), dict(nm=True, bm=True), f'R2 two-state law ({q}) + NM + BM'),
                 (f'{tag}_ic', 2, b(qb), dict(ic=True), f'R2 two-state law ({q}) + IC'),
                 (f'{tag}_nm_ic', 2, b(qb), dict(nm=True, ic=True), f'R2 two-state law ({q}) + NM + IC'),
                 (f'{tag}_nm_bm_ic', 2, b(qb), dict(nm=True, bm=True, ic=True), f'R2 two-state law ({q}) + NM + BM + IC')]
    spec += [('b4_nm_bm_ic_btab', 2, b(V.QB4), dict(nm=True, bm=True, ic=True, btab=[True]), 'R2 two-state (T_rec/4) + NM + BM + IC + logged bias table'),
             ('b2_nm_bm_ic_btab', 2, b(V.QB2), dict(nm=True, bm=True, ic=True, btab=[True]), 'R2 two-state (T_rec/2) + NM + BM + IC + logged bias table')]
    spec = [s for s in spec if names is None or s[0] in names]
    ests, cols, D = [], [], {}
    for name, kind, rule, opts, desc in spec:
        if kind == 1: ests.append(Law1R([rule], **opts))
        else: ests.append(Law2SSR([rule], **opts))
        cols.append(name); D[name] = desc
    ref = ['KF_rival', 'LAW_c1', 'b_ss', 'b_ss_qb2']
    if names is not None: ref = [r for r in ref if r in names]
    if ref:
        e1, c1, d1, _ = V.build(ref); ests += e1; cols += c1; D.update(d1)
    return ests, cols, D, None
