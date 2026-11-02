"""Estimator variants (vectorised, plug into common.run), 2 Nov 2026, exploratory. Every constant comes from the
records' noise facts (results/facts2.json, results/facts_allan.json) or from the published engine; none from truth.

Law1   one-state law, the published engine exactly (start-up cap tau <= max(0.5 s, t), 1/tau travel at 128 ms, floor
       DTA/0.45, ceiling 30 s unless disabled, warm-up 1 s, bias loop off), with a pluggable target rule tau(s2, regime).
Law2SS two-state law, steady-state Kalman-Bucy gains k2 = sqrt(q_b/r), k1 = sqrt(q_theta/r + 2 k2) (rule-dependent r),
       pole-matched discretisation at the accelerometer rate (equal to the discrete ARE gains, facts2.py), start-up cap on
       the angle loop (k1 >= 1/max(0.5 s, t)), k1 and sqrt(k2) travel at 128 ms, floor, optional ceiling; bias loop from
       arming with b = 0 (the arming calibration is the bias estimate at that instant).
Law2TV two-state Kalman with the Riccati recursion run online: R = s2_live / DTA (live), Q_theta, Q_b from the idle record,
       P0 by the bench's convention (p11 = R at the first update); the bias state starts at arming (b = 0,
       p22 = q_theta / T_pre: the calibration's own variance, from the gyro's angle random walk over the 5 s pre-roll).
Optional bias table (all classes): after arming the gyro is corrected by the logged residual bias of the throttle level
(per-record mean minus the idle mean), a fact measured once per throttle on the held frame.
"""
import os, json, math
import numpy as np
import common as C
F, G = C.F, C.G
D2R, DT, DTA, T_PRE = C.D2R, C.DT, C.DTA, C.T_PRE
FLOOR_A, TAU_MAX, T_TR, WARMUP = F.FLOOR_A, F.TAU_MAX, F.T_TR, G.WARMUP
F2 = json.load(open(os.path.join(C.OUT, 'facts2.json')))
TG = np.array(F2['TG'])
QT, QB4, QB2 = F2['q_theta'], F2['q_b'], F2['q_b_half']
LAM = {r: C.FX['lam'][r] for r in ('idle', 'low', 'half')}
B_RES = {r: F2['curves']['raw'][r]['b_res'] for r in ('idle', 'low', 'half')}     # rad/s
REGS = ('idle', 'low', 'half')
S2G = 10.0 ** np.arange(-9.0, 0.5, 0.01)            # lookup grid of live s2 values (rad^2 s)


# ------------------------------------------------------------------------------------------------ rule tables
def _curve(key, reg):
    cv = F2['curves'][key][reg]; return np.array(cv['r']), np.array(cv['lam']), np.array(cv['sz']), cv['s2_rec'], cv['b_res']


def _interp_loglog(T, x, y, b_res=None, kind='flat'):
    """interpolate y(x) log-log at T; beyond the grid: flat (r, sz-like) or lam-like (last + b^2 (T - T_last))"""
    lx, ly = np.log(x), np.log(np.maximum(y, 1e-300))
    v = np.exp(np.interp(np.log(T), lx, ly))
    if kind == 'lam':
        hi = T > x[-1]; v = np.where(hi, y[-1] + (b_res ** 2) * (T - x[-1]), v)
    return v


def table_a_fp(key):
    """(a) fixed point tau = sqrt(s2 rho(2 tau) / lam(2 tau)), rho = r / s2_rec; returns {reg: tau(S2G)}"""
    out = {}
    taus = np.exp(np.linspace(np.log(FLOOR_A), np.log(3000.0), 3000))
    for reg in REGS:
        r, lam, _, s2r, b = _curve(key, reg)
        rho = _interp_loglog(2 * taus, TG, r) / s2r; lm = _interp_loglog(2 * taus, TG, lam, b, 'lam')
        res = np.empty(len(S2G))
        for i, s2 in enumerate(S2G):
            g = taus ** 2 * lm - s2 * rho
            j = np.flatnonzero(g >= 0)
            res[i] = taus[j[0]] if len(j) else taus[-1]
        out[reg] = res
    return out


def table_a_min(key):
    """(a) classical crossover: T* = argmin_T sigma_z^2(T; s2), sigma_z^2(T; s2) = sz_rec(T) + (s2/s2_rec - 1) r(T)/T;
    tau = T* / sqrt(3) (canonical class: sigma_z^2 = r/T + q T/3 -> T* = sqrt(3 r / q), tau_opt = sqrt(r / q))"""
    out = {}
    Tf = np.exp(np.linspace(np.log(TG[0]), np.log(TG[-1]), 1200))
    for reg in REGS:
        r, _, sz, s2r, b = _curve(key, reg)
        szf = _interp_loglog(Tf, TG, sz); rf = _interp_loglog(Tf, TG, r)
        res = np.empty(len(S2G))
        for i, s2 in enumerate(S2G):
            v = szf + (s2 / s2r - 1.0) * rf / Tf; res[i] = Tf[int(np.argmin(v))] / math.sqrt(3.0)
        out[reg] = res
    return out


def gains_closed(r, qt, qb):
    k2 = np.sqrt(qb / r); return np.sqrt(qt / r + 2 * k2), k2


def table_c(key, qb):
    """(c) two-state with r at the operating scale: r_op = s2 rho(2/k1), k1 = k1(r_op) fixed point (bisection in log k1:
    h(k1) = k1 - k1(s2 rho(2/k1)) is increasing when rho does not increase with scale faster than ~T^2)"""
    out = {}
    for reg in REGS:
        r, _, _, s2r, _ = _curve(key, reg)
        K1 = np.empty(len(S2G)); K2 = np.empty(len(S2G))
        for i, s2 in enumerate(S2G):
            lo, hi = np.log(1e-5), np.log(1.0 / FLOOR_A)
            for _ in range(60):
                mid = 0.5 * (lo + hi); k = math.exp(mid)
                rop = s2 * float(_interp_loglog(np.array([2.0 / k]), TG, r)[0]) / s2r
                h = k - gains_closed(rop, QT, qb)[0]
                if h > 0: hi = mid
                else: lo = mid
            k = math.exp(0.5 * (lo + hi)); rop = s2 * float(_interp_loglog(np.array([2.0 / k]), TG, r)[0]) / s2r
            K1[i], K2[i] = gains_closed(rop, QT, qb)
        out[reg] = (K1, K2)
    return out


def lookup(tab, reg, s2):
    return np.exp(np.interp(np.log(np.maximum(s2, 1e-300)), np.log(S2G), np.log(tab[reg])))


# ------------------------------------------------------------------------------------------------ estimators
class Law1:
    """one-state law. rules: list (one per column) of callables rule(s2_col, reg) -> target window (s)."""
    def __init__(self, rules, ceiling=None, btab=None):
        self.n = n = len(rules); self.rules = rules
        self.ceil = np.array([TAU_MAX if c is None else c for c in (ceiling or [None] * n)], float)
        self.btab = np.array(btab if btab is not None else [False] * n, bool)
        self.th = np.zeros(n); self.thr = np.zeros(n); self.tau = np.full(n, 0.5); self.taul = np.full(n, 0.5)
        self.had = G.HadNoise(n); self.armed = False
    def arm(self, bias_cal): self.armed = True
    def step(self, g, a, reg_lam, t, reg, k):
        gg = np.full(self.n, g)
        if self.armed: gg = gg - np.where(self.btab, B_RES[reg], 0.0)
        self.th = self.th + gg * DT; self.thr = self.thr + gg * DT
        cap = max(0.5, t)
        tgt = np.clip(np.minimum(self.taul, cap), FLOOR_A, self.ceil)
        rr = 1.0 / self.tau; rr = rr + (DT / T_TR) * (1.0 / tgt - rr); self.tau = 1.0 / rr
        if a is not None:
            e = a - self.th
            s2 = self.had.update(a - self.thr)
            if s2 is not None and t > WARMUP:
                self.taul = np.array([self.rules[i](s2[i], reg) for i in range(self.n)])
            self.th = self.th + np.minimum(DTA / self.tau, 1.0) * e
        return self.th


class Law2SS:
    """two-state law, steady-state gains. rules: callables rule(s2_col, reg) -> (k1, k2) continuous gains."""
    def __init__(self, rules, ceiling=None, btab=None):
        self.n = n = len(rules); self.rules = rules
        self.kmin = np.array([0.0 if c is None else 1.0 / c for c in (ceiling or [None] * n)], float)
        self.btab = np.array(btab if btab is not None else [False] * n, bool)
        self.th = np.zeros(n); self.b = np.zeros(n); self.thr = np.zeros(n)
        self.k1 = np.full(n, 2.0); self.w2 = np.zeros(n); self.k1l = np.full(n, 2.0); self.w2l = np.zeros(n)
        self.had = G.HadNoise(n); self.armed = False
    def arm(self, bias_cal):
        self.armed = True; self.b = np.zeros(self.n)
    def step(self, g, a, reg_lam, t, reg, k):
        gg = np.full(self.n, g)
        if self.armed: gg = gg - np.where(self.btab, B_RES[reg], 0.0)
        self.th = self.th + (gg - self.b) * DT; self.thr = self.thr + gg * DT
        cap = max(0.5, t)
        k1t = np.clip(np.maximum(self.k1l, 1.0 / cap), self.kmin, 1.0 / FLOOR_A)
        self.k1 = self.k1 + (DT / T_TR) * (k1t - self.k1)
        w2t = self.w2l if self.armed else np.zeros(self.n)
        self.w2 = self.w2 + (DT / T_TR) * (w2t - self.w2)
        if a is not None:
            e = a - self.th
            s2 = self.had.update(a - self.thr)
            if s2 is not None and t > WARMUP:
                kk = np.array([self.rules[i](s2[i], reg) for i in range(self.n)])
                self.k1l = kk[:, 0]; self.w2l = np.sqrt(kk[:, 1])
            k1, k2 = self.k1, self.w2 ** 2
            sq = np.sqrt((k1 * k1 - 4 * k2).astype(complex)); z1 = np.exp((-k1 + sq) / 2 * DTA); z2 = np.exp((-k1 - sq) / 2 * DTA)
            K1 = np.real(1 - z1 * z2); K2 = np.real((1 - z1) * (1 - z2)) / DTA
            self.th = self.th + K1 * e
            if self.armed: self.b = self.b - K2 * e
        return self.th


class Law2TV:
    """two-state Kalman, Riccati recursion online, live R. qt, qb per column."""
    def __init__(self, qt, qb, btab=None):
        self.n = n = len(qt); self.qt = np.asarray(qt, float); self.qb = np.asarray(qb, float)
        self.btab = np.array(btab if btab is not None else [False] * n, bool)
        self.th = np.zeros(n); self.b = np.zeros(n); self.thr = np.zeros(n)
        self.p11 = np.zeros(n); self.p12 = np.zeros(n); self.p22 = np.zeros(n); self.R = None
        self.had = G.HadNoise(n); self.armed = False
    def arm(self, bias_cal):
        self.armed = True; self.b = np.zeros(self.n); self.p22 = self.qt / T_PRE; self.p12 = np.zeros(self.n)
    def step(self, g, a, reg_lam, t, reg, k):
        gg = np.full(self.n, g)
        if self.armed: gg = gg - np.where(self.btab, B_RES[reg], 0.0)
        self.th = self.th + (gg - self.b) * DT; self.thr = self.thr + gg * DT
        qb = self.qb if self.armed else 0.0
        self.p11 = self.p11 - 2 * DT * self.p12 + DT * DT * self.p22 + self.qt * DT
        self.p12 = self.p12 - DT * self.p22; self.p22 = self.p22 + qb * DT
        if a is not None:
            s2 = self.had.update(a - self.thr)
            if s2 is not None and t > WARMUP:
                Rn = s2 / DTA
                if self.R is None: self.p11 = Rn.copy()
                self.R = Rn
            if self.R is not None:
                e = a - self.th; S = self.p11 + self.R; k1 = self.p11 / S; k2 = self.p12 / S; p12 = self.p12
                self.th = self.th + k1 * e; self.b = self.b + k2 * e
                self.p11 = (1 - k1) * self.p11; self.p12 = (1 - k1) * p12; self.p22 = self.p22 - k2 * p12
        return self.th


class KF(C.BenchEst):
    """the bench's two-state Kalman (F.Est kind 1) with an arming variant: 'bench' (as F.Est), 'reset' (b <- 0)."""
    def __init__(self, kR, kQ, qB, arming='bench'):
        super().__init__(np.ones(len(kR), int), kR=kR, kQ=kQ, qB=qB); self.arming = arming
    def arm(self, bias_cal):
        if self.arming == 'reset': self.e.b = np.zeros(self.n)


# ------------------------------------------------------------------------------------------------ registry
_TABS = {}
def build(names=None):
    """returns (list of estimator objects, list of column names, dict name -> description)"""
    rv = json.load(open(os.path.join(C.OUT, 'rival.json')))['rival_params']; TU = C.TUNED_A2; VT = C.FX['var_tilt']
    pub = lambda s2, reg: math.sqrt(max(s2, 1e-300) / LAM[reg])
    lamidle = lambda s2, reg: math.sqrt(max(s2, 1e-300) / LAM['idle'])
    if not _TABS:
        for key in ('raw', 'loop'):
            _TABS['a_fp_' + key] = table_a_fp(key); _TABS['a_min_' + key] = table_a_min(key); _TABS['c_' + key] = table_c(key, QB4)
    tabs = _TABS
    T = lambda nm: (lambda s2, reg, _t=tabs[nm]: float(lookup(_t, reg, s2)))
    def Tc(nm):
        t = tabs[nm]
        def f(s2, reg):
            K1, K2 = t[reg]; ls = np.log(max(s2, 1e-300)); lg = np.log(S2G)
            return float(np.exp(np.interp(ls, lg, np.log(K1)))), float(np.exp(np.interp(ls, lg, np.log(K2))))
        return f
    b = lambda qb: (lambda s2, reg: gains_closed(max(s2, 1e-300), QT, qb))
    D = {}
    # references
    kf = [('KF_rival', rv['R'], rv['QG'], rv['QB'], 'reset', 'rival: bench Kalman tuned against truth (training C3, R/Q_G/Q_B grid), reset arming'),
          ('KF_bench', TU['kf_R'], TU['kf_QG'], F.Q_B, 'bench', 'reference: the bench published tuned Kalman (R x8192, Q_G ~0, Q_B fixed)'),
          ('KF_ridge', VT * 2.0 ** 16, 0.0, F.Q_B * 2.0 ** 18, 'reset', 'reference: another point of the rival flat ridge (R x2^16, Q_B x2^18)'),
          ('KF_facts', VT, LAM['idle'], F.Q_B, 'bench', 'reference: Kalman set from facts (R = var_tilt, Q_G = lam_idle), bench arming'),
          ('KF_facts_reset', VT, LAM['idle'], F.Q_B, 'reset', 'reference: the same with the bias state reset at arming')]
    l1 = [('LAW_c1', pub, None, False, 'reference: published law kappa = 1'),
          ('a_fp_raw', T('a_fp_raw'), None, False, '(a) fixed point tau = sqrt(r(2tau)/lam(2tau)), raw-record curves'),
          ('a_min_raw', T('a_min_raw'), None, False, '(a) residual Allan crossover tau = T*/sqrt(3), raw-record curves'),
          ('a_fp_loop', T('a_fp_loop'), None, False, '(a) fixed point, looped-record (bench-world) curves'),
          ('a_min_loop', T('a_min_loop'), None, False, '(a) crossover, looped-record (bench-world) curves'),
          ('d_law_noceil', pub, 1e9, False, '(d) published law without the 30 s ceiling'),
          ('d_law_lamidle', lamidle, None, False, '(d) published law with lambda = idle ARW in every regime'),
          ('d_law_btab', pub, None, True, '(d) published law + logged per-throttle bias table')]
    l2 = [('b_ss', b(QB4), None, False, '(b) two-state law, r = live s2, q_theta = lam_idle, q_b = idle RRW bound (T_rec/4)'),
          ('b_ss_qb2', b(QB2), None, False, '(b) as b_ss with q_b at T_rec/2 (sensitivity)'),
          ('b_ss_ceil', b(QB4), TAU_MAX, False, '(b) as b_ss with the 30 s ceiling on 1/k1 (sensitivity)'),
          ('c_ss_raw', Tc('c_raw'), None, False, '(c) two-state law, r at the operating scale 2/k1, raw-record curves'),
          ('c_ss_loop', Tc('c_loop'), None, False, '(c) two-state law, r at the operating scale, looped-record curves'),
          ('d_b_ss_btab', b(QB4), None, True, '(d) b_ss + logged per-throttle bias table')]
    tv = [('b_tv', QT, QB4, False, '(b/d) two-state Kalman, Riccati online, R = live s2/DTA, Q from idle facts'),
          ('b_tv_qb2', QT, QB2, False, '(b/d) as b_tv with q_b at T_rec/2'),
          ('d_b_tv_btab', QT, QB4, True, '(d) b_tv + logged per-throttle bias table')]
    sel = lambda L: [x for x in L if names is None or x[0] in names]
    ests, cols = [], []
    k = sel(kf)
    for arming in ('reset', 'bench'):
        kk = [x for x in k if x[4] == arming]
        if kk:
            ests.append(KF(np.array([x[1] for x in kk]), np.array([x[2] for x in kk]), np.array([x[3] for x in kk]), arming)); cols += [x[0] for x in kk]
    for x in k: D[x[0]] = x[5]
    x1 = sel(l1)
    if x1:
        ests.append(Law1([x[1] for x in x1], ceiling=[x[2] for x in x1], btab=[x[3] for x in x1])); cols += [x[0] for x in x1]
        for x in x1: D[x[0]] = x[4]
    x2 = sel(l2)
    if x2:
        ests.append(Law2SS([x[1] for x in x2], ceiling=[x[2] for x in x2], btab=[x[3] for x in x2])); cols += [x[0] for x in x2]
        for x in x2: D[x[0]] = x[4]
    x3 = sel(tv)
    if x3:
        ests.append(Law2TV([x[1] for x in x3], [x[2] for x in x3], btab=[x[3] for x in x3])); cols += [x[0] for x in x3]
        for x in x3: D[x[0]] = x[4]
    if names is None or 'LAW_c2.8' in names:
        ests.append(C.BenchEst(np.array([2]), cfac=np.array([2.8]))); cols.append('LAW_c2.8'); D['LAW_c2.8'] = 'reference: published law c = 2.8 (calibrated against truth)'
    return ests, cols, D, tabs
