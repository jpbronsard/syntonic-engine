"""STEP 6, pre-registered. The attitude estimator of a quadrotor on a gimbal, in software.

Question. Step 5 showed that the rate loop of a quadrotor gives the horizon law nothing to decide; the place where a
window costs on both sides is the attitude estimator (roll/pitch fusion of gyro and accelerometer). Before any
hardware: placed there, does the law land on the oracle window, does it hold parity with a tuned Kalman filter of
the EKF type, and does it keep the optimum when the vibration level changes with throttle?

Plant and sensors (one axis, 1 kHz, gimbal: no translation). Values are representative, not measured; the gimbal
bench will measure them.
  truth      theta(t): Ornstein-Uhlenbeck, sd 5 deg, correlation time 2 s (attitude wander under the PID).
  gyro       w = dtheta/dt + b(t) + white, angle-random-walk PSD Q_G = (0.01 deg/s/sqrt(Hz))^2;
             bias b(t) = b0 + random walk of PSD Q_B = (0.005 deg/s)^2 per 100 s.
  accel tilt a_i = theta + AA80(v + e_i), i = 1, 2 (two IMUs on one board):
             v   = common-mode vibration: white of low-frequency PSD S_V = V Q_G (1 s)^2, plus a rotor line of
                   amplitude 1 deg sqrt(V) at 92 Hz V^(1/4);
             e_i = independent sensor noise, PSD (7e-5 rad/sqrt(Hz))^2 per IMU;
             AA80 = 2nd-order Butterworth anti-alias at 80 Hz (on the noise only; the truth is not delayed).
  Vibration level V per 120 s segment, three segments, 360 s per stream, seeds 600..609:
             A_steady    V = 1, 1, 1         b0 = 0
             B_throttle  V = 1, 16, 1/16     b0 = 0
             C_gyrobias  V = 1, 1, 1         b0 = 0.3 deg/s, not calibrated at start
  Low-frequency accel-tilt PSD of one IMU: S0(V) = S_V + (7e-5)^2. The optimal first-order window of this fusion
  (Kalman steady state, kappa = 1, H = 1/2) is sqrt(S0 / Q_G): 1.08 s, 4.02 s, 0.47 s.

Branchement, declared before the run. The two sensors reduce to one signal z = a - Theta, Theta the open-loop
integral of the bias-corrected gyro: its slow part is the gyro angle random walk (lambda), its fast part the accel
noise (sigma^2); the truth cancels exactly. The window tau sets the proportional gain of a Mahony filter, Kp = 1/tau;
the gyro-bias integral is tied to it by the damping of the flown firmware (Kp = 2, Ki = 0.3, zeta = 1.826):
Ki = Kp^2 / (4 zeta^2). Never the raw innovation a - theta_hat (the 22 May trap).

Candidates:
  ORACLE   Mahony, fixed tau, Ki tied, 40 windows 0.02 .. 50 s; best per segment. Yardstick of the pass criterion.
  SYN-1S   Mahony + law. Identification: engine 0.2.0 block-mean structure function of z on scales 16 ms .. 32.8 s
           (powers of two). One stream, no second read: the noise floor is A / l with A = min over l in
           {16, 32, 64, 128} ms of D(l) l. Then the engine's fit (snr gate 2, six finest scales clearing the floor)
           gives H, lambda_H, and tau = kappa_H (A/2 / lambda_H)^(1/(2H+1)). Floor dt/0.45, ceiling 30 s, rate limit
           0.4 s/s, averaging 1 s, warm-up 1 s, start 0.5 s. Engine constants unchanged except scale range and ceiling.
  SYN-2R   Mahony + the engine 0.2.0 identification UNCHANGED (scales 2 ms .. 32.8 s): read noise from the two IMUs,
           n2 from (a1 - a2)/sqrt 2, z on the pair mean. Common-mode vibration cancels in the difference.
  SYN-K    Mahony + law with the true S0 and Q_G (clairvoyant inputs): the law and the branchement, without the
           identification.
  SYN-noI  the law (SYN-1S identification) without the bias integral, z on the raw gyro: the witness.
  KF-OPT   2-state Kalman (angle, gyro bias), Q_G and Q_B true, measurement variance R swept (33 values, factor
           sqrt 2 around S0(1)/dt); the single best R per stream (an expert tuning the EKF once on the rig).
  KF-FIX   same Kalman, R = S0(1)/dt: tuned at nominal vibration and kept (EKF2's fixed noise parameters).
  KF-CLAIR same Kalman, R = S0(V)/dt switched at each segment (clairvoyant adaptive Kalman, the floor).
  MADG     one-axis reduction of Madgwick's IMU filter, theta_hat += (w - 2 beta sign(theta_hat - a)) dt, beta swept
           (best per stream); beta = 0.04 (dRehmFlight default) reported.
  CF-0.5   the flown attitude_update: fixed 0.5 s, no bias integral.
Scores: RMS angle error (deg) on the steady part of each segment (its last 80 s); window: median over the steady part.
Transients (first 40 s of a segment) are reported, not scored.

Hypotheses and decision rules (cell = scenario x segment; seed medians unless stated):
  H6.1 law, clairvoyant: SYN-K window within x1.5 of the oracle window and regret <= 10 %, all 9 cells.
  H6.2 main, the master protocol's criterion: SYN-1S within x1.5 and regret <= 10 %, all 6 cells of A and B.
  H6.3 exponent: in B, d log tau_1S / d log S0 between segments (seg 1 vs 0, seg 2 vs 0) = 0.5 +/- 0.1, seed medians;
       identified H in [0.4, 0.6] in all A and B cells.
  H6.4 parity with a tuned Kalman, steady vibration: in A, mean RMS(SYN-1S) / mean RMS(KF-OPT) <= 1.10, each segment.
  H6.5 adaptation: in B, mean RMS(SYN-1S) < mean RMS(KF-FIX) in segments 1 and 2, and mean RMS(SYN-1S) /
       mean RMS(KF-CLAIR) <= 1.10 in all three segments.
  H6.6 predicted failure (control of the branchement): SYN-2R fails the H6.2 criterion in B segments 0 and 1 (window
       too short: the inter-IMU difference does not see common-mode vibration). If it passes, that mechanism is
       wrong for this model.
  H6.7 bias: in C, SYN-1S passes the H6.2 criterion in all three cells, and mean RMS(SYN-noI) >= 2 mean RMS(SYN-1S)
       in each.
  Madgwick: no hypothesis; reported.
Not modelled, declared: accelerometer DC bias from rectified vibration (invisible to every candidate here without
aiding), lever-arm accelerations of an IMU off the gimbal pivot, large-angle nonlinearity, gyro noise rising with
throttle, temperature drift. Identification runs whenever a block completes (every 2 ticks for SYN-2R, every 16 for
SYN-1S); the window's rate limiter runs every tick. A failure is published as such.

Smoke test policy: one short run (2 seeds, 40 s segments) to check the code runs and the oracle optimum sits inside
its grid; it prints only oracle windows against sqrt(S0/Q_G), never candidate scores. Any change after it is
written here as an amendment before the full run.
"""
import hashlib, json, os, sys, time
import numpy as np
from scipy.signal import butter, lfilter
from scipy.special import gamma as G

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'results'); os.makedirs(OUT, exist_ok=True)

DT = 1e-3
D2R = np.pi / 180.0
Q_G = (0.01 * D2R) ** 2                       # rad^2/s
Q_B = (0.005 * D2R) ** 2 / 100.0              # rad^2/s^3
S_SENS = (7e-5) ** 2                          # rad^2 s, per IMU
LINE_AMP, LINE_F = 1.0 * D2R, 92.0
TRUTH_SD, TRUTH_TC = 5.0 * D2R, 2.0
ZETA = 2.0 / (2.0 * np.sqrt(0.3))             # 1.826, firmware Kp = 2, Ki = 0.3
SCEN = {'A_steady': dict(V=[1.0, 1.0, 1.0], b0=0.0),
        'B_throttle': dict(V=[1.0, 16.0, 1.0 / 16.0], b0=0.0),
        'C_gyrobias': dict(V=[1.0, 1.0, 1.0], b0=0.3 * D2R)}
# engine constants (engine.py 0.2.0)
SNR_GATE, EMA, WARMUP, TAU_RATE, FINEST = 2.0, 1.0, 1.0, 0.4, 6
FLOOR_TAU, TAU_MAX, TAU0 = DT / 0.45, 30.0, 0.5

SEG_S, SKIP_S, SEEDS = 120.0, 40.0, list(range(600, 610))
if os.environ.get('STEP6_SMOKE'):
    SEG_S, SKIP_S, SEEDS = 40.0, 15.0, [600, 601]


def S0(V):
    return V * Q_G * 1.0 ** 2 + S_SENS


def kappa_H(H):
    return (2 * H * G(2 * H + 1)) ** (-1.0 / (2 * H + 1))


def C_H(H):
    return (2.0 ** (2 * H + 2) - 4.0) / ((2 * H + 1) * (2 * H + 2))


# ---------------------------------------------------------------------------------------------------------- data
def generate(scn, seeds):
    V = np.array(SCEN[scn]['V']); b0 = SCEN[scn]['b0']
    nseg = int(round(SEG_S / DT)); N = 3 * nseg
    Vk = np.repeat(V, nseg)
    bb, ab = butter(2, 80.0 / (0.5 / DT))
    rho = np.exp(-DT / TRUTH_TC); c = TRUTH_SD * np.sqrt(1 - rho ** 2)
    TH, W, A1, A2 = (np.empty((N, len(seeds))) for _ in range(4))
    TH1 = np.empty((N + 1, len(seeds)))
    for j, s in enumerate(seeds):
        r = np.random.default_rng(s)
        th = lfilter([c], [1, -rho], r.standard_normal(N + 1)); TH1[:, j] = th
        bias = b0 + np.cumsum(np.sqrt(Q_B * DT) * r.standard_normal(N))
        W[:, j] = np.diff(th) / DT + bias + np.sqrt(Q_G / DT) * r.standard_normal(N)
        ph = 2 * np.pi * np.cumsum(LINE_F * Vk ** 0.25) * DT + r.uniform(0, 2 * np.pi)
        vib = np.sqrt(V.repeat(nseg) * Q_G / DT) * r.standard_normal(N) + LINE_AMP * np.sqrt(Vk) * np.sin(ph)
        e1 = np.sqrt(S_SENS / DT) * r.standard_normal(N); e2 = np.sqrt(S_SENS / DT) * r.standard_normal(N)
        A1[:, j] = th[1:] + lfilter(bb, ab, vib + e1)
        A2[:, j] = th[1:] + lfilter(bb, ab, vib + e2)
        TH[:, j] = th[1:]
    return dict(TH=TH, W=W, A1=A1, A2=A2, Vk=Vk, nseg=nseg, N=N, th0=TH1[0])


# ------------------------------------------------------------------------------------------------------- filters
class Mahony:
    """attitude_update() with the MAIDEN11 bias integral, one axis, vectorized over columns."""

    def __init__(self, a0, withI=True):
        self.th = a0.copy(); self.b = np.zeros_like(a0); self.Th = a0.copy(); self.withI = withI

    def step(self, w, a, tau):
        kp = 1.0 / tau
        wb = w - self.b
        self.th = self.th + wb * DT
        self.Th = self.Th + wb * DT                                  # open loop, for the identification
        e = a - self.th
        self.th = self.th + np.minimum(DT * kp, 1.0) * e
        if self.withI:
            ki = kp * kp / (4 * ZETA ** 2)
            self.b = np.clip(self.b - ki * e * DT, -0.15, 0.15)
        return self.th


class SynID:
    """engine.py 0.2.0 block-mean identification and window law. mode 'two': engine unchanged (read noise from the
    pair difference); mode 'one': one stream, noise floor read from the scales where noise dominates."""

    def __init__(self, n, mode):
        self.mode = mode
        lo = 2 if mode == 'two' else 16
        self.scales = np.array([2 ** i for i in range(16) if lo <= 2 ** i <= 32768])
        self.logl = np.log(self.scales * DT); m = len(self.scales)
        self.Aidx = np.array([i for i, s in enumerate(self.scales) if s in (16, 32, 64, 128)])
        self.acc = np.zeros((m, n)); self.cnt = np.zeros(m, dtype=int)
        self.prev = np.full((m, n), np.nan); self.DA = np.zeros((m, n)); self.nblk = np.zeros(m, dtype=int)
        self.a_blk = np.array([min(1.0, DT * s / max(EMA, 8 * DT * s)) for s in self.scales])
        self.a_w = DT / EMA; self.warm = int(round(WARMUP / DT)); self.k = 0
        self.n2 = np.zeros(n); self.tau = np.full(n, TAU0); self.tau_t = np.full(n, TAU0); self.H = np.full(n, 0.5)
        self.s2 = np.full(n, np.nan); self.lam = np.full(n, np.nan)

    def _fit(self, D, floor):
        ok_g = D > SNR_GATE * floor; ok_u = D > 0
        use_g = ok_g.sum(0) >= 2
        ok = np.where(use_g[None], ok_g, ok_u); ok &= np.cumsum(ok, axis=0) <= FINEST
        nok = ok.sum(0)
        logD = np.where(ok, np.log(np.maximum(D, 1e-300)), 0.0)
        W = ok.astype(float); L = self.logl[:, None]
        Sw = W.sum(0); Sx = (W * L).sum(0); Sy = (W * logD).sum(0)
        Sxx = (W * L ** 2).sum(0); Sxy = (W * L * logD).sum(0)
        den = Sw * Sxx - Sx ** 2
        slope = np.where(den > 1e-12, (Sw * Sxy - Sx * Sy) / np.maximum(den, 1e-12), 1.0)
        Hn = np.clip(0.5 * slope, 0.05, 1.0)
        lam = np.exp((W * (logD - 2 * Hn[None] * L)).sum(0) / np.maximum(Sw, 1)) / C_H(Hn)
        return Hn, lam, nok >= 2

    def update(self, z, d=None):
        if self.mode == 'two':
            self.n2 += self.a_w * (d * d - self.n2)
        self.acc += z[None, :]; self.cnt += 1
        done = self.cnt == self.scales; changed = False
        if done.any():
            for i in np.flatnonzero(done):
                mean = self.acc[i] / self.scales[i]
                if not np.isnan(self.prev[i, 0]):
                    dd = mean - self.prev[i]
                    aw = self.a_blk[i] if self.nblk[i] > 0 else 1.0
                    self.DA[i] += aw * (dd * dd - self.DA[i]); self.nblk[i] += 1
                self.prev[i] = mean; self.acc[i] = 0.0; self.cnt[i] = 0
            changed = True
        if self.k > self.warm:
            if changed:
                l = (self.scales * DT)[:, None]
                if self.mode == 'two':
                    floor = self.n2[None, :] * DT / l; s2 = 0.5 * self.n2 * DT; okA = np.ones(z.shape, bool)
                else:
                    good = self.nblk[self.Aidx] >= 2
                    if good.any():
                        A = np.min(self.DA[self.Aidx[good]] * l[self.Aidx[good]], axis=0); okA = np.ones(z.shape, bool)
                    else:
                        A = np.full(z.shape, np.nan); okA = np.zeros(z.shape, bool)
                    floor = A[None, :] / l; s2 = 0.5 * A
                D = np.where((self.nblk >= 2)[:, None], self.DA - floor, -1.0)
                Hn, lam, valid = self._fit(D, floor); valid &= okA
                ts = kappa_H(Hn) * (np.maximum(s2, 1e-300) / np.maximum(lam, 1e-300)) ** (1.0 / (2 * Hn + 1))
                ts = np.clip(np.nan_to_num(ts, nan=TAU0), FLOOR_TAU, TAU_MAX)
                self.tau_t = np.where(valid, ts, self.tau_t); self.H = np.where(valid, Hn, self.H)
                self.s2 = np.where(valid, s2, self.s2); self.lam = np.where(valid, lam, self.lam)
            rr = 1.0 / self.tau; dr = TAU_RATE * DT / self.tau ** 2
            rr = rr + np.clip(1.0 / self.tau_t - rr, -dr, dr); self.tau = 1.0 / rr
        self.k += 1
        return self.tau


class Kalman2:
    """angle + gyro bias, discrete, vectorized; R per column (may change between calls)."""

    def __init__(self, a0, R):
        n = a0.shape[0]
        self.x1 = a0.copy(); self.x2 = np.zeros(n)
        self.p11 = np.asarray(R, float).copy() * np.ones(n); self.p12 = np.zeros(n); self.p22 = np.full(n, (0.5 * D2R) ** 2)

    def step(self, w, a, R):
        self.x1 = self.x1 + (w - self.x2) * DT
        p11 = self.p11 - 2 * DT * self.p12 + DT * DT * self.p22 + Q_G * DT
        p12 = self.p12 - DT * self.p22; p22 = self.p22 + Q_B * DT
        y = a - self.x1; S = p11 + R; k1 = p11 / S; k2 = p12 / S
        self.x1 = self.x1 + k1 * y; self.x2 = self.x2 + k2 * y
        self.p11 = (1 - k1) * p11; self.p12 = (1 - k1) * p12; self.p22 = p22 - k2 * p12
        return self.x1


class Madg1:
    def __init__(self, a0):
        self.th = a0.copy()

    def step(self, w, a, two_beta):
        self.th = self.th + (w - two_beta * np.sign(self.th - a)) * DT
        return self.th


# ----------------------------------------------------------------------------------------------------- scoring
class Score:
    def __init__(self, ncol, nseg, skip):
        self.nseg, self.skip = nseg, skip
        self.ss = np.zeros((3, ncol)); self.tr = np.zeros((3, ncol)); self.ns = np.zeros(3); self.nt = np.zeros(3)

    def add(self, k, err):
        sg = k // self.nseg
        if (k % self.nseg) >= self.skip:
            self.ss[sg] += err * err; self.ns[sg] += 1
        else:
            self.tr[sg] += err * err; self.nt[sg] += 1

    def rms(self):
        return np.sqrt(self.ss / self.ns[:, None]) / D2R, np.sqrt(self.tr / self.nt[:, None]) / D2R


TAU_GRID = np.logspace(np.log10(0.02), np.log10(50.0), 40)
R_MULT = 2.0 ** np.arange(-8, 8.5, 0.5)                                     # 33 values
BETA2_GRID = np.concatenate([np.logspace(-4, 0, 30), [0.08]])               # 2 beta; last = dRehmFlight 0.04


def run(scn):
    d = generate(scn, SEEDS); ns = len(SEEDS); N, nseg = d['N'], d['nseg']; skip = int(round(SKIP_S / DT))
    TH, W, A1, A2, Vk = d['TH'], d['W'], d['A1'], d['A2'], d['Vk']
    Am = 0.5 * (A1 + A2); Dd = (A1 - A2) / np.sqrt(2.0)
    nt, nr, nb = len(TAU_GRID), len(R_MULT), len(BETA2_GRID)
    rep = lambda x, m: np.tile(x, m)                                        # column layout: [grid value][seed]
    a0 = A1[0]
    orc = Mahony(rep(a0, nt)); tau_or = np.repeat(TAU_GRID, ns)
    cf05 = Mahony(a0, withI=False)
    s1 = Mahony(a0); id1 = SynID(ns, 'one')
    s2r = Mahony(Am[0]); id2 = SynID(ns, 'two')
    sk = Mahony(a0)
    snoi = Mahony(a0, withI=False); idn = SynID(ns, 'one')
    R_nom = S0(1.0) / DT
    kopt = Kalman2(rep(a0, nr), np.repeat(R_MULT * R_nom, ns)); Ropt = np.repeat(R_MULT * R_nom, ns)
    kfix = Kalman2(a0, R_nom); kcl = Kalman2(a0, S0(Vk[0]) / DT)
    mad = Madg1(rep(a0, nb)); b2 = np.repeat(BETA2_GRID, ns)
    sc = {k: Score(c, nseg, skip) for k, c in [('orc', nt * ns), ('cf05', ns), ('s1', ns), ('s2r', ns), ('sk', ns),
                                               ('snoi', ns), ('kopt', nr * ns), ('kfix', ns), ('kcl', ns), ('mad', nb * ns)]}
    taulog = {k: [] for k in ('s1', 's2r', 'snoi')}; Hlog = {k: [] for k in ('s1', 's2r', 'snoi')}; klog = []
    t0 = time.time()
    for k in range(N):
        w, a, th, V = W[k], A1[k], TH[k], Vk[k]
        tk = np.full(ns, np.sqrt(S0(V) / Q_G))
        sc['orc'].add(k, orc.step(rep(w, nt), rep(a, nt), tau_or) - rep(th, nt))
        sc['cf05'].add(k, cf05.step(w, a, 0.5) - th)
        t1 = id1.update(a - s1.Th); sc['s1'].add(k, s1.step(w, a, t1) - th)
        t2 = id2.update(Am[k] - s2r.Th, Dd[k]); sc['s2r'].add(k, s2r.step(w, Am[k], t2) - th)
        sc['sk'].add(k, sk.step(w, a, tk) - th)
        tn = idn.update(a - snoi.Th); sc['snoi'].add(k, snoi.step(w, a, tn) - th)
        sc['kopt'].add(k, kopt.step(rep(w, nr), rep(a, nr), Ropt) - rep(th, nr))
        sc['kfix'].add(k, kfix.step(w, a, R_nom) - th)
        sc['kcl'].add(k, kcl.step(w, a, S0(V) / DT) - th)
        sc['mad'].add(k, mad.step(rep(w, nb), rep(a, nb), b2) - rep(th, nb))
        if k % 100 == 0:
            for key, idr in (('s1', id1), ('s2r', id2), ('snoi', idn)):
                taulog[key].append(idr.tau.copy()); Hlog[key].append(idr.H.copy())
            klog.append(k)
        if k % 60000 == 0 and k:
            print(f'  {scn} t = {k * DT:.0f} s  ({time.time() - t0:.0f} s)', flush=True)
    klog = np.array(klog)
    res = {'scenario': scn, 'V': SCEN[scn]['V'], 'tau_star_true': [float(np.sqrt(S0(v) / Q_G)) for v in SCEN[scn]['V']]}
    for key, s in sc.items():
        st, tr = s.rms(); res[key] = {'steady': st.tolist(), 'trans': tr.tolist()}
    for key in taulog:
        T = np.array(taulog[key]); Hh = np.array(Hlog[key])
        med_t, med_H = [], []
        for g in range(3):
            m = (klog // nseg == g) & ((klog % nseg) >= skip)
            med_t.append(np.median(T[m], axis=0).tolist()); med_H.append(np.median(Hh[m], axis=0).tolist())
        res[key]['tau_med'] = med_t; res[key]['H_med'] = med_H
    res['grids'] = {'tau': TAU_GRID.tolist(), 'R_mult': R_MULT.tolist(), 'beta2': BETA2_GRID.tolist(), 'R_nom': R_nom}
    res['seeds'] = SEEDS; res['runtime_s'] = time.time() - t0
    return res


def oracle_window(res, ns):
    """per segment, per seed: best fixed window (parabolic refinement in log tau) and its RMS."""
    st = np.array(res['orc']['steady']).reshape(3, len(TAU_GRID), ns)
    out_t, out_r = np.empty((3, ns)), np.empty((3, ns))
    lg = np.log(TAU_GRID)
    for g in range(3):
        for j in range(ns):
            y = st[g, :, j]; i = int(np.argmin(y)); out_r[g, j] = y[i]
            if 0 < i < len(y) - 1:
                x0, x1, x2 = lg[i - 1:i + 2]; y0, y1, y2 = y[i - 1:i + 2]
                den = (y0 - 2 * y1 + y2)
                off = 0.5 * (y0 - y2) / den if den > 0 else 0.0
                out_t[g, j] = np.exp(x1 + off * (x1 - x0))
            else:
                out_t[g, j] = TAU_GRID[i]
    return out_t, out_r


def main():
    sha = hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest()
    print('step6_quad_attitude_estimator.py  sha256', sha, flush=True)
    smoke = bool(os.environ.get('STEP6_SMOKE'))
    allres = {'sha256': sha, 'smoke': smoke}
    for scn in SCEN:
        res = run(scn); ns = len(SEEDS)
        tor, ror = oracle_window(res, ns)
        res['oracle'] = {'tau': tor.tolist(), 'rms': ror.tolist()}
        allres[scn] = res
        print(f'{scn}: oracle window per segment (seed median) {np.round(np.median(tor, 1), 3).tolist()}'
              f'  vs sqrt(S0/Q_G) {np.round(res["tau_star_true"], 3).tolist()}  runtime {res["runtime_s"]:.0f} s', flush=True)
    fn = os.path.join(OUT, 'step6_smoke.json' if smoke else 'step6_rows.json')
    json.dump(allres, open(fn, 'w'))
    print('written', fn)


if __name__ == '__main__':
    main()
