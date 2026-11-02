# Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
# Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
"""A syntonic flight controller (one axis, roll) against a PID and against a PID with a tuned Kalman filter, in the same
conditions, on the X500's own logged vibration noise. Written 2 Nov 2026 before any run of this file except the smoke test;
request and go from JP ("monte moi le meilleur controleur de vol syntonique avec toutes nos connaissances, k = 1, compare
le a un PID dans les memes conditions et a un Kalman ajuste"; scope chosen: one axis, real X500 noise, three full
controllers; then "peut-etre k dynamique ici").

PLANT AND CONDITIONS (identical for every controller)
  One axis: motor state u follows the command u_c with time constant TAU_M = 30 ms; omega_dot = B u - DAMP omega + d;
  theta_dot = omega. B = 4000 deg/s^2 per unit, DAMP = 2 /s, |u_c| <= 0.25 (step-5 plant, X500-like). Loop 1 kHz.
  Disturbance d: Ornstein-Uhlenbeck torque (std 60 deg/s^2, correlation 1 s) plus 4 gusts per flight (+-300 deg/s^2 for
  0.25 s, random times). Pre-roll 5 s on the ground (theta clamped at 0, idle noise): every estimator starts there; the
  gyro bias is calibrated at the end of the pre-roll by the pre-roll mean (arming calibration, same for all).
  Flight 60 s, throttle regimes for the noise: half 0-15 s, low 15-30, half 30-45, low 45-60 s.
  Setpoint: hover (0 deg) 0-10 s and 50-60 s; manoeuvres 10-50 s, square +-10 deg switching every 2.5 s, shaped by a
  critically damped second-order prefilter (omega 10 rad/s); the shaped setpoint is the reference for every controller.
  Sensors: gyro = omega + real gyro noise (1 kHz, the four sub-samples of the logger); accelerometer tilt = theta + real
  tilt noise (250 Hz, held 4 ticks). Real noise, mean removed, from bench/software_2026-10-26/data: idle = S6-idle,
  low = S6-low (40 s), half = S5-half (8 s, the only half-throttle record: looped, declared). Each seed draws random
  circular offsets into the records, its own disturbance and gust times. Linear tilt model (as step 6g).

CONTROLLERS (same cascade: angle P -> rate PID -> motor; gyro low-passed at the engine's floor, DT_ACC/0.45 ... see code)
  PID   attitude by a complementary (Mahony) filter with the FIXED window 0.5 s of the flown firmware (ATT_TAU_S) and its
        bias rule Ki = Kp^2 / (4 zeta^2), zeta 1.826; no feed-forward. Gains tuned (grid below).
  KF    attitude by a two-state Kalman filter (angle, gyro bias), R and Q_G tuned once on training seeds (estimation
        error against truth; R in var(tilt) x 2^(-8..8), Q_G in lambda_hat x 2^(-2..8)); no feed-forward. Gains tuned.
  SYN1  the syntonic controller with what the benches taught: window tau = kappa sqrt(sigma^2 / lambda), KAPPA = 1
        (canonical class: gyro angle random walk against white accelerometer noise; at this class the steady-state
        Kalman window is the same square root); sigma^2 read live (block-mean structure function of accel - open-loop
        gyro integral at 32, 64, 128 ms, engine EMA 1 s, step 6g); lambda = the gyro's drift fact at the operating
        throttle (Allan variance at 1, 2, 4 s, floored at the idle fact; step 6h), measured once per throttle level;
        start-up schedule tau <= max(0.5 s, t) (6h); bias loop off, arming calibration (6h exploratory); the window
        travels in 1/tau (e-geodesic) with T = 128 ms, the longest identification scale (1-2 Nov transport benches);
        floor DT_ACC/0.45, ceiling 30 s; rate loop at the floor (step 5); setpoint feed-forward (rate and torque from
        the shaped setpoint) with the motor lead TAU_M on the feed-forward torque (MAESTRO lead, 1 Nov). Gains tuned
        like the PID's: the law decides when, not how much.
  SYND  SYN1 with an in-flight dynamic kappa: a bank of five windows c x tau, c in {0.5, 0.7, 1, 1.4, 2}, each its own
        filter; cost = EMA (10 s) of the squared accelerometer innovation (accel - prior estimate); the output is the
        filter with the lowest cost. The only in-flight signal; declared constants.
  SYNB  SYN1 with kappa calibrated once against truth on the training seeds (the gimbal bench analogue: encoders as
        truth), c in {0.25, 0.35, 0.5, 0.7, 1, 1.4, 2}, then fixed in flight (the corpus's kappa* = kappa0 sqrt(theta0)
        for a biased measurement pipeline).
  Gain grid, each arm tuned on the training seeds 100-102 (mean attitude RMS): Ka in {2, 4, 8, 16} /s; rate Kp in
  {10, 20, 40, 80} / B; rate Kd in {0, 2e-4, 5e-4}; rate Ki in {0, 2 Kp}. SYND and SYNB use SYN1's gains.
Diagnostics, described: PIDFF (PID estimator and gains plus SYN's feed-forward), SYN1noFF (SYN1 without feed-forward).

SCORES (test seeds 200-209, never used in tuning): J_att = RMS(theta - setpoint) over the 60 s flight (deg); J_est =
RMS(estimate - theta) (deg); also by regime and phase; motor jitter RMS(diff u_c). |theta| > 90 deg = unstable (inf).

PREDICTIONS (judged; pooled = rms over the 10 test seeds of the per-seed RMS):
F1  SYN1's estimator beats the flown fixed window: pooled J_est(SYN1) <= 0.5 x pooled J_est(PID).
F2  SYN1's estimator matches or beats the tuned Kalman: pooled J_est(SYN1) <= pooled J_est(KF).
F3  SYN1 flies better than the PID: pooled J_att(SYN1) <= 0.8 x pooled J_att(PID), below in >= 8 of 10 seeds.
F4  SYN1 flies at least as well as PID + tuned Kalman: pooled J_att(SYN1) <= pooled J_att(KF), below in >= 7 of 10.
K1  (assistant) the in-flight dynamic kappa hears nothing: pooled J_est(SYND) >= 0.95 x pooled J_est(SYN1).
K2  (assistant) a pipeline bias is real here (6h: windows 1.7 to 2.1x too long): the calibrated factor c_bench <= 0.7,
    and pooled J_est(SYNB) <= 0.9 x pooled J_est(SYN1).
JPq (author, "kappa tends to 1") the calibrated factor c_bench lies in [0.8, 1.25].

Run: python flight_bench.py <repo_root> smoke | facts | calibrate | tune <arm> | run | summarize
"""
import os, sys, json, time, hashlib, math
import numpy as np

REPO = os.path.abspath(sys.argv[1]); MODE = sys.argv[2]
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, 'results'); os.makedirs(OUT, exist_ok=True)
sys.path.insert(0, os.path.join(REPO, 'bench', 'software_2026-10-26'))
import step6g_real_noise as G
SHA = hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest()

D2R = math.pi / 180
DT, DTA, NSUB = 1e-3, 4e-3, 4
T_PRE, T_FLY = 5.0, 60.0
REGS = [('half', 0, 15), ('low', 15, 30), ('half', 30, 45), ('low', 45, 60)]
B, DAMP, UMAX, TAU_M = 4000.0, 2.0, 0.25, 0.030
OU_STD, OU_TAU, N_GUST, GUST_A, GUST_T = 60.0, 1.0, 4, 300.0, 0.25
SP_AMP, SP_HALF, SP_W = 10.0, 2.5, 10.0
ZETA = 1.826; TAU_CF = 0.5
FLOOR_A = DTA / 0.45; TAU_MAX = 30.0; T_TR = 0.128
TAU_G = DT / 0.45                                     # gyro low-pass at the 1 kHz floor, dt/0.45
TAU_D = 0.008
BANK = np.array([0.5, 0.7, 1.0, 1.4, 2.0]); BANK_EMA = 10.0
CAL = np.array([0.25, 0.35, 0.5, 0.7, 1.0, 1.4, 2.0])
TRAIN, TEST = [100, 101, 102], list(range(200, 210))
KAS, KPS, KDS, KIM = [2.0, 4.0, 8.0, 16.0], [10.0 / B, 20.0 / B, 40.0 / B, 80.0 / B], [0.0, 2e-4, 5e-4], [0.0, 2.0]
GRID = [(ka, kp, kd, km) for ka in KAS for kp in KPS for kd in KDS for km in KIM]
Q_B = (0.005 * D2R) ** 2 / 100.0

# ----------------------------------------------------------------------------------------------------- noise records
_REC = {}
def records():
    if not _REC:
        for name, seg in (('idle', G.SEGS[0]), ('low', G.SEGS[1]), ('half', G.SEGS[2])):
            tilt, _, g1k = G.load(*seg[:4]); _REC[name] = (tilt - tilt.mean(), g1k - g1k.mean())
    return _REC

def facts():
    R = records(); _, _, gi = G.load(*G.SEGS[0][:4]); lam_hat, _ = G.allan_lambda(gi)
    lam = {'idle': lam_hat}
    for name in ('low', 'half'):
        g = R[name][1]
        def av(l):
            k = int(round(l / 1e-3)); mb = g[: len(g) // k * k].reshape(-1, k).mean(1); return 0.5 * np.mean(np.diff(mb) ** 2)
        ls = [l for l in (1.0, 2.0, 4.0) if len(g) > 6 * int(l / 1e-3)]
        lam[name] = max(float(np.mean([av(l) * l for l in ls])), lam_hat)
    var_tilt = float(np.mean(np.concatenate([R[n][0] for n in ('low', 'half')]) ** 2))
    return dict(lam=lam, lam_hat=lam_hat, var_tilt=var_tilt)

def build_noise(seed):
    """gyro noise per 1 ms tick and tilt noise per 4 ms sample over pre-roll + flight; regime label per tick."""
    R = records(); rng = np.random.default_rng(seed)
    segs = [('idle', 0.0, T_PRE)] + [(n, T_PRE + a, T_PRE + b) for n, a, b in REGS]
    gn, tn, reg = [], [], []
    for name, a, b in segs:
        tilt, g = R[name]; ng = int(round((b - a) / DT)); na = ng // NSUB
        o = int(rng.integers(0, len(tilt))); idx_a = (o + np.arange(na)) % len(tilt)
        idx_g = (o * NSUB + np.arange(ng)) % len(g)
        gn.append(g[idx_g]); tn.append(tilt[idx_a]); reg += [name] * ng
    return np.concatenate(gn), np.concatenate(tn), np.array(reg)

def build_world(seed):
    rng = np.random.default_rng(10_000 + seed)
    n = int(round((T_PRE + T_FLY) / DT)); t = np.arange(n) * DT - T_PRE
    a = math.exp(-DT / OU_TAU); s = OU_STD * math.sqrt(1 - a * a); w = rng.standard_normal(n); d = np.empty(n); x = 0.0
    for k in range(n): x = a * x + s * w[k]; d[k] = x
    for tg in rng.uniform(2.0, T_FLY - 2.0, N_GUST):
        d[(t >= tg) & (t < tg + GUST_T)] += GUST_A * rng.choice([-1.0, 1.0])
    raw = np.where((t >= 10) & (t < 50), np.where(((t - 10) // SP_HALF) % 2 == 0, SP_AMP, -SP_AMP), 0.0)
    sp = np.zeros(n); spd = np.zeros(n); spdd = np.zeros(n); x = v = 0.0
    for k in range(n):
        acc = SP_W ** 2 * (raw[k] - x) - 2 * SP_W * v; v += acc * DT; x += v * DT
        sp[k], spd[k], spdd[k] = x, v, acc
    return t, d, sp, spd, spdd

# ----------------------------------------------------------------------------------------------------- estimators
class Est:
    """Vectorised attitude estimators. kind per column: 0 Mahony fixed tau (+ flown Ki rule), 1 Kalman 2-state, 2 law
    (bank_mask: law columns whose output comes from the five-window bank)."""
    def __init__(self, n, kind, kf_R=None, kf_QG=None, cfac=None, bank_mask=None):
        self.n = n; self.kind = np.asarray(kind)
        self.th = np.zeros(n); self.b = np.zeros(n); self.thr = np.zeros(n)
        self.kR = np.ones(n) if kf_R is None else np.asarray(kf_R, float); self.kQ = np.ones(n) if kf_QG is None else np.asarray(kf_QG, float)
        self.p11 = self.kR.copy(); self.p12 = np.zeros(n); self.p22 = np.full(n, (0.5 * D2R) ** 2)
        self.cfac = np.ones(n) if cfac is None else np.asarray(cfac, float)
        self.had = G.HadNoise(n); self.tau = np.full(n, 0.5); self.taul = np.full(n, 0.5)
        self.bmask = np.zeros(n, bool) if bank_mask is None else np.asarray(bank_mask, bool)
        self.bank = bool(self.bmask.any())
        if self.bank:
            self.bth = np.zeros((len(BANK), n)); self.bcost = np.zeros((len(BANK), n)); self.bsel = np.full(n, 2)
            self.btau = np.full((len(BANK), n), 0.5)

    def step(self, g, a, reg_lam, t_since):
        """g: calibrated gyro (rad/s) per column; a: tilt (rad) per column, or None between accelerometer samples."""
        mah = self.kind == 0; kal = self.kind == 1; law = self.kind == 2
        wb = g - self.b
        self.th = self.th + wb * DT; self.thr = self.thr + wb * DT
        if self.bank: self.bth = self.bth + wb[None, :] * DT
        p11 = self.p11 - 2 * DT * self.p12 + DT * DT * self.p22 + self.kQ * DT
        p12 = self.p12 - DT * self.p22; p22 = self.p22 + Q_B * DT
        self.p11 = np.where(kal, p11, self.p11); self.p12 = np.where(kal, p12, self.p12); self.p22 = np.where(kal, p22, self.p22)
        if law.any():                                   # the window travels in 1/tau toward its scheduled target
            cap = max(0.5, t_since)
            tgt = np.clip(np.minimum(self.cfac * self.taul, cap), FLOOR_A, TAU_MAX)
            rr = 1.0 / self.tau; rr = rr + (DT / T_TR) * (1.0 / tgt - rr); self.tau = np.where(law, 1.0 / rr, self.tau)
            if self.bank:
                btgt = np.clip(np.minimum(BANK[:, None] * self.taul[None, :], cap), FLOOR_A, TAU_MAX)
                brr = 1.0 / self.btau; brr = brr + (DT / T_TR) * (1.0 / btgt - brr); self.btau = 1.0 / brr
        if a is not None:
            e = a - self.th
            kp = 1.0 / TAU_CF; ki = kp * kp / (4 * ZETA ** 2)
            th_m = self.th + min(DTA * kp, 1.0) * e; b_m = np.clip(self.b - ki * e * DTA, -0.15, 0.15)
            S = self.p11 + self.kR; k1 = self.p11 / S; k2 = self.p12 / S
            th_k = self.th + k1 * e; b_k = self.b + k2 * e
            p12_pred = self.p12
            self.p11 = np.where(kal, (1 - k1) * self.p11, self.p11); self.p12 = np.where(kal, (1 - k1) * p12_pred, self.p12)
            self.p22 = np.where(kal, self.p22 - k2 * p12_pred, self.p22)
            th_l = self.th
            if law.any():
                s2 = self.had.update(a - self.thr)
                if s2 is not None and t_since > G.WARMUP:
                    self.taul = np.where(law, np.sqrt(np.maximum(s2, 1e-300) / reg_lam), self.taul)
                th_l = self.th + np.minimum(DTA / self.tau, 1.0) * e
                if self.bank:
                    eb = a[None, :] - self.bth
                    self.bcost = self.bcost + (DTA / BANK_EMA) * (eb * eb - self.bcost)
                    self.bth = self.bth + np.minimum(DTA / self.btau, 1.0) * eb
                    self.bsel = np.argmin(self.bcost, axis=0)
                    th_l = np.where(self.bmask, self.bth[self.bsel, np.arange(self.n)], th_l)
            self.th = np.where(mah, th_m, np.where(kal, th_k, th_l))
            self.b = np.where(mah, b_m, np.where(kal, b_k, self.b))
        return self.th

# ----------------------------------------------------------------------------------------------------- closed loop
def arm_spec(name, fx, gains, tuned):
    """returns (estimator kind, cfac, bank, feed-forward on, gains)"""
    if name == 'PID': return 0, 1.0, False, False, gains
    if name == 'PIDFF': return 0, 1.0, False, True, gains
    if name == 'KF': return 1, 1.0, False, False, gains
    if name == 'SYN1': return 2, 1.0, False, True, gains
    if name == 'SYN1noFF': return 2, 1.0, False, False, gains
    if name == 'SYND': return 2, 1.0, True, True, gains
    if name == 'SYNB': return 2, tuned.get('c_bench', 1.0), False, True, gains
    raise ValueError(name)

def fly(seed, cols, fx, tuned, plant=True):
    """cols: list of (arm name, gains tuple). Returns per-column scores."""
    gn, tn, reg = build_noise(seed); t, d, sp, spd, spdd = build_world(seed)
    n = len(cols); N = len(t); kpre = int(round(T_PRE / DT))
    specs = [arm_spec(nm, fx, gn_, tuned) for nm, gn_ in cols]
    kind = np.array([s[0] for s in specs]); cfac = np.array([s[1] for s in specs]); ff = np.array([s[3] for s in specs])
    bank_cols = [i for i, s in enumerate(specs) if s[2]]
    Ka = np.array([s[4][0] for s in specs]); Kp = np.array([s[4][1] for s in specs]); Kd = np.array([s[4][2] for s in specs])
    Ki = np.array([s[4][1] * s[4][3] for s in specs])
    kR = np.full(n, tuned.get('kf_R', 1.0)); kQ = np.full(n, tuned.get('kf_QG', 1.0))
    bmask = np.zeros(n, bool); bmask[bank_cols] = True
    est = Est(n, kind, kf_R=kR, kf_QG=kQ, cfac=cfac, bank_mask=bmask)
    lam_of = {r: fx['lam'][r] for r in ('idle', 'low', 'half')}
    th = np.zeros(n); om = np.zeros(n); u = np.zeros(n); wf = np.zeros(n); dwf = np.zeros(n); I = np.zeros(n); uc_prev = np.zeros(n)
    gsum = 0.0; gcnt = 0; bias_cal = 0.0
    e_att = np.zeros(n); e_est = np.zeros(n); jit = np.zeros(n); m = 0; blown = np.zeros(n, bool)
    regs_att = {r: np.zeros(n) for r in ('low', 'half')}; regs_m = {r: 0 for r in ('low', 'half')}
    a_ = math.exp(-DT / TAU_G); a_d = math.exp(-DT / TAU_D)
    for k in range(N):
        flying = k >= kpre
        graw = om * D2R + gn[k]; g = graw - bias_cal
        a = (th * D2R + tn[k // NSUB]) if k % NSUB == 0 else None
        thh = est.step(g, a, lam_of[reg[k]], t[k] + T_PRE)
        if not flying:
            gsum += gn[k]; gcnt += 1
            if k == kpre - 1: bias_cal = gsum / gcnt          # arming calibration: pre-roll mean of the gyro
            continue
        # rate measurement: gyro minus the estimator's bias, low-passed at the floor
        wm = (g - est.b) / D2R
        wf_prev = wf.copy(); wf = a_ * wf + (1 - a_) * wm
        dwf = a_d * dwf + (1 - a_d) * (wf - wf_prev) / DT
        ffr = np.where(ff, spd[k], 0.0)
        wsp = Ka * (sp[k] - thh / D2R) + ffr
        e = wsp - wf
        uff = np.where(ff, (spdd[k] + DAMP * spd[k]) / B, 0.0)
        uff_dot = np.where(ff, ((spdd[k] + DAMP * spd[k]) - (spdd[k - 1] + DAMP * spd[k - 1])) / (B * DT), 0.0)
        ud = Kp * e + Ki * I - Kd * dwf + uff + TAU_M * uff_dot
        uc = np.clip(ud, -UMAX, UMAX)
        I = np.where(np.abs(ud) < UMAX, I + e * DT, I)                  # anti-windup
        if plant:
            u = u + (DT / TAU_M) * (uc - u)
            om = om + DT * (B * u - DAMP * om + d[k]); th = th + DT * om
        blown |= ~np.isfinite(th) | (np.abs(th) > 90)
        th = np.where(blown, 0.0, th); om = np.where(blown, 0.0, om)
        err = th - sp[k]; e_att += err ** 2; e_est += (thh / D2R - th) ** 2; jit += (uc - uc_prev) ** 2; uc_prev = uc; m += 1
        r = reg[k]; regs_att[r] += err ** 2; regs_m[r] += 1
    J_att = np.where(blown, np.inf, np.sqrt(e_att / m)); J_est = np.sqrt(e_est / m)
    return dict(J_att=J_att.tolist(), J_est=J_est.tolist(), jit=np.sqrt(jit / m).tolist(),
                att_by_reg={r: np.where(blown, np.inf, np.sqrt(v / max(regs_m[r], 1))).tolist() for r, v in regs_att.items()})

def est_only(seed, fx, kind, cfac=None, kR=None, kQ=None):
    """estimation error with the plant at rest (the filters are linear: the error does not depend on the motion)."""
    n = len(kind); cols = [('EST', None)] * n
    gn, tn, reg = build_noise(seed); N = len(gn); kpre = int(round(T_PRE / DT))
    est = Est(n, np.asarray(kind), kf_R=kR, kf_QG=kQ, cfac=cfac)
    lam_of = {r: fx['lam'][r] for r in ('idle', 'low', 'half')}; gsum = 0.0; gcnt = 0; ss = np.zeros(n); m = 0; bias_cal = 0.0
    for k in range(N):
        g = np.full(n, gn[k] - bias_cal); a = np.full(n, tn[k // NSUB]) if k % NSUB == 0 else None
        thh = est.step(g, a, lam_of[reg[k]], k * DT)
        if k < kpre:
            gsum += gn[k]; gcnt += 1
            if k == kpre - 1: bias_cal = gsum / gcnt
            continue
        ss += (thh / D2R) ** 2; m += 1
    return np.sqrt(ss / m)

# ----------------------------------------------------------------------------------------------------- modes
def load_json(name, default=None):
    p = os.path.join(OUT, name)
    return json.load(open(p)) if os.path.exists(p) else default

if __name__ == '__main__':
    if MODE == 'smoke':
        fx = facts(); tuned = dict(kf_R=fx['var_tilt'], kf_QG=fx['lam']['half'], c_bench=1.0)
        t0 = time.time(); g = (8.0, 40.0 / B, 2e-4, 2.0)
        r = fly(100, [('PID', g), ('KF', g), ('SYN1', g), ('SYND', g), ('SYNB', g)], fx, tuned)
        print('smoke ok', {k: np.round(v, 3).tolist() if isinstance(v, list) else v for k, v in r.items() if k != 'att_by_reg'}, round(time.time() - t0, 1), 's')
        sys.exit(0)
    if MODE == 'facts':
        fx = facts(); fx['script_sha256'] = SHA; json.dump(fx, open(os.path.join(OUT, 'facts.json'), 'w'), indent=1); print(json.dumps(fx)); sys.exit(0)
    if MODE == 'calibrate':
        fx = load_json('facts.json')
        # Kalman: R x Q_G grid, estimation error on the training seeds
        Rg = fx['var_tilt'] * 2.0 ** np.arange(-8, 9); Qg = fx['lam_hat'] * 2.0 ** np.arange(-2, 9)
        kR = np.repeat(Rg, len(Qg)); kQ = np.tile(Qg, len(Rg)); nk = len(kR)
        Jk = np.mean([est_only(s, fx, np.ones(nk, int), kR=kR, kQ=kQ) ** 2 for s in TRAIN], 0) ** 0.5
        i = int(np.argmin(Jk))
        # bench calibration of the law's factor against truth
        Jc = np.mean([est_only(s, fx, np.full(len(CAL), 2), cfac=CAL) ** 2 for s in TRAIN], 0) ** 0.5
        j = int(np.argmin(Jc))
        tuned = dict(kf_R=float(kR[i]), kf_QG=float(kQ[i]), kf_R_mult=float(kR[i] / fx['var_tilt']), kf_QG_mult=float(kQ[i] / fx['lam_hat']),
                     kf_J_est=float(Jk[i]), c_bench=float(CAL[j]), cal_curve={str(c): float(v) for c, v in zip(CAL, Jc)}, script_sha256=SHA)
        json.dump(tuned, open(os.path.join(OUT, 'tuned_estimators.json'), 'w'), indent=1); print(json.dumps(tuned)); sys.exit(0)
    if MODE == 'tune':
        arm = sys.argv[3]; fx = load_json('facts.json'); tuned = load_json('tuned_estimators.json')
        cols = [(arm, g) for g in GRID]
        J = np.mean([np.square(fly(s, cols, fx, tuned)['J_att']) for s in TRAIN], 0) ** 0.5
        i = int(np.argmin(J)); res = dict(arm=arm, gains=list(GRID[i]), J_train=float(J[i]), script_sha256=SHA,
                                          top5=[(list(GRID[j]), float(J[j])) for j in np.argsort(J)[:5]])
        json.dump(res, open(os.path.join(OUT, f'gains_{arm}.json'), 'w'), indent=1); print(json.dumps(res)); sys.exit(0)
    if MODE == 'run':
        fx = load_json('facts.json'); tuned = load_json('tuned_estimators.json')
        gP = tuple(load_json('gains_PID.json')['gains']); gK = tuple(load_json('gains_KF.json')['gains']); gS = tuple(load_json('gains_SYN1.json')['gains'])
        cols = [('PID', gP), ('KF', gK), ('SYN1', gS), ('SYND', gS), ('SYNB', gS), ('PIDFF', gP), ('SYN1noFF', gS)]
        seeds = [int(s) for s in sys.argv[3:]] if len(sys.argv) > 3 else TEST
        for s in seeds:
            t0 = time.time(); r = fly(s, cols, fx, tuned)
            row = dict(seed=s, sha=SHA, arms=[c[0] for c in cols], **r)
            with open(os.path.join(OUT, 'test_runs.jsonl'), 'a') as f: f.write(json.dumps(row) + '\n')
            print(s, {c[0]: (round(r['J_att'][i], 3), round(r['J_est'][i], 3)) for i, c in enumerate(cols)}, round(time.time() - t0), 's', flush=True)
        sys.exit(0)
    if MODE == 'summarize':
        rows = [json.loads(l) for l in open(os.path.join(OUT, 'test_runs.jsonl'))]; assert len(rows) == 10, len(rows)
        arms = rows[0]['arms']; tuned = load_json('tuned_estimators.json')
        A = {a: np.array([r['J_att'][i] for r in rows]) for i, a in enumerate(arms)}
        E = {a: np.array([r['J_est'][i] for r in rows]) for i, a in enumerate(arms)}
        P = lambda v: float(np.sqrt(np.mean(np.square(v))))
        PA = {a: P(v) for a, v in A.items()}; PE = {a: P(v) for a, v in E.items()}
        V = dict(F1=PE['SYN1'] <= 0.5 * PE['PID'], F2=PE['SYN1'] <= PE['KF'],
                 F3=PA['SYN1'] <= 0.8 * PA['PID'] and int(np.sum(A['SYN1'] < A['PID'])) >= 8,
                 F4=PA['SYN1'] <= PA['KF'] and int(np.sum(A['SYN1'] < A['KF'])) >= 7,
                 K1=PE['SYND'] >= 0.95 * PE['SYN1'],
                 K2=tuned['c_bench'] <= 0.7 and PE['SYNB'] <= 0.9 * PE['SYN1'],
                 JPq=0.8 <= tuned['c_bench'] <= 1.25)
        reg = {a: {r: P([row['att_by_reg'][r][i] for row in rows]) for r in ('low', 'half')} for i, a in enumerate(arms)}
        S = dict(script_sha256=SHA, pooled_att=PA, pooled_est=PE, wins=dict(SYN1_vs_PID=int(np.sum(A['SYN1'] < A['PID'])), SYN1_vs_KF=int(np.sum(A['SYN1'] < A['KF']))),
                 by_regime=reg, jitter={a: P([r['jit'][i] for r in rows]) for i, a in enumerate(arms)}, verdict={k: bool(v) for k, v in V.items()},
                 tuned=tuned, gains={a: load_json(f'gains_{a}.json')['gains'] for a in ('PID', 'KF', 'SYN1')})
        json.dump(S, open(os.path.join(OUT, 'summary.json'), 'w'), indent=1); print(json.dumps(S, indent=1))
