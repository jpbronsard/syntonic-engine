"""STEP 6g, pre-registered (30 October 2026). The attitude estimator on REAL vibration noise, synthetic motion.

Origin. Every step-6 bench used simulated noise. The one real-noise test so far (the Gauge, 29 Oct) flagged the raw
1 kHz gyro out of class on every segment, idle included. Before the gimbal (cradle) bench: does the law, reading the
accel noise live and the gyro's angle random walk as a fact, land on the oracle window when the noise is the X500's
own vibration, and is a pre-filter needed?

Measured facts (looked at before writing this, no score computed): the accelerometer at the logger's 250 Hz has no
line to notch; its vibration is broadband (variance spread evenly from 20 Hz to Nyquist, aliased), per-sample tilt
noise 5.7 deg at idle, 5.6 deg at low throttle, 50 deg at half throttle. The gyro's vibration sits in a cluster at
250 to 500 Hz, not at one line. So the planned "notch" arm is replaced by a gyro low-pass arm, and the central
question becomes whether broadband vibration is white at the engine's scales (32, 64, 128 ms block means).

Data. jpbronsard/syntonic-engine bench/software_2026-10-26/data (raw SD logs, 250 Hz rows, 4 gyro sub-samples per row).
  S6-idle   sess6.csv.gz   65 to 100 s  (35 s, motors idling ~0.56 A, ballasted frame still)   reported, and the gyro fact
  S6-low    sess6.csv.gz    5 to  45 s  (40 s, ~0.92 A, throttle 0.11, held)                   judged
  S5-half   session 5       15 to  23 s  ( 8 s, ~9.8 A, throttle 0.44, held)                    judged where the record allows
  S7-flight crash7.csv.gz   19 to 24.5 s (5.5 s airborne; real motion in the streams)          reported only
Streams, roll axis: tilt_real = atan2(ay, az); w_real = mean of the 4 gyro sub-samples per row (exact angle increment),
gx / 16.384 deg/s. DT = 4 ms. Both treated as noise (the held frames' own small motion counts as noise).

Synthetic truth, as in step 6: theta = Ornstein-Uhlenbeck, sd 5 deg, correlation time 2 s, seeds 600..609 (seeds vary the
truth only; the noise is one real record). gyro = dtheta/dt + w_real; accel = theta + tilt_real; truth = theta + mean(tilt_real).
The atan2 tilt is bounded and nonlinear at half throttle; adding theta in the angle domain is a linearization (declared).

Candidate SYN-6g = SYN-6B unchanged, at DT = 4 ms: Mahony with Ki = Kp^2/(4 zeta^2), zeta = 1.826; sigma^2 live from
second differences of block means of z_raw = accel - open-loop raw gyro integral at 32, 64, 128 ms (engine block EMA,
s2 = min over scales of D_H(l) l / 6); lambda_hat = the real gyro's angle random walk from S6-idle, Allan variance at
l in {0.128, 0.256, 0.512, 1.024} s, lambda_hat = mean(sigma_A^2(l) l); tau = sqrt(s2 / lambda_hat), floor DT/0.45,
ceiling 30 s, rate limit 0.4 s/s, warm-up 1 s, start 0.5 s. Variants lambda_hat x2 and /2 reported (the fact is measured
with motors idling, not off).
Comparators, same streams: ORACLE = Mahony with a fixed window on a 40-point grid 0.02 to 50 s (parabolic refinement),
best per seed; KF-OPT = two-state Kalman (angle, gyro bias) tuned per seed on the scored data over Q_G in lambda_hat x
{1/4, 1/2, 1, 2, 4} and R in var(tilt_real) x 2^(-8..8 by 0.5) (clairvoyant: the ceiling).
Scoring: skip = min(10 s, 30 % of the record); steady RMS error in degrees after the skip; window = median tau after the skip.
Reading rule for short records: a window above 25 % of the record length is "at the record limit"; the x1.5 window
criterion is then not judged (reported), only the regret is.

Hypotheses (cells judged: S6-low and S5-half; S6-idle and S7-flight reported):
  H6g.1 class at engine scales: on tilt_real alone, ratio max/min of D_H(l) l / 6 across l in {32, 64, 128} ms <= 1.3
        (white at these scales) in S6-low and S5-half. Prediction: holds (broadband vibration averages white).
  H6g.2 the law lands: SYN-6g regret = RMS(SYN) / RMS(ORACLE) - 1 <= 10 % (seed median), and, where the oracle is not at
        the record limit, window within x1.5 of the oracle, in S6-low and S5-half.
  H6g.3 parity: mean RMS(SYN-6g) / mean RMS(KF-OPT) <= 1.10 in S6-low and S5-half.
  H6g.4 gyro pre-filter: with a causal 2nd-order Butterworth low-pass at 150 Hz on the 1 kHz gyro before the row
        average, SYN-6g's steady RMS changes by less than 5 % in both judged cells. Prediction: holds (the 250-500 Hz
        cluster averages out in the 4 ms integration). If it improves by more than 5 %, the engine should carry it.
Reported, not judged: lambda_hat against the IMU data sheet; Allan slope; S6-idle and S7-flight metrics; transients;
window traces; the x2 and /2 variants; the KF-OPT tuning found.
Declared limits: one axis; held frames, not flight dynamics; 8 s at half throttle; the gyro fact measured with motors
idling; one noise realization per segment. The cradle supersedes all of this. A failure is published as such.
Smoke policy: STEP6G_SMOKE=1 runs 2 seeds and prints lambda_hat and the runtime only, never a score.
"""
import hashlib, json, os, sys, time
import numpy as np, pandas as pd
from scipy.signal import butter, lfilter

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, 'data')   # pre-registered copy (sha 54f0c7fc) pointed at the author's scratchpad layout; this is the only change
OUT = os.path.join(HERE, 'results'); os.makedirs(OUT, exist_ok=True)
SMOKE = bool(os.environ.get('STEP6G_SMOKE'))
DT = 4e-3; D2R = np.pi / 180; LSB = 16.384
SEEDS = [600, 601] if SMOKE else list(range(600, 610))
TRUTH_SD, TRUTH_TC = 5.0 * D2R, 2.0
ZETA = 2.0 / (2.0 * np.sqrt(0.3)); Q_B = (0.005 * D2R) ** 2 / 100.0
EMA, WARMUP, TAU_RATE = 1.0, 1.0, 0.4
FLOOR_TAU, TAU_MAX, TAU0 = DT / 0.45, 30.0, 0.5
HSCALES = np.array([8, 16, 32])                      # samples at 4 ms = 32, 64, 128 ms
ALLAN_L = np.array([0.128, 0.256, 0.512, 1.024])
TAU_GRID = np.logspace(np.log10(0.02), np.log10(50.0), 40)
R_MULT = 2.0 ** np.arange(-8, 8.5, 0.5); QG_MULT = np.array([0.25, 0.5, 1.0, 2.0, 4.0])
SEGS = [('S6-idle', 'sess6.csv.gz', 65, 100, 'reported'), ('S6-low', 'sess6.csv.gz', 5, 45, 'judged'),
        ('S5-half', 'syntonic_00000000_000000_00005.csv.gz', 15, 23, 'judged'), ('S7-flight', 'crash7.csv.gz', 19, 24.5, 'reported')]
LOG = []
def P(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); LOG.append(s)


def load(name, f, t0, t1, gyro_lp=False):
    df = pd.read_csv(os.path.join(DATA, f)); t = (df.t_us.values - df.t_us.values[0]) * 1e-6; m = (t >= t0) & (t < t1); d = df[m]
    ay, az = d.ay_raw.values.astype(float), d.az_raw.values.astype(float)
    tilt = np.arctan2(ay, az)
    g1k = np.stack([d[f'gx_s{k}'].values for k in range(4)], 1).reshape(-1) / LSB * D2R      # 1 kHz, rad/s
    if gyro_lp:
        b, a = butter(2, 150.0 / 500.0); g1k = lfilter(b, a, g1k)
    w = g1k.reshape(-1, 4).mean(1)
    return tilt, w, g1k


def allan_lambda(g1k):
    N = len(g1k)
    def av(l):
        k = int(round(l / 1e-3)); mb = g1k[: N // k * k].reshape(-1, k).mean(1); return 0.5 * np.mean(np.diff(mb) ** 2)
    lam = np.mean([av(l) * l for l in ALLAN_L])
    Lfit = np.array([0.016, 0.032, 0.064, 0.128, 0.256, 0.512, 1.024, 2.048])
    slope = np.polyfit(np.log(Lfit), np.log([av(l) for l in Lfit]), 1)[0] / 2
    return lam, slope


def class_ratio(x):
    """ratio max/min of D_H(l) l / 6 across the engine scales on one real stream (white: near 1)"""
    s2 = []
    for s in HSCALES:
        mb = x[: len(x) // s * s].reshape(-1, s).mean(1); dd = mb[2:] - 2 * mb[1:-1] + mb[:-2]
        s2.append(np.mean(dd ** 2) * (s * DT) / 6.0)
    return max(s2) / min(s2), s2


class HadNoise:
    def __init__(self, n):
        m = len(HSCALES); self.acc = np.zeros((m, n)); self.cnt = np.zeros(m, dtype=int)
        self.p1 = np.full((m, n), np.nan); self.p2 = np.full((m, n), np.nan); self.DH = np.zeros((m, n)); self.nb = np.zeros(m, dtype=int)
        self.a_blk = np.array([min(1.0, DT * s / max(EMA, 8 * DT * s)) for s in HSCALES])
    def update(self, z):
        self.acc += z[None, :]; self.cnt += 1
        for i in np.flatnonzero(self.cnt == HSCALES):
            mean = self.acc[i] / HSCALES[i]
            if not np.isnan(self.p2[i, 0]):
                dd = mean - 2 * self.p1[i] + self.p2[i]; aw = self.a_blk[i] if self.nb[i] > 0 else 1.0
                self.DH[i] += aw * (dd * dd - self.DH[i]); self.nb[i] += 1
            self.p2[i] = self.p1[i]; self.p1[i] = mean; self.acc[i] = 0.0; self.cnt[i] = 0
        if (self.nb >= 2).all(): return np.min(self.DH * (HSCALES * DT)[:, None] / 6.0, axis=0)
        return None


class LawTau:
    def __init__(self, n, lam):
        self.lam = lam; self.tau = np.full(n, TAU0); self.tau_t = np.full(n, TAU0); self.k = 0; self.warm = int(round(WARMUP / DT))
    def update(self, s2):
        if self.k > self.warm:
            if s2 is not None: self.tau_t = np.clip(np.sqrt(np.maximum(s2, 1e-300) / self.lam), FLOOR_TAU, TAU_MAX)
            rr = 1.0 / self.tau; dr = TAU_RATE * DT / self.tau ** 2
            rr = rr + np.clip(1.0 / self.tau_t - rr, -dr, dr); self.tau = 1.0 / rr
        self.k += 1
        return self.tau


class Mahony:
    def __init__(self, a0):
        self.th = a0.copy(); self.b = np.zeros_like(a0); self.Th = a0.copy()
    def step(self, w, a, tau):
        kp = 1.0 / tau; wb = w - self.b
        self.th = self.th + wb * DT; self.Th = self.Th + wb * DT
        e = a - self.th; self.th = self.th + np.minimum(DT * kp, 1.0) * e
        ki = kp * kp / (4 * ZETA ** 2); self.b = np.clip(self.b - ki * e * DT, -0.15, 0.15)
        return self.th


class Kalman2:
    def __init__(self, a0, R, QG):
        n = a0.shape[0]; self.x1 = a0.copy(); self.x2 = np.zeros(n); self.QG = QG
        self.p11 = np.asarray(R, float) * np.ones(n); self.p12 = np.zeros(n); self.p22 = np.full(n, (0.5 * D2R) ** 2)
    def step(self, w, a, R):
        self.x1 = self.x1 + (w - self.x2) * DT
        p11 = self.p11 - 2 * DT * self.p12 + DT * DT * self.p22 + self.QG * DT
        p12 = self.p12 - DT * self.p22; p22 = self.p22 + Q_B * DT
        y = a - self.x1; S = p11 + R; k1 = p11 / S; k2 = p12 / S
        self.x1 = self.x1 + k1 * y; self.x2 = self.x2 + k2 * y
        self.p11 = (1 - k1) * p11; self.p12 = (1 - k1) * p12; self.p22 = p22 - k2 * p12
        return self.x1


def truth(N, seeds):
    rho = np.exp(-DT / TRUTH_TC); c = TRUTH_SD * np.sqrt(1 - rho ** 2); TH = np.empty((N + 1, len(seeds)))
    for j, s in enumerate(seeds):
        r = np.random.default_rng(s); TH[:, j] = lfilter([c], [1, -rho], r.standard_normal(N + 1))
    return TH


def refine(lg, y):
    i = int(np.argmin(y))
    if 0 < i < len(y) - 1:
        y0, y1, y2 = y[i - 1:i + 2]; den = y0 - 2 * y1 + y2; off = 0.5 * (y0 - y2) / den if den > 0 else 0.0
        return np.exp(lg[i] + off * (lg[i] - lg[i - 1])), y[i]
    return np.exp(lg[i]), y[i]


def run_segment(name, f, t0, t1, lam_hat, gyro_lp=False):
    tilt, w_real, _ = load(name, f, t0, t1, gyro_lp); N = len(tilt); ns = len(SEEDS)
    skip = int(round(min(10.0, 0.3 * N * DT) / DT)); TH = truth(N, SEEDS)
    W = np.diff(TH, axis=0) / DT + w_real[:, None]; A = TH[1:] + tilt[:, None]; TRU = TH[1:] + tilt.mean()
    # law, three lambda variants
    var = [('6g', 1.0), ('6g_x2', 2.0), ('6g_d2', 0.5)]; nv = len(var); rep = lambda x: np.tile(x, nv)
    lam = np.concatenate([np.full(ns, lam_hat * fct) for _, fct in var])
    mah = Mahony(rep(A[0])); had = HadNoise(nv * ns); law = LawTau(nv * ns, lam); Th_raw = A[0].copy()
    # oracle grid and KF grid, vectorized
    nt = len(TAU_GRID); orc = Mahony(np.tile(A[0], nt)); tau_orc = np.repeat(TAU_GRID, ns)
    Rg = tilt.var() * R_MULT; nR, nQ = len(Rg), len(QG_MULT); nk = nR * nQ
    kf = Kalman2(np.tile(A[0], nk), np.repeat(Rg, nQ * ns), np.repeat(np.repeat(QG_MULT * lam_hat, ns)[None, :], nR, 0).reshape(-1))
    R_kf = np.repeat(Rg, nQ * ns)
    ss = {k: 0.0 for k in ('law', 'orc', 'kf')}; n_sc = 0; tr = {k: 0.0 for k in ss}; n_tr = 0; tlog = []
    for k in range(N):
        w, a, tru = W[k], A[k], TRU[k]
        Th_raw = Th_raw + w * DT; s2 = had.update(rep(a - Th_raw)); tau = law.update(s2)
        e_law = mah.step(rep(w), rep(a), tau) - rep(tru)
        e_orc = orc.step(np.tile(w, nt), np.tile(a, nt), tau_orc) - np.tile(tru, nt)
        e_kf = kf.step(np.tile(w, nk), np.tile(a, nk), R_kf) - np.tile(tru, nk)
        if k >= skip:
            ss['law'] = ss['law'] + e_law ** 2; ss['orc'] = ss['orc'] + e_orc ** 2; ss['kf'] = ss['kf'] + e_kf ** 2; n_sc += 1
            if k % 25 == 0: tlog.append(law.tau[:ns].copy())
        else:
            tr['law'] = tr['law'] + e_law ** 2; n_tr += 1
    rms = {k: np.sqrt(v / n_sc) / D2R for k, v in ss.items()}
    law_rms = rms['law'].reshape(nv, ns); orc_rms = rms['orc'].reshape(nt, ns); kf_rms = rms['kf'].reshape(nR, nQ, ns)
    lg = np.log(TAU_GRID); o_tau = np.empty(ns); o_rms = np.empty(ns)
    for j in range(ns): o_tau[j], o_rms[j] = refine(lg, orc_rms[:, j])
    kf_best = kf_rms.reshape(-1, ns).min(0); kf_arg = kf_rms.reshape(-1, ns).argmin(0)
    tau_med = np.median(np.array(tlog), axis=0) if tlog else np.full(ns, np.nan)
    cr, s2l = class_ratio(tilt)
    rec = N * DT
    out = dict(name=name, record_s=rec, skip_s=skip * DT, tilt_sd_deg=float(tilt.std() / D2R), gyro_sd_dps=float(w_real.std() / D2R),
               class_ratio=float(cr), s2_scales=[float(x) for x in s2l],
               law_rms=law_rms.tolist(), law_tau_med=tau_med.tolist(), orc_tau=o_tau.tolist(), orc_rms=o_rms.tolist(),
               kf_rms=kf_best.tolist(), kf_R_mult=[float(R_MULT[i // nQ]) for i in kf_arg], kf_QG_mult=[float(QG_MULT[i % nQ]) for i in kf_arg],
               trans_rms=(np.sqrt(tr['law'] / max(n_tr, 1)) / D2R).reshape(nv, ns).tolist(), gyro_lp=gyro_lp)
    return out


def main():
    sha = hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest()
    P('step6g_real_noise.py sha256', sha, '| smoke' if SMOKE else '| FULL', '| seeds', SEEDS)
    t0 = time.time()
    _, _, g_idle = load(*SEGS[0][:4]); lam_hat, slope = allan_lambda(g_idle)
    ds = (0.0028 * D2R) ** 2   # ICM-42688-P data-sheet ARW 0.0028 deg/s/sqrt(Hz) (reference only)
    P(f'gyro fact (S6-idle, motors idling): lambda_hat = {lam_hat:.3e} rad^2/s = ({np.sqrt(lam_hat)/D2R:.4f} deg/s/sqrt(Hz))^2; '
      f'data-sheet ratio {lam_hat/ds:.1f}; Allan-deviation slope {slope:.2f} (ARW: -0.5)')
    res = dict(sha256=sha, smoke=SMOKE, seeds=SEEDS, lambda_hat=lam_hat, allan_slope=slope, segments={})
    for name, f, a, b, role in SEGS:
        for lp in (False, True):
            r = run_segment(name, f, a, b, lam_hat, gyro_lp=lp); r['role'] = role; res['segments'][name + ('|lp' if lp else '')] = r
            P(f'{name}{" +gyro-lp" if lp else ""}: done ({time.time()-t0:.0f} s)')
    P(f'runtime {time.time()-t0:.0f} s')
    fn = os.path.join(OUT, 'step6g_smoke.json' if SMOKE else 'step6g_rows.json'); json.dump(res, open(fn, 'w'), indent=1); P('written', fn)
    if SMOKE: return
    # ------------------------------------------------------------------ verdicts
    P('\n=== STEP 6g verdicts (seed medians for windows and regret, means for RMS ratios) ===')
    P(f'{"segment":10s} {"rec":>5s} {"tiltSD":>6s} {"class":>6s} {"law tau":>8s} {"orc tau":>8s} {"regret":>7s} {"RMS law":>7s} {"RMS orc":>7s} {"RMS KF":>7s} {"law/KF":>6s} {"lp dRMS":>7s}')
    H1 = H2 = H3 = H4 = True
    for name, f, a, b, role in SEGS:
        r = res['segments'][name]; rl = res['segments'][name + '|lp']
        law = np.array(r['law_rms'][0]); orc = np.array(r['orc_rms']); kf = np.array(r['kf_rms'])
        regret = np.median(law / orc - 1); ltau = np.median(r['law_tau_med']); otau = np.median(r['orc_tau'])
        ratio_kf = law.mean() / kf.mean(); d_lp = np.array(rl['law_rms'][0]).mean() / law.mean() - 1
        at_limit = otau > 0.25 * r['record_s']
        P(f'{name:10s} {r["record_s"]:5.1f} {r["tilt_sd_deg"]:6.1f} {r["class_ratio"]:6.2f} {ltau:8.2f} {otau:8.2f}{"*" if at_limit else " "} {regret*100:6.1f}% {law.mean():7.3f} {orc.mean():7.3f} {kf.mean():7.3f} {ratio_kf:6.3f} {d_lp*100:+6.1f}%   [{role}]')
        if role == 'judged':
            H1 &= r['class_ratio'] <= 1.3
            H2 &= (regret <= 0.10) and (at_limit or (1 / 1.5 <= ltau / otau <= 1.5))
            H3 &= ratio_kf <= 1.10
            H4 &= abs(d_lp) < 0.05
    P('  * oracle window above 25 % of the record: window criterion not judged, regret only')
    P(f'\nH6g.1 white at engine scales (class ratio <= 1.3, judged cells): {"PASS" if H1 else "FAIL"}')
    P(f'H6g.2 the law lands (regret <= 10 %, window x1.5 where judgeable): {"PASS" if H2 else "FAIL"}')
    P(f'H6g.3 parity with tuned Kalman (<= 1.10): {"PASS" if H3 else "FAIL"}')
    P(f'H6g.4 gyro low-pass changes RMS by < 5 %: {"PASS" if H4 else "FAIL"}')
    P('\nReported: lambda variants (RMS mean, x2 / /2 relative to 6g):')
    for name, *_ in SEGS:
        r = res['segments'][name]; L = np.array(r['law_rms']); P(f'  {name:10s} x2 {L[1].mean()/L[0].mean():.3f}  /2 {L[2].mean()/L[0].mean():.3f}   KF-OPT tuning (R mult, QG mult) medians {np.median(r["kf_R_mult"]):.3g}, {np.median(r["kf_QG_mult"]):.3g}   s2 by scale {np.round(r["s2_scales"], 9).tolist()}')
    open(os.path.join(OUT, 'step6g_results.txt'), 'w').write('\n'.join(LOG) + '\n')


if __name__ == '__main__':
    main()
