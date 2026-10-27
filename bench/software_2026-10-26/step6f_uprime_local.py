"""STEP 6f, pre-registered. Theorem U'-Local against the square root and against U' read with a raw local index: three
predictions frozen before the evaluation, the exponent of the optimal window as the main test, the window and its
regret as the second, and noise levels that keep the noise's spectral shape.

Origin. The note Theorem U'-Local (v1.0, deposit pending) gives the exact local exponent of the optimal window,
alpha_loc = 1 / (s[A'] - s[-beta']) at the optimum, for every separable loss whose balance increases from 0 to infinity;
for the exponential window under white noise the optimum is the root of tau E[(E - 1) D(tau E)] = sigma^2 (its eq. 12).
Its Section 11 pre-registers P1 (exact local law) and P2 (slope bias) for a bench of this estimator with a two-regime
gyro drift; this is that bench. Design by J.-P. Bronsard (2026-10-27): three predictors frozen before any result, the
third computing its slopes from calibration statistics without seeing the evaluation's optimal windows; the prediction
error of the exponent on held-out data as the main test; the optimal window and its regret as the second; changes of
sigma^2 that keep the noise's spectral shape, so that the measured exponent is the partial derivative the theorem
defines.

Setting. One axis; the attitude estimator of steps 6d and 6e inside the note's class: Mahony without bias integral, a
first-order window of gain DT / tau on the accelerometer, the gyro carrying the drift (UL-H3) and the accelerometer the
noise (UL-H4). The filter is linear and the attitude truth cancels, so the error at noise level c is e_d + sqrt(c) e_n
(e_d from the drift alone, e_n from the unit noise alone) and the steady mean-square error per seed is
A + c B + 2 sqrt(c) X, exactly: every level comes from one drift path and one unit noise path per seed, and the noise
keeps its shape across levels by construction. DT = 1 ms; paths of 600 s; steady part after max(10 tau, 1 s);
windows: 40 values from 0.03 to 30 s, log-spaced. Q = Q_G = (0.01 deg/s/sqrt(Hz))^2, step 6's angle random walk.

Arms (drift structure functions D, accelerometer noise):
  C50  D = Q s                           angle random walk.                               Control: all predictors agree.
  C90  D = Q s^1.8                       fractional drift, H = 0.9.        Control: U' agrees with any reading; sqrt not.
  X9   D = Q (s + s^1.8)                 angle random walk plus a fractional component of index 0.9, crossover at 1 s:
                                         the note's P1-P2 bench.                                            Discriminant.
  X10  D = Q s + Q s^2                   angle random walk plus a constant random bias of standard deviation sqrt(Q)
                                         rad/s (0.01 deg/s), crossover at 1 s; wider index spread.         Discriminant.
  J    D = Q s + (100 Q / 3)(1 - cos 3s) angle random walk plus an angle oscillation of random phase at 3 rad/s: the
                                         note's Remark 6 (balance not monotone; the optimum jumps).
  These five: white noise of intensity sigma^2 = c Q, the same unit path scaled at every level. Levels c = 10^k, k from
  -2 in steps of 1/2 to 1.5 (C50) or 3 (C90, X9, X10): optima between 0.08 and 8 s. J: sigma^2 = s2c x (0.6, 0.75,
  0.85, 0.92, 0.96, 0.98, 1.02, 1.04, 1.08, 1.15, 1.3, 1.6), s2c the jump level predicted by stage 1.
  S50, S90  step 6e's drifts (fBm, H = 0.5 and 0.9, lambda_H such that U''s window at V = 1 is 1.077 s) with step 6e's
  accelerometer noise at constant V = 1, 256, 1/16 (white vibration V Q, a line of amplitude 1 deg x sqrt(V) at
  92 V^(1/4) Hz, white sensor noise, 80 Hz Butterworth; the same draws at every V). Its shape changes with V: the
  control's negative, what a bench measures when the shape is not kept.

Predictors. All receive the same inputs: D measured at rest on calibration seeds (400 records of 1200 s, seeds
3000..3399, pooled; by FFT, all lags to 400 s) and sigma^2 at each level (known: c Q, or S0(V)). They differ only in
the law; nothing from the evaluation enters.
  LOCAL   U'-Local, Theorem 1 on the discrete filter: the window minimising A(tau) + sigma^2 / (2 tau - DT), with
          A(tau) = a (1 - a) / (2 - a) sum_(j >= 1) (1 - a)^j D(j DT), a = DT / tau, the exact discrete counterpart of
          the note's Lemma 4 (the noise term is the discrete one too). Exponents: secants of its windows between
          adjacent levels. For the tests (H6f.1, H6f.2, H6f.3, H6f.5) its predicted curve is read the way the bench reads
          the measured one, on the evaluation's 40-point grid with the same parabola (LOCAL-G; AWARE-G likewise): read
          on the fine grid, it would differ from the grid fit by up to 0.003 in exponent, a fixed bias that exceeds 3 SE
          at large N; read the same way, the two agree to 3e-5 on noise-free curves (checked).
  RAWS    U' with the raw local index (the design's middle row): exponent 1 / (2 H_D + 1) for each pair of adjacent
          levels, H_D half the log-slope of D at the scale where it is read, the geometric mean of RAWF's two windows.
  SQRT    the square root: exponent 1/2; window sqrt(sigma^2 / lambda), lambda the mean of D(l) / l, l = 0.128 .. 1.024 s.
  Reported, not judged:
  RAWF    U' applied level by level with the power law fitted to D at the window itself, 2 H Gamma(2H + 1) tau D(tau)
          = sigma^2 with H = H_D(tau). It agrees with LOCAL to first order in the spread of the indices (within 0.010 on
          X9 and 0.016 on X10, computed below): this bench cannot tell them apart and does not try.
  ENGINE  the engine's U' law with H and lambda_H from one block-mean fit at rest over 0.128 .. 8.192 s (step 6e's
          procedure, pooled): one exponent for all levels.
  LOCAL-C the root of the continuous law (eq. 12) by quadrature on the measured D, the theorem as written for a continuous
          window, its first root; in C50, C90, X9 and X10 at most 0.7 % below LOCAL's window (discretisation); in J the
          first root sits on the lower branch above the jump.
  Shape arms: LOCAL is shape-blind (white noise of intensity S0(V)); AWARE minimises A(tau) + N_V(tau), N_V the
  noise-only error curve of the filter measured on calibration noise records at each V (seeds 3000..3099, 600 s).

Expected values from the true structure functions, computed before any run (the frozen predictions come from the
calibration and will differ by its sampling):
  C50  LOCAL exponents 0.498 .. 0.500; all predictors within 0.002.     C90  LOCAL 0.356 .. 0.357; SQRT 0.5.
  X9   LOCAL 0.438 0.425 0.412 0.400 0.390 0.382 0.376 0.371 0.367 0.364
       RAWS - LOCAL 0.034 0.037 0.038 0.037 0.034 0.030 0.025 0.020 0.016 0.012 (RMS 0.030); SQRT - LOCAL RMS 0.110;
       ENGINE - LOCAL -0.032 .. +0.042 (H = 0.73); RAWF - LOCAL +0.010 .. -0.006.
  X10  LOCAL 0.431 0.412 0.394 0.378 0.366 0.356 0.349 0.344 0.341 0.339
       RAWS - LOCAL 0.044 0.050 0.053 0.051 0.046 0.039 0.031 0.024 0.018 0.013 (RMS 0.040); SQRT - LOCAL RMS 0.133;
       ENGINE - LOCAL -0.044 .. +0.048 (H = 0.79); RAWF - LOCAL +0.016 .. -0.010 (-0.0105).
  J    s2c = 7.26 Q; LOCAL windows 0.270 .. 0.456 s below the jump, 1.943 .. 3.033 s above, a ratio of 4.26 between the
       levels 0.98 and 1.02; RAWF and SQRT: no adjacent ratio above 1.12.
  S50  down exponent LOCAL (blind) 0.500, AWARE 0.471; up 0.500 and 0.499.   S90  down 0.357 and 0.341; up 0.357, 0.356.

Stages. (1) calibrate: the predictions and their calibration bootstrap (2000 draws over calibration seeds), frozen in
results/step6f_predictions.json and step6f_pred_draws.npz and committed before stage 2. (2) pilot: 100 seeds
(2900..2999) through the evaluation; it prints only standard errors, grid-edge flags, the identity check and N, the
number of evaluation seeds: the smallest multiple of 100 in [200, 2000] such that three standard errors (calibration
part, plus the evaluation part scaled to N) fit inside the predicted gap of H6f.2 (X9 and X10, the RMS of RAWS - LOCAL)
and of H6f.5 (S50, LOCAL - AWARE). A hypothesis whose gap is not reachable at N = 2000 is declared non-discriminant
there: reported, not judged. With the computed gaps (0.030 to 0.040, and 0.029) against the standard errors seen in
toy runs, N is expected at the floor of 200. Committed before stage 3; the judge checks that its N and the predictions'
hash match the pilot file. (3) evaluate: N seeds from 4000; rows in
results/step6f_rows.npz; then the judge.

Measurement. Population optimum at each level: the seed-mean mean-square error on the window grid, minimised by a
parabola in log tau through the grid minimum and its neighbours (step 6e). Exponents: secants of log tau_opt between
adjacent levels; shape arms 1 -> 256 (up) and 1 -> 1/16 (down) against S0(V), as in step 6e. Regret of a window:
sqrt(J(tau) / J_min) - 1, both read on a cubic spline in log tau of the population curve. Uncertainty: joint bootstrap, each draw pairing an evaluation draw
(over evaluation seeds, all levels together) with a calibration draw; 2000 draws; 95 % percentile intervals.

Hypotheses and decision rules (judge: evaluate_step6f.py, written with this script, before any run):
  H6f.0 validity: the decomposition equals step 6's Mahony (no bias integral, truth included) to 1e-9 rad, and no
        population optimum lies at a grid edge. A failure voids the run; nothing else is judged.
  H6f.1 U'-Local exact (the note's P1): in C50, C90, X9, X10, every secant has |e - e_LOCAL-G| <= 3 SE. With 37
        secants, the family's false-alarm rate is about 10 % if the law is exact.
  H6f.2 U'-Local beats the raw local index (the main comparison; the note's P2): in X9 and X10,
        d = RMS(e - e_RAWS) - RMS(e - e_LOCAL-G) > 0 with its 95 % interval above 0.
  H6f.3 the square root is rejected: in C90, X9, X10, RMS(e - 1/2) - RMS(e - e_LOCAL-G) has its interval above 0.
  H6f.4 the window and its regret (the second test): in C50, C90, X9, X10, at every level, LOCAL's window within x1.1
        of tau_opt and its regret at most 1 %.
  H6f.5 the noise shape (the control): in S50 the down exponent lies within 3 SE of AWARE-G and more than 3 SE below
        LOCAL-G (shape-blind). This also tests the exploratory reading of 2026-10-27 that the shape change accounts for
        part of step 6e's low down exponents: 0.471 predicted here against 6e's 0.454 at H = 0.5 (at most 0.029 of the
        0.046 shortfall), 0.341 against 0.340 at H = 0.9; step 6e switched V inside 180 s runs, this bench holds V
        constant. S90 reported.
  H6f.6 the jump (Remark 6): in J, exactly one adjacent ratio of tau_opt is at least 2 and no other reaches 1.5; it lies
        between the relative levels 0.92 and 1.08; every tau_opt is within x1.1 of LOCAL's minimiser on the branch it
        is on (the lower or upper minimum of A + sigma^2 beta, below or above 1 s).
  Reported, not judged: RAWF, ENGINE, LOCAL-C and the window and regret of every predictor (regret read on a cubic
  spline of the population curve); the white-noise check of B; the shape arms' noise curves against S0(V); the up
  exponents; S90; the smooth predictors' largest ratio in J; the number of non-finite bootstrap draws.
Reading rules, declared. H6f.1 says U'-Local, with D measured once at rest, predicts the exponents of the optimal window
within resolution, crossovers included, in this generator; if it fails only by deviations under 0.005 (the judge prints
whether all are), the failure is reported as a bench bias of that size (the verdict stays FAILS). H6f.2 with H6f.1 says the theorem's correction over the
raw slope reading is real at the resolved level. H6f.3 rejects the square root wherever the drift is not Brownian.
H6f.4 says the predicted window is usable. H6f.5 says the control matters: when the noise changes shape with its level,
the theorem's exponent is not what a bench measures, and the shape-aware prediction is. H6f.6 says the optimum jumps
where the balance stops being monotone, at the predicted level and onto the predicted branches.
Declared limits. One simulated axis; drifts generated from the declared laws, so the bench tests the mathematics, the
calibration chain and the code, not a real gyroscope; the filter is the class's (no bias integral), not the flown
Mahony; sigma^2 is known at each level (the live noise estimate, and the note's P3, are not tested here); RAWF is not
judged. A failure is published as such.

Smoke policy: STEP6F_SMOKE=1 runs every stage at toy size (2 calibration seeds of 240 s, 4 seeds of 120 s, windows to
10 s) to check that the code runs, the identity of the decomposition and the FFT structure function against direct
sums; its outputs (suffix _smoke) are not results and are not kept, apart from results/step6f_smoke.txt. Any change
after the rules are committed is written here as an amendment before stage 1.
"""
import hashlib, json, math, os, sys, time
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from scipy.integrate import quad
from scipy.signal import butter, lfilter
from scipy.special import gamma as G

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import step6_quad_attitude_estimator as S6
import step6d_uprime_mechanism as S6D

SMOKE = bool(os.environ.get('STEP6F_SMOKE'))
DT, Q = S6.DT, S6.Q_G
RES = os.path.join(HERE, 'results')

# ------------------------------------------------------------------------------------------------ arms and constants
# drift components: ('bm', lam) D = lam s; ('fbm', H, lam) D = lam s^2H; ('ramp', v) D = v s^2 (a constant random bias of
# variance v); ('sin', a2, w) D = a2 (1 - cos w s) (an angle oscillation of random phase)
ARMS = {
    'C50': dict(kind='scaled', comps=[('bm', Q)]),
    'C90': dict(kind='scaled', comps=[('fbm', 0.9, Q)]),
    'X9': dict(kind='scaled', comps=[('bm', Q), ('fbm', 0.9, Q)]),
    'X10': dict(kind='scaled', comps=[('bm', Q), ('ramp', Q)]),
    'J': dict(kind='scaled', comps=[('bm', Q), ('sin', 100.0 * Q / 3.0, 3.0)]),
    'S50': dict(kind='shape', comps=[('bm', S6D.lam_for(0.5))]),
    'S90': dict(kind='shape', comps=[('fbm', 0.9, S6D.lam_for(0.9))]),
}
ARM_ID = {a: i for i, a in enumerate(ARMS)}
# noise levels. Scaled arms: white noise of intensity sigma^2 = c Q, c = 10^k, the same unit path scaled at every level.
# J: sigma^2 = s2c x J_REL, s2c the jump level predicted from the calibration. Shape arms: step 6e's V levels.
LADDER_K = {'C50': np.arange(-4, 4) / 2, 'C90': np.arange(-4, 7) / 2, 'X9': np.arange(-4, 7) / 2,
            'X10': np.arange(-4, 7) / 2}
J_REL = np.array([0.6, 0.75, 0.85, 0.92, 0.96, 0.98, 1.02, 1.04, 1.08, 1.15, 1.3, 1.6])
V_S = [1.0, 256.0, 1.0 / 16.0]
GRID = np.logspace(np.log10(0.03), np.log10(30.0), 40)                       # evaluation windows (s)
TF = np.logspace(np.log10(0.03), np.log10(40.0), 260)                        # fine grid of the predictions (s)
T_EV, T_CAL, T_NCAL = 600.0, 1200.0, 600.0                                   # path lengths (s)
L_MAX = int(round(400.0 / DT))                                               # structure function to 400 s
L_BM = S6D.L_REST_U                                                          # block means 0.128 .. 8.192 s
L_SQ = 0.128 * 2.0 ** np.arange(4)                                           # square-root lambda, 0.128 .. 1.024 s
CAL_SEEDS = list(range(3000, 3400)); NCAL_SEEDS = set(range(3000, 3100))
PILOT_SEEDS = list(range(2900, 3000)); EV_SEED0 = 4000
B_BOOT = 2000; BATCH = 25; N_MIN, N_MAX = 200, 2000
if SMOKE:
    CAL_SEEDS, NCAL_SEEDS, PILOT_SEEDS = [3000, 3001], {3000, 3001}, [2900, 2901, 2902, 2903]
    T_EV, T_CAL, T_NCAL, B_BOOT, BATCH = 120.0, 240.0, 60.0, 20, 2
    L_MAX = int(round(60.0 / DT)); GRID = GRID[GRID <= 10.0]
BB, AB = butter(2, 80.0 / (0.5 / DT))
TAG = '_smoke' if SMOKE else ''


def rng(seed, arm, purpose, k=0):
    """independent stream per (seed, arm, purpose, component); purpose 0 calibration drift, 1 evaluation drift,
    2 evaluation noise, 3 calibration noise, 4 identity check."""
    return np.random.default_rng([seed, ARM_ID[arm], purpose, k])


def sha_of(path):
    return hashlib.sha256(open(path, 'rb').read()).hexdigest()


# ------------------------------------------------------------------------------------------------ drift and structure
def D_true(arm, s):
    s = np.asarray(s, float); out = np.zeros_like(s)
    for c in ARMS[arm]['comps']:
        if c[0] == 'bm': out = out + c[1] * s
        elif c[0] == 'fbm': out = out + c[2] * s ** (2 * c[1])
        elif c[0] == 'ramp': out = out + c[1] * s * s
        elif c[0] == 'sin': out = out + c[1] * (1 - np.cos(c[2] * s))
    return out


def drift_increments(arm, n, seed, purpose):
    """increments of the gyro drift angle over n steps of DT."""
    dx = np.zeros(n)
    for k, c in enumerate(ARMS[arm]['comps']):
        r = rng(seed, arm, purpose, k)
        if c[0] == 'bm':
            dx += math.sqrt(c[1] * DT) * r.standard_normal(n)
        elif c[0] == 'fbm':
            dx += math.sqrt(c[2]) * DT ** c[1] * S6D.fgn(n, c[1], r)
        elif c[0] == 'ramp':
            dx += math.sqrt(c[1]) * r.standard_normal() * DT
        elif c[0] == 'sin':
            t = np.arange(n + 1) * DT; dx += np.diff(math.sqrt(c[1]) * np.cos(c[2] * t + r.uniform(0, 2 * np.pi)))
    return dx


def struct_fn(x, L):
    """D(j DT) = mean over t of (x(t + j) - x(t))^2 for j = 0..L, by FFT."""
    n = len(x); x = x - x.mean(); m = 1 << int(np.ceil(np.log2(2 * n)))
    F = np.fft.rfft(x, m); ac = np.fft.irfft(F * np.conj(F), m)[: L + 1]
    c2 = np.concatenate([[0.0], np.cumsum(x * x)]); j = np.arange(L + 1)
    return np.maximum((c2[n] - c2[j]) + c2[n - j] - 2 * ac, 0.0) / (n - j)


def a_disc(Dfull, taus):
    """tracking mean-square error of the discrete filter e_k = (1 - a)(e_(k-1) + dx_k), a = DT / tau, from the structure
    function at all lags: A = a (1 - a) / (2 - a) sum_(j >= 1) (1 - a)^j D_j (the discrete Laplace representation)."""
    out = np.empty(len(taus)); L = len(Dfull) - 1
    for i, t in enumerate(taus):
        a = DT / t; jm = min(L, int(40 * t / DT)); j = np.arange(1, jm + 1)
        out[i] = a * (1 - a) / (2 - a) * (np.exp(j * math.log1p(-a)) @ Dfull[1: jm + 1])
    return out


def beta_disc(taus):
    """noise mean-square error of the discrete filter per unit of white intensity: a / ((2 - a) DT) = 1 / (2 tau - DT)."""
    return 1.0 / (2 * np.asarray(taus, float) - DT)


def filt_err(dx, nz, tau):
    """drift part and noise part of the error of the Mahony filter without bias integral (they add: the filter is
    linear and the attitude truth cancels; checked against step6_quad_attitude_estimator.Mahony by identity_check)."""
    a = DT / tau
    ed = None if dx is None else lfilter([1 - a], [1, -(1 - a)], dx)
    en = None if nz is None else lfilter([a], [1, -(1 - a)], nz)
    return ed, en


def shape_noise(n, seed, arm, purpose):
    """step 6e's accelerometer noise at constant V, the same draws at every V: white vibration of PSD V Q, a line of
    amplitude LINE_AMP sqrt(V) at 92 V^(1/4) Hz, white sensor noise S_SENS, through the 80 Hz Butterworth."""
    r = rng(seed, arm, purpose, 9); z1 = r.standard_normal(n); z2 = r.standard_normal(n); ph = r.uniform(0, 2 * np.pi)
    t = np.arange(1, n + 1) * DT; out = []
    for V in V_S:
        vib = math.sqrt(V * Q / DT) * z1 + S6.LINE_AMP * math.sqrt(V) * np.sin(2 * np.pi * S6.LINE_F * V ** 0.25 * t + ph)
        out.append(lfilter(BB, AB, vib + math.sqrt(S6.S_SENS / DT) * z2))
    return out


def unit_noise(n, seed, arm, purpose):
    """white noise of intensity Q (level c = 1); the level c uses sqrt(c) times this path."""
    return math.sqrt(Q / DT) * rng(seed, arm, purpose, 9).standard_normal(n)


# ------------------------------------------------------------------------------------------------ predictors
def argmin_log(J, grid):
    """minimiser and minimum of a curve on a log grid: parabola in log tau through the grid minimum and its neighbours;
    (nan, nan) at an edge."""
    i = int(np.argmin(J))
    if i in (0, len(J) - 1):
        return float('nan'), float('nan')
    x = np.log(grid[i - 1: i + 2]); y = J[i - 1: i + 2]; c = np.polyfit(x, y, 2)
    if c[0] <= 0:
        return float(grid[i]), float(J[i])
    xm = -c[1] / (2 * c[0]); return float(np.exp(xm)), float(np.polyval(c, xm))


def local_tau(A, s2=None, Nf=None):
    """U'-Local, Theorem 1 on the discrete filter: argmin over TF of A + s2 beta_disc (white) or A + Nf (a noise curve)."""
    return argmin_log(A + (s2 * beta_disc(TF) if Nf is None else Nf), TF)[0]


def grid_read(A, s2=None, Nf=None):
    """the predicted curve read the way the bench reads the measured one: A (and the noise curve) interpolated in log-log
    from TF to the evaluation grid GRID, then the same parabola fit on GRID."""
    lg, lt = np.log(GRID), np.log(TF); AG = np.exp(np.interp(lg, lt, np.log(A)))
    NG = s2 * beta_disc(GRID) if Nf is None else np.exp(np.interp(lg, lt, np.log(Nf)))
    return argmin_log(AG + NG, GRID)[0]


def make_Dl(Dm, D0, Dp):
    """D at a lag in [TF[0] / 1.1, TF[-1] * 1.1], log-log interpolated between the stored lags."""
    lt = np.log(TF); xs = np.concatenate([lt - math.log(1.1), lt, lt + math.log(1.1)])
    ys = np.log(np.maximum(np.concatenate([Dm, D0, Dp]), 1e-300)); o = np.argsort(xs, kind='stable'); xs, ys = xs[o], ys[o]
    return lambda l: float(math.exp(np.interp(math.log(l), xs, ys)))


def hd_at(Dl, t):
    """half the log-slope of D at lag t, from D at t / 1.1 and t * 1.1."""
    return 0.5 * math.log(Dl(t * 1.1) / Dl(t / 1.1)) / math.log(1.21)


def first_crossing(R, s2):
    for i in range(len(TF) - 1):
        if np.isfinite(R[i]) and np.isfinite(R[i + 1]) and 0 < R[i] < s2 <= R[i + 1]:
            f = (math.log(s2) - math.log(R[i])) / (math.log(R[i + 1]) - math.log(R[i]))
            return float(math.exp(math.log(TF[i]) + f * math.log(TF[i + 1] / TF[i])))
    return float('nan')


def rawf_curve(Dl):
    R = []
    for t in TF:
        h = hd_at(Dl, t); R.append(2 * h * G(2 * h + 1) * t * Dl(t) if h > 0 else float('nan'))
    return np.array(R)


def engine_fit(Dbm):
    sl, ic = np.polyfit(np.log(L_BM), np.log(Dbm), 1); h = sl / 2
    return (h, math.exp(ic) / S6.C_H(h)) if h > 0 else (h, float('nan'))


def localc_curve(Dfull):
    lags = np.arange(len(Dfull)) * DT
    return np.array([t * quad(lambda u: math.exp(-u) * (u - 1) * np.interp(t * u, lags, Dfull), 0, 60, limit=400)[0]
                     for t in TF])


def jump_level(A):
    """J arm: the noise intensity at which the global minimiser of A + s2 beta_disc switches branch (bisection)."""
    br = lambda s2: local_tau(A, s2) > 1.0
    lo, hi = 1e-9, 4e-7
    if br(lo) or not br(hi):
        return float('nan')
    for _ in range(100):
        mid = math.sqrt(lo * hi)
        if br(mid): hi = mid
        else: lo = mid
    return math.sqrt(lo * hi)


def levels(arm, s2c=None):
    if ARMS[arm]['kind'] == 'shape':
        return np.array([S6.S0(v) for v in V_S])
    if arm == 'J':
        return s2c * J_REL
    return Q * 10.0 ** LADDER_K[arm]


def secants(tau, lev):
    tau = np.asarray(tau, float); L = np.log(np.asarray(lev, float))
    return np.diff(np.log(tau)) / np.diff(L)


def updown(tau, lev):
    """shape arms (levels V = 1, 256, 1/16): the exponents 1 -> 256 (up) and 1 -> 1/16 (down), as in step 6e."""
    tau = np.asarray(tau, float); lev = np.asarray(lev, float)
    return np.array([math.log(tau[1] / tau[0]) / math.log(lev[1] / lev[0]), math.log(tau[2] / tau[0]) / math.log(lev[2] / lev[0])])


def predict(arm, A, Dm, D0, Dp, Dsq, Dbm, NV=None, s2c=None):
    """every predictor at the arm's levels, from pooled calibration statistics (nothing from the evaluation)."""
    Dl = make_Dl(Dm, D0, Dp)
    if arm == 'J' and s2c is None:
        s2c = jump_level(A)
    lev = levels(arm, s2c); P = {'levels': lev}
    P['LOCAL'] = np.array([local_tau(A, s2) for s2 in lev])                    # white noise (shape-blind in S arms)
    P['LOCAL_G'] = np.array([grid_read(A, s2) for s2 in lev])
    if ARMS[arm]['kind'] == 'shape':
        P['AWARE'] = np.array([local_tau(A, Nf=np.asarray(NV[i])) for i in range(len(V_S))])
        P['AWARE_G'] = np.array([grid_read(A, Nf=np.asarray(NV[i])) for i in range(len(V_S))])
    Rf = rawf_curve(Dl); P['RAWF'] = np.array([first_crossing(Rf, s2) for s2 in lev])
    h, lam = engine_fit(Dbm); P['H_engine'] = h
    P['ENGINE'] = np.array([S6.kappa_H(h) * (s2 / lam) ** (1 / (2 * h + 1)) for s2 in lev])
    lsq = float(np.mean(np.asarray(Dsq) / L_SQ)); P['SQRT'] = np.sqrt(lev / lsq)
    rs = []
    for i in range(len(lev) - 1):
        ta, tb = P['RAWF'][i], P['RAWF'][i + 1]
        rs.append(1.0 / (2 * hd_at(Dl, math.sqrt(ta * tb)) + 1) if np.isfinite(ta) and np.isfinite(tb) else float('nan'))
    P['RAWS_exp'] = np.array(rs)
    if arm == 'J':
        P['s2c'] = s2c
    return P


HORIZON_PREDICTORS = ('LOCAL', 'LOCAL_G', 'AWARE', 'AWARE_G', 'RAWF', 'ENGINE', 'SQRT', 'LOCALC')


# ------------------------------------------------------------------------------------------------ stage 1: calibration
def cal_batch(args):
    """per calibration seed (a rest record of the drift): the tracking error A on TF from the structure function at all
    lags, D at TF / 1.1, TF, TF * 1.1 and at L_SQ, the block-mean structure function at L_BM, and the batch sum of D
    at all lags (pooled for LOCAL-C); shape arms, NCAL_SEEDS: the noise-only error curve at each V on TF."""
    arm, seeds = args; n = int(round(T_CAL / DT)); rows = []; Dsum = np.zeros(L_MAX + 1); fft_chk = 0.0
    for s in seeds:
        x = np.concatenate([[0.0], np.cumsum(drift_increments(arm, n, s, 0))])
        Df = struct_fn(x, L_MAX); Dsum += Df; lags = np.arange(L_MAX + 1) * DT
        Dl = lambda l: float(np.interp(l, lags, Df))
        if SMOKE:                                                   # FFT structure function against direct sums
            for j in (1, 10, 1000, L_MAX // 2):
                fft_chk = max(fft_chk, abs(Df[j] / np.mean((x[j:] - x[:-j]) ** 2) - 1))
        bm = []
        for l in L_BM:
            k = int(round(l / DT)); mb = x[: (n + 1) // k * k].reshape(-1, k).mean(1); bm.append(float(np.mean(np.diff(mb) ** 2)))
        o = dict(seed=s, A=a_disc(Df, TF), Dm=np.array([Dl(t / 1.1) for t in TF]), D0=np.array([Dl(t) for t in TF]),
                 Dp=np.array([Dl(t * 1.1) for t in TF]), Dsq=np.array([Dl(l) for l in L_SQ]), Dbm=np.array(bm))
        if ARMS[arm]['kind'] == 'shape' and s in NCAL_SEEDS:
            nn = int(round(T_NCAL / DT)); NV = []
            for nz in shape_noise(nn, s, arm, 3):
                row = []
                for t in TF:
                    _, en = filt_err(None, nz, t); sk = int(max(10 * t, 1.0) / DT); row.append(float(np.mean(en[sk:] ** 2)))
                NV.append(row)
            o['NV'] = np.array(NV)
        rows.append(o)
    return arm, rows, Dsum, fft_chk


def stage_calibrate():
    sha = sha_of(os.path.abspath(__file__))
    print('step6f_uprime_local.py  sha256', sha, '  stage: calibrate', flush=True)
    t0 = time.time(); tasks = [(arm, CAL_SEEDS[i:i + BATCH]) for arm in ARMS for i in range(0, len(CAL_SEEDS), BATCH)]
    parts = {arm: [] for arm in ARMS}; Dsum = {arm: np.zeros(L_MAX + 1) for arm in ARMS}; fchk = 0.0
    with ProcessPoolExecutor(max_workers=2) as ex:
        for arm, rows, ds, fc in ex.map(cal_batch, tasks):
            parts[arm] += rows; Dsum[arm] += ds; fchk = max(fchk, fc)
            print(f'  calibration {arm} seeds {rows[0]["seed"]}..{rows[-1]["seed"]} ({time.time() - t0:.0f} s)', flush=True)
    out = {'sha256': sha, 'smoke': SMOKE, 'TF': TF.tolist(), 'GRID': GRID.tolist(), 'cal_seeds': CAL_SEEDS,
           'B': B_BOOT}
    if SMOKE:
        out['fft_check_max_rel'] = fchk
    draws = {}; rb = np.random.default_rng(2025); keys = ('A', 'Dm', 'D0', 'Dp', 'Dsq', 'Dbm')
    for arm in ARMS:
        rows = sorted(parts[arm], key=lambda o: o['seed']); ns = len(rows)
        M = {k: np.array([o[k] for o in rows]) for k in keys}
        NVs = np.array([o['NV'] for o in rows if 'NV' in o]) if ARMS[arm]['kind'] == 'shape' else None
        P = predict(arm, *(M[k].mean(0) for k in keys), NV=None if NVs is None else NVs.mean(0))
        assert arm != 'J' or np.isfinite(P['s2c']), 'J: the jump level was not bracketed'
        assert SMOKE or np.all(np.isfinite(P['LOCAL_G'])), f'{arm}: a predicted window lies at a grid edge'
        P['LOCALC'] = np.array([first_crossing(localc_curve(Dsum[arm] / ns), s2) for s2 in P['levels']]) \
            if ARMS[arm]['kind'] == 'scaled' else np.full(len(P['levels']), np.nan)
        rec = {k: (v.tolist() if isinstance(v, np.ndarray) else v) for k, v in P.items()}
        for k in HORIZON_PREDICTORS:
            if k in P:
                rec[k + '_exp'] = secants(P[k], P['levels']).tolist() if ARMS[arm]['kind'] == 'scaled' else \
                    updown(P[k], P['levels']).tolist()
        rec['A_pooled'] = M['A'].mean(0).tolist()
        if NVs is not None:
            rec['NV_pooled'] = NVs.mean(0).tolist()
        out[arm] = rec
        Bd = {k: [] for k in ('LOCAL', 'LOCAL_G', 'AWARE', 'AWARE_G', 'RAWF', 'ENGINE', 'SQRT', 'RAWS_exp')}
        for b in range(B_BOOT):
            idx = rb.integers(0, ns, ns)
            nv = None if NVs is None else NVs[rb.integers(0, len(NVs), len(NVs))].mean(0)
            Pb = predict(arm, *(M[k][idx].mean(0) for k in keys), NV=nv, s2c=P.get('s2c'))
            for k in Bd:
                if k in Pb:
                    Bd[k].append(Pb[k])
        draws[arm] = {k: np.array(v, dtype=float) for k, v in Bd.items() if len(v)}
        print(f'  predictions {arm} frozen ({time.time() - t0:.0f} s)', flush=True)
    out['runtime_s'] = time.time() - t0
    json.dump(out, open(os.path.join(RES, f'step6f_predictions{TAG}.json'), 'w'), indent=1)
    np.savez_compressed(os.path.join(RES, f'step6f_pred_draws{TAG}.npz'),
                        **{f'{arm}__{k}': v for arm in draws for k, v in draws[arm].items()})
    print(f'written results/step6f_predictions{TAG}.json and results/step6f_pred_draws{TAG}.npz', flush=True)


# ------------------------------------------------------------------------------------------------ stages 2 and 3
def identity_check(arm, seed):
    """the error of step6_quad_attitude_estimator.Mahony (no bias integral), attitude truth included, against the sum
    of the two filtered parts used here; 20 s, three windows; returns the largest absolute difference (rad)."""
    n = int(round(20.0 / DT)); r = rng(seed, arm, 4, 8)
    rho = math.exp(-DT / S6.TRUTH_TC); c = S6.TRUTH_SD * math.sqrt(1 - rho ** 2)
    th = lfilter([c], [1, -rho], r.standard_normal(n + 1))
    dx = drift_increments(arm, n, seed, 4)
    nz = (unit_noise(n, seed, arm, 4) * math.sqrt(10.0)) if ARMS[arm]['kind'] == 'scaled' else shape_noise(n, seed, arm, 4)[0]
    W = np.diff(th) / DT + dx / DT; Acc = th[1:] + nz; worst = 0.0
    for tau in (0.05, 0.7, 9.0):
        m = S6.Mahony(np.array([th[0]]), withI=False); e = np.empty(n)
        for k in range(n):
            e[k] = m.step(np.array([W[k]]), np.array([Acc[k]]), np.array([tau]))[0] - th[k + 1]
        ed, en = filt_err(dx, nz, tau); worst = max(worst, float(np.max(np.abs(e - ed - en))))
    return worst


def ev_batch(args):
    """per evaluation seed and window: the steady mean squares of the drift part A, of the unit noise part B (one row
    per V in the shape arms) and of their product X. The error at level c is e_d + sqrt(c) e_n, so the steady
    mean-square error is A + c B + 2 sqrt(c) X, exactly; steady = after max(10 tau, 1 s)."""
    arm, seeds = args; n = int(round(T_EV / DT)); ng = len(GRID); out = {'seeds': seeds, 'A': [], 'B': [], 'X': []}
    for s in seeds:
        dx = drift_increments(arm, n, s, 1)
        nzs = [unit_noise(n, s, arm, 2)] if ARMS[arm]['kind'] == 'scaled' else shape_noise(n, s, arm, 2)
        A = np.empty(ng); B = np.empty((len(nzs), ng)); X = np.empty((len(nzs), ng))
        for g, t in enumerate(GRID):
            sk = int(max(10 * t, 1.0) / DT); ed, _ = filt_err(dx, None, t); ed = ed[sk:]; A[g] = ed @ ed / len(ed)
            for v, nz in enumerate(nzs):
                _, en = filt_err(None, nz, t); en = en[sk:]; B[v, g] = en @ en / len(en); X[v, g] = ed @ en / len(en)
        out['A'].append(A); out['B'].append(B); out['X'].append(X)
    out['ident'] = identity_check(arm, seeds[0]) if seeds[0] in (PILOT_SEEDS[0], EV_SEED0) else None
    return arm, out


def run_eval(seeds):
    tasks = [(arm, seeds[i:i + BATCH]) for arm in ARMS for i in range(0, len(seeds), BATCH)]
    parts = {arm: [] for arm in ARMS}; t0 = time.time()
    with ProcessPoolExecutor(max_workers=2) as ex:
        for arm, o in ex.map(ev_batch, tasks):
            parts[arm].append(o); print(f'  {arm} seeds {o["seeds"][0]}..{o["seeds"][-1]} ({time.time() - t0:.0f} s)', flush=True)
    rows = {}
    for arm in ARMS:
        ps = sorted(parts[arm], key=lambda o: o['seeds'][0])
        rows[arm] = {'seeds': np.array(sum((o['seeds'] for o in ps), [])),
                     'A': np.concatenate([np.array(o['A']) for o in ps]), 'B': np.concatenate([np.array(o['B']) for o in ps]),
                     'X': np.concatenate([np.array(o['X']) for o in ps]),
                     'ident': max([o['ident'] for o in ps if o['ident'] is not None], default=float('nan'))}
    return rows, time.time() - t0


def stage_pilot():
    """100 pilot seeds through the evaluation machinery. Prints and writes only standard errors, grid-edge flags, the
    identity check and the number of evaluation seeds N given by the pre-registered rule; never a window, an exponent,
    a difference or a verdict."""
    import evaluate_step6f as EV
    sha = sha_of(os.path.abspath(__file__))
    print('step6f_uprime_local.py  sha256', sha, '  stage: pilot', flush=True)
    rows, rt = run_eval(PILOT_SEEDS)
    txt = EV.pilot_report(rows, sha, rt)
    print(txt); open(os.path.join(RES, f'step6f_pilot{TAG}.txt'), 'w').write(txt + '\n')


def stage_evaluate():
    sha = sha_of(os.path.abspath(__file__))
    print('step6f_uprime_local.py  sha256', sha, '  stage: evaluate', flush=True)
    if SMOKE:
        N = len(PILOT_SEEDS)
    else:
        line = [l for l in open(os.path.join(RES, 'step6f_pilot.txt')) if l.startswith('N_EVAL =')][0]
        N = int(line.split('=')[1].split()[0])
    rows, rt = run_eval(list(range(EV_SEED0, EV_SEED0 + N)))
    np.savez_compressed(os.path.join(RES, f'step6f_rows{TAG}.npz'),
                        **{f'{arm}__{k}': np.asarray(v) for arm in rows for k, v in rows[arm].items()},
                        meta=np.array(json.dumps({'sha256': sha, 'smoke': SMOKE, 'N': N, 'runtime_s': rt})))
    print(f'written results/step6f_rows{TAG}.npz ({N} seeds per arm, {rt:.0f} s)', flush=True)


if __name__ == '__main__':
    stage = sys.argv[1] if len(sys.argv) > 1 else ''
    {'calibrate': stage_calibrate, 'pilot': stage_pilot, 'evaluate': stage_evaluate}[stage]()
