"""STEP 0 of the software bench: does the published syntonic engine land on the swept optimum for H != 1/2?

PRE-REGISTRATION (written before any run; the SHA-256 of this file is printed at the top of the results).

Object under test
  The identification arithmetic of jpbronsard/syntonic-engine (python/syntonic_engine/engine.py, v0.1.1):
  noise from two independent reads, structure function of the stream at lags 2..64 ms, H from its log-log slope,
  lambda_H from its intercept, window tau* = kappa_H (s2 / lambda_H)^(1/(2H+1)) (Theorem U'), floor dt/0.45,
  ceiling tau_max, rate limit tau_rate. Here it is applied as a SMOOTHER (an exponential window on the pair-mean of
  the two reads), the natural estimation-side counterpart of the joint correction loop; the identification lines
  are copied from engine.py without change.

Targets
  Fractional Brownian motion of Hurst index H in {0.30, 0.50, 0.70, 0.85, 0.95}, structure function
  D(s) = lambda_H s^(2H) with lambda_H = 1 (rad^2 / s^2H), sampled at dt = 1 ms, 2^16 samples (65.5 s),
  Davies-Harte exact synthesis. Two independent white Gaussian reads per tick, read noise chosen so that the U'
  optimum tau_U' = kappa_H (s2/lambda_H)^(1/(2H+1)) sits at 7 log-spaced targets between 5 ms and 160 ms
  (1.5 decades of tau, i.e. 1.5(2H+1) decades of s2/lambda_H: 3 decades at H = 1/2). s2 = noise PSD of the pair
  mean = sigma_read^2 dt / 2, exactly as engine.py defines it.

Measurements, per (H, level, seed)
  tau_opt : the window of the fixed exponential smoother with the smallest MSE against the true x, from a sweep
            of 40 windows between 2.5 ms and 1 s, refined by a parabola in log tau on the three best points.
  tau_eng : the median window of the syntonic smoother after warm-up (t > 5 s).
  regret  : MSE_eng / MSE_opt - 1, both MSEs on t > 5 s.
  H_eng   : the median H the engine estimated after warm-up.

Hypotheses (fixed now)
  H1 (Theorem U' in this class): for each H, the slope of log tau_opt against log(s2/lambda_H) over the 7 levels
      lies within +/- 0.05 of alpha = 1/(2H+1). (alpha = 0.625, 0.500, 0.417, 0.370, 0.345.)
  H2 (window landing): for every (H, level), the geometric mean over seeds of tau_eng / tau_opt lies in [1/1.5, 1.5].
  H3 (regret): for every (H, level), the mean regret over seeds is at most 0.10.
  H4 (H estimate): for every (H, level), the mean H_eng over seeds is within +/- 0.10 of the true H.

Failure reading
  H1 fails  -> the U' exponent does not describe an exponential smoother on fBm (would contradict the theorem;
               check the synthesis first).
  H2/H3 fail at some H -> the engine's identification is biased there; fix the engine before any physical run.
  H4 fails but H3 holds -> the H estimate is biased but the window is insensitive to it in that region.

Seeds
  Seeds 0..9 are used here. Seeds >= 100 are reserved, untouched, for the main bench.
  Within one (H, seed), the SAME fBm path is used at all 7 noise levels and only the read noise changes, so the
  H1 slope isolates the s2 scaling from the path-to-path amplitude scatter of long-memory fBm (at H = 0.95 the
  single-path structure function varies by about +/- 40 % across paths; checked on 30 paths before this run,
  mean ratio to theory 1.02).

Synthesis check (done before this run, not a result): Davies-Harte circulant eigenvalues are all positive for
  H <= 0.95 at n = 2^16..2^20; the mean sample structure function over 30 paths matches s^(2H) within 2 %.

Nothing in this script is tuned after seeing results. One pass.
"""
import hashlib, json, os, sys, time
from math import gamma as _gamma
import numpy as np
from scipy.signal import lfilter

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'step0')
os.makedirs(OUT, exist_ok=True)

# ----------------------------------------------------------------------------------------------- fixed settings
DT = 1e-3
N = 2 ** 16
H_LIST = [0.30, 0.50, 0.70, 0.85, 0.95]
LAMBDA_H = 1.0
TAU_TARGETS = np.geomspace(5e-3, 160e-3, 7)          # where the U' optimum is placed, per level
SEEDS = list(range(10))
SWEEP = np.geomspace(2.5e-3, 1.0, 40)                # fixed windows swept
T_SKIP = 5.0                                         # seconds discarded (engine warm-up and EMA settling)
ENGINE = dict(floor=DT / 0.45, tau_max=1.0, tau_rate=0.4, lags=(0.002, 0.004, 0.008, 0.016, 0.032, 0.064),
              snr_gate=2.0, ema=1.0, warmup=1.0)     # engine.py defaults, except tau_max raised from 0.5 to 1.0
                                                     # so the ceiling never masks a 160 ms optimum, and no servo
                                                     # floor (there is no servo: floor = dt / 0.45 of TH41-A)


def kappa_H(H):
    """Copied from engine.py."""
    H = np.asarray(H, dtype=float)
    return (2 * H * np.vectorize(_gamma)(2 * H + 1)) ** (-1.0 / (2 * H + 1))


def fbm(n, H, rng, dt):
    """Exact fractional Brownian motion on n points (Davies-Harte), B(0) = 0, D(s) = s^(2H) in seconds."""
    k = np.arange(n + 1, dtype=float)
    g = 0.5 * (np.abs(k + 1) ** (2 * H) - 2 * np.abs(k) ** (2 * H) + np.abs(k - 1) ** (2 * H))
    r = np.concatenate([g, g[-2:0:-1]])              # circulant row, length 2n
    lam = np.maximum(np.fft.fft(r).real, 0.0)
    M = r.size
    W = rng.standard_normal(M) + 1j * rng.standard_normal(M)
    Y = np.fft.ifft(np.sqrt(lam) * W) * np.sqrt(M)
    fgn = Y.real[:n]
    return dt ** H * np.concatenate([[0.0], np.cumsum(fgn)[:-1]])


class SyntonicSmoother:
    """engine.py's identification, vectorized over n streams, driving an exponential window on the pair mean."""

    def __init__(self, n, dt, floor, tau_max, tau_rate, lags, snr_gate, ema, warmup):
        self.n, self.dt = n, dt
        self.floor, self.tau_max, self.tau_rate, self.snr_gate = floor, tau_max, tau_rate, snr_gate
        self.lags = np.unique(np.maximum(1, np.round(np.asarray(lags) / dt))).astype(int)
        self.logl = np.log(self.lags * dt)
        self.a_w = dt / ema
        self.warm = int(round(warmup / dt))
        self.k = 0
        self.tau = np.full(n, 0.5); self.H = np.full(n, 0.5); self.tau_target = np.full(n, 0.5)
        self.n2 = np.zeros(n)
        self.r_hist = np.zeros((self.lags.max() + 1, n)); self.Dk = np.zeros((len(self.lags), n))
        self.xhat = None

    def step(self, a, b):
        dt = self.dt
        meas = 0.5 * (a + b)
        d = (a - b) / np.sqrt(2.0)
        self.n2 += self.a_w * (d * d - self.n2)                              # per-read noise variance   [engine.py]
        r = meas                                                             # the stream itself (no plan, no servo)
        self.r_hist = np.roll(self.r_hist, 1, axis=0); self.r_hist[0] = r
        if self.k >= self.lags.max():
            for i, lg in enumerate(self.lags):
                dd = r - self.r_hist[lg]; self.Dk[i] += self.a_w * (dd * dd - self.Dk[i])
        if self.k > self.warm:                                               # ---- identification, verbatim engine.py
            D = self.Dk - self.n2[None, :]
            ok_g = D > self.snr_gate * self.n2[None, :]
            ok_u = D > 0
            use_g = ok_g.sum(0) >= 2
            ok = np.where(use_g[None], ok_g, ok_u)
            nok = ok.sum(0)
            logD = np.where(ok, np.log(np.maximum(D, 1e-30)), 0.0)
            W = ok.astype(float); L = self.logl[:, None]
            Sw = W.sum(0); Sx = (W * L).sum(0); Sy = (W * logD).sum(0)
            Sxx = (W * L ** 2).sum(0); Sxy = (W * L * logD).sum(0)
            den = Sw * Sxx - Sx ** 2
            slope = np.where(den > 1e-12, (Sw * Sxy - Sx * Sy) / np.maximum(den, 1e-12), 1.0)
            Hn = np.clip(0.5 * slope, 0.05, 1.0)
            valid = nok >= 2
            s2 = 0.5 * self.n2 * dt
            lam = np.exp((W * (logD - 2 * Hn[None] * L)).sum(0) / np.maximum(Sw, 1))
            ts = kappa_H(Hn) * (np.maximum(s2, 1e-30) / np.maximum(lam, 1e-30)) ** (1.0 / (2 * Hn + 1))
            ts = np.clip(ts, self.floor, self.tau_max)
            self.tau_target = np.where(valid, ts, self.tau_target); self.H = np.where(valid, Hn, self.H)
            rr = 1.0 / self.tau; dr_max = self.tau_rate * dt / self.tau ** 2
            rr = rr + np.clip(1.0 / self.tau_target - rr, -dr_max, dr_max); self.tau = 1.0 / rr
        g = 1.0 - np.exp(-dt / self.tau)                                     # ---- the window, applied as a smoother
        self.xhat = meas.copy() if self.xhat is None else self.xhat + g * (meas - self.xhat)
        self.k += 1
        return self.xhat


def sweep_opt(meas, x, keep):
    """Best fixed exponential window on this stream: coarse sweep, parabolic refinement in log tau."""
    mse = np.empty(SWEEP.size)
    for i, tau in enumerate(SWEEP):
        g = 1.0 - np.exp(-DT / tau)
        xh, _ = lfilter([g], [1.0, -(1.0 - g)], meas, zi=[(1.0 - g) * x[0]])
        mse[i] = np.mean((xh[keep] - x[keep]) ** 2)
    j = int(np.argmin(mse))
    if 0 < j < SWEEP.size - 1:
        l = np.log(SWEEP[j - 1:j + 2]); y = mse[j - 1:j + 2]
        A = np.vstack([l ** 2, l, np.ones(3)]).T; c = np.linalg.solve(A, y)
        lt = -c[1] / (2 * c[0]) if c[0] > 0 else l[1]
        lt = float(np.clip(lt, l[0], l[2]))
        return float(np.exp(lt)), float(c[0] * lt ** 2 + c[1] * lt + c[2]) if c[0] > 0 else float(mse[j]), j
    return float(SWEEP[j]), float(mse[j]), j


def main():
    src = open(os.path.abspath(__file__), 'rb').read()
    sha = hashlib.sha256(src).hexdigest()
    print('step0_uprime_engine.py  sha256', sha, flush=True)
    t = np.arange(N) * DT; keep = t > T_SKIP
    rows = []
    t0 = time.time()
    for H in H_LIST:
        kH = float(kappa_H(H)); alpha = 1.0 / (2 * H + 1)
        # read noise per level: tau_U' = kH (s2/lambda)^alpha  ->  s2 = lambda (tau/kH)^(1/alpha);  s2 = sig^2 dt / 2
        s2_levels = LAMBDA_H * (TAU_TARGETS / kH) ** (2 * H + 1)
        sig_levels = np.sqrt(2 * s2_levels / DT)
        # one vectorized engine run per H: n = levels x seeds streams
        X, A, B, cells = [], [], [], []
        for li, sig in enumerate(sig_levels):
            for seed in SEEDS:
                rng_path = np.random.default_rng(10_000 * seed + int(round(H * 100)))        # same path at every level
                rng = np.random.default_rng(10_000 * seed + 100 * (li + 1) + int(round(H * 100)))  # fresh read noise
                x = np.sqrt(LAMBDA_H) * fbm(N, H, rng_path, DT)
                a = x + sig * rng.standard_normal(N); b = x + sig * rng.standard_normal(N)
                X.append(x); A.append(a); B.append(b); cells.append((li, seed))
        X = np.array(X).T; A = np.array(A).T; B = np.array(B).T           # (N, n)
        n = X.shape[1]
        eng = SyntonicSmoother(n, DT, **ENGINE)
        XH = np.empty_like(X); TAU = np.empty_like(X); HE = np.empty_like(X)
        for k in range(N):
            XH[k] = eng.step(A[k], B[k]); TAU[k] = eng.tau; HE[k] = eng.H
        mse_eng = np.mean((XH[keep] - X[keep]) ** 2, axis=0)
        tau_eng = np.median(TAU[keep], axis=0); h_eng = np.median(HE[keep], axis=0)
        for c, (li, seed) in enumerate(cells):
            meas = 0.5 * (A[:, c] + B[:, c])
            tau_opt, mse_opt, j = sweep_opt(meas, X[:, c], keep)
            rows.append(dict(H=H, level=li, seed=seed, s2=float(s2_levels[li]), tau_uprime=float(TAU_TARGETS[li]),
                             tau_opt=tau_opt, mse_opt=mse_opt, sweep_idx=j, tau_eng=float(tau_eng[c]),
                             mse_eng=float(mse_eng[c]), H_eng=float(h_eng[c]),
                             regret=float(mse_eng[c] / mse_opt - 1.0)))
        print(f'H={H:.2f} done, {time.time() - t0:.0f} s', flush=True)
    json.dump(dict(sha256=sha, rows=rows, settings=dict(DT=DT, N=N, H_LIST=H_LIST, LAMBDA_H=LAMBDA_H,
              TAU_TARGETS=TAU_TARGETS.tolist(), SEEDS=SEEDS, SWEEP=SWEEP.tolist(), T_SKIP=T_SKIP, ENGINE=ENGINE)),
              open(os.path.join(OUT, 'step0_results.json'), 'w'), indent=1)

    # ---------------------------------------------------------------------------------------- evaluation of H1-H4
    print('\n=== STEP 0 results (sha256 %s) ===' % sha[:12])
    verdict = {}
    for H in H_LIST:
        R = [r for r in rows if r['H'] == H]
        alpha = 1.0 / (2 * H + 1)
        # H1: slope of log tau_opt vs log s2 (lambda fixed), per-seed fit then mean
        slopes = []
        for seed in SEEDS:
            rs = sorted([r for r in R if r['seed'] == seed], key=lambda r: r['level'])
            slopes.append(np.polyfit(np.log([r['s2'] for r in rs]), np.log([r['tau_opt'] for r in rs]), 1)[0])
        sl, sls = float(np.mean(slopes)), float(np.std(slopes, ddof=1) / np.sqrt(len(slopes)))
        h1 = abs(sl - alpha) <= 0.05
        print(f'\nH = {H:.2f}   alpha_U\' = {alpha:.3f}   measured slope = {sl:.3f} +/- {sls:.3f}   H1 {"PASS" if h1 else "FAIL"}')
        print('  level  tau_U\'    tau_opt    tau_eng/tau_opt   regret    H_eng      H2    H3    H4')
        h2 = h3 = h4 = True
        for li in range(len(TAU_TARGETS)):
            rl = [r for r in R if r['level'] == li]
            ratio = float(np.exp(np.mean([np.log(r['tau_eng'] / r['tau_opt']) for r in rl])))
            reg = float(np.mean([r['regret'] for r in rl])); he = float(np.mean([r['H_eng'] for r in rl]))
            topt = float(np.exp(np.mean([np.log(r['tau_opt']) for r in rl])))
            p2 = 1 / 1.5 <= ratio <= 1.5; p3 = reg <= 0.10; p4 = abs(he - H) <= 0.10
            h2 &= p2; h3 &= p3; h4 &= p4
            print(f'  {li:5d}  {TAU_TARGETS[li]*1e3:6.1f} ms  {topt*1e3:6.1f} ms   {ratio:6.2f}         {reg:6.3f}   {he:5.2f}     '
                  f'{"ok" if p2 else "NO":3s}   {"ok" if p3 else "NO":3s}   {"ok" if p4 else "NO":3s}')
        verdict[H] = dict(H1=h1, H2=h2, H3=h3, H4=h4, slope=sl, slope_se=sls)
    print('\nSummary:', json.dumps({f'{H:.2f}': v for H, v in verdict.items()}, indent=0))
    json.dump(verdict, open(os.path.join(OUT, 'step0_verdict.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
