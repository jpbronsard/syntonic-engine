"""STEP 0b, EXPLORATORY (not pre-registered): the same grid as step0, with the identification changed in one way.

Finding of step 0 (sha 7d1d419d...): Theorem U' holds (H1 passes at every H, slopes within 0.01 of 1/(2H+1)), but
the engine's identification fails as soon as the read noise is large enough that the signal's increments over
its longest structure-function lag (64 ms) fall below the per-read noise floor. The lag at which the signal
clears the floor grows from 0.1 s (level 0) to 100 s (level 6, H = 1/2), so no fixed raw-read lag set can work.

Change tested here: the structure function is taken on BLOCK MEANS (adjacent non-overlapping blocks of length l,
l = 2 ms .. 4.096 s in powers of two), the Allan-variance construction. For white read noise the floor of the
block-mean difference is exactly n2 dt / l (n2 = per-read variance measured from the two reads, as in engine.py),
so it is subtracted analytically; for fBm the signal term is lambda C(H) l^(2H) with
C(H) = (2^(2H+2) - 4) / ((2H+1)(2H+2))  (= 2/3 at H = 1/2, the Allan variance of a random walk).
Slope -> H; intercept / C(H) -> lambda_H; the window law and all the rest unchanged.

Everything else (targets, levels, seeds, sweep, criteria H2-H4) is identical to step 0 so the tables compare.
This is an engineering fix found after seeing step-0 results; it is reported as exploratory. The main bench will be
pre-registered on the fixed engine.
"""
import hashlib, json, os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import step0_uprime_engine as S0                      # targets, sweep, kappa_H, settings
OUT = os.path.join(HERE, 'step0b'); os.makedirs(OUT, exist_ok=True)

SCALES = [2 ** i for i in range(1, 13)]               # block lengths in samples: 2 .. 4096 (2 ms .. 4.096 s)


def C_H(H):
    return (2.0 ** (2 * H + 2) - 4.0) / ((2 * H + 1) * (2 * H + 2))


class AllanSmoother:
    """engine.py's law with a block-mean (Allan) structure function in place of the raw-read one."""

    def __init__(self, n, dt, floor, tau_max, tau_rate, snr_gate=2.0, ema=1.0, warmup=1.0, scales=SCALES):
        self.n, self.dt = n, dt
        self.floor, self.tau_max, self.tau_rate, self.snr_gate = floor, tau_max, tau_rate, snr_gate
        self.scales = np.asarray(scales); self.logl = np.log(self.scales * dt)
        self.a_w = dt / ema; self.warm = int(round(warmup / dt))
        m = len(scales)
        self.acc = np.zeros((m, n)); self.cnt = np.zeros(m, dtype=int)
        self.prev = np.full((m, n), np.nan); self.DA = np.zeros((m, n)); self.nblk = np.zeros(m, dtype=int)
        # per-scale averaging of the squared block differences: EMA over blocks with time constant max(ema, 8 l)
        self.a_blk = np.array([min(1.0, dt * s / max(ema, 8 * dt * s)) for s in scales])
        self.k = 0; self.n2 = np.zeros(n)
        self.tau = np.full(n, 0.5); self.H = np.full(n, 0.5); self.tau_target = np.full(n, 0.5)
        self.xhat = None

    def step(self, a, b):
        dt = self.dt
        meas = 0.5 * (a + b); d = (a - b) / np.sqrt(2.0)
        self.n2 += self.a_w * (d * d - self.n2)
        self.acc += meas[None, :]; self.cnt += 1
        done = self.cnt == self.scales
        if done.any():
            for i in np.flatnonzero(done):
                mean = self.acc[i] / self.scales[i]
                if not np.isnan(self.prev[i, 0]):
                    dd = mean - self.prev[i]
                    aw = self.a_blk[i] if self.nblk[i] > 0 else 1.0
                    self.DA[i] += aw * (dd * dd - self.DA[i]); self.nblk[i] += 1
                self.prev[i] = mean; self.acc[i] = 0.0; self.cnt[i] = 0
        if self.k > self.warm:
            floor = self.n2[None, :] * dt / (self.scales * dt)[:, None]           # n2 dt / l
            D = self.DA - floor
            have = (self.nblk >= 2)[:, None] & np.ones((1, self.n), bool)
            ok_g = have & (D > self.snr_gate * floor)
            ok_u = have & (D > 0)
            use_g = ok_g.sum(0) >= 2
            ok = np.where(use_g[None], ok_g, ok_u); nok = ok.sum(0)
            logD = np.where(ok, np.log(np.maximum(D, 1e-30)), 0.0)
            W = ok.astype(float); L = self.logl[:, None]
            Sw = W.sum(0); Sx = (W * L).sum(0); Sy = (W * logD).sum(0)
            Sxx = (W * L ** 2).sum(0); Sxy = (W * L * logD).sum(0)
            den = Sw * Sxx - Sx ** 2
            slope = np.where(den > 1e-12, (Sw * Sxy - Sx * Sy) / np.maximum(den, 1e-12), 1.0)
            Hn = np.clip(0.5 * slope, 0.05, 1.0)
            valid = nok >= 2
            s2 = 0.5 * self.n2 * dt
            lam = np.exp((W * (logD - 2 * Hn[None] * L)).sum(0) / np.maximum(Sw, 1)) / C_H(Hn)
            ts = S0.kappa_H(Hn) * (np.maximum(s2, 1e-30) / np.maximum(lam, 1e-30)) ** (1.0 / (2 * Hn + 1))
            ts = np.clip(ts, self.floor, self.tau_max)
            self.tau_target = np.where(valid, ts, self.tau_target); self.H = np.where(valid, Hn, self.H)
            rr = 1.0 / self.tau; dr_max = self.tau_rate * dt / self.tau ** 2
            rr = rr + np.clip(1.0 / self.tau_target - rr, -dr_max, dr_max); self.tau = 1.0 / rr
        g = 1.0 - np.exp(-dt / self.tau)
        self.xhat = meas.copy() if self.xhat is None else self.xhat + g * (meas - self.xhat)
        self.k += 1
        return self.xhat


def main():
    sha = hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest()
    print('step0b_allan_identification.py  sha256', sha, flush=True)
    DT, N = S0.DT, S0.N; t = np.arange(N) * DT; keep = t > S0.T_SKIP
    rows = []; t0 = time.time()
    for H in S0.H_LIST:
        kH = float(S0.kappa_H(H))
        s2_levels = S0.LAMBDA_H * (S0.TAU_TARGETS / kH) ** (2 * H + 1); sig_levels = np.sqrt(2 * s2_levels / DT)
        X, A, B, cells = [], [], [], []
        for li, sig in enumerate(sig_levels):
            for seed in S0.SEEDS:
                rng_path = np.random.default_rng(10_000 * seed + int(round(H * 100)))
                rng = np.random.default_rng(10_000 * seed + 100 * (li + 1) + int(round(H * 100)))
                x = np.sqrt(S0.LAMBDA_H) * S0.fbm(N, H, rng_path, DT)
                a = x + sig * rng.standard_normal(N); b = x + sig * rng.standard_normal(N)
                X.append(x); A.append(a); B.append(b); cells.append((li, seed))
        X = np.array(X).T; A = np.array(A).T; B = np.array(B).T; n = X.shape[1]
        eng = AllanSmoother(n, DT, S0.ENGINE['floor'], S0.ENGINE['tau_max'], S0.ENGINE['tau_rate'])
        XH = np.empty_like(X); TAU = np.empty_like(X); HE = np.empty_like(X)
        for k in range(N):
            XH[k] = eng.step(A[k], B[k]); TAU[k] = eng.tau; HE[k] = eng.H
        mse_eng = np.mean((XH[keep] - X[keep]) ** 2, axis=0)
        tau_eng = np.median(TAU[keep], axis=0); h_eng = np.median(HE[keep], axis=0)
        for c, (li, seed) in enumerate(cells):
            meas = 0.5 * (A[:, c] + B[:, c])
            tau_opt, mse_opt, j = S0.sweep_opt(meas, X[:, c], keep)
            rows.append(dict(H=H, level=li, seed=seed, tau_opt=tau_opt, mse_opt=mse_opt, tau_eng=float(tau_eng[c]),
                             mse_eng=float(mse_eng[c]), H_eng=float(h_eng[c]), regret=float(mse_eng[c] / mse_opt - 1)))
        print(f'H={H:.2f} done, {time.time() - t0:.0f} s', flush=True)
    json.dump(dict(sha256=sha, rows=rows), open(os.path.join(OUT, 'step0b_results.json'), 'w'), indent=1)
    print('\n=== STEP 0b (exploratory, Allan identification) ===')
    for H in S0.H_LIST:
        R = [r for r in rows if r['H'] == H]
        print(f'\nH = {H:.2f}')
        print('  level  tau_U\'    tau_opt    tau_eng/tau_opt   regret    H_eng      H2    H3    H4')
        for li in range(len(S0.TAU_TARGETS)):
            rl = [r for r in R if r['level'] == li]
            ratio = float(np.exp(np.mean([np.log(r['tau_eng'] / r['tau_opt']) for r in rl])))
            reg = float(np.mean([r['regret'] for r in rl])); he = float(np.mean([r['H_eng'] for r in rl]))
            topt = float(np.exp(np.mean([np.log(r['tau_opt']) for r in rl])))
            p2 = 1 / 1.5 <= ratio <= 1.5; p3 = reg <= 0.10; p4 = abs(he - H) <= 0.10
            print(f'  {li:5d}  {S0.TAU_TARGETS[li]*1e3:6.1f} ms  {topt*1e3:6.1f} ms   {ratio:6.2f}         {reg:6.3f}   {he:5.2f}     '
                  f'{"ok" if p2 else "NO":3s}   {"ok" if p3 else "NO":3s}   {"ok" if p4 else "NO":3s}')


if __name__ == '__main__':
    main()
