"""STEP 1 of the software bench: the syntonic engine (Allan identification, step 0b) against adaptive Kalman filters.

PRE-REGISTRATION (written before the run; SHA-256 printed at the top of the results).

The claim under test (JP, 26 Oct 2026, evening): an adaptive Kalman filter adapts the two intensities (process
noise q, read noise r) but keeps its model class; the syntonic engine also reads the memory class H from the stream.
Hence: parity with the right-model Kalman at H = 1/2, better than a wrong-model Kalman elsewhere.

Competitors, all fed the SAME two reads per tick as the engine, nothing fitted to the truth:
  KF-RW  : bank of scalar Kalman filters with the random-walk model x_{k+1} = x_k + w, w ~ N(0, q), q on a log grid
           (q/r from 1e-8 to 1e2, 6 points per decade), r from the two reads (n2/2, the pair-mean variance, exactly
           the engine's statistic). Bayesian multiple-model adaptive estimation (Magill 1965): posterior over the
           grid from the cumulative innovation likelihood, estimate = posterior mixture. No forgetting in this
           stationary step, so this is the maximum-likelihood adaptive Kalman of its class: the strongest fair
           version. (Amendment before the run: innovation-covariance matching with a 1 s window, first drafted,
           was found noise-dominated in a sanity check at H = 1/2 (q underestimated by 37 %, +44 % MSE); it was
           replaced by the bank. No result of this step existed when the change was made.)
  KF-CV  : the same bank with the constant-velocity (white-acceleration) model, Q = q [[dt^3/3, dt^2/2],[dt^2/2, dt]],
           q dt^3 / r on a log grid from 1e-12 to 1e4, 6 points per decade.
  MM-UNION : one posterior over the union of both banks (RW grid + CV grid): the Bayesian choice of model class
           AND intensity. For stationary targets this is what an IMM converges to; the switching IMM is kept for
           the regime-change step. (Named IMM in the hypotheses below for continuity; it is this union.)
  ENGINE : the syntonic smoother with Allan identification (step 0b), zero configuration, 1 s averaging.
  Yardstick: the best fixed exponential window on the same stream (sweep of step 0), MSE_opt. Every result is
  reported as the ratio MSE / MSE_opt - 1 ("regret"; negative means better than any fixed exponential window).

Targets, levels, seeds: identical to step 0 (fBm H in {0.30, 0.50, 0.70, 0.85, 0.95}; 7 noise levels placing the
U' optimum from 5 ms to 160 ms; seeds 0..9; same path across levels within a seed; 65.5 s at 1 kHz; t > 5 s scored).

Hypotheses (fixed now)
  C1 (right model, parity): at H = 0.50, KF-RW regret <= 0.05 and ENGINE regret <= 0.10 at every level.
  C2 (wrong model costs): at H in {0.85, 0.95}, KF-RW regret >= 0.20 at every level with tau_opt >= 16 ms;
      at H in {0.30, 0.50}, KF-CV regret >= 0.20 at every level with tau_opt >= 16 ms.
  C3 (engine vs single-model adaptive Kalman): at every (H, level), ENGINE MSE <= 1.10 x min(KF-RW, KF-CV) MSE.
      Named risk, stated before the run: at H >= 0.85 the constant-velocity filter can extrapolate a locally
      linear path, which no first-order exponential window can. If KF-CV beats the engine there by more than
      10 %, the limitation is the engine's estimator class (one first-order window), not its identification;
      the law would then need a second-order (vector) window to compete. C3 may fail for that reason.
  C4 (IMM): at every (H, level), IMM MSE <= 1.10 x min(KF-RW, KF-CV) MSE (it should find the better model).
  Q5 (open question, no prediction): ENGINE against IMM at the intermediate H = 0.70, where neither IMM model
      is right.

Failure reading
  C1 fails -> a bug in the Kalman implementation or in the engine; check before anything else.
  C2 fails -> model mis-specification costs less than assumed in this noise range; the "wrong model" argument
      is weaker than claimed.
  C3 fails only at high H -> the estimator-class limitation above; report it as such.
  C3 fails at low or mid H -> the engine's identification is beaten by covariance matching; the claim is wrong.

One pass. Seeds 0..9. Nothing changed after seeing results.

Run note: the first launch (same code, sha 2a1c32bd...) was killed by the session after H = 0.85, before any
table was printed or saved; no result of it was seen. Per-H checkpointing was added (I/O only, no change of
method) and the run restarted from scratch.
"""
import hashlib, json, os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import step0_uprime_engine as S0
from step0b_allan_identification import AllanSmoother
OUT = os.path.join(HERE, 'step1'); os.makedirs(OUT, exist_ok=True)

DT, N = S0.DT, S0.N
EMA = 1.0                                   # averaging time of every adaptive statistic (engine's default)


class BankKF:
    """n streams x m Kalman filters of one model class, q on a log grid; returns per-filter estimates and
    per-filter log-likelihood increments. kind = 'rw' (scalar) or 'cv' (constant velocity, 2 states)."""

    def __init__(self, n, kind, x0, grid):
        self.n, self.kind, self.m = n, kind, grid.size
        self.grid = grid                                      # q/r (rw) or q dt^3/r (cv), dimensionless
        self.a_w = DT / EMA; self.n2 = np.zeros(n); self.k = 0
        if kind == 'rw':
            self.x = np.tile(x0[:, None], (1, self.m)); self.P = np.ones((n, self.m))
        else:
            self.x = np.zeros((n, self.m, 2)); self.x[:, :, 0] = x0[:, None]
            self.P = np.tile(np.diag([1.0, 1.0]), (n, self.m, 1, 1))
            self.F = np.array([[1.0, DT], [0.0, 1.0]]); self.Qs = np.array([[DT ** 3 / 3, DT ** 2 / 2], [DT ** 2 / 2, DT]])

    def step(self, meas, r):
        """r: (n,) read variance of the pair mean. Returns est (n,m), loglik increment (n,m)."""
        if self.kind == 'rw':
            q = r[:, None] * self.grid[None, :]
            Pm = self.P + q; S = Pm + r[:, None]; innov = meas[:, None] - self.x; K = Pm / S
            self.x = self.x + K * innov; self.P = (1 - K) * Pm
            est = self.x
        else:
            q = r[:, None] * self.grid[None, :] / DT ** 3
            xm = self.x @ self.F.T
            Pm = np.einsum('ij,nmjk,lk->nmil', self.F, self.P, self.F) + q[:, :, None, None] * self.Qs[None, None]
            S = Pm[:, :, 0, 0] + r[:, None]; innov = meas[:, None] - xm[:, :, 0]; K = Pm[:, :, :, 0] / S[:, :, None]
            self.x = xm + K * innov[:, :, None]
            IKH = np.tile(np.eye(2), (self.n, self.m, 1, 1)); IKH[:, :, :, 0] -= K
            self.P = np.einsum('nmij,nmjk->nmik', IKH, Pm)
            est = self.x[:, :, 0]
        ll = -0.5 * (innov * innov / S + np.log(S))
        return est, ll


class Competitors:
    """KF-RW bank, KF-CV bank and their union, sharing the read-noise estimate; posteriors by cumulative likelihood."""

    def __init__(self, n, x0):
        self.n = n; self.a_w = DT / EMA; self.n2 = np.zeros(n); self.k = 0
        self.rw = BankKF(n, 'rw', x0, np.geomspace(1e-8, 1e2, 61))
        self.cv = BankKF(n, 'cv', x0, np.geomspace(1e-12, 1e4, 97))
        self.L_rw = np.zeros((n, self.rw.m)); self.L_cv = np.zeros((n, self.cv.m))

    @staticmethod
    def _mix(est, L):
        w = np.exp(L - L.max(1, keepdims=True)); w /= w.sum(1, keepdims=True)
        return (w * est).sum(1)

    def step(self, a, b):
        meas = 0.5 * (a + b); d = (a - b) / np.sqrt(2.0)
        self.n2 += self.a_w * (d * d - self.n2); r = np.maximum(self.n2 / 2, 1e-30)
        e_rw, ll_rw = self.rw.step(meas, r); e_cv, ll_cv = self.cv.step(meas, r)
        if self.k > int(EMA / DT):                       # likelihood counted once the read-noise estimate has settled
            self.L_rw += ll_rw; self.L_cv += ll_cv
        self.k += 1
        est_u = self._mix(np.concatenate([e_rw, e_cv], 1), np.concatenate([self.L_rw, self.L_cv], 1))
        return dict(KF_RW=self._mix(e_rw, self.L_rw), KF_CV=self._mix(e_cv, self.L_cv), IMM=est_u)


def main():
    sha = hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest()
    print('step1_vs_kalman.py  sha256', sha, flush=True)
    t = np.arange(N) * DT; keep = t > S0.T_SKIP
    ck = os.path.join(OUT, 'step1_checkpoint.json')
    rows = json.load(open(ck))['rows'] if os.path.exists(ck) else []
    done_H = sorted({r['H'] for r in rows}); t0 = time.time()
    if done_H: print('checkpoint: H already done', done_H, flush=True)
    for H in S0.H_LIST:
        if H in done_H: continue
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
        x0 = 0.5 * (A[0] + B[0])
        eng = AllanSmoother(n, DT, S0.ENGINE['floor'], S0.ENGINE['tau_max'], S0.ENGINE['tau_rate'])
        comp = Competitors(n, x0)
        err = {k: np.zeros(n) for k in ('ENGINE', 'KF_RW', 'KF_CV', 'IMM')}; cnt = 0
        for k in range(N):
            a, b = A[k], B[k]
            xh = eng.step(a, b); out = comp.step(a, b)
            if keep[k]:
                err['ENGINE'] += (xh - X[k]) ** 2
                for name in ('KF_RW', 'KF_CV', 'IMM'): err[name] += (out[name] - X[k]) ** 2
                cnt += 1
        mse = {k: v / cnt for k, v in err.items()}
        for c, (li, seed) in enumerate(cells):
            meas = 0.5 * (A[:, c] + B[:, c])
            tau_opt, mse_opt, j = S0.sweep_opt(meas, X[:, c], keep)
            rows.append(dict(H=H, level=li, seed=seed, tau_opt=tau_opt, mse_opt=mse_opt,
                             **{f'regret_{k}': float(mse[k][c] / mse_opt - 1) for k in mse},
                             **{f'mse_{k}': float(mse[k][c]) for k in mse}))
        json.dump(dict(sha256=sha, rows=rows), open(ck, 'w'))
        print(f'H={H:.2f} done, {time.time() - t0:.0f} s', flush=True)
    json.dump(dict(sha256=sha, rows=rows), open(os.path.join(OUT, 'step1_results.json'), 'w'), indent=1)

    names = ['ENGINE', 'KF_RW', 'KF_CV', 'IMM']
    print('\n=== STEP 1: regret = MSE / MSE(best fixed exponential window) - 1, mean over 10 seeds ===')
    verdict = dict(C1=True, C2=True, C3=True, C3_fail_cells=[], C4=True)
    for H in S0.H_LIST:
        R = [r for r in rows if r['H'] == H]
        print(f'\nH = {H:.2f}')
        print('  level  tau_opt    ENGINE    KF_RW     KF_CV     IMM      | ENGINE/minKF  IMM/minKF  ENGINE/IMM')
        for li in range(len(S0.TAU_TARGETS)):
            rl = [r for r in R if r['level'] == li]
            topt = float(np.exp(np.mean([np.log(r['tau_opt']) for r in rl])))
            reg = {k: float(np.mean([r[f'regret_{k}'] for r in rl])) for k in names}
            m = {k: float(np.exp(np.mean([np.log(r[f'mse_{k}']) for r in rl]))) for k in names}
            minkf = min(m['KF_RW'], m['KF_CV'])
            e_kf, i_kf, e_i = m['ENGINE'] / minkf, m['IMM'] / minkf, m['ENGINE'] / m['IMM']
            print(f'  {li:5d}  {topt*1e3:6.1f} ms  {reg["ENGINE"]:+7.3f}  {reg["KF_RW"]:+7.3f}  {reg["KF_CV"]:+7.3f}  {reg["IMM"]:+7.3f}  |'
                  f'   {e_kf:6.3f}       {i_kf:6.3f}      {e_i:6.3f}')
            if H == 0.50: verdict['C1'] &= reg['KF_RW'] <= 0.05 and reg['ENGINE'] <= 0.10
            if topt >= 16e-3:
                if H in (0.85, 0.95): verdict['C2'] &= reg['KF_RW'] >= 0.20
                if H in (0.30, 0.50): verdict['C2'] &= reg['KF_CV'] >= 0.20
            if e_kf > 1.10: verdict['C3'] = False; verdict['C3_fail_cells'].append((H, li, round(e_kf, 3)))
            if i_kf > 1.10: verdict['C4'] = False
    print('\nVerdict:', json.dumps(verdict))
    json.dump(verdict, open(os.path.join(OUT, 'step1_verdict.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
