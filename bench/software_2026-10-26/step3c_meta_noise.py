"""STEP 3c: the symmetric question, done without parameters. The noise estimator obeys the law itself (META level).

PRE-REGISTRATION (written before the run; SHA-256 printed at the top of the results).

Context. Step 3 (sha 9086d5a3): against banks hardened with a thresholded change detector, the engine keeps the
four noise-INCREASE cells and loses the four noise-DECREASE cells. Mechanism: the engine estimates the read noise
n2 with a fixed 1 s exponential average, which decays slowly after a decrease. A detector with thresholds was
drafted for the engine (step 3b) and withdrawn before running: thresholds are parameters, and the engine's
point is to have none. The corpus's own answer is the Identification Duality (Theorem 2): the identification
window must itself follow the law, N_e* = sqrt(v / Lambda), v the relative variance of the moment stream and
Lambda its drift. This step implements it as a second syntonic engine, identical in construction, tracking the
noise level: the META level of the hierarchy, with no new constant.

ENGINE_META: the Allan smoother of step 0b, in which n2 is supplied by a second Allan smoother that tracks the
  stream y = log d^2 (d = (a - b)/sqrt(2), one value per tick):
    - y has a known additive scatter: Var(log chi^2_1) = pi^2/2, mean offset E[log chi^2_1] = -1.2704. The
      meta-engine therefore needs no second read: its per-read noise variance is set to the constant pi^2/2
      (expressed in the class's convention, n2_meta = pi^2 so that the pair-mean variance is pi^2/2).
    - its window comes from the same law (Allan structure function of y at the same scales, H_meta,
      lambda_meta, kappa_H, U' exponent); n2_hat = exp(xhat_meta + 1.2704).
    - during the meta-engine's own warm-up (1 s, the engine's existing constant) n2 falls back to the 1 s EMA.
  Nothing else changes. No threshold, no new parameter: the scales, block averaging, warm-up and rate limit are
  the engine's existing design constants, reused as they are.
Opponents: the hardened banks of step 3, from its checkpoints (same seeds 100..109, same streams).
FIXED-A, ORACLE, measures, cells and the 14-cell rule: as in steps 2 and 3. ENGINE (step 2) shown for reference.

Prediction, written now (Claude, 26 Oct 2026): ENGINE_META recovers a good part of the four noise-decrease cells,
probably not to the level of a 110 ms thresholded detector; keeps the noise-increase cells; J3 and J4 unchanged.
Expected outcome on the 14-cell rule: "matches" (at least as many wins as losses against each bank), not "beats".

Hypotheses (fixed now)
  M1 (steady state): ENGINE_META steady regret <= ENGINE's + 0.02 in every segment (the meta window must not
      make the noise estimate worse in steady state).
  T3m: 14-cell rule against HARD-1s and HARD-4s; outcomes beats / matches / loses as in step 3b's rule
      (beats: >= 8 wins against both; matches: wins >= losses against each; else loses).
Reading: this answers "engine with a law-driven identification against bank with a detector": two systems,
each hardened in its own idiom. The law-alone result remains that of step 3.
"""
import hashlib, json, os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import step0_uprime_engine as S0
import step0b_allan_identification as S0b
import step2_regime_jumps as S2
OUT = os.path.join(HERE, 'step3c'); os.makedirs(OUT, exist_ok=True)
SEEDS = list(range(100, 110)); DT, N, t = S2.DT, S2.N, S2.t
LOGCHI2_VAR, LOGCHI2_MEAN = np.pi ** 2 / 2, -1.2704


class FixedNoiseAllan(S0b.AllanSmoother):
    """Allan smoother for a single-read stream whose additive scatter variance is known (no second read)."""

    def __init__(self, n, dt, floor, tau_max, tau_rate, var_known):
        super().__init__(n, dt, floor, tau_max, tau_rate)
        self.n2_fixed = 2.0 * var_known                 # class convention: pair-mean variance = n2 / 2

    def step(self, y):
        self.n2[:] = self.n2_fixed
        return super().step(y, y)                       # a = b = y: meas = y, d = 0; n2 reset above before use


class MetaNoiseAllan(S0b.AllanSmoother):
    """Allan smoother whose n2 comes from a second (meta) Allan smoother on log d^2."""

    def __init__(self, n, dt, floor, tau_max, tau_rate):
        super().__init__(n, dt, floor, tau_max, tau_rate)
        self.meta = FixedNoiseAllan(n, dt, floor, tau_max, tau_rate, LOGCHI2_VAR)
        self.n2_ema = np.zeros(n)

    def step(self, a, b):
        d = (a - b) / np.sqrt(2.0); d2 = np.maximum(d * d, 1e-300)
        self.n2_ema += self.a_w * (d2 - self.n2_ema)
        ym = self.meta.step(np.log(d2))
        n2_meta = np.exp(ym + (-LOGCHI2_MEAN))
        use_meta = self.meta.k > self.meta.warm
        out = super().step(a, b)                        # updates self.n2 by its own 1 s EMA ...
        self.n2 = np.where(use_meta, n2_meta, self.n2_ema)   # ... then overridden by the meta estimate
        return out


def main():
    sha = hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest()
    print('step3c_meta_noise.py  sha256', sha, flush=True)
    X, A, B, meta = [], [], [], []
    for kind, H in S2.CELLS:
        for seed in SEEDS:
            x, a, b, _ = S2.make_stream(kind, H, seed); X.append(x); A.append(a); B.append(b); meta.append((kind, H, seed))
    X = np.array(X).T; A = np.array(A).T; B = np.array(B).T; n = X.shape[1]; t0 = time.time()
    MEAS = 0.5 * (A + B); orc = np.zeros((n, 3)); orc_tau = np.zeros((n, 3)); err2 = {'FIXED_A': np.empty_like(X)}
    for c in range(n):
        for i in range(3): orc_tau[c, i], orc[c, i] = S2.oracle_segment(MEAS[:, c], X[:, c], i)
        err2['FIXED_A'][:, c] = (S2.ema_series(MEAS[:, c], orc_tau[c, 0], X[0, c]) - X[:, c]) ** 2
    eng = MetaNoiseAllan(n, DT, S0.ENGINE['floor'], S0.ENGINE['tau_max'], S0.ENGINE['tau_rate']); E = np.empty_like(X)
    for k in range(N): E[k] = (eng.step(A[k], B[k]) - X[k]) ** 2
    err2['ENGINE_META'] = E
    err2['ENGINE'] = np.load(os.path.join(HERE, 'step2', 'ck_engine.npy'))
    for mem in S2.MEMS:
        tag = f'{mem:g}s'; z = np.load(os.path.join(HERE, 'step3', f'ck_hard_{tag}.npz'), allow_pickle=True)
        for k in ('KF_RW', 'KF_CV', 'IMM'): err2[f'{k}_{tag}'] = z[k]
    print(f'done computing, {time.time() - t0:.0f} s', flush=True)
    rows = []
    for c, (kind, H, seed) in enumerate(meta):
        for name, Em in err2.items():
            steady = [float(np.mean(Em[S2.seg_mask(*S2.STEADY[i]), c])) for i in range(3)]
            rec = dict(kind=kind, H=H, seed=seed, method=name, steady_regret=[steady[i] / orc[c, i] - 1 for i in range(3)])
            for j, tj in enumerate((S2.T1, S2.T2)):
                tr5, rc = S2.transient_measures(Em[:, c], tj, steady[j + 1])
                rec[f'excess{j+1}'] = tr5 / orc[c, j + 1] - 1; rec[f'adapt{j+1}'] = tr5 / steady[j + 1] - 1; rec[f'recovery{j+1}'] = rc
            rows.append(rec)
    json.dump(dict(sha256=sha, rows=rows), open(os.path.join(OUT, 'step3c_results.json'), 'w'), indent=1)
    methods = ['ENGINE_META', 'ENGINE', 'IMM_1s', 'IMM_4s', 'FIXED_A']
    def M(kind, H, name, key, j=None):
        rs = [r for r in rows if r['kind'] == kind and r['H'] == H and r['method'] == name]
        return float(np.mean([r[key][j] for r in rs])) if key == 'steady_regret' else float(np.mean([r[key] for r in rs]))
    print(f'\n=== STEP 3c (sha {sha[:12]}): engine with META noise estimator vs hardened banks; means over 10 seeds ===')
    cells = []; M1 = True
    for kind, H in S2.CELLS:
        hl = f'H={H:.2f}' if H is not None else 'H 0.50->0.85->0.50'
        print(f'\n{kind}  {hl}'); print('  method       steady regret A / B / A\'      jump1: excess  adapt  recov(s)     jump2: excess  adapt  recov(s)')
        for name in methods:
            sr = [M(kind, H, name, 'steady_regret', i) for i in range(3)]
            v = [M(kind, H, name, k) for k in ('excess1', 'adapt1', 'recovery1', 'excess2', 'adapt2', 'recovery2')]
            print(f'  {name:11s}  {sr[0]:+6.3f} / {sr[1]:+6.3f} / {sr[2]:+6.3f}        {v[0]:+7.3f} {v[1]:+7.3f}  {v[2]:5.2f}          {v[3]:+7.3f} {v[4]:+7.3f}  {v[5]:5.2f}')
        srm = [M(kind, H, 'ENGINE_META', 'steady_regret', i) for i in range(3)]; sre = [M(kind, H, 'ENGINE', 'steady_regret', i) for i in range(3)]
        M1 &= all(srm[i] <= sre[i] + 0.02 for i in range(3))
        for j in (1, 2):
            aE, rE = M(kind, H, 'ENGINE_META', f'adapt{j}'), M(kind, H, 'ENGINE_META', f'recovery{j}'); res = {}
            for bank in ('IMM_1s', 'IMM_4s'):
                aB, rB = M(kind, H, bank, f'adapt{j}'), M(kind, H, bank, f'recovery{j}')
                res[bank] = 'ENGINE_META' if (aE < 0.9 * aB and rE <= rB) else (bank if (aB < 0.9 * aE and rB <= rE) else 'tie')
            cells.append((kind, H, j, res))
    nE = {b: sum(1 for c in cells if c[3][b] == 'ENGINE_META') for b in ('IMM_1s', 'IMM_4s')}
    nB = {b: sum(1 for c in cells if c[3][b] == b) for b in ('IMM_1s', 'IMM_4s')}
    outcome = 'BEATS' if all(nE[b] >= 8 for b in nE) else ('MATCHES' if all(nE[b] >= nB[b] for b in nE) else 'LOSES')
    print(f'\nM1 steady state not worse than ENGINE (+0.02): {"PASS" if M1 else "FAIL"}')
    print(f'T3m cells (14): ENGINE_META wins {nE}, hardened bank wins {nB}, rest ties  ->  outcome: {outcome}')
    for c in cells: print('   ', c)
    json.dump(dict(M1=bool(M1), outcome=outcome, engine_wins=nE, bank_wins=nB), open(os.path.join(OUT, 'step3c_verdict.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
