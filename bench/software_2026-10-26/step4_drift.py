"""STEP 4, the last of the engineering branch: slow DRIFT of the regime, the case the identification law is built for.

PRE-REGISTRATION (written before the run; SHA-256 printed at the top of the results).

Context. Steps 2-3c tested STEP changes of regime. Against a bank hardened with a change detector, the engine
ties; its parameter-free META level does not replace a detector for steps. Steps are outside the class of the
Identification Duality (Theorem 2 assumes a log random walk of the parameter). This step tests that class.

Streams. dt = 1 ms, 2^16 samples. Base regime = level 3 of step 0 at the given H. Both log sigma^2 (read-noise
PSD) and log lambda_H (signal intensity) follow independent log random walks with per-sample step variance
Lambda, so that the e-folding time of the parameter is T_e = 1 / (Lambda / dt):
  slow drift : T_e = 10 s      fast drift : T_e = 2 s
H in {0.50, 0.85}. Signal: fGn increments scaled by sqrt(lambda_k) (exact for H = 1/2, local-scale approximation
otherwise); reads: a, b = x + sigma_k * N(0,1). The random walks are reflected softly by a weak mean reversion
(rate 1/60 s) to keep the regime inside the engine's window range over 65 s. 4 cells x 10 seeds (100..109).

Methods, all fed the same two reads: ENGINE (step 0b), ENGINE_META (step 3c, law-driven noise estimator),
HARD-1s and HARD-4s (step 3 banks with detector), FIXED-BEST (the single best fixed exponential window over the
whole run, chosen with hindsight: the strongest static baseline), and LOCAL ORACLE: the best fixed window within
each 5 s block (t > 5 s), the time-varying optimum a practitioner could not know.

Measures: regret = total MSE over t > 5 s divided by the LOCAL ORACLE's total - 1, mean over seeds; the
detector's false-alarm count per stream (it should not fire on drift).

Hypotheses (fixed now)
  R1 (drift is not a step): the hardened banks' detector fires fewer than 2 times per stream on average.
  R2 (the law's home ground), three statements, each on the mean regret over 10 seeds, per cell:
      a. ENGINE_META <= ENGINE in the two fast-drift cells (the META window helps when the regime moves);
      b. the better of ENGINE / ENGINE_META is within 10 % (relative, in 1 + regret) of the BETTER of the two
         bank memories in every cell: parity with a bank whose memory happens to fit;
      c. the better of ENGINE / ENGINE_META beats the WORSE of the two bank memories by more than 10 % in at
         least 3 of the 4 cells: the law chooses its memory, a bank user must guess.
  Reading: a and b and c -> "the law's territory": parity with the right memory, no memory to choose.
           b fails -> the law is beaten even on its own ground; the branch closes on a loss.
           c fails -> memory choice does not matter in this range; the parsimony argument stands alone.

Prediction, written now (Claude): a holds in fast drift; b holds (parity); c holds in at least 3 cells.
One pass. Bank checkpoints are I/O only. Nothing changed after seeing results.
"""
import hashlib, json, os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import step0_uprime_engine as S0
from step0b_allan_identification import AllanSmoother
import step2_regime_jumps as S2
import step3_hardened_bank as S3
import step3c_meta_noise as S3c
OUT = os.path.join(HERE, 'step4'); os.makedirs(OUT, exist_ok=True)
SEEDS = list(range(100, 110)); DT, N, t = S2.DT, S2.N, S2.t
CELLS = [(10.0, 0.50), (10.0, 0.85), (2.0, 0.50), (2.0, 0.85)]        # (T_e seconds, H)
BLOCK = 5.0; REVERT = 1.0 / 60.0


def fgn(n, H, rng):
    """Exact fractional Gaussian noise (unit variance per step) via Davies-Harte."""
    k = np.arange(n + 1, dtype=float)
    g = 0.5 * (np.abs(k + 1) ** (2 * H) - 2 * np.abs(k) ** (2 * H) + np.abs(k - 1) ** (2 * H))
    r = np.concatenate([g, g[-2:0:-1]]); lam = np.maximum(np.fft.fft(r).real, 0.0); M = r.size
    W = rng.standard_normal(M) + 1j * rng.standard_normal(M)
    return (np.fft.ifft(np.sqrt(lam) * W) * np.sqrt(M)).real[:n]


def log_rw(n, Te, rng):
    """Log random walk with e-folding time Te and weak mean reversion; returns log-factor series (starts at 0)."""
    Lam = DT / Te; step = np.sqrt(Lam) * rng.standard_normal(n); z = np.empty(n); v = 0.0
    for k in range(n):
        v += step[k] - REVERT * DT * v; z[k] = v
    return z


def make_stream(Te, H, seed):
    rng_p = np.random.default_rng(seed); rng_d = np.random.default_rng(seed + 70_000); rng_n = np.random.default_rng(seed + 50_000)
    zl = log_rw(N, Te, rng_d); zs = log_rw(N, Te, np.random.default_rng(seed + 90_000))
    inc = DT ** H * fgn(N, H, rng_p) * np.exp(0.5 * zl)            # increments scaled by sqrt(lambda_k)
    x = np.concatenate([[0.0], np.cumsum(inc)[:-1]])
    s2 = S2.base_s2(H) * np.exp(zs); sig = np.sqrt(2 * s2 / DT)
    a = x + sig * rng_n.standard_normal(N); b = x + sig * rng_n.standard_normal(N)
    return x, a, b


def local_oracle_and_best(meas, x):
    """Per 5 s block (t > 5 s): best fixed-window MSE; and the single best window over the whole scored run."""
    scored = t > 5.0; nb = int((N * DT - 5.0) // BLOCK); blocks = [(5.0 + i * BLOCK, 5.0 + (i + 1) * BLOCK) for i in range(nb)]
    per_block = np.full((S0.SWEEP.size, nb), np.inf); total = np.empty(S0.SWEEP.size)
    for i, tau in enumerate(S0.SWEEP):
        e2 = (S2.ema_series(meas, tau, x[0]) - x) ** 2; total[i] = e2[scored].sum()
        for j, (b0, b1) in enumerate(blocks): per_block[i, j] = e2[(t >= b0) & (t < b1)].sum()
    return per_block.min(0).sum(), total.min(), scored


def main():
    sha = hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest()
    print('step4_drift.py  sha256', sha, flush=True)
    X, A, B, meta = [], [], [], []
    for Te, H in CELLS:
        for seed in SEEDS:
            x, a, b = make_stream(Te, H, seed); X.append(x); A.append(a); B.append(b); meta.append((Te, H, seed))
    X = np.array(X).T; A = np.array(A).T; B = np.array(B).T; n = X.shape[1]; x0 = 0.5 * (A[0] + B[0]); t0 = time.time()
    MEAS = 0.5 * (A + B); orc = np.zeros(n); fixed = np.zeros(n)
    for c in range(n): orc[c], fixed[c], scored = local_oracle_and_best(MEAS[:, c], X[:, c])
    print(f'oracles done, {time.time() - t0:.0f} s', flush=True)
    tot = {'FIXED_BEST': fixed}
    eng = AllanSmoother(n, DT, S0.ENGINE['floor'], S0.ENGINE['tau_max'], S0.ENGINE['tau_rate'])
    engm = S3c.MetaNoiseAllan(n, DT, S0.ENGINE['floor'], S0.ENGINE['tau_max'], S0.ENGINE['tau_rate'])
    E1 = np.zeros(n); E2 = np.zeros(n)
    for k in range(N):
        e1 = (eng.step(A[k], B[k]) - X[k]) ** 2; e2 = (engm.step(A[k], B[k]) - X[k]) ** 2
        if scored[k]: E1 += e1; E2 += e2
    tot['ENGINE'] = E1; tot['ENGINE_META'] = E2; print(f'engines done, {time.time() - t0:.0f} s', flush=True)
    fires = {}
    for mem in S2.MEMS:
        tag = f'{mem:g}s'; ck = os.path.join(OUT, f'ck_hard_{tag}.npz')
        if os.path.exists(ck):
            z = np.load(ck, allow_pickle=True); tot[f'HARD_{tag}'] = z['IMM']; fires[tag] = list(z['fires'])
        else:
            hb = S3.HardenedBank(n, x0, mem); Eb = np.zeros(n)
            for k in range(N):
                out = hb.step(A[k], B[k])
                if scored[k]: Eb += (out['IMM'] - X[k]) ** 2
            tot[f'HARD_{tag}'] = Eb; fires[tag] = hb.fires; np.savez(ck, IMM=Eb, fires=np.array(hb.fires, dtype=object))
        print(f'hardened bank {tag} done, {time.time() - t0:.0f} s', flush=True)
    rows = []
    for c, (Te, H, seed) in enumerate(meta):
        rows.append(dict(Te=Te, H=H, seed=seed, fires_1s=len(fires['1s'][c]), fires_4s=len(fires['4s'][c]),
                         **{f'regret_{k}': float(v[c] / orc[c] - 1) for k, v in tot.items()}))
    json.dump(dict(sha256=sha, rows=rows), open(os.path.join(OUT, 'step4_results.json'), 'w'), indent=1)
    names = ['ENGINE', 'ENGINE_META', 'HARD_1s', 'HARD_4s', 'FIXED_BEST']
    print(f'\n=== STEP 4 (sha {sha[:12]}): drift; regret vs the local 5 s oracle, mean over 10 seeds ===')
    print('  T_e     H      ENGINE   ENGINE_META   HARD_1s   HARD_4s   FIXED_BEST   fires/stream(1s, 4s)')
    R1 = True; a_ok = True; b_ok = True; c_cnt = 0
    for Te, H in CELLS:
        rs = [r for r in rows if r['Te'] == Te and r['H'] == H]
        m = {k: float(np.mean([r[f'regret_{k}'] for r in rs])) for k in names}
        f1 = float(np.mean([r['fires_1s'] for r in rs])); f4 = float(np.mean([r['fires_4s'] for r in rs]))
        print(f'  {Te:4.0f} s  {H:.2f}   {m["ENGINE"]:+7.3f}   {m["ENGINE_META"]:+7.3f}      {m["HARD_1s"]:+7.3f}   {m["HARD_4s"]:+7.3f}   {m["FIXED_BEST"]:+7.3f}        {f1:4.1f}  {f4:4.1f}')
        R1 &= (f1 < 2) and (f4 < 2)
        best_e = min(m['ENGINE'], m['ENGINE_META']); bank_best = min(m['HARD_1s'], m['HARD_4s']); bank_worst = max(m['HARD_1s'], m['HARD_4s'])
        if Te == 2.0: a_ok &= m['ENGINE_META'] <= m['ENGINE']
        b_ok &= (1 + best_e) <= 1.10 * (1 + bank_best)
        c_cnt += (1 + bank_worst) > 1.10 * (1 + best_e)
    print(f'\nR1 detector quiet under drift: {"PASS" if R1 else "FAIL"}')
    print(f'R2a META helps in fast drift: {"PASS" if a_ok else "FAIL"}   R2b parity with the better memory: {"PASS" if b_ok else "FAIL"}   R2c beats the worse memory in {c_cnt}/4 cells: {"PASS" if c_cnt >= 3 else "FAIL"}')
    json.dump(dict(R1=bool(R1), R2a=bool(a_ok), R2b=bool(b_ok), R2c_cells=int(c_cnt)), open(os.path.join(OUT, 'step4_verdict.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
