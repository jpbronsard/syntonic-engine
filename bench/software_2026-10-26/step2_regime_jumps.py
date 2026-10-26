"""STEP 2 of the software bench: regime jumps. Does the syntonic engine re-adapt faster than Bayesian Kalman banks?

PRE-REGISTRATION (written before the run; SHA-256 printed at the top of the results).

Context. Step 1 (sha fc13474b) settled the stationary case: parity. A likelihood-selected Kalman bank finds the
optimal window at every H without reading H, and on smooth targets (H >= 0.85) a constant-velocity Kalman beats
every first-order window. The last place an advantage of the law could appear is the transient after a change of
regime, where a bank needs a forgetting factor (a tuned memory) and the engine has its own identification law.
No directional claim is made here: this step measures, and the decision rule below says what would count.

Streams. dt = 1 ms, 2^16 samples (65.5 s), three segments: A = [0, 22 s), B = [22, 44 s), A' = [44, 65.5 s).
Jumps at t1 = 22 s and t2 = 44 s. Base regime = level 3 of step 0 (U' optimum ~ 28 ms) at the given H.
  J1  read-noise PSD s2 x4 in B (sigma x2), back in A'.          H in {0.50, 0.85}
  J2  s2 x16 in B (sigma x4), back in A'.                         H in {0.50, 0.85}
  J3  lambda_H x16 in B (increment amplitude x4), continuous path. H in {0.50, 0.85}
  J4  H switch 0.50 -> 0.85 -> 0.50 at fixed read noise, continuous path.
10 seeds (100..109, reserved for this purpose and untouched so far). 70 streams in all.

Methods, all fed the same two reads per tick.
  ENGINE    : syntonic smoother, Allan identification (step 0b), unchanged: 1 s averaging, >= 8 blocks per scale,
              window rate limit 0.4 s/s. Zero configuration.
  BANK-RW, BANK-CV, UNION (= IMM of step 1): the Bayesian banks of step 1 with exponential forgetting of the
              log-likelihood, memory M. Two settings, both reported: M = 1 s (the engine's averaging time) and
              M = 4 s. The memory is the bank's one tuning knob; two values are run so that no single choice
              decides the outcome.
  FIXED-A   : the best fixed exponential window of segment A, held fixed. What a practitioner who tuned once gets.
  ORACLE    : the best fixed exponential window per segment (sweep on that segment's steady window).

Measures, per stream, per jump (t_j = t1, t2), per method:
  trans5        mean squared error over [t_j, t_j + 5 s]
  own_steady    the method's MSE over the last 10 s of the post-jump segment
  oracle_steady ORACLE's MSE over the same window
  steady_regret own_steady / oracle_steady - 1
  excess        trans5 / oracle_steady - 1   (total transient cost against the best achievable steady state)
  adapt_excess  trans5 / own_steady - 1      (cost of adapting, separated from any steady-state class deficit)
  recovery      first t >= t_j + 0.5 s at which the running 1 s MSE <= 1.2 x own_steady (cap 12 s)

Hypotheses (fixed now)
  T1 (re-convergence): ENGINE, BANK-RW, UNION have steady_regret <= 0.10 in every segment of every cell
      (BANK-CV is exempt where its model is wrong: J1-J3 at H = 0.50, and the H = 0.50 segments of J4).
  T2 (stakes): FIXED-A has steady_regret >= 0.20 on segment B for J2 and J3 at both H.
  T3 (the question): ENGINE against UNION-1s and UNION-4s, cell by cell (7 jump-type x H cells, x 2 jumps = 14).
      A method wins a cell when its mean adapt_excess is lower by more than 10 % (relative) AND its mean recovery
      is not longer. The claim "the engine re-adapts faster without a tuned memory" is supported only if ENGINE
      wins a majority of the 14 cells against BOTH bank memories. Otherwise it is not supported.
  T4 (class switch, J4): UNION is expected to move to the CV class in the H = 0.85 segment (steady_regret < 0,
      as in step 1); ENGINE is expected to lose there in steady state (known). Only adaptation is scored for J4.

Failure reading
  T1 fails for a bank -> its forgetting is too short for that regime (report, do not retune).
  T1 fails for ENGINE -> the Allan identification does not re-converge within 12 s after that jump: a defect.
  T3 not supported -> the engineering branch is closed: the engine is a cheap window optimizer, at parity with
      banks in every regime tested.

One pass. Per-setting checkpoints (I/O only). Nothing changed after seeing results.

Replication (26 Oct 2026, after the first run, sha 397c22b5, supported T3 with 8/14 and 13/14): the same script
is re-run on fresh seeds 200..209 into step2_rep/ (environment variables STEP2_SEED0 and STEP2_OUT; nothing else
changed). The same decision rule applies; a replication that fails it is reported as a failure.
"""
import hashlib, json, os, sys, time
import numpy as np
from scipy.signal import lfilter

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import step0_uprime_engine as S0
from step0b_allan_identification import AllanSmoother
from step1_vs_kalman import BankKF, EMA
OUT = os.path.join(HERE, os.environ.get('STEP2_OUT', 'step2')); os.makedirs(OUT, exist_ok=True)

DT, N = S0.DT, S0.N
T1, T2 = 22.0, 44.0
SEEDS = list(range(int(os.environ.get('STEP2_SEED0', 100)), int(os.environ.get('STEP2_SEED0', 100)) + 10))
CELLS = [('J1', 0.50), ('J1', 0.85), ('J2', 0.50), ('J2', 0.85), ('J3', 0.50), ('J3', 0.85), ('J4', None)]
MEMS = [1.0, 4.0]
BASE_LEVEL = 3
t = np.arange(N) * DT
SEG = [(0.0, T1), (T1, T2), (T2, N * DT)]
STEADY = [(s1 - 10.0, s1) for (s0, s1) in SEG]


def base_s2(H):
    return S0.LAMBDA_H * (S0.TAU_TARGETS[BASE_LEVEL] / float(S0.kappa_H(H))) ** (2 * H + 1)


def make_stream(kind, H, seed):
    """Returns x (N,), a, b (N,), the read-noise PSD per segment and the path H per segment (for the record)."""
    rng_p = np.random.default_rng(seed); rng_n = np.random.default_rng(seed + 50_000)
    k1, k2 = int(T1 / DT), int(T2 / DT)
    if kind == 'J4':
        xa = S0.fbm(N, 0.50, rng_p, DT); xb = S0.fbm(N, 0.85, np.random.default_rng(seed + 1), DT)
        x = xa.copy(); x[k1:k2] = xa[k1] + (xb[k1:k2] - xb[k1]); x[k2:] = x[k2 - 1] + (xa[k2:] - xa[k2 - 1])
        s2 = np.full(N, base_s2(0.50)); Hs = (0.50, 0.85, 0.50)
    else:
        x = S0.fbm(N, H, rng_p, DT); s2 = np.full(N, base_s2(H)); Hs = (H, H, H)
        if kind == 'J1': s2[k1:k2] *= 4
        if kind == 'J2': s2[k1:k2] *= 16
        if kind == 'J3':
            xo = x.copy(); x[k1:k2] = xo[k1] + 4 * (xo[k1:k2] - xo[k1]); x[k2:] = x[k2 - 1] + (xo[k2:] - xo[k2 - 1])
    sig = np.sqrt(2 * s2 / DT)
    a = x + sig * rng_n.standard_normal(N); b = x + sig * rng_n.standard_normal(N)
    return x, a, b, Hs


class Competitors2:
    """Step-1 banks with exponential forgetting of the log-likelihood (memory mem, seconds)."""

    def __init__(self, n, x0, mem):
        self.n = n; self.a_w = DT / EMA; self.n2 = np.zeros(n); self.k = 0; self.rho = 1.0 - DT / mem
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
        if self.k > int(EMA / DT):
            self.L_rw = self.rho * self.L_rw + ll_rw; self.L_cv = self.rho * self.L_cv + ll_cv
        self.k += 1
        est_u = self._mix(np.concatenate([e_rw, e_cv], 1), np.concatenate([self.L_rw, self.L_cv], 1))
        return dict(KF_RW=self._mix(e_rw, self.L_rw), KF_CV=self._mix(e_cv, self.L_cv), IMM=est_u)


def ema_series(meas, tau, x0):
    g = 1.0 - np.exp(-DT / tau)
    return lfilter([g], [1.0, -(1.0 - g)], meas, zi=[(1.0 - g) * x0])[0]


def seg_mask(t0, t1):
    return (t >= t0) & (t < t1)


def oracle_segment(meas, x, i):
    """Best fixed window on segment i's steady window (EMA run from the start)."""
    m = seg_mask(*STEADY[i]); best = (None, np.inf)
    for tau in S0.SWEEP:
        e = ema_series(meas, tau, x[0]); v = np.mean((e[m] - x[m]) ** 2)
        if v < best[1]: best = (tau, v)
    return best


def transient_measures(err2, tj, own_steady):
    """err2: squared error series of one method. Returns trans5, recovery."""
    m5 = (t >= tj) & (t < tj + 5.0); trans5 = float(np.mean(err2[m5]))
    w = int(1.0 / DT); c = np.cumsum(np.concatenate([[0.0], err2]))
    run = (c[w:] - c[:-w]) / w                           # running 1 s mean, value at index i covers [i-w+1 .. i]
    i0 = int((tj + 0.5) / DT); i1 = min(int((tj + 12.0) / DT), run.size - 1)
    idx = np.flatnonzero(run[i0:i1] <= 1.2 * own_steady)
    rec = (i0 + idx[0]) * DT - tj if idx.size else 12.0
    return trans5, float(rec)


def main():
    sha = hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest()
    print('step2_regime_jumps.py  sha256', sha, flush=True)
    # ---- streams
    X, A, B, meta = [], [], [], []
    for kind, H in CELLS:
        for seed in SEEDS:
            x, a, b, Hs = make_stream(kind, H, seed); X.append(x); A.append(a); B.append(b); meta.append((kind, H, seed, Hs))
    X = np.array(X).T; A = np.array(A).T; B = np.array(B).T; n = X.shape[1]
    x0 = 0.5 * (A[0] + B[0]); t0 = time.time()
    # ---- oracle and FIXED-A
    MEAS = 0.5 * (A + B)
    orc = np.zeros((n, 3)); orc_tau = np.zeros((n, 3)); fixedA = np.zeros((n, 3)); err2 = {}
    err2['FIXED_A'] = np.empty_like(X)
    for c in range(n):
        for i in range(3):
            orc_tau[c, i], orc[c, i] = oracle_segment(MEAS[:, c], X[:, c], i)
        e = ema_series(MEAS[:, c], orc_tau[c, 0], X[0, c]); err2['FIXED_A'][:, c] = (e - X[:, c]) ** 2
    print(f'oracle + fixed done, {time.time() - t0:.0f} s', flush=True)
    # ---- engine
    ck = os.path.join(OUT, 'ck_engine.npy')
    if os.path.exists(ck): err2['ENGINE'] = np.load(ck)
    else:
        eng = AllanSmoother(n, DT, S0.ENGINE['floor'], S0.ENGINE['tau_max'], S0.ENGINE['tau_rate'])
        E = np.empty_like(X)
        for k in range(N): E[k] = (eng.step(A[k], B[k]) - X[k]) ** 2
        err2['ENGINE'] = E; np.save(ck, E)
    print(f'engine done, {time.time() - t0:.0f} s', flush=True)
    # ---- banks
    for mem in MEMS:
        tag = f'{mem:g}s'; ck = os.path.join(OUT, f'ck_bank_{tag}.npz')
        if os.path.exists(ck):
            z = np.load(ck); [err2.__setitem__(f'{k}_{tag}', z[k]) for k in ('KF_RW', 'KF_CV', 'IMM')]
        else:
            comp = Competitors2(n, x0, mem); E = {k: np.empty_like(X) for k in ('KF_RW', 'KF_CV', 'IMM')}
            for k in range(N):
                out = comp.step(A[k], B[k])
                for name in E: E[name][k] = (out[name] - X[k]) ** 2
            np.savez(ck, **E); [err2.__setitem__(f'{k}_{tag}', E[k]) for k in E]
        print(f'bank {tag} done, {time.time() - t0:.0f} s', flush=True)
    # ---- measures
    rows = []
    for c, (kind, H, seed, Hs) in enumerate(meta):
        for name, E in err2.items():
            steady = [float(np.mean(E[seg_mask(*STEADY[i]), c])) for i in range(3)]
            rec = dict(kind=kind, H=H, seed=seed, method=name,
                       steady_regret=[steady[i] / orc[c, i] - 1 for i in range(3)], oracle_tau=orc_tau[c].tolist())
            for j, tj in enumerate((T1, T2)):
                tr5, rc = transient_measures(E[:, c], tj, steady[j + 1])
                rec[f'excess{j+1}'] = tr5 / orc[c, j + 1] - 1; rec[f'adapt{j+1}'] = tr5 / steady[j + 1] - 1
                rec[f'recovery{j+1}'] = rc
            rows.append(rec)
    json.dump(dict(sha256=sha, rows=rows), open(os.path.join(OUT, 'step2_results.json'), 'w'), indent=1)

    # ---- report
    methods = ['ENGINE', 'IMM_1s', 'IMM_4s', 'KF_RW_1s', 'KF_RW_4s', 'KF_CV_1s', 'KF_CV_4s', 'FIXED_A']
    def M(kind, H, name, key, j=None):
        rs = [r for r in rows if r['kind'] == kind and r['H'] == H and r['method'] == name]
        if key == 'steady_regret': return float(np.mean([r[key][j] for r in rs]))
        return float(np.mean([r[key] for r in rs]))
    print('\n=== STEP 2 (sha %s): means over 10 seeds ===' % sha[:12])
    T1ok, T2ok = True, True; wins = {m: 0 for m in ('ENGINE', 'IMM_1s', 'IMM_4s')}; cells14 = []
    for kind, H in CELLS:
        hl = f'H={H:.2f}' if H is not None else 'H 0.50->0.85->0.50'
        print(f'\n{kind}  {hl}   oracle windows A/B/A\': ' + ' / '.join(
            f'{1e3*np.exp(np.mean([np.log(r["oracle_tau"][i]) for r in rows if r["kind"]==kind and r["H"]==H and r["method"]=="ENGINE"])):.0f} ms' for i in range(3)))
        print('  method      steady regret A / B / A\'      jump1: excess  adapt  recov(s)     jump2: excess  adapt  recov(s)')
        for name in methods:
            sr = [M(kind, H, name, 'steady_regret', i) for i in range(3)]
            e1, a1, r1 = M(kind, H, name, 'excess1'), M(kind, H, name, 'adapt1'), M(kind, H, name, 'recovery1')
            e2, a2, r2 = M(kind, H, name, 'excess2'), M(kind, H, name, 'adapt2'), M(kind, H, name, 'recovery2')
            print(f'  {name:10s}  {sr[0]:+6.3f} / {sr[1]:+6.3f} / {sr[2]:+6.3f}        {e1:+7.3f} {a1:+7.3f}  {r1:5.2f}          {e2:+7.3f} {a2:+7.3f}  {r2:5.2f}')
            if name in ('ENGINE', 'KF_RW_1s', 'KF_RW_4s', 'IMM_1s', 'IMM_4s') and kind != 'J4':
                T1ok &= all(s <= 0.10 for s in sr)
            if name == 'FIXED_A' and kind in ('J2', 'J3'): T2ok &= sr[1] >= 0.20
        for j in (1, 2):
            aE, rE = M(kind, H, 'ENGINE', f'adapt{j}'), M(kind, H, 'ENGINE', f'recovery{j}')
            res = {}
            for bank in ('IMM_1s', 'IMM_4s'):
                aB, rB = M(kind, H, bank, f'adapt{j}'), M(kind, H, bank, f'recovery{j}')
                if aE < 0.9 * aB and rE <= rB: res[bank] = 'ENGINE'
                elif aB < 0.9 * aE and rB <= rE: res[bank] = bank
                else: res[bank] = 'tie'
            cells14.append((kind, H, j, res))
    nE = {b: sum(1 for c in cells14 if c[3][b] == 'ENGINE') for b in ('IMM_1s', 'IMM_4s')}
    nB = {b: sum(1 for c in cells14 if c[3][b] == b) for b in ('IMM_1s', 'IMM_4s')}
    T3 = all(nE[b] >= 8 for b in nE)
    print(f'\nT1 re-convergence: {"PASS" if T1ok else "FAIL"}   T2 stakes: {"PASS" if T2ok else "FAIL"}')
    print(f'T3 cells (14): ENGINE wins {nE}, bank wins {nB}, rest ties  ->  claim supported: {T3}')
    for c in cells14: print('   ', c)
    json.dump(dict(T1=T1ok, T2=T2ok, T3=T3, engine_wins=nE, bank_wins=nB, cells=[list(map(str, c)) for c in cells14]),
              open(os.path.join(OUT, 'step2_verdict.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
