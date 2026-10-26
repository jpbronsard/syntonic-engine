"""STEP 3 of the software bench: the hardened opponent. Does the step-2 result survive a bank with change detection?

PRE-REGISTRATION (written before the run; SHA-256 printed at the top of the results).

Context. Step 2 (sha 397c22b5, replicated 4e34029f) found: the engine beats Bayesian Kalman banks with forgetting
after sensor-noise jumps (8/8 cells, both passes), loses to the 1 s bank after signal-amplitude jumps (J3),
splits on the H switch. The banks' failure on noise jumps has a known cause: for the second during which the
read-noise estimate catches up, the innovations are out of scale and the posterior runs to fast models. The
standard engineering defence is a change detector with a reset. This step adds it, with conventional settings
fixed here, and asks whether the engine's advantage survives.

Hardened bank (HARD-M, M = 1 s and 4 s): the step-2 bank (forgetting memory M) plus ONE detector with two inputs
and one action, all settings fixed now:
  noise input   : ratio of a short (50 ms) EMA of the squared read difference d^2 to the long (1 s) EMA n2.
                  Change declared when the ratio stays > 2 or < 1/2 for 100 ms.
  signal input  : the mixture's normalized innovation squared, nu^2/S, averaged over 100 ms (NIS; ~1 for a
                  consistent filter). Change declared when the average stays > 3 for 100 ms.
  action        : on a noise change, n2 is re-initialised to the short EMA (the bank gets the new noise level at
                  once); on either change, the log-likelihoods of all models are reset to zero (uniform posterior,
                  immediate re-identification with memory M). Refractory period 1 s after a reset.
  The thresholds (2, 1/2, 3) and windows (50 ms, 100 ms) are textbook values chosen before any run; they are
  not tuned on this bench.
ENGINE: unchanged (Allan identification, step 0b). FIXED-A and ORACLE as in step 2.

Hypotheses (fixed now)
  T1, T2 as in step 2 (re-convergence of the adaptive methods; stakes of the fixed window).
  T3h (the question): ENGINE against HARD-1s and HARD-4s, the 14-cell rule of step 2 (win = adapt_excess lower by
      more than 10 % relative AND recovery not longer). The step-2 claim ("re-adapts faster without a tuned
      memory") survives only if ENGINE still wins a majority of the 14 cells against BOTH hardened banks.
      Expectation stated for the record, not as a prediction: the detector should close most of the gap on noise
      jumps (J1, J2); the bank's J3 advantage should persist; the net outcome is unknown.
  D1 (detector sanity, noise jumps): in the J1 and J2 streams the detector fires within 1 s after each of the two
      noise jumps in at least 90 % of the streams; across all streams it fires fewer than 3 false alarms per
      stream outside the 2 s following a jump. (Pre-run check on three streams: the noise input fires ~110 ms
      after each noise jump; the NIS input does not fire on the amplitude or H changes, whose one-step innovation
      statistics change by only ~25 % under this read noise. D1 is therefore stated for the noise jumps, which
      are the cells where the engine won in step 2 and where the hardening matters for fairness.)

Failure reading
  T3h not supported -> the step-2 advantage was an artefact of an un-hardened opponent; report both.
  D1 fails -> the detector settings are not doing their job; the comparison is not fair and is reported as such.

Run note: the first run (sha 9086d5a3) completed and printed its table, then crashed on the final JSON write
(numpy bool); the cast was added and the script re-run from its bank checkpoints, which reprints the same table.

Seeds 100..109, then, if T3h is supported, 200..209 (environment variables as in step 2). One pass each.
"""
import hashlib, json, os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import step0_uprime_engine as S0
from step0b_allan_identification import AllanSmoother
from step1_vs_kalman import BankKF, EMA
import step2_regime_jumps as S2
OUT = os.path.join(HERE, os.environ.get('STEP3_OUT', 'step3')); os.makedirs(OUT, exist_ok=True)
SEEDS = list(range(int(os.environ.get('STEP3_SEED0', 100)), int(os.environ.get('STEP3_SEED0', 100)) + 10))
DT, N, t = S2.DT, S2.N, S2.t

SHORT, HOLD, NIS_WIN = 0.050, 0.100, 0.100
R_HI, R_LO, NIS_HI, REFRACT = 2.0, 0.5, 3.0, 1.0


class HardenedBank:
    def __init__(self, n, x0, mem):
        self.n = n; self.a_w = DT / EMA; self.a_s = DT / SHORT; self.a_nis = DT / NIS_WIN
        self.n2 = np.zeros(n); self.n2s = np.zeros(n); self.nis = np.ones(n); self.k = 0; self.rho = 1.0 - DT / mem
        self.rw = BankKF(n, 'rw', x0, np.geomspace(1e-8, 1e2, 61)); self.cv = BankKF(n, 'cv', x0, np.geomspace(1e-12, 1e4, 97))
        self.L_rw = np.zeros((n, self.rw.m)); self.L_cv = np.zeros((n, self.cv.m))
        self.hold_noise = np.zeros(n, dtype=int); self.hold_nis = np.zeros(n, dtype=int); self.last_reset = np.full(n, -1e9)
        self.fires = [[] for _ in range(n)]

    @staticmethod
    def _w(L):
        w = np.exp(L - L.max(1, keepdims=True)); return w / w.sum(1, keepdims=True)

    def step(self, a, b):
        meas = 0.5 * (a + b); d = (a - b) / np.sqrt(2.0); d2 = d * d
        self.n2 += self.a_w * (d2 - self.n2); self.n2s += self.a_s * (d2 - self.n2s)
        r = np.maximum(self.n2 / 2, 1e-30)
        e_rw, ll_rw = self.rw.step(meas, r); e_cv, ll_cv = self.cv.step(meas, r)
        est_all = np.concatenate([e_rw, e_cv], 1); L_all = np.concatenate([self.L_rw, self.L_cv], 1)
        w = self._w(L_all); est_u = (w * est_all).sum(1)
        # mixture NIS from the per-model innovations: nu_j^2/S_j weighted by the posterior
        S_rw = self.rw.P * 0 + 1  # placeholder, replaced below
        if self.k > int(EMA / DT):
            self.L_rw = self.rho * self.L_rw + ll_rw; self.L_cv = self.rho * self.L_cv + ll_cv
            # posterior-weighted normalized innovation squared: recover nu^2/S from ll = -0.5(nu^2/S + log S)
            # we do not have S here directly; use the mixture innovation against the mixture's predicted scale:
            # approximate NIS by (meas - est_u)^2 / (r + Pmix), Pmix = posterior-weighted RW predicted variance.
            Pmix = (self._w(self.L_rw) * self.rw.P).sum(1)
            nis_now = (meas - est_u) ** 2 / (r + Pmix)
            self.nis += self.a_nis * (nis_now - self.nis)
            # detector
            ratio = self.n2s / np.maximum(self.n2, 1e-30)
            self.hold_noise = np.where((ratio > R_HI) | (ratio < R_LO), self.hold_noise + 1, 0)
            self.hold_nis = np.where(self.nis > NIS_HI, self.hold_nis + 1, 0)
            tnow = self.k * DT
            fire_noise = self.hold_noise >= int(HOLD / DT); fire_nis = self.hold_nis >= int(HOLD / DT)
            fire = (fire_noise | fire_nis) & (tnow - self.last_reset > REFRACT)
            if fire.any():
                idx = np.flatnonzero(fire)
                self.n2[idx] = np.where(fire_noise[idx], self.n2s[idx], self.n2[idx])
                self.L_rw[idx] = 0.0; self.L_cv[idx] = 0.0; self.last_reset[idx] = tnow
                self.hold_noise[idx] = 0; self.hold_nis[idx] = 0; self.nis[idx] = 1.0
                for i in idx: self.fires[i].append((tnow, 'noise' if fire_noise[i] else 'nis'))
        self.k += 1
        return dict(KF_RW=(self._w(self.L_rw) * e_rw).sum(1), KF_CV=(self._w(self.L_cv) * e_cv).sum(1), IMM=est_u)


def main():
    sha = hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest()
    print('step3_hardened_bank.py  sha256', sha, '  seeds', SEEDS[0], '..', SEEDS[-1], flush=True)
    X, A, B, meta = [], [], [], []
    for kind, H in S2.CELLS:
        for seed in SEEDS:
            x, a, b, Hs = S2.make_stream(kind, H, seed); X.append(x); A.append(a); B.append(b); meta.append((kind, H, seed))
    X = np.array(X).T; A = np.array(A).T; B = np.array(B).T; n = X.shape[1]; x0 = 0.5 * (A[0] + B[0]); t0 = time.time()
    MEAS = 0.5 * (A + B); orc = np.zeros((n, 3)); orc_tau = np.zeros((n, 3)); err2 = {'FIXED_A': np.empty_like(X)}
    for c in range(n):
        for i in range(3): orc_tau[c, i], orc[c, i] = S2.oracle_segment(MEAS[:, c], X[:, c], i)
        err2['FIXED_A'][:, c] = (S2.ema_series(MEAS[:, c], orc_tau[c, 0], X[0, c]) - X[:, c]) ** 2
    eng = AllanSmoother(n, DT, S0.ENGINE['floor'], S0.ENGINE['tau_max'], S0.ENGINE['tau_rate']); E = np.empty_like(X)
    for k in range(N): E[k] = (eng.step(A[k], B[k]) - X[k]) ** 2
    err2['ENGINE'] = E; print(f'oracle, fixed, engine done, {time.time() - t0:.0f} s', flush=True)
    fires = {}
    for mem in S2.MEMS:
        tag = f'{mem:g}s'; ck = os.path.join(OUT, f'ck_hard_{tag}.npz')
        if os.path.exists(ck):
            z = np.load(ck, allow_pickle=True); [err2.__setitem__(f'{k}_{tag}', z[k]) for k in ('KF_RW', 'KF_CV', 'IMM')]
            fires[tag] = list(z['fires'])
        else:
            hb = HardenedBank(n, x0, mem); Eb = {k: np.empty_like(X) for k in ('KF_RW', 'KF_CV', 'IMM')}
            for k in range(N):
                out = hb.step(A[k], B[k])
                for name in Eb: Eb[name][k] = (out[name] - X[k]) ** 2
            fires[tag] = hb.fires
            np.savez(ck, fires=np.array(hb.fires, dtype=object), **Eb); [err2.__setitem__(f'{k}_{tag}', Eb[k]) for k in Eb]
        print(f'hardened bank {tag} done, {time.time() - t0:.0f} s', flush=True)
    # ---- detector sanity D1 (on the 1 s bank's detector; the 4 s one has the same detector)
    hit = 0; nn = 0; false = []
    for c, (kind, H, seed) in enumerate(meta):
        f = [tt for tt, _ in fires['1s'][c]]
        if kind in ('J1', 'J2'):
            nn += 1; hit += all(any(tj <= tt <= tj + 1.0 for tt in f) for tj in (S2.T1, S2.T2))
        false.append(sum(1 for tt in f if not any(tj <= tt <= tj + 2.0 for tj in (S2.T1, S2.T2))))
    D1 = (hit / nn >= 0.9) and (np.mean(false) < 3)
    # ---- measures and report (as step 2)
    rows = []
    for c, (kind, H, seed) in enumerate(meta):
        for name, Em in err2.items():
            steady = [float(np.mean(Em[S2.seg_mask(*S2.STEADY[i]), c])) for i in range(3)]
            rec = dict(kind=kind, H=H, seed=seed, method=name, steady_regret=[steady[i] / orc[c, i] - 1 for i in range(3)])
            for j, tj in enumerate((S2.T1, S2.T2)):
                tr5, rc = S2.transient_measures(Em[:, c], tj, steady[j + 1])
                rec[f'excess{j+1}'] = tr5 / orc[c, j + 1] - 1; rec[f'adapt{j+1}'] = tr5 / steady[j + 1] - 1; rec[f'recovery{j+1}'] = rc
            rows.append(rec)
    json.dump(dict(sha256=sha, seeds=SEEDS, rows=rows), open(os.path.join(OUT, 'step3_results.json'), 'w'), indent=1)
    methods = ['ENGINE', 'IMM_1s', 'IMM_4s', 'KF_RW_1s', 'KF_RW_4s', 'KF_CV_1s', 'KF_CV_4s', 'FIXED_A']
    def M(kind, H, name, key, j=None):
        rs = [r for r in rows if r['kind'] == kind and r['H'] == H and r['method'] == name]
        return float(np.mean([r[key][j] for r in rs])) if key == 'steady_regret' else float(np.mean([r[key] for r in rs]))
    print(f'\n=== STEP 3 (sha {sha[:12]}, seeds {SEEDS[0]}..{SEEDS[-1]}): hardened banks; means over 10 seeds ===')
    print(f'D1 detector: fired within 1 s after both noise jumps in {100*hit/nn:.0f} % of J1/J2 streams; false alarms per stream {np.mean(false):.2f}  -> {"PASS" if D1 else "FAIL"}')
    T1ok, T2ok = True, True; cells = []
    for kind, H in S2.CELLS:
        hl = f'H={H:.2f}' if H is not None else 'H 0.50->0.85->0.50'
        print(f'\n{kind}  {hl}'); print('  method      steady regret A / B / A\'      jump1: excess  adapt  recov(s)     jump2: excess  adapt  recov(s)')
        for name in methods:
            sr = [M(kind, H, name, 'steady_regret', i) for i in range(3)]
            v = [M(kind, H, name, k) for k in ('excess1', 'adapt1', 'recovery1', 'excess2', 'adapt2', 'recovery2')]
            print(f'  {name:10s}  {sr[0]:+6.3f} / {sr[1]:+6.3f} / {sr[2]:+6.3f}        {v[0]:+7.3f} {v[1]:+7.3f}  {v[2]:5.2f}          {v[3]:+7.3f} {v[4]:+7.3f}  {v[5]:5.2f}')
            if name in ('ENGINE', 'KF_RW_1s', 'KF_RW_4s', 'IMM_1s', 'IMM_4s') and kind != 'J4': T1ok &= all(s <= 0.10 for s in sr)
            if name == 'FIXED_A' and kind in ('J2', 'J3'): T2ok &= sr[1] >= 0.20
        for j in (1, 2):
            aE, rE = M(kind, H, 'ENGINE', f'adapt{j}'), M(kind, H, 'ENGINE', f'recovery{j}'); res = {}
            for bank in ('IMM_1s', 'IMM_4s'):
                aB, rB = M(kind, H, bank, f'adapt{j}'), M(kind, H, bank, f'recovery{j}')
                res[bank] = 'ENGINE' if (aE < 0.9 * aB and rE <= rB) else (bank if (aB < 0.9 * aE and rB <= rE) else 'tie')
            cells.append((kind, H, j, res))
    nE = {b: sum(1 for c in cells if c[3][b] == 'ENGINE') for b in ('IMM_1s', 'IMM_4s')}
    nB = {b: sum(1 for c in cells if c[3][b] == b) for b in ('IMM_1s', 'IMM_4s')}
    T3h = all(nE[b] >= 8 for b in nE)
    print(f'\nT1: {"PASS" if T1ok else "FAIL"}   T2: {"PASS" if T2ok else "FAIL"}   D1: {"PASS" if D1 else "FAIL"}')
    print(f'T3h cells (14): ENGINE wins {nE}, hardened bank wins {nB}, rest ties  ->  step-2 claim survives: {T3h}')
    for c in cells: print('   ', c)
    json.dump(dict(T1=bool(T1ok), T2=bool(T2ok), D1=bool(D1), T3h=bool(T3h), engine_wins=nE, bank_wins=nB), open(os.path.join(OUT, 'step3_verdict.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
