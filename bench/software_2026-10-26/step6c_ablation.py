"""STEP 6c, pre-registered. Ablation of the attitude-estimator result of step 6b: does the square-root law contribute
anything beyond the new noise estimator and a well-chosen fixed horizon?

Origin. Step 6b (sha aa082341) put the window law in a Mahony attitude filter with the gyro drift as a fact measured at
rest and the accelerometer noise measured live, and found it within 2 to 4 % of the swept optimum and of a tuned
two-state Kalman filter. A critique asks three things of that result; this script answers the one that software can
answer: an ablation. Two parts of SYN-6B are separable, the live noise estimate s2 and the map from s2 to the window
(the square-root law). The ablation keeps the first and replaces the second.

Rows. Same generator, filters and scoring as steps 6 and 6b, imported from step6_quad_attitude_estimator.py and
step6b_gyro_fact.py unchanged; seeds 600..609; the drift fact lambda_hat measured at rest per seed exactly as in 6b.
  B_throttle  V = 1, 16, 1/16 (step 6), segments 120 s, steady part the last 80 s. The oracle sweep and the Kalman
              rows are read from results/step6_rows.json unchanged.
  W_wide      V = 1, 256, 1/16 (new): the low-frequency accel-tilt PSD S0 = V Q_G + S_SENS rises x221 in segment 1, so
              the optimal window (sqrt(S0/Q_G) = 1.08, 16.0, 0.47 s) moves by x15, where exponents separate. Segments
              180 s, steady part the last 120 s (the window needs about 40 s at the rate limit to climb to 16 s). Computed
              here: the oracle sweep (Mahony with 40 fixed windows 0.02 to 50 s, as in step 6), KF-OPT (33 R values,
              best per stream), KF-FIX (R at nominal, EKF2's fixed parameters), KF-CLAIR (R switched with the true level).

Variants (the ablation). One live noise estimate s2 per seed, exactly step 6b's (second differences of block means of
z_raw = accel tilt - open-loop integral of the raw gyro, l in {32, 64, 128} ms, block EMA max(1 s, 8 l), the minimum
over the three scales of D_H(l) l / 6). One drift fact lambda_hat per seed (6b). The map from s2 to the window:
  POW-p   tau = sqrt(S0_nom / lambda_hat) * (s2 / S0_nom) ^ p,   S0_nom = S0(V = 1),   p in {0, 0.25, 0.375, 0.5,
          0.625, 0.75, 1}.
All variants pass through the same point at nominal vibration (calibrated at the true nominal level, which favours the
alternatives); POW-0.5 is the square-root law, algebraically identical to SYN-6B; POW-0 is the law's nominal window
held fixed. Same floor dt/0.45, ceiling 30 s, rate limit 0.4 s/s, warm-up 1 s, start 0.5 s, Mahony branchement
(Kp = 1/tau, Ki = Kp^2 / (4 zeta^2), zeta = 1.826) as step 6b. Declared: at V = 256 the exponents p >= 0.625 reach the
30 s ceiling (targets 31.5, 62 and 239 s) and climb at the rate limit; the ablation measures each variant as it would
run in this machinery.
  BEST-FIX  the best single fixed window over a whole stream, chosen with hindsight from the oracle sweep (minimum of the
          total steady mean-square error J, below): the strongest fixed horizon a designer could have picked.

Score. RMS angle error per seed on the steady part of each segment (deg), as in step 6. Total J = sum over the three
segments of the seed-mean of RMS^2. Windows: seed medians of the steady-part median.

Hypotheses and decision rules (judge: evaluate_step6c.py, written with this script, before any run):
  H6c.1 identity: in B, POW-0.5 reproduces step 6b's SYN-6B steady RMS rows to 1e-9 deg, or the run is void.
  H6c.2 beyond a fixed horizon: in B and in W, J(POW-0.5) < J(BEST-FIX) and J(POW-0.5) < J(POW-0).
  H6c.3 beyond the noise estimator, the exponent: in W, the exponent that minimises J, p* (parabola through the grid
        minimum and its two neighbours), lies in [0.4, 0.6]; and J(POW-0.25) / J(POW-0.5) >= 1.05 and
        J(POW-0.75) / J(POW-0.5) >= 1.05.
  H6c.4 predicted against optimal horizon, kappa not refitted: in W, POW-0.5's window within x1.5 of the oracle window
        and regret <= 10 % in all three segments (seed medians, the 6b criterion).
  H6c.5 near-Kalman across the wide range: in W, mean RMS(POW-0.5) / mean RMS(KF-CLAIR) <= 1.10 in all three segments,
        and mean RMS(POW-0.5) < mean RMS(KF-FIX) in segments 1 and 2.
  Reported, not judged: the J curve over p in B and W, per-segment RMS of every row, the oracle's own exponent in W,
  BEST-FIX's window, KF-OPT, transients (first part of each segment).
Declared limits: one simulated axis, representative parameters, the generator of steps 6 and 6b (so this is not an
independent test of the horizons, which needs hardware); the anchor at the true nominal level is given to every
variant; the noise estimator is the one written after step 6. A failure is published as such.

Smoke policy: one short run of W only (2 seeds, 90 s segments, 30 s skipped) that prints only lambda_hat / Q_G, the
oracle windows against sqrt(S0/Q_G), and the runtime; never a variant or Kalman score. Any change after it is written
here as an amendment before the full run.
"""
import hashlib, json, os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import step6_quad_attitude_estimator as S6
import step6b_gyro_fact as S6B

SMOKE = bool(os.environ.get('STEP6C_SMOKE'))
DT = S6.DT
P_GRID = [0.0, 0.25, 0.375, 0.5, 0.625, 0.75, 1.0]
S0_NOM = S6.S0(1.0)
S6.SCEN['W_wide'] = dict(V=[1.0, 256.0, 1.0 / 16.0], b0=0.0)
SEG = {'B_throttle': (120.0, 40.0), 'W_wide': (180.0, 60.0)}
SEEDS = S6.SEEDS
if SMOKE:
    SEG = {'W_wide': (90.0, 30.0)}; SEEDS = [600, 601]


class LawPow:
    """the step-6b window machinery with the map tau = anchor * (s2 / S0_nom) ^ p, per column."""

    def __init__(self, anchor, p):
        n = len(anchor)
        self.anchor = anchor; self.p = p
        self.tau = np.full(n, S6.TAU0); self.tau_t = np.full(n, S6.TAU0)
        self.k = 0; self.warm = int(round(S6.WARMUP / DT))

    def update(self, s2):
        if self.k > self.warm:
            if s2 is not None:
                self.tau_t = np.clip(self.anchor * (np.maximum(s2, 1e-300) / S0_NOM) ** self.p, S6.FLOOR_TAU, S6.TAU_MAX)
            rr = 1.0 / self.tau; dr = S6.TAU_RATE * DT / self.tau ** 2
            rr = rr + np.clip(1.0 / self.tau_t - rr, -dr, dr); self.tau = 1.0 / rr
        self.k += 1
        return self.tau


def run(scn, lam_hat, full_competitors):
    S6.SEG_S, S6.SKIP_S = SEG[scn]
    d = S6.generate(scn, SEEDS); ns = len(SEEDS); N, nseg = d['N'], d['nseg']; skip = int(round(S6.SKIP_S / DT))
    TH, W, A1, Vk = d['TH'], d['W'], d['A1'], d['Vk']
    nP = len(P_GRID); rep = lambda x, m: np.tile(x, m)
    anchor = rep(np.sqrt(S0_NOM / lam_hat), nP); pcol = np.repeat(P_GRID, ns)
    mah = S6.Mahony(rep(A1[0], nP)); had = S6B.HadNoise(ns); law = LawPow(anchor, pcol)
    Th_raw = A1[0].copy()
    sc = S6.Score(nP * ns, nseg, skip)
    if full_competitors:
        nt, nr = len(S6.TAU_GRID), len(S6.R_MULT); R_nom = S6.S0(1.0) / DT
        orc = S6.Mahony(rep(A1[0], nt)); tau_or = np.repeat(S6.TAU_GRID, ns)
        Ropt = np.repeat(S6.R_MULT * R_nom, ns); kopt = S6.Kalman2(rep(A1[0], nr), Ropt)
        kfix = S6.Kalman2(A1[0], R_nom); kcl = S6.Kalman2(A1[0], S6.S0(Vk[0]) / DT)
        sco = {k: S6.Score(c, nseg, skip) for k, c in (('orc', nt * ns), ('kopt', nr * ns), ('kfix', ns), ('kcl', ns))}
    tlog, klog = [], []
    t0 = time.time()
    for k in range(N):
        w, a, th, V = W[k], A1[k], TH[k], Vk[k]
        Th_raw = Th_raw + w * DT
        s2 = had.update(a - Th_raw)
        tau = law.update(None if s2 is None else rep(s2, nP))
        if not SMOKE:
            sc.add(k, mah.step(rep(w, nP), rep(a, nP), tau) - rep(th, nP))
        if full_competitors:
            sco['orc'].add(k, orc.step(rep(w, nt), rep(a, nt), tau_or) - rep(th, nt))
            if not SMOKE:
                sco['kopt'].add(k, kopt.step(rep(w, nr), rep(a, nr), Ropt) - rep(th, nr))
                sco['kfix'].add(k, kfix.step(w, a, R_nom) - th)
                sco['kcl'].add(k, kcl.step(w, a, S6.S0(V) / DT) - th)
        if k % 100 == 0:
            tlog.append(law.tau.copy()); klog.append(k)
        if k % 90000 == 0 and k:
            print(f'  {scn} t = {k * DT:.0f} s  ({time.time() - t0:.0f} s)', flush=True)
    klog = np.array(klog); T = np.array(tlog)
    res = {'scenario': scn, 'V': S6.SCEN[scn]['V'], 'seg_s': S6.SEG_S, 'skip_s': S6.SKIP_S, 'seeds': SEEDS,
           'tau_star_true': [float(np.sqrt(S6.S0(v) / S6.Q_G)) for v in S6.SCEN[scn]['V']], 'runtime_s': time.time() - t0}
    if not SMOKE:
        st, tr = sc.rms()
        for j, p in enumerate(P_GRID):
            cols = slice(j * ns, (j + 1) * ns)
            med = [np.median(T[(klog // nseg == g) & ((klog % nseg) >= skip)][:, cols], axis=0).tolist() for g in range(3)]
            res[f'pow_{p}'] = {'p': p, 'steady': st[:, cols].tolist(), 'trans': tr[:, cols].tolist(), 'tau_med': med}
    if full_competitors:
        for key, s in sco.items():
            stc, trc = s.rms(); res[key] = {'steady': stc.tolist(), 'trans': trc.tolist()}
        tor, ror = S6.oracle_window(res, ns)
        res['oracle'] = {'tau': tor.tolist(), 'rms': ror.tolist()}
        res['grids'] = {'tau': S6.TAU_GRID.tolist(), 'R_mult': S6.R_MULT.tolist()}
    return res


def main():
    sha = hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest()
    print('step6c_ablation.py  sha256', sha, flush=True)
    lam_hat, slope = S6B.allan_lambda(SEEDS)
    print(f'Allan at rest: lambda_hat / Q_G per seed {np.round(lam_hat / S6.Q_G, 3).tolist()}', flush=True)
    out = {'sha256': sha, 'smoke': SMOKE, 'P_grid': P_GRID, 'lambda_hat': lam_hat.tolist(), 'allan_slope': slope.tolist()}
    for scn in SEG:
        r = run(scn, lam_hat, full_competitors=(scn == 'W_wide')); out[scn] = r
        msg = f'{scn}: runtime {r["runtime_s"]:.0f} s'
        if 'oracle' in r:
            msg += (f'; oracle window per segment (seed median) {np.round(np.median(np.array(r["oracle"]["tau"]), 1), 3).tolist()}'
                    f' vs sqrt(S0/Q_G) {np.round(r["tau_star_true"], 3).tolist()}')
        print(msg, flush=True)
    fn = os.path.join(S6.OUT, 'step6c_smoke.json' if SMOKE else 'step6c_rows.json')
    json.dump(out, open(fn, 'w')); print('written', fn)


if __name__ == '__main__':
    main()
