"""STEP 6b, pre-registered. The attitude estimator with the gyro drift as a fact and the accel noise measured live.

Origin. Step 6 (sha 035647c9) held H6.1 (the law with the true noise levels lands on the oracle and sits at parity
with a tuned two-state Kalman) and failed H6.2 to H6.5, H6.7: the engine's identification of the gyro drift lambda
ran away through the Mahony bias integral, and lambda is visible only at scales of several seconds. The exploratory
diagnosis (step-6 note, section 4) found the accel noise level measured within 12 % from second differences of block
means at 32 to 128 ms. Step 6b tests the configuration that follows: lambda is a fact of the gyro, the noise level is
measured live, the law is unchanged.

Rows. Same generator, same scenarios (A_steady, B_throttle, C_gyrobias), same seeds 600..609, same scoring (steady
parts: last 80 s of each 120 s segment) as step 6, imported from step6_quad_attitude_estimator.py. Oracle, Kalman
filters, Madgwick and CF-0.5 rows are read from results/step6_rows.json unchanged. Identity check: CF-0.5 is
recomputed here and must match the stored rows to 1e-9 deg, or the run is void.

Candidate SYN-6B (Mahony, Ki = Kp^2 / (4 zeta^2), zeta = 1.826, as in step 6):
  lambda  the gyro's angle random walk, measured once at rest: a 300 s static gyro record per seed (rng seed + 10000,
          same gyro model, motors off), Allan variance of the rate at l in {0.128, 0.256, 0.512, 1.024} s,
          lambda_hat = mean of sigma_A^2(l) l. H = 1/2: an angle random walk is a Wiener process (Allan slope -1/2,
          reported). This is the sensor fact; nothing is fitted in flight.
  sigma^2 measured live on z_raw = accel tilt - open-loop integral of the RAW gyro (no bias feedback): second
          differences of block means (Hadamard), l in {32, 64, 128} ms, each averaged with the engine's block EMA
          (time constant max(1 s, 8 l)); s2 = min over the three scales of D_H(l) l / 6 (white noise: D_H = 6 v dt / l;
          a constant or slowly varying gyro bias, a ramp in z_raw, cancels in second differences).
  law     tau = kappa_{1/2} sqrt(s2 / lambda_hat) = sqrt(s2 / lambda_hat); floor dt/0.45, ceiling 30 s, rate limit
          0.4 s/s, warm-up 1 s, start 0.5 s (engine constants).
Sensitivity variants (the fact off by a factor two, as a data sheet can be): SYN-6B x2 (lambda_hat x 2) and
SYN-6B /2 (lambda_hat / 2).

Hypotheses and decision rules (cell = scenario x segment; seed medians for windows and regret; means for RMS ratios):
  H6b.1 main: SYN-6B window within x1.5 of the step-6 oracle window and regret <= 10 %, all 6 cells of A and B.
  H6b.2 exponent: in B, d log tau / d log S0 (seg 1 vs 0, seg 2 vs 0) = 0.5 +/- 0.1.
  H6b.3 parity with a tuned Kalman: in A, mean RMS(SYN-6B) / mean RMS(KF-OPT) <= 1.10, each segment.
  H6b.4 adaptation: in B, mean RMS(SYN-6B) < mean RMS(KF-FIX) in segments 1 and 2, and mean RMS(SYN-6B) /
        mean RMS(KF-CLAIR) <= 1.10 in all three segments.
  H6b.5 bias: in C, SYN-6B passes the H6b.1 criterion in all three cells.
  H6b.6 robustness of the fact: SYN-6B x2 and SYN-6B /2 each pass the H6b.1 criterion in all 6 cells of A and B.
  Reported, not judged: lambda_hat / Q_G and the Allan slope; transients (first 40 s of each segment) against KF-FIX
  and KF-CLAIR; window traces.
Declared limits: H is fixed at 1/2 by the sensor physics, so the U' exponent measurement is not exercised here; in
flight the gyro's low-frequency noise may rise with throttle, which a fact measured at rest would miss (the gimbal
measures it); parameters are representative, as in step 6. A failure is published as such.

Smoke policy: one short run (2 seeds, 40 s segments) that prints only lambda_hat / Q_G and the runtime, never a score.
Any change after it is written here as an amendment before the full run.
"""
import hashlib, json, os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import step6_quad_attitude_estimator as S6

SMOKE = bool(os.environ.get('STEP6_SMOKE'))
DT = S6.DT
HSCALES = np.array([32, 64, 128])
ALLAN_L = np.array([0.128, 0.256, 0.512, 1.024]); REST_S = 300.0


def allan_lambda(seeds):
    """lambda_hat and Allan slope from a static gyro record per seed (same gyro model, motors off)."""
    lam, slope = [], []
    Lfit = np.array([0.016, 0.032, 0.064, 0.128, 0.256, 0.512, 1.024, 2.048, 4.096])
    for s in seeds:
        r = np.random.default_rng(s + 10000); N = int(REST_S / DT)
        w = S6.SCEN['A_steady']['b0'] + np.cumsum(np.sqrt(S6.Q_B * DT) * r.standard_normal(N)) \
            + np.sqrt(S6.Q_G / DT) * r.standard_normal(N)
        def av(l):
            k = int(round(l / DT)); m = w[: N // k * k].reshape(-1, k).mean(1); return 0.5 * np.mean(np.diff(m) ** 2)
        lam.append(np.mean([av(l) * l for l in ALLAN_L]))
        slope.append(np.polyfit(np.log(Lfit), np.log([av(l) for l in Lfit]), 1)[0] / 2)   # Allan deviation slope
    return np.array(lam), np.array(slope)


class HadNoise:
    """live noise level of z_raw from second differences of block means at 32, 64, 128 ms (engine block EMA)."""

    def __init__(self, n):
        m = len(HSCALES)
        self.acc = np.zeros((m, n)); self.cnt = np.zeros(m, dtype=int)
        self.p1 = np.full((m, n), np.nan); self.p2 = np.full((m, n), np.nan)
        self.DH = np.zeros((m, n)); self.nb = np.zeros(m, dtype=int)
        self.a_blk = np.array([min(1.0, DT * s / max(S6.EMA, 8 * DT * s)) for s in HSCALES])

    def update(self, z):
        self.acc += z[None, :]; self.cnt += 1
        for i in np.flatnonzero(self.cnt == HSCALES):
            mean = self.acc[i] / HSCALES[i]
            if not np.isnan(self.p2[i, 0]):
                dd = mean - 2 * self.p1[i] + self.p2[i]
                aw = self.a_blk[i] if self.nb[i] > 0 else 1.0
                self.DH[i] += aw * (dd * dd - self.DH[i]); self.nb[i] += 1
            self.p2[i] = self.p1[i]; self.p1[i] = mean; self.acc[i] = 0.0; self.cnt[i] = 0
        if (self.nb >= 2).all():
            return np.min(self.DH * (HSCALES * DT)[:, None] / 6.0, axis=0)
        return None


class LawTau:
    def __init__(self, n, lam):
        self.lam = lam; self.tau = np.full(n, S6.TAU0); self.tau_t = np.full(n, S6.TAU0)
        self.k = 0; self.warm = int(round(S6.WARMUP / DT))

    def update(self, s2):
        if self.k > self.warm:
            if s2 is not None:
                self.tau_t = np.clip(np.sqrt(np.maximum(s2, 1e-300) / self.lam), S6.FLOOR_TAU, S6.TAU_MAX)
            rr = 1.0 / self.tau; dr = S6.TAU_RATE * DT / self.tau ** 2
            rr = rr + np.clip(1.0 / self.tau_t - rr, -dr, dr); self.tau = 1.0 / rr
        self.k += 1
        return self.tau


def run(scn, lam_hat):
    d = S6.generate(scn, S6.SEEDS); ns = len(S6.SEEDS); N, nseg = d['N'], d['nseg']; skip = int(round(S6.SKIP_S / DT))
    TH, W, A1 = d['TH'], d['W'], d['A1']
    var = [('6b', 1.0), ('6b_x2', 2.0), ('6b_d2', 0.5)]; nv = len(var)
    lam = np.concatenate([lam_hat * f for _, f in var])
    rep = lambda x: np.tile(x, nv)
    mah = S6.Mahony(rep(A1[0])); had = HadNoise(nv * ns); law = LawTau(nv * ns, lam)
    Th_raw = A1[0].copy()
    cf05 = S6.Mahony(A1[0], withI=False)
    sc = S6.Score(nv * ns, nseg, skip); sc_cf = S6.Score(ns, nseg, skip)
    tlog, klog, slog = [], [], []
    t0 = time.time()
    for k in range(N):
        w, a, th = W[k], A1[k], TH[k]
        Th_raw = Th_raw + w * DT
        s2 = had.update(rep(a - Th_raw))
        tau = law.update(s2)
        sc.add(k, mah.step(rep(w), rep(a), tau) - rep(th))
        sc_cf.add(k, cf05.step(w, a, 0.5) - th)
        if k % 100 == 0:
            tlog.append(law.tau.copy()); slog.append(law.tau_t.copy()); klog.append(k)
    klog = np.array(klog); T = np.array(tlog)
    st, tr = sc.rms(); cst, _ = sc_cf.rms()
    res = {'scenario': scn, 'runtime_s': time.time() - t0, 'cf05_steady': cst.tolist()}
    for j, (name, _) in enumerate(var):
        cols = slice(j * ns, (j + 1) * ns)
        med = []
        for g in range(3):
            m = (klog // nseg == g) & ((klog % nseg) >= skip); med.append(np.median(T[m][:, cols], axis=0).tolist())
        res[name] = {'steady': st[:, cols].tolist(), 'trans': tr[:, cols].tolist(), 'tau_med': med}
    res['tau_trace_6b'] = {'k': klog[::10].tolist(), 'tau': T[::10, :ns].tolist()}
    return res


def main():
    sha = hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest()
    print('step6b_gyro_fact.py  sha256', sha, flush=True)
    lam_hat, slope = allan_lambda(S6.SEEDS)
    print(f'Allan at rest: lambda_hat / Q_G per seed {np.round(lam_hat / S6.Q_G, 3).tolist()};'
          f' Allan-deviation slope median {np.median(slope):.3f} (ARW: -0.5)', flush=True)
    out = {'sha256': sha, 'smoke': SMOKE, 'lambda_hat': lam_hat.tolist(), 'allan_slope': slope.tolist()}
    for scn in S6.SCEN:
        r = run(scn, lam_hat); out[scn] = r
        print(f'{scn}: runtime {r["runtime_s"]:.0f} s', flush=True)
    fn = os.path.join(S6.OUT, 'step6b_smoke.json' if SMOKE else 'step6b_rows.json')
    json.dump(out, open(fn, 'w')); print('written', fn)


if __name__ == '__main__':
    main()
