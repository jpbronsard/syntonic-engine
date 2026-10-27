"""STEP 6d, pre-registered. Two questions left by step 6c (sha d7c8b23b), on the attitude estimator in software.

(M) Mechanism. In step 6c the square-root law set the window 1.65 times too long when the vibration rose x256
    (13.3 s against 8.0 s; H6c.4 failed), and the best exponent there was 0.375 to 0.4. Attributed, untested: the
    gyro's bias random walk acting through the Mahony bias integral tied to the window. Arm M tests the attribution.
(U) Does the generalized law U' predict a change of exponent? Theorem U' gives tau* = kappa_H (sigma^2 / lambda_H)
    ^(1 / (2H + 1)) for a first-order window tracking a drift with structure function lambda_H s^(2H), H in (0, 1),
    through white noise: the exponent is 1/(2H + 1), below 1/2 when H > 1/2. Step 6c's setting lies outside U' (a
    bias random walk gives an H = 3/2 component, and Mahony with a tied bias integral has two coupled timescales).
    Arm U puts the attitude estimator inside U': a gyro whose drift is fractional (fBm, H = 0.6 and 0.75), no bias,
    a first-order window (Mahony without bias integral), and asks whether the optimal window's exponent follows
    1/(2H + 1), and whether the law with H and lambda_H measured once at rest lands on the optimum.

Common to both arms: the generator, Mahony filter, scoring and constants of steps 6, 6b and 6c (imported from
step6_quad_attitude_estimator.py and step6b_gyro_fact.py unchanged), scenario W of step 6c: vibration V = 1, 256, 1/16
over three segments of 180 s, steady part the last 120 s, seeds 600..609; accelerometer noise and attitude truth drawn
exactly as in step 6c (same random streams). Oracle: Mahony with fixed windows, 40 values 0.02 to 50 s, best per
segment (step 6). Window criterion (as H6b.1, H6c.4): window within x1.5 of the oracle window and regret <= 10 %,
seed medians, per segment. Exponent of a window between two segments: d log tau / d log S0, S0 = V Q_G + S_SENS.

Arm M (gyro of step 6: angle random walk Q_G, bias random walk Q_B; the law of step 6b: tau = sqrt(s2 / lambda_hat),
lambda_hat the Allan measurement at rest of 6b, s2 the live noise estimate of 6b, engine constants of 6b):
  REPLAY  the step-6c data and filter (Mahony with the bias integral tied to the window): must reproduce step 6c's
          POW-0.5 rows (identity check).
  M0      the same data with the bias random walk switched off (Q_B = 0, every other random draw unchanged), Mahony
          with the tied bias integral, the 6b law.
  M1      the step-6c data (bias random walk present); the filter is given the true gyro bias and has no bias
          integral (Mahony without integral on the bias-corrected gyro), the 6b law. The noise estimate s2 is built on
          the raw gyro, as in 6b and 6c, so it is identical to step 6c's.
  Oracles: M0 with the tied integral; M1 without integral, bias-corrected.

Arm U (gyro drift = fractional Brownian motion of index H and intensity lambda_H, D(s) = lambda_H s^(2H); no angle
random walk, no bias; Mahony without bias integral, which is exactly a first-order window on z = accel - integral of
the gyro, the class of U'). lambda_H is set so that the U' window at nominal vibration equals step 6's 1.08 s.
Scenarios U60 (H = 0.6) and U75 (H = 0.75). Predicted optimal windows sqrt-free: tau* = kappa_H (S0 / lambda_H)
^(1/(2H+1)), exponent 0.455 (U60) and 0.400 (U75).
  ULAW-K  U' with the true H, lambda_H and S0 (clairvoyant inputs, no rate limit, as step 6's SYN-K).
  ULAW    U' with H and lambda_H measured once at rest (a 300 s static record per seed: block means of the integrated
          gyro at 0.128 to 8.192 s, first differences, log-log fit; lambda_H = intercept / C(H), the engine's C(H))
          and s2 measured live (6b's estimator); 6b's machinery (floor, ceiling 30 s, rate limit, warm-up, start).
  SQRT    the square-root law of 6b on the same gyro: lambda from the Allan variance of the rest record at 0.128 to
          1.024 s (6b's formula), tau = sqrt(s2 / lambda), same machinery.

Hypotheses and decision rules (judge: evaluate_step6d.py, written with this script, before any run):
  H6d.0 identity: REPLAY reproduces step 6c's POW-0.5 steady rows in W to 1e-9 deg, or the run is void.
  H6d.1 (M0, no bias walk): the 6b law passes the window criterion in all three segments, and the oracle exponent from
        V = 1 to 256 lies in [0.45, 0.55].
  H6d.2 (M1, bias walk present, true bias given, no integral): the same two criteria.
  H6d.3 (U', the formula): in U60 and U75, ULAW-K passes the window criterion in all three segments.
  H6d.4 (U', the exponent): in U60 and U75, the oracle exponent from V = 1 to 256 lies within 0.05 of 1/(2H+1).
  H6d.5 (U', the practical law): in U60 and U75, ULAW passes the window criterion in all three segments.
  H6d.6 (control): in U75 the square-root law SQRT fails the window criterion in the V = 256 segment. If it passes,
        this arm cannot tell the exponents apart.
  Reported, not judged: H and lambda_H measured at rest against the truth; every window, regret and exponent; SQRT in
  U60; the 6c rows for reference.
Reading rules, declared: if H6d.1 and H6d.2 hold, the 6c failure is attributed in software to the gyro's bias random
walk acting through the bias channel; if H6d.1 fails, its cause lies elsewhere (the noise estimate or the vibration
model at x256). If H6d.4 holds, U' predicts the change of exponent on this substrate when the drift is in its class;
if H6d.5 holds as well, the engine can use it with facts measured once at rest.
Declared limits: one simulated axis, representative parameters, the generator of steps 6 to 6c (not an independent
test); fBm is a stand-in for gyro low-frequency noise (real gyros show angle random walk, bias instability and rate
random walk, not a pure power law); arm U has no bias channel; M1 gives the filter the true bias, a diagnostic, not
a filter one can build. A failure is published as such.

Smoke policy: one short run of U75 only (2 seeds, 90 s segments, 30 s skipped) that prints only the fGn variance
check, H and lambda_H measured at rest against the truth, whether each oracle optimum lies strictly inside its grid
(yes or no), and the runtime; never a window, an oracle value or a score. Any change after it is written here as an
amendment before the full run.
"""
import hashlib, json, os, sys, time
import numpy as np
from scipy.signal import butter, lfilter

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import step6_quad_attitude_estimator as S6
import step6b_gyro_fact as S6B

SMOKE = bool(os.environ.get('STEP6D_SMOKE'))
DT = S6.DT
V_W = [1.0, 256.0, 1.0 / 16.0]
SEG_S, SKIP_S, SEEDS = 180.0, 60.0, list(range(600, 610))
REST_S = 300.0
L_REST_U = 0.128 * 2.0 ** np.arange(7)                   # 0.128 .. 8.192 s
TAU_NOM = float(np.sqrt(S6.S0(1.0) / S6.Q_G))            # 1.077 s, step 6's nominal window
H_U = {'U60': 0.6, 'U75': 0.75}
ARMS = ['REPLAY', 'M0', 'M1', 'U60', 'U75']
if SMOKE:
    SEG_S, SKIP_S, SEEDS = 90.0, 30.0, [600, 601]; ARMS = ['U75']


def lam_for(H):
    """lambda_H such that the U' window at nominal vibration equals TAU_NOM."""
    return S6.S0(1.0) * S6.kappa_H(H) ** (2 * H + 1) / TAU_NOM ** (2 * H + 1)


def fgn(n, H, rng):
    """unit-variance fractional Gaussian noise of index H (Davies-Harte, circulant embedding)."""
    k = np.arange(n + 1, dtype=float)
    g = 0.5 * ((k + 1) ** (2 * H) - 2 * k ** (2 * H) + np.abs(k - 1) ** (2 * H))
    c = np.concatenate([g, g[-2:0:-1]]); m = len(c)
    lam = np.fft.fft(c).real
    assert lam.min() > -1e-8 * lam.max()
    Z = rng.standard_normal(m) + 1j * rng.standard_normal(m)
    return np.fft.fft(np.sqrt(np.maximum(lam, 0) / m) * Z).real[:n]


def generate(seeds, q_b=S6.Q_B, fbm=None):
    """step 6's generator (same draws in the same order), with the bias walk scaled by q_b and, if fbm = (H, lam_H),
    the angle random walk and bias replaced by a fractional drift (own stream, seed + 20000)."""
    V = np.array(V_W); nseg = int(round(SEG_S / DT)); N = 3 * nseg
    Vk = np.repeat(V, nseg)
    bb, ab = butter(2, 80.0 / (0.5 / DT))
    rho = np.exp(-DT / S6.TRUTH_TC); c = S6.TRUTH_SD * np.sqrt(1 - rho ** 2)
    TH, W, A1, BIAS = (np.empty((N, len(seeds))) for _ in range(4))
    for j, s in enumerate(seeds):
        r = np.random.default_rng(s)
        th = lfilter([c], [1, -rho], r.standard_normal(N + 1))
        bias = 0.0 + np.cumsum(np.sqrt(q_b * DT) * r.standard_normal(N))
        arw = np.sqrt(S6.Q_G / DT) * r.standard_normal(N)
        if fbm is not None:
            H, lamH = fbm
            gerr = np.sqrt(lamH) * DT ** H * fgn(N, H, np.random.default_rng(s + 20000)) / DT
            bias = np.zeros(N)
        else:
            gerr = bias + arw
        W[:, j] = np.diff(th) / DT + gerr
        ph = 2 * np.pi * np.cumsum(S6.LINE_F * Vk ** 0.25) * DT + r.uniform(0, 2 * np.pi)
        vib = np.sqrt(V.repeat(nseg) * S6.Q_G / DT) * r.standard_normal(N) + S6.LINE_AMP * np.sqrt(Vk) * np.sin(ph)
        e1 = np.sqrt(S6.S_SENS / DT) * r.standard_normal(N)
        A1[:, j] = th[1:] + lfilter(bb, ab, vib + e1)
        TH[:, j] = th[1:]; BIAS[:, j] = bias
    return dict(TH=TH, W=W, A1=A1, BIAS=BIAS, Vk=Vk, nseg=nseg, N=N)


def rest_facts_fbm(seeds, H, lamH):
    """300 s at rest per seed (stream seed + 30000): H and lambda_H from block means of the integrated gyro (first
    differences, log-log fit over 0.128 .. 8.192 s), and the 6b Allan lambda of the rate at 0.128 .. 1.024 s."""
    Hh, Lh, Lsq = [], [], []
    for s in seeds:
        n = int(REST_S / DT)
        x = np.sqrt(lamH) * DT ** H * fgn(n, H, np.random.default_rng(s + 30000))
        ang = np.cumsum(x); rate = x / DT
        D = []
        for l in L_REST_U:
            k = int(round(l / DT)); mb = ang[: n // k * k].reshape(-1, k).mean(1); D.append(np.mean(np.diff(mb) ** 2))
        sl, ic = np.polyfit(np.log(L_REST_U), np.log(D), 1)
        h = sl / 2; Hh.append(h); Lh.append(np.exp(ic) / S6.C_H(h))
        av = []
        for l in S6B.ALLAN_L:
            k = int(round(l / DT)); mb = rate[: n // k * k].reshape(-1, k).mean(1); av.append(0.5 * np.mean(np.diff(mb) ** 2) * l)
        Lsq.append(np.mean(av))
    return np.array(Hh), np.array(Lh), np.array(Lsq)


class Law:
    """6b's window machinery with tau_target = kappa * (s2 / lam) ^ e per column (kappa, lam, e arrays)."""

    def __init__(self, kappa, lam, e):
        n = len(lam); self.kappa, self.lam, self.e = kappa, lam, e
        self.tau = np.full(n, S6.TAU0); self.tau_t = np.full(n, S6.TAU0)
        self.k = 0; self.warm = int(round(S6.WARMUP / DT))

    def update(self, s2):
        if self.k > self.warm:
            if s2 is not None:
                self.tau_t = np.clip(self.kappa * (np.maximum(s2, 1e-300) / self.lam) ** self.e, S6.FLOOR_TAU, S6.TAU_MAX)
            rr = 1.0 / self.tau; dr = S6.TAU_RATE * DT / self.tau ** 2
            rr = rr + np.clip(1.0 / self.tau_t - rr, -dr, dr); self.tau = 1.0 / rr
        self.k += 1
        return self.tau


def run(arm, lam6b):
    ns = len(SEEDS); rep = lambda x, m: np.tile(x, m)
    if arm in ('REPLAY', 'M1'):
        d = generate(SEEDS)
    elif arm == 'M0':
        d = generate(SEEDS, q_b=0.0)
    else:
        H = H_U[arm]; lamH = lam_for(H); d = generate(SEEDS, fbm=(H, lamH))
    TH, W, A1, BIAS, Vk, N, nseg = d['TH'], d['W'], d['A1'], d['BIAS'], d['Vk'], d['N'], d['nseg']
    skip = int(round(SKIP_S / DT)); nt = len(S6.TAU_GRID)
    res = {'arm': arm, 'V': V_W, 'seeds': SEEDS, 'seg_s': SEG_S, 'skip_s': SKIP_S}
    withI = arm in ('REPLAY', 'M0')
    corr = arm == 'M1'
    if arm in ('U60', 'U75'):
        H = H_U[arm]; lamH = lam_for(H)
        Hh, Lh, Lsq = rest_facts_fbm(SEEDS, H, lamH)
        res.update(H_true=H, lamH_true=lamH, H_rest=Hh.tolist(), lamH_rest=Lh.tolist(), lam_sqrt_rest=Lsq.tolist(),
                   tau_pred=[float(S6.kappa_H(H) * (S6.S0(v) / lamH) ** (1 / (2 * H + 1))) for v in V_W])
        if SMOKE:
            print(f'  {arm}: fGn variance check {np.var(fgn(2 ** 16, H, np.random.default_rng(1))):.3f};'
                  f' H at rest {np.round(Hh, 3).tolist()} (true {H}); lambda_H at rest / true {np.round(Lh / lamH, 3).tolist()}', flush=True)
        names = ['ULAW', 'SQRT']
        law = Law(np.concatenate([S6.kappa_H(Hh), np.ones(ns)]), np.concatenate([Lh, Lsq]),
                  np.concatenate([1 / (2 * Hh + 1), np.full(ns, 0.5)]))
        nl = 2
    else:
        names = ['LAW6B']; nl = 1
        law = Law(np.ones(ns), lam6b, np.full(ns, 0.5))
    mah = S6.Mahony(rep(A1[0], nl), withI=withI); had = S6B.HadNoise(ns)
    orc = S6.Mahony(rep(A1[0], nt), withI=withI); tau_or = np.repeat(S6.TAU_GRID, ns)
    sc = S6.Score(nl * ns, nseg, skip); sco = S6.Score(nt * ns, nseg, skip)
    if arm in ('U60', 'U75'):
        mk = S6.Mahony(A1[0], withI=False); sck = S6.Score(ns, nseg, skip)
    Th_raw = A1[0].copy(); tlog, klog = [], []
    t0 = time.time()
    for k in range(N):
        w, a, th, V = W[k], A1[k], TH[k], Vk[k]
        Th_raw = Th_raw + w * DT
        s2 = had.update(a - Th_raw)
        tau = law.update(None if s2 is None else rep(s2, nl))
        wf = w - BIAS[k] if corr else w
        if not SMOKE:
            sc.add(k, mah.step(rep(wf, nl), rep(a, nl), tau) - rep(th, nl))
            if arm in ('U60', 'U75'):
                tk = np.full(ns, S6.kappa_H(H) * (S6.S0(V) / lamH) ** (1 / (2 * H + 1)))
                sck.add(k, mk.step(wf, a, tk) - th)
        sco.add(k, orc.step(rep(wf, nt), rep(a, nt), tau_or) - rep(th, nt))
        if k % 100 == 0:
            tlog.append(law.tau.copy()); klog.append(k)
        if k % 90000 == 0 and k:
            print(f'  {arm} t = {k * DT:.0f} s  ({time.time() - t0:.0f} s)', flush=True)
    res['runtime_s'] = time.time() - t0
    st, _ = sco.rms(); res['orc'] = {'steady': st.tolist()}
    tor, ror = S6.oracle_window(res, ns)
    if SMOKE:
        stg = np.array(st).reshape(3, nt, ns); inside = [bool(0 < np.argmin(stg[g].mean(1)) < nt - 1) for g in range(3)]
        print(f'  {arm}: oracle optimum strictly inside the grid per segment: {inside}; runtime {res["runtime_s"]:.0f} s', flush=True)
        return res
    res['oracle'] = {'tau': tor.tolist(), 'rms': ror.tolist()}
    klog = np.array(klog); T = np.array(tlog); st, tr = sc.rms()
    for j, nm in enumerate(names):
        cols = slice(j * ns, (j + 1) * ns)
        med = [np.median(T[(klog // nseg == g) & ((klog % nseg) >= skip)][:, cols], axis=0).tolist() for g in range(3)]
        res[nm] = {'steady': st[:, cols].tolist(), 'trans': tr[:, cols].tolist(), 'tau_med': med}
    if arm in ('U60', 'U75'):
        stk, trk = sck.rms()
        res['ULAW_K'] = {'steady': stk.tolist(), 'trans': trk.tolist(),
                         'tau_med': [[res['tau_pred'][g]] * ns for g in range(3)]}
    res['grids'] = {'tau': S6.TAU_GRID.tolist()}
    return res


def main():
    sha = hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest()
    print('step6d_uprime_mechanism.py  sha256', sha, flush=True)
    lam6b, _ = S6B.allan_lambda(SEEDS)
    out = {'sha256': sha, 'smoke': SMOKE, 'lambda_6b': lam6b.tolist(), 'tau_nom': TAU_NOM}
    for arm in ARMS:
        r = run(arm, lam6b); out[arm] = r
        if not SMOKE:
            print(f'{arm}: runtime {r["runtime_s"]:.0f} s', flush=True)
    fn = os.path.join(S6.OUT, 'step6d_smoke.json' if SMOKE else 'step6d_rows.json')
    if SMOKE:                                                # the smoke file keeps no oracle values
        for arm in ARMS:
            out[arm].pop('orc', None)
    json.dump(out, open(fn, 'w')); print('written', fn)


if __name__ == '__main__':
    main()
