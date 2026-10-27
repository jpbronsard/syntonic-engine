"""STEP 6e, pre-registered. Does the optimal window's exponent follow Theorem U' (1/(2H+1)) as the gyro drift's memory
index H changes, with enough statistical power to tell it from the square root?

Origin. Step 6d (sha c6e49019) put the attitude estimator inside U''s class (gyro drift = fractional Brownian motion of
index H, no bias, Mahony without bias integral = a first-order window on z) and found optimal exponents 0.484 (H = 0.6)
and 0.441 (H = 0.75) against 0.455 and 0.400 predicted. An independent check showed the test had little power: ten
seeds, per-seed exponents with a standard deviation of about 0.14, 95 % intervals containing 0.5; and its square-root
control failed mostly on its level (a drift fact miscalibrated on fractional drift), not on its exponent. Step 6e is the
powered version.

Design. The generator, filters and constants of step 6d (imported from step6d_uprime_mechanism.py unchanged: same
accelerometer noise, attitude truth and fractional drift construction, lambda_H set so that U''s window at nominal
vibration is 1.077 s for every H), scenario W (V = 1, 256, 1/16; 180 s segments; steady part the last 120 s), and:
  seeds   1000..1099 (100 new seeds), the same seeds for every H (paired design: the accelerometer noise and truth of a
          seed are identical across H; only the drift differs).
  H       0.5, 0.6, 0.75, 0.9 (arms U50, U60, U75, U90); predicted exponents p_H = 1/(2H+1) = 0.500, 0.455, 0.400, 0.357.
  ORACLE  Mahony without bias integral, fixed windows, 36 values 0.05 to 100 s (log-spaced; wider than step 6d's grid,
          whose 50 s edge censored some optima).
  POW-U   tau = tau_nom (s2 / S0_nom) ^ p_H, tau_nom = 1.077 s, S0_nom = S0(V = 1): U''s exponent, calibrated at nominal.
  POW-S   tau = tau_nom (s2 / S0_nom) ^ 0.5: the square root, calibrated at the same nominal point. POW-U and POW-S differ
          only by their exponent (identical at H = 0.5).
  ULAW    U' with H and lambda_H measured once at rest (step 6d's procedure: 300 s static record, block means at 0.128 to
          8.192 s) and s2 measured live: tau = kappa_H (s2 / lambda_H) ^ (1 / (2H + 1)).
  ULAW-K  U' with the true H, lambda_H and S0 (clairvoyant; no rate limit, as step 6's SYN-K).
  s2 is step 6b's live noise estimate; POW-U, POW-S and ULAW use 6b's machinery (floor, ceiling 30 s, rate limit, warm-up,
  start). All filters are Mahony without bias integral.

Population optimum (new here; step 6d used the seed median of per-seed optima). For each arm and segment, the steady
mean-square error of every fixed window averaged over the 100 seeds, J_g(tau); tau_opt(g) its minimiser (parabola in
log tau through the grid minimum and its neighbours), J_opt(g) the parabola's minimum. Exponent upward
e_H = log(tau_opt(1) / tau_opt(0)) / log(S0(256) / S0(1)); downward the same with segment 2. Regret of an estimator in a
segment: sqrt(seed-mean of its squared RMS / J_opt) - 1. Uncertainty: bootstrap over seeds (2000 draws, the same seed
draw for all arms, so the pairing is kept), 95 % percentile intervals.

Hypotheses and decision rules (judge: evaluate_step6e.py, written with this script, before any run):
  H6e.1 the exponent changes with H as U' says: the least-squares slope b of the four upward exponents e_H against p_H
        lies in [0.7, 1.3], and its 95 % interval excludes 0 (the square root everywhere predicts b = 0).
  H6e.2 the exponent's value: |e_H - p_H| <= 0.05 for each of the four H.
  H6e.3 power: the 95 % intervals of e_0.75 and e_0.9 exclude 0.5.
  H6e.4 the exponent matters in performance: at V = 256, J(POW-U) / J(POW-S) < 1 with its 95 % interval entirely below 1,
        for H = 0.75 and H = 0.9.
  H6e.5 the practical law: ULAW's regret <= 10 % and its window (seed median) within x1.5 of tau_opt, in all three segments,
        for all four H.
  Reported, not judged: downward exponents and their slope; ULAW-K; POW-U and POW-S regrets; H and lambda_H measured at
  rest; per-seed spreads; whether any population optimum lies at a grid edge.
Reading rules, declared: H6e.1 and H6e.3 together say the exponent follows H as U' predicts and reject the square root
everywhere, in this generator. H6e.1 holding with H6e.2 failing means U''s dependence on H is right with a systematic
offset (of the generator or the oracle construction), to be reported as such. H6e.4 says the exponent changes the error,
not only the window. H6e.5 says the engine can use it with facts measured once at rest.
Declared limits: one simulated axis, representative parameters, the generator family of steps 6 to 6d (not an
independent test); fractional Brownian drift is a stand-in for real gyro low-frequency noise (no angle random walk, no
bias instability, no rate random walk) and H is identified on the same pure power law; no bias channel; POW-U and POW-S
are calibrated at the true nominal level. A failure is published as such.

Smoke policy: one short run of U90 only (4 seeds, 90 s segments, 30 s skipped) that prints only the fractional-noise
variance check, H measured at rest, whether each population optimum lies strictly inside the grid (yes or no), and the
runtime; never a window, an exponent or a score. Any change after it is written here as an amendment before the full run.
"""
import hashlib, json, os, sys, time
from concurrent.futures import ProcessPoolExecutor
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import step6_quad_attitude_estimator as S6
import step6b_gyro_fact as S6B
import step6d_uprime_mechanism as S6D

SMOKE = bool(os.environ.get('STEP6E_SMOKE'))
DT = S6.DT
H_ARMS = {'U50': 0.5, 'U60': 0.6, 'U75': 0.75, 'U90': 0.9}
SEEDS = list(range(1000, 1100)); BATCH = 25
GRID = np.logspace(np.log10(0.05), np.log10(100.0), 36)
SEG_S, SKIP_S = 180.0, 60.0
if SMOKE:
    H_ARMS = {'U90': 0.9}; SEEDS = [1000, 1001, 1002, 1003]; BATCH = 2; SEG_S, SKIP_S = 90.0, 30.0
S0_NOM = S6.S0(1.0); TAU_NOM = S6D.TAU_NOM


def run_batch(args):
    arm, seeds = args
    S6D.SEG_S = SEG_S
    H = H_ARMS[arm]; lamH = S6D.lam_for(H); pH = 1.0 / (2 * H + 1)
    d = S6D.generate(seeds, fbm=(H, lamH))
    TH, W, A1, Vk, N, nseg = d['TH'], d['W'], d['A1'], d['Vk'], d['N'], d['nseg']
    ns = len(seeds); nt = len(GRID); skip = int(round(SKIP_S / DT)); rep = lambda x, m: np.tile(x, m)
    Hh, Lh, _ = S6D.rest_facts_fbm(seeds, H, lamH)
    out = {'seeds': seeds, 'H_rest': Hh.tolist(), 'lamH_rest': (Lh / lamH).tolist()}
    if SMOKE:
        out['fgn_var'] = float(np.var(S6D.fgn(2 ** 16, H, np.random.default_rng(1))))
    # laws: POW-U, POW-S, ULAW (columns [law][seed])
    law = S6D.Law(np.concatenate([np.full(ns, TAU_NOM), np.full(ns, TAU_NOM), S6.kappa_H(Hh)]),
                  np.concatenate([np.full(ns, S0_NOM), np.full(ns, S0_NOM), Lh]),
                  np.concatenate([np.full(ns, pH), np.full(ns, 0.5), 1.0 / (2 * Hh + 1)]))
    nl = 3
    mah = S6.Mahony(rep(A1[0], nl), withI=False); had = S6B.HadNoise(ns)
    mk = S6.Mahony(A1[0], withI=False)
    orc = S6.Mahony(rep(A1[0], nt), withI=False); tau_or = np.repeat(GRID, ns)
    sc = S6.Score(nl * ns, nseg, skip); sck = S6.Score(ns, nseg, skip); sco = S6.Score(nt * ns, nseg, skip)
    Th_raw = A1[0].copy(); tlog, klog = [], []
    for k in range(N):
        w, a, th, V = W[k], A1[k], TH[k], Vk[k]
        Th_raw = Th_raw + w * DT
        s2 = had.update(a - Th_raw)
        tau = law.update(None if s2 is None else rep(s2, nl))
        if not SMOKE:
            sc.add(k, mah.step(rep(w, nl), rep(a, nl), tau) - rep(th, nl))
            tk = np.full(ns, S6.kappa_H(H) * (S6.S0(V) / lamH) ** (1 / (2 * H + 1)))
            sck.add(k, mk.step(w, a, tk) - th)
        sco.add(k, orc.step(rep(w, nt), rep(a, nt), tau_or) - rep(th, nt))
        if k % 100 == 0:
            tlog.append(law.tau.copy()); klog.append(k)
    so, _ = sco.rms(); out['orc'] = np.array(so).reshape(3, nt, ns).tolist()
    if not SMOKE:
        st, tr = sc.rms(); stk, trk = sck.rms()
        klog = np.array(klog); T = np.array(tlog)
        for j, nm in enumerate(('POW_U', 'POW_S', 'ULAW')):
            cols = slice(j * ns, (j + 1) * ns)
            med = [np.median(T[(klog // nseg == g) & ((klog % nseg) >= skip)][:, cols], axis=0).tolist() for g in range(3)]
            out[nm] = {'steady': st[:, cols].tolist(), 'trans': tr[:, cols].tolist(), 'tau_med': med}
        out['ULAW_K'] = {'steady': stk.tolist(), 'trans': trk.tolist()}
    return arm, out


def main():
    sha = hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest()
    print('step6e_uprime_powered.py  sha256', sha, flush=True)
    tasks = [(arm, SEEDS[i:i + BATCH]) for arm in H_ARMS for i in range(0, len(SEEDS), BATCH)]
    t0 = time.time(); parts = {arm: [] for arm in H_ARMS}
    with ProcessPoolExecutor(max_workers=2) as ex:
        for arm, o in ex.map(run_batch, tasks):
            parts[arm].append(o); print(f'  {arm} seeds {o["seeds"][0]}..{o["seeds"][-1]} done ({time.time() - t0:.0f} s)', flush=True)
    res = {'sha256': sha, 'smoke': SMOKE, 'grid': GRID.tolist(), 'V': S6D.V_W, 'seg_s': SEG_S, 'skip_s': SKIP_S,
           'tau_nom': TAU_NOM, 'runtime_s': time.time() - t0}
    for arm, H in H_ARMS.items():
        ps = sorted(parts[arm], key=lambda o: o['seeds'][0])
        r = {'H': H, 'p_H': 1 / (2 * H + 1), 'lamH': S6D.lam_for(H), 'seeds': sum((o['seeds'] for o in ps), []),
             'H_rest': sum((o['H_rest'] for o in ps), []), 'lamH_rest_ratio': sum((o['lamH_rest'] for o in ps), []),
             'orc': np.concatenate([np.array(o['orc']) for o in ps], axis=2).tolist()}
        if SMOKE:
            orc = np.array(r['orc']); J = (orc ** 2).mean(2)
            inside = [bool(0 < int(np.argmin(J[g])) < len(GRID) - 1) for g in range(3)]
            print(f'{arm}: fGn variance check {ps[0]["fgn_var"]:.3f}; H at rest {np.round(r["H_rest"], 3).tolist()} (true {H});'
                  f' population optimum strictly inside the grid per segment: {inside}', flush=True)
            r.pop('orc')
        else:
            for nm in ('POW_U', 'POW_S', 'ULAW', 'ULAW_K'):
                r[nm] = {k: np.concatenate([np.array(o[nm][k]) for o in ps], axis=1).tolist() for k in ps[0][nm]}
        res[arm] = r
    print(f'runtime {res["runtime_s"]:.0f} s', flush=True)
    fn = os.path.join(S6.OUT, 'step6e_smoke.json' if SMOKE else 'step6e_rows.json')
    json.dump(res, open(fn, 'w')); print('written', fn)


if __name__ == '__main__':
    main()
