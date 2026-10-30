"""STEP 6h, pre-registered (30 October 2026, after the step-6g diagnosis).

AMENDMENT 1 (after the first run, sha 83eabb57, which is void and recorded): (a) the schedule started at the window floor (9 ms),
where the Ki rule tied to the window explodes; the engine itself starts at 0.5 s, so the schedule is tau(t) = min(tau_law,
max(0.5 s, t)). (b) H6h.1 failed on its own terms in that run: D_H of the fused residual at 1 to 4 s lies below the white
extrapolation in both judged cells, so the accelerometer's wander is not a drift source at those scales, and a 40 s record cannot
resolve a lambda of 1e-7 there (a 4 % effect at 4 s on nine blocks). The drift the gyro fact misses is the gyro's own low-frequency
noise at throttle (Allan deviation at 3 to 10 s: 3x higher at low throttle than at idle), the limit 6b declared. So lambda is now
read FROM THE GYRO ALONE at the operating condition: Allan variance of the raw rate at l in {1, 2, 4} s over the record,
lambda_op = mean of sigma_A^2(l) l, floored at the idle fact. The fused reading is kept and reported, not used. (c) The gyro bias
is calibrated by the record mean (an arming calibration; optimistic for held frames, declared). Hypotheses unchanged except H6h.1,
which now reads: lambda_op / lambda_hat >= 2 at low throttle (idle is the reference and must give about 1). Same records as 6g: a pass is consistency with the
diagnosis, not an independent test; the cradle is.

Two changes to SYN-6B, both from the corpus, decided before this run:
 1. lambda read on the fused residual at the scale of seconds. z_raw = accel tilt - open-loop raw gyro integral (the truth cancels).
    Second differences of block means of z_raw at L in {32, 64, 128, 250, 500, 1000, 2000, 4000} ms. Model, white + drift:
    D_H(L) = 6 s2 / L + c lambda L, c calibrated here on a simulated Wiener process (printed). s2 = engine reading unchanged
    (min over 32, 64, 128 ms of D_H L / 6). lambda_eff = median over L in {1, 2, 4} s of (D_H(L) - 6 s2 / L) / (c L), floored at
    the gyro fact lambda_hat (a drift cannot be less than the gyro's own angle random walk). tau_law = sqrt(s2 / lambda_eff).
    Read over the whole record (declared: offline; the engine would carry it as a block EMA, which 40 s records cannot exercise).
 2. start-up schedule tau(t) = min(tau_law, max(DT/0.45, t)): the window cannot exceed the elapsed time (running mean until the
    law's window is reached). Replaces the 0.5 s start and the 0.4 s/s rate limit. Nothing else changes: same Mahony, same Ki rule.
Filters start on the first accel sample (the real start), gyro bias NOT calibrated (raw stream), 10 seeds of synthetic truth.
Comparators: ORACLE-S = fixed window from the first accel sample (same start as the law), grid and refinement as in 6g;
ORACLE-T = fixed window from the true angle (the clean window, as in 6g exploratory 3); KF-OPT from the first sample, tuned.
Scoring as in 6g (skip min(10 s, 30 %), steady RMS). Judged cells: S6-idle and S6-low (S5-half 8 s and S7-flight reported).
Hypotheses:
  H6h.1 the accelerometer's wander is a drift source the gyro fact misses: lambda_eff / lambda_hat >= 2 in both judged cells.
  H6h.2 the window lands: tau_law within x1.5 of the ORACLE-T window in both judged cells.
  H6h.3 the filter lands: regret vs ORACLE-S (same start) <= 10 %, seed median, in both judged cells.
  H6h.4 parity: mean RMS(law) / mean RMS(KF-OPT) <= 1.10 in both judged cells.
Reported: S5-half, S7-flight; the fitted D_H curve; the schedule's effect (law with the old 0.5 s start + rate limit, same lambda_eff).
Declared limits as in 6g, plus: lambda_eff read on 35 to 40 s records at 4 s scale uses about 9 blocks; its sampling error is large.
"""
import os, sys, json, numpy as np, hashlib, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import step6g_real_noise as G
G.LOG.clear(); P = G.P
LSCALES_MS = np.array([32, 64, 128, 250, 500, 1000, 2000, 4000]); LONG = [1000, 2000, 4000]

def hadamard(z, Lms):
    s = int(round(Lms / 1000 / G.DT)); mb = z[: len(z) // s * s].reshape(-1, s).mean(1); dd = mb[2:] - 2 * mb[1:-1] + mb[:-2]
    return np.mean(dd ** 2), len(dd)

def calibrate_c(n=400000, seed=1):
    r = np.random.default_rng(seed); wnr = np.cumsum(np.sqrt(G.DT) * r.standard_normal(n))     # Wiener, lambda = 1
    return {L: hadamard(wnr, L)[0] / (L / 1000) for L in LSCALES_MS}

def read_lambda(z, lam_hat, c):
    D = {L: hadamard(z, L)[0] for L in LSCALES_MS}
    s2 = min(D[L] * (L / 1000) / 6 for L in (32, 64, 128))
    lam_eff = np.median([(D[L] - 6 * s2 / (L / 1000)) / (c[L] * L / 1000) for L in LONG])
    return s2, max(lam_eff, lam_hat), D

class LawSched:
    def __init__(self, n, tau_law): self.tau_law = tau_law * np.ones(n); self.k = 0
    def update(self):
        t = self.k * G.DT; self.k += 1
        return np.minimum(self.tau_law, np.maximum(G.TAU0, t))

def run(name, f, a, b, lam_hat, c, old_start=False):
    tilt, w, g1k = G.load(name, f, a, b); w = w - w.mean(); N = len(tilt); ns = len(G.SEEDS); seeds = G.SEEDS
    def av(l):
        k = int(round(l / 1e-3)); x = g1k - g1k.mean(); mb = x[: len(x) // k * k].reshape(-1, k).mean(1); return 0.5 * np.mean(np.diff(mb) ** 2)
    lam_op = max(np.mean([av(l) * l for l in (1.0, 2.0, 4.0) if len(g1k) > 6 * int(l / 1e-3)]), lam_hat)
    TH = G.truth(N, seeds); W = np.diff(TH, axis=0) / G.DT + w[:, None]; A = TH[1:] + tilt[:, None]; TRU = TH[1:] + tilt.mean()
    Th = np.cumsum(w) * G.DT; s2, lam_fused, D = read_lambda(tilt - Th, lam_hat, c); lam_eff = lam_op; tau_law = float(np.sqrt(s2 / lam_eff))
    skip = int(round(min(10.0, 0.3 * N * G.DT) / G.DT))
    mah = G.Mahony(A[0].copy()); sch = LawSched(ns, tau_law)
    if old_start: law_old = G.LawTau(ns, np.full(ns, lam_eff)); law_old.tau_t = np.full(ns, min(tau_law, G.TAU_MAX))
    nt = len(G.TAU_GRID); orcS = G.Mahony(np.tile(A[0], nt)); orcT = G.Mahony(np.tile(TRU[0], nt)); tau_o = np.repeat(G.TAU_GRID, ns)
    Rg = tilt.var() * G.R_MULT; nR, nQ = len(Rg), len(G.QG_MULT); nk = nR * nQ
    kf = G.Kalman2(np.tile(A[0], nk), np.repeat(Rg, nQ * ns), np.repeat(np.repeat(G.QG_MULT * lam_hat, ns)[None, :], nR, 0).reshape(-1)); R_kf = np.repeat(Rg, nQ * ns)
    ss = dict(law=0.0, orcS=0.0, orcT=0.0, kf=0.0, old=0.0); n = 0
    for k in range(N):
        wk, ak, tk = W[k], A[k], TRU[k]
        if old_start:
            s2v = None; tau_old = law_old.update(np.full(ns, s2))   # old engine start: 0.5 s start, 0.4 s/s ramp toward the same window
        e = dict(law=mah.step(wk, ak, sch.update()) - tk,
                 orcS=orcS.step(np.tile(wk, nt), np.tile(ak, nt), tau_o) - np.tile(tk, nt),
                 orcT=orcT.step(np.tile(wk, nt), np.tile(ak, nt), tau_o) - np.tile(tk, nt),
                 kf=kf.step(np.tile(wk, nk), np.tile(ak, nk), R_kf) - np.tile(tk, nk))
        if k >= skip:
            for key in e: ss[key] = ss[key] + e[key] ** 2
            n += 1
    rms = {key: np.sqrt(v / n) / G.D2R for key, v in ss.items() if key != 'old'}
    lg = np.log(G.TAU_GRID)
    oS = np.array([G.refine(lg, rms['orcS'].reshape(nt, ns)[:, j]) for j in range(ns)]); oT = np.array([G.refine(lg, rms['orcT'].reshape(nt, ns)[:, j]) for j in range(ns)])
    kfb = rms['kf'].reshape(-1, ns).min(0)
    return dict(name=name, record_s=N * G.DT, s2=s2, lam_eff=lam_eff, lam_fused=lam_fused, lam_ratio=lam_eff / lam_hat, tau_law=tau_law, D={int(k): float(v) for k, v in D.items()},
                law_rms=rms['law'].tolist(), orcS_tau=oS[:, 0].tolist(), orcS_rms=oS[:, 1].tolist(), orcT_tau=oT[:, 0].tolist(), orcT_rms=oT[:, 1].tolist(), kf_rms=kfb.tolist())

if __name__ == '__main__':
    sha = hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest(); P('step6h_fused_lambda.py sha256', sha, '| seeds', G.SEEDS)
    t0 = time.time(); c = calibrate_c(); P('Wiener calibration c(L) = D_H / (lambda L):', {int(k): round(float(v), 4) for k, v in c.items()})
    _, _, g_idle = G.load(*G.SEGS[0][:4]); lam_hat, _ = G.allan_lambda(g_idle); P(f'gyro fact lambda_hat = {lam_hat:.3e} rad^2/s (idle, as in 6g)')
    res = {}; H = [True] * 4
    P(f'\n{"segment":10s} {"rec":>5s} {"lam_op/hat":>11s} {"tau law":>8s} {"orc-T tau":>9s} {"orc-S tau":>9s} {"RMS law":>8s} {"RMS orc-S":>9s} {"RMS orc-T":>9s} {"RMS KF":>7s} {"regret(S)":>9s} {"law/KF":>6s}')
    for name, f, a, b, role in G.SEGS:
        r = run(name, f, a, b, lam_hat, c); res[name] = r
        law = np.array(r['law_rms']); oS = np.array(r['orcS_rms']); kf = np.array(r['kf_rms'])
        regret = np.median(law / oS - 1); tT = np.median(r['orcT_tau']); tS = np.median(r['orcS_tau']); ratio = law.mean() / kf.mean(); win = r['tau_law'] / tT
        judged = name in ('S6-idle', 'S6-low')
        P(f'{name:10s} {r["record_s"]:5.1f} {r["lam_ratio"]:11.1f} {r["tau_law"]:8.1f} {tT:9.1f} {tS:9.1f} {law.mean():8.3f} {oS.mean():9.3f} {np.mean(r["orcT_rms"]):9.3f} {kf.mean():7.3f} {regret*100:8.1f}% {ratio:6.3f}   [{"judged" if judged else "reported"}]')
        if judged:
            H[0] &= (r['lam_ratio'] >= 2) if name == 'S6-low' else True; H[1] &= 1 / 1.5 <= win <= 1.5; H[2] &= regret <= 0.10; H[3] &= ratio <= 1.10
    P(f'\nH6h.1 (amended) the gyro drift at throttle exceeds the idle fact (lambda_op/lambda_hat >= 2 at low throttle): {"PASS" if H[0] else "FAIL"}')
    P(f'H6h.2 window within x1.5 of the clean oracle: {"PASS" if H[1] else "FAIL"}')
    P(f'H6h.3 filter regret vs same-start oracle <= 10 %: {"PASS" if H[2] else "FAIL"}')
    P(f'H6h.4 parity with tuned Kalman <= 1.10: {"PASS" if H[3] else "FAIL"}')
    P('\nReported: D_H(L) of z_raw (rad^2) by scale, and the white-only prediction 6 s2 / L')
    for name, r in res.items():
        P(f'  {name:10s} ' + ' '.join(f'{L}:{r["D"][L]:.2e}/{6*r["s2"]/(L/1000):.2e}' for L in LSCALES_MS))
    for r in res.values(): r['D'] = {str(k): v for k, v in r['D'].items()}
    P(f'runtime {time.time()-t0:.0f} s'); json.dump(res, open(os.path.join(G.OUT, 'step6h_rows.json'), 'w'), indent=1, default=float)
    open(os.path.join(G.OUT, 'step6h_results.txt'), 'w').write('\n'.join(G.LOG) + '\n')
