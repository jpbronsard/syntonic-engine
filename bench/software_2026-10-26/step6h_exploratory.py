"""STEP 6h EXPLORATORY (after 6h-A1, sha 12751351, which passed H6h.1 and failed H6h.2 to H6h.4). Not judged.
Mechanism to test: the Ki rule slaved to the window (Ki = Kp^2 / 4 zeta^2) poisons the bias estimate during the start (error 5 deg, Ki = 0.3 at
tau = 0.5 s), then cannot repair it once the window is 25 to 40 s (bias-loop time constant 4 zeta^2 tau = 5 to 9 min).
Arms, all with lambda_op from the gyro at the operating condition, the 0.5 s schedule, bias calibrated by the record mean:
  FROZEN  Ki = 0 while t < tau_law (schedule phase), then the Ki rule at the law's window;
  NOKI    Ki = 0 always (bias calibrated at arming, no bias loop);
  SLAVED  the 6h-A1 filter, for reference.
Same comparators and scoring as 6h."""
import os, sys, json, numpy as np, hashlib
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import step6g_real_noise as G, step6h_fused_lambda as H
G.LOG.clear(); P = G.P
class MahonyKi(G.Mahony):
    def step(self, w, a, tau, ki_on):
        kp = 1.0 / tau; wb = w - self.b; self.th = self.th + wb * G.DT; self.Th = self.Th + wb * G.DT
        e = a - self.th; self.th = self.th + np.minimum(G.DT * kp, 1.0) * e
        ki = kp * kp / (4 * G.ZETA ** 2) * ki_on; self.b = np.clip(self.b - ki * e * G.DT, -0.15, 0.15); return self.th
def run(name, f, a, b, lam_hat):
    tilt, w, g1k = G.load(name, f, a, b); w = w - w.mean(); N = len(tilt); ns = len(G.SEEDS)
    def av(l):
        k = int(round(l / 1e-3)); x = g1k - g1k.mean(); mb = x[: len(x) // k * k].reshape(-1, k).mean(1); return 0.5 * np.mean(np.diff(mb) ** 2)
    lam_op = max(np.mean([av(l) * l for l in (1.0, 2.0, 4.0) if len(g1k) > 6 * int(l / 1e-3)]), lam_hat)
    TH = G.truth(N, G.SEEDS); W = np.diff(TH, axis=0) / G.DT + w[:, None]; A = TH[1:] + tilt[:, None]; TRU = TH[1:] + tilt.mean()
    Th = np.cumsum(w) * G.DT; s2, _, _ = H.read_lambda(tilt - Th, lam_hat, {L: 1.0 for L in H.LSCALES_MS}); tau_law = float(np.sqrt(s2 / lam_op))
    skip = int(round(min(10.0, 0.3 * N * G.DT) / G.DT))
    arms = {k: MahonyKi(A[0].copy()) for k in ('FROZEN', 'NOKI', 'SLAVED')}; sch = H.LawSched(ns, tau_law)
    nt = len(G.TAU_GRID); orcS = G.Mahony(np.tile(A[0], nt)); orcT = G.Mahony(np.tile(TRU[0], nt)); tau_o = np.repeat(G.TAU_GRID, ns)
    ss = {k: 0.0 for k in list(arms) + ['orcS', 'orcT']}; n = 0
    for k in range(N):
        wk, ak, tk = W[k], A[k], TRU[k]; tau = sch.update(); t = k * G.DT
        e = {'FROZEN': arms['FROZEN'].step(wk, ak, tau, 1.0 if t >= tau_law else 0.0) - tk, 'NOKI': arms['NOKI'].step(wk, ak, tau, 0.0) - tk,
             'SLAVED': arms['SLAVED'].step(wk, ak, tau, 1.0) - tk,
             'orcS': orcS.step(np.tile(wk, nt), np.tile(ak, nt), tau_o) - np.tile(tk, nt), 'orcT': orcT.step(np.tile(wk, nt), np.tile(ak, nt), tau_o) - np.tile(tk, nt)}
        if k >= skip:
            for key in e: ss[key] = ss[key] + e[key] ** 2
            n += 1
    rms = {key: np.sqrt(v / n) / G.D2R for key, v in ss.items()}; lg = np.log(G.TAU_GRID)
    oS = np.array([G.refine(lg, rms['orcS'].reshape(nt, ns)[:, j]) for j in range(ns)]); oT = np.array([G.refine(lg, rms['orcT'].reshape(nt, ns)[:, j]) for j in range(ns)])
    return dict(tau_law=tau_law, lam_ratio=lam_op / lam_hat, orcS_tau=np.median(oS[:, 0]), orcS_rms=oS[:, 1], orcT_tau=np.median(oT[:, 0]), orcT_rms=oT[:, 1], **{k: rms[k] for k in arms})
if __name__ == '__main__':
    P('step6h_exploratory.py sha256', hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest()[:12])
    _, _, g_idle = G.load(*G.SEGS[0][:4]); lam_hat, _ = G.allan_lambda(g_idle - g_idle.mean())
    P(f'{"segment":10s} {"tau law":>8s} {"orc-S tau":>9s} {"orc-T tau":>9s} | {"SLAVED":>7s} {"FROZEN":>7s} {"NOKI":>7s} | {"orc-S":>7s} {"orc-T":>7s} | regret FROZEN / NOKI vs orc-S')
    for name, f, a, b, role in G.SEGS[:3]:
        r = run(name, f, a, b, lam_hat)
        rf = np.median(r['FROZEN'] / r['orcS_rms'] - 1); rn = np.median(r['NOKI'] / r['orcS_rms'] - 1)
        P(f'{name:10s} {r["tau_law"]:8.1f} {r["orcS_tau"]:9.1f} {r["orcT_tau"]:9.1f} | {r["SLAVED"].mean():7.3f} {r["FROZEN"].mean():7.3f} {r["NOKI"].mean():7.3f} | {r["orcS_rms"].mean():7.3f} {r["orcT_rms"].mean():7.3f} | {rf*100:+6.1f}% / {rn*100:+6.1f}%')
    P('  RMS in deg, steady part; orc-S = fixed window from the first accel sample, orc-T = from the true angle')
    open(os.path.join(G.OUT, 'step6h_exploratory.txt'), 'w').write('\n'.join(G.LOG) + '\n')
