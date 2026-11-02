"""Diagnostic, no truth used: the published law's internal window (tau actually used, and the law's target taul =
sqrt(s2/lambda)) and the live s2 over a flight, per circumstance, on dev seed 300. Prints medians per regime/phase."""
import sys, numpy as np
import common as C
F = C.F
circ = sys.argv[1] if len(sys.argv) > 1 else 'C3'
seed = 300
gn, tn, reg = C.build_noise_c(seed, circ); N = len(gn)
est = F.Est(2, np.array([2, 2]), cfac=np.array([1.0, 2.8]))
lam_of = {r: C.FX['lam'][r] for r in ('idle', 'low', 'half')}
bc = 0.0; gs = 0.0; rows = []
for k in range(N):
    g = np.full(2, gn[k] - bc); a = np.full(2, tn[k // 4]) if k % 4 == 0 else None
    est.step(g, a, lam_of[reg[k]], k * C.DT)
    if k < C.KPRE:
        gs += gn[k]
        if k == C.KPRE - 1: bc = gs / C.KPRE
    if k % 100 == 0:
        s2 = est.had.DH * (G_ := C.G.HSCALES * C.G.DT)[:, None] / 6.0
        rows.append((k * C.DT - C.T_PRE, reg[k], est.tau[0], est.tau[1], est.taul[0], s2[:, 0].min()))
t = np.array([r[0] for r in rows]); rg = np.array([r[1] for r in rows]); tau1 = np.array([r[2] for r in rows]); tau28 = np.array([r[3] for r in rows])
taul = np.array([r[4] for r in rows]); s2 = np.array([r[5] for r in rows])
print(circ, 'seed', seed)
for name in ('idle', 'low', 'half'):
    m = (rg == name) & (t >= 0)
    if m.any():
        print(f"  {name:5s} flight: tau(c=1) median {np.median(tau1[m]):6.2f} [p10 {np.percentile(tau1[m], 10):6.2f}, p90 {np.percentile(tau1[m], 90):6.2f}]  "
              f"tau(c=2.8) median {np.median(tau28[m]):6.2f}  target sqrt(s2/lam) median {np.median(taul[m]):7.2f}  s2 median {np.median(s2[m]):.2e} [p10 {np.percentile(s2[m], 10):.2e}, p90 {np.percentile(s2[m], 90):.2e}]")
for t0 in range(-5, 60, 5):
    m = (t >= t0) & (t < t0 + 5)
    print(f"  t {t0:3d}..{t0 + 5:3d} {rg[m][0]:5s} tau1 {np.median(tau1[m]):6.2f}  tau2.8 {np.median(tau28[m]):6.2f}  taul {np.median(taul[m]):7.2f}")
