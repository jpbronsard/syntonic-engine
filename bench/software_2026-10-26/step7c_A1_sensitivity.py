"""STEP 7c A1, exploratory sensitivity: flight window without the departure second (19-23 s, 20-23 s), and the
V5.1 matrix window (T_sym) against the 3-axis Kalman window on the flight matrices."""
import os, sys, numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import step7c_rate3_th44 as S
from scipy.linalg import inv
out = []; P_ = lambda *a: out.append(' '.join(str(x) for x in a))
d7 = pd.read_csv(S.UP + 'crash7.csv'); t7 = (d7.t_us.values - d7.t_us.values[0]) * 1e-6
rng = np.random.default_rng(79)
for a, b in ((19.0, 24.0), (19.0, 23.0), (20.0, 23.0)):
    fn = lambda tt: (tt >= a) & (tt < b); m = fn(t7)
    runs = S.segs(d7, t7, fn); tau = [float(np.median(d7[f'tau_{x}'].values[m])) for x in 'RPY']
    D = {l: S.dseq(runs, l) for l in S.LS + (128,)}
    def est(D):
        Sn = np.mean([D[l].T @ D[l] / len(D[l]) * (l * 1e-3) / 6 for l in S.LS], 0)
        return Sn, D[128].T @ D[128] / len(D[128]) * 3 / (2 * 0.128)
    Sn, Lam = est(D); g, Td, Tk = S.d2r(Lam, Sn); Ts = S.gmean(inv(Lam), Sn)
    gs = S.trP(Ts, Lam, Sn) / S.trP(Tk, Lam, Sn) - 1
    bg = []; bc = []; bf = []; bs = []
    for _ in range(S.NBOOT):
        Snb, Lamb = est(S.boot(D, rng))
        if np.all(np.linalg.eigvalsh(Lamb) > 0):
            bg.append(S.d2r(Lamb, Snb)[0]); bc.append(S.corr(Lamb)[0, 1])
            cf, ca, fac, sh, _ = S.th44(tau, Lamb, 4e-3); bf.append(fac); bs.append(sh.max())
    cf, ca, fac, sh, _ = S.th44(tau, Lam, 4e-3)
    q = lambda x: f'{np.quantile(x,.05):+.3f}..{np.quantile(x,.95):+.3f}'
    P_(f'{a:.0f}-{b:.0f} s: n128 {len(D[128])}, gyro mean dps {np.round(np.concatenate(runs).mean(0)*180/np.pi,1).tolist()},'
       f' corr(Lam) RP {S.corr(Lam)[0,1]:+.2f} ({q(bc)}) RY {S.corr(Lam)[0,2]:+.2f} PY {S.corr(Lam)[1,2]:+.2f};'
       f' D2r {100*g:.2f} % ({q(100*np.array(bg))}); T_sym vs Kalman {100*gs:.3f} %;'
       f' D4 (4 ms) factor {100*fac:+.2f} % ({q(100*np.array(bf))}), v_c R/P/Y {np.round(sh,3).tolist()} (largest share {q(bs)})')
txt = '\n'.join(out); print(txt)
open(os.path.join(HERE, 'results', 'step7c_A1_sensitivity.txt'), 'w').write(txt + '\n')
