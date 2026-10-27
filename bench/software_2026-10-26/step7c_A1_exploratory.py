"""STEP 7c, amendment A1, EXPLORATORY (written after the pre-registered run, sha d1f43d53, found Lam not positive
definite in all three windows). Cause seen in that output: the short-scale noise is not white (rotor vibration is
narrowband), so subtracting 6 Sn / l over-corrects at 32 ms and drives diagonal terms negative.
A1: Lam from the 128 ms scale alone, C(128) * 3 / (2 l), no noise subtraction; the share of the white-noise term
6 Sn / l in C(128) is printed per axis as the bias bound. 256 ms reported beside it. D2r and D4 then exactly as
written in step7c_rate3_th44.py. These numbers are exploratory, not pre-registered decisions."""
import json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import step7c_rate3_th44 as S
import pandas as pd

out = []; P_ = lambda *a: out.append(' '.join(str(x) for x in a))
rng = np.random.default_rng(78)
W = {}
df = pd.read_csv(S.UP + 'sess6.csv'); t = (df.t_us.values - df.t_us.values[0]) * 1e-6
bins = np.floor(t / 5).astype(int); med = pd.Series(df.amps_A.values).groupby(bins).transform('median').values
W['W1 bench'] = (df, t, lambda tt: med > 0.8)
d5 = pd.read_csv(S.UP + 'syntonic_00000000_000000_00005.csv'); t5 = (d5.t_us.values - d5.t_us.values[0]) * 1e-6
W['W2 thrust'] = (d5, t5, lambda tt: (tt >= 15.0) & (tt < 23.0))
d7 = pd.read_csv(S.UP + 'crash7.csv'); t7 = (d7.t_us.values - d7.t_us.values[0]) * 1e-6
W['W3 flight'] = (d7, t7, lambda tt: (tt >= 19.0) & (tt < 24.0))

def est(D):
    Sn = np.mean([D[l].T @ D[l] / len(D[l]) * (l * 1e-3) / 6 for l in S.LS], 0)
    C = D[128].T @ D[128] / len(D[128]); Lam = C * 3 / (2 * 0.128)
    return Sn, Lam, C

res = {}
for wn, (d, tt, fn) in W.items():
    runs = S.segs(d, tt, fn); m = fn(tt)
    tau = [float(np.median(d[f'tau_{a}'].values[m])) for a in 'RPY']
    D = {l: S.dseq(runs, l) for l in S.LS + (128, 256)}
    Sn, Lam, C = est(D)
    C256 = D[256].T @ D[256] / len(D[256]); L256 = C256 * 3 / (2 * 0.256)
    share = 6 * np.diag(Sn) / 0.128 / np.diag(C)
    pd_ok = bool(np.all(np.linalg.eigvalsh(Lam) > 0))
    P_(f'\n== {wn}: n(128 ms) = {len(D[128])}, n(256 ms) = {len(D[256])}; tau R/P/Y {np.round(np.array(tau)*1e3,1).tolist()} ms')
    P_(f'  white-noise share of C(128) per axis: {np.round(share,3).tolist()}; Lam PD: {pd_ok}')
    P_('  corr(Lam, 128 ms) =', np.array2string(S.corr(Lam), precision=3).replace('\n', ''),
       ' | corr(256 ms) =', np.array2string(S.corr(L256), precision=3).replace('\n', ''))
    P_(f'  Lam diag 128 ms {np.array2string(np.diag(Lam), precision=4)}, 256 ms {np.array2string(np.diag(L256), precision=4)}')
    r = dict(n128=len(D[128]), tau=tau, Lam=Lam.tolist(), corr_Lam=S.corr(Lam).tolist(), corr_Lam256=S.corr(L256).tolist(),
             noise_share=share.tolist())
    if not pd_ok: res[wn] = r; continue
    comm = np.linalg.norm(Lam @ Sn - Sn @ Lam) / (np.linalg.norm(Lam) * np.linalg.norm(Sn))
    g, Td, Tk = S.d2r(Lam, Sn)
    bg = []; bt = {T: [] for T in S.TENV}; bc = []
    for _ in range(S.NBOOT):
        Db = S.boot(D, rng); Snb, Lamb, _ = est(Db)
        if np.all(np.linalg.eigvalsh(Lamb) > 0):
            bg.append(S.d2r(Lamb, Snb)[0]); bc.append(S.corr(Lamb)[np.triu_indices(3, 1)])
            for T in S.TENV:
                cf, ca, fac, sh, _ = S.th44(tau, Lamb, T); bt[T].append([fac, sh.max()])
    bg = np.array(bg); bc = np.array(bc)
    P_(f'  corr(Lam) off-diagonals RP, RY, PY 90 % boot: ' + ', '.join(f'{np.quantile(bc[:,k],0.05):+.2f}..{np.quantile(bc[:,k],0.95):+.2f}' for k in range(3)))
    P_(f'  commutator {comm:.4f}; per-axis windows {np.round(np.diag(Td)*1e3,2).tolist()} ms')
    P_(f'  D2r: per-axis / 3-axis Kalman - 1 = {100*g:.3f} % (90 % {100*np.quantile(bg,0.05):.3f} .. {100*np.quantile(bg,0.95):.3f} %, {len(bg)} draws)')
    r.update(comm=comm, D2r=g, D2r_ci=[np.quantile(bg, .05), np.quantile(bg, .95)], D4={})
    for T in S.TENV:
        cf, ca, fac, sh, sha_ = S.th44(tau, Lam, T); b = np.array(bt[T])
        kind = 'single-axis' if sh.max() >= 0.90 else 'mixed'
        P_(f'  D4 T_env {T*1e3:.0f} ms: factor {100*fac:+.3f} % (90 % {100*np.quantile(b[:,0],.05):+.3f} .. {100*np.quantile(b[:,0],.95):+.3f}),'
           f' v_c shares R/P/Y {np.round(sh,3).tolist()} -> {kind} (largest share 90 % {np.quantile(b[:,1],.05):.3f}..{np.quantile(b[:,1],.95):.3f}); per-axis v_c {np.round(sha_,3).tolist()}')
        r['D4'][str(T)] = dict(factor=fac, shares=sh.tolist(), kind=kind, factor_ci=[np.quantile(b[:, 0], .05), np.quantile(b[:, 0], .95)])
    res[wn] = r
txt = '\n'.join(out); print(txt)
open(os.path.join(HERE, 'results', 'step7c_A1_results.txt'), 'w').write(txt + '\n')
json.dump(res, open(os.path.join(HERE, 'results', 'step7c_A1_rows.json'), 'w'), default=float)
