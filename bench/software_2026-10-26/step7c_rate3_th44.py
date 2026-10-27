"""STEP 7c, rules written before any covariance of these windows was computed. Roll, pitch AND yaw: are the three
rate channels coupled enough to need the matrix law, and is the TH44 collapse direction single-axis or mixed?

Why a second check. Step 7a answered for the attitude estimator, roll and pitch only (yaw has no absolute sensor in the
SD log: no magnetometer values). The briefing's open diagnostic is about the 3x3 covariances of the rate channels,
where the M5/M6 collapses happened. The M5/M6 telemetry (analysis/maiden6lucky.csv, about 50 Hz) logs only per-axis
sig2 and lam, no cross-axis terms, so the certificate cannot be computed on those flights. The raw 1 kHz gyro of the
SD logger exists for three sessions; this script uses them.

Windows (chosen from the throttle, current, accelerometer and gyro-mean profile, before any covariance):
  W1 bench    tools/sess6.csv, the 5 s bins with median current > 0.8 A (as in step 7a): ballasted, still, about 1 A.
  W2 thrust   tools/syntonic_00000000_000000_00005.csv, 15.0 to 23.0 s: throttle 0.436, 9.8 A, accelerometer at 1 g and
              gyro mean near zero, so the frame is held; rotor vibration at about half throttle.
  W3 flight   tools/crash7.csv (session 7), 19.0 to 24.0 s: airborne (accelerometer off 1 g, non-zero gyro means),
              before the throttle cut at about 24.5 s. Five seconds, so its numbers carry wide intervals.
Gyro: four 1 kHz snapshots per record, chip frame mapped to body (x, -y, -z), rad/s.

Quantities, 3x3, per window (Hadamard: second differences of block means, d, C(l) = E[d d^T]):
  Sn      white-noise intensity of the rate measurement: mean over l = 2, 4, 8 ms of C(l) l / 6.
  Lam     random-walk intensity of the rate (what the rate path tracks): mean over l = 32, 64, 128 ms of
          (C(l) - 6 Sn / l) * 3 / (2 l). (For white noise s and a random walk q: E[d d^T] = 6 s / l + 2 q l / 3.)
          If Lam is not positive definite the window is reported and D2r is not computed for it.
  tau_i   the engine's own per-axis horizons, medians of tau_R, tau_P, tau_Y over the window (logged).

Decision rules:
  D1r (reported): correlation matrices of Sn and Lam, commutator ratio ||[Lam, Sn]|| / (||Lam|| ||Sn||).
  D2r (build decision, same rule as 7a D2): the rate path as a filter on a random walk Lam measured with white noise
      Sn: tr P(T_diag) / tr P(T_K) - 1, T_diag = diag(sqrt(Sn_ii / Lam_ii)), T_K the 3-axis steady-state Kalman window.
      Above 5 % in W2 or W3: a matrix-law rate bench is worth building. At or under 5 %: per-axis enough.
  D4  (TH44/VM-2 certificate): generator Q = diag(1 / tau_i), environmental covariance Sigma = R = Lam T_env
      (independent innovations, H21 with R = Sigma), T_env in {1, 4, 16} ms. Reported: chi_full =
      lambda_max(Sigma^-1/2 P* Sigma^-1/2), P* = sum_k A^k R A^k, A = exp(-Q T_env); chi_axes = the same with the
      off-diagonal terms of Lam set to zero (the per-axis picture, which is the slowest axis); coupling factor
      chi_full / chi_axes - 1; the collapse direction v_c (dominant generalized eigenvector of (P*, Sigma)), its squared
      component shares. Rule: v_c is "single-axis" if its largest share >= 0.90, else "mixed"; the coupling matters
      for the certificate if the factor exceeds 5 %. The absolute chi is NOT interpreted: T_env (the interval between
      independent environmental events) is not identified on this machine. With equal tau_i the factor is zero
      identically (isotropic generator), so any effect comes from unequal horizons and correlated innovations.
Uncertainty: moving-block bootstrap of the second-difference sequences (blocks of 4, each scale resampled on its
own, 400 draws), 90 % intervals for D2r and the D4 factor and shares.
Limits declared: W1 and W2 are on a held frame (Lam there is vibration and drift, not motion); W3 is one short
flight that ended in a crash; T_env is not identified; yaw attitude (magnetometer) is still out of reach.
"""
import hashlib, json, os, sys
import numpy as np
import pandas as pd
from scipy.linalg import sqrtm, solve_continuous_lyapunov, inv, expm, eigh

HERE = os.path.dirname(os.path.abspath(__file__))
UP = '/mnt/user-data/uploads/SyntonicQuad_NuttX/tools/'
LSB_DPS = 16.384; D2R = np.pi / 180
LS = (2, 4, 8); LL = (32, 64, 128); TENV = (1e-3, 4e-3, 16e-3)
NBOOT = 400; BLK = 4


def gmean(A, B):
    Ah = np.real(sqrtm(A)); Aih = inv(Ah)
    return np.real(Ah @ sqrtm(Aih @ B @ Aih) @ Ah)


def trP(T, Q, S):
    K = inv(T)
    return float(np.trace(solve_continuous_lyapunov(-K, -(Q + K @ S @ K.T))))


def d2r(Lam, Sn):
    Td = np.diag(np.sqrt(np.diag(Sn) / np.diag(Lam)))
    P = gmean(Sn, Lam); Tk = inv(P @ inv(Sn))
    return trP(Td, Lam, Sn) / trP(Tk, Lam, Sn) - 1, Td, Tk


def th44(tau, Lam, Tenv):
    q = 1.0 / np.asarray(tau)
    def chi(L):
        Sig = L * Tenv
        H = 1.0 / np.expm1((q[:, None] + q[None, :]) * Tenv)
        P = H * Sig                                   # Q is diagonal: U = I, P*_ij = Sig_ij / (e^{(qi+qj)T} - 1)
        w, V = eigh(P, Sig)                           # generalized: P v = chi Sigma v
        v = V[:, -1]; v = v / np.linalg.norm(v)
        return w[-1], v ** 2
    cf, sh = chi(Lam); ca, sha = chi(np.diag(np.diag(Lam)))
    return cf, ca, cf / ca - 1, sh, sha


def segs(df, t, mask_fn):
    """contiguous record runs where mask is true -> list of 1 kHz gyro arrays (n x 3, rad/s)"""
    m = mask_fn(t)
    G = np.stack([np.stack([df[f'gx_s{k}'].values, -df[f'gy_s{k}'].values, -df[f'gz_s{k}'].values], 1)
                  for k in range(4)], 1) / LSB_DPS * D2R                     # records x 4 x 3
    out = []; i = 0; n = len(m)
    while i < n:
        if not m[i]: i += 1; continue
        j = i
        while j + 1 < n and m[j + 1]: j += 1
        out.append(G[i:j + 1].reshape(-1, 3)); i = j + 1
    return out


def dseq(runs, l):
    ds = []
    for X in runs:
        if len(X) < 4 * l: continue
        mb = X[: len(X) // l * l].reshape(-1, l, 3).mean(1)
        ds.append(mb[2:] - 2 * mb[1:-1] + mb[:-2])
    return np.concatenate(ds) if ds else np.zeros((0, 3))


def estimate(D):
    Sn = np.mean([D[l].T @ D[l] / len(D[l]) * (l * 1e-3) / 6 for l in LS], 0)
    Lams = {l: (D[l].T @ D[l] / len(D[l]) - 6 * Sn / (l * 1e-3)) * 3 / (2 * l * 1e-3) for l in LL}
    return Sn, np.mean(list(Lams.values()), 0), Lams


def corr(M):
    s = np.sqrt(np.diag(M)); return M / np.outer(s, s)


def boot(D, rng):
    Db = {}
    for l, d in D.items():
        n = len(d); nb = int(np.ceil(n / BLK)); st = rng.integers(0, max(1, n - BLK + 1), nb)
        idx = (st[:, None] + np.arange(BLK)[None, :]).ravel()[:n]
        Db[l] = d[np.minimum(idx, n - 1)]
    return Db


def main():
    sha = hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest()
    out = [f'step7c_rate3_th44.py sha256 {sha}']
    P_ = lambda *a: out.append(' '.join(str(x) for x in a))
    W = {}
    df = pd.read_csv(UP + 'sess6.csv'); t = (df.t_us.values - df.t_us.values[0]) * 1e-6
    bins = np.floor(t / 5).astype(int); med = pd.Series(df.amps_A.values).groupby(bins).transform('median').values
    W['W1 bench'] = (df, t, lambda tt: med > 0.8)
    df5 = pd.read_csv(UP + 'syntonic_00000000_000000_00005.csv'); t5 = (df5.t_us.values - df5.t_us.values[0]) * 1e-6
    W['W2 thrust'] = (df5, t5, lambda tt: (tt >= 15.0) & (tt < 23.0))
    df7 = pd.read_csv(UP + 'crash7.csv'); t7 = (df7.t_us.values - df7.t_us.values[0]) * 1e-6
    W['W3 flight'] = (df7, t7, lambda tt: (tt >= 19.0) & (tt < 24.0))
    rng = np.random.default_rng(77); res = {}
    names = ['R', 'P', 'Y']
    for wn, (d, tt, fn) in W.items():
        runs = segs(d, tt, fn); m = fn(tt)
        gaps = int(np.sum(np.diff(tt[m]) > 0.0055))
        tau = [float(np.median(d[f'tau_{a}'].values[m])) for a in names]
        D = {l: dseq(runs, l) for l in LS + LL}
        Sn, Lam, Lams = estimate(D)
        P_(f'\n== {wn}: {sum(len(r) for r in runs)/1e3:.1f} s of 1 kHz gyro, {len(runs)} run(s), record gaps > 5.5 ms: {gaps};'
           f' n second differences at 128 ms: {len(D[128])}; engine tau medians R/P/Y = {np.round(np.array(tau)*1e3,1).tolist()} ms;'
           f' throttle median {np.median(d.thr.values[m]):.3f}, current {np.median(d.amps_A.values[m]):.2f} A')
        P_('  Sn diag (rad^2/s):', np.array2string(np.diag(Sn), precision=3), ' corr(Sn) =', np.array2string(corr(Sn), precision=3).replace('\n', ''))
        P_('  Lam diag (rad^2/s^3) per scale:', {l: np.round(np.diag(v), 6).tolist() for l, v in Lams.items()})
        pd_ok = bool(np.all(np.linalg.eigvalsh(Lam) > 0))
        P_('  Lam diag:', np.array2string(np.diag(Lam), precision=4), ' positive definite:', pd_ok)
        r = dict(seconds=sum(len(x) for x in runs) / 1e3, gaps=gaps, tau=tau, Sn=Sn.tolist(), Lam=Lam.tolist(), Lam_pd=pd_ok)
        if not pd_ok:
            res[wn] = r; P_('  Lam not positive definite: D2r and D4 not computed'); continue
        P_('  corr(Lam) =', np.array2string(corr(Lam), precision=3).replace('\n', ''))
        comm = np.linalg.norm(Lam @ Sn - Sn @ Lam) / (np.linalg.norm(Lam) * np.linalg.norm(Sn))
        g, Td, Tk = d2r(Lam, Sn)
        P_(f'  commutator ratio {comm:.4f}; per-axis windows sqrt(Sn/Lam) = {np.round(np.diag(Td)*1e3, 2).tolist()} ms')
        bg = []; bt = {T: [] for T in TENV}
        for _ in range(NBOOT):
            Snb, Lamb, _ = estimate(boot(D, rng))
            if np.all(np.linalg.eigvalsh(Lamb) > 0) and np.all(np.linalg.eigvalsh(Snb) > 0):
                bg.append(d2r(Lamb, Snb)[0])
                for T in TENV:
                    cf, ca, fac, sh, _ = th44(tau, Lamb, T); bt[T].append([fac, sh.max(), *sh])
        bg = np.array(bg)
        P_(f'  D2r: tr P(per-axis) / tr P(3-axis Kalman) - 1 = {100*g:.3f} %   (90 % boot {100*np.quantile(bg,0.05):.3f} .. {100*np.quantile(bg,0.95):.3f} %, {len(bg)} valid draws)')
        r.update(corr_Sn=corr(Sn).tolist(), corr_Lam=corr(Lam).tolist(), comm=comm, D2r=g, D2r_ci=[np.quantile(bg, 0.05), np.quantile(bg, 0.95)], D4={})
        for T in TENV:
            cf, ca, fac, sh, sha_ = th44(tau, Lam, T); b = np.array(bt[T])
            kind = 'single-axis' if sh.max() >= 0.90 else 'mixed'
            P_(f'  D4 T_env {T*1e3:.0f} ms: chi_full {cf:.4g}, chi_axes {ca:.4g} -> coupling factor {100*fac:+.3f} %'
               f' (90 % {100*np.quantile(b[:,0],0.05):+.3f} .. {100*np.quantile(b[:,0],0.95):+.3f} %);'
               f' v_c shares R/P/Y {np.round(sh,3).tolist()} -> {kind} (largest share 90 % {np.quantile(b[:,1],0.05):.3f} .. {np.quantile(b[:,1],0.95):.3f});'
               f' per-axis v_c {np.round(sha_,3).tolist()}')
            r['D4'][str(T)] = dict(chi_full=cf, chi_axes=ca, factor=fac, shares=sh.tolist(), kind=kind,
                                   factor_ci=[np.quantile(b[:, 0], 0.05), np.quantile(b[:, 0], 0.95)],
                                   share_max_ci=[np.quantile(b[:, 1], 0.05), np.quantile(b[:, 1], 0.95)])
        res[wn] = r
    P_('\n== Decisions')
    for wn in ('W2 thrust', 'W3 flight'):
        r = res[wn]
        if 'D2r' in r:
            P_(f'D2r {wn}: {100*r["D2r"]:.3f} % -> ' + ('BUILD a matrix-law rate bench (above 5 %)' if r['D2r'] > 0.05 else 'per-axis enough (at or under 5 %)'))
    for wn, r in res.items():
        if 'D4' in r:
            P_(f'D4 {wn}: ' + '; '.join(f'T_env {float(T)*1e3:.0f} ms factor {100*v["factor"]:+.2f} %, v_c {v["kind"]}' for T, v in r['D4'].items()))
    txt = '\n'.join(out); print(txt)
    os.makedirs(os.path.join(HERE, 'results'), exist_ok=True)
    open(os.path.join(HERE, 'results', 'step7c_results.txt'), 'w').write(txt + '\n')
    json.dump({'sha256': sha, **res}, open(os.path.join(HERE, 'results', 'step7c_rows.json'), 'w'), default=float)


if __name__ == '__main__':
    main()
