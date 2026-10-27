"""STEP 7a, decision rules written before any coupling was computed. Does the attitude estimator of the X500 need
the matrix law (V5.1, Theorem 35.1), or are per-axis windows enough?

Data. Session 6 of the SD binary logger (tools/sess6.csv, decoded schema v3): motor test of the X500 with 2 kg of
ballast per leg, the frame provably immobile, 135 s, current 0.5 to 1.2 A (idle and low throttle only). Raw IMU:
accelerometer at 250 Hz, gyroscope at 1 kHz (four snapshots per record), chip frame, mapped to the body frame of
the firmware (x, -y, -z; the further board remap is a sign/permutation applied to both sensors and leaves every
quantity below invariant). Because the frame is still, the accelerometer tilt and the gyroscope rate are noise
(plus vibration and slow drift).

Quantities (per regime: idle, current < 0.65 A; throttle, current > 0.8 A; 5 s bins assigned by their median
current, contiguous runs of bins analysed separately and pooled):
  Sigma  2x2 low-frequency PSD matrix of the accelerometer tilt noise (roll = atan2(-ay, -az), pitch =
         atan2(ax, sqrt(ay^2 + az^2)), firmware definitions), from second differences of block means (Hadamard,
         removes slow drift): C(l) = E[d d^T], Sigma = C(l) l / 6, l = 64, 128, 256 ms; the three scales averaged,
         their spread reported.
  Q      2x2 low-frequency PSD matrix of the gyroscope rate noise (roll and pitch rates, rad/s), same
         construction at l = 128, 256, 512, 1024 ms. (Q here is the drift intensity lambda of the fusion.)
  Windows (kappa = 1):
         T_diag = diag(sqrt(Sigma_ii / Q_ii))                     the per-axis law (what step 6b runs, per axis);
         T_sym  = Q^-1 # Sigma                                     the matrix law of V5.1 Theorem 35.1;
         T_K    = (K_opt)^-1, K_opt = P Sigma^-1, P = Sigma # Q   the steady-state two-axis Kalman gain
                                                                   (random-walk model, both axes measured).
  Predicted steady-state error covariance of the filter xhat' = T^-1 (y - xhat) on a random walk of intensity Q
  with white noise Sigma: -K P - P K^T + Q + K Sigma K^T = 0, K = T^-1; score tr(P).

Decision rules:
  D1 (reported, no threshold): noise correlations rho_Sigma, rho_Q, and the commutator ratio
     ||[Q, Sigma]||_F / (||Q||_F ||Sigma||_F), per regime.
  D2 (build decision): if tr P(T_diag) / tr P(T_K) - 1 > 5 % in the throttle regime, a matrix-law bench (step 7b)
     is worth building; otherwise the per-axis law is enough on this machine in the regimes measured.
  D3 (a point of theory checked numerically, prediction written here): T_sym and T_K both satisfy T Q T^T = Sigma;
     they coincide, and score the same, when Q and Sigma commute; when they do not, T_K (the MSE optimum by
     construction) scores strictly better. Checked on the measured matrices and on 1000 random SPD pairs.
Limits declared: low throttle only, ballasted and still (flight vibration is larger and differently shaped); yaw
is left out (its absolute sensor is the magnetometer, not in this log); 2x2 only.
"""
import hashlib, json, os, sys
import numpy as np
import pandas as pd
from scipy.linalg import sqrtm, solve_continuous_lyapunov, inv

HERE = os.path.dirname(os.path.abspath(__file__))
CSV = sys.argv[1] if len(sys.argv) > 1 else '/mnt/user-data/uploads/SyntonicQuad_NuttX/tools/sess6.csv'
LSB_DPS = 16.384; D2R = np.pi / 180


def gmean(A, B):
    """A # B, the Pusz-Woronowicz geometric mean of SPD matrices."""
    Ah = np.real(sqrtm(A)); Aih = inv(Ah)
    return np.real(Ah @ sqrtm(Aih @ B @ Aih) @ Ah)


def trP(T, Q, S):
    K = inv(T)
    P = solve_continuous_lyapunov(-K, -(Q + K @ S @ K.T))
    return float(np.trace(P)), P


def windows(Q, S):
    Td = np.diag(np.sqrt(np.diag(S) / np.diag(Q)))
    Ts = gmean(inv(Q), S)
    P = gmean(S, Q); Kopt = P @ inv(S); Tk = inv(Kopt)
    return Td, Ts, Tk


def had_cov(X, n):
    """second differences of block means of the columns of X (blocks of n samples): E[d d^T]."""
    m = X[: len(X) // n * n].reshape(-1, n, X.shape[1]).mean(1)
    d = m[2:] - 2 * m[1:-1] + m[:-2]
    return d.T @ d, len(d)


def main():
    sha = hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest()
    out = [f'step7_axis_coupling.py sha256 {sha}']
    P_ = lambda *a: out.append(' '.join(str(x) for x in a))
    df = pd.read_csv(CSV, usecols=['t_us', 'armed', 'amps_A', 'ax_raw', 'ay_raw', 'az_raw'] +
                     [f'g{c}_s{k}' for k in range(4) for c in 'xyz'])
    t = (df.t_us.values - df.t_us.values[0]) * 1e-6
    ax, ay, az = df.ax_raw.values, -df.ay_raw.values, -df.az_raw.values          # body frame (x, -y, -z)
    roll = np.arctan2(-ay, -az); pitch = np.arctan2(ax, np.sqrt(ay ** 2 + az ** 2))
    G = np.stack([np.stack([df[f'gx_s{k}'].values, -df[f'gy_s{k}'].values, -df[f'gz_s{k}'].values], 1)
                  for k in range(4)], 1).reshape(-1, 3) / LSB_DPS * D2R                # 1 kHz, rad/s
    tg = np.repeat(t, 4) + np.tile(np.arange(4) * 1e-3 - 3e-3, len(t))
    dt_rec = np.median(np.diff(t)); P_(f'records {len(t)}, duration {t[-1]:.1f} s, record period {1e3*dt_rec:.2f} ms, armed {df.armed.min()}..{df.armed.max()}')
    # regimes from 5 s bins
    nb = int(t[-1] // 5) + 1; lab = []
    for b in range(nb):
        a = df.amps_A.values[(t >= 5 * b) & (t < 5 * b + 5)]
        lab.append('idle' if np.median(a) < 0.65 else ('throttle' if np.median(a) > 0.8 else 'mixed'))
    runs = []; b = 0
    while b < nb:
        e = b
        while e + 1 < nb and lab[e + 1] == lab[b]: e += 1
        if lab[b] != 'mixed': runs.append((lab[b], 5 * b, 5 * (e + 1)))
        b = e + 1
    P_('runs:', [(r, a, z) for r, a, z in runs])
    res = {}
    for reg in ('idle', 'throttle'):
        Cs = {l: [np.zeros((2, 2)), 0] for l in (64, 128, 256)}; Cq = {l: [np.zeros((2, 2)), 0] for l in (128, 256, 512, 1024)}
        dur = 0.0
        for r, a, z in runs:
            if r != reg: continue
            m = (t >= a) & (t < z); mg = (tg >= a) & (tg < z); dur += z - a
            X = np.stack([roll[m], pitch[m]], 1); Y = G[mg][:, :2]
            for l in Cs:
                n = int(round(l * 1e-3 / dt_rec))
                if len(X) >= 4 * n: c, k = had_cov(X, n); Cs[l][0] += c; Cs[l][1] += k
            for l in Cq:
                n = int(round(l * 1e-3 / 1e-3))
                if len(Y) >= 4 * n: c, k = had_cov(Y, n); Cq[l][0] += c; Cq[l][1] += k
        Ss = {l: Cs[l][0] / Cs[l][1] * (l * 1e-3) / 6 for l in Cs if Cs[l][1] > 0}
        Qs = {l: Cq[l][0] / Cq[l][1] * (l * 1e-3) / 6 for l in Cq if Cq[l][1] > 0}
        S = np.mean(list(Ss.values()), 0); Q = np.mean(list(Qs.values()), 0)
        rS = S[0, 1] / np.sqrt(S[0, 0] * S[1, 1]); rQ = Q[0, 1] / np.sqrt(Q[0, 0] * Q[1, 1])
        comm = np.linalg.norm(Q @ S - S @ Q) / (np.linalg.norm(Q) * np.linalg.norm(S))
        Td, Ts, Tk = windows(Q, S)
        jd, js, jk = trP(Td, Q, S)[0], trP(Ts, Q, S)[0], trP(Tk, Q, S)[0]
        P_(f'\n== {reg}: {dur:.0f} s')
        P_('  Sigma per scale (diag, rad^2 s):', {l: np.round(np.diag(v), 12).tolist() for l, v in Ss.items()})
        P_('  Q per scale (diag, rad^2/s):', {l: np.round(np.diag(v), 12).tolist() for l, v in Qs.items()})
        P_(f'  Sigma = {np.array2string(S, precision=3)}  rho_Sigma = {rS:+.3f}')
        P_(f'  Q     = {np.array2string(Q, precision=3)}  rho_Q = {rQ:+.3f}')
        P_(f'  commutator ratio ||[Q,S]|| / (||Q|| ||S||) = {comm:.4f}')
        P_(f'  per-axis windows sqrt(S_ii/Q_ii) = {np.round(np.diag(Td), 3).tolist()} s')
        P_(f'  T_sym = {np.array2string(Ts, precision=3)}   T_K = {np.array2string(Tk, precision=3)}')
        P_(f'  tr P: diag {jd:.4e}  sym {js:.4e}  Kalman {jk:.4e}  ->  diag/K - 1 = {100*(jd/jk-1):.3f} %,  sym/K - 1 = {100*(js/jk-1):.3f} %')
        res[reg] = dict(Sigma=S.tolist(), Q=Q.tolist(), rho_Sigma=rS, rho_Q=rQ, comm=comm, T_diag=Td.tolist(), T_sym=Ts.tolist(),
                        T_K=Tk.tolist(), trP_diag=jd, trP_sym=js, trP_K=jk, Sigma_scales={str(k): v.tolist() for k, v in Ss.items()},
                        Q_scales={str(k): v.tolist() for k, v in Qs.items()}, duration_s=dur)
    # D3 on random SPD pairs
    rng = np.random.default_rng(7); gaps = []; comms = []
    for _ in range(1000):
        A = rng.standard_normal((2, 2)); B = rng.standard_normal((2, 2))
        Qr = A @ A.T + 0.1 * np.eye(2); Sr = B @ B.T + 0.1 * np.eye(2)
        Td, Ts, Tk = windows(Qr, Sr)
        gaps.append(trP(Ts, Qr, Sr)[0] / trP(Tk, Qr, Sr)[0] - 1)
        comms.append(np.linalg.norm(Qr @ Sr - Sr @ Qr) / (np.linalg.norm(Qr) * np.linalg.norm(Sr)))
    gaps = np.array(gaps); comms = np.array(comms)
    Qc = np.diag([1.0, 3.0]); Sc = np.diag([2.0, 0.5])
    Td, Ts, Tk = windows(Qc, Sc); g0 = trP(Ts, Qc, Sc)[0] / trP(Tk, Qc, Sc)[0] - 1
    P_('\n== D3, random SPD pairs (1000): sym/K - 1  min %.2e  median %.3e  max %.3e ; fraction > 1e-9: %.3f ; commuting case: %.1e'
       % (gaps.min(), np.median(gaps), gaps.max(), np.mean(gaps > 1e-9), g0))
    P_('   corr(log gap, log commutator) = %.2f' % np.corrcoef(np.log(gaps[gaps > 1e-12]), np.log(comms[gaps > 1e-12]))[0, 1])
    thr = res['throttle']; g = thr['trP_diag'] / thr['trP_K'] - 1
    P_('\n== Decisions')
    P_(f'D2: throttle regime, per-axis excess over the two-axis Kalman = {100*g:.3f} % -> '
       + ('BUILD step 7b (above 5 %)' if g > 0.05 else 'per-axis law enough here (at or under 5 %)'))
    P_(f'D3: measured, sym/K - 1 = idle {100*(res["idle"]["trP_sym"]/res["idle"]["trP_K"]-1):.4f} %, throttle {100*(thr["trP_sym"]/thr["trP_K"]-1):.4f} %;'
       f' random pairs: T_K never worse (min gap {gaps.min():.1e}), gap > 0 in {100*np.mean(gaps > 1e-9):.1f} % of pairs')
    txt = '\n'.join(out); print(txt)
    os.makedirs(os.path.join(HERE, 'results'), exist_ok=True)
    open(os.path.join(HERE, 'results', 'step7_results.txt'), 'w').write(txt + '\n')
    json.dump({'sha256': sha, **res, 'D3_random_gap': {'min': gaps.min(), 'median': float(np.median(gaps)), 'max': gaps.max()}},
              open(os.path.join(HERE, 'results', 'step7_rows.json'), 'w'), default=float)


if __name__ == '__main__':
    main()
