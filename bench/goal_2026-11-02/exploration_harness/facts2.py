"""Facts for the new rules, from the logged records only (no truth, no score), 2 Nov 2026, exploratory.

1. Two-state facts from the idle record (at rest, once):
   q_theta = the bench's lam_hat (gyro angle random walk, Allan 0.128-1.024 s, idle) [rad^2/s];
   q_b = largest rate random walk consistent with the idle Allan curve: the +1/2-slope asymptote sigma_A^2(T) = q_b T / 3
         through the longest scale with >= 4 non-overlapping blocks (T <= T_rec/4): q_b = 3 sigma_A^2(T) / T [rad^2/s^3];
         sensitivity: the same at T <= T_rec/2.
2. Check of the closed-form steady-state gains against scipy's continuous ARE and the discrete ARE at the accelerometer
   rate (pole-matched discretisation K1 = 1 - z1 z2, K2 = (1 - z1)(1 - z2) / T, z = exp(s T)).
3. Operating-scale curves per throttle level (record 'raw' = as logged, and 'loop' = circularly looped as the bench
   replays it), on a log grid of scales T (2^(k/4) x 4 ms):
   r(T)   = sigma_tilt^2(T) T          (overlapping Allan variance of the tilt record; white -> r);
   lam(T) = E[mean_T(omega_c)^2] T     (omega_c = gyro minus the idle mean, the arming calibration in expectation; white
                                        rate noise q -> q; a residual bias b adds b^2 T);
   sz(T)  = sigma_z^2(T)               (Allan variance of the residual z = tilt - open-loop integral of omega_c);
   s2_rec = the bench's live reading computed offline on the record (min over 32/64/128 ms of D_H(l) l / 6).
   'raw' curves are measured up to T_lim = T_rec/4 and extrapolated beyond with the canonical forms
   (r flat at its top-octave mean; lam = top-octave Allan part + b_res^2 T; sz = r/T + q T/3 + b_res^2 T^2/2).
Output: results/facts2.json"""
import numpy as np
from scipy.linalg import solve_continuous_are, solve_discrete_are
import common as C
F, D2R, DT, DTA = C.F, C.D2R, C.DT, C.DTA
TG = DTA * 2.0 ** (np.arange(0, 4 * 15) / 4.0)          # 4 ms .. ~130 s, quarter-octave


def blockmeans_overlap(x, m, wrap):
    N = len(x)
    if wrap:
        xx = np.concatenate([x] * (int(np.ceil(2 * m / N)) + 2)); X = np.concatenate([[0.0], np.cumsum(xx)]); i = np.arange(N)
    else:
        X = np.concatenate([[0.0], np.cumsum(x)]); i = np.arange(max(N - 2 * m + 1, 0))
    return (X[i + m] - X[i]) / m, (X[i + 2 * m] - X[i + m]) / m


def allan(x, m, wrap=False):
    if not wrap and 2 * m > len(x): return np.nan
    a, b = blockmeans_overlap(x, m, wrap); return 0.5 * np.mean((b - a) ** 2)


def msq_mean(x, m, wrap=False):
    if not wrap and 2 * m > len(x): return np.nan
    a, _ = blockmeans_overlap(x, m, wrap); return np.mean(a ** 2)


def s2_bench(z):
    s = []
    for l in (8, 16, 32):
        mb = z[: len(z) // l * l].reshape(-1, l).mean(1); dd = mb[2:] - 2 * mb[1:-1] + mb[:-2]; s.append(np.mean(dd ** 2) * l * DTA / 6.0)
    return min(s)


def gains_closed(r, qt, qb):
    k2 = np.sqrt(qb / r); k1 = np.sqrt(qt / r + 2 * k2); return k1, k2


def disc_gains(k1, k2, T=DTA):
    sq = np.sqrt(np.asarray(k1 * k1 - 4 * k2, complex)); s1 = (-k1 + sq) / 2; s2 = (-k1 - sq) / 2
    z1, z2 = np.exp(s1 * T), np.exp(s2 * T)
    return np.real(1 - z1 * z2), np.real((1 - z1) * (1 - z2)) / T


def check_are(r, qt, qb):
    """returns closed-form (k1, k2), CARE gains, and the DARE closed-loop poles vs pole-matched poles"""
    A = np.array([[0.0, -1.0], [0.0, 0.0]]); Cm = np.array([[1.0, 0.0]]); Q = np.diag([qt, qb])
    P = solve_continuous_are(A.T, Cm.T, Q, np.array([[r]])); K = P @ Cm.T / r
    k1, k2 = gains_closed(r, qt, qb)
    T = DTA; Fm = np.array([[1.0, -T], [0.0, 1.0]])
    Qd = np.array([[qt * T + qb * T ** 3 / 3, -qb * T ** 2 / 2], [-qb * T ** 2 / 2, qb * T]]); Rd = r / T
    Pd = solve_discrete_are(Fm.T, Cm.T, Qd, np.array([[Rd]]))       # a priori covariance
    Kd = Pd @ Cm.T / (Cm @ Pd @ Cm.T + Rd)                               # update gain
    M = (np.eye(2) - Kd @ Cm) @ Fm; zd = np.sort_complex(np.linalg.eigvals(M))
    K1, K2 = disc_gains(k1, k2)
    Mp = (np.eye(2) - np.array([[K1], [-K2]]) @ Cm) @ Fm; zp = np.sort_complex(np.linalg.eigvals(Mp))
    return dict(k1=k1, k2=k2, care_K=[K[0, 0], K[1, 0]], dare_K=[Kd[0, 0], Kd[1, 0]], pm_K=[K1, -K2],
                dare_poles=zd.tolist(), pm_poles=zp.tolist())


def curves(name, R, b_idle, wrap):
    tilt, g = R[name]; wc = g - b_idle; N = len(tilt)
    z = tilt - DT * np.cumsum(wc)[np.arange(N) * C.NSUB]
    if wrap:   # looped residual: loop tilt and gyro together; the integral keeps its ramp -> use 8 loops, plain overlapping
        L = max(8, int(np.ceil(2.2 * TG[-1] / (N * DTA))))
        tt = np.concatenate([tilt] * L); ww = np.concatenate([wc] * L)
        zz = tt - DT * np.cumsum(ww)[np.arange(len(tt)) * C.NSUB]
    rT, lT, sz, aT = [], [], [], []
    for T in TG:
        ma = int(round(T / DTA)); mg = int(round(T / DT))
        rT.append(allan(tilt, ma, wrap) * T)
        lT.append(msq_mean(wc, mg, wrap) * T)
        aT.append(allan(wc, mg, wrap) * T)
        sz.append(allan(zz, ma) if wrap else allan(z, ma))
    rT, lT, sz, aT = map(np.array, (rT, lT, sz, aT))
    T_rec = N * DTA; b_res = float(wc.mean())
    out = dict(T_rec=T_rec, b_res=b_res, s2_rec=float(s2_bench(z)), r=rT, lam=lT, lamA=aT, sz=sz)
    if not wrap:      # canonical extrapolation beyond T_lim = T_rec / 4
        T_lim = T_rec / 4; ok = TG <= T_lim * 1.0001; top = ok & (TG >= T_lim / 2)
        r_top = float(np.mean(rT[top])); q_top = float(np.mean(aT[top]))
        ext = ~ok
        rT2, lT2, sz2 = rT.copy(), lT.copy(), sz.copy()
        rT2[ext] = r_top; lT2[ext] = q_top + b_res ** 2 * TG[ext]
        sz2[ext] = r_top / TG[ext] + q_top * TG[ext] / 3 + b_res ** 2 * TG[ext] ** 2 / 2
        out.update(T_lim=T_lim, r_top=r_top, q_top=q_top, r=rT2, lam=lT2, sz=sz2, measured=ok)
    return out


def main():
    R = F.records(); b_idle = float(R['idle'][1].mean())
    # ---- 1. two-state facts from idle
    gi = R['idle'][1]; Nrec = len(gi) * DT
    def qb_at(Tmax):
        Ts = [T for T in DT * 2.0 ** np.arange(0, 20) if T <= Tmax]; T = Ts[-1]
        return 3 * allan(gi, int(round(T / DT))) / T, T
    qb4, T4 = qb_at(Nrec / 4); qb2, T2 = qb_at(Nrec / 2)
    qt = C.FX['lam_hat']
    print(f'idle record {Nrec:.1f} s: q_theta = lam_hat = {qt:.3e} rad^2/s ({np.sqrt(qt) / D2R * 60:.3f} deg/sqrt(h)); '
          f'q_b (T_rec/4, T = {T4:.3f} s) = {qb4:.3e} rad^2/s^3; sensitivity q_b (T_rec/2, T = {T2:.3f} s) = {qb2:.3e}; bench Q_B0 = {F.Q_B:.3e}')
    # ---- 2. ARE check at representative r values
    chk = {}
    for nm, r in (('low', 2.2e-5), ('half', 3.5e-3), ('idle', 1.9e-5)):
        for qb in (qb4, qb2):
            c = check_are(r, qt, qb); chk[f'{nm}_{qb:.1e}'] = c
            print(f'  r={r:.1e} q_b={qb:.1e}: closed k1={c["k1"]:.5f} k2={c["k2"]:.3e} | CARE K={c["care_K"][0]:.5f},{c["care_K"][1]:.3e} | '
                  f'1/k1={1 / c["k1"]:.1f} s, zeta={c["k1"] / (2 * np.sqrt(c["k2"])):.3f} | DARE K={c["dare_K"][0]:.3e},{c["dare_K"][1]:.3e} '
                  f'pole-matched K={c["pm_K"][0]:.3e},{c["pm_K"][1]:.3e} | DARE poles {np.round(c["dare_poles"], 7)} PM poles {np.round(c["pm_poles"], 7)}')
    # ---- 3. operating-scale curves
    out = dict(q_theta=qt, q_b=qb4, q_b_T=T4, q_b_half=qb2, q_b_half_T=T2, b_idle=b_idle, TG=TG, are_check=chk, curves={})
    for wrap in (False, True):
        key = 'loop' if wrap else 'raw'; out['curves'][key] = {}
        for name in ('idle', 'low', 'half'):
            cv = curves(name, R, b_idle, wrap); out['curves'][key][name] = cv
    C.dump('facts2.json', out)
    for key in ('raw', 'loop'):
        print(f'\n== {key} curves: T, then per throttle r(T), lam(T), sigma_z(T) [deg]  (s2_rec, b_res deg/s)')
        for name in ('idle', 'low', 'half'):
            cv = out['curves'][key][name]
            print(f'  {name}: s2_rec {cv["s2_rec"]:.3e}, b_res {cv["b_res"] / D2R:+.4f} deg/s' + (f', T_lim {cv["T_lim"]:.1f} s, r_top {cv["r_top"]:.2e}, q_top {cv["q_top"]:.2e}' if key == 'raw' else ''))
        for i in range(0, len(TG), 4):
            row = f'   T {TG[i]:7.3f}'
            for name in ('idle', 'low', 'half'):
                cv = out['curves'][key][name]
                row += f' | {cv["r"][i]:.2e} {cv["lam"][i]:.2e} {np.sqrt(cv["sz"][i]) / D2R:7.4f}'
            print(row)


if __name__ == '__main__':
    main()
