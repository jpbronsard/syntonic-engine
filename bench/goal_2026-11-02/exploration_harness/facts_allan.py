"""Noise facts from the logged records only (no truth, no score): Allan variances, 2 Nov 2026, exploratory.

For each record (idle 35 s, low 40 s, half 8 s):
  gyro rate Allan variance (1 kHz, overlapping) -> ARW region, bias, long-scale slope;
  tilt Allan variance (250 Hz, overlapping);
  residual z = tilt - open-loop integral of the arming-calibrated gyro (gyro minus the idle record's mean, which is what
  the arming calibration subtracts in expectation), Allan variance of z (an angle signal: white tilt -> r/T, gyro angle
  random walk q -> q T / 3, residual bias b -> b^2 T^2 / 2);
  the bench's live reading s2 = min over 32/64/128 ms of D_H(l) l / 6 (Hadamard), and the law's window sqrt(s2/lambda).
Two versions: 'raw' (the record as logged, scales up to N/2) and 'wrap' (the record looped circularly, as build_noise
replays it; scales up to 32 s). Output: results/facts_allan.json, printed table.
"""
import numpy as np
import common as C
F, G, D2R, DT, DTA = C.F, C.G, C.D2R, C.DT, C.DTA


def avar(x, dt, ms, wrap=False):
    """overlapping Allan variance of x (sampled at dt) at block lengths ms (samples). wrap: circular record."""
    out = []
    N = len(x)
    for m in ms:
        if wrap:
            xx = np.concatenate([x] * (int(np.ceil(2 * m / N)) + 2))
            X = np.concatenate([[0.0], np.cumsum(xx)])
            i = np.arange(N); mb0 = (X[i + m] - X[i]) / m; mb1 = (X[i + 2 * m] - X[i + m]) / m
        else:
            if 2 * m > N: out.append(np.nan); continue
            X = np.concatenate([[0.0], np.cumsum(x)]); i = np.arange(N - 2 * m + 1)
            mb0 = (X[i + m] - X[i]) / m; mb1 = (X[i + 2 * m] - X[i + m]) / m
        out.append(0.5 * np.mean((mb1 - mb0) ** 2))
    return np.array(out)


def hvar_s2(z, ls=(8, 16, 32)):
    """the bench's live reading offline: D_H(l) l / 6 on non-overlapping block means, per scale, and the min"""
    s = []
    for l in ls:
        mb = z[: len(z) // l * l].reshape(-1, l).mean(1); dd = mb[2:] - 2 * mb[1:-1] + mb[:-2]
        s.append(np.mean(dd ** 2) * l * DTA / 6.0)
    return np.array(s)


def main():
    R = F.records()
    b_idle = float(R['idle'][1].mean())
    Ta = DTA * 2.0 ** np.arange(0, 13)                     # 4 ms .. 16.4 s
    Ta_w = DTA * 2.0 ** np.arange(0, 14)                   # .. 32.8 s (wrapped)
    Tg = DT * 2.0 ** np.arange(0, 15)                      # 1 ms .. 16.4 s
    out = dict(b_idle_dps=b_idle / D2R, scales_tilt=Ta_w.tolist(), scales_gyro=Tg.tolist(), rec={})
    for name in ('idle', 'low', 'half'):
        tilt, g = R[name]; N = len(tilt)
        gc = g - b_idle                                    # arming-calibrated gyro (expected calibration)
        th_g = DT * np.cumsum(gc)[NSUB_IDX(N)]             # open-loop integral at each tilt sample
        z = tilt - th_g
        d = dict(T_s=N * DTA, bias_dps=float(g.mean() / D2R), resid_bias_dps=float(gc.mean() / D2R),
                 bias_halves_dps=[float(g[: len(g) // 2].mean() / D2R), float(g[len(g) // 2:].mean() / D2R)],
                 tilt_sd_deg=float(tilt.std() / D2R), gyro_sd_dps=float(g.std() / D2R))
        d['avar_gyro_raw'] = avar(g, DT, (Tg / DT).round().astype(int)).tolist()
        d['avar_tilt_raw'] = avar(tilt, DTA, (Ta_w / DTA).round().astype(int)).tolist()
        d['avar_tilt_wrap'] = avar(tilt, DTA, (Ta_w / DTA).round().astype(int), wrap=True).tolist()
        d['avar_z_raw'] = avar(z, DTA, (Ta_w / DTA).round().astype(int)).tolist()
        # wrapped residual: loop the calibrated gyro and the tilt together, the integral keeps its ramp
        zw = loop_resid(tilt, gc, 4)                        # 4 loops of the record, then circular Allan is not valid:
        d['avar_z_loop4'] = avar(zw, DTA, (Ta_w / DTA).round().astype(int)).tolist()   # plain overlapping on 4 loops
        d['s2_bench_scales'] = hvar_s2(z).tolist(); d['s2_bench'] = float(min(d['s2_bench_scales']))
        d['s2_tilt_only'] = float(min(hvar_s2(tilt)))
        d['r_white_persample'] = float(tilt.var() * DTA)   # per-sample variance x dt (white-equivalent if white)
        lam = C.FX['lam'][name]; d['lam_fact'] = lam
        d['tau_law'] = float(np.sqrt(d['s2_bench'] / lam))
        out['rec'][name] = d
    C.dump('facts_allan.json', out)
    # ---- print
    print(f"idle bias (arming calibration in expectation) {b_idle / D2R:.4f} deg/s")
    for name, d in out['rec'].items():
        print(f"\n== {name}: {d['T_s']:.1f} s, bias {d['bias_dps']:.4f} deg/s (halves {d['bias_halves_dps'][0]:.4f} / {d['bias_halves_dps'][1]:.4f}), "
              f"residual after idle calibration {d['resid_bias_dps']:+.4f} deg/s; tilt sd {d['tilt_sd_deg']:.2f} deg, gyro sd {d['gyro_sd_dps']:.2f} deg/s")
        print(f"   s2_bench (z, Hadamard 32/64/128 ms) {d['s2_bench']:.3e} rad^2 s  (per scale {np.array(d['s2_bench_scales'])}), "
              f"tilt alone {d['s2_tilt_only']:.3e}, white-equivalent per-sample {d['r_white_persample']:.3e}")
        print(f"   lambda fact {d['lam_fact']:.3e} rad^2/s -> tau_law = sqrt(s2/lambda) = {d['tau_law']:.2f} s")
        print('   gyro rate Allan: T, sigma2*T (rad^2/s), slope of adev')
        ag = np.array(d['avar_gyro_raw'])
        for i, T in enumerate(Tg):
            sl = 0.5 * np.log(ag[i] / ag[i - 1]) / np.log(2) if i > 0 and np.isfinite(ag[i]) else np.nan
            if np.isfinite(ag[i]): print(f"     {T:8.3f}  adev {np.sqrt(ag[i]) / D2R:9.5f} deg/s  s2T {ag[i] * T:.3e}  slope {sl:+.2f}")
        print('   tilt / z Allan (angle signals): T, s2*T raw tilt, wrap tilt, raw z, loop4 z (rad^2 s), z adev deg')
        at, aw, az, a4 = (np.array(d[k]) for k in ('avar_tilt_raw', 'avar_tilt_wrap', 'avar_z_raw', 'avar_z_loop4'))
        for i, T in enumerate(Ta_w):
            print(f"     {T:8.3f}  {at[i] * T:.3e}  {aw[i] * T:.3e}  {az[i] * T:.3e}  {a4[i] * T:.3e}   adev z raw {np.sqrt(az[i]) / D2R:8.4f}  loop4 {np.sqrt(a4[i]) / D2R:8.4f}")


def NSUB_IDX(N):
    return np.arange(N) * C.NSUB             # tick index of tilt sample j (the estimator integrates ticks 0..4j)


def loop_resid(tilt, gc, nloop):
    t = np.concatenate([tilt] * nloop); g = np.concatenate([gc] * nloop)
    return t - DT * np.cumsum(g)[np.arange(len(t)) * C.NSUB]


if __name__ == '__main__':
    main()
