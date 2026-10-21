"""v10: U' against PID over ten noise draws per condition. Written 21 Oct 2026, before any run of this file.

Question: is U' on par with the tuned PID at 1 kHz (servo frame 3 ms), beyond the single noise draw of v5c?

Design (fixed before the run):
- Engines: PID + interpolation with the v5 tunings (tuning_v5_pair_A/B.json, chosen on the rest of the movement with
  clean noise), and U' derived with the 2.22 ms floor of the v4 video engine (TH41-A limit, above the 2.1 ms loop
  delay; nothing tuned). Same plant, plan, excerpts, conditions and scores as v5c.
- Ten fresh noise draws per condition (A, B; steady noise, noise x10 from 10 to 20 s): one generator per condition
  (seeds 1000, 2000, 3000, 4000), one noise row per draw, the same rows for both engines (paired design).
- Statistic: per draw the ratio rho = tip rms (U') / tip rms (PID); per condition the geometric mean of rho and its
  95 % interval from the t distribution on log rho (n = 10). Jitter ratio reported the same way.
- Equivalence margin, declared: "on par" holds for an excerpt if the 95 % interval of rho with steady noise lies
  inside [0.90, 1.10].
Predictions:
E1  A steady: interval contains 1 and lies inside [0.90, 1.10].
E2  B steady: PID slightly ahead (rho about 1.05-1.08), interval inside [0.90, 1.10].
E3  Noise step: interval entirely below 1 on both excerpts (U' better).
Check before the draws: one draw with the v5c generator (seed 5, A steady) reproduces v5c for PID (0.73 mm) and U'
with the 2 ms floor (0.73 mm).
A failed line is reported as such."""
import os, sys, json, time, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
os.environ['BENCH_HZ'] = '1000'
import flow_v7_ctrl as V7
F = V7.F
DT = F.DT; DEG = F.DEG; LIM = F.LIM; J = 12
LAGS, LOGL, SNR_GATE, TAU_MAX, RATE, kappa_H = V7.LAGS, V7.LOGL, V7.SNR_GATE, V7.TAU_MAX, V7.RATE, V7.kappa_H
FLOOR = 1e-3 / 0.45; NDRAW = 10

def sim(kind, tt, q, tau, sd, seed, C, prm=None, floor=FLOOR, plant=F.PLANT):
    """kind 'pid' | 'uprime'; C noise draws in parallel, one noise row each (C = 1 reproduces the v5 generator)."""
    rng = np.random.default_rng(seed)
    Ts = plant['Ts']; k_inv = (plant['stiff_deg'] * DEG) / F.TAU_S; hb = 0.5 * plant['play_deg'] * DEG
    th_g = np.tile(q[0], (C, 1)); cmd = th_g.copy(); corr = np.zeros((C, J)); u = np.zeros((C, J))
    I = np.zeros((C, J)); ef = np.zeros((C, J)); efp = np.zeros((C, J))
    th_m = np.tile(q[0], (C, 1)); th_ff = q[0].copy()
    tau_c = np.full((C, J), 0.5); H = np.full((C, J), 0.5); tau_target = np.full((C, J), 0.5); a_w = DT / 1.0
    n2 = np.zeros((C, J)); r_hist = np.zeros((LAGS.max() + 1, C, J)); Dk = np.zeros((len(LAGS), C, J))
    n = len(tt); keep = []
    for k in range(n):
        load = tau[k]
        th_out = th_g - load * k_inv - hb * np.tanh(load / (0.02 * F.TAU_S))
        A_ = np.round((th_out + rng.normal(0.0, 1.0, (C, J)) * sd[k]) / F.ENC_LSB) * F.ENC_LSB
        B_ = np.round((th_out + rng.normal(0.0, 1.0, (C, J)) * sd[k]) / F.ENC_LSB) * F.ENC_LSB
        meas = 0.5 * (A_ + B_); e = q[k] - meas
        if kind == 'pid':
            ef += (DT / 0.005) * (e - ef)
            raw = prm['kp'] * e + prm['ki'] * I + prm['kd'] * (ef - efp) / DT
            sat = np.abs(raw) >= LIM
            I += np.where(sat, 0.0, e * DT); efp = ef.copy()
            corr = np.clip(raw, -LIM, LIM)
        else:
            d = (A_ - B_) / np.sqrt(2.0)
            n2 += a_w * (d * d - n2)
            r = (q[k] - th_ff)[None, :] + (meas - th_m)
            r_hist = np.roll(r_hist, 1, axis=0); r_hist[0] = r
            if k >= LAGS.max():
                for i_, lg in enumerate(LAGS):
                    dd = r - r_hist[lg]; Dk[i_] += a_w * (dd * dd - Dk[i_])
            if k > 1000:
                D = Dk - n2[None]
                ok_g = D > SNR_GATE * n2[None]
                ok_u = D > 0
                use_g = ok_g.sum(0) >= 2
                ok = np.where(use_g[None], ok_g, ok_u)
                nok = ok.sum(0)
                logD = np.where(ok, np.log(np.maximum(D, 1e-30)), 0.0)
                W = ok.astype(float); Sw = W.sum(0); Sx = (W * LOGL[:, None, None]).sum(0); Sy = (W * logD).sum(0)
                Sxx = (W * LOGL[:, None, None] ** 2).sum(0); Sxy = (W * LOGL[:, None, None] * logD).sum(0)
                den = Sw * Sxx - Sx ** 2
                slope = np.where(den > 1e-12, (Sw * Sxy - Sx * Sy) / np.maximum(den, 1e-12), 1.0)
                Hn = np.clip(0.5 * slope, 0.05, 1.0)
                valid = nok >= 2
                s2 = 0.5 * n2 * DT
                loglam = (W * (logD - 2 * Hn[None] * LOGL[:, None, None])).sum(0) / np.maximum(Sw, 1)
                lam = np.exp(loglam)
                ts = kappa_H(Hn) * (np.maximum(s2, 1e-30) / np.maximum(lam, 1e-30)) ** (1.0 / (2 * Hn + 1))
                ts = np.clip(ts, floor, TAU_MAX)
                tau_target = np.where(valid, ts, tau_target); H = np.where(valid, Hn, H)
                rr = 1.0 / tau_c; dr_max = RATE * DT / tau_c ** 2
                rr = rr + np.clip(1.0 / tau_target - rr, -dr_max, dr_max); tau_c = 1.0 / rr
            g = 1.0 - np.exp(-DT / tau_c)
            u = np.clip(u + g * e, -LIM, LIM); corr = u
        if k % F.CMD_EVERY == 0: cmd = q[k] + corr
        w_eff = F.W0 * np.maximum(0.05, 1.0 - np.abs(load) / F.TAU_S)
        th_g = th_g + DT * np.clip((cmd - th_g) / Ts, -w_eff, w_eff)
        if kind == 'uprime':
            th_m = th_m + DT * np.clip((cmd - th_m) / 0.025, -F.W0, F.W0)
            th_ff = th_ff + DT * np.clip((q[k] - th_ff) / 0.025, -F.W0, F.W0)
        if k % 10 == 0: keep.append(th_out.copy())
    return dict(th=np.array(keep), t=tt[::10][:len(keep)])

def ci(x):
    lx = np.log(np.asarray(x)); m = lx.mean(); h = 2.262 * lx.std(ddof=1) / np.sqrt(len(lx))
    return dict(gmean=float(np.exp(m)), lo=float(np.exp(m - h)), hi=float(np.exp(m + h)))

if __name__ == '__main__':
    t_start = time.time(); out = dict(floor_ms=FLOOR * 1e3, ndraw=NDRAW, check={}, results={})
    TUN = {nm: json.load(open(f'{F.OUT}/tuning_v5_pair_{nm}.json')) for nm in 'AB'}
    def gains(nm, C): return {kk: np.repeat(np.array([[TUN[nm]['pid'][j][i] for j in F.JN]]), C, 0) for i, kk in enumerate(('kp', 'ki', 'kd'))}
    # check: the v5c generator, A steady
    tt, q, tau = F.plan(6.0 - 3.0, 36.0); sd = np.full(len(tt), F.ENC_SD)
    for lab, kind, kw in (('PID', 'pid', dict(prm=gains('A', 1))), ('U floor 2 ms', 'uprime', dict(floor=0.002))):
        r = sim(kind, tt, q, tau, sd, 5, 1, **kw); k = r['t'] >= 6.0 - 1e-9
        s_, _ = F.score(r['th'][k][:, 0, :], r['t'][k], 6.0); out['check'][lab] = (s_['tip_rms_mm'], s_['jitter_mm'])
    print('CHECK A steady, seed 5:', {k: tuple(round(x, 2) for x in v) for k, v in out['check'].items()}, round(time.time() - t_start), 's', flush=True)
    for (nm, t0), cond, seed in ((('A', 6.0), 'real', 1000), (('A', 6.0), 'step', 2000), (('B', 217.5), 'real', 3000), (('B', 217.5), 'step', 4000)):
        tt, q, tau = F.plan(t0 - 3.0, t0 + 30.0)
        sd = np.full(len(tt), F.ENC_SD)
        if cond == 'step': sd[(tt >= t0 + 10) & (tt < t0 + 20)] *= np.sqrt(10.0)
        rp = sim('pid', tt, q, tau, sd, seed, NDRAW, prm=gains(nm, NDRAW)); ru = sim('uprime', tt, q, tau, sd, seed, NDRAW)
        k = rp['t'] >= t0 - 1e-9; rows = []
        for c in range(NDRAW):
            sp, _ = F.score(rp['th'][k][:, c, :], rp['t'][k], t0); su, _ = F.score(ru['th'][k][:, c, :], ru['t'][k], t0)
            rows.append(dict(pid=sp['tip_rms_mm'], u=su['tip_rms_mm'], pid_jit=sp['jitter_mm'], u_jit=su['jitter_mm']))
        rho = [x['u'] / x['pid'] for x in rows]; rj = [x['u_jit'] / x['pid_jit'] for x in rows]
        res = dict(draws=rows, rms_ratio=ci(rho), jitter_ratio=ci(rj),
                   pid_mean=float(np.mean([x['pid'] for x in rows])), u_mean=float(np.mean([x['u'] for x in rows])),
                   u_wins=int(sum(x['u'] < x['pid'] for x in rows)))
        out['results'][f'{nm} {cond}'] = res
        print(f'RESULT {nm} {cond}: PID {res["pid_mean"]:.3f} mm, U {res["u_mean"]:.3f} mm, ratio {res["rms_ratio"]["gmean"]:.3f} '
              f'[{res["rms_ratio"]["lo"]:.3f}, {res["rms_ratio"]["hi"]:.3f}], jitter ratio {res["jitter_ratio"]["gmean"]:.3f} '
              f'[{res["jitter_ratio"]["lo"]:.3f}, {res["jitter_ratio"]["hi"]:.3f}], U better in {res["u_wins"]}/10', round(time.time() - t_start), 's', flush=True)
        json.dump(out, open(f'{F.OUT}/summary_v10_seeds.json', 'w'), indent=1)
    print('done', round(time.time() - t_start), 's', flush=True)
