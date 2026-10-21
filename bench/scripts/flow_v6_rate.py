"""v6: the v5 joint bench (summary_v5c) run at another loop rate. Written 21 Oct 2026, before any run of this file.

Question: does a 2 kHz loop (encoder pair reads and engines every 0.5 ms) change U' against a PID re-tuned at 2 kHz?

What changes with the rate, and only that:
- dt = 1 / BENCH_HZ for the encoder reads, every engine and the plant integration;
- the servo command stays at 333 Hz (one command every 3 ms: CMD_EVERY = 3 ms / dt);
- the encoder noise per read is unchanged (0.06 deg rms, 0.0219 deg steps): the AS5048A refreshes its angle at
  10.2-12.4 kHz, so pair reads 0.25 ms apart still carry independent noise (data sheet DS000298);
- every time constant stays in seconds, as v5 declares it: lags {2, 4, 8, 16, 32, 64} ms, EMA 1 s, first use after 1 s,
  tau path 0.4 s per second, PID derivative filter 5 ms, servo 25 ms, scores at 100 Hz;
- the window floor follows the bench's own rule at 1 kHz, tau_min = 2 dt (2 ms at 1 kHz, 1 ms at 2 kHz); a second U'
  row keeps the floor at 2 ms, to separate the effect of the floor from the effect of the noise per unit time.
Baselines: PID + interpolation and best fixed window re-tuned at the new rate on the v5 grids, same spans and seeds.
U' derived: the v5 code, copied below with only dt-expressed constants; nothing tuned (kappa_m = 1).
Check: at 1000 Hz with the cached v5 tunings this file must reproduce summary_v5c to the digit.

Predictions, written before the first 2 kHz run:
Q1  Every engine improves a little at 2 kHz (half the noise per unit time).
Q2  U' gains most on excerpt B, real noise, where it sat 54 % of the time on the 2 ms floor at 1 kHz and lost to PID
    (1.28 against 1.19 mm); the floor-2 ms row recovers less of that gain than the floor-1 ms row.
Q3  The U' window median falls by about 2^(-1/3) = 0.79 where it is not on a bound.
A failed line is reported as such."""
import os, sys, json, time, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
HZ = int(os.environ.get('BENCH_HZ', '2000'))
import flow_vs_pid as F
import flow_v3 as V3
import flow_v5_uprime as V5

DT = 1.0 / HZ
F.DT = DT; V3.DT = DT; V5.DT = DT
F.CMD_EVERY = int(round(0.003 / DT))
DEG = F.DEG; LIM = F.LIM; J = 12; OUT = F.OUT
LAGS = np.round(np.array([2, 4, 8, 16, 32, 64]) * 1e-3 / DT).astype(int); LOGL = np.log(LAGS * DT)
WARM = int(round(1.0 / DT)); SE = int(round(0.01 / DT)); DEC = int(round(SE / 10))
SNR_GATE = V5.SNR_GATE; TAU_MAX, RATE = 0.5, 0.4
kappa_H = V5.kappa_H

def simulate_u(km, tt, q, tau, sd, seed, tau_min, plant=F.PLANT):
    """v5 simulate, U' chart, store=True; constants expressed in dt."""
    rng = np.random.default_rng(seed); C = km.shape[0]
    Ts = plant['Ts']; k_inv = (plant['stiff_deg'] * DEG) / F.TAU_S; hb = 0.5 * plant['play_deg'] * DEG
    th_g = np.tile(q[0], (C, 1)); cmd = th_g.copy(); u = np.zeros((C, J)); th_m = np.tile(q[0], (C, 1)); th_ff = q[0].copy()
    tau_c = np.full((C, J), 0.5); H = np.full((C, J), 0.5); tau_target = np.full((C, J), 0.5)
    a_w = DT / 1.0
    n2 = np.zeros(J); r_hist = np.zeros((LAGS.max() + 1, C, J)); Dk = np.zeros((len(LAGS), C, J))
    err2 = np.zeros((C, J)); n = len(tt); keep, ktau, kH = [], [], []
    for k in range(n):
        load = tau[k]
        th_out = th_g - load * k_inv - hb * np.tanh(load / (0.02 * F.TAU_S))
        A_ = np.round((th_out + rng.normal(0.0, 1.0, J) * sd[k]) / F.ENC_LSB) * F.ENC_LSB
        B_ = np.round((th_out + rng.normal(0.0, 1.0, J) * sd[k]) / F.ENC_LSB) * F.ENC_LSB
        meas = 0.5 * (A_ + B_); e = q[k] - meas
        d = (A_[0] - B_[0]) / np.sqrt(2.0)
        n2 += a_w * (d * d - n2)
        r = (q[k] - th_ff)[None, :] + (meas - th_m)
        r_hist = np.roll(r_hist, 1, axis=0); r_hist[0] = r
        if k >= LAGS.max():
            for i_, lg in enumerate(LAGS):
                dd = r - r_hist[lg]; Dk[i_] += a_w * (dd * dd - Dk[i_])
        if k > WARM:
            D = Dk - n2[None, None, :]
            ok_g = D > SNR_GATE * n2[None, None, :]
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
            ts = np.clip(km * ts, tau_min, TAU_MAX)
            tau_target = np.where(valid, ts, tau_target); H = np.where(valid, Hn, H)
            rr = 1.0 / tau_c; dr_max = RATE * DT / tau_c ** 2
            rr = rr + np.clip(1.0 / tau_target - rr, -dr_max, dr_max); tau_c = 1.0 / rr
        g = 1.0 - np.exp(-DT / tau_c)
        u = np.clip(u + g * e, -LIM, LIM)
        if k % F.CMD_EVERY == 0: cmd = q[k] + u
        w_eff = F.W0 * np.maximum(0.05, 1.0 - np.abs(load) / F.TAU_S)
        th_g = th_g + DT * np.clip((cmd - th_g) / Ts, -w_eff, w_eff)
        th_m = th_m + DT * np.clip((cmd - th_m) / 0.025, -F.W0, F.W0)
        th_ff = th_ff + DT * np.clip((q[k] - th_ff) / 0.025, -F.W0, F.W0)
        err2 += (q[k] - th_out) ** 2
        if k % SE == 0:
            keep.append(th_out[0].copy()); ktau.append(tau_c[0].copy()); kH.append(H[0].copy())
    return dict(mse=err2 / n, th=np.array(keep), t=tt[::SE][:len(keep)], tau=np.array(ktau), H=np.array(kH))

def base(kind, prm, tt, q, tau, sd, seed):
    r = V3.simulate(kind, prm, tt, q, tau, sd, seed, F.PLANT, True)
    r['th'] = r['th'][::DEC]; r['t'] = r['t'][::DEC]; return r

def run(t0, tun, cond, tag):
    tt, q, tau = F.plan(t0 - 3.0, t0 + 30.0)
    sd = np.full(len(tt), F.ENC_SD)
    if cond == 'step': sd[(tt >= t0 + 10) & (tt < t0 + 20)] *= np.sqrt(10.0)
    seed = 7 if cond == 'step' else 5
    pid = {kk: np.array([[tun['pid'][j][i] for j in F.JN]]) for i, kk in enumerate(('kp', 'ki', 'kd'))}
    b0 = dict(tau=np.array([[tun['b0'][j] for j in F.JN]]))
    one = np.ones((1, J))
    res = {'PID + interpolation': base('pid', pid, tt, q, tau, sd, seed),
           'best fixed window': base('b0', b0, tt, q, tau, sd, seed),
           'SYN U-prime derived': simulate_u(one, tt, q, tau, sd, seed, 2 * DT)}
    floors = {'SYN U-prime derived': 2 * DT}
    if HZ != 1000:
        res['SYN U-prime derived, floor 2 ms'] = simulate_u(one, tt, q, tau, sd, seed, 0.002); floors['SYN U-prime derived, floor 2 ms'] = 0.002
    sc = {}
    for name, r in res.items():
        k = r['t'] >= t0 - 1e-9
        s_, e = F.score(r['th'][k], r['t'][k], t0); sc[name] = s_
        if 'tau' in r:
            t_ = r['t'][k]; tau_ = r['tau'][k][:, :6]
            parts = {p: (t_ >= t0 + 10 * i) & (t_ < t0 + 10 * (i + 1)) for i, p in enumerate(('0-10', '10-20', '20-30'))}
            s_['tau_median_ms'] = {p: float(1e3 * np.median(tau_[m_])) for p, m_ in parts.items()}
            s_['tau_step_factor'] = s_['tau_median_ms']['10-20'] / s_['tau_median_ms']['0-10']
            s_['tau_floor_frac'] = float(np.mean(tau_ <= floors[name] * 1.001))
            s_['tau_edge_frac'] = float(np.mean((tau_ <= floors[name] * 1.001) | (tau_ >= TAU_MAX * 0.999)))
            s_['th41a_inside_frac'] = float(np.mean(DT / tau_ <= 0.45))
            H_ = r['H'][k][:, :6]; s_['H_median_arm'] = float(np.median(H_))
            np.savez_compressed(f'{OUT}/v6_{HZ}Hz_{tag}_{cond}_{name.replace(" ", "_").replace(",", "")}.npz', th=r['th'][k], t=t_, err=e, tau=r['tau'][k], H=r['H'][k])
    return sc

if __name__ == '__main__':
    t_start = time.time(); A, B = 6.0, 217.5
    F.PID_GRID = [(kp, ki, kd) for kp in (0, 0.5, 1, 2, 3, 4, 6, 8) for ki in (0, 80, 160, 320, 640, 1280, 2560) for kd in (0, 0.0005, 0.001, 0.002, 0.005)]
    F.B0_GRID = [0.001, 0.0015, 0.002, 0.003, 0.005, 0.007, 0.01, 0.015, 0.02, 0.03, 0.05, 0.07, 0.1, 0.15, 0.2, 0.3, 0.5]
    G = np.array(F.PID_GRID)
    summary = dict(hz=HZ, dt=DT, cmd_every=F.CMD_EVERY, lags_samples=LAGS.tolist(), tau_min_ms=2e3 * DT, results={}, tuning={})
    print(f'{HZ} Hz: dt {DT * 1e3:.3f} ms, command every {F.CMD_EVERY} steps, lags {LAGS.tolist()} samples', flush=True)
    for name, t0 in (('A', A), ('B', B)):
        tag = f'v5_pair_{name}' if HZ == 1000 else f'v6_{HZ}Hz_pair_{name}'
        tun = V3.tune(t0, tag=tag)
        edges = {}
        for j in F.JN:
            p = tun['pid'][j]; e_ = [ax for ax, v in zip(('kp', 'ki', 'kd'), p) if v == G[:, ['kp', 'ki', 'kd'].index(ax)].max()]
            if tun['b0'][j] in (F.B0_GRID[0], F.B0_GRID[-1]): e_.append('b0')
            if e_: edges[j] = e_
        summary['tuning'][name] = dict(tun, top_edges=edges)
        print(name, 'baselines tuned', round(time.time() - t_start), 's; at top of grid:', edges, flush=True)
        for cond in ('real', 'step'):
            sc = run(t0, tun, cond, name); summary['results'][f'{name} {cond}'] = sc
            print(name, cond, json.dumps({kk: dict(rms=round(v['tip_rms_mm'], 2), jit=round(v['jitter_mm'], 2),
                                                    tau=v.get('tau_median_ms'), floor=round(v['tau_floor_frac'], 3) if 'tau_floor_frac' in v else None,
                                                    th41a=round(v['th41a_inside_frac'], 3) if 'th41a_inside_frac' in v else None,
                                                    H=round(v['H_median_arm'], 2) if 'H_median_arm' in v else None)
                                          for kk, v in sc.items()}, indent=0), round(time.time() - t_start), 's', flush=True)
            json.dump(summary, open(f'{OUT}/summary_v6_{HZ}Hz.json', 'w'), indent=1)
    print('done', round(time.time() - t_start), 's')
