"""v7: the v5 joint bench with the controller slower than the world. Written 21 Oct 2026, before any run of this file.

Question: what changes if the loop (encoder pair reads + engine) runs at 500 Hz, or 333 Hz, instead of 1 kHz?

Model. The plant (servo lag 25 ms, speed derating, compliance, play, loads) is integrated at 1 kHz exactly as in v5. The
controller reads the encoder pair and updates its engine every R plant steps (R = 1000 / HZ: 2 at 500 Hz, 3 at 333 Hz).
The servo PWM frame stays at 333 Hz (every 3 ms) and carries the plan at the frame time plus the latest correction
(asynchronous PWM, as on the Teensy). Engine constants in controller time, dt_c = R ms: EMA 1 s; lags {2, 4, 8, 16, 32,
64} ms rounded to whole ticks (at least one tick, duplicates dropped); first use after 1 s; tau path 0.4 s per second;
PID derivative filter 5 ms; the engine's own servo models integrated at the controller ticks. Noise per read unchanged.
Floors: the bench rule tau_min = 2 dt_c (4 ms at 500 Hz, 6 ms at 333 Hz), and a row with the 2 ms floor of the 1 kHz bench.
PID + interpolation and best fixed window re-tuned at the controller rate on the v5 grids, same spans and seeds.
Check: HZ = 1000 (R = 1) with the cached v5 tunings must reproduce summary_v5c.

Predictions, written after v6 (best loop window 1.5-2 ms, set by the servo rather than by the noise) and before this run:
S1  PID and the fixed window get a little worse at 500 Hz (fewer reads, older corrections).
S2  U' with the 2 dt_c floor (4 ms) sits on that floor most of the time and loses more than PID on steady noise.
S3  U' with the 2 ms floor stays close to its 1 kHz result; its window lengthens by about 2^(1/3) = 1.26 off the bounds.
A failed line is reported as such."""
import os, sys, json, time, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
HZ = int(os.environ.get('BENCH_HZ', '500'))
import flow_vs_pid as F
import flow_v5_uprime as V5
assert F.DT == 1e-3 and F.CMD_EVERY == 3
R = int(round(1000 / HZ)); DTP = F.DT; DTC = R * DTP
DEG = F.DEG; LIM = F.LIM; J = 12; OUT = F.OUT
LAGS = np.unique(np.maximum(1, np.round(np.array([2, 4, 8, 16, 32, 64]) * 1e-3 / DTC))).astype(int); LOGL = np.log(LAGS * DTC)
SNR_GATE = V5.SNR_GATE; TAU_MAX, RATE = 0.5, 0.4; kappa_H = V5.kappa_H

def simulate(kind, prm, tt, q, tau, sd, seed, tau_min=None, plant=F.PLANT, store=False):
    """kind 'pid' | 'b0' (vectorized over configs) | 'uprime' (prm = dict(km=(C, 12)))."""
    rng = np.random.default_rng(seed)
    C = next(iter(prm.values())).shape[0]
    Ts = plant['Ts']; k_inv = (plant['stiff_deg'] * DEG) / F.TAU_S; hb = 0.5 * plant['play_deg'] * DEG
    th_g = np.tile(q[0], (C, 1)); cmd = th_g.copy(); corr = np.zeros((C, J))
    I = np.zeros((C, J)); ef = np.zeros((C, J)); efp = np.zeros((C, J)); u = np.zeros((C, J))
    if kind == 'b0': g_b0 = 1.0 - np.exp(-DTC / prm['tau'])
    if kind == 'uprime':
        km = prm['km']; th_m = np.tile(q[0], (C, 1)); th_ff = q[0].copy()
        tau_c = np.full((C, J), 0.5); H = np.full((C, J), 0.5); tau_target = np.full((C, J), 0.5); a_w = DTC / 1.0
        n2 = np.zeros(J); r_hist = np.zeros((LAGS.max() + 1, C, J)); Dk = np.zeros((len(LAGS), C, J)); kc = 0
    err2 = np.zeros((C, J)); n = len(tt); keep, ktau, kH = [], [], []
    for k in range(n):
        load = tau[k]
        th_out = th_g - load * k_inv - hb * np.tanh(load / (0.02 * F.TAU_S))
        tick = (k % R == 0)
        if tick:
            A_ = np.round((th_out + rng.normal(0.0, 1.0, J) * sd[k]) / F.ENC_LSB) * F.ENC_LSB
            B_ = np.round((th_out + rng.normal(0.0, 1.0, J) * sd[k]) / F.ENC_LSB) * F.ENC_LSB
            meas = 0.5 * (A_ + B_); e = q[k] - meas
            if kind == 'pid':
                ef += (DTC / 0.005) * (e - ef)
                raw = prm['kp'] * e + prm['ki'] * I + prm['kd'] * (ef - efp) / DTC
                sat = np.abs(raw) >= LIM
                I += np.where(sat, 0.0, e * DTC); efp = ef.copy()
                corr = np.clip(raw, -LIM, LIM)
            elif kind == 'b0':
                u = np.clip(u + g_b0 * e, -LIM, LIM); corr = u
            else:
                d = (A_[0] - B_[0]) / np.sqrt(2.0)
                n2 += a_w * (d * d - n2)
                r = (q[k] - th_ff)[None, :] + (meas - th_m)
                r_hist = np.roll(r_hist, 1, axis=0); r_hist[0] = r
                if kc >= LAGS.max():
                    for i_, lg in enumerate(LAGS):
                        dd = r - r_hist[lg]; Dk[i_] += a_w * (dd * dd - Dk[i_])
                if k > 1000:
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
                    s2 = 0.5 * n2 * DTC
                    loglam = (W * (logD - 2 * Hn[None] * LOGL[:, None, None])).sum(0) / np.maximum(Sw, 1)
                    lam = np.exp(loglam)
                    ts = kappa_H(Hn) * (np.maximum(s2, 1e-30) / np.maximum(lam, 1e-30)) ** (1.0 / (2 * Hn + 1))
                    ts = np.clip(km * ts, tau_min, TAU_MAX)
                    tau_target = np.where(valid, ts, tau_target); H = np.where(valid, Hn, H)
                    rr = 1.0 / tau_c; dr_max = RATE * DTC / tau_c ** 2
                    rr = rr + np.clip(1.0 / tau_target - rr, -dr_max, dr_max); tau_c = 1.0 / rr
                g = 1.0 - np.exp(-DTC / tau_c)
                u = np.clip(u + g * e, -LIM, LIM); corr = u
                kc += 1
        if k % F.CMD_EVERY == 0: cmd = q[k] + corr
        w_eff = F.W0 * np.maximum(0.05, 1.0 - np.abs(load) / F.TAU_S)
        th_g = th_g + DTP * np.clip((cmd - th_g) / Ts, -w_eff, w_eff)
        if kind == 'uprime' and tick:
            th_m = th_m + DTC * np.clip((cmd - th_m) / 0.025, -F.W0, F.W0)
            th_ff = th_ff + DTC * np.clip((q[k] - th_ff) / 0.025, -F.W0, F.W0)
        err2 += (q[k] - th_out) ** 2
        if store and k % 10 == 0:
            keep.append(th_out[0].copy())
            if kind == 'uprime': ktau.append(tau_c[0].copy()); kH.append(H[0].copy())
    out = dict(mse=err2 / n)
    if store:
        out.update(th=np.array(keep), t=tt[::10][:len(keep)])
        if kind == 'uprime': out.update(tau=np.array(ktau), H=np.array(kH))
    return out

def tune(t0, tag):
    f = f'{OUT}/tuning_{tag}.json'
    if os.path.exists(f): return json.load(open(f))
    G = np.array(F.PID_GRID); C = len(G)
    prm = dict(kp=np.repeat(G[:, :1], 12, 1), ki=np.repeat(G[:, 1:2], 12, 1), kd=np.repeat(G[:, 2:], 12, 1))
    tb = np.array(F.B0_GRID); prm_b = dict(tau=np.repeat(tb[:, None], 12, 1))
    m_pid = np.zeros((C, 12)); m_b0 = np.zeros((len(tb), 12)); ntot = 0
    for s_, (a, b) in enumerate(F.training_spans(t0)):
        if b - a < 5: continue
        tt, q, tau = F.plan(a, b); sd = np.full(len(tt), F.ENC_SD); w = len(tt)
        m_pid += w * simulate('pid', prm, tt, q, tau, sd, 100 + s_)['mse']
        m_b0 += w * simulate('b0', prm_b, tt, q, tau, sd, 100 + s_)['mse']; ntot += w
    ip = np.argmin(m_pid / ntot, 0); ib = np.argmin(m_b0 / ntot, 0)
    res = dict(pid={F.JN[j]: list(map(float, G[ip[j]])) for j in range(12)}, b0={F.JN[j]: float(tb[ib[j]]) for j in range(12)})
    json.dump(res, open(f, 'w'), indent=1); return res

def run(t0, tun, cond, tag):
    tt, q, tau = F.plan(t0 - 3.0, t0 + 30.0)
    sd = np.full(len(tt), F.ENC_SD)
    if cond == 'step': sd[(tt >= t0 + 10) & (tt < t0 + 20)] *= np.sqrt(10.0)
    seed = 7 if cond == 'step' else 5
    pid = {kk: np.array([[tun['pid'][j][i] for j in F.JN]]) for i, kk in enumerate(('kp', 'ki', 'kd'))}
    b0 = dict(tau=np.array([[tun['b0'][j] for j in F.JN]]))
    one = dict(km=np.ones((1, J)))
    res = {'PID + interpolation': simulate('pid', pid, tt, q, tau, sd, seed, store=True),
           'best fixed window': simulate('b0', b0, tt, q, tau, sd, seed, store=True),
           'SYN U-prime derived': simulate('uprime', one, tt, q, tau, sd, seed, 2 * DTC, store=True)}
    floors = {'SYN U-prime derived': 2 * DTC}
    if R != 1:
        res['SYN U-prime derived, floor 2 ms'] = simulate('uprime', one, tt, q, tau, sd, seed, 0.002, store=True); floors['SYN U-prime derived, floor 2 ms'] = 0.002
    sc = {}
    for name, r in res.items():
        k = r['t'] >= t0 - 1e-9
        s_, e = F.score(r['th'][k], r['t'][k], t0); sc[name] = s_
        if 'tau' in r:
            t_ = r['t'][k]; tau_ = r['tau'][k][:, :6]
            parts = {p: (t_ >= t0 + 10 * i) & (t_ < t0 + 10 * (i + 1)) for i, p in enumerate(('0-10', '10-20', '20-30'))}
            s_['tau_median_ms'] = {p: float(1e3 * np.median(tau_[m_])) for p, m_ in parts.items()}
            s_['tau_floor_frac'] = float(np.mean(tau_ <= floors[name] * 1.001))
            s_['th41a_inside_frac'] = float(np.mean(DTC / tau_ <= 0.45))
            s_['H_median_arm'] = float(np.median(r['H'][k][:, :6]))
            np.savez_compressed(f'{OUT}/v7_{HZ}Hz_{tag}_{cond}_{name.replace(" ", "_").replace(",", "")}.npz', th=r['th'][k], t=t_, err=e, tau=r['tau'][k], H=r['H'][k])
    return sc

if __name__ == '__main__':
    t_start = time.time(); A, B = 6.0, 217.5
    # BENCH_WIDE=1 (amendment after the first 500 / 333 Hz pass, where kd sat at the top and the fixed window at the
    # bottom of the v5 grids): kd extended to 0.05 s, fixed window down to 0.2 ms; U' rows unchanged
    WIDE = os.environ.get('BENCH_WIDE') == '1'; W_ = 'w' if WIDE else ''
    KD = (0, 0.0005, 0.001, 0.002, 0.005) + ((0.01, 0.02, 0.05) if WIDE else ())
    F.PID_GRID = [(kp, ki, kd) for kp in (0, 0.5, 1, 2, 3, 4, 6, 8) for ki in (0, 80, 160, 320, 640, 1280, 2560) for kd in KD]
    F.B0_GRID = ([0.0002, 0.0005] if WIDE else []) + [0.001, 0.0015, 0.002, 0.003, 0.005, 0.007, 0.01, 0.015, 0.02, 0.03, 0.05, 0.07, 0.1, 0.15, 0.2, 0.3, 0.5]
    G = np.array(F.PID_GRID)
    summary = dict(hz=HZ, wide=WIDE, ctrl_every_ms=DTC * 1e3, lags_ticks=LAGS.tolist(), tau_min_ms=2e3 * DTC, results={}, tuning={})
    print(f'{HZ} Hz controller: tick {DTC * 1e3:.0f} ms, PWM frame 3 ms, lags {LAGS.tolist()} ticks', flush=True)
    for name, t0 in (('A', A), ('B', B)):
        tag = f'v5_pair_{name}' if R == 1 and not WIDE else f'v7{W_}_{HZ}Hz_pair_{name}'
        tun = tune(t0, tag)
        edges = {}
        for j in F.JN:
            p = tun['pid'][j]; e_ = [ax for i, (ax, v) in enumerate(zip(('kp', 'ki', 'kd'), p)) if v == G[:, i].max()]
            if tun['b0'][j] == F.B0_GRID[-1]: e_.append('b0 top')
            if tun['b0'][j] == F.B0_GRID[0]: e_.append('b0 bottom')
            if e_: edges[j] = e_
        summary['tuning'][name] = dict(tun, top_edges=edges)
        print(name, 'baselines tuned', round(time.time() - t_start), 's; at top of grid:', edges, flush=True)
        for cond in ('real', 'step'):
            sc = run(t0, tun, cond, name); summary['results'][f'{name} {cond}'] = sc
            print(f'RESULT {name} {cond}', {kk: (round(v['tip_rms_mm'], 2), round(v['jitter_mm'], 2)) for kk, v in sc.items()}, round(time.time() - t_start), 's', flush=True)
            json.dump(summary, open(f'{OUT}/summary_v7{W_}_{HZ}Hz.json', 'w'), indent=1)
    print('done', round(time.time() - t_start), 's')
