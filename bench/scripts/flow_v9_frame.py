"""v9: does the best window floor follow the servo frame? Written 21 Oct 2026, before any run of this file.

After v8, a post-hoc candidate: never correct more than the full error within one servo frame,
  g = 1 - exp(-dt / tau) <= dt / T_frame,  i.e.  tau_min = -dt / ln(1 - dt / T_frame).
Test: change the servo PWM frame T_frame in {2, 3, 5, 8} ms (controller and plant at 1 kHz; 2 ms is faster than the
BLS3355's 330 Hz, a simulation probe only), sweep the U' floor over {1, 1.5, 2, 3, 4, 5, 6, 8, 10} ms, and add the rules:
  per-frame rule: 1.44, 2.47, 4.48, 7.49 ms;   delay rule (v8): T_frame / 2 + dt / 2 + 0.1 ms = 1.6, 2.1, 3.1, 4.6 ms.
U' derived as in v5 / v7 (nothing tuned); every floor run in one vectorized pass, each configuration with its own reads.
Score per frame: root of the mean square of the tip rms over the four conditions (A, B; steady noise, noise x10).

Predictions:
G1  The best floor of the sweep rises with the frame (monotonic from 2 to 8 ms).
G2  The per-frame rule is within 2 % of the best floor of the sweep at every frame.
G3  A fixed 2 ms floor is not within 2 % at 5 and 8 ms (if it is, the frame does not matter and the fixed floor stands).
Check: at T_frame = 3 ms the 2 ms row reproduces v5c / v8 (0.73, 0.92, 1.28, 1.45 mm).
A failed line is reported as such."""
import os, sys, json, time, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
os.environ['BENCH_HZ'] = '1000'
import flow_v7_ctrl as V7
F = V7.F
DT = F.DT; DEG = F.DEG; LIM = F.LIM; J = 12
LAGS, LOGL, SNR_GATE, TAU_MAX, RATE, kappa_H = V7.LAGS, V7.LOGL, V7.SNR_GATE, V7.TAU_MAX, V7.RATE, V7.kappa_H
SWEEP = [1, 1.5, 2, 3, 4, 5, 6, 8, 10]

def sim(tt, q, tau, sd, seed, floors, cmd_every, plant=F.PLANT):
    rng = np.random.default_rng(seed); C = len(floors); fl = np.asarray(floors, float)[:, None]
    Ts = plant['Ts']; k_inv = (plant['stiff_deg'] * DEG) / F.TAU_S; hb = 0.5 * plant['play_deg'] * DEG
    th_g = np.tile(q[0], (C, 1)); cmd = th_g.copy(); u = np.zeros((C, J)); th_m = np.tile(q[0], (C, 1)); th_ff = q[0].copy()
    tau_c = np.full((C, J), 0.5); H = np.full((C, J), 0.5); tau_target = np.full((C, J), 0.5); a_w = DT / 1.0
    n2 = np.zeros((C, J)); r_hist = np.zeros((LAGS.max() + 1, C, J)); Dk = np.zeros((len(LAGS), C, J))
    n = len(tt); keep, ktau = [], []
    for k in range(n):
        load = tau[k]
        th_out = th_g - load * k_inv - hb * np.tanh(load / (0.02 * F.TAU_S))
        A_ = np.round((th_out + rng.normal(0.0, 1.0, J) * sd[k]) / F.ENC_LSB) * F.ENC_LSB
        B_ = np.round((th_out + rng.normal(0.0, 1.0, J) * sd[k]) / F.ENC_LSB) * F.ENC_LSB
        meas = 0.5 * (A_ + B_); e = q[k] - meas
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
            ts = np.clip(ts, fl, TAU_MAX)
            tau_target = np.where(valid, ts, tau_target); H = np.where(valid, Hn, H)
            rr = 1.0 / tau_c; dr_max = RATE * DT / tau_c ** 2
            rr = rr + np.clip(1.0 / tau_target - rr, -dr_max, dr_max); tau_c = 1.0 / rr
        g = 1.0 - np.exp(-DT / tau_c)
        u = np.clip(u + g * e, -LIM, LIM)
        if k % cmd_every == 0: cmd = q[k] + u
        w_eff = F.W0 * np.maximum(0.05, 1.0 - np.abs(load) / F.TAU_S)
        th_g = th_g + DT * np.clip((cmd - th_g) / Ts, -w_eff, w_eff)
        th_m = th_m + DT * np.clip((cmd - th_m) / 0.025, -F.W0, F.W0)
        th_ff = th_ff + DT * np.clip((q[k] - th_ff) / 0.025, -F.W0, F.W0)
        if k % 10 == 0:
            keep.append(th_out.copy()); ktau.append(tau_c.copy())
    return dict(th=np.array(keep), t=tt[::10][:len(keep)], tau=np.array(ktau))

if __name__ == '__main__':
    t_start = time.time(); out = dict(sweep_ms=SWEEP, frames={})
    for T in (3, 2, 5, 8):
        rule_f = -DT / np.log(1 - DT / (T * 1e-3)); rule_d = T * 1e-3 / 2 + DT / 2 + 1e-4
        labs = [f'{x:g} ms' for x in SWEEP] + ['frame rule', 'delay rule']
        floors = [x * 1e-3 for x in SWEEP] + [rule_f, rule_d]
        res = {l: {} for l in labs}
        for name, t0 in (('A', 6.0), ('B', 217.5)):
            for cond in ('real', 'step'):
                tt, q, tau = F.plan(t0 - 3.0, t0 + 30.0)
                sd = np.full(len(tt), F.ENC_SD)
                if cond == 'step': sd[(tt >= t0 + 10) & (tt < t0 + 20)] *= np.sqrt(10.0)
                seed = 7 if cond == 'step' else 5
                r = sim(tt, q, tau, sd, seed, floors, T); k = r['t'] >= t0 - 1e-9
                for c, l in enumerate(labs):
                    s_, _ = F.score(r['th'][k][:, c, :], r['t'][k], t0)
                    res[l][f'{name} {cond}'] = dict(rms=s_['tip_rms_mm'], jit=s_['jitter_mm'],
                                                     on_floor=float(np.mean(r['tau'][k][:, c, :6] <= floors[c] * 1.001)))
                print(f'T {T} ms {name} {cond} done', round(time.time() - t_start), 's', flush=True)
        score = {l: float(np.sqrt(np.mean([v['rms'] ** 2 for v in res[l].values()]))) for l in labs}
        best = min(SWEEP, key=lambda x: score[f'{x:g} ms'])
        out['frames'][str(T)] = dict(rule_frame_ms=rule_f * 1e3, rule_delay_ms=rule_d * 1e3, score=score, best_sweep_ms=best, results=res)
        print(f'RESULT T {T} ms: best sweep floor {best} ms ({score[f"{best:g} ms"]:.3f}); frame rule {rule_f * 1e3:.2f} ms ({score["frame rule"]:.3f}); '
              f'delay rule {rule_d * 1e3:.2f} ms ({score["delay rule"]:.3f}); fixed 2 ms ({score["2 ms"]:.3f}); 2 ms row by condition',
              {c: round(v['rms'], 2) for c, v in res['2 ms'].items()}, flush=True)
        json.dump(out, open(f'{F.OUT}/summary_v9_frame.json', 'w'), indent=1)
    print('done', round(time.time() - t_start), 's', flush=True)
