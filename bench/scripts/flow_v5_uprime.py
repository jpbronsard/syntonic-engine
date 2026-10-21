"""v5: the syntonic window at MAESTRO's joints in two charts, canonical (Theorem A) and U' (measured H), against the
tuned baselines. Written 19 Oct 2026, after the re-validation that separated the exponent (set by the smoothness of the
tracked signal) from the kappa drift (an artefact of reading a smooth-signal loop in the canonical chart).

Everything below was fixed before any run of this file.

Plant, plan, pairing, excerpts, conditions, scores: as flow_v3 (pair reads each 1 ms, every engine on the pair mean;
excerpts A = 0:06-0:36 and B = 3:37-4:07 held out; real noise and the noise-variance x10 step from 10 to 20 s; baton
tip rms and jitter above 5 Hz from forward kinematics). The plan is regenerated here from the delivered solver
(conduct_big.py limit) because the lean archive carried no video/; loads and tunings are recomputed on it.

One measurement layer, two charts. Per joint, from the pair reads only (nothing from the engine's own state):
- noise: n2 = EMA_1s(d^2), d = (A - B)/sqrt(2), the per-read encoder noise variance; the pair-mean channel has
  variance n2/2 and noise PSD sigma2 = (n2/2) dt;
- level: the offset the integral loop has to cancel, taken from exogenous channels: L = (q - theta_FF) + r, where
  theta_FF is the nominal servo model (25 ms lag, unloaded speed limit) driven by the PLAN alone (no correction), so
  q - theta_FF is the lag offset the plan imposes (exogenous, as DIFF-EST's actuator-model residual on the quad), and
  r = pair mean - theta_model is the measured residual of the same model driven by the actual commands (compliance,
  play, speed derating: the unmodelled part). The smoke test of the first version, which used r alone as in flow_v3,
  read a level a thousand times too quiet for the loop (tau 14-20 ms where the best window is 2 ms): the lag offset
  is what the correction cancels, and a residual that excludes it by construction cannot set the window.
  Structure function at lags k in {2, 4, 8, 16, 32} ms: D(k) = EMA_1s((L_t - L_{t-k})^2) - n2 (noise floor of a
  difference of two pair means);
- memory index: 2H = slope of log D(k) against log(k dt) over the lags with D(k) > 0 (at least two), H clipped to
  [0.05, 1]; intensity lambda_H = geometric mean over those lags of D(k) / (k dt)^(2H), in deg^2 / s^(2H).
Canonical chart (Theorem A):  tau = kappa_m sqrt(sigma2 / lambda_half), lambda_half = mean D(k) / (k dt)  (H forced to 1/2).
U' chart (Theorem U', white noise, gamma = 0):  tau = kappa_m kappa_H (sigma2 / lambda_H)^(1 / (2H + 1)),
   with the kernel prefactor DERIVED for the exponential kernel on a signal of structure function lambda_H u^(2H):
   tracking term (Gamma(2H+1)/2) lambda_H tau^(2H), filtering term sigma2 / (2 tau), hence
   kappa_H = (2H Gamma(2H+1))^(-1/(2H+1)); kappa_H = 1 at H = 1/2 and 0.63 at H = 1. No constant is tuned.
Both charts: tau bounded to [2 ms, 500 ms] (H4), moved along the rate-linear path at most 0.4 s per second, first used
after 1 s of data; correction law u += (1 - exp(-dt/tau)) e on the pair-mean error, limited to +-5 deg, command at 330 Hz.
kappa_m = 1 for the two DERIVED engines. The oracle variants choose kappa_m per joint from
{0.3, 0.45, 0.6, 0.8, 1, 1.25, 1.6, 2} on the rest of the movement with clean noise (same spans and seeds as every
baseline). They are TUNED and labelled so; they measure how much a loop-lag correction, not derived in the corpus, is worth.

Predictions, written before the first run:
P1  Under the x10 noise step, the U' engine lengthens its window by about 10^(1/(2H+1)), i.e. 2.15x at H = 1; the
    canonical engine by 3.16x. (Median tau over the step part against the part before it, baton arm.)
P2  Tip error on the excerpt: the derived U' engine within 1.10x of the best fixed window in rms (1.21x in mean square),
    real noise, both excerpts; the derived canonical engine stays at 1.2-1.4x in rms as this morning.
P3  The measured H on the model residual is above 0.75 on the baton arm's joints (the plan-driven offset is smooth).
P4  The oracle kappa for the U' chart sits in [0.6, 1.2] on most joints; for the canonical chart below 0.6 (the drift).
A failed line is reported as such.

Amended after the first pass on excerpt A (so NOT fixed before seeing results; the first-pass summary is kept as
summary_v5_pass1.json). First pass: real noise, both derived engines 0.74 mm against 0.86 mm for the best fixed window
and 0.73 mm for the PID; noise step, canonical 0.93 mm and U' 1.09 mm against 1.45 mm (fixed) and 1.33 mm (PID). But
the window response to the step came out 2.2x for the canonical chart and 2.7x for U' (predicted 3.16x and 2.15x): the
H instrument broke under the x10 noise (H fell to 0.3-0.5 on RJ3 and RJ5 during the step, and already flipped between
0.5 and 1.0 on the slow joints RJ2 and RJ3 at x1 noise depending on the noise draw), because the short lags fall under
the noise floor and were still used in the slope fit. On the joints where H held near 1 (RJ1, RJ4, RJ6) the U' factor
was 2.0-2.4, as predicted. Amendment, declared: a lag is used in the H and lambda fits only when its structure
function stands above the noise floor, D(k) - n2 > 2 n2 (signal-to-floor ratio 2), and the lag set is extended to
{2, 4, 8, 16, 32, 64} ms so that at least two lags survive at x10 noise on the slow joints; with fewer than two, the
previous H and target are held. Nothing else changed; the baselines and the oracle sweeps are recomputed on the
amended instrument.

Second amendment (after the gated instrument's first run on excerpt A, real noise; that run is kept as
summary_v5_pass2_A_real_only.json). With the gate alone, every syntonic engine lost 5.3 mm in the first ten seconds of
excerpt A (0.6 mm before): in the quiet opening bars no lag passes the gate on several joints, the target is held at
its initial 500 ms, and when the loud motif arrives the rate-linear path (0.4 s per second) needs over a second to
bring the window down. Amendment, declared: the gated fit is primary; when fewer than two lags pass the gate, the
ungated fit of the first pass (lags with D(k) > 0) is used as a fallback; only when that too has fewer than two lags
is the previous target held. The initial target is still 500 ms, as in every earlier version."""
import os, sys, json, time, numpy as np
from math import gamma as Gamma
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import flow_vs_pid as F
import flow_v3 as V3

DT = F.DT; OUT = F.OUT; DEG = F.DEG; LIM = F.LIM; J = 12
LAGS = np.array([2, 4, 8, 16, 32, 64]); LOGL = np.log(LAGS * DT)
SNR_GATE = 2.0
KM_GRID = [0.3, 0.45, 0.6, 0.8, 1.0, 1.25, 1.6, 2.0]
TAU_MIN, TAU_MAX, RATE = 0.002, 0.5, 0.4

def kappa_H(H):
    return (2 * H * np.vectorize(Gamma)(2 * H + 1)) ** (-1.0 / (2 * H + 1))

def simulate(chart, km, tt, q, tau, sd, seed, plant=F.PLANT, store=False):
    """chart: 'canon' | 'uprime'. km: kappa multipliers, shape (C, 12). Vectorized over C."""
    rng = np.random.default_rng(seed); C = km.shape[0]
    Ts = plant['Ts']; k_inv = (plant['stiff_deg'] * DEG) / F.TAU_S; hb = 0.5 * plant['play_deg'] * DEG
    th_g = np.tile(q[0], (C, 1)); cmd = th_g.copy(); u = np.zeros((C, J)); th_m = np.tile(q[0], (C, 1)); th_ff = q[0].copy()
    tau_c = np.full((C, J), 0.5); H = np.full((C, J), 0.5); tau_target = np.full((C, J), 0.5)
    a_w = DT / 1.0
    n2 = np.zeros(J); r_hist = np.zeros((LAGS.max() + 1, C, J)); Dk = np.zeros((len(LAGS), C, J))
    err2 = np.zeros((C, J)); n = len(tt); keep, ktau, kH, klam, ks2 = [], [], [], [], []
    for k in range(n):
        load = tau[k]
        th_out = th_g - load * k_inv - hb * np.tanh(load / (0.02 * F.TAU_S))
        A_ = np.round((th_out + rng.normal(0.0, 1.0, J) * sd[k]) / F.ENC_LSB) * F.ENC_LSB
        B_ = np.round((th_out + rng.normal(0.0, 1.0, J) * sd[k]) / F.ENC_LSB) * F.ENC_LSB
        meas = 0.5 * (A_ + B_); e = q[k] - meas
        d = (A_[0] - B_[0]) / np.sqrt(2.0)                       # noise draws are identical across configs (th_out differs, d does not)
        n2 += a_w * (d * d - n2)
        r = (q[k] - th_ff)[None, :] + (meas - th_m)
        r_hist = np.roll(r_hist, 1, axis=0); r_hist[0] = r
        if k >= LAGS.max():
            for i_, lg in enumerate(LAGS):
                dd = r - r_hist[lg]; Dk[i_] += a_w * (dd * dd - Dk[i_])
        if k > 1000:
            D = Dk - n2[None, None, :]                           # noise floor removed
            ok_g = D > SNR_GATE * n2[None, None, :]               # amended: lag used only above the floor
            ok_u = D > 0                                          # second amendment: ungated fallback
            use_g = ok_g.sum(0) >= 2
            ok = np.where(use_g[None], ok_g, ok_u)
            nok = ok.sum(0)
            logD = np.where(ok, np.log(np.maximum(D, 1e-30)), 0.0)
            # per (C, J) least-squares slope of log D on log lag over the valid lags
            W = ok.astype(float); Sw = W.sum(0); Sx = (W * LOGL[:, None, None]).sum(0); Sy = (W * logD).sum(0)
            Sxx = (W * LOGL[:, None, None] ** 2).sum(0); Sxy = (W * LOGL[:, None, None] * logD).sum(0)
            den = Sw * Sxx - Sx ** 2
            slope = np.where(den > 1e-12, (Sw * Sxy - Sx * Sy) / np.maximum(den, 1e-12), 1.0)
            Hn = np.clip(0.5 * slope, 0.05, 1.0)
            valid = nok >= 2
            s2 = 0.5 * n2 * DT                                   # pair-mean noise PSD
            if chart == 'canon':
                lam = (W * np.where(ok, D, 0.0) / (LAGS * DT)[:, None, None]).sum(0) / np.maximum(Sw, 1)
                ts = np.sqrt(np.maximum(s2, 1e-30) / np.maximum(lam, 1e-30))
                Hn = np.full_like(Hn, 0.5)
            else:
                loglam = (W * (logD - 2 * Hn[None] * LOGL[:, None, None])).sum(0) / np.maximum(Sw, 1)
                lam = np.exp(loglam)
                ts = kappa_H(Hn) * (np.maximum(s2, 1e-30) / np.maximum(lam, 1e-30)) ** (1.0 / (2 * Hn + 1))
            ts = np.clip(km * ts, TAU_MIN, TAU_MAX)
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
        if store and k % 10 == 0:
            keep.append(th_out[0].copy()); ktau.append(tau_c[0].copy()); kH.append(H[0].copy())
    out = dict(mse=err2 / n)
    if store:
        out.update(th=np.array(keep), t=tt[::10][:len(keep)], tau=np.array(ktau), H=np.array(kH))
    return out

def tune_kappa(t0, chart, tag):
    f = f'{OUT}/tuning_v5c_kappa_{chart}_{tag}.json'
    if os.path.exists(f): return json.load(open(f))
    km = np.repeat(np.array(KM_GRID)[:, None], J, 1); m = np.zeros((len(KM_GRID), J)); ntot = 0
    for s_, (a, b) in enumerate(F.training_spans(t0)):
        if b - a < 5: continue
        tt, q, tau = F.plan(a, b); sd = np.full(len(tt), F.ENC_SD); w = len(tt)
        m += w * simulate(chart, km, tt, q, tau, sd, 100 + s_)['mse']; ntot += w
    ib = np.argmin(m / ntot, 0)
    res = dict(kappa_m={F.JN[j]: float(KM_GRID[ib[j]]) for j in range(J)}, edge={F.JN[j]: bool(ib[j] in (0, len(KM_GRID) - 1)) for j in range(J)},
               rms_deg_by_kappa={str(KM_GRID[i]): float(np.sqrt(np.mean(m[i] / ntot)) / DEG) for i in range(len(KM_GRID))})
    json.dump(res, open(f, 'w'), indent=1); return res

def run(t0, tun, tk_c, tk_u, cond, tag):
    tt, q, tau = F.plan(t0 - 3.0, t0 + 30.0)
    sd = np.full(len(tt), F.ENC_SD)
    if cond == 'step': sd[(tt >= t0 + 10) & (tt < t0 + 20)] *= np.sqrt(10.0)
    seed = 7 if cond == 'step' else 5
    pid = {kk: np.array([[tun['pid'][j][i] for j in F.JN]]) for i, kk in enumerate(('kp', 'ki', 'kd'))}
    b0 = dict(tau=np.array([[tun['b0'][j] for j in F.JN]]))
    one = np.ones((1, J))
    res = {'PID + interpolation': V3.simulate('pid', pid, tt, q, tau, sd, seed, F.PLANT, True),
           'best fixed window': V3.simulate('b0', b0, tt, q, tau, sd, seed, F.PLANT, True),
           'SYN canonical derived': simulate('canon', one, tt, q, tau, sd, seed, store=True),
           'SYN U-prime derived': simulate('uprime', one, tt, q, tau, sd, seed, store=True),
           'SYN canonical oracle kappa': simulate('canon', np.array([[tk_c['kappa_m'][j] for j in F.JN]]), tt, q, tau, sd, seed, store=True),
           'SYN U-prime oracle kappa': simulate('uprime', np.array([[tk_u['kappa_m'][j] for j in F.JN]]), tt, q, tau, sd, seed, store=True)}
    sc = {}
    for name, r in res.items():
        k = r['t'] >= t0 - 1e-9
        s_, e = F.score(r['th'][k], r['t'][k], t0); sc[name] = s_
        if 'tau' in r:
            t_ = r['t'][k]; tau_ = r['tau'][k][:, :6]
            parts = {p: (t_ >= t0 + 10 * i) & (t_ < t0 + 10 * (i + 1)) for i, p in enumerate(('0-10', '10-20', '20-30'))}
            s_['tau_median_ms'] = {p: float(1e3 * np.median(tau_[m_])) for p, m_ in parts.items()}
            s_['tau_step_factor'] = s_['tau_median_ms']['10-20'] / s_['tau_median_ms']['0-10']
            s_['tau_edge_frac'] = float(np.mean((tau_ <= TAU_MIN * 1.001) | (tau_ >= TAU_MAX * 0.999)))
            if 'H' in r:
                H_ = r['H'][k][:, :6]
                s_['H_median_by_joint'] = {F.JN[j]: float(np.median(H_[:, j])) for j in range(6)}
                s_['H_median_arm'] = float(np.median(H_))
            np.savez_compressed(f'{OUT}/v5c_{tag}_{cond}_{name.replace(" ", "_")}.npz', th=r['th'][k], t=t_, err=e, tau=r['tau'][k], H=r.get('H', np.zeros(0)))
    m = {kk: v['tip_rms_mm'] ** 2 for kk, v in sc.items()}
    sc['ratios_rms_over_fixed'] = {kk: float(np.sqrt(m[kk] / m['best fixed window'])) for kk in m}
    return sc

if __name__ == '__main__':
    t_start = time.time(); A, B = 6.0, 217.5
    # baselines on the widened v4 grids (no grid edge this morning), retuned on the regenerated plan, pair reads
    F.PID_GRID = [(kp, ki, kd) for kp in (0, 0.5, 1, 2, 3, 4, 6, 8) for ki in (0, 80, 160, 320, 640, 1280, 2560) for kd in (0, 0.0005, 0.001, 0.002, 0.005)]
    F.B0_GRID = [0.001, 0.0015, 0.002, 0.003, 0.005, 0.007, 0.01, 0.015, 0.02, 0.03, 0.05, 0.07, 0.1, 0.15, 0.2, 0.3, 0.5]
    summary = dict(A=A, B=B, kappa_H_at={h: float(kappa_H(np.array(h))) for h in (0.5, 0.75, 1.0)}, results={}, kappa_tuning={})
    for name, t0 in (('A', A), ('B', B)):
        tun = V3.tune(t0, tag=f'v5_pair_{name}')
        print(name, 'baselines tuned', round(time.time() - t_start), 's', flush=True)
        tk_c = tune_kappa(t0, 'canon', name); tk_u = tune_kappa(t0, 'uprime', name)
        summary['kappa_tuning'][name] = dict(canon=tk_c, uprime=tk_u)
        print(name, 'kappa oracle', {c: {j: v for j, v in t['kappa_m'].items()} for c, t in (('canon', tk_c), ('uprime', tk_u))}, round(time.time() - t_start), 's', flush=True)
        for cond in ('real', 'step'):
            sc = run(t0, tun, tk_c, tk_u, cond, name); summary['results'][f'{name} {cond}'] = sc
            print(name, cond, json.dumps({kk: dict(rms=round(v['tip_rms_mm'], 2), jit=round(v['jitter_mm'], 2),
                                                    tau=v.get('tau_median_ms'), step=round(v['tau_step_factor'], 2) if 'tau_step_factor' in v else None,
                                                    H=round(v['H_median_arm'], 2) if 'H_median_arm' in v else None, edge=round(v['tau_edge_frac'], 3) if 'tau_edge_frac' in v else None)
                                          for kk, v in sc.items() if not kk.startswith('ratios')}, indent=0),
                  'rms/fixed', {kk: round(v, 2) for kk, v in sc['ratios_rms_over_fixed'].items()}, round(time.time() - t_start), 's', flush=True)
            json.dump(summary, open(f'{OUT}/summary_v5c.json', 'w'), indent=1)
    print('done', round(time.time() - t_start), 's')
