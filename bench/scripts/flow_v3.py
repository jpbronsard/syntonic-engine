"""FLOW with the conforming estimator of the financial and deep-learning branches (conforming_est.py from the Fin V3.0
package, SHA256 9a5f3d7f..., used unchanged), at MAESTRO's joints.

Fixed before any run of this file.

Pairing (the analogue of the day pairs). Each 1 ms step the Teensy reads every AS5048A twice, at least one internal
sample apart (the chip samples at 11.25 kHz), so the two reads carry independent noise (assumed) of the same joint
angle (the joint moves under 0.02 deg in 0.1 ms). Every engine, baselines included, then works on the pair mean, so all
see the same information (noise 0.06/sqrt(2) deg); PID and the fixed window are retuned on it.

Engines:
- PID + interpolation and best fixed window: as before, retuned on the pair mean (rest of the movement, clean noise).
- FLOW conforming: per joint, one ConformingTau instance (n = 1, its own stride), fed the pair of servo-model
  residuals r_A = theta_A - theta_model, r_B = theta_B - theta_model (theta_model: nominal 25 ms, unloaded speed
  limit, driven by the actual commands). So: disagreement d = (theta_A - theta_B)/sqrt(2), pure encoder noise; level =
  the model residual, which excludes the engine's own corrections (DIFF-EST's exogenous channel). Identity transform,
  not log: the encoder noise is additive (H16's log map is for multiplicative noise). The instrument's tau is in
  steps; the controller uses tau = tau_inst x 1 ms, bounded to [2 ms, 500 ms] and moved along the rate-linear path at
  most 0.4 s per second, with v1's correction law u += (1 - exp(-dt/tau)) e on the pair-mean error.
- FLOW conforming + meta-kappa: the same, times a dynamic kappa per joint, exactly as in flow_v2 (Phi from parallel
  EMAs on the same residual, noise = the instrument's n2/2, every 0.2 s kappa += clamp(eta (Phi - 0.95), +-0.02),
  eta = 0.2 s / tau, kappa in [0.5, 2]).
Same excerpts (A = 0:06-0:36, B = 3:37-4:07 held out), conditions (real noise; noise variance x10 from 10 to 20 s),
scores and robustness plants as flow_v2."""
import os, sys, json, time, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'ext'))
import flow_vs_pid as F
from conforming_est import ConformingTau

DT = F.DT; OUT = F.OUT

def reads(th_out, sdk, rng, J):
    a = np.round((th_out + rng.normal(0.0, 1.0, J) * sdk) / F.ENC_LSB) * F.ENC_LSB
    b = np.round((th_out + rng.normal(0.0, 1.0, J) * sdk) / F.ENC_LSB) * F.ENC_LSB
    return a, b

def simulate(kind, prm, tt, q, tau, sd, seed, plant=F.PLANT, store=False):
    """kind: 'pid' | 'b0' (vectorized over configs, as flow_vs_pid) | 'conf' | 'conf_meta' (one config)."""
    rng = np.random.default_rng(seed); J = q.shape[1]
    C = next(iter(prm.values())).shape[0] if prm else 1
    Ts = plant['Ts']; k_inv = (plant['stiff_deg'] * F.DEG) / F.TAU_S; hb = 0.5 * plant['play_deg'] * F.DEG
    th_g = np.tile(q[0], (C, 1)); cmd = th_g.copy(); corr = np.zeros((C, J))
    I = np.zeros((C, J)); ef = np.zeros((C, J)); efp = np.zeros((C, J)); u = np.zeros((C, J))
    if kind == 'b0': g_b0 = 1.0 - np.exp(-DT / prm['tau'])
    if kind.startswith('conf'):
        est = [ConformingTau(n=1) for _ in range(J)]
        th_m = q[0].copy(); tau_c = np.full(J, 0.5); kappa = np.ones(J); phi = np.ones(J)
        s = np.zeros(J); er2 = np.zeros(J); ee = np.zeros(J); ee2 = np.zeros(J); a_w = DT / 1.0
    err2 = np.zeros((C, J)); n = len(tt); keep, ktau, kk, kp, kin = [], [], [], [], []
    for k in range(n):
        load = tau[k]
        th_out = th_g - load * k_inv - hb * np.tanh(load / (0.02 * F.TAU_S))
        if kind.startswith('conf'):
            A_, B_ = reads(th_out[0], sd[k], rng, J); meas = 0.5 * (A_ + B_)
            e = q[k] - meas
            tin = np.array([est[j].update(A_[j] - th_m[j], B_[j] - th_m[j])[0] for j in range(J)])
            if k > 1000:
                ts = np.clip(kappa * tin * DT, 0.002, 0.5)
                r = 1.0 / tau_c; dr_max = 0.4 * DT / tau_c ** 2
                r = r + np.clip(1.0 / ts - r, -dr_max, dr_max); tau_c = 1.0 / r
            g = 1.0 - np.exp(-DT / tau_c)
            res = e - s; s = s + g * (e - s)
            er2 += a_w * (res * res - er2); ee += a_w * (e - ee); ee2 += a_w * (e * e - ee2)
            if kind == 'conf_meta' and k > 1000 and k % 200 == 0:
                s2n = np.array([es.n2_hat[0] for es in est]) / 2.0
                var = np.maximum(ee2 - ee * ee, 1e-14)
                phi = np.clip(1.0 - np.maximum(er2 - s2n, 0.0) / var, 0.0, 1.0)
                eta = 0.2 / np.maximum(tau_c, 0.002)
                kappa = np.clip(kappa + np.clip(eta * (phi - 0.95), -0.02, 0.02), 0.5, 2.0)
            u = np.clip(u + g * e, -F.LIM, F.LIM); corr = u
        else:
            A_ = np.round((th_out + rng.normal(0.0, 1.0, J) * sd[k]) / F.ENC_LSB) * F.ENC_LSB
            B_ = np.round((th_out + rng.normal(0.0, 1.0, J) * sd[k]) / F.ENC_LSB) * F.ENC_LSB
            e = q[k] - 0.5 * (A_ + B_)
            if kind == 'pid':
                ef += (DT / 0.005) * (e - ef)
                raw = prm['kp'] * e + prm['ki'] * I + prm['kd'] * (ef - efp) / DT
                sat = np.abs(raw) >= F.LIM
                I += np.where(sat, 0.0, e * DT); efp = ef.copy()
                corr = np.clip(raw, -F.LIM, F.LIM)
            else:
                u = np.clip(u + g_b0 * e, -F.LIM, F.LIM); corr = u
        if k % F.CMD_EVERY == 0: cmd = q[k] + corr
        w_eff = F.W0 * np.maximum(0.05, 1.0 - np.abs(load) / F.TAU_S)
        th_g = th_g + DT * np.clip((cmd - th_g) / Ts, -w_eff, w_eff)
        if kind.startswith('conf'):
            th_m = th_m + DT * np.clip((cmd[0] - th_m) / 0.025, -F.W0, F.W0)
        err2 += (q[k] - th_out) ** 2
        if store and k % 10 == 0:
            keep.append(th_out[0].copy())
            if kind.startswith('conf'):
                ktau.append(tau_c.copy()); kk.append(kappa.copy()); kp.append(phi.copy()); kin.append(tin * DT)
    out = dict(mse=err2 / n)
    if store:
        out.update(th=np.array(keep), t=tt[::10][:len(keep)])
        if kind.startswith('conf'): out.update(tau=np.array(ktau), kappa=np.array(kk), phi=np.array(kp), tau_inst=np.array(kin))
    return out

def tune(t0, plant=F.PLANT, tag='pair'):
    f = f'{OUT}/tuning_{tag}.json'
    if os.path.exists(f): return json.load(open(f))
    G = np.array(F.PID_GRID); C = len(G)
    prm = dict(kp=np.repeat(G[:, :1], 12, 1), ki=np.repeat(G[:, 1:2], 12, 1), kd=np.repeat(G[:, 2:], 12, 1))
    tb = np.array(F.B0_GRID); prm_b = dict(tau=np.repeat(tb[:, None], 12, 1))
    m_pid = np.zeros((C, 12)); m_b0 = np.zeros((len(tb), 12)); ntot = 0
    for s_, (a, b) in enumerate(F.training_spans(t0)):
        if b - a < 5: continue
        tt, q, tau = F.plan(a, b); sd = np.full(len(tt), F.ENC_SD); w = len(tt)
        m_pid += w * simulate('pid', prm, tt, q, tau, sd, 100 + s_, plant)['mse']
        m_b0 += w * simulate('b0', prm_b, tt, q, tau, sd, 100 + s_, plant)['mse']; ntot += w
    ip = np.argmin(m_pid / ntot, 0); ib = np.argmin(m_b0 / ntot, 0)
    res = dict(pid={F.JN[j]: list(map(float, G[ip[j]])) for j in range(12)}, b0={F.JN[j]: float(tb[ib[j]]) for j in range(12)})
    json.dump(res, open(f, 'w'), indent=1); return res

def run(t0, tun, cond, plant=F.PLANT, tag=''):
    tt, q, tau = F.plan(t0 - 3.0, t0 + 30.0)
    sd = np.full(len(tt), F.ENC_SD)
    if cond == 'step': sd[(tt >= t0 + 10) & (tt < t0 + 20)] *= np.sqrt(10.0)
    seed = 7 if cond == 'step' else 5
    pid = {kk: np.array([[tun['pid'][j][i] for j in F.JN]]) for i, kk in enumerate(('kp', 'ki', 'kd'))}
    b0 = dict(tau=np.array([[tun['b0'][j] for j in F.JN]]))
    res = {'PID + interpolation': simulate('pid', pid, tt, q, tau, sd, seed, plant, True),
           'best fixed window': simulate('b0', b0, tt, q, tau, sd, seed, plant, True),
           'FLOW conforming': simulate('conf', {}, tt, q, tau, sd, seed, plant, True),
           'FLOW conforming + meta-kappa': simulate('conf_meta', {}, tt, q, tau, sd, seed, plant, True)}
    sc = {}
    for name, r in res.items():
        k = r['t'] >= t0 - 1e-9
        s_, e = F.score(r['th'][k], r['t'][k], t0); sc[name] = s_
        if tag:
            extra = {x: r[x][k] for x in ('tau', 'kappa', 'phi', 'tau_inst') if x in r}
            np.savez_compressed(f'{OUT}/v3_{tag}_{cond}_{name.replace(" ", "_").replace("+", "plus")}.npz', th=r['th'][k], t=r['t'][k], err=e, **extra)
    m = {kk: v['tip_rms_mm'] ** 2 for kk, v in sc.items()}
    sc['ratios'] = {f'{kk} / PID': m[kk] / m['PID + interpolation'] for kk in m} | {f'{kk} / fixed': m[kk] / m['best fixed window'] for kk in m}
    return sc

if __name__ == '__main__':
    t_start = time.time(); A, B = 6.0, 217.5
    summary = dict(A=A, B=B, nominal={}, robustness={})
    for name, t0, tag in (('A', A, 'pair'), ('B', B, 'pair_heldout')):
        tun = tune(t0, tag=tag)
        for cond in ('real', 'step'):
            sc = run(t0, tun, cond, tag=name); summary['nominal'][f'{name} {cond}'] = sc
            print(name, cond, json.dumps({kk: (round(v['tip_rms_mm'], 2), round(v['jitter_mm'], 2), round(v['tip_10-20 s'], 2), round(v['jitter_10-20 s'], 2))
                                          for kk, v in sc.items() if kk != 'ratios'}), flush=True)
    for ptag, plant in (('Ts15', dict(F.PLANT, Ts=0.015)), ('Ts40', dict(F.PLANT, Ts=0.040)),
                        ('play0', dict(F.PLANT, play_deg=0.0)), ('play0.6', dict(F.PLANT, play_deg=0.6))):
        tun = tune(B, plant, f'pair_heldout_{ptag}')
        for cond in ('real', 'step'):
            sc = run(B, tun, cond, plant); summary['robustness'][f'B {ptag} {cond}'] = sc
            print('B', ptag, cond, json.dumps({kk: round(v['tip_rms_mm'], 2) for kk, v in sc.items() if kk != 'ratios'}), flush=True)
    json.dump(summary, open(f'{OUT}/summary_v3.json', 'w'), indent=1)
    print('done', round(time.time() - t_start), 's')
