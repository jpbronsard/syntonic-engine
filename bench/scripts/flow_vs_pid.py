"""FLOW against PID with interpolation, at the joints, on 30 s of the Karajan-scale gestures.

Everything below was fixed before any engine was run.

Plan and loads. The plan is the delivered Karajan run (video/big_solution.npz, BLS3355 at J3, BLS3375 at J4,
7.4 V): joint angles at 30 Hz, the same interpolation for every engine (cubic spline to 1 kHz). Load torque on each
joint is the inverse dynamics of the planned motion (the servo check's own computation). The turntable J0 is ideal.

Servo and encoder (one model per joint, both arms; values assumed where no data sheet gives them):
- command refreshed at 330 Hz (DSSERVO maximum PWM rate), held between refreshes;
- internal position loop: first-order, time constant 25 ms, speed limited to the no-load speed derated by load
  (w0 (1 - |tau|/tau_s), the line used by the servo check, floor 5 %);
- compliance: 1.0 deg of deflection at stall torque, proportional to load;
- gear play 0.3 deg total, taken up on the loaded side (tanh over 2 % of stall);
- encoder AS5048A on the joint output: 0.06 deg rms white noise and 0.0219 deg steps (data sheet), read at 1 kHz,
  assumed calibrated (its +-0.8 deg non-linearity is not modelled).
Robustness (reported, never used to choose anything): servo time constant 15 and 40 ms; gear play 0 and 0.6 deg.

Engines (all add a correction to the interpolated plan; correction limited to +-5 deg, no integration while limited):
- Open loop: the plan alone (reference row).
- PID + interpolation: P, I and D on the encoder error (D on the error low-passed at 5 ms). Gains per joint from a
  grid (Kp 0 .. 4, Ki 0 .. 80 /s, Kd 0 .. 0.01 s), chosen once on the rest of the movement with clean noise.
- Best fixed window (B0): first-order flow of the correction, u += (1 - exp(-dt/tau)) e, tau per joint from a grid
  of 2 ms to 500 ms, chosen the same way.
- FLOW: the same first-order flow with tau = sqrt(sigma2 / lambda) per joint (kappa = 1, MAESTRO-FLOW V1.0). sigma2 and
  lambda come from the structure function of the measured error, V(q) = q lambda dt + 2 sigma2 fitted over
  q = 2 .. 10 samples (V1.0: fitted from q = 2), averaged over the last 1 s. tau is bounded to [2 ms, 500 ms] (H4) and
  moves along the rate-linear path (the e-geodesic of Dual Geodesy) at most 0.4 s per second (H4: below 1/(4 tau_max)).
  Nothing in FLOW is tuned.

Excerpt. The continuous 30 s, starting on a 0.5 s grid, with at least three soft (< -33 dB) and three loud (> -22 dB)
bars, that maximizes the spread (p90 / p10) of the baton arm's per-bar joint speed. Chosen from the plan only.

Every engine runs from 3 s before the excerpt (run-in, not scored), so FLOW's estimate and the integrators settle.

Conditions on the excerpt. Real noise throughout; and a synthetic step (Phase -1 style): encoder noise variance x10
(rms x3.16) from 10 s to 20 s into the excerpt, on every joint. Same noise draws for every engine.

Scores. Baton-tip position error (actual against planned, from forward kinematics at 100 Hz): rms over the excerpt
(primary), by part (0-10, 10-20, 20-30 s) and in soft and loud bars; jitter = rms of the tip error above 5 Hz. Ratios
of mean square errors: FLOW / PID and FLOW / best fixed window (below 1 means FLOW does better)."""
import os, sys, json, time, numpy as np, mujoco
from scipy.interpolate import CubicSpline
from scipy.ndimage import gaussian_filter1d
from scipy.signal import butter, sosfiltfilt
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); os.chdir(ROOT); sys.path.insert(0, ROOT)
import presto as P, sim
P.mm.SERVO[4] = ('DSSERVO BLS3375 7.4 V', 75 * 0.0980665, np.radians(60) / 0.12)
assert P.mm.SERVO[3][0].startswith('DSSERVO BLS3355')

DT = 1e-3; CMD_EVERY = 3
DEG = np.pi / 180
PLANT = dict(Ts=0.025, stiff_deg=1.0, play_deg=0.3)
ENC_SD = 0.06 * DEG; ENC_LSB = 0.0219 * DEG
LIM = 5 * DEG
OUT = 'flowcmp/out'; os.makedirs(OUT, exist_ok=True)

# ------------------------------------------------------------------ plan, loads, servo ratings
MY = mujoco.MjModel.from_xml_string(P.with_yaw(P.mm.make_xml(P.P_DES)))
JN = [f'{s}J{k}' for s in 'RL' for k in range(1, 7)]
QADR = np.array([MY.jnt_qposadr[MY.joint(n).id] for n in JN]); VADR = np.array([MY.jnt_dofadr[MY.joint(n).id] for n in JN])
Z = np.load('video/big_solution.npz'); T30 = Z['t']; QY = Z['Q']
QS = gaussian_filter1d(QY, 1.0, axis=0, mode='nearest')
_arm = sim.Arm(MY, 'R')
TAU_S = np.tile(_arm.tau_s, 2); W0 = np.tile(_arm.w0, 2)

def loads():
    f = f'{OUT}/loads.npy'
    if os.path.exists(f): return np.load(f)
    d = mujoco.MjData(MY); qd = np.gradient(QS, P.DT, axis=0); qdd = np.gradient(qd, P.DT, axis=0)
    MY.opt.disableflags |= mujoco.mjtDisableBit.mjDSBL_CONTACT
    tau = np.zeros((len(QS), len(JN)))
    for i in range(len(QS)):
        d.qpos[:] = QS[i]; d.qvel[:] = qd[i]; d.qacc[:] = qdd[i]; mujoco.mj_inverse(MY, d); tau[i] = d.qfrc_inverse[VADR]
    np.save(f, tau); return tau
TAU30 = loads()
SPL_Q = CubicSpline(T30, QS[:, QADR], axis=0); SPL_T = CubicSpline(T30, TAU30, axis=0)

def plan(t0, t1):
    tt = np.arange(t0, t1, DT); return tt, SPL_Q(tt), SPL_T(tt)

# ------------------------------------------------------------------ excerpt
D = json.load(open('fifth/bars_full.json')); LT, LDB = np.array(D['loud']['t']), np.array(D['loud']['db'])
BEATS = [b['t'] for b in D['bars'] if b['kind'] == 'beat']
def bar_table():
    v = np.gradient(QS[:, QADR[:6]], P.DT, axis=0)
    rows = []
    for a, b in zip(BEATS[:-1], BEATS[1:]):
        if b - a > 1.0: continue
        k = (T30 >= a) & (T30 < b)
        rows.append((a, b, float(LDB[(LT >= a) & (LT < b)].max()), float(np.sqrt(np.mean(np.sum(v[k] ** 2, 1))))))
    return np.array(rows)
BT = bar_table()
def pick_excerpt():
    best = None
    for s in np.arange(5.0, T30[-1] - 31.0, 0.5):
        k = (BT[:, 0] >= s) & (BT[:, 1] <= s + 30)
        if (BT[k, 2] < -33).sum() < 3 or (BT[k, 2] > -22).sum() < 3: continue
        sp = np.percentile(BT[k, 3], 90) / np.percentile(BT[k, 3], 10)
        if best is None or sp > best[1]: best = (float(s), float(sp))
    return best

# ------------------------------------------------------------------ simulation
def simulate(kind, prm, tt, q, tau, sd, seed, plant=PLANT, store=False):
    """Vectorized over configurations C (leading axis of the parameter arrays, shape (C, 12)).
    Returns the mean square joint error (C, 12) against the plan, and stored trajectories if asked."""
    rng = np.random.default_rng(seed)
    C = next(iter(prm.values())).shape[0] if prm else 1
    J = q.shape[1]
    Ts = plant['Ts']; k_inv = (plant['stiff_deg'] * DEG) / TAU_S; hb = 0.5 * plant['play_deg'] * DEG
    th_g = np.tile(q[0], (C, 1)); cmd = th_g.copy(); corr = np.zeros((C, J))
    I = np.zeros((C, J)); ef = np.zeros((C, J)); efp = np.zeros((C, J)); u = np.zeros((C, J))
    if kind == 'flow':
        QS_ = np.arange(2, 11); xm = QS_.mean(); xd = QS_ - xm; sxx = float(np.sum(xd ** 2))
        hist = np.zeros((11, C, J)); V = np.zeros((len(QS_), C, J)); a_w = DT / 1.0
        tau_c = np.full((C, J), 0.5); nfill = 0
    if kind == 'b0': g_b0 = 1.0 - np.exp(-DT / prm['tau'])
    err2 = np.zeros((C, J)); n = len(tt)
    keep = []; ktau = []
    for k in range(n):
        load = tau[k]
        th_out = th_g - load * k_inv - hb * np.tanh(load / (0.02 * TAU_S))
        meas = th_out + rng.normal(0.0, 1.0, J) * sd[k]
        meas = np.round(meas / ENC_LSB) * ENC_LSB
        e = q[k] - meas
        if kind == 'pid':
            ef += (DT / 0.005) * (e - ef)
            raw = prm['kp'] * e + prm['ki'] * I + prm['kd'] * (ef - efp) / DT
            sat = np.abs(raw) >= LIM
            I += np.where(sat, 0.0, e * DT); efp = ef.copy()
            corr = np.clip(raw, -LIM, LIM)
        elif kind == 'b0':
            u = np.clip(u + g_b0 * e, -LIM, LIM); corr = u
        elif kind == 'flow':
            hist = np.roll(hist, 1, axis=0); hist[0] = e; nfill = min(nfill + 1, 11)
            if nfill > 10:
                for i_, qq in enumerate(QS_):
                    V[i_] += a_w * ((e - hist[qq]) ** 2 - V[i_])
                if k > 1000:                               # one window of data before the estimate is used
                    Vm = V.mean(0); slope = np.tensordot(xd, V, 1) / sxx; icpt = Vm - slope * xm
                    lam = slope / DT; s2 = 0.5 * icpt * DT
                    ts = np.where(slope <= 0, 0.5, np.where(icpt <= 0, 0.002,
                                  np.sqrt(np.maximum(s2, 1e-30) / np.maximum(lam, 1e-30))))
                    ts = np.clip(ts, 0.002, 0.5)
                    r = 1.0 / tau_c; dr_max = 0.4 * DT / tau_c ** 2       # rate-linear move, |dtau/dt| <= 0.4
                    r = r + np.clip(1.0 / ts - r, -dr_max, dr_max); tau_c = 1.0 / r
            u = np.clip(u + (1.0 - np.exp(-DT / tau_c)) * e, -LIM, LIM); corr = u
        if k % CMD_EVERY == 0: cmd = q[k] + corr
        w_eff = W0 * np.maximum(0.05, 1.0 - np.abs(load) / TAU_S)
        th_g = th_g + DT * np.clip((cmd - th_g) / Ts, -w_eff, w_eff)
        err2 += (q[k] - th_out) ** 2
        if store and k % 10 == 0:
            keep.append(th_out.copy())
            if kind == 'flow': ktau.append(tau_c.copy())
    out = dict(mse=err2 / n)
    if store:
        out['th'] = np.array(keep); out['t'] = tt[::10][:len(keep)]
        if kind == 'flow': out['tau'] = np.array(ktau)
    return out

# ------------------------------------------------------------------ tuning on the rest of the movement
PID_GRID = [(kp, ki, kd) for kp in (0, 0.5, 1, 2, 4) for ki in (0, 5, 10, 20, 40, 80) for kd in (0, 0.002, 0.005, 0.01)]
B0_GRID = [0.002, 0.003, 0.005, 0.007, 0.01, 0.015, 0.02, 0.03, 0.05, 0.07, 0.1, 0.15, 0.2, 0.3, 0.5]

def training_spans(t0):
    return [(1.0, t0 - 2.0), (t0 + 32.0, T30[-1] - 1.0)]

def tune(t0, plant=PLANT, tag='nominal'):
    f = f'{OUT}/tuning_{tag}.json'
    if os.path.exists(f): return json.load(open(f))
    G = np.array(PID_GRID); C = len(G)
    prm = dict(kp=np.repeat(G[:, :1], 12, 1), ki=np.repeat(G[:, 1:2], 12, 1), kd=np.repeat(G[:, 2:], 12, 1))
    tb = np.array(B0_GRID); prm_b = dict(tau=np.repeat(tb[:, None], 12, 1))
    m_pid = np.zeros((C, 12)); m_b0 = np.zeros((len(tb), 12)); ntot = 0
    for s, (a, b) in enumerate(training_spans(t0)):
        if b - a < 5: continue
        tt, q, tau = plan(a, b); sd = np.full(len(tt), ENC_SD); w = len(tt)
        m_pid += w * simulate('pid', prm, tt, q, tau, sd, 100 + s, plant)['mse']
        m_b0 += w * simulate('b0', prm_b, tt, q, tau, sd, 100 + s, plant)['mse']; ntot += w
    m_pid /= ntot; m_b0 /= ntot
    ip = np.argmin(m_pid, 0); ib = np.argmin(m_b0, 0)
    res = dict(pid={JN[j]: list(map(float, G[ip[j]])) for j in range(12)}, b0={JN[j]: float(tb[ib[j]]) for j in range(12)},
               pid_rms_deg={JN[j]: float(np.sqrt(m_pid[ip[j], j]) / DEG) for j in range(12)},
               b0_rms_deg={JN[j]: float(np.sqrt(m_b0[ib[j], j]) / DEG) for j in range(12)})
    json.dump(res, open(f, 'w'), indent=1); return res

# ------------------------------------------------------------------ scoring
def tips(th, t):
    d = mujoco.MjData(MY); out = np.zeros((len(t), 3)); ref = np.zeros((len(t), 3))
    qp = CubicSpline(T30, QS, axis=0)(t)
    for i in range(len(t)):
        d.qpos[:] = qp[i]; mujoco.mj_kinematics(MY, d); ref[i] = d.site('Rtip').xpos
        d.qpos[QADR] = th[i]; mujoco.mj_kinematics(MY, d); out[i] = d.site('Rtip').xpos
    return out, ref

def score(th, t, t0):
    a, r = tips(th, t); e = a - r
    hp = sosfiltfilt(butter(4, 5.0, 'highpass', fs=100.0, output='sos'), e, axis=0)
    def rms(x, k): return float(1e3 * np.sqrt(np.mean(np.sum(x[k] ** 2, 1))))
    allk = np.ones(len(t), bool)
    parts = {f'{p}': (t >= t0 + 10 * i) & (t < t0 + 10 * (i + 1)) for i, p in enumerate(('0-10 s', '10-20 s', '20-30 s'))}
    soft = np.zeros(len(t), bool); loud = np.zeros(len(t), bool)
    for a_, b_, db, _ in BT:
        k = (t >= a_) & (t < b_)
        if db < -33: soft |= k
        if db > -22: loud |= k
    return dict(tip_rms_mm=rms(e, allk), jitter_mm=rms(hp, allk), **{f'tip_{p}': rms(e, k) for p, k in parts.items()},
                **{f'jitter_{p}': rms(hp, k) for p, k in parts.items()},
                tip_soft=rms(e, soft) if soft.any() else None, tip_loud=rms(e, loud) if loud.any() else None), e

def run_excerpt(t0, tun, cond, plant=PLANT, store=True):
    tt, q, tau = plan(t0 - 3.0, t0 + 30.0)
    sd = np.full(len(tt), ENC_SD)
    if cond == 'step': sd[(tt >= t0 + 10) & (tt < t0 + 20)] *= np.sqrt(10.0)
    seed = 7 if cond == 'step' else 5
    pid = {k: np.array([[tun['pid'][j][i] for j in JN]]) for i, k in enumerate(('kp', 'ki', 'kd'))}
    b0 = dict(tau=np.array([[tun['b0'][j] for j in JN]]))
    res = {}
    for name, kind, prm in (('open loop', 'pid', dict(kp=np.zeros((1, 12)), ki=np.zeros((1, 12)), kd=np.zeros((1, 12)))),
                            ('PID + interpolation', 'pid', pid), ('best fixed window', 'b0', b0), ('FLOW', 'flow', {})):
        r = simulate(kind, prm, tt, q, tau, sd, seed, plant, store=store)
        res[name] = r
    return tt, res

if __name__ == '__main__':
    t_start = time.time()
    ex = pick_excerpt(); t0 = ex[0]
    print('excerpt', ex, flush=True)
    tun = tune(t0); print('tuned', round(time.time() - t_start), 's', json.dumps(tun), flush=True)
    summary = dict(excerpt=dict(start=t0, end=t0 + 30, spread=ex[1]), tuning=tun, conditions={})
    for cond in ('real', 'step'):
        tt, res = run_excerpt(t0, tun, cond)
        sc = {}
        for name, r in res.items():
            k = r['t'] >= t0 - 1e-9
            r['th'], r['t'] = r['th'][k], r['t'][k]
            if 'tau' in r: r['tau'] = r['tau'][k]
            s, e = score(r['th'][:, :6], r['t'], t0); sc[name] = s
            np.savez_compressed(f'{OUT}/{cond}_{name.replace(" ", "_").replace("+", "plus")}.npz', th=r['th'], t=r['t'],
                                tau=r.get('tau', np.zeros(0)), err=e)
        mse = {k: v['tip_rms_mm'] ** 2 for k, v in sc.items()}
        sc['ratios'] = dict(FLOW_over_PID=mse['FLOW'] / mse['PID + interpolation'], FLOW_over_fixed=mse['FLOW'] / mse['best fixed window'])
        summary['conditions'][cond] = sc
        print(cond, json.dumps(sc, indent=1), flush=True)
    json.dump(summary, open(f'{OUT}/summary.json', 'w'), indent=1)
    print('done', round(time.time() - t_start), 's')
