"""Log the syntonic engine's window tau* against the gesture's acceleration on the MAESTRO 1 servo bench
(examples/servo_bench.py arithmetic, published engine 'raw', steady noise), plus the same engine with its window
frozen at its own mean. Excerpt and seed fixed before looking: A, steady noise, seed 5 (the published ones).
Usage: python accel_tau.py <repo_root> <excerpt> <out.npz>"""
import os, sys
import numpy as np
from scipy.interpolate import CubicSpline

ROOT = sys.argv[1]; NAME = sys.argv[2]; OUT = sys.argv[3]
sys.path.insert(0, os.path.join(ROOT, 'python')); sys.path.insert(0, os.path.join(ROOT, 'examples'))
from syntonic_engine import SyntonicEngine
from arm_fk import tip

DT = float(os.environ.get('DT', 1e-3)); FRAME = float(os.environ.get('FRAME', 3e-3)); EVERY = int(round(0.01 / DT))
DEG = np.pi / 180
Z = np.load(os.path.join(ROOT, 'examples', 'data', 'maestro_excerpts.npz'))
TAU_S = Z['stall_torque']; W0 = Z['no_load_speed']
ENC_SD, ENC_LSB = float(Z['enc_sd_rad']), float(Z['enc_lsb_rad'])
TS, K_INV, HB = 0.025, (1.0 * DEG) / TAU_S, 0.5 * 0.3 * DEG
SEED = int(sys.argv[4]) if len(sys.argv) > 4 else 5


def run(tau_fix=None):
    t0 = float(Z[f'{NAME}_t0']); tt = np.arange(t0 - 3.0, t0 + 30.0, DT)
    S = CubicSpline(Z[f'{NAME}_t30'], Z[f'{NAME}_q30'], axis=0)
    q = S(tt); tau = CubicSpline(Z[f'{NAME}_t30'], Z[f'{NAME}_tau30'], axis=0)(tt)
    rng = np.random.default_rng(SEED); J = q.shape[1]
    th_g = q[0].copy()
    eng = SyntonicEngine(q[0], DT, FRAME, W0, identification='raw')
    if tau_fix is not None:
        tf = np.asarray(tau_fix, float).copy()
        eng.floor = tf.copy(); eng.tau_max = tf.copy(); eng.tau = tf.copy(); eng.tau_target = tf.copy()
    T, WIN, TH, QP = [], [], [], []
    for k in range(len(tt)):
        load = tau[k]
        th_out = th_g - load * K_INV - HB * np.tanh(load / (0.02 * TAU_S))
        a = np.round((th_out + rng.normal(0.0, 1.0, J) * ENC_SD) / ENC_LSB) * ENC_LSB
        b = np.round((th_out + rng.normal(0.0, 1.0, J) * ENC_SD) / ENC_LSB) * ENC_LSB
        cmd = eng.step(q[k], a, b)
        w_eff = W0 * np.maximum(0.05, 1.0 - np.abs(load) / TAU_S)
        th_g = th_g + DT * np.clip((cmd - th_g) / TS, -w_eff, w_eff)
        if k % EVERY == 0 and tt[k] >= t0 - 1e-9:
            T.append(tt[k] - t0); WIN.append(np.broadcast_to(eng.tau, (J,)).copy())
            TH.append(th_out[:6].copy()); QP.append(q[k, :6].copy())
    T = np.array(T); WIN = np.array(WIN); TH = np.array(TH); QP = np.array(QP)
    qdd = S(T + t0, 2)[:, :6]
    tip_act = np.array([tip(x) for x in TH]); tip_plan = np.array([tip(x) for x in QP])
    return T, WIN, qdd, tip_act, tip_plan, float(eng.floor if np.ndim(eng.floor) == 0 else np.mean(eng.floor))


T, WIN, qdd, tip_a, tip_p, floor = run()
fix = WIN.mean(0)
T2, WIN2, _, tip_a2, _, _ = run(tau_fix=fix)
np.savez(OUT, T=T, win=WIN[:, :6], win_frozen=WIN2[:, :6], qdd=qdd, tip_act=tip_a, tip_act_frozen=tip_a2, tip_plan=tip_p,
         floor=floor, tau_fix=fix)
e1 = 1e3 * np.sqrt(np.mean(np.sum((tip_a - tip_p) ** 2, 1))); e2 = 1e3 * np.sqrt(np.mean(np.sum((tip_a2 - tip_p) ** 2, 1)))
print(NAME, 'floor ms', round(floor * 1e3, 3), 'tau mean ms by joint', np.round(fix[:6] * 1e3, 2))
print(NAME, 'tip rms mm adaptive', round(e1, 3), 'frozen', round(e2, 3))
