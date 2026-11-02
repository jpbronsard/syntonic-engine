"""One cell of the rate grid: servo command rate x loop rate, one excerpt, one seed. Runs the published engine
(adaptive) then the same engine frozen at its own mean window; prints one JSON line.
Usage: python mc_cell.py <repo_root> <servo_hz> <loop_hz> <excerpt> <seed>"""
import os, sys, json
import numpy as np
from scipy.interpolate import CubicSpline
ROOT = sys.argv[1]; SERVO_HZ, LOOP_HZ = float(sys.argv[2]), float(sys.argv[3]); NAME = sys.argv[4]; SEED = int(sys.argv[5])
sys.path.insert(0, os.path.join(ROOT, 'python')); sys.path.insert(0, os.path.join(ROOT, 'examples'))
from syntonic_engine import SyntonicEngine
from arm_fk import tip
DT, FRAME = 1.0 / LOOP_HZ, 1.0 / SERVO_HZ; EVERY = int(round(0.01 / DT)); DEG = np.pi / 180
Z = np.load(os.path.join(ROOT, 'examples', 'data', 'maestro_excerpts.npz'))
TAU_S = Z['stall_torque']; W0 = Z['no_load_speed']
ENC_SD, ENC_LSB = float(Z['enc_sd_rad']), float(Z['enc_lsb_rad'])
TS, K_INV, HB = 0.025, (1.0 * DEG) / TAU_S, 0.5 * 0.3 * DEG

def run(tau_fix=None):
    t0 = float(Z[f'{NAME}_t0']); tt = np.arange(t0 - 3.0, t0 + 30.0, DT)
    q = CubicSpline(Z[f'{NAME}_t30'], Z[f'{NAME}_q30'], axis=0)(tt)
    tau = CubicSpline(Z[f'{NAME}_t30'], Z[f'{NAME}_tau30'], axis=0)(tt)
    rng = np.random.default_rng(SEED); J = q.shape[1]; th_g = q[0].copy()
    eng = SyntonicEngine(q[0], DT, FRAME, W0, identification='raw')
    if tau_fix is not None:
        tf = np.asarray(tau_fix, float).copy(); eng.floor = tf.copy(); eng.tau_max = tf.copy(); eng.tau = tf.copy(); eng.tau_target = tf.copy()
    keep = []; wsum = np.zeros(J); nw = 0; atf = np.zeros(J); fprev = -1
    for k in range(len(tt)):
        load = tau[k]
        th_out = th_g - load * K_INV - HB * np.tanh(load / (0.02 * TAU_S))
        a = np.round((th_out + rng.normal(0.0, 1.0, J) * ENC_SD) / ENC_LSB) * ENC_LSB
        b = np.round((th_out + rng.normal(0.0, 1.0, J) * ENC_SD) / ENC_LSB) * ENC_LSB
        fi = int(k * DT / FRAME + 1e-9); fr = (k == 0) or (fi != fprev); fprev = fi
        cmd = eng.step(q[k], a, b, frame=fr)
        w_eff = W0 * np.maximum(0.05, 1.0 - np.abs(load) / TAU_S)
        th_g = th_g + DT * np.clip((cmd - th_g) / TS, -w_eff, w_eff)
        if tt[k] >= t0:
            tv = np.broadcast_to(eng.tau, (J,)); wsum += tv; nw += 1
            if tau_fix is None: atf += tv < float(eng.floor) * 1.02
        if k % EVERY == 0 and tt[k] >= t0 - 1e-9: keep.append((th_out[:6].copy(), q[k, :6].copy()))
    e = np.array([tip(x) - tip(y) for x, y in keep])
    return 1e3 * np.sqrt(np.mean(np.sum(e ** 2, 1))), wsum / nw, atf / max(nw, 1), float(np.min(np.atleast_1d(eng.floor)))

ta, wmean, atf, floor = run()
tf_, _, _, _ = run(tau_fix=wmean)
print(json.dumps(dict(servo=SERVO_HZ, loop=LOOP_HZ, excerpt=NAME, seed=SEED, floor_ms=floor * 1e3, tip_adaptive=ta,
                      tip_frozen=tf_, gain=1 - ta / tf_, at_floor=atf[:6].tolist())), flush=True)
