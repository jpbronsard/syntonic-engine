"""The bench's tuned PID (gains from maestro_excerpts.npz, authority clip), with and without the same servo-lag lead.
Usage: python pid_ff.py <repo_root> <excerpt> <steady|step> <seed>"""
import os, sys
import numpy as np
from scipy.interpolate import CubicSpline
ROOT, NAME, COND, SEED = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
sys.path.insert(0, os.path.join(ROOT, 'examples'))
from arm_fk import tip
DT, FRAME = 1e-3, 3e-3; DEG = np.pi / 180
Z = np.load(os.path.join(ROOT, 'examples', 'data', 'maestro_excerpts.npz'))
TAU_S = Z['stall_torque']; W0 = Z['no_load_speed']
ENC_SD, ENC_LSB = float(Z['enc_sd_rad']), float(Z['enc_lsb_rad'])
TS, K_INV, HB = 0.025, (1.0 * DEG) / TAU_S, 0.5 * 0.3 * DEG
LIM = TS * W0

def run(ff):
    t0 = float(Z[f'{NAME}_t0']); tt = np.arange(t0 - 3.0, t0 + 30.0, DT)
    S = CubicSpline(Z[f'{NAME}_t30'], Z[f'{NAME}_q30'], axis=0)
    q = S(tt); qd = S(tt, 1); tau = CubicSpline(Z[f'{NAME}_t30'], Z[f'{NAME}_tau30'], axis=0)(tt)
    sd = np.full(len(tt), ENC_SD)
    if COND == 'step': sd[(tt >= t0 + 10) & (tt < t0 + 20)] *= np.sqrt(10.0)
    rng = np.random.default_rng(SEED); J = q.shape[1]
    th_g = q[0].copy(); cmd = q[0].copy()
    kp, ki, kd = (Z[f'{NAME}_pid'][:, i] for i in range(3))
    I = np.zeros(J); ef = np.zeros(J); efp = np.zeros(J)
    keep = []; err2 = np.zeros(J); m = 0
    for k in range(len(tt)):
        load = tau[k]
        th_out = th_g - load * K_INV - HB * np.tanh(load / (0.02 * TAU_S))
        a = np.round((th_out + rng.normal(0.0, 1.0, J) * sd[k]) / ENC_LSB) * ENC_LSB
        b = np.round((th_out + rng.normal(0.0, 1.0, J) * sd[k]) / ENC_LSB) * ENC_LSB
        e = q[k] - 0.5 * (a + b)
        ef += (DT / 0.005) * (e - ef)
        raw = kp * e + ki * I + kd * (ef - efp) / DT
        I += np.where(np.abs(raw) >= LIM, 0.0, e * DT); efp = ef.copy()
        if k % 3 == 0: cmd = q[k] + np.clip(raw, -LIM, LIM) + (TS * qd[k] if ff else 0.0)
        w_eff = W0 * np.maximum(0.05, 1.0 - np.abs(load) / TAU_S)
        th_g = th_g + DT * np.clip((cmd - th_g) / TS, -w_eff, w_eff)
        if tt[k] >= t0: err2 += (q[k] - th_out) ** 2; m += 1
        if k % 10 == 0 and tt[k] >= t0 - 1e-9: keep.append((th_out[:6].copy(), q[k, :6].copy()))
    e = np.array([tip(x) - tip(y) for x, y in keep])
    return 1e3 * np.sqrt(np.mean(np.sum(e ** 2, 1)))

p0 = run(False); p1 = run(True)
print(f'{NAME} {COND:6s} seed {SEED} | PID {p0:.3f} mm, PID+lead {p1:.3f} mm ({1 - p1 / p0:+.0%})')
