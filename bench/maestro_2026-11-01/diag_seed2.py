"""Exploratory, after the verdicts: why the lead arms lose in A steady seed 2. Per-joint rms error and time with the
correction at its authority limit, F15 against F15L, seeds 2 and 3."""
import sys, numpy as np
import os
R = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')); sys.argv = ['x', R, 'none']
import eg732_lead as M
from scipy.interpolate import CubicSpline
def go(seed):
    name, cond, dt, joints, arms = 'A', 'steady', 5e-4, list(range(6)), ['F15', 'F15L']
    t0 = float(M.Z[f'{name}_t0']); tt = np.arange(t0 - 3.0, t0 + 30.0, dt); J = 6; A = 2
    S = CubicSpline(M.Z[f'{name}_t30'], M.Z[f'{name}_q30'][:, joints], axis=0); qJ = S(tt); qdJ = S(tt, 1)
    tauJ = CubicSpline(M.Z[f'{name}_t30'], M.Z[f'{name}_tau30'][:, joints], axis=0)(tt)
    rng = np.random.default_rng(seed); tauS = np.tile(M.TAU_S[joints], A); w0 = np.tile(M.W0[joints], A)
    k_inv = (1.0 * M.DEG) / tauS; hb = 0.5 * 0.3 * M.DEG
    allowed = np.zeros((A * J, len(M.FR)), bool); allowed[:, M.FR.index(1.5)] = True
    eng = M.EGLEngine(np.tile(qJ[0], A), dt, w0, allowed, np.zeros(A * J, int), np.r_[np.zeros(J, bool), np.ones(J, bool)])
    th_g = np.tile(qJ[0], A); err2 = np.zeros(A * J); atl = np.zeros(A * J); n = 0; worst = np.zeros(A * J); tw = np.zeros(A * J)
    for k in range(len(tt)):
        load = np.tile(tauJ[k], A); q = np.tile(qJ[k], A); qd = np.tile(qdJ[k], A)
        th_out = th_g - load * k_inv - hb * np.tanh(load / (0.02 * tauS))
        na = np.tile(rng.normal(0, 1, J), A); nb = np.tile(rng.normal(0, 1, J), A)
        ra = np.round((th_out + na * M.ENC_SD) / M.ENC_LSB) * M.ENC_LSB; rb = np.round((th_out + nb * M.ENC_SD) / M.ENC_LSB) * M.ENC_LSB
        cmd = eng.step(q, qd, ra, rb)
        w_eff = w0 * np.maximum(0.05, 1.0 - np.abs(load) / tauS)
        th_g = th_g + dt * np.clip((cmd - th_g) / M.TS, -w_eff, w_eff)
        if tt[k] >= t0:
            e = q - th_out; err2 += e ** 2; n += 1; atl += np.abs(eng.u) >= 0.999 * eng.limit
            big = np.abs(e) > worst; worst = np.where(big, np.abs(e), worst); tw = np.where(big, tt[k] - t0, tw)
    rms = np.sqrt(err2 / n) / M.DEG
    for i, a in enumerate(arms):
        s = slice(i * J, (i + 1) * J)
        print(f'seed {seed} {a:5s} rms deg', np.round(rms[s], 4), '| at limit', np.round(atl[s] / n, 3), '| worst deg', np.round(worst[s] / M.DEG, 2), 'at s', np.round(tw[s], 1))
for s in (2, 3): go(s)
