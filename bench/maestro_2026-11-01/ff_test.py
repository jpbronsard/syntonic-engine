"""Servo-lag feedforward on the MAESTRO 1 servo bench: command = plan + servo_tau x planned velocity + syntonic
correction. servo_tau (25 ms) is the authority fact the engine already takes; nothing new is tuned.
Usage: python ff_test.py <repo_root> <excerpt> <cond: steady|step> <seed> [DT FRAME]"""
import os, sys
import numpy as np
from scipy.interpolate import CubicSpline
ROOT, NAME, COND, SEED = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
DT = float(sys.argv[5]) if len(sys.argv) > 5 else 1e-3
FRAME = float(sys.argv[6]) if len(sys.argv) > 6 else 3e-3
sys.path.insert(0, os.path.join(ROOT, 'python')); sys.path.insert(0, os.path.join(ROOT, 'examples'))
from syntonic_engine import SyntonicEngine
from syntonic_engine.engine import kappa_H
from arm_fk import tip

DEG = np.pi / 180
Z = np.load(os.path.join(ROOT, 'examples', 'data', 'maestro_excerpts.npz'))
TAU_S = Z['stall_torque']; W0 = Z['no_load_speed']
ENC_SD, ENC_LSB = float(Z['enc_sd_rad']), float(Z['enc_lsb_rad'])
TS, K_INV, HB = 0.025, (1.0 * DEG) / TAU_S, 0.5 * 0.3 * DEG
EVERY = int(round(0.01 / DT))


class FFEngine(SyntonicEngine):
    """The published engine, unchanged, except that the command sent at each frame adds the servo-lag lead
    servo_tau x planned velocity (the plan's own derivative). The internal servo model sees the same command."""
    def step_ff(self, q_plan, qd_plan, read_a, read_b):
        dt = self.dt
        q = np.asarray(q_plan, float); a = np.asarray(read_a, float); b = np.asarray(read_b, float)
        meas = 0.5 * (a + b); e = q - meas; d = (a - b) / np.sqrt(2.0)
        self.n2 += self.a_w * (d * d - self.n2)
        r = (q - self.th_ff) + (meas - self.th_m)
        fit = self._identify_raw(r) if self.identification == 'raw' else self._identify_allan(r)
        if fit is not None:
            Hn, lam, valid = fit
            s2 = 0.5 * self.n2 * dt
            ts = kappa_H(Hn) * (np.maximum(s2, 1e-30) / np.maximum(lam, 1e-30)) ** (1.0 / (2 * Hn + 1))
            ts = np.clip(ts, self.floor, self.tau_max)
            self.tau_target = np.where(valid, ts, self.tau_target); self.H = np.where(valid, Hn, self.H)
            rr = 1.0 / self.tau; dr_max = self.tau_rate * dt / self.tau ** 2
            rr = rr + np.clip(1.0 / self.tau_target - rr, -dr_max, dr_max); self.tau = 1.0 / rr
        g = 1.0 - np.exp(-dt / self.tau)
        self.u = np.clip(self.u + g * e, -self.limit, self.limit)
        lead = self.servo_tau * np.asarray(qd_plan, float)
        if self.k % self.frame_every == 0:
            self.cmd = q + lead + self.u
        # the servo model and the plan model both see the lead, so the residual stays exogenous
        self.th_m = self.th_m + dt * np.clip((self.cmd - self.th_m) / self.servo_tau, -self.w0, self.w0)
        self.th_ff = self.th_ff + dt * np.clip((q + lead - self.th_ff) / self.servo_tau, -self.w0, self.w0)
        self.k += 1
        return self.cmd


def run(ff):
    t0 = float(Z[f'{NAME}_t0']); tt = np.arange(t0 - 3.0, t0 + 30.0, DT)
    S = CubicSpline(Z[f'{NAME}_t30'], Z[f'{NAME}_q30'], axis=0)
    q = S(tt); qd = S(tt, 1); tau = CubicSpline(Z[f'{NAME}_t30'], Z[f'{NAME}_tau30'], axis=0)(tt)
    sd = np.full(len(tt), ENC_SD)
    if COND == 'step': sd[(tt >= t0 + 10) & (tt < t0 + 20)] *= np.sqrt(10.0)
    rng = np.random.default_rng(SEED); J = q.shape[1]
    th_g = q[0].copy()
    eng = FFEngine(q[0], DT, FRAME, W0, identification='raw')
    keep = []; err2 = np.zeros(J); m = 0
    for k in range(len(tt)):
        load = tau[k]
        th_out = th_g - load * K_INV - HB * np.tanh(load / (0.02 * TAU_S))
        a = np.round((th_out + rng.normal(0.0, 1.0, J) * sd[k]) / ENC_LSB) * ENC_LSB
        b = np.round((th_out + rng.normal(0.0, 1.0, J) * sd[k]) / ENC_LSB) * ENC_LSB
        cmd = eng.step_ff(q[k], qd[k] if ff else np.zeros(J), a, b)
        w_eff = W0 * np.maximum(0.05, 1.0 - np.abs(load) / TAU_S)
        th_g = th_g + DT * np.clip((cmd - th_g) / TS, -w_eff, w_eff)
        if tt[k] >= t0: err2 += (q[k] - th_out) ** 2; m += 1
        if k % EVERY == 0 and tt[k] >= t0 - 1e-9: keep.append((th_out[:6].copy(), q[k, :6].copy()))
    e = np.array([tip(x) - tip(y) for x, y in keep])
    return 1e3 * np.sqrt(np.mean(np.sum(e ** 2, 1))), float(np.mean(np.sqrt(err2 / m)[:6]) / DEG)


t0, j0 = run(False); t1, j1 = run(True)
print(f'{NAME} {COND:6s} seed {SEED} frame {FRAME*1e3:.0f} ms loop {1/DT:.0f} Hz | tip: engine {t0:.3f} mm, '
      f'engine+lead {t1:.3f} mm ({1 - t1 / t0:+.0%}) | joint: {j0:.4f} -> {j1:.4f} deg')
