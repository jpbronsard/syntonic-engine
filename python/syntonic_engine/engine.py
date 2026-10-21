# Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
# Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
"""Syntonic engine: a joint correction loop whose time window follows Theorem U' of the Syntony Principle.

Each loop tick, for every joint:
  1. read the encoder twice (two independent reads of the same angle);
  2. measure the encoder noise from the difference of the two reads;
  3. measure how the tracking residual moves, over lags of 2 to 64 ms (its structure function), and fit its
     smoothness H and intensity lambda_H;
  4. set the window  tau* = kappa_H (sigma^2 / lambda_H)^(1 / (2H + 1)),  kappa_H = (2H Gamma(2H + 1))^(-1/(2H + 1)),
     bounded below by the loop delay (never react faster than the loop can respond);
  5. integrate the error with that window and add it to the planned position when the servo frame is sent.

Nothing is tuned. The engine needs three facts about the hardware: the loop period, the servo frame period, and the
servo's nominal response (reaction time and no-load speed, from its data sheet). This is the engine of the MAESTRO 1
joint bench (flow_v5_uprime, 19-21 Oct 2026), extracted without change to the arithmetic. The same engine is in
cpp/SyntonicEngine.h and in the Rust crate rust/; both agree with this module to about 1e-14 rad.
"""
from math import gamma as _gamma
import numpy as np

__all__ = ['SyntonicEngine', 'kappa_H', 'loop_delay', 'window_floor']


def kappa_H(H):
    """Kernel prefactor of Theorem U' for an exponential window on a signal of structure function lambda_H u^(2H)."""
    H = np.asarray(H, dtype=float)
    return (2 * H * np.vectorize(_gamma)(2 * H + 1)) ** (-1.0 / (2 * H + 1))


def loop_delay(dt, frame, sensor_delay=1e-4):
    """Delay from an encoder read to the servo acting on it: half the servo frame (zero-order hold), half a loop
    tick (age of the latest correction when an asynchronous frame takes it), and the sensor's own delay
    (AS5048: about 0.1 ms)."""
    return frame / 2 + dt / 2 + sensor_delay


def window_floor(dt, frame, sensor_delay=1e-4):
    """Lower bound of the window: the larger of the loop delay and the discretisation limit dt / tau <= 0.45
    (TH41-A). 2.22 ms for a 1 kHz loop driving 333 Hz servos."""
    return max(loop_delay(dt, frame, sensor_delay), dt / 0.45)


class SyntonicEngine:
    """Per-joint syntonic correction for n joints, vectorized.

    Parameters
    ----------
    q0 : array (n,)          initial joint positions (rad)
    dt : float               loop period (s), e.g. 0.001
    frame : float            servo command frame period (s), e.g. 0.003 for 333 Hz PWM
    w0 : array (n,) or float no-load servo speed (rad/s), from the data sheet
    servo_tau : float        nominal servo reaction time (s); 0.025 when the data sheet gives none
    floor : float or None    window lower bound (s); None = window_floor(dt, frame, sensor_delay)
    sensor_delay : float     encoder delay (s)
    tau_max : float          window upper bound (s)
    tau_rate : float         largest change of the window, in seconds per second (H4: below 1 / (4 tau_max))
    lags : tuple             lags of the structure function (s)
    snr_gate : float         a lag is used only where its structure function stands this many noise floors above it
    ema : float              averaging time of the noise and structure-function estimates (s)
    warmup : float           data collected before the window is first set (s)
    limit : float            correction limit (rad)
    """

    def __init__(self, q0, dt, frame, w0, servo_tau=0.025, floor=None, sensor_delay=1e-4, tau_max=0.5,
                 tau_rate=0.4, lags=(0.002, 0.004, 0.008, 0.016, 0.032, 0.064), snr_gate=2.0, ema=1.0,
                 warmup=1.0, limit=np.radians(5.0)):
        q0 = np.asarray(q0, dtype=float); n = q0.shape[0]
        self.n, self.dt = n, float(dt)
        self.frame_every = max(1, int(round(frame / dt)))
        self.w0 = np.broadcast_to(np.asarray(w0, dtype=float), (n,)).copy()
        self.servo_tau = float(servo_tau)
        self.floor = window_floor(dt, frame, sensor_delay) if floor is None else float(floor)
        self.tau_max, self.tau_rate, self.snr_gate, self.limit = float(tau_max), float(tau_rate), float(snr_gate), float(limit)
        self.lags = np.unique(np.maximum(1, np.round(np.asarray(lags) / dt))).astype(int)
        self.logl = np.log(self.lags * dt)
        self.a_w = dt / ema
        self.warm = int(round(warmup / dt))
        # state
        self.k = 0
        self.u = np.zeros(n); self.cmd = q0.copy()
        self.th_m = q0.copy(); self.th_ff = q0.copy()
        self.tau = np.full(n, 0.5); self.H = np.full(n, 0.5); self.tau_target = np.full(n, 0.5)
        self.n2 = np.zeros(n)
        self.r_hist = np.zeros((self.lags.max() + 1, n)); self.Dk = np.zeros((len(self.lags), n))

    def step(self, q_plan, read_a, read_b, frame=None):
        """One loop tick. q_plan: planned positions now (rad); read_a, read_b: two independent encoder reads of the
        joint angles this tick (rad); frame: True when a servo frame is sent this tick (None: every frame period,
        counted from the first tick). Returns the command to send (rad), held between frames."""
        dt = self.dt
        q = np.asarray(q_plan, dtype=float)
        a = np.asarray(read_a, dtype=float); b = np.asarray(read_b, dtype=float)
        meas = 0.5 * (a + b); e = q - meas
        d = (a - b) / np.sqrt(2.0)
        self.n2 += self.a_w * (d * d - self.n2)                              # per-read noise variance
        r = (q - self.th_ff) + (meas - self.th_m)                            # residual from exogenous channels
        self.r_hist = np.roll(self.r_hist, 1, axis=0); self.r_hist[0] = r
        if self.k >= self.lags.max():
            for i, lg in enumerate(self.lags):
                dd = r - self.r_hist[lg]; self.Dk[i] += self.a_w * (dd * dd - self.Dk[i])
        if self.k > self.warm:
            D = self.Dk - self.n2[None, :]                                    # noise floor removed
            ok_g = D > self.snr_gate * self.n2[None, :]                       # lags above the floor
            ok_u = D > 0                                                      # fallback
            use_g = ok_g.sum(0) >= 2
            ok = np.where(use_g[None], ok_g, ok_u)
            nok = ok.sum(0)
            logD = np.where(ok, np.log(np.maximum(D, 1e-30)), 0.0)
            W = ok.astype(float); L = self.logl[:, None]
            Sw = W.sum(0); Sx = (W * L).sum(0); Sy = (W * logD).sum(0)
            Sxx = (W * L ** 2).sum(0); Sxy = (W * L * logD).sum(0)
            den = Sw * Sxx - Sx ** 2
            slope = np.where(den > 1e-12, (Sw * Sxy - Sx * Sy) / np.maximum(den, 1e-12), 1.0)
            Hn = np.clip(0.5 * slope, 0.05, 1.0)                               # smoothness H
            valid = nok >= 2
            s2 = 0.5 * self.n2 * dt                                            # noise PSD of the pair mean
            lam = np.exp((W * (logD - 2 * Hn[None] * L)).sum(0) / np.maximum(Sw, 1))
            ts = kappa_H(Hn) * (np.maximum(s2, 1e-30) / np.maximum(lam, 1e-30)) ** (1.0 / (2 * Hn + 1))
            ts = np.clip(ts, self.floor, self.tau_max)
            self.tau_target = np.where(valid, ts, self.tau_target); self.H = np.where(valid, Hn, self.H)
            rr = 1.0 / self.tau; dr_max = self.tau_rate * dt / self.tau ** 2
            rr = rr + np.clip(1.0 / self.tau_target - rr, -dr_max, dr_max); self.tau = 1.0 / rr
        g = 1.0 - np.exp(-dt / self.tau)
        self.u = np.clip(self.u + g * e, -self.limit, self.limit)
        if (self.k % self.frame_every == 0) if frame is None else frame:
            self.cmd = q + self.u
        self.th_m = self.th_m + dt * np.clip((self.cmd - self.th_m) / self.servo_tau, -self.w0, self.w0)
        self.th_ff = self.th_ff + dt * np.clip((q - self.th_ff) / self.servo_tau, -self.w0, self.w0)
        self.k += 1
        return self.cmd
