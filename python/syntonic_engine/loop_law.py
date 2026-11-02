# Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
# Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
"""The syntonic loop law: a joint correction whose PI gains are not tuned but computed, every 0.5 s, as the gains that
minimise the loop error predicted from statistics measured live. Written 2 Nov 2026 (exploratory development on seeds
91-93 of the published plant, /home/claude/work/goal/core/proto.py, classes SPIB..SPIO8); frozen here for the
pre-registered confirmation.

FACTS (measured once, as the published engine): loop period dt, servo frame period, servo time constant servo_tau
(data sheet; 25 ms), no-load speed w0 (data sheet), and the servo's internal delay `dead` (measured once at commissioning,
as the engine's README asks for the full loop delay). Design constants, declared: averaging time 1 s, warm-up 1 s,
noise gate 2, block lags 2 ms to 256 ms, at most 3 fitted lags, gain update every 0.5 s, phase margin >= 30 deg and
max |S| <= 2 (textbook robustness), candidate grid of (Kp, Ki).

EACH TICK, per joint:
  1. two encoder reads a, b: y = (a + b)/2, read-pair noise variance R = EMA_1s((a - b)^2 / 2) / 2 (the engine's n2 / 2);
  2. nominal servo models (the engine's): th_m driven by the commands actually sent, th_ff by the plan feed-forward
     alone, both delayed by `dead`; the error the loop must cancel is r = (q - th_ff) + (th_m - y)
     (plan-following lag of the nominal servo plus the output disturbance, same sign);
  3. identification (block-mean, Hadamard): second differences of adjacent block means of r at lags u = 2..256 ms,
     EMA'd, white-noise floor 6 R dt / u removed; fit S(u) = A u^gamma on the finest lags (<= 3) clear of the floor
     (gate 2) where S still grows; gamma clipped to [1, 3]; corner Oc = 1 / (largest lag up to which S grows);
  4. every 0.5 s: disturbance spectrum Phi(W) = A / K(gamma) x (W^2 + Oc^2)^(-(1+gamma)/2) (two-sided; K the
     block-mean constant, K(3) = 11/20 exact), read noise white with intensity r = R dt; among the candidate PI gains
     that keep the loop G = exp(-s D) / (servo_tau s + 1) robust (D = loop_delay(dt, frame) + dead), choose the one
     minimising J = int(|S|^2 Phi + |T|^2 r) dW / 2 pi;
  5. PI: the integral accumulates Ki e dt every tick (e = q - y; frozen while |Kp e + I| >= the authority limit
     servo_tau x w0); at each frame the proportional term uses the mean error over the ticks since the last frame;
     command = q + L x qd + clip(Kp e_mean + I, +-limit), L = the lead coefficient (shared with the rival PID).

In the package since 0.3.0 (in C++ and Rust since 0.3.1): this module is the frozen file
bench/goal_2026-11-02/servo/syntonic_loop_law.py (sha256 a7a234ccefb2c9bb963b87b9a3268493644127bfd43d44c387e1411b3910b2b3)
with this paragraph added; the code is unchanged. Its pre-registered confirmation against tuned PIDs on a realistic
MAESTRO 1 bench (servo frames 333 Hz, nine circumstances, fresh seeds) is bench/goal_2026-11-02/servo/confirm_servo.py,
with the verdicts in bench/goal_2026-11-02/NOTE.md. The development prototype named above stayed in the working
session and is not in the repository.
"""
import numpy as np

DEG = np.pi / 180
_KG = np.array([[1.0, 0.99999999987], [1.5, 0.74370490905], [2.0, 0.59538447860], [2.5, 0.52771124824], [3.0, 0.55]])


def loop_delay(dt, frame, sensor_delay=1e-4):
    return frame / 2 + dt / 2 + sensor_delay


def K_gamma(g):
    return np.interp(g, _KG[:, 0], _KG[:, 1])


class LoopLaw:
    """Frequency responses of the PI loop over the candidate grid, and the per-joint choice of the gains."""
    def __init__(self, servo_tau, D, pm_min=30.0, ms_max=2.0):
        self.W = np.logspace(-1, 4, 1500); dW = np.gradient(self.W)
        kp = np.concatenate([[0.0], np.geomspace(0.05, 12, 40)]); ki = np.geomspace(5, 5000, 60)
        KP, KI = np.meshgrid(kp, ki, indexing='ij'); self.KP, self.KI = KP.ravel(), KI.ravel()
        s = 1j * self.W[None, :]
        G = np.exp(-s * D) / (1 + s * servo_tau)
        L = (self.KP[:, None] + self.KI[:, None] / s) * G
        S = 1 / (1 + L); T = L / (1 + L)
        aL = np.abs(L); idx = np.argmax(aL < 1, axis=1)
        pm = 180 + np.degrees(np.angle(L[np.arange(len(L)), idx])); pm = np.where(pm > 180, pm - 360, pm)
        self.ok = (pm >= pm_min) & (np.abs(S).max(1) <= ms_max) & (aL[:, 0] > 1)
        self.MS = (np.abs(S) ** 2) * dW[None, :] / (2 * np.pi) * 2
        self.MTr = ((np.abs(T) ** 2) * dW[None, :] / (2 * np.pi) * 2).sum(1)

    def gains(self, A, gam, Oc, r):
        C = A / K_gamma(gam)
        Phi = C[:, None] * (self.W[None, :] ** 2 + Oc[:, None] ** 2) ** (-(1 + gam[:, None]) / 2)
        J = Phi @ self.MS.T + r[:, None] * self.MTr[None, :]
        J = np.where(self.ok[None, :], J, np.inf); b = np.argmin(J, axis=1)
        return self.KP[b], self.KI[b]


_CACHE = {}


class SyntonicLoopLaw:
    LAGS = (2, 4, 8, 16, 32, 64, 128, 256)

    def __init__(self, q0, w0, dt=1e-3, frame=3e-3, servo_tau=0.025, dead=0.0, lead=None, ema=1.0, warmup=1.0, gate=2.0,
                 nfit=3, every=500):
        q0 = np.asarray(q0, float); n = len(q0); self.n, self.dt = n, dt
        self.w0 = np.broadcast_to(np.asarray(w0, float), (n,)).copy(); self.ts = servo_tau; self.lim = servo_tau * self.w0
        D = loop_delay(dt, frame) + dead
        self.L = servo_tau + loop_delay(dt, frame) if lead is None else float(lead)
        key = (servo_tau, round(D, 12)); _CACHE.setdefault(key, LoopLaw(servo_tau, D)); self.LL = _CACHE[key]
        self.aw = dt / ema; self.warm = int(round(warmup / dt)); self.gate = gate; self.nfit = nfit; self.every = every
        nd = dead / dt; self.nd = int(np.floor(nd + 1e-9)); self.fd = nd - self.nd          # whole and fractional ticks
        self.pc = [q0.copy()] * (self.nd + 2); self.pf = [q0.copy()] * (self.nd + 2)
        self.th_m = q0.copy(); self.th_ff = q0.copy(); self.cmd = q0.copy(); self.cff = q0.copy(); self.n2 = np.zeros(n)
        self.R = np.full(n, 1e-14); self.k = 0; self.started = False
        self.Kp = np.zeros(n); self.Ki = np.full(n, 50.0); self.I = np.zeros(n); self.esum = np.zeros(n); self.ecnt = 0
        self.lags = np.array(self.LAGS); m = len(self.lags)
        self.acc = np.zeros((m, n)); self.cnt = np.zeros(m, int); self.m1 = np.full((m, n), np.nan); self.m2 = np.full((m, n), np.nan)
        self.S2 = np.zeros((m, n)); self.nb = np.zeros(m, int)
        self.ab = np.array([min(1.0, dt * L / max(ema, 8 * dt * L)) for L in self.lags])
        self.gam = np.full(n, 3.0); self.A = np.zeros(n); self.Oc = np.full(n, 1.0)

    def _identify(self, z):
        self.acc += z[None, :]; self.cnt += 1
        for i in np.flatnonzero(self.cnt == self.lags):
            mean = self.acc[i] / self.lags[i]
            if not np.isnan(self.m2[i, 0]):
                dd = mean - 2 * self.m1[i] + self.m2[i]
                a = self.ab[i] if self.nb[i] > 0 else 1.0
                self.S2[i] += a * (dd * dd - self.S2[i]); self.nb[i] += 1
            self.m2[i] = self.m1[i]; self.m1[i] = mean; self.acc[i] = 0.0; self.cnt[i] = 0
        if self.k > self.warm:
            u = self.lags * self.dt; fl = 6.0 * self.R[None, :] * self.dt / u[:, None]
            sig = self.S2 - fl
            valid = (sig > self.gate * fl) & (self.nb >= 2)[:, None]
            lsig = np.log(np.maximum(sig, 1e-30)); grow = np.vstack([np.ones((1, self.n), bool), lsig[1:] > lsig[:-1]])
            ok = valid & grow; ok = ok & (np.cumsum(ok, axis=0) <= self.nfit)
            lu = np.log(u)[:, None]; W = ok.astype(float)
            Sw = W.sum(0); Sx = (W * lu).sum(0); Sy = (W * lsig).sum(0); Sxx = (W * lu * lu).sum(0); Sxy = (W * lu * lsig).sum(0)
            den = Sw * Sxx - Sx * Sx
            g = np.where(den > 1e-12, (Sw * Sxy - Sx * Sy) / np.maximum(den, 1e-12), 3.0); g = np.clip(g, 1.0, 3.0)
            lA = (Sy - g * Sx) / np.maximum(Sw, 1)
            good = Sw >= 2
            self.gam = np.where(good, g, self.gam); self.A = np.where(good, np.exp(lA), self.A)
            self.started = self.started or bool(good.any())
            ug = np.full(self.n, u[-1])
            for j in range(self.n):
                v = np.flatnonzero(valid[:, j])
                if len(v) == 0: continue
                i = v[0]
                while i + 1 < len(u) and valid[i + 1, j] and lsig[i + 1, j] > lsig[i, j]: i += 1
                ug[j] = u[i]
            self.Oc = 1.0 / ug

    def _delayed(self, pipe, x):
        pipe.append(x.copy()); pipe.pop(0)
        if self.fd == 0: return pipe[-1 - self.nd]
        return (1 - self.fd) * pipe[-1 - self.nd] + self.fd * pipe[-2 - self.nd]

    def step(self, q, qd, a, b, send):
        dt = self.dt
        y = 0.5 * (a + b); dn = (a - b) / np.sqrt(2.0)
        self.n2 += self.aw * (dn * dn - self.n2); self.R = np.maximum(0.5 * self.n2, 1e-14)
        self._identify((q - self.th_ff) + (self.th_m - y))
        if self.started and (self.k % self.every == 0):
            self.Kp, self.Ki = self.LL.gains(self.A, self.gam, self.Oc, self.R * dt)
        e = q - y; self.esum += e; self.ecnt += 1
        self.I += np.where(np.abs(self.Kp * e + self.I) >= self.lim, 0.0, self.Ki * e * dt)
        if send:
            eb = self.esum / self.ecnt; self.esum[:] = 0.0; self.ecnt = 0
            self.cmd = q + self.L * qd + np.clip(self.Kp * eb + self.I, -self.lim, self.lim); self.cff = q + self.L * qd
        cm = self._delayed(self.pc, self.cmd); cf = self._delayed(self.pf, self.cff)
        self.th_m = self.th_m + dt * np.clip((cm - self.th_m) / self.ts, -self.w0, self.w0)
        self.th_ff = self.th_ff + dt * np.clip((cf - self.th_ff) / self.ts, -self.w0, self.w0)
        self.k += 1
        return self.cmd
