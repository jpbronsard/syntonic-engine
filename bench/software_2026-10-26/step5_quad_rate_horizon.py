#!/usr/bin/env python3
# Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche). Apache-2.0.
"""Step 5, 27 Oct 2026: the rate-loop horizon of the dRehmFlight_syntonic cascade, raw proxy against measurement layer.

PRE-REGISTERED. Hypotheses, cells, criteria and decision rule are written here before the first run; the SHA-256 of
this file is printed first. Nothing below the hypotheses was changed after seeing results.

Question. The graft (Syntonic.h, 2026-05) feeds the square-root law tau* = kappa sqrt(sigma^2 / lambda) with a proxy:
sigma = std of the raw loop error, lambda = RMS of its time derivative. Under gyro noise this proxy collapses to the
floor (for white noise sigma/lambda = dt/sqrt 2). The measurement layer of engine 0.2.0 feeds the same law with the
read-noise variance (known) and the innovation read from the structure function of block means with the noise floor
subtracted. Does the measurement layer keep tau* alive and near the best fixed horizon, where the proxy sits at the
floor? And does a horizon law placed inside a feedback loop need a cap derived from the loop?

Plant (sim_compare.c of the graft folder): one axis, omega_dot = B u - D omega + d, B = 4000 deg/s^2 per unit of
differential, D = 2 /s, angle integrates omega. Loop 2 kHz (dt = 0.5 ms). u = AUTHORITY * c, |u| <= AUTHORITY = 0.25.
Controller: the graft cascade. Outer angle loop e_a = (theta_sp - theta_m)/30 deg through a first-order filter with
a FIXED 2 ms horizon (controlled in every variant: the question is the inner loop), rate_sp = 0.9 * c_a (unit = 250
deg/s). Inner rate loop e_r = rate_sp - gyro/250 through a first-order filter with the horizon under test, c_r its
output, u = 0.25 c_r. This is a cascade of two proportional loops when both horizons sit at their floors: inner pole
B*0.25/250 = 4 /s, outer pole 0.9*250/30 = 7.5 /s.
Signals: setpoint square wave +-10 deg, step every 5 s; disturbance torque d: Ornstein-Uhlenbeck, correlation 1 s,
std 60 deg/s^2 (about 2 deg of static droop); gyro noise white, std SIGMA_G in {1, 3, 10} deg/s per sample; angle
measurement noise white 0.05 deg. Run 40 s, scored from 5 s. Seeds 0..4 (same noise draws across variants).
Scores: J1 = RMS(theta - theta_sp) over the scored window (deg); J2 = RMS of the per-tick increment of u (motor
jitter, dimensionless). A run whose |theta| exceeds 90 deg is unstable (J1 = inf).

Inner-horizon variants:
  FIXED(tau)    tau in geomspace(2 ms, 0.5 s, 13). tau_opt = argmin J1 (median J1 over seeds).
  RAW           Syntonic.h as written: EMA windows max(tau*, 50 ms), floor 2 ms, cap 2 s, meta-kappa at 5 Hz on this
                axis (eta 0.05, step 5e-4, [0.9, 1.1], Phi target 0.95), kappa 1 at start.
  ML-half       measurement layer, H fixed at 1/2 (the square-root law): sigma^2 = known noise variance of e_r times
                dt; lambda from the structure function of adjacent block means (scales 2 ms .. 4.096 s), noise floor
                n2 dt / l subtracted, six finest scales clearing the SNR gate 2, lambda = exp mean(log D - log l)/C(1/2);
                tau clipped to [2 ms, cap], rate limit 0.4 as in the engine; identification refreshed every 20 ticks.
                cap = 0.5 s (the servo engine's default).
  ML-half-cap   same, cap = 1/(3 * 7.5 /s) = 44 ms, a cap derived from the loop (inner lag at most a third of the
                outer time constant).
  ML-U-cap      same as ML-half-cap with H fitted (Theorem U': tau = kappa_H (sigma^2/lambda_H)^(1/(2H+1))).
  ML-half-cap with the noise variance mis-stated by x0.5 and x2 (SIGMA_G = 3 only).

Hypotheses (median over seeds unless stated):
  R1  RAW's inner tau* is within 1.5x of the floor for at least 90 % of the scored time at every SIGMA_G.
  R2  ML-half's median tau* is at least 3x the floor at every SIGMA_G, and its J1 regret against the FIXED optimum
      is at most 10 % at every SIGMA_G.
  R3  On J2, ML-half-cap <= 0.7 x RAW at SIGMA_G = 3 and 10.
  R4  Without the loop cap (ML-half, cap 0.5 s), at least one cell shows a J1 regret above 25 % or an unstable run,
      and ML-half-cap shows none: a horizon law inside a loop needs a cap from the loop.
  R5  (exploratory) ML-U-cap within 10 % of ML-half-cap on J1: the exponent does not matter here (H near 1/2).
  R6  (exploratory) a 2x error on the noise variance moves ML-half-cap's median tau* by less than 1.5x (the law
      predicts sqrt 2).
Decision: the measurement-layer claim for the quad stands if R1, R2 (with the cap variant allowed to be the one that
meets the regret) and R3 hold; it fails if ML-half-cap misses the 10 % regret at any SIGMA_G. R4 decides whether the
cap goes into the protocol as a fourth hardware-derived quantity.

AMENDMENT BEFORE THE FULL RUN (after a smoke test: one seed, SIGMA_G = 3, B = 4000, five single runs, no sweep).
The smoke test showed the J1 optimum of the FIXED sweep sitting at the floor for B = 4000 (inner pole 4 /s, slower
than the outer 7.5 /s: any inner lag costs), with the measurement layer also at the floor. A floor optimum cannot
tell a law that sits at the floor for the right reason from a proxy that sits there because of noise. Two changes,
written before the first full run:
  - B becomes a cell: B in {4000, 30000} deg/s^2 per unit (the second is a stiff small quad, inner pole 30 /s).
  - R2 is restated in the form of the servo bench: ML-half-cap's J1 regret <= 10 % in every cell, and in every cell
    whose FIXED optimum is >= 3x the floor, ML-half-cap's median tau* within 1.5x of tau_opt. The "tau_med >= 3x
    floor everywhere" clause is dropped (it presumed where the optimum is).
  - A steady-state score J1h (RMS angle error over the last 3 s of each 5 s segment) is reported, not decided on.
Everything else unchanged. The decision rule reads R2 in its restated form.
"""
import hashlib, json, math, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
print('SHA-256 of this file:', hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest()[:16])

DT, FS = 0.0005, 2000.0
T_RUN, T_SCORE = 40.0, 5.0
BS, DAMP = [4000.0, 30000.0], 2.0
AUTH, ANGLE_RANGE, GYRO_RANGE, ANGLE_TO_RATE = 0.25, 30.0, 250.0, 0.9
SP_AMP, SP_HALF = 10.0, 5.0
OU_TAU, OU_STD = 1.0, 60.0
ANGLE_NOISE = 0.05
SIGMAS = [1.0, 3.0, 10.0]
SEEDS = [0, 1, 2, 3, 4]
FLOOR = 0.002
TAUS = np.geomspace(0.002, 0.5, 13)
OUTER_TAU = 0.002
K_OUTER = ANGLE_TO_RATE * GYRO_RANGE / ANGLE_RANGE          # 7.5 /s
CAP_LOOP = 1.0 / (3.0 * K_OUTER)                              # 44 ms
OUT = os.path.join(HERE, 'step5_quad_rate_horizon_results.txt')

def C_H(H):
    return (2.0 ** (2 * H + 2) - 4.0) / ((2 * H + 1) * (2 * H + 2))

def kappa_H(H):
    return (2 * H * math.gamma(2 * H + 1)) ** (-1.0 / (2 * H + 1))

# ------------------------------------------------------------------------------------------------------ signals
def make_signals(sigma_g, seed):
    rng = np.random.default_rng(1000 + seed)
    n = int(T_RUN / DT)
    t = np.arange(n) * DT
    sp = np.where((t // SP_HALF) % 2 == 0, SP_AMP, -SP_AMP)
    a = math.exp(-DT / OU_TAU); s = OU_STD * math.sqrt(1 - a * a)
    d = np.empty(n); x = 0.0
    w = rng.standard_normal(n)
    for k in range(n):
        x = a * x + s * w[k]; d[k] = x
    gn = rng.standard_normal(n) * sigma_g
    an = rng.standard_normal(n) * ANGLE_NOISE
    return t, sp, d, gn, an

# ------------------------------------------------------------------------------------------------- horizon laws
class Fixed:
    def __init__(self, tau): self.tau = tau
    def update(self, e): return self.tau

class Raw:
    """Syntonic.h, syn_update + syn_meta_kappa on one axis."""
    def __init__(self, cap=2.0):
        self.cap = cap; self.kappa = 1.0; self.tau = FLOOR; self.seeded = False; self.n = 0
        self.mu = self.e2 = self.de2 = self.mse = 0.0; self.prev = 0.0; self.integ = 0.0; self.phi = 0.5
        self.meta_ticks = int(round(0.2 / DT)); self.k = 0
    def update(self, e):
        self.n += 1
        if not self.seeded:
            self.mu, self.e2, self.de2, self.mse, self.prev, self.seeded = e, e * e, 0.0, 0.0, e, True
            return self.tau
        de = (e - self.prev) / DT; self.prev = e
        tau_obs = max(self.tau, 0.05)
        bs = math.exp(-DT / (2 * tau_obs)); bl = math.exp(-DT / tau_obs)
        self.mu = bs * self.mu + (1 - bs) * e; self.e2 = bs * self.e2 + (1 - bs) * e * e
        self.de2 = bl * self.de2 + (1 - bl) * de * de
        var = max(self.e2 - self.mu * self.mu, 1e-9); sig = math.sqrt(var); lam = math.sqrt(max(self.de2, 1e-9))
        self.tau = min(max(self.kappa * sig / lam, FLOOR), self.cap)
        # the graft's own first-order relaxation, only to form Phi (the loop filter is applied by the caller
        # with the same tau, so the two integrals coincide)
        al = 1 - math.exp(-DT / self.tau); self.integ += al * (e - self.integ)
        r = e - self.integ; self.mse = bl * self.mse + (1 - bl) * r * r
        if self.n >= 10:
            self.phi = min(max(1 - self.mse / max(var, 1e-4), 0.05), 1.0)
        self.k += 1
        if self.k % self.meta_ticks == 0:
            dk = min(max(0.05 * (self.phi - 0.95) * 0.2, -5e-4), 5e-4)
            self.kappa = min(max(self.kappa + dk, 0.9), 1.1)
        return self.tau

class MeasurementLayer:
    """Engine 0.2.0 identification on block means, fed a known read-noise variance; H fixed at 1/2 or fitted."""
    SCALES_S = [2 ** i * 0.001 for i in range(1, 13)]
    def __init__(self, noise_var, cap, fitH=False, refresh=20, snr_gate=2.0, ema=1.0, warm=1.0, tau_rate=0.4, finest=6):
        self.n2 = 2.0 * noise_var                     # single read: n2 is defined for the pair mean (var n2/2)
        self.cap, self.fitH, self.refresh = cap, fitH, refresh
        self.scales = np.unique(np.maximum(1, np.round(np.asarray(self.SCALES_S) / DT))).astype(int)
        self.logl = np.log(self.scales * DT); m = len(self.scales)
        self.acc = np.zeros(m); self.cnt = np.zeros(m, dtype=int); self.prev = np.full(m, np.nan)
        self.DA = np.zeros(m); self.nblk = np.zeros(m, dtype=int)
        self.a_blk = np.array([min(1.0, DT * s / max(ema, 8 * DT * s)) for s in self.scales])
        self.snr_gate, self.warm, self.tau_rate, self.finest = snr_gate, int(round(warm / DT)), tau_rate, finest
        self.k = 0; self.tau = cap; self.tau_target = cap; self.H = 0.5
    def _fit(self):
        floor = self.n2 * DT / (self.scales * DT)
        D = np.where(self.nblk >= 2, self.DA - floor, -1.0)
        ok_g = D > self.snr_gate * floor; ok_u = D > 0
        ok = ok_g if ok_g.sum() >= 2 else ok_u
        ok &= np.cumsum(ok) <= self.finest
        if ok.sum() < 2: return None
        L = self.logl[ok]; y = np.log(D[ok])
        if self.fitH:
            slope = np.polyfit(L, y, 1)[0]; H = min(max(0.5 * slope, 0.05), 1.0)
        else:
            H = 0.5
        lam = math.exp(np.mean(y - 2 * H * L)) / C_H(H)
        return H, lam
    def update(self, e):
        self.acc += e; self.cnt += 1
        done = self.cnt == self.scales
        if done.any():
            for i in np.flatnonzero(done):
                mean = self.acc[i] / self.scales[i]
                if not np.isnan(self.prev[i]):
                    dd = mean - self.prev[i]; aw = self.a_blk[i] if self.nblk[i] > 0 else 1.0
                    self.DA[i] += aw * (dd * dd - self.DA[i]); self.nblk[i] += 1
                self.prev[i] = mean; self.acc[i] = 0.0; self.cnt[i] = 0
        self.k += 1
        if self.k > self.warm and self.k % self.refresh == 0:
            fit = self._fit()
            if fit is not None:
                H, lam = fit; s2 = 0.5 * self.n2 * DT
                ts = kappa_H(H) * (max(s2, 1e-30) / max(lam, 1e-30)) ** (1.0 / (2 * H + 1))
                self.tau_target = min(max(ts, FLOOR), self.cap); self.H = H
        # smooth move of the window, as in the engine (rate limit on 1/tau)
        rr = 1.0 / self.tau; dr_max = self.tau_rate * DT / self.tau ** 2
        rr += min(max(1.0 / self.tau_target - rr, -dr_max), dr_max); self.tau = 1.0 / rr
        return self.tau

# ---------------------------------------------------------------------------------------------------- the loop
def run(law, sigma_g, seed, B):
    t, sp, d, gn, an = make_signals(sigma_g, seed)
    n = len(t); theta = 0.0; omega = 0.0; ca = 0.0; cr = 0.0; u_prev = 0.0
    al_a = 1 - math.exp(-DT / OUTER_TAU)
    k0 = int(T_SCORE / DT)
    err2 = 0.0; du2 = 0.0; cnt = 0; taus = []; errh2 = 0.0; cnth = 0
    for k in range(n):
        th_m = theta + an[k]; gy = omega + gn[k]
        ea = (sp[k] - th_m) / ANGLE_RANGE
        ca += al_a * (ea - ca); ca = min(max(ca, -1.0), 1.0)
        rate_sp = ANGLE_TO_RATE * ca
        er = rate_sp - gy / GYRO_RANGE
        tau = law.update(er)
        al = 1 - math.exp(-DT / tau)
        cr += al * (er - cr); cr = min(max(cr, -1.0), 1.0)
        u = AUTH * cr
        omega += DT * (B * u - DAMP * omega + d[k]); theta += DT * omega
        if abs(theta) > 90.0:
            return dict(J1=float('inf'), J1h=float('inf'), J2=float('inf'), tau_med=float('nan'), floor_frac=float('nan'), unstable=True)
        if k >= k0:
            e = theta - sp[k]; err2 += e * e; du2 += (u - u_prev) ** 2; cnt += 1
            if (t[k] % SP_HALF) >= 2.0: errh2 += e * e; cnth += 1
            if k % 20 == 0: taus.append(tau)
        u_prev = u
    taus = np.asarray(taus)
    return dict(J1=math.sqrt(err2 / cnt), J1h=math.sqrt(errh2 / cnth), J2=math.sqrt(du2 / cnt), tau_med=float(np.median(taus)),
                floor_frac=float(np.mean(taus <= 1.5 * FLOOR)), unstable=False)

def med(xs): return float(np.median(xs))

def main():
    rows = {}
    lines = []
    def out(s):
        print(s); lines.append(s)
    opt = {}; summary = {}
    for B in BS:
        out(f'\n==== B = {B:.0f} deg/s^2 per unit (inner pole {B*AUTH/GYRO_RANGE:.0f} /s, outer {K_OUTER:.1f} /s) ====')
        out('== FIXED sweep: median J1 (deg), J1h, J2 over seeds ==')
        out('sigma_g  ' + ' '.join(f'{tau*1e3:7.1f}' for tau in TAUS) + '   (tau, ms)')
        for sg in SIGMAS:
            J1 = []; J1h = []; J2 = []
            for tau in TAUS:
                r = [run(Fixed(tau), sg, s, B) for s in SEEDS]
                rows[f'FIXED,{B:.0f},{sg},{tau:.5f}'] = r
                J1.append(med([x['J1'] for x in r])); J1h.append(med([x['J1h'] for x in r])); J2.append(med([x['J2'] for x in r]))
            i = int(np.argmin(J1)); opt[(B, sg)] = (TAUS[i], J1[i])
            out(f'{sg:6.1f}  ' + ' '.join(f'{v:7.3f}' for v in J1) + f'   J1; opt tau {TAUS[i]*1e3:.1f} ms, J1 {J1[i]:.3f}')
            out(f'        ' + ' '.join(f'{v:7.3f}' for v in J1h) + f'   J1h (best at {TAUS[int(np.argmin(J1h))]*1e3:.1f} ms)')
            out(f'        ' + ' '.join(f'{v:7.5f}' for v in J2) + '   J2')
        out('== Laws: median over seeds ==')
        out(f"{'variant':14s} {'sigma_g':>7s} {'J1':>7s} {'regret':>7s} {'J1h':>7s} {'J2':>8s} {'J2/RAW':>7s} {'tau_med':>8s} {'x floor':>7s} {'x opt':>6s} {'at floor':>8s} {'unst':>4s} {'H':>5s}")
        for sg in SIGMAS:
            res = {}
            variants = [('RAW', lambda nv: Raw()),
                        ('ML-half', lambda nv: MeasurementLayer(nv, cap=0.5)),
                        ('ML-half-cap', lambda nv: MeasurementLayer(nv, cap=CAP_LOOP)),
                        ('ML-U-cap', lambda nv: MeasurementLayer(nv, cap=CAP_LOOP, fitH=True))]
            if sg == 3.0:
                variants += [('ML-cap-n2x0.5', lambda nv: MeasurementLayer(0.5 * nv, cap=CAP_LOOP)),
                             ('ML-cap-n2x2', lambda nv: MeasurementLayer(2.0 * nv, cap=CAP_LOOP))]
            nv = (sg / GYRO_RANGE) ** 2
            for name, mk in variants:
                laws = [mk(nv) for s in SEEDS]
                r = [run(law, sg, s, B) for law, s in zip(laws, SEEDS)]
                rows[f'{name},{B:.0f},{sg}'] = r; res[name] = r
                J1 = med([x['J1'] for x in r]); J1h = med([x['J1h'] for x in r]); J2 = med([x['J2'] for x in r])
                tm = med([x['tau_med'] for x in r]); ff = med([x['floor_frac'] for x in r])
                uns = sum(x['unstable'] for x in r)
                Hm = med([getattr(l, 'H', float('nan')) for l in laws])
                regret = J1 / opt[(B, sg)][1] - 1 if math.isfinite(J1) else float('inf')
                j2raw = J2 / med([x['J2'] for x in res['RAW']]) if name != 'RAW' else 1.0
                summary[(name, B, sg)] = dict(J1=J1, J1h=J1h, regret=regret, J2=J2, J2raw=j2raw, tau_med=tm, floor_frac=ff, unstable=uns, xopt=tm / opt[(B, sg)][0])
                out(f'{name:14s} {sg:7.1f} {J1:7.3f} {regret:+7.1%} {J1h:7.3f} {J2:8.5f} {j2raw:7.2f} {tm*1e3:6.1f}ms {tm/FLOOR:7.1f} {tm/opt[(B, sg)][0]:6.2f} {ff:8.0%} {uns:4d} {Hm:5.2f}')
    cells = [(B, sg) for B in BS for sg in SIGMAS]
    out('\n== Hypotheses ==')
    r1 = all(summary[('RAW', B, sg)]['floor_frac'] >= 0.9 for B, sg in cells)
    out(f"R1 RAW at the floor >= 90 % of the time in every cell: {'held' if r1 else 'failed'}")
    r2a = all(summary[('ML-half-cap', B, sg)]['regret'] <= 0.10 for B, sg in cells)
    off = [(B, sg) for B, sg in cells if opt[(B, sg)][0] >= 3 * FLOOR]
    r2b = all(1 / 1.5 <= summary[('ML-half-cap', B, sg)]['xopt'] <= 1.5 for B, sg in off)
    out(f"R2 ML-half-cap regret <= 10 % in every cell: {'held' if r2a else 'failed'}; cells with optimum >= 3x floor: "
        f"{[(int(B), sg) for B, sg in off]}; tau_med within 1.5x of tau_opt there: {'held' if r2b else 'failed'} -> "
        f"{'held' if (r2a and r2b) else 'failed'}")
    r3 = all(summary[('ML-half-cap', B, sg)]['J2raw'] <= 0.7 for B in BS for sg in (3.0, 10.0))
    out(f"R3 ML-half-cap J2 <= 0.7 x RAW at 3 and 10 deg/s, both B: {'held' if r3 else 'failed'}")
    r4a = any(summary[('ML-half', B, sg)]['regret'] > 0.25 or summary[('ML-half', B, sg)]['unstable'] > 0 for B, sg in cells)
    r4b = all(summary[('ML-half-cap', B, sg)]['regret'] <= 0.25 and summary[('ML-half-cap', B, sg)]['unstable'] == 0 for B, sg in cells)
    out(f"R4 no-cap shows regret > 25 % or instability somewhere: {'yes' if r4a else 'no'}; cap variant clean everywhere: "
        f"{'yes' if r4b else 'no'} -> {'held' if (r4a and r4b) else 'failed'}")
    r5 = all(abs(summary[('ML-U-cap', B, sg)]['J1'] / summary[('ML-half-cap', B, sg)]['J1'] - 1) <= 0.10 for B, sg in cells)
    out(f"R5 (expl.) ML-U-cap within 10 % of ML-half-cap on J1 in every cell: {'held' if r5 else 'failed'}")
    for B in BS:
        t0 = summary[('ML-half-cap', B, 3.0)]['tau_med']
        rat = [summary[('ML-cap-n2x0.5', B, 3.0)]['tau_med'] / t0, summary[('ML-cap-n2x2', B, 3.0)]['tau_med'] / t0]
        out(f"R6 (expl.) B={B:.0f}: tau_med ratio for n2 x0.5, x2: {rat[0]:.2f}, {rat[1]:.2f} (law: 0.71, 1.41); within 1.5x: "
            f"{'held' if all(1/1.5 <= x <= 1.5 for x in rat) else 'failed'}")
    verdict = r1 and r2a and r2b and r3
    out(f"\nDecision: measurement-layer claim for the quad rate loop {'STANDS' if verdict else 'FAILS'}; "
        f"loop cap {'goes into the protocol' if (r4a and r4b) else 'not supported by this bench'}.")
    open(OUT, 'w').write('\n'.join(lines) + '\n')
    json.dump(rows, open(OUT.replace('.txt', '_rows.json'), 'w'))
    print('\nwritten', OUT)

if __name__ == '__main__':
    main()
