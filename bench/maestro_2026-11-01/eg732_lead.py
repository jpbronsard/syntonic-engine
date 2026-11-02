# Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
# Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
"""The e-geodesic cadence at the servo's own ceiling, with the servo lead. Written 1 Nov 2026, before any run of this
file except the smoke test; idea and go from JP ("On suit l'e-geodesique", then "Oui").

What changes from egeodesic_cadence.py (31 Oct):
  1. Engine 0.2.8: correction limit from the authority fact (servo_tau x no-load speed per joint), not 5 degrees.
  2. Ceiling: the servos accept at most 732 Hz (JP, 1 Nov). On the 2 kHz loop's 0.5 ms grid the fastest frame at or
     under that rate is 1.5 ms (667 Hz). Allowed frames: 1.5, 2, 2.5, 3 ms.
  3. Servo lead, per column on or off: the command adds servo_tau x planned velocity (servo_tau = 25 ms, the authority
     fact the engine already takes; nothing new is fitted). The engine's servo model and plan model see the same
     lead, so its residual stays exogenous (exploratory test ff_test.py, 1 Nov: tip error down 12 to 24 % at 333 Hz).
The window law, its identification, its floor rule and the e-geodesic carrier (1/tau stepped toward the raw law at
the engine's own bounded speed, frame = slowest allowed frame whose floor stays at or below tau_e) are unchanged.

Arms, all in one pass as columns sharing every encoder read (baton arm RJ1-RJ6, 2 kHz loop, excerpts A and B, steady
and x10 noise step from 10 to 20 s, seeds 1 to 5: 20 paired runs):
  F3, F2, F15   fixed frames 3, 2, 1.5 ms, no lead          F3L, F2L, F15L   the same with the lead
  EG            e-geodesic cadence over {1.5, 2, 2.5, 3}    EGL              the same with the lead
Gates, before the runs:
  G1  this engine with the lead off, 1 kHz, 12 joints, fixed 3 ms, reproduces examples/servo_bench.py (0.2.8) to
      1e-9 mm on A steady seed 5 and B step seed 7.
  G2  with the lead on, the same setting reproduces ff_test.py (1 Nov) on A steady seed 5: 0.599 mm as printed (to
      5e-4 mm).
Efficiency, per family (no lead / lead): frontier(rate) = pooled tip error of that family's fixed arms, linear in the
frame rate between 333, 500 and 667 Hz; gain(EG) = 1 - pooled(EG) / frontier at EG's mean rate.

Predictions (judged; pooled = rms over the 20 runs of the per-run tip rms):
L1  The lead helps the fastest fixed frame: pooled F15L <= 0.92 x pooled F15, and F15L below F15 in >= 18 of 20 runs.
L2  The lead helps the e-geodesic cadence as much: pooled EGL <= 0.92 x pooled EG, and EGL below EG in >= 18 of 20.
E2  The e-geodesic's efficiency holds under 0.2.8 at this ceiling: gain(EG) >= 5 % (31 Oct: 6.7 %).
E2L The same with the lead: gain(EGL) >= 5 %.
A1  Form, not accuracy (31 Oct reading): the fastest fixed frame stays at least as accurate, pooled F15L <= pooled EGL.
Described, not judged: per-cell values, mean frame rates, switches per joint per second, time on the floor, joint
error.

Run: python eg732_lead.py <repo_root> gate | run <excerpt> <cond> | summarize
"""
import os, sys, json, time, hashlib
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.signal import butter, sosfiltfilt

ROOT = os.path.abspath(sys.argv[1]); HERE = os.path.join(ROOT, 'bench', 'cadence_2026-10-31')
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(ROOT, 'examples')); sys.path.insert(0, os.path.join(ROOT, 'python'))
from cadence_bench import SEEDS, Z, TAU_S, W0, ENC_SD, ENC_LSB, TS, DEG, tip
from syntonic_engine import SyntonicEngine
from syntonic_engine.engine import kappa_H, window_floor

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results'); os.makedirs(OUT, exist_ok=True)
SHA = hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest()
FR = [1.5, 2.0, 2.5, 3.0]
ARMS = {'F3': ([3.0], 0, False), 'F2': ([2.0], 0, False), 'F15': ([1.5], 0, False), 'EG': (FR, 2, False),
        'F3L': ([3.0], 0, True), 'F2L': ([2.0], 0, True), 'F15L': ([1.5], 0, True), 'EGL': (FR, 2, True)}


class EGLEngine(SyntonicEngine):
    """egeodesic_cadence.EGEngine with the frame grid FR and a per-column servo lead. source: 0 fixed, 2 e-geodesic."""

    def __init__(self, q0, dt, w0, allowed, source, lead):
        super().__init__(q0, dt, 3e-3, w0, identification='raw')                 # limit: 'authority' (0.2.8)
        n = self.n
        self.Fms = np.array(FR)
        self.F_ticks = np.array([max(1, int(round(f * 1e-3 / dt))) for f in FR])
        self.F_floor = np.array([window_floor(dt, f * 1e-3) for f in FR])
        self.allowed = np.asarray(allowed, bool); self.source = np.asarray(source, int); self.lead_on = np.asarray(lead, bool)
        i3 = FR.index(3.0)
        self.fi = np.array([i3 if (self.source[c] > 0 or self.allowed[c, i3]) else int(np.flatnonzero(self.allowed[c])[0])
                            for c in range(n)])
        self.floor = self.F_floor[self.fi].copy()
        self.next_send = np.zeros(n, dtype=int)
        self.law = np.full(n, np.nan); self.tau_e = np.full(n, self.tau_max)
        self.sent = np.zeros(n); self.switches = np.zeros(n)

    def step(self, q_plan, qd_plan, read_a, read_b):
        dt = self.dt
        q = np.asarray(q_plan, dtype=float)
        a = np.asarray(read_a, dtype=float); b = np.asarray(read_b, dtype=float)
        meas = 0.5 * (a + b); e = q - meas
        d = (a - b) / np.sqrt(2.0)
        self.n2 += self.a_w * (d * d - self.n2)
        r = (q - self.th_ff) + (meas - self.th_m)
        fit = self._identify_raw(r)
        if fit is not None:
            Hn, lam, valid = fit
            s2 = 0.5 * self.n2 * dt
            ts = kappa_H(Hn) * (np.maximum(s2, 1e-30) / np.maximum(lam, 1e-30)) ** (1.0 / (2 * Hn + 1))
            raw = np.minimum(ts, self.tau_max)
            self.law = np.where(valid, raw, self.law)
            re = 1.0 / self.tau_e; dre = self.tau_rate * dt / self.tau_e ** 2          # e-geodesic: rate-linear step
            re = re + np.clip(1.0 / raw - re, -dre, dre)
            self.tau_e = np.where(valid, 1.0 / re, self.tau_e)
            ts = np.clip(ts, self.floor, self.tau_max)
            self.tau_target = np.where(valid, ts, self.tau_target); self.H = np.where(valid, Hn, self.H)
            rr = 1.0 / self.tau; dr_max = self.tau_rate * dt / self.tau ** 2
            rr = rr + np.clip(1.0 / self.tau_target - rr, -dr_max, dr_max); self.tau = 1.0 / rr
        g = 1.0 - np.exp(-dt / self.tau)
        self.u = np.clip(self.u + g * e, -self.limit, self.limit)
        lead = np.where(self.lead_on, self.servo_tau * np.asarray(qd_plan, float), 0.0)
        send = self.k >= self.next_send
        if send.any():
            self.cmd = np.where(send, q + lead + self.u, self.cmd)
            self.sent += send
            src = self.tau_e
            pick = send & (self.source == 2) & ~np.isnan(self.law)
            if pick.any():
                ok = self.allowed & (self.F_floor[None, :] <= src[:, None])
                slowest = np.where(ok.any(1), len(FR) - 1 - np.argmax(ok[:, ::-1], axis=1), np.argmax(self.allowed, axis=1))
                new = np.where(pick, slowest, self.fi)
                self.switches += new != self.fi
                self.fi = new; self.floor = self.F_floor[self.fi]
            self.next_send = np.where(send, self.k + self.F_ticks[self.fi], self.next_send)
        self.th_m = self.th_m + dt * np.clip((self.cmd - self.th_m) / self.servo_tau, -self.w0, self.w0)
        self.th_ff = self.th_ff + dt * np.clip((q + lead - self.th_ff) / self.servo_tau, -self.w0, self.w0)
        self.k += 1
        return self.cmd


def simulate(name, cond, seed, dt, joints, arms, t_end=30.0):
    t0 = float(Z[f'{name}_t0']); tt = np.arange(t0 - 3.0, t0 + t_end, dt)
    J = len(joints); A = len(arms)
    S = CubicSpline(Z[f'{name}_t30'], Z[f'{name}_q30'][:, joints], axis=0)
    qJ = S(tt); qdJ = S(tt, 1)
    tauJ = CubicSpline(Z[f'{name}_t30'], Z[f'{name}_tau30'][:, joints], axis=0)(tt)
    sd = np.full(len(tt), ENC_SD)
    if cond == 'step': sd[(tt >= t0 + 10) & (tt < t0 + 20)] *= np.sqrt(10.0)
    rng = np.random.default_rng(seed)
    tauS = np.tile(TAU_S[joints], A); w0 = np.tile(W0[joints], A)
    k_inv = (1.0 * DEG) / tauS; hb = 0.5 * 0.3 * DEG
    allowed = np.zeros((A * J, len(FR)), bool); source = np.zeros(A * J, int); lead = np.zeros(A * J, bool)
    for i, lab in enumerate(arms):
        fr, src, ld = ARMS[lab]
        for f in fr: allowed[i * J:(i + 1) * J, FR.index(f)] = True
        source[i * J:(i + 1) * J] = src; lead[i * J:(i + 1) * J] = ld
    eng = EGLEngine(np.tile(qJ[0], A), dt, w0, allowed, source, lead)
    th_g = np.tile(qJ[0], A); every = int(round(0.01 / dt)); keep = []
    err2 = np.zeros(A * J); m = 0; below = np.zeros(A * J); sent0 = sw0 = None
    for k in range(len(tt)):
        load = np.tile(tauJ[k], A); q = np.tile(qJ[k], A); qd = np.tile(qdJ[k], A)
        th_out = th_g - load * k_inv - hb * np.tanh(load / (0.02 * tauS))
        na = np.tile(rng.normal(0.0, 1.0, J), A); nb = np.tile(rng.normal(0.0, 1.0, J), A)
        ra = np.round((th_out + na * sd[k]) / ENC_LSB) * ENC_LSB
        rb = np.round((th_out + nb * sd[k]) / ENC_LSB) * ENC_LSB
        cmd = eng.step(q, qd, ra, rb)
        w_eff = w0 * np.maximum(0.05, 1.0 - np.abs(load) / tauS)
        th_g = th_g + dt * np.clip((cmd - th_g) / TS, -w_eff, w_eff)
        if tt[k] >= t0 - 1e-9:
            if sent0 is None: sent0, sw0 = eng.sent.copy(), eng.switches.copy()
            err2 += (q - th_out) ** 2; m += 1; below += eng.law < eng.floor
            if k % every == 0: keep.append((th_out.copy(), q.copy()))
    T = len(tt) * dt - 3.0; out = {}
    for i, lab in enumerate(arms):
        s = slice(i * J, (i + 1) * J)
        e = np.array([tip(a[s][:6]) - tip(b[s][:6]) for a, b in keep])
        hp = sosfiltfilt(butter(4, 5.0, 'highpass', fs=100.0, output='sos'), e, axis=0)
        out[lab] = dict(tip_mm=float(1e3 * np.sqrt(np.mean(np.sum(e ** 2, 1)))),
                        jitter_mm=float(1e3 * np.sqrt(np.mean(np.sum(hp ** 2, 1)))),
                        joint_deg=float(np.mean(np.sqrt(err2[s][:6] / m)) / DEG),
                        rate_hz=float(np.mean(eng.sent[s] - sent0[s]) / T),
                        switches_per_s=float(np.mean(eng.switches[s] - sw0[s]) / T),
                        below_floor=float(np.mean(below[s] / m)))
    return out


if __name__ == '__main__':
    mode = sys.argv[2]
    if mode == 'smoke':
        t = time.time(); r = simulate('A', 'steady', 1, 5e-4, list(range(6)), list(ARMS), t_end=1.5)
        print('smoke ok', sorted(r), round(time.time() - t, 1), 's'); sys.exit(0)
    if mode == 'gate':
        import servo_bench as SB
        g = {}
        ref1 = SB.run('syntonic', 'A', 'steady noise', 5)[0]; ref2 = SB.run('syntonic', 'B', 'noise step', 7)[0]
        m1 = simulate('A', 'steady', 5, 1e-3, list(range(12)), ['F3', 'F3L'])
        m2 = simulate('B', 'step', 7, 1e-3, list(range(12)), ['F3'])
        g['G1'] = dict(A_steady_5=[ref1, m1['F3']['tip_mm']], B_step_7=[ref2, m2['F3']['tip_mm']],
                       ok=bool(abs(ref1 - m1['F3']['tip_mm']) <= 1e-9 and abs(ref2 - m2['F3']['tip_mm']) <= 1e-9))
        g['G2'] = dict(ff_test_printed=0.599, this=m1['F3L']['tip_mm'], ok=bool(abs(m1['F3L']['tip_mm'] - 0.599) <= 5e-4))
        g['script_sha256'] = SHA
        json.dump(g, open(os.path.join(OUT, 'gates.json'), 'w'), indent=1); print(json.dumps(g)); sys.exit(0)
    if mode == 'run':
        name, cond = sys.argv[3], sys.argv[4]
        for seed in SEEDS:
            t = time.time(); r = simulate(name, cond, seed, 5e-4, list(range(6)), list(ARMS))
            with open(os.path.join(OUT, f'runs_{name}_{cond}.jsonl'), 'a') as f:
                f.write(json.dumps(dict(key=f'{name} {cond} {seed}', sha=SHA, **r)) + '\n')
            print(name, cond, seed, {a: round(r[a]['tip_mm'], 3) for a in ARMS}, round(time.time() - t), 's', flush=True)
        sys.exit(0)
    if mode == 'summarize':
        runs = {}
        for fn in sorted(os.listdir(OUT)):
            if fn.startswith('runs_'):
                for line in open(os.path.join(OUT, fn)):
                    d = json.loads(line); runs[d.pop('key')] = d
        assert len(runs) == 20, len(runs)
        G = json.load(open(os.path.join(OUT, 'gates.json'))); gates_ok = G['G1']['ok'] and G['G2']['ok']
        P = {a: float(np.sqrt(np.mean([runs[k][a]['tip_mm'] ** 2 for k in runs]))) for a in ARMS}
        R = {a: float(np.mean([runs[k][a]['rate_hz'] for k in runs])) for a in ARMS}
        def frontier(rate, fam):
            return float(np.interp(rate, [R[f'F3{fam}'], R[f'F2{fam}'], R[f'F15{fam}']], [P[f'F3{fam}'], P[f'F2{fam}'], P[f'F15{fam}']]))
        gain = {'EG': 1 - P['EG'] / frontier(R['EG'], ''), 'EGL': 1 - P['EGL'] / frontier(R['EGL'], 'L')}
        win = lambda a, b: int(sum(runs[k][a]['tip_mm'] < runs[k][b]['tip_mm'] for k in runs))
        V = dict(L1=bool(gates_ok and P['F15L'] <= 0.92 * P['F15'] and win('F15L', 'F15') >= 18),
                 L2=bool(gates_ok and P['EGL'] <= 0.92 * P['EG'] and win('EGL', 'EG') >= 18),
                 E2=bool(gates_ok and gain['EG'] >= 0.05), E2L=bool(gates_ok and gain['EGL'] >= 0.05),
                 A1=bool(gates_ok and P['F15L'] <= P['EGL']))
        cells = {a: {c: float(np.sqrt(np.mean([runs[k][a]['tip_mm'] ** 2 for k in runs if k.startswith(c)])))
                     for c in ('A steady', 'A step', 'B steady', 'B step')} for a in ARMS}
        desc = {a: dict(rate_hz=R[a], switches_per_s=float(np.mean([runs[k][a]['switches_per_s'] for k in runs])),
                        below_floor=float(np.mean([runs[k][a]['below_floor'] for k in runs])),
                        joint_deg=float(np.mean([runs[k][a]['joint_deg'] for k in runs]))) for a in ARMS}
        S = dict(script_sha256=SHA, gates=G, pooled_tip_mm=P, cells=cells, gain=gain,
                 wins=dict(F15L_vs_F15=win('F15L', 'F15'), EGL_vs_EG=win('EGL', 'EG'), F15L_vs_EGL=win('F15L', 'EGL')),
                 ratios=dict(F15L_F15=P['F15L'] / P['F15'], EGL_EG=P['EGL'] / P['EG'], EGL_F15L=P['EGL'] / P['F15L']),
                 verdict=V, described=desc)
        json.dump(S, open(os.path.join(OUT, 'summary.json'), 'w'), indent=1); print(json.dumps(S, indent=1))
