# Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
# Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
"""The servo lead with the loop delay: does leading the frame's hold close the frame-rate gap, and does the law's
cadence then regain a use? Written 1 Nov 2026 after eg732_lead.py (sha256 507791a9...), before any run of this file
except the smoke test; go from JP ("Ok").

eg732_lead.py found: with the lead servo_tau x planned velocity, the fastest fixed frame (1.5 ms) is 11 % more accurate
than the law's e-geodesic cadence, which falls to the slowest frame because its window no longer sees a need for fast
frames. Reading to test: what fast frames still buy is the hold of the frame (a command held for a frame is on
average half a frame old), which that lead does not cover.

Change, per column: lead mode 1 = servo_tau x planned velocity (as eg732_lead.py); lead mode 2 = (servo_tau + D) x
planned velocity, with D = loop_delay(dt, frame) = frame/2 + dt/2 + sensor delay (0.1 ms), the delay the engine already
declares and uses for its floor, for the frame being held after each send (for the e-geodesic arm, the frame it has
just chosen). No new constant. The engine's servo and plan models see the same lead. Everything else as eg732_lead.py:
engine 0.2.8, 2 kHz loop, baton arm, frames 1.5, 2, 2.5, 3 ms (1.5 ms = fastest grid frame under the servos' 732 Hz),
excerpts A and B, steady and x10 noise step, seeds 1 to 5: 20 paired runs.

Arms: F3L, F2L, F15L, EGL (mode 1, the reference, identical to eg732_lead.py), F3D, F2D, F15D, EGD (mode 2).
Gates: G1, lead off, 1 kHz, 12 joints, fixed 3 ms, reproduces examples/servo_bench.py (0.2.8) on A steady seed 5 and
B step seed 7 to 1e-9 mm. G3, the four mode-1 arms reproduce eg732_lead.py's results run by run to 1e-9 mm.

Predictions (judged; pooled = rms over the 20 runs of the per-run tip rms):
D1  The frame gap closes: pooled F3D <= 1.05 x pooled F15D (with mode 1 it was 1.134).
D2  The hold lead helps the slowest frame: pooled F3D <= 0.92 x pooled F3L, and F3D below F3L in >= 18 of 20 runs.
D3  It costs nothing at the fastest frame: pooled F15D <= pooled F15L.
D4  The law's cadence regains a use: pooled EGD <= 1.03 x pooled F15D, at a mean frame rate <= 500 Hz.
Described, not judged: per-cell values, mean frame rates, the law below the floor, joint error, A steady seed 2.

Run: python eg732_hold.py <repo_root> smoke | gate | run <excerpt> <cond> | summarize <eg732_lead results dir>
"""
import os, sys, json, time, hashlib
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.signal import butter, sosfiltfilt

ROOT = os.path.abspath(sys.argv[1]); HERE = os.path.join(ROOT, 'bench', 'cadence_2026-10-31')
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(ROOT, 'examples')); sys.path.insert(0, os.path.join(ROOT, 'python'))
from cadence_bench import SEEDS, Z, TAU_S, W0, ENC_SD, ENC_LSB, TS, DEG, tip
from syntonic_engine import SyntonicEngine
from syntonic_engine.engine import kappa_H, window_floor, loop_delay

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results_hold'); os.makedirs(OUT, exist_ok=True)
SHA = hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest()
FR = [1.5, 2.0, 2.5, 3.0]
ARMS = {'F3L': ([3.0], 0, 1), 'F2L': ([2.0], 0, 1), 'F15L': ([1.5], 0, 1), 'EGL': (FR, 2, 1),
        'F3D': ([3.0], 0, 2), 'F2D': ([2.0], 0, 2), 'F15D': ([1.5], 0, 2), 'EGD': (FR, 2, 2)}
GATE_ARMS = {'F3': ([3.0], 0, 0)}
MAIN = list(ARMS)


class HoldEngine(SyntonicEngine):
    """eg732_lead.EGLEngine with a lead mode per column: 0 none, 1 servo_tau, 2 servo_tau + loop delay of the held frame.
    At a send the frame is chosen first, then the command is formed (modes 0 and 1 are unaffected by the order)."""

    def __init__(self, q0, dt, w0, allowed, source, mode):
        super().__init__(q0, dt, 3e-3, w0, identification='raw')                 # limit: 'authority' (0.2.8)
        n = self.n
        self.Fms = np.array(FR)
        self.F_ticks = np.array([max(1, int(round(f * 1e-3 / dt))) for f in FR])
        self.F_floor = np.array([window_floor(dt, f * 1e-3) for f in FR])
        self.F_D = np.array([loop_delay(dt, f * 1e-3) for f in FR])
        self.allowed = np.asarray(allowed, bool); self.source = np.asarray(source, int); self.mode = np.asarray(mode, int)
        i3 = FR.index(3.0)
        self.fi = np.array([i3 if (self.source[c] > 0 or self.allowed[c, i3]) else int(np.flatnonzero(self.allowed[c])[0])
                            for c in range(n)])
        self.floor = self.F_floor[self.fi].copy()
        self.next_send = np.zeros(n, dtype=int)
        self.law = np.full(n, np.nan); self.tau_e = np.full(n, self.tau_max)
        self.sent = np.zeros(n); self.switches = np.zeros(n)

    def _coef(self):
        return np.where(self.mode == 0, 0.0, np.where(self.mode == 1, self.servo_tau, self.servo_tau + self.F_D[self.fi]))

    def step(self, q_plan, qd_plan, read_a, read_b):
        dt = self.dt
        q = np.asarray(q_plan, dtype=float); qd = np.asarray(qd_plan, dtype=float)
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
            re = 1.0 / self.tau_e; dre = self.tau_rate * dt / self.tau_e ** 2
            re = re + np.clip(1.0 / raw - re, -dre, dre)
            self.tau_e = np.where(valid, 1.0 / re, self.tau_e)
            ts = np.clip(ts, self.floor, self.tau_max)
            self.tau_target = np.where(valid, ts, self.tau_target); self.H = np.where(valid, Hn, self.H)
            rr = 1.0 / self.tau; dr_max = self.tau_rate * dt / self.tau ** 2
            rr = rr + np.clip(1.0 / self.tau_target - rr, -dr_max, dr_max); self.tau = 1.0 / rr
        g = 1.0 - np.exp(-dt / self.tau)
        self.u = np.clip(self.u + g * e, -self.limit, self.limit)
        send = self.k >= self.next_send
        if send.any():
            pick = send & (self.source == 2) & ~np.isnan(self.law)
            if pick.any():
                ok = self.allowed & (self.F_floor[None, :] <= self.tau_e[:, None])
                slowest = np.where(ok.any(1), len(FR) - 1 - np.argmax(ok[:, ::-1], axis=1), np.argmax(self.allowed, axis=1))
                new = np.where(pick, slowest, self.fi)
                self.switches += new != self.fi
                self.fi = new; self.floor = self.F_floor[self.fi]
            lead = self._coef() * qd
            self.cmd = np.where(send, q + lead + self.u, self.cmd)
            self.sent += send
            self.next_send = np.where(send, self.k + self.F_ticks[self.fi], self.next_send)
        lead = self._coef() * qd
        self.th_m = self.th_m + dt * np.clip((self.cmd - self.th_m) / self.servo_tau, -self.w0, self.w0)
        self.th_ff = self.th_ff + dt * np.clip((q + lead - self.th_ff) / self.servo_tau, -self.w0, self.w0)
        self.k += 1
        return self.cmd


def simulate(name, cond, seed, dt, joints, arms, table, t_end=30.0):
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
    allowed = np.zeros((A * J, len(FR)), bool); source = np.zeros(A * J, int); mode = np.zeros(A * J, int)
    for i, lab in enumerate(arms):
        fr, src, md = table[lab]
        for f in fr: allowed[i * J:(i + 1) * J, FR.index(f)] = True
        source[i * J:(i + 1) * J] = src; mode[i * J:(i + 1) * J] = md
    eng = HoldEngine(np.tile(qJ[0], A), dt, w0, allowed, source, mode)
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
        t = time.time(); r = simulate('A', 'steady', 1, 5e-4, list(range(6)), MAIN, ARMS, t_end=1.5)
        print('smoke ok', sorted(r), round(time.time() - t, 1), 's'); sys.exit(0)
    if mode == 'gate':
        import servo_bench as SB
        ref1 = SB.run('syntonic', 'A', 'steady noise', 5)[0]; ref2 = SB.run('syntonic', 'B', 'noise step', 7)[0]
        m1 = simulate('A', 'steady', 5, 1e-3, list(range(12)), ['F3'], GATE_ARMS)['F3']['tip_mm']
        m2 = simulate('B', 'step', 7, 1e-3, list(range(12)), ['F3'], GATE_ARMS)['F3']['tip_mm']
        g = dict(G1=dict(A_steady_5=[ref1, m1], B_step_7=[ref2, m2], ok=bool(abs(ref1 - m1) <= 1e-9 and abs(ref2 - m2) <= 1e-9)),
                 script_sha256=SHA)
        json.dump(g, open(os.path.join(OUT, 'gates.json'), 'w'), indent=1); print(json.dumps(g)); sys.exit(0)
    if mode == 'run':
        name, cond = sys.argv[3], sys.argv[4]
        for seed in SEEDS:
            t = time.time(); r = simulate(name, cond, seed, 5e-4, list(range(6)), MAIN, ARMS)
            with open(os.path.join(OUT, f'runs_{name}_{cond}.jsonl'), 'a') as f:
                f.write(json.dumps(dict(key=f'{name} {cond} {seed}', sha=SHA, **r)) + '\n')
            print(name, cond, seed, {a: round(r[a]['tip_mm'], 3) for a in MAIN}, round(time.time() - t), 's', flush=True)
        sys.exit(0)
    if mode == 'summarize':
        prev_dir = sys.argv[3]
        def load(d):
            runs = {}
            for fn in sorted(os.listdir(d)):
                if fn.startswith('runs_'):
                    for line in open(os.path.join(d, fn)):
                        x = json.loads(line); runs[x.pop('key')] = x
            return runs
        runs = load(OUT); prev = load(prev_dir)
        assert len(runs) == 20 and len(prev) == 20, (len(runs), len(prev))
        g3 = max(abs(runs[k][a]['tip_mm'] - prev[k][a]['tip_mm']) for k in runs for a in ('F3L', 'F2L', 'F15L', 'EGL'))
        G = json.load(open(os.path.join(OUT, 'gates.json'))); G['G3'] = dict(max_diff_mm=g3, ok=bool(g3 <= 1e-9))
        gates_ok = G['G1']['ok'] and G['G3']['ok']
        P = {a: float(np.sqrt(np.mean([runs[k][a]['tip_mm'] ** 2 for k in runs]))) for a in MAIN}
        R = {a: float(np.mean([runs[k][a]['rate_hz'] for k in runs])) for a in MAIN}
        win = lambda a, b: int(sum(runs[k][a]['tip_mm'] < runs[k][b]['tip_mm'] for k in runs))
        V = dict(D1=bool(gates_ok and P['F3D'] <= 1.05 * P['F15D']),
                 D2=bool(gates_ok and P['F3D'] <= 0.92 * P['F3L'] and win('F3D', 'F3L') >= 18),
                 D3=bool(gates_ok and P['F15D'] <= P['F15L']),
                 D4=bool(gates_ok and P['EGD'] <= 1.03 * P['F15D'] and R['EGD'] <= 500.0))
        cells = {a: {c: float(np.sqrt(np.mean([runs[k][a]['tip_mm'] ** 2 for k in runs if k.startswith(c)])))
                     for c in ('A steady', 'A step', 'B steady', 'B step')} for a in MAIN}
        desc = {a: dict(rate_hz=R[a], below_floor=float(np.mean([runs[k][a]['below_floor'] for k in runs])),
                        joint_deg=float(np.mean([runs[k][a]['joint_deg'] for k in runs])),
                        switches_per_s=float(np.mean([runs[k][a]['switches_per_s'] for k in runs]))) for a in MAIN}
        S = dict(script_sha256=SHA, gates=G, pooled_tip_mm=P, cells=cells, verdict=V,
                 ratios=dict(F3D_F15D=P['F3D'] / P['F15D'], F3L_F15L=P['F3L'] / P['F15L'], F3D_F3L=P['F3D'] / P['F3L'],
                             F15D_F15L=P['F15D'] / P['F15L'], EGD_F15D=P['EGD'] / P['F15D']),
                 wins=dict(F3D_vs_F3L=win('F3D', 'F3L'), F15D_vs_F15L=win('F15D', 'F15L'), EGD_vs_F15D=win('EGD', 'F15D')),
                 a_steady_2={a: runs['A steady 2'][a]['tip_mm'] for a in MAIN}, described=desc)
        json.dump(S, open(os.path.join(OUT, 'summary.json'), 'w'), indent=1); print(json.dumps(S, indent=1))
