# Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
# Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
"""Where the e-geodesic can make a difference: how the window travels between two targets. Written 1 Nov 2026, before
any run of this file except the smoke test; request from JP ("trouve moi la place de l'e-geodesique ici et qui fait une
difference").

Map so far: as a cadence path (31 Oct, E1-E3) the e-geodesic kept the law's timing but was accuracy-neutral (EG within
1 % of the raw law, its speed bound almost never binding); with the servo lead the cadence falls to the slowest frame
(eg732_lead, eg732_hold). What remains is the window's own travel. The adaptive window just beat the best fixed window
by 4.5 % held out, most when the noise jumps (window_vs_frozen). Dual Geodesy (Theorem 1; Register V1.3, repaired P2):
the canonical fiber is an e-geodesic, and the open question is the time parameterization, which coordinate relaxes
exponentially, theta = 1/tau (replicator, e) or tau (Newton). This test puts that choice in the engine and scores it.

Arms (F15D setting of window_vs_frozen.py: engine 0.2.8, 2 kHz loop, fixed 1.5 ms frames, lead (servo_tau + D) x
planned velocity; the law's target tau_target unchanged in every arm):
  BASE   the published transport: 1/tau stepped toward 1/tau_target, clipped at tau_rate (0.4 s/s); rarely binds.
  E1     theta = 1/tau relaxes exponentially toward 1/tau_target with time constant 1 s (the engine's ema).
  N1     tau relaxes exponentially toward tau_target with the same 1 s.
  E01    as E1 with 0.1 s.        N01    as N1 with 0.1 s.
Expected shape, before any run: from the same start and target, theta-relaxation shortens the window faster and
lengthens it slower than tau-relaxation (at t = T/10 from 10.8 to 5 ms: 9.7 against 10.25 ms; from 5 to 10.8 ms: 5.27
against 5.55 ms). On this bench lag costs more than noise, so a transport quick to shorten should track better.
Start-up: E and N arms take the law's target as their window at the first fit (after the 1 s warm-up), then
relax; a 0.5 s start relaxed by tau at 1 s would still be at about 70 ms when scoring begins, a start-up artifact.
Gate: G6, BASE reproduces eg732_hold.py's F15D run by run to 1e-9 mm.

Predictions (judged; 20 paired runs: excerpts A and B, steady and x10 noise step, seeds 1 to 5):
P2a  The e coordinate transports better at 1 s: pooled E1 <= pooled N1, and E1 below N1 in >= 15 of 20 runs.
P2b  The same at 0.1 s: pooled E01 <= pooled N01, and E01 below N01 in >= 15 of 20 runs.
Null clause, described: if |pooled E / pooled N - 1| < 1 % at both time constants, the coordinate makes no measurable
difference on this bench, whatever the signs.
Described, not judged: the transient windows 10 to 12 s and 20 to 22 s of the noise-step runs; BASE against E and N.

Run: python window_transport.py <repo_root> smoke | run <excerpt> <cond> | summarize <eg732_hold results>
"""
import os, sys, json, time, hashlib
import numpy as np
from scipy.interpolate import CubicSpline

ROOT = os.path.abspath(sys.argv[1]); MODE = sys.argv[2]
HERE_W = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE_W)
_argv = sys.argv; sys.argv = [_argv[0], ROOT, 'import']
import eg732_hold as H
from syntonic_engine.engine import kappa_H
sys.argv = _argv
OUT = os.path.join(HERE_W, 'results_transport'); os.makedirs(OUT, exist_ok=True)
SHA = hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest()
ARMS = {'BASE': (0, 0.0), 'E1': (1, 1.0), 'N1': (2, 1.0), 'E01': (1, 0.1), 'N01': (2, 0.1)}   # transport, T (s)


class TransportEngine(H.HoldEngine):
    """HoldEngine (fixed frames, lead mode 2) with a per-column window transport: 0 published, 1 theta, 2 tau."""

    def __init__(self, *a, transport, T, **k):
        super().__init__(*a, **k)
        self.tr = np.asarray(transport, int); self.T = np.asarray(T, float)

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
            self.law = np.where(valid, np.minimum(ts, self.tau_max), self.law)
            ts = np.clip(ts, self.floor, self.tau_max)
            self.tau_target = np.where(valid, ts, self.tau_target); self.H = np.where(valid, Hn, self.H)
            rr = 1.0 / self.tau; dr_max = self.tau_rate * dt / self.tau ** 2
            base = 1.0 / (rr + np.clip(1.0 / self.tau_target - rr, -dr_max, dr_max))
            a_r = np.where(self.T > 0, dt / np.maximum(self.T, 1e-12), 0.0)
            th = rr + a_r * (1.0 / self.tau_target - rr)
            e_tau = 1.0 / th
            n_tau = self.tau + a_r * (self.tau_target - self.tau)
            if not getattr(self, '_started', False):              # E and N start on the target, not at 0.5 s
                e_tau = n_tau = self.tau_target.copy(); self._started = True
            self.tau = np.where(self.tr == 0, base, np.where(self.tr == 1, e_tau, n_tau))
        g = 1.0 - np.exp(-dt / self.tau)
        self.u = np.clip(self.u + g * e, -self.limit, self.limit)
        send = self.k >= self.next_send
        if send.any():
            lead = self._coef() * qd
            self.cmd = np.where(send, q + lead + self.u, self.cmd)
            self.sent += send
            self.next_send = np.where(send, self.k + self.F_ticks[self.fi], self.next_send)
        lead = self._coef() * qd
        self.th_m = self.th_m + dt * np.clip((self.cmd - self.th_m) / self.servo_tau, -self.w0, self.w0)
        self.th_ff = self.th_ff + dt * np.clip((q + lead - self.th_ff) / self.servo_tau, -self.w0, self.w0)
        self.k += 1
        return self.cmd


def simulate(name, cond, seed, arms, t_end=30.0):
    dt = 5e-4; joints = list(range(6)); J = 6; A = len(arms)
    t0 = float(H.Z[f'{name}_t0']); tt = np.arange(t0 - 3.0, t0 + t_end, dt)
    S = CubicSpline(H.Z[f'{name}_t30'], H.Z[f'{name}_q30'][:, joints], axis=0); qJ = S(tt); qdJ = S(tt, 1)
    tauJ = CubicSpline(H.Z[f'{name}_t30'], H.Z[f'{name}_tau30'][:, joints], axis=0)(tt)
    sd = np.full(len(tt), H.ENC_SD)
    if cond == 'step': sd[(tt >= t0 + 10) & (tt < t0 + 20)] *= np.sqrt(10.0)
    rng = np.random.default_rng(seed)
    tauS = np.tile(H.TAU_S[joints], A); w0 = np.tile(H.W0[joints], A)
    k_inv = (1.0 * H.DEG) / tauS; hb = 0.5 * 0.3 * H.DEG
    allowed = np.zeros((A * J, len(H.FR)), bool); allowed[:, H.FR.index(1.5)] = True
    tr = np.repeat([ARMS[a][0] for a in arms], J); T = np.repeat([ARMS[a][1] for a in arms], J)
    eng = TransportEngine(np.tile(qJ[0], A), dt, w0, allowed, np.zeros(A * J, int), np.full(A * J, 2), transport=tr, T=T)
    th_g = np.tile(qJ[0], A); every = int(round(0.01 / dt)); keep = []; ts = []
    for k in range(len(tt)):
        load = np.tile(tauJ[k], A); q = np.tile(qJ[k], A); qd = np.tile(qdJ[k], A)
        th_out = th_g - load * k_inv - hb * np.tanh(load / (0.02 * tauS))
        na = np.tile(rng.normal(0.0, 1.0, J), A); nb = np.tile(rng.normal(0.0, 1.0, J), A)
        ra = np.round((th_out + na * sd[k]) / H.ENC_LSB) * H.ENC_LSB
        rb = np.round((th_out + nb * sd[k]) / H.ENC_LSB) * H.ENC_LSB
        cmd = eng.step(q, qd, ra, rb)
        w_eff = w0 * np.maximum(0.05, 1.0 - np.abs(load) / tauS)
        th_g = th_g + dt * np.clip((cmd - th_g) / H.TS, -w_eff, w_eff)
        if tt[k] >= t0 - 1e-9 and k % every == 0: keep.append((th_out.copy(), q.copy())); ts.append(tt[k] - t0)
    ts = np.array(ts); out = {}
    for i, lab in enumerate(arms):
        s = slice(i * J, (i + 1) * J)
        e2 = np.array([np.sum((H.tip(a[s]) - H.tip(b[s])) ** 2) for a, b in keep])
        tr_m = ((ts >= 10) & (ts < 12)) | ((ts >= 20) & (ts < 22))
        out[lab] = dict(tip_mm=float(1e3 * np.sqrt(np.mean(e2))), transient_mm=float(1e3 * np.sqrt(np.mean(e2[tr_m]))))
    return out


if __name__ == '__main__':
    arms = list(ARMS)
    if MODE == 'smoke':
        t = time.time(); r = simulate('A', 'step', 1, arms, t_end=1.5); print('smoke ok', sorted(r), round(time.time() - t, 1), 's'); sys.exit(0)
    if MODE == 'run':
        name, cond = sys.argv[3], sys.argv[4]
        for seed in H.SEEDS:
            t = time.time(); r = simulate(name, cond, seed, arms)
            with open(os.path.join(OUT, f'runs_{name}_{cond}.jsonl'), 'a') as f: f.write(json.dumps(dict(key=f'{name} {cond} {seed}', sha=SHA, **r)) + '\n')
            print(name, cond, seed, {a: round(r[a]['tip_mm'], 3) for a in arms}, round(time.time() - t), 's', flush=True)
        sys.exit(0)
    if MODE == 'summarize':
        def load(d):
            runs = {}
            for fn in sorted(os.listdir(d)):
                if fn.startswith('runs_'):
                    for line in open(os.path.join(d, fn)):
                        x = json.loads(line); runs[x.pop('key')] = x
            return runs
        runs = load(OUT); hold = load(sys.argv[3]); assert len(runs) == 20
        g6 = max(abs(runs[k]['BASE']['tip_mm'] - hold[k]['F15D']['tip_mm']) for k in runs)
        P = {a: float(np.sqrt(np.mean([runs[k][a]['tip_mm'] ** 2 for k in runs]))) for a in arms}
        steps = [k for k in runs if ' step ' in k]
        PT = {a: float(np.sqrt(np.mean([runs[k][a]['transient_mm'] ** 2 for k in steps]))) for a in arms}
        win = lambda x, y: int(sum(runs[k][x]['tip_mm'] < runs[k][y]['tip_mm'] for k in runs))
        ok = g6 <= 1e-9
        V = dict(P2a=bool(ok and P['E1'] <= P['N1'] and win('E1', 'N1') >= 15),
                 P2b=bool(ok and P['E01'] <= P['N01'] and win('E01', 'N01') >= 15),
                 null_clause=bool(abs(P['E1'] / P['N1'] - 1) < 0.01 and abs(P['E01'] / P['N01'] - 1) < 0.01))
        S = dict(script_sha256=SHA, G6=dict(max_diff_mm=g6, ok=ok), pooled_tip_mm=P, transient_step_mm=PT, verdict=V,
                 ratios=dict(E1_N1=P['E1'] / P['N1'], E01_N01=P['E01'] / P['N01'], E1_BASE=P['E1'] / P['BASE'], E01_BASE=P['E01'] / P['BASE']),
                 wins=dict(E1_vs_N1=win('E1', 'N1'), E01_vs_N01=win('E01', 'N01'), E01_vs_BASE=win('E01', 'BASE'), E1_vs_BASE=win('E1', 'BASE')),
                 cells={c: {a: float(np.sqrt(np.mean([runs[k][a]['tip_mm'] ** 2 for k in runs if k.startswith(c)]))) for a in arms}
                        for c in ('A steady', 'A step', 'B steady', 'B step')})
        json.dump(S, open(os.path.join(OUT, 'summary.json'), 'w'), indent=1); print(json.dumps(S, indent=1))
