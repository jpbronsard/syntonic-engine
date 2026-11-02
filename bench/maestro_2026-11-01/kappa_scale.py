# Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
# Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
"""Is kappa_H placing the window where it should on MAESTRO 1? Written 2 Nov 2026, before any run of this file except the
smoke test; go from JP ("go").

The engine does not fit kappa: kappa_H = (2H Gamma(2H+1))^(-1/(2H+1)) from Theorem U', with H identified live (0.95 to
0.99 on the baton arm, so kappa_H about 0.64). The oracle test of the 27 Oct handoff: sweep the window by hand, let the
law predict its own, compare (within x1.5, regret under 10 %). Here the law's target is multiplied by a factor c before
the floor and ceiling clip (equivalent to replacing kappa_H by c x kappa_H), c in {0.5, 0.7, 1, 1.4, 2}, all as columns
sharing every encoder read.

Settings (baton arm RJ1-RJ6, excerpts A and B, steady and x10 noise step from 10 to 20 s, seeds 31 to 35, never used:
20 paired runs per setting):
  S1  the published bench: 1 kHz loop, fixed 3 ms frames (333 Hz), no lead (eg732_hold.HoldEngine, lead mode 0).
  S2  tonight's best configuration: 2 kHz loop, fixed 1.5 ms frames, lead (servo_tau + D) x planned velocity (mode 2).
Engine 0.2.8 (authority limit), identification 'raw', published transport.
Gates: G9, S1 with c = 1 on 12 joints, A steady seed 5, reproduces examples/servo_bench.py (0.683910 mm) to 1e-9 mm;
G10, S2 with c = 1, A steady seed 1, reproduces eg732_hold.py's F15D value to 1e-9 mm.

Predictions (judged; pooled = rms over the 20 runs of the per-run tip rms; c* = the factor with the lowest pooled):
Q1  S1: c* in {0.7, 1, 1.4} (the x1.5 band on this grid), and regret(c = 1) = pooled(1) / pooled(c*) - 1 < 10 %.
Q2  S2: the same.
Described, not judged: the optimum by a parabola in log c through the three lowest points; pooled by cell; the share
of time the law's window sits on the floor at each c.

Run: python kappa_scale.py <repo_root> smoke | gate | run <S1|S2> <excerpt> <cond> | summarize
"""
import os, sys, json, time, hashlib
import numpy as np
from scipy.interpolate import CubicSpline

ROOT = os.path.abspath(sys.argv[1]); MODE = sys.argv[2]; _args = list(sys.argv)
HERE_W = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE_W)
sys.argv = [_args[0], ROOT, 'import']
import eg732_hold as H
from syntonic_engine.engine import kappa_H
sys.argv = _args
OUT = os.path.join(HERE_W, 'results_kappa'); os.makedirs(OUT, exist_ok=True)
SHA = hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest()
CS = [0.5, 0.7, 1.0, 1.4, 2.0]
SEEDS = [31, 32, 33, 34, 35]
SETTINGS = {'S1': dict(dt=1e-3, frame=3.0, mode=0), 'S2': dict(dt=5e-4, frame=1.5, mode=2)}


class ScaledEngine(H.HoldEngine):
    """HoldEngine with fixed frames; the law's target is multiplied by a per-column factor c before the clip."""

    def __init__(self, *a, c, **k):
        super().__init__(*a, **k)
        self.c = np.asarray(c, float); self.on_floor = np.zeros(self.n)

    def step(self, q_plan, qd_plan, read_a, read_b):
        dt = self.dt
        q = np.asarray(q_plan, dtype=float); qd = np.asarray(qd_plan, dtype=float)
        a = np.asarray(read_a, dtype=float); b = np.asarray(read_b, dtype=float)
        meas = 0.5 * (a + b); e = q - meas
        d = (a - b) / np.sqrt(2.0)
        self.n2 += self.a_w * (d * d - self.n2)
        r = (q - self.th_ff) + (meas - self.th_m)
        fit = self._identify_raw(r)
        self.on_floor = np.zeros(self.n)
        if fit is not None:
            Hn, lam, valid = fit
            s2 = 0.5 * self.n2 * dt
            ts = self.c * kappa_H(Hn) * (np.maximum(s2, 1e-30) / np.maximum(lam, 1e-30)) ** (1.0 / (2 * Hn + 1))
            self.law = np.where(valid, np.minimum(ts, self.tau_max), self.law)
            self.on_floor = (ts <= self.floor).astype(float)
            ts = np.clip(ts, self.floor, self.tau_max)
            self.tau_target = np.where(valid, ts, self.tau_target); self.H = np.where(valid, Hn, self.H)
            rr = 1.0 / self.tau; dr_max = self.tau_rate * dt / self.tau ** 2
            rr = rr + np.clip(1.0 / self.tau_target - rr, -dr_max, dr_max); self.tau = 1.0 / rr
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


def simulate(setting, name, cond, seed, cs, joints=None, t_end=30.0):
    st = SETTINGS[setting]; dt = st['dt']
    joints = list(range(6)) if joints is None else joints; J = len(joints); A = len(cs)
    t0 = float(H.Z[f'{name}_t0']); tt = np.arange(t0 - 3.0, t0 + t_end, dt)
    S = CubicSpline(H.Z[f'{name}_t30'], H.Z[f'{name}_q30'][:, joints], axis=0); qJ = S(tt); qdJ = S(tt, 1)
    tauJ = CubicSpline(H.Z[f'{name}_t30'], H.Z[f'{name}_tau30'][:, joints], axis=0)(tt)
    sd = np.full(len(tt), H.ENC_SD)
    if cond == 'step': sd[(tt >= t0 + 10) & (tt < t0 + 20)] *= np.sqrt(10.0)
    rng = np.random.default_rng(seed)
    tauS = np.tile(H.TAU_S[joints], A); w0 = np.tile(H.W0[joints], A)
    k_inv = (1.0 * H.DEG) / tauS; hb = 0.5 * 0.3 * H.DEG
    allowed = np.zeros((A * J, len(H.FR)), bool); allowed[:, H.FR.index(st['frame'])] = True
    eng = ScaledEngine(np.tile(qJ[0], A), dt, w0, allowed, np.zeros(A * J, int), np.full(A * J, st['mode']), c=np.repeat(cs, J))
    th_g = np.tile(qJ[0], A); every = int(round(0.01 / dt)); keep = []; flo = np.zeros(A * J); m = 0
    for k in range(len(tt)):
        load = np.tile(tauJ[k], A); q = np.tile(qJ[k], A); qd = np.tile(qdJ[k], A)
        th_out = th_g - load * k_inv - hb * np.tanh(load / (0.02 * tauS))
        na = np.tile(rng.normal(0.0, 1.0, J), A); nb = np.tile(rng.normal(0.0, 1.0, J), A)
        ra = np.round((th_out + na * sd[k]) / H.ENC_LSB) * H.ENC_LSB
        rb = np.round((th_out + nb * sd[k]) / H.ENC_LSB) * H.ENC_LSB
        cmd = eng.step(q, qd, ra, rb)
        w_eff = w0 * np.maximum(0.05, 1.0 - np.abs(load) / tauS)
        th_g = th_g + dt * np.clip((cmd - th_g) / H.TS, -w_eff, w_eff)
        if tt[k] >= t0 - 1e-9:
            flo += eng.on_floor; m += 1
            if k % every == 0: keep.append((th_out.copy(), q.copy()))
    out = {}
    for i, c in enumerate(cs):
        s = slice(i * J, (i + 1) * J)
        e = np.array([H.tip(a[s][:6]) - H.tip(b[s][:6]) for a, b in keep])
        out[str(c)] = dict(tip_mm=float(1e3 * np.sqrt(np.mean(np.sum(e ** 2, 1)))), on_floor=float(np.mean(flo[s][:6] / m)))
    return out


if __name__ == '__main__':
    if MODE == 'smoke':
        t = time.time(); r1 = simulate('S1', 'A', 'step', 31, CS, t_end=1.5); r2 = simulate('S2', 'A', 'step', 31, CS, t_end=1.5)
        print('smoke ok', sorted(r1), sorted(r2), round(time.time() - t, 1), 's'); sys.exit(0)
    if MODE == 'gate':
        g9 = simulate('S1', 'A', 'steady', 5, [1.0], joints=list(range(12)))['1.0']['tip_mm']
        ref10 = None
        for line in open(os.path.join(HERE_W, 'results_hold', 'runs_A_steady.jsonl')):
            x = json.loads(line)
            if x['key'] == 'A steady 1': ref10 = x['F15D']['tip_mm']
        g10 = simulate('S2', 'A', 'steady', 1, [1.0])['1.0']['tip_mm']
        G = dict(G9=dict(ref=0.6839103668504339, this=g9, ok=bool(abs(g9 - 0.6839103668504339) <= 1e-9)),
                 G10=dict(ref=ref10, this=g10, ok=bool(abs(g10 - ref10) <= 1e-9)), script_sha256=SHA)
        json.dump(G, open(os.path.join(OUT, 'gates.json'), 'w'), indent=1); print(json.dumps(G)); sys.exit(0)
    if MODE == 'run':
        setting, name, cond = sys.argv[3], sys.argv[4], sys.argv[5]
        for seed in SEEDS:
            t = time.time(); r = simulate(setting, name, cond, seed, CS)
            with open(os.path.join(OUT, 'runs.jsonl'), 'a') as f: f.write(json.dumps(dict(setting=setting, key=f'{name} {cond} {seed}', sha=SHA, **r)) + '\n')
            print(setting, name, cond, seed, {c: round(v['tip_mm'], 3) for c, v in r.items()}, round(time.time() - t), 's', flush=True)
        sys.exit(0)
    if MODE == 'summarize':
        rows = [json.loads(l) for l in open(os.path.join(OUT, 'runs.jsonl'))]
        G = json.load(open(os.path.join(OUT, 'gates.json'))); ok = G['G9']['ok'] and G['G10']['ok']
        S = dict(script_sha256=SHA, gates=G)
        for st in ('S1', 'S2'):
            R = [r for r in rows if r['setting'] == st]; assert len(R) == 20, (st, len(R))
            P = {str(c): float(np.sqrt(np.mean([r[str(c)]['tip_mm'] ** 2 for r in R]))) for c in CS}
            cstar = min(P, key=P.get); regret = P['1.0'] / P[cstar] - 1
            idx = np.argsort([P[str(c)] for c in CS])[:3]; x = np.log([CS[i] for i in idx]); y = np.array([P[str(CS[i])] for i in idx])
            co = np.polyfit(x, y, 2); c_par = float(np.exp(-co[1] / (2 * co[0]))) if co[0] > 0 else None
            cells = {cl: {str(c): float(np.sqrt(np.mean([r[str(c)]['tip_mm'] ** 2 for r in R if r['key'].startswith(cl)]))) for c in CS}
                     for cl in ('A steady', 'A step', 'B steady', 'B step')}
            S[st] = dict(pooled=P, c_star=float(cstar), regret_c1=regret, c_parabola=c_par, cells=cells,
                         on_floor={str(c): float(np.mean([r[str(c)]['on_floor'] for r in R])) for c in CS},
                         wins_c1_vs_cstar=int(sum(r['1.0']['tip_mm'] <= r[cstar]['tip_mm'] for r in R)))
            S[st]['verdict'] = bool(ok and float(cstar) in (0.7, 1.0, 1.4) and regret < 0.10)
        S['verdict'] = dict(Q1=S['S1']['verdict'], Q2=S['S2']['verdict'])
        json.dump(S, open(os.path.join(OUT, 'summary.json'), 'w'), indent=1); print(json.dumps(S, indent=1))
