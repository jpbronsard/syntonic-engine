# Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
# Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
"""kappa on MAESTRO 1: a static map, and a dynamic kappa that learns online from the encoders. Written 2 Nov 2026,
before any run of this file except the smoke test; go from JP ("k dynamique", then "go"). Supersedes kappa_scale.py
(sha256 235ebf89..., frozen, never run past its smoke test).

The engine takes kappa_H = (2H Gamma(2H+1))^(-1/(2H+1)) from Theorem U', with H identified live (0.95 to 0.99 on the
baton arm: kappa_H 0.62 to 0.66). Dynamic kappa in the corpus is an operational state kappa_op (META hierarchy). The
Register V1.3 warns that regulating kappa_op on kappa_obs = tau sqrt(lambda / sigma^2), computed from the same moment
streams as the window, is circular (kappa_obs = kappa_op algebraically). Here the dynamic kappa learns from an outcome
the window does not define: the encoder tracking error e = plan - measured angle.

Arms, all as columns sharing every encoder read; the law's target is multiplied by a factor c (kappa -> c x kappa_H):
  C0.5, C0.7, C1.0, C1.4, C2.0   static factors (C1.0 is the published engine)
  DYN   per joint, c starts at 1 and is moved by extremum seeking on the encoder error:
        dither: c_eff = c x exp(+a) for 1 s, then c x exp(-a) for 1 s (period 2 s = twice the engine's 1 s averaging
        time), a = ln(1.5) / 2 (half the x1.5 oracle band of the 27 Oct protocol);
        cost: mean of e^2 over each 1 s half; g = (J+ - J-) / (J+ + J-);
        update once per period: ln c -= (a / 4) x sign(g); c kept in [0.25, 4]; starts after the 1 s warm-up.
        These four constants are declared here, not fitted.
Settings (baton arm, excerpts A and B, steady and x10 noise step 10 to 20 s, seeds 31 to 35 never used: 20 paired runs
each): S1 published bench (1 kHz loop, fixed 3 ms frames, no lead); S2 tonight's best (2 kHz loop, fixed 1.5 ms frames,
lead (servo_tau + D) x planned velocity). Engine 0.2.8, identification 'raw', published transport.
Gates: G9, S1 C1.0 on 12 joints, A steady seed 5, reproduces examples/servo_bench.py (0.683910 mm) to 1e-9 mm; G10, S2
C1.0, A steady seed 1, reproduces eg732_hold.py's F15D to 1e-9 mm.

Predictions (judged; pooled = rms over the 20 runs of the per-run tip rms; c* = the static factor with the lowest pooled):
Q1  S1: c* in {0.7, 1, 1.4} and regret(C1.0) = pooled(C1.0) / pooled(c*) - 1 < 10 %.
Q2  S2: the same.
D1  S1 and S2: the dynamic kappa does no harm: pooled DYN <= pooled C1.0 in each setting.
D2  S1 and S2: the dynamic kappa finds the best static factor: pooled DYN <= 1.02 x pooled c* in each setting.
Described, not judged: the parabola optimum in ln c; DYN's final factor per joint (mean over the last 10 s); per cell.

Run: python kappa_dyn.py <repo_root> smoke | gate | run <S1|S2> <excerpt> <cond> | summarize
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
OUT = os.path.join(HERE_W, 'results_kdyn'); os.makedirs(OUT, exist_ok=True)
SHA = hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest()
CS = [0.5, 0.7, 1.0, 1.4, 2.0]
ARMS = [f'C{c}' for c in CS] + ['DYN']
SEEDS = [31, 32, 33, 34, 35]
SETTINGS = {'S1': dict(dt=1e-3, frame=3.0, mode=0), 'S2': dict(dt=5e-4, frame=1.5, mode=2)}
PERIOD, AMP = 2.0, np.log(1.5) / 2
STEP, CMIN, CMAX = AMP / 4, 0.25, 4.0


class KEngine(H.HoldEngine):
    """HoldEngine with fixed frames; the law's target is multiplied by c per column; dynamic columns learn c online."""

    def __init__(self, *a, c, dyn, **k):
        super().__init__(*a, **k)
        n = self.n
        self.lnc = np.log(np.asarray(c, float)); self.dyn = np.asarray(dyn, bool)
        self.jp = np.zeros(n); self.np_ = 0; self.jm = np.zeros(n); self.nm = 0; self.half_prev = 0

    def c_eff(self):
        t = self.k * self.dt
        if self.k <= self.warm: return np.exp(self.lnc)
        half = int(t / (PERIOD / 2)) % 2
        sgn = 1.0 if half == 0 else -1.0
        return np.exp(self.lnc + np.where(self.dyn, sgn * AMP, 0.0))

    def es_update(self, e):
        if self.k <= self.warm: return
        t = self.k * self.dt; half = int(t / (PERIOD / 2)) % 2
        if half == 0 and self.half_prev == 1 and self.np_ > 0 and self.nm > 0:      # a full period (+ then -) has ended
            Jp = self.jp / self.np_; Jm = self.jm / self.nm
            g = (Jp - Jm) / np.maximum(Jp + Jm, 1e-30)
            self.lnc = np.where(self.dyn, np.clip(self.lnc - STEP * np.sign(g), np.log(CMIN), np.log(CMAX)), self.lnc)
            self.jp = np.zeros(self.n); self.np_ = 0; self.jm = np.zeros(self.n); self.nm = 0
        if half == 0: self.jp += e * e; self.np_ += 1
        else: self.jm += e * e; self.nm += 1
        self.half_prev = half

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
            ts = self.c_eff() * kappa_H(Hn) * (np.maximum(s2, 1e-30) / np.maximum(lam, 1e-30)) ** (1.0 / (2 * Hn + 1))
            self.law = np.where(valid, np.minimum(ts, self.tau_max), self.law)
            ts = np.clip(ts, self.floor, self.tau_max)
            self.tau_target = np.where(valid, ts, self.tau_target); self.H = np.where(valid, Hn, self.H)
            rr = 1.0 / self.tau; dr_max = self.tau_rate * dt / self.tau ** 2
            rr = rr + np.clip(1.0 / self.tau_target - rr, -dr_max, dr_max); self.tau = 1.0 / rr
        self.es_update(e)
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


def simulate(setting, name, cond, seed, arms, joints=None, t_end=30.0):
    st = SETTINGS[setting]; dt = st['dt']
    joints = list(range(6)) if joints is None else joints; J = len(joints); A = len(arms)
    t0 = float(H.Z[f'{name}_t0']); tt = np.arange(t0 - 3.0, t0 + t_end, dt)
    S = CubicSpline(H.Z[f'{name}_t30'], H.Z[f'{name}_q30'][:, joints], axis=0); qJ = S(tt); qdJ = S(tt, 1)
    tauJ = CubicSpline(H.Z[f'{name}_t30'], H.Z[f'{name}_tau30'][:, joints], axis=0)(tt)
    sd = np.full(len(tt), H.ENC_SD)
    if cond == 'step': sd[(tt >= t0 + 10) & (tt < t0 + 20)] *= np.sqrt(10.0)
    rng = np.random.default_rng(seed)
    tauS = np.tile(H.TAU_S[joints], A); w0 = np.tile(H.W0[joints], A)
    k_inv = (1.0 * H.DEG) / tauS; hb = 0.5 * 0.3 * H.DEG
    allowed = np.zeros((A * J, len(H.FR)), bool); allowed[:, H.FR.index(st['frame'])] = True
    c0 = np.repeat([1.0 if a == 'DYN' else float(a[1:]) for a in arms], J); dyn = np.repeat([a == 'DYN' for a in arms], J)
    eng = KEngine(np.tile(qJ[0], A), dt, w0, allowed, np.zeros(A * J, int), np.full(A * J, st['mode']), c=c0, dyn=dyn)
    th_g = np.tile(qJ[0], A); every = int(round(0.01 / dt)); keep = []; clog = []
    for k in range(len(tt)):
        load = np.tile(tauJ[k], A); q = np.tile(qJ[k], A); qd = np.tile(qdJ[k], A)
        th_out = th_g - load * k_inv - hb * np.tanh(load / (0.02 * tauS))
        na = np.tile(rng.normal(0.0, 1.0, J), A); nb = np.tile(rng.normal(0.0, 1.0, J), A)
        ra = np.round((th_out + na * sd[k]) / H.ENC_LSB) * H.ENC_LSB
        rb = np.round((th_out + nb * sd[k]) / H.ENC_LSB) * H.ENC_LSB
        cmd = eng.step(q, qd, ra, rb)
        w_eff = w0 * np.maximum(0.05, 1.0 - np.abs(load) / tauS)
        th_g = th_g + dt * np.clip((cmd - th_g) / H.TS, -w_eff, w_eff)
        if tt[k] >= t0 - 1e-9 and k % every == 0:
            keep.append((th_out.copy(), q.copy())); clog.append(np.exp(eng.lnc).copy())
    clog = np.array(clog); out = {}
    for i, lab in enumerate(arms):
        s = slice(i * J, (i + 1) * J)
        e = np.array([H.tip(a[s][:6]) - H.tip(b[s][:6]) for a, b in keep])
        out[lab] = dict(tip_mm=float(1e3 * np.sqrt(np.mean(np.sum(e ** 2, 1)))))
        if lab == 'DYN': out[lab]['c_final'] = np.exp(np.mean(np.log(clog[-1000:, s][:, :6]), 0)).tolist()
    return out


if __name__ == '__main__':
    if MODE == 'smoke':
        t = time.time(); r1 = simulate('S1', 'A', 'step', 31, ARMS, t_end=4.0); r2 = simulate('S2', 'A', 'step', 31, ARMS, t_end=1.5)
        print('smoke ok', sorted(r1), r1['DYN'].get('c_final'), round(time.time() - t, 1), 's'); sys.exit(0)
    if MODE == 'gate':
        g9 = simulate('S1', 'A', 'steady', 5, ['C1.0'], joints=list(range(12)))['C1.0']['tip_mm']
        ref10 = None
        for line in open(os.path.join(HERE_W, 'results_hold', 'runs_A_steady.jsonl')):
            x = json.loads(line)
            if x['key'] == 'A steady 1': ref10 = x['F15D']['tip_mm']
        g10 = simulate('S2', 'A', 'steady', 1, ['C1.0'])['C1.0']['tip_mm']
        G = dict(G9=dict(ref=0.6839103668504339, this=g9, ok=bool(abs(g9 - 0.6839103668504339) <= 1e-9)),
                 G10=dict(ref=ref10, this=g10, ok=bool(abs(g10 - ref10) <= 1e-9)), script_sha256=SHA)
        json.dump(G, open(os.path.join(OUT, 'gates.json'), 'w'), indent=1); print(json.dumps(G)); sys.exit(0)
    if MODE == 'run':
        setting, name, cond = sys.argv[3], sys.argv[4], sys.argv[5]
        for seed in SEEDS:
            t = time.time(); r = simulate(setting, name, cond, seed, ARMS)
            with open(os.path.join(OUT, 'runs.jsonl'), 'a') as f: f.write(json.dumps(dict(setting=setting, key=f'{name} {cond} {seed}', sha=SHA, **r)) + '\n')
            print(setting, name, cond, seed, {a: round(v['tip_mm'], 3) for a, v in r.items()}, 'c_dyn', np.round(r['DYN']['c_final'], 2).tolist(), round(time.time() - t), 's', flush=True)
        sys.exit(0)
    if MODE == 'summarize':
        rows = [json.loads(l) for l in open(os.path.join(OUT, 'runs.jsonl'))]
        G = json.load(open(os.path.join(OUT, 'gates.json'))); ok = G['G9']['ok'] and G['G10']['ok']
        S = dict(script_sha256=SHA, gates=G); V = {}
        for st in ('S1', 'S2'):
            R = [r for r in rows if r['setting'] == st]; assert len(R) == 20, (st, len(R))
            P = {a: float(np.sqrt(np.mean([r[a]['tip_mm'] ** 2 for r in R]))) for a in ARMS}
            stat = {a: P[a] for a in ARMS if a != 'DYN'}; cstar = min(stat, key=stat.get); regret = P['C1.0'] / P[cstar] - 1
            idx = np.argsort([P[f'C{c}'] for c in CS])[:3]; x = np.log([CS[i] for i in idx]); y = np.array([P[f'C{CS[i]}'] for i in idx])
            co = np.polyfit(x, y, 2); c_par = float(np.exp(-co[1] / (2 * co[0]))) if co[0] > 0 else None
            cfin = np.array([r['DYN']['c_final'] for r in R])
            S[st] = dict(pooled=P, c_star=cstar, regret_C1=regret, c_parabola=c_par,
                         dyn_final_c_geomean_by_joint=np.exp(np.mean(np.log(cfin), 0)).tolist(),
                         wins_DYN_vs_C1=int(sum(r['DYN']['tip_mm'] < r['C1.0']['tip_mm'] for r in R)),
                         cells={cl: {a: float(np.sqrt(np.mean([r[a]['tip_mm'] ** 2 for r in R if r['key'].startswith(cl)]))) for a in ARMS}
                                for cl in ('A steady', 'A step', 'B steady', 'B step')})
            V[f'Q_{st}'] = bool(ok and cstar in ('C0.7', 'C1.0', 'C1.4') and regret < 0.10)
            V[f'D1_{st}'] = bool(ok and P['DYN'] <= P['C1.0'])
            V[f'D2_{st}'] = bool(ok and P['DYN'] <= 1.02 * P[cstar])
        S['verdict'] = dict(Q1=V['Q_S1'], Q2=V['Q_S2'], D1=V['D1_S1'] and V['D1_S2'], D2=V['D2_S1'] and V['D2_S2'], by_setting=V)
        json.dump(S, open(os.path.join(OUT, 'summary.json'), 'w'), indent=1); print(json.dumps(S, indent=1))
