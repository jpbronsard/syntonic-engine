# Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
# Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
"""Does the law's adaptive window still beat a frozen window once the servo lead is on? Written 1 Nov 2026 after
eg732_hold.py (sha256 64134c25...), before any run of this file except the smoke test; go from JP ("Oui lance le").

Why: tonight's 30 % gain came from the servo lead and faster frames, classical engineering. What the Syntony law itself
contributes on this bench is its adaptive window. Earlier tonight (exploratory, accel_window_graph note) the adaptive
window beat the same engine with its window frozen at its own mean by 7 to 12 % (333 Hz frames, 1 kHz loop) and 10 to 14 %
(1 kHz frames, 2 kHz loop, one seed), without any lead. This test asks whether that survives in the best configuration.

Setting: eg732_hold.HoldEngine unchanged (engine 0.2.8, authority limit), 2 kHz loop, fixed 1.5 ms frames (the fastest
grid frame under the servos' 732 Hz), baton arm RJ1-RJ6, excerpts A and B, steady and x10 noise step (10 to 20 s),
seeds 1 to 5: 20 paired runs.
Arms: F15 (no lead) and F15D (lead (servo_tau + D) x planned velocity), each adaptive; then F15z and F15Dz, the same
engines with the window frozen per joint at the mean of their own adaptive run over the scored 30 s (floor = tau_max =
tau = tau_target = that mean), the method of plan_intervention.py. Frozen arms share the run's encoder noise draws.
Gates: G4, F15D adaptive reproduces eg732_hold.py run by run to 1e-9 mm; G5, F15 adaptive reproduces eg732_lead.py
run by run to 1e-9 mm.

Predictions (judged; pooled = rms over the 20 runs of the per-run tip rms; gap = 1 - pooled adaptive / pooled frozen):
W1  Without the lead, the window still pays: gap(F15) >= 5 %, adaptive below frozen in >= 18 of 20 runs.
W2  With the lead, the window still pays: gap(F15D) >= 5 %, adaptive below frozen in >= 18 of 20 runs.
W3  The lead does not eat most of the window's share: gap(F15D) >= 0.5 x gap(F15).
Described, not judged: gaps by cell (steady against noise step), frozen window values, joint error.

Run: python window_vs_frozen.py <repo_root> smoke | run <excerpt> <cond> | summarize <eg732_lead results> <eg732_hold results>
"""
import os, sys, json, time, hashlib
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.signal import butter, sosfiltfilt

ROOT = os.path.abspath(sys.argv[1]); MODE = sys.argv[2]
HERE_W = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE_W)
_argv = sys.argv; sys.argv = [_argv[0], ROOT, 'import']
import eg732_hold as H                                                   # HoldEngine, FR, bench constants
sys.argv = _argv
OUT = os.path.join(HERE_W, 'results_window'); os.makedirs(OUT, exist_ok=True)
SHA = hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest()
TABLE = {'F15': ([1.5], 0, 0), 'F15D': ([1.5], 0, 2)}


def simulate(name, cond, seed, arms, freeze=None, t_end=30.0):
    dt = 5e-4; joints = list(range(6))
    t0 = float(H.Z[f'{name}_t0']); tt = np.arange(t0 - 3.0, t0 + t_end, dt)
    J = 6; A = len(arms)
    S = CubicSpline(H.Z[f'{name}_t30'], H.Z[f'{name}_q30'][:, joints], axis=0)
    qJ = S(tt); qdJ = S(tt, 1)
    tauJ = CubicSpline(H.Z[f'{name}_t30'], H.Z[f'{name}_tau30'][:, joints], axis=0)(tt)
    sd = np.full(len(tt), H.ENC_SD)
    if cond == 'step': sd[(tt >= t0 + 10) & (tt < t0 + 20)] *= np.sqrt(10.0)
    rng = np.random.default_rng(seed)
    tauS = np.tile(H.TAU_S[joints], A); w0 = np.tile(H.W0[joints], A)
    k_inv = (1.0 * H.DEG) / tauS; hb = 0.5 * 0.3 * H.DEG
    allowed = np.zeros((A * J, len(H.FR)), bool); source = np.zeros(A * J, int); mode = np.zeros(A * J, int)
    for i, lab in enumerate(arms):
        fr, src, md = TABLE[lab]
        for f in fr: allowed[i * J:(i + 1) * J, H.FR.index(f)] = True
        source[i * J:(i + 1) * J] = src; mode[i * J:(i + 1) * J] = md
    eng = H.HoldEngine(np.tile(qJ[0], A), dt, w0, allowed, source, mode)
    if freeze is not None:
        tf = np.concatenate([np.asarray(freeze[lab], float) for lab in arms])
        eng.floor = tf.copy(); eng.tau_max = tf.copy(); eng.tau = tf.copy(); eng.tau_target = tf.copy()
    th_g = np.tile(qJ[0], A); every = int(round(0.01 / dt)); keep = []
    err2 = np.zeros(A * J); m = 0; wsum = np.zeros(A * J)
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
            err2 += (q - th_out) ** 2; m += 1; wsum += np.broadcast_to(eng.tau, (A * J,))
            if k % every == 0: keep.append((th_out.copy(), q.copy()))
    out = {}
    for i, lab in enumerate(arms):
        s = slice(i * J, (i + 1) * J)
        e = np.array([H.tip(a[s]) - H.tip(b[s]) for a, b in keep])
        out[lab] = dict(tip_mm=float(1e3 * np.sqrt(np.mean(np.sum(e ** 2, 1)))),
                        joint_deg=float(np.mean(np.sqrt(err2[s] / m)) / H.DEG),
                        window_mean_ms=(wsum[s] / m * 1e3).tolist())
    return out


if __name__ == '__main__':
    arms = ['F15', 'F15D']
    if MODE == 'smoke':
        t = time.time(); r = simulate('A', 'steady', 1, arms, t_end=1.5)
        r2 = simulate('A', 'steady', 1, arms, freeze={a: np.array(r[a]['window_mean_ms']) * 1e-3 for a in arms}, t_end=1.5)
        print('smoke ok', sorted(r), sorted(r2), round(time.time() - t, 1), 's'); sys.exit(0)
    if MODE == 'run':
        name, cond = sys.argv[3], sys.argv[4]
        for seed in H.SEEDS:
            t = time.time()
            ad = simulate(name, cond, seed, arms)
            fz = simulate(name, cond, seed, arms, freeze={a: np.array(ad[a]['window_mean_ms']) * 1e-3 for a in arms})
            row = dict(key=f'{name} {cond} {seed}', sha=SHA, F15=ad['F15'], F15D=ad['F15D'], F15z=fz['F15'], F15Dz=fz['F15D'])
            with open(os.path.join(OUT, f'runs_{name}_{cond}.jsonl'), 'a') as f: f.write(json.dumps(row) + '\n')
            print(name, cond, seed, {a: round(row[a]['tip_mm'], 3) for a in ('F15', 'F15z', 'F15D', 'F15Dz')}, round(time.time() - t), 's', flush=True)
        sys.exit(0)
    if MODE == 'summarize':
        def load(d):
            runs = {}
            for fn in sorted(os.listdir(d)):
                if fn.startswith('runs_'):
                    for line in open(os.path.join(d, fn)):
                        x = json.loads(line); runs[x.pop('key')] = x
            return runs
        runs = load(OUT); lead = load(sys.argv[3]); hold = load(sys.argv[4])
        assert len(runs) == 20 and len(lead) == 20 and len(hold) == 20
        g4 = max(abs(runs[k]['F15D']['tip_mm'] - hold[k]['F15D']['tip_mm']) for k in runs)
        g5 = max(abs(runs[k]['F15']['tip_mm'] - lead[k]['F15']['tip_mm']) for k in runs)
        G = dict(G4=dict(max_diff_mm=g4, ok=bool(g4 <= 1e-9)), G5=dict(max_diff_mm=g5, ok=bool(g5 <= 1e-9)))
        ok = G['G4']['ok'] and G['G5']['ok']
        AR = ('F15', 'F15z', 'F15D', 'F15Dz')
        P = {a: float(np.sqrt(np.mean([runs[k][a]['tip_mm'] ** 2 for k in runs]))) for a in AR}
        gap = {'F15': 1 - P['F15'] / P['F15z'], 'F15D': 1 - P['F15D'] / P['F15Dz']}
        win = lambda a, b: int(sum(runs[k][a]['tip_mm'] < runs[k][b]['tip_mm'] for k in runs))
        V = dict(W1=bool(ok and gap['F15'] >= 0.05 and win('F15', 'F15z') >= 18),
                 W2=bool(ok and gap['F15D'] >= 0.05 and win('F15D', 'F15Dz') >= 18),
                 W3=bool(ok and gap['F15D'] >= 0.5 * gap['F15']))
        cells = {}
        for c in ('A steady', 'A step', 'B steady', 'B step'):
            pc = {a: float(np.sqrt(np.mean([runs[k][a]['tip_mm'] ** 2 for k in runs if k.startswith(c)]))) for a in AR}
            cells[c] = dict(pooled=pc, gap_F15=1 - pc['F15'] / pc['F15z'], gap_F15D=1 - pc['F15D'] / pc['F15Dz'])
        frozen_ms = {a: np.mean([runs[k][a]['window_mean_ms'] for k in runs], 0).tolist() for a in ('F15', 'F15D')}
        S = dict(script_sha256=SHA, gates=G, pooled_tip_mm=P, gap=gap, wins=dict(F15_vs_F15z=win('F15', 'F15z'),
                 F15D_vs_F15Dz=win('F15D', 'F15Dz')), verdict=V, cells=cells, mean_window_ms_by_joint=frozen_ms,
                 joint_deg={a: float(np.mean([runs[k][a]['joint_deg'] for k in runs])) for a in AR})
        json.dump(S, open(os.path.join(OUT, 'summary.json'), 'w'), indent=1); print(json.dumps(S, indent=1))
