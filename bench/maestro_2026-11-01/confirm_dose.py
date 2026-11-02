# Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
# Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
"""A more marked effect: large noise changes and frequent noise bursts, on fresh seeds. Written 1 Nov 2026, after the
exploratory explore_dose.py / explore_dose2.py passes and before any run of this file except the smoke test; goal set
by JP ("comment obtenir un effet plus marque ?").

Exploratory findings this is built on (seeds 21 to 23, F15D setting, not judged):
- The formula's lead over a fixed window grows with the size of the noise change: against the 5 ms window, +4.4 %
  steady, +8.5 % at variance x10, +32 % at x100; against the best fixed window for the x100 profile (13 ms, chosen on
  seed 21), +13.5 % on seeds 22-23 (4 of 4); with the e-geodesic travel (E64), +15.8 % (4 of 4).
- Under frequent bursts (variance x30, 0.5 s every 2.5 s) the published transport loses to the best fixed window (5 ms,
  -5.2 %), and the e-geodesic travel turns it around: E64 +4.2 % over the fixed window, +9.0 % over the published
  transport (4 of 4).
The fixed rivals are the best of 3, 5, 8, 13, 21 ms for each profile on seed 21; they are fixed here, before any run.

Setting: window_transport.TransportEngine (engine 0.2.8, 2 kHz loop, fixed 1.5 ms frames, lead (servo_tau + D) x
planned velocity), baton arm, excerpts A and B, seeds 26 to 30 (never used): 10 paired runs per profile.
  step100  encoder read variance x100 from 10 to 20 s.     Arms: BASE, E64, FIX13.
  bursts   variance x30 for 0.5 s every 2.5 s, 5 to 25 s.   Arms: BASE, E64, FIX5.
BASE = published transport; E64 = 1/tau relaxing with T = 64 ms (longest identification lag); FIXv = window fixed at
v ms on every joint (floor = tau_max = tau = tau_target).
Gate: G8, on step100 A seed 21 this script reproduces explore_dose.py's BASE and E64 and explore_dose2.py's FIX13 to
1e-9 mm.

Predictions (judged; pooled = rms over the 10 runs of the per-run tip rms):
K1  Large change, formula alone: pooled BASE <= 0.92 x pooled FIX13, BASE below FIX13 in >= 9 of 10.
K2  Large change, with the e-geodesic travel: pooled E64 <= 0.90 x pooled FIX13, below in >= 9 of 10.
K3  Bursts, the e-geodesic travel against the published transport: pooled E64 <= 0.95 x pooled BASE, below in >= 9 of 10.
K4  Bursts, the e-geodesic travel against the best fixed window: pooled E64 <= pooled FIX5, below in >= 7 of 10.

Run: python confirm_dose.py <repo_root> smoke | gate | run <profile> <excerpt> | summarize
"""
import os, sys, json, time, hashlib
import numpy as np
from scipy.interpolate import CubicSpline

ROOT = os.path.abspath(sys.argv[1]); MODE = sys.argv[2]; _args = list(sys.argv)
HERE_W = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE_W)
sys.argv = [_args[0], ROOT, 'import']
import window_transport as WT
sys.argv = _args
H = WT.H
OUT = os.path.join(HERE_W, 'results_confirm_dose'); os.makedirs(OUT, exist_ok=True)
SHA = hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest()
SEEDS = [26, 27, 28, 29, 30]
PROFILES = {'step100': ['BASE', 'E64', 'FIX13'], 'bursts': ['BASE', 'E64', 'FIX5']}


def profile(name, trel):
    m = np.ones_like(trel)
    if name == 'step100': m[(trel >= 10) & (trel < 20)] = 10.0
    if name == 'bursts':
        ph = (trel - 5.0) % 2.5
        m[(trel >= 5) & (trel < 25) & (ph < 0.5)] = np.sqrt(30.0)
    return m


def arm_spec(lab):
    if lab == 'BASE': return 0, 0.0, None
    if lab == 'E64': return 1, 0.064, None
    return 0, 0.0, float(lab[3:]) * 1e-3


def simulate(prof, name, seed, arms, t_end=30.0):
    dt = 5e-4; joints = list(range(6)); J = 6; A = len(arms)
    t0 = float(H.Z[f'{name}_t0']); tt = np.arange(t0 - 3.0, t0 + t_end, dt)
    S = CubicSpline(H.Z[f'{name}_t30'], H.Z[f'{name}_q30'][:, joints], axis=0); qJ = S(tt); qdJ = S(tt, 1)
    tauJ = CubicSpline(H.Z[f'{name}_t30'], H.Z[f'{name}_tau30'][:, joints], axis=0)(tt)
    sd = H.ENC_SD * profile(prof, tt - t0)
    rng = np.random.default_rng(seed)
    tauS = np.tile(H.TAU_S[joints], A); w0 = np.tile(H.W0[joints], A)
    k_inv = (1.0 * H.DEG) / tauS; hb = 0.5 * 0.3 * H.DEG
    allowed = np.zeros((A * J, len(H.FR)), bool); allowed[:, H.FR.index(1.5)] = True
    specs = [arm_spec(a) for a in arms]
    tr = np.repeat([s[0] for s in specs], J); T = np.repeat([s[1] for s in specs], J)
    eng = WT.TransportEngine(np.tile(qJ[0], A), dt, w0, allowed, np.zeros(A * J, int), np.full(A * J, 2), transport=tr, T=T)
    eng.tau_max = np.full(A * J, float(eng.tau_max)); eng.floor = np.asarray(eng.floor, float).copy()
    for i, sp in enumerate(specs):
        if sp[2] is not None:
            s = slice(i * J, (i + 1) * J); v = sp[2]
            eng.floor[s] = v; eng.tau_max[s] = v; eng.tau[s] = v; eng.tau_target[s] = v
    th_g = np.tile(qJ[0], A); every = int(round(0.01 / dt)); keep = []
    for k in range(len(tt)):
        load = np.tile(tauJ[k], A); q = np.tile(qJ[k], A); qd = np.tile(qdJ[k], A)
        th_out = th_g - load * k_inv - hb * np.tanh(load / (0.02 * tauS))
        na = np.tile(rng.normal(0.0, 1.0, J), A); nb = np.tile(rng.normal(0.0, 1.0, J), A)
        ra = np.round((th_out + na * sd[k]) / H.ENC_LSB) * H.ENC_LSB
        rb = np.round((th_out + nb * sd[k]) / H.ENC_LSB) * H.ENC_LSB
        cmd = eng.step(q, qd, ra, rb)
        w_eff = w0 * np.maximum(0.05, 1.0 - np.abs(load) / tauS)
        th_g = th_g + dt * np.clip((cmd - th_g) / H.TS, -w_eff, w_eff)
        if tt[k] >= t0 - 1e-9 and k % every == 0: keep.append((th_out.copy(), q.copy()))
    out = {}
    for i, lab in enumerate(arms):
        s = slice(i * J, (i + 1) * J)
        e = np.array([H.tip(a[s]) - H.tip(b[s]) for a, b in keep])
        out[lab] = float(1e3 * np.sqrt(np.mean(np.sum(e ** 2, 1))))
    return out


if __name__ == '__main__':
    if MODE == 'smoke':
        t = time.time(); r = simulate('bursts', 'A', 26, PROFILES['bursts'], t_end=1.5); print('smoke ok', sorted(r), round(time.time() - t, 1), 's'); sys.exit(0)
    if MODE == 'gate':
        def find(fn, key):
            for line in open(fn):
                x = json.loads(line)
                if x['profile'] == 'step100' and x['excerpt'] == 'A' and x['seed'] == 21: return x[key]
        refB = find(os.path.join(HERE_W, 'results_dose', 'runs.jsonl'), 'BASE'); refE = find(os.path.join(HERE_W, 'results_dose', 'runs.jsonl'), 'E64')
        refF = find(os.path.join(HERE_W, 'results_dose_sweep', 'runs.jsonl'), 'FIX13')
        m = simulate('step100', 'A', 21, ['BASE', 'E64', 'FIX13'])
        d = max(abs(m['BASE'] - refB), abs(m['E64'] - refE), abs(m['FIX13'] - refF))
        g = dict(G8=dict(ref=[refB, refE, refF], this=[m['BASE'], m['E64'], m['FIX13']], max_diff_mm=d, ok=bool(d <= 1e-9)), script_sha256=SHA)
        json.dump(g, open(os.path.join(OUT, 'gates.json'), 'w'), indent=1); print(json.dumps(g)); sys.exit(0)
    if MODE == 'run':
        prof, name = sys.argv[3], sys.argv[4]
        for seed in SEEDS:
            t = time.time(); r = simulate(prof, name, seed, PROFILES[prof])
            with open(os.path.join(OUT, 'runs.jsonl'), 'a') as f: f.write(json.dumps(dict(profile=prof, excerpt=name, seed=seed, sha=SHA, **r)) + '\n')
            print(prof, name, seed, {a: round(v, 3) for a, v in r.items()}, round(time.time() - t), 's', flush=True)
        sys.exit(0)
    if MODE == 'summarize':
        rows = [json.loads(l) for l in open(os.path.join(OUT, 'runs.jsonl'))]
        G = json.load(open(os.path.join(OUT, 'gates.json'))); ok = G['G8']['ok']
        S = dict(script_sha256=SHA, gates=G)
        for prof, arms in PROFILES.items():
            R = [r for r in rows if r['profile'] == prof]; assert len(R) == 10, (prof, len(R))
            P = {a: float(np.sqrt(np.mean([r[a] ** 2 for r in R]))) for a in arms}
            W = {f'{a}_vs_{b}': int(sum(r[a] < r[b] for r in R)) for a in arms for b in arms if a != b}
            S[prof] = dict(pooled=P, wins=W, by_excerpt={n: {a: float(np.sqrt(np.mean([r[a] ** 2 for r in R if r['excerpt'] == n]))) for a in arms} for n in 'AB'})
        s1, s2 = S['step100'], S['bursts']
        S['verdict'] = dict(
            K1=bool(ok and s1['pooled']['BASE'] <= 0.92 * s1['pooled']['FIX13'] and s1['wins']['BASE_vs_FIX13'] >= 9),
            K2=bool(ok and s1['pooled']['E64'] <= 0.90 * s1['pooled']['FIX13'] and s1['wins']['E64_vs_FIX13'] >= 9),
            K3=bool(ok and s2['pooled']['E64'] <= 0.95 * s2['pooled']['BASE'] and s2['wins']['E64_vs_BASE'] >= 9),
            K4=bool(ok and s2['pooled']['E64'] <= s2['pooled']['FIX5'] and s2['wins']['E64_vs_FIX5'] >= 7))
        S['ratios'] = dict(step100_BASE_FIX13=s1['pooled']['BASE'] / s1['pooled']['FIX13'], step100_E64_FIX13=s1['pooled']['E64'] / s1['pooled']['FIX13'],
                           bursts_E64_BASE=s2['pooled']['E64'] / s2['pooled']['BASE'], bursts_E64_FIX5=s2['pooled']['E64'] / s2['pooled']['FIX5'],
                           bursts_BASE_FIX5=s2['pooled']['BASE'] / s2['pooled']['FIX5'])
        json.dump(S, open(os.path.join(OUT, 'summary.json'), 'w'), indent=1); print(json.dumps(S, indent=1))
