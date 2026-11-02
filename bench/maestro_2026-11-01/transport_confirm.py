# Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
# Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
"""Confirmation on fresh seeds: the window's travel along the e-geodesic, with its time constant taken from a declared
engine constant instead of a value chosen by hand. Written 1 Nov 2026 after window_transport.py (sha256 101655fa...),
before any run of this file except the smoke test; go from JP ("roule 1").

window_transport.py found the e coordinate (theta = 1/tau relaxing) better than the tau coordinate at 1 s (-19 %) and
0.1 s (-3 %), and, exploratory, E at 0.1 s 2.9 % better than the published transport pooled but only in 10 of 20 runs;
0.1 s was chosen by hand. Here the time constant is declared before the run, from the engine's own design constants:

  T = the longest lag of the 'raw' identification, max(lags) = 64 ms (engine default lags 2 to 64 ms).
  Rationale: the law's target is identified from the residual's structure function at lags up to 64 ms; the window
  travels on the slowest time scale the identification resolves. Stated openly: this rule was chosen after seeing that
  0.1 s helped and 1 s did not, so only the fresh-seed result below counts.

Setting: window_transport.TransportEngine unchanged (F15D: engine 0.2.8, 2 kHz loop, fixed 1.5 ms frames, lead
(servo_tau + D) x planned velocity), baton arm, excerpts A and B, steady and x10 noise step, seeds 21 to 25 (never
used tonight): 20 paired runs. Arms: BASE (published transport), E64 and N64 (theta / tau relaxation, T = 64 ms), E01
(theta relaxation, T = 0.1 s, the exploratory value, for replication).
Gate: G7, BASE on A steady seed 1 reproduces window_transport.py's BASE value to 1e-9 mm.

Predictions (judged; pooled = rms over the 20 runs of the per-run tip rms):
C1  The declared rule beats the published transport: pooled E64 <= 0.98 x pooled BASE, and E64 below BASE in >= 14
    of 20 runs.
C2  The coordinate effect replicates at the declared constant: pooled E64 <= pooled N64, E64 below N64 in >= 15 of 20.
C3  The exploratory 0.1 s result replicates on fresh seeds: pooled E01 <= 0.98 x pooled BASE.
Described, not judged: by cell, the noise-step transients (10 to 12 s and 20 to 22 s), win counts for E01.

Run: python transport_confirm.py <repo_root> smoke | gate | run <excerpt> <cond> | summarize
"""
import os, sys, json, time, hashlib
import numpy as np

ROOT = os.path.abspath(sys.argv[1]); MODE = sys.argv[2]; _args = list(sys.argv)
HERE_W = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE_W)
sys.argv = [_args[0], ROOT, 'import']
import window_transport as WT
sys.argv = _args
OUT = os.path.join(HERE_W, 'results_confirm'); os.makedirs(OUT, exist_ok=True)
SHA = hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest()
T_DECL = 0.064
WT.ARMS.update({'E64': (1, T_DECL), 'N64': (2, T_DECL)})
ARMS = ['BASE', 'E64', 'N64', 'E01']
SEEDS = [21, 22, 23, 24, 25]

if __name__ == '__main__':
    if MODE == 'smoke':
        t = time.time(); r = WT.simulate('A', 'step', 21, ARMS, t_end=1.5); print('smoke ok', sorted(r), round(time.time() - t, 1), 's'); sys.exit(0)
    if MODE == 'gate':
        ref = None
        for line in open(os.path.join(HERE_W, 'results_transport', 'runs_A_steady.jsonl')):
            x = json.loads(line)
            if x['key'] == 'A steady 1': ref = x['BASE']['tip_mm']
        mine = WT.simulate('A', 'steady', 1, ['BASE'])['BASE']['tip_mm']
        g = dict(G7=dict(ref=ref, this=mine, ok=bool(abs(ref - mine) <= 1e-9)), script_sha256=SHA)
        json.dump(g, open(os.path.join(OUT, 'gates.json'), 'w'), indent=1); print(json.dumps(g)); sys.exit(0)
    if MODE == 'run':
        name, cond = sys.argv[3], sys.argv[4]
        for seed in SEEDS:
            t = time.time(); r = WT.simulate(name, cond, seed, ARMS)
            with open(os.path.join(OUT, f'runs_{name}_{cond}.jsonl'), 'a') as f: f.write(json.dumps(dict(key=f'{name} {cond} {seed}', sha=SHA, **r)) + '\n')
            print(name, cond, seed, {a: round(r[a]['tip_mm'], 3) for a in ARMS}, round(time.time() - t), 's', flush=True)
        sys.exit(0)
    if MODE == 'summarize':
        runs = {}
        for fn in sorted(os.listdir(OUT)):
            if fn.startswith('runs_'):
                for line in open(os.path.join(OUT, fn)):
                    x = json.loads(line); runs[x.pop('key')] = x
        assert len(runs) == 20
        G = json.load(open(os.path.join(OUT, 'gates.json'))); ok = G['G7']['ok']
        P = {a: float(np.sqrt(np.mean([runs[k][a]['tip_mm'] ** 2 for k in runs]))) for a in ARMS}
        steps = [k for k in runs if ' step ' in k]
        PT = {a: float(np.sqrt(np.mean([runs[k][a]['transient_mm'] ** 2 for k in steps]))) for a in ARMS}
        win = lambda x, y: int(sum(runs[k][x]['tip_mm'] < runs[k][y]['tip_mm'] for k in runs))
        V = dict(C1=bool(ok and P['E64'] <= 0.98 * P['BASE'] and win('E64', 'BASE') >= 14),
                 C2=bool(ok and P['E64'] <= P['N64'] and win('E64', 'N64') >= 15),
                 C3=bool(ok and P['E01'] <= 0.98 * P['BASE']))
        S = dict(script_sha256=SHA, gates=G, pooled_tip_mm=P, transient_step_mm=PT, verdict=V,
                 ratios=dict(E64_BASE=P['E64'] / P['BASE'], E64_N64=P['E64'] / P['N64'], E01_BASE=P['E01'] / P['BASE']),
                 wins=dict(E64_vs_BASE=win('E64', 'BASE'), E64_vs_N64=win('E64', 'N64'), E01_vs_BASE=win('E01', 'BASE')),
                 cells={c: {a: float(np.sqrt(np.mean([runs[k][a]['tip_mm'] ** 2 for k in runs if k.startswith(c)]))) for a in ARMS}
                        for c in ('A steady', 'A step', 'B steady', 'B step')})
        json.dump(S, open(os.path.join(OUT, 'summary.json'), 'w'), indent=1); print(json.dumps(S, indent=1))
