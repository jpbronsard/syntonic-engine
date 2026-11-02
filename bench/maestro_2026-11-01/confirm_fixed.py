"""Exploratory confirmation: the best single uniform fixed window chosen on seed 1 (pooled over the four cells), run on
held-out seeds 2 to 5, against the adaptive F15D values already recorded in results_window."""
import sys, json, numpy as np
V = float(sys.argv[-3]); n = sys.argv[-2]; cond = sys.argv[-1]
import os
R = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')); sys.argv = ['x', R, 'import']
import window_vs_frozen as W
W.TABLE['U'] = ([1.5], 0, 2)
for seed in (2, 3, 4, 5):
    r = W.simulate(n, cond, seed, ['U'], freeze={'U': np.full(6, V * 1e-3)})['U']['tip_mm']
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results_window', 'confirm_fixed.jsonl'), 'a') as f: f.write(json.dumps(dict(key=f'{n} {cond} {seed}', value_ms=V, tip_mm=r)) + '\n')
    print(n, cond, seed, V, round(r, 3), flush=True)
