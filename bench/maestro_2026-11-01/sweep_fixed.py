"""Exploratory, after the W1-W3 verdicts: is 'frozen at its own mean' a fair rival once the lead is on? Fixed windows,
the same value on every joint, swept, against the adaptive window; F15D setting (lead servo_tau + D, 1.5 ms frames,
2 kHz loop). Seed 1, one run per cell; then the best uniform value on seeds 2 to 5."""
import os, sys, json, numpy as np
import os
R = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')); sys.argv = ['x', R, 'import']
import window_vs_frozen as W
VALS = [1.111, 1.5, 2.0, 3.0, 5.0, 8.0, 13.0]
def run(name, cond, seed, vals):
    W.TABLE.update({f'U{v}': ([1.5], 0, 2) for v in vals})
    arms = ['F15D'] + [f'U{v}' for v in vals]
    ad = W.simulate(name, cond, seed, ['F15D'])
    fz = W.simulate(name, cond, seed, [f'U{v}' for v in vals], freeze={f'U{v}': np.full(6, v * 1e-3) for v in vals})
    return {'adaptive': ad['F15D']['tip_mm'], **{f'{v} ms': fz[f'U{v}']['tip_mm'] for v in vals}}
mode = sys.argv[-1] if len(sys.argv) > 3 else None
out = {}
for c in ('A steady', 'A step', 'B steady', 'B step'):
    n, cond = c.split(); r = run(n, cond, 1, VALS); out[c] = r
    print(c, {k: round(v, 3) for k, v in r.items()}, flush=True)
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results_window', 'sweep_seed1.json'), 'w'), indent=1)
