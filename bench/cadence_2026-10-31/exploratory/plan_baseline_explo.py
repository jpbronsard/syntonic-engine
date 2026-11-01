# Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
# Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
"""Exploratory, before plan_intervention.py's predictions were written: baseline plan only, seed 1 (already used by
the exploratory passes), excerpts A and B. Checks the runner against omega_phase_rates.run at identity, the phase maps,
and reads at baseline the acceleration gain c_j against TS x the window actually used, and the acceleration part's
share of the tip mean square. No intervened plan is run here.   Usage: python plan_baseline_explo.py <bench dir>"""
import sys, numpy as np
sys.path.insert(0, sys.argv[1])
import plan_intervention as P
from omega_phase_rates import run as run_ref
for arm in ('EG', 'PL', 'U50', 'U70'):
    T, phi, om, omd = P.tempo(arm)
    t = np.array([-P.PRE, -1e-9, 0.0, T]); print(arm, 'T %.2f' % T, 'phi', np.round(phi(t), 4), 'omega', np.round(om(t), 4))
for name in ('A', 'B'):
    r = P.run(name, 1, 'syntonic', 'B'); ref = run_ref(name, 'steady', 1, 'syntonic', 6)
    print(name, 'identity: tip diff', abs(r['tip_mm'] - ref['tip_mm']), 'th diff', np.max(np.abs(r['th'] - ref['th'])))
    Js = P.jacobians(r['q']); f = P.fit(r, Js); tb = r['tau_used'].mean(0)
    print(' c (1e-5 s2)', np.round(f['c'] * 1e5, 2)); print(' TS x mean window (1e-5 s2)', np.round(P.TS * tb * 1e5, 2))
    print(' window used ms mean', np.round(tb * 1e3, 2), 'law below floor share', np.round(np.mean(r['law'] < 2.2e-3, 0), 2))
    print(' A2 / tip MS', round(f['A2'] / r['tip_ms'], 3), 'R2', np.round(f['R2'], 2))
    rf = P.run(name, 1, 'frozen', 'B', tau_fix=tb); ff = P.fit(rf, P.jacobians(rf['q']))
    print(' frozen: tip mm %.3f vs %.3f' % (rf['tip_mm'], r['tip_mm']), 'c', np.round(ff['c'] * 1e5, 2), 'A2/tipMS', round(ff['A2'] / rf['tip_ms'], 3))
