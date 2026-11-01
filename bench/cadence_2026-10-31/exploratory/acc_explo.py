"""Exploratory (31 Oct 2026, after the lag test failed L1 to L3; no prediction): is the syntonic engine's joint error
an acceleration term of the planned gesture, and does that term carry the tip cancellation? Steady noise, 4 runs."""
import sys, numpy as np
sys.path.insert(0, sys.argv[1])
from omega_phase_rates import run, Z
from lag_cancellation import jacobians, canc
from scipy.interpolate import CubicSpline
for name in ('A', 'B'):
    for seed in (1, 2):
        for eng in ('syntonic', 'pid'):
            r = run(name, 'steady', seed, eng, 6); n = len(r['q'])
            t0 = float(Z[f'{name}_t0']); ts = t0 + 0.01 * np.arange(n)
            sp = CubicSpline(Z[f'{name}_t30'], Z[f'{name}_q30'][:, :6], axis=0); qdd = sp(ts, 2)
            e = r['th'] - r['q']; Js = jacobians(r['q']); R2 = []; eacc = np.zeros_like(e)
            for j in range(6):
                X = np.column_stack([np.ones(n), qdd[:, j]]); c, *_ = np.linalg.lstsq(X, e[:, j], rcond=None)
                R2.append(1 - np.var(e[:, j] - X @ c) / np.var(e[:, j])); eacc[:, j] = c[1] * qdd[:, j]
            ix = list(range(n))
            print(name, seed, eng, 'R2(e ~ qdd) by joint', np.round(R2, 2), '| C_err %.3f C_acc_part %.3f C_plan_acc %.3f'
                  % (canc(Js, e, ix), canc(Js, eacc, ix), canc(Js, qdd, ix)), flush=True)
