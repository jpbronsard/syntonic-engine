# Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
# Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
"""The acceleration reading, on fresh seeds. Written 31 Oct 2026 after the exploratory pass exploratory/acc_explo.py
(seeds 1 and 2), before any run of this file; go from JP ("Peut-on verifier si cette lecture tient ?").

Reading under test. The syntonic engine's joint error is largely the planned acceleration seen through the servo's
dynamics; that part inherits the gesture's own cancellation of joint accelerations at the baton tip, which is why the
engine's joint errors cancel there and the PID's do not.

Same plant, engines, sampling and definitions as lag_cancellation.py and omega_phase_rates.py (the published bench
at 1 kHz with a 3 ms frame, baton arm, 100 Hz samples over the scored 30 s), with the planned acceleration qddot from
the plan's spline in place of the velocity. Held-out seeds 6 to 10 (the exploratory pass used 1 and 2), excerpts A
and B, steady noise: 10 runs. Per joint and run: least squares of e_j on [1, qddot_j]; R2_j; acceleration part
e_acc_j = (fitted qddot coefficient) x qddot_j. Cancellation C(x) as before, over the whole run.

Predictions:
A1  The planned acceleration explains the syntonic engine's joint error: median R2 over 6 joints x 10 runs >= 0.35,
    and the PID's median R2 at most half of the syntonic engine's.
A2  That part cancels at the tip like the gesture's own accelerations: median over the runs of
    |C(e_acc) - C(qddot)| <= 0.10 for the syntonic engine.
A3  The cancellation replicates on fresh seeds: the syntonic engine's run-level C(e) is at least 0.2 below the PID's
    in every run.
Described: R2 by joint, C(e), C(e_acc), C(qddot), and the same for the PID. A failed line is reported as such. Before
this file was committed it was run only as a smoke test (imports and one 2 s run; no statistic printed).

Run (from the repository root, about 4 minutes):   python bench/cadence_2026-10-31/accel_reading.py
"""
import os, sys, json, time, hashlib
import numpy as np
from scipy.interpolate import CubicSpline
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from omega_phase_rates import run, Z
from lag_cancellation import jacobians, canc

SEEDS_HELD_OUT = [6, 7, 8, 9, 10]


def analyse(r, qdd, Js):
    e = r['th'] - r['q']; n = len(e); R2 = np.zeros(6); eacc = np.zeros_like(e)
    for j in range(6):
        X = np.column_stack([np.ones(n), qdd[:, j]]); c, *_ = np.linalg.lstsq(X, e[:, j], rcond=None)
        R2[j] = 1.0 - np.var(e[:, j] - X @ c) / np.var(e[:, j]); eacc[:, j] = c[1] * qdd[:, j]
    ix = list(range(n))
    return dict(R2=R2.tolist(), C_err=canc(Js, e, ix), C_acc_part=canc(Js, eacc, ix), C_plan_acc=canc(Js, qdd, ix))


if __name__ == '__main__':
    sha = hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest()
    if '--smoke' in sys.argv:
        r = run('A', 'steady', 6, 'syntonic', 6, t_end=2.0); print('smoke ok', r['th'].shape); sys.exit(0)
    t_start = time.time(); out = dict(script_sha256=sha, runs={})
    for name in ('A', 'B'):
        for seed in SEEDS_HELD_OUT:
            rs = run(name, 'steady', seed, 'syntonic', 6); rp = run(name, 'steady', seed, 'pid', 6)
            t0 = float(Z[f'{name}_t0']); ts = t0 + 0.01 * np.arange(len(rs['q']))
            qdd = CubicSpline(Z[f'{name}_t30'], Z[f'{name}_q30'][:, :6], axis=0)(ts, 2)
            Js = jacobians(rs['q'])
            s, p = analyse(rs, qdd, Js), analyse(rp, qdd, Js)
            out['runs'][f'{name} steady {seed}'] = dict(syntonic=s, pid=p)
            print(name, seed, 'syn R2', np.round(s['R2'], 2), 'C %.3f Cacc %.3f Cplan %.3f | pid R2' % (s['C_err'], s['C_acc_part'], s['C_plan_acc']),
                  np.round(p['R2'], 2), 'C %.3f' % p['C_err'], round(time.time() - t_start), 's', flush=True)
    R = out['runs'].values()
    R2s = float(np.median([x for v in R for x in v['syntonic']['R2']])); R2p = float(np.median([x for v in R for x in v['pid']['R2']]))
    d2 = float(np.median([abs(v['syntonic']['C_acc_part'] - v['syntonic']['C_plan_acc']) for v in R]))
    a3 = [v['pid']['C_err'] - v['syntonic']['C_err'] for v in R]
    out['verdict'] = dict(A1=bool(R2s >= 0.35 and R2p <= 0.5 * R2s), A1_R2_syn=R2s, A1_R2_pid=R2p,
                          A2=bool(d2 <= 0.10), A2_median_absdiff=d2, A3=bool(min(a3) >= 0.2), A3_min_gap=float(min(a3)))
    out['minutes'] = round((time.time() - t_start) / 60, 1)
    json.dump(out, open(os.path.join(HERE, 'results', 'accel_reading.json'), 'w'), indent=1)
    print('verdict', out['verdict'])
