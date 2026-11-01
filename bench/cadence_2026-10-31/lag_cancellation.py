# Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
# Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
"""Why do the syntonic engine's joint errors cancel at the baton tip? Written 31 Oct 2026 after the Theorem Omega test
(W2 failed; described: cancellation ratio 0.58 to 0.72 for the syntonic engine, 0.94 to 0.98 for the tuned PID),
before any run of this file; go from JP.

Hypothesis. The syntonic engine's error is mostly a lag along the planned motion, e_j(t) ~ -L_j qdot_j(t), with lags
close across the joints; its tip error is then the planned gesture slightly late, and it inherits the gesture's own
kinematic cancellation at the tip (shoulder and elbow working against each other to hold the baton's path). The PID
adds errors of each joint's own, which reach the tip as independent errors.

Same runs, plant, engines and scores as omega_phase_rates.py (the published bench at 1 kHz with a 3 ms frame,
baton arm, 100 Hz samples over the scored 30 s; excerpts A and B, steady noise seeds 1 to 5 judged, noise step seeds 1
to 5 described). Planned velocity qdot from the same spline as the plan. Tip Jacobian J(q) at the planned pose.
- Cancellation of a joint-space signal x over samples W: C_W(x) = sum |J x|^2 / sum_t sum_j |J_j x_j|^2.
- The gesture's own cancellation: C_vel = C(qdot). The error's: C_err = C(e), e = actual - planned angle.
- Lag fit, per joint and run: least squares of e_j on [1, qdot_j]; L_j = minus the qdot coefficient; R2_j the share
  of e_j's variance the fit explains. Lag part e_lag_j = -L_j qdot_j (the offset is not lag and stays in the rest).
Gate: every run reproduces the tip errors recorded in results/omega_phase_rates.json to 1e-9 mm (syntonic and PID),
or nothing is judged.

Predictions (10 steady runs):
L1  The error's cancellation follows the gesture's: within a run, the Spearman correlation over the 30 one-second
    windows between C_err and C_vel has a median of at least 0.5 for the syntonic engine, and the PID's median is at
    least 0.2 lower.
L2  The lag explains at least half of the syntonic engine's joint error: median R2 over 6 joints x 10 runs >= 0.5,
    and the PID's median R2 is lower.
L3  The lag part alone reproduces the syntonic engine's cancellation: median over the runs of |C(e_lag) - C_err|,
    computed over the whole run, <= 0.05.
Described, not judged: the lags L_j (ms) and their spread across joints; C_vel per run; the same for the noise step.
If L1 and L3 hold, the cancellation is the gesture's own, carried by a coherent lag; if L1 holds and L3 fails, it
follows the gesture by another path; if L1 fails, this lead is closed. A failed line is reported as such. Before this
file was committed it was run only as a smoke test (gate path on one run; no statistic printed).

Run (from the repository root, about 6 minutes):   python bench/cadence_2026-10-31/lag_cancellation.py
"""
import os, sys, json, time, hashlib
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.stats import spearmanr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from omega_phase_rates import run, jacobian, Z, SEEDS


def jacobians(q):
    return np.array([jacobian(qq) for qq in q])                               # (n, 3, 6)


def canc(Js, x, idx):
    contrib = Js[idx] * x[idx][:, None, :]                                     # (m, 3, 6)
    return float(np.sum(contrib.sum(2) ** 2) / np.sum(contrib ** 2))


def analyse(r, qd, Js):
    e = r['th'] - r['q']; n = len(e)
    L = np.zeros(6); R2 = np.zeros(6); elag = np.zeros_like(e)
    for j in range(6):
        X = np.column_stack([np.ones(n), qd[:, j]]); coef, *_ = np.linalg.lstsq(X, e[:, j], rcond=None)
        fit = X @ coef; L[j] = -coef[1]
        R2[j] = 1.0 - np.var(e[:, j] - fit) / np.var(e[:, j]); elag[:, j] = -L[j] * qd[:, j]
    W = 100; nw = n // W; win = [range(w * W, (w + 1) * W) for w in range(nw)]
    Cerr_w = [canc(Js, e, list(ix)) for ix in win]; Cvel_w = [canc(Js, qd, list(ix)) for ix in win]
    allix = list(range(n))
    return dict(lag_ms=(1e3 * L).tolist(), R2=R2.tolist(), C_err=canc(Js, e, allix), C_vel=canc(Js, qd, allix),
                C_lag=canc(Js, elag, allix), rho_err_vel=float(spearmanr(Cerr_w, Cvel_w).correlation))


if __name__ == '__main__':
    sha = hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest()
    ref = json.load(open(os.path.join(HERE, 'results', 'omega_phase_rates.json')))['runs']
    smoke = '--smoke' in sys.argv
    t_start = time.time(); out = dict(script_sha256=sha, runs={}); g = 0.0
    for name in ('A', 'B'):
        for cond in ('steady', 'step'):
            for seed in SEEDS:
                key = f'{name} {cond} {seed}'
                rs = run(name, cond, seed, 'syntonic', 6); rp = run(name, cond, seed, 'pid', 6)
                g = max(g, abs(rs['tip_mm'] - ref[key]['tip_mm']), abs(rp['tip_mm'] - ref[key]['tip_mm_pid']))
                if smoke: print('smoke ok, gate diff', g); sys.exit(0)
                t0 = float(Z[f'{name}_t0']); ts = t0 + 0.01 * np.arange(len(rs['q']))
                qd = CubicSpline(Z[f'{name}_t30'], Z[f'{name}_q30'][:, :6], axis=0)(ts, 1)
                Js = jacobians(rs['q'])
                out['runs'][key] = dict(syntonic=analyse(rs, qd, Js), pid=analyse(rp, qd, Js))
                s, p = out['runs'][key]['syntonic'], out['runs'][key]['pid']
                print(key, 'C_vel %.3f | syn C_err %.3f C_lag %.3f rho %.3f R2 %.2f | pid C_err %.3f rho %.3f R2 %.2f'
                      % (s['C_vel'], s['C_err'], s['C_lag'], s['rho_err_vel'], np.median(s['R2']), p['C_err'], p['rho_err_vel'],
                         np.median(p['R2'])), 'gate', g, round(time.time() - t_start), 's', flush=True)
    out['gate_max_diff_mm'] = g; out['gate_ok'] = g <= 1e-9
    st = [v for k, v in out['runs'].items() if ' steady ' in k]
    rs_ = [v['syntonic']['rho_err_vel'] for v in st]; rp_ = [v['pid']['rho_err_vel'] for v in st]
    R2s = np.median([x for v in st for x in v['syntonic']['R2']]); R2p = np.median([x for v in st for x in v['pid']['R2']])
    d3 = np.median([abs(v['syntonic']['C_lag'] - v['syntonic']['C_err']) for v in st])
    ok = out['gate_ok']
    out['verdict'] = dict(L1=bool(ok and np.median(rs_) >= 0.5 and np.median(rp_) <= np.median(rs_) - 0.2),
                          L1_median_syn=float(np.median(rs_)), L1_median_pid=float(np.median(rp_)),
                          L2=bool(ok and R2s >= 0.5 and R2p < R2s), L2_median_R2_syn=float(R2s), L2_median_R2_pid=float(R2p),
                          L3=bool(ok and d3 <= 0.05), L3_median_absdiff=float(d3))
    out['minutes'] = round((time.time() - t_start) / 60, 1)
    json.dump(out, open(os.path.join(HERE, 'results', 'lag_cancellation.json'), 'w'), indent=1)
    print('gate', g, out['gate_ok']); print('verdict', out['verdict'])
