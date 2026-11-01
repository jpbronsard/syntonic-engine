# Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
# Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
"""The correction limit from the authority fact. Written 31 Oct 2026 after the exploratory passes of 024b90d, before
any run of this file; go from JP ("Go", for the option recommended: the limit derived from the authority fact).

Background. The engine clips its accumulated correction at a fixed 5 degrees, one of the design constants the main
statement of 26 October keeps identical on every machine. The correction must carry the servo's lag, about the servo
time constant times the joint speed; at the published tempo RJ4 and RJ6 need more than 5 degrees and sit at the limit
(noise-free, frozen window: 0.6 to 0.7 % of the time in excerpt A, 3.7 % in B; never at tempo 0.7 or 0.5), and that is
the extra tip error found by the causal test. The option tested: limit = 'authority', per joint servo_tau x w0, the
largest lag the servo itself can produce, from the two data-sheet values the engine already takes (RJ1 and RJ2 7.2
degrees, RJ3 and RJ5 15, RJ4 12.5, RJ6 25). Nothing new is fitted. The option was added to the Python engine in the
same commit as this file; the default stays 5 degrees until a verdict.

Setup: the plant, plans and arms of plan_intervention.py (published tempo B; uniform tempos U70 and U50), baton arm,
steady noise, fresh seeds 16 to 20 (1 to 15 were used before), excerpts A and B. The adaptive engine (as published)
with the 5 degree limit and with the authority limit, paired on the same seeds; the frozen-window engine (windows of
results/plan_intervention.json) with the authority limit, for k. Gate: with the 5 degree limit this runner reproduces
plan_intervention.run exactly (seed 16, excerpt A, arm B, adaptive), or nothing is judged.

Predictions (adaptive engine):
L1  The authority limit does not bind: the share of time a joint's correction sits at its limit is at most 0.1 % for
    every joint in every run (3 arms x 5 seeds x 2 excerpts).
L2  The tip error falls where the 5 degree limit bound: at the published tempo the tip error with the authority limit
    is below the 5 degree run of the same seed in every run of both excerpts, and at least 20 % below it in every run
    of excerpt B.
L3  Nothing changes where it did not bind: at tempos 0.7 and 0.5 the two tip errors of the same seed differ by at most
    1 % in every run.
Described, not judged: the share of time at the 5 degree limit by joint and tempo; the tip errors; with the frozen
window and the authority limit, k as in plan_intervention.py (exploratory noise-free values with a 30 degree limit:
0.87 and 0.81 in A, 0.82 and 0.76 in B); a blocked joint (RJ4 held still from 10 to 12 s, published tempo, seed 16):
the largest command offset from the plan and the time to come back within 0.2 degree after release, for both limits,
since a wider limit lets the command move further from the plan when a joint cannot follow. A failed line is reported
as such. Before this file was committed it was run only as a smoke test of the gate (no statistic printed).

Run (from the repository root, about 25 minutes):   python bench/cadence_2026-10-31/authority_limit.py
"""
import os, sys, json, time, hashlib
import numpy as np
from scipy.interpolate import CubicSpline
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import plan_intervention as P

SEEDS = [16, 17, 18, 19, 20]; ARMS = ('B', 'U70', 'U50'); J = 6; LIM5 = np.radians(5.0)


def run(name, seed, arm, limit, tau_fix=None, block=False):
    t0 = float(P.Z[f'{name}_t0']); T, phi, om, omd = P.tempo(arm)
    tt = np.arange(t0 - P.PRE, t0 + T, P.DT)
    ph = tt if arm == 'B' else t0 + phi(tt - t0)
    q = CubicSpline(P.Z[f'{name}_t30'], P.Z[f'{name}_q30'][:, :J], axis=0)(ph)
    tau = CubicSpline(P.Z[f'{name}_t30'], P.Z[f'{name}_tau30'][:, :J], axis=0)(ph)
    sd = np.full(len(tt), P.ENC_SD); rng = np.random.default_rng(seed)
    tauS, w0 = P.TAU_S[:J], P.W0[:J]; k_inv = (1.0 * P.DEG) / tauS; hb = 0.5 * 0.3 * P.DEG
    th_g = q[0].copy()
    allowed = np.zeros((J, 5), bool); allowed[:, 2] = True
    eng = P.CadenceEngine(q[0], P.DT, w0, allowed, np.zeros(J, bool), limit=limit)
    if tau_fix is not None:
        tf = np.asarray(tau_fix, float); eng.floor = tf.copy(); eng.tau_max = tf.copy(); eng.tau = tf.copy(); eng.tau_target = tf.copy()
    keep, ks = [], []; at_lim = np.zeros(J); n_sc = 0; off_max = 0.0; e3 = []
    for k in range(len(tt)):
        load = tau[k]
        th_out = th_g - load * k_inv - hb * np.tanh(load / (0.02 * tauS))
        a = np.round((th_out + rng.normal(0.0, 1.0, J) * sd[k]) / P.ENC_LSB) * P.ENC_LSB
        b = np.round((th_out + rng.normal(0.0, 1.0, J) * sd[k]) / P.ENC_LSB) * P.ENC_LSB
        cmd = eng.step(q[k], a, b)
        w_eff = w0 * np.maximum(0.05, 1.0 - np.abs(load) / tauS)
        new = th_g + P.DT * np.clip((cmd - th_g) / P.TS, -w_eff, w_eff)
        if block and t0 + 10.0 <= tt[k] < t0 + 12.0:
            new[3] = th_g[3]; off_max = max(off_max, abs(cmd[3] - q[k, 3]))
        th_g = new
        if tt[k] >= t0 - 1e-9:
            at_lim += np.abs(eng.u) >= eng.limit - 1e-12; n_sc += 1
            if block and t0 + 12.0 <= tt[k] < t0 + 15.0: e3.append(abs(q[k, 3] - th_out[3]))
        if k % 10 == 0 and tt[k] >= t0 - 1e-9:
            keep.append((th_out[:6].copy(), q[k, :6].copy())); ks.append(k)
    ks = np.array(ks); e = np.array([P.tip(x) - P.tip(y) for x, y in keep])
    trel = tt[ks] - t0; w = np.ones_like(trel) if arm == 'B' else om(trel)
    out = dict(etip=e, phi=ph[ks] - t0, omega=w, tip_mm=float(1e3 * np.sqrt(np.mean(np.sum(e ** 2, 1)))), at_limit=(at_lim / n_sc).tolist())
    if block:
        e3 = np.array(e3); ok = e3 < np.radians(0.2); rec = np.nan
        for i in range(len(e3) - 200):
            if ok[i:i + 200].all(): rec = i * P.DT; break
        out.update(block_cmd_offset_max_deg=float(np.degrees(off_max)), block_recovery_s=float(rec))
    return out


if __name__ == '__main__':
    sha = hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest()
    ref = P.run('A', SEEDS[0], 'syntonic', 'B'); mine = run('A', SEEDS[0], 'B', LIM5)
    gate = dict(tip_diff=abs(ref['tip_mm'] - mine['tip_mm']), th_equal=bool(np.array_equal(ref['etip'], mine['etip'])))
    gate_ok = gate['tip_diff'] <= 1e-9 and gate['th_equal']; print('gate', gate, flush=True)
    if '--smoke' in sys.argv: sys.exit(0)
    t_start = time.time(); tfs = json.load(open(os.path.join(HERE, 'results', 'plan_intervention.json')))['tau_fix']
    out = dict(script_sha256=sha, gate=gate, gate_ok=bool(gate_ok), runs={}, cells={}, blocked={})
    for name in ('A', 'B'):
        for arm in ARMS:
            for seed in SEEDS:
                r5 = run(name, seed, arm, LIM5); ra = run(name, seed, arm, 'authority')
                out['runs'][f'{name} {arm} {seed}'] = dict(tip_5=r5['tip_mm'], tip_auth=ra['tip_mm'], at_limit_5=r5['at_limit'], at_limit_auth=ra['at_limit'])
                print(name, arm, seed, 'tip 5deg %.4f auth %.4f' % (r5['tip_mm'], ra['tip_mm']), 'at 5deg', np.round(r5['at_limit'], 4),
                      'at auth', np.round(ra['at_limit'], 4), round(time.time() - t_start), 's', flush=True)
        G = {arm: [P.on_grid(run(name, s, arm, 'authority', tau_fix=tfs[name])) for s in SEEDS] for arm in ARMS}
        W = {arm: P.on_grid(dict(etip=np.zeros((2, 3)), phi=np.array([0.0, 30.0]), omega=np.full(2, dict(B=1.0, U70=0.7, U50=0.5)[arm])), 'omega') for arm in ARMS}
        ebB, NB, PB = P.cell_stats(G['B'])
        for arm in ('U70', 'U50'):
            eb, N, Pp = P.cell_stats(G[arm]); w2 = W[arm] ** 2; m = w2[:, None] * ebB
            out['cells'][f'{name} frozen authority {arm}'] = dict(k=float(np.sum(eb * m) / (np.sum(m ** 2) - (NB / 5.0) * np.sum(w2 ** 2))))
            print(name, 'frozen authority', arm, 'k %.3f' % out['cells'][f'{name} frozen authority {arm}']['k'], flush=True)
        for lim, lab in ((LIM5, '5deg'), ('authority', 'authority')):
            rb = run(name, SEEDS[0], 'B', lim, block=True)
            out['blocked'][f'{name} {lab}'] = dict(cmd_offset_max_deg=rb['block_cmd_offset_max_deg'], recovery_s=rb['block_recovery_s'], tip_mm=rb['tip_mm'])
            print(name, 'blocked RJ4', lab, out['blocked'][f'{name} {lab}'], flush=True)
    R = out['runs']
    l1 = max(max(v['at_limit_auth']) for v in R.values())
    l2 = all(v['tip_auth'] < v['tip_5'] for k, v in R.items() if ' B ' in k) and \
        all(v['tip_auth'] <= 0.8 * v['tip_5'] for k, v in R.items() if k.startswith('B B '))
    l3 = max(abs(v['tip_auth'] - v['tip_5']) / v['tip_5'] for k, v in R.items() if (' U70 ' in k or ' U50 ' in k))
    out['verdict'] = dict(L1=bool(gate_ok and l1 <= 0.001), L1_max_share=float(l1), L2=bool(gate_ok and l2),
                          L2_ratios={k: v['tip_auth'] / v['tip_5'] for k, v in R.items() if ' B ' in k},
                          L3=bool(gate_ok and l3 <= 0.01), L3_max_rel_diff=float(l3))
    out['minutes'] = round((time.time() - t_start) / 60, 1)
    json.dump(out, open(os.path.join(HERE, 'results', 'authority_limit.json'), 'w'), indent=1)
    print('verdict', {k: v for k, v in out['verdict'].items() if k != 'L2_ratios'})
