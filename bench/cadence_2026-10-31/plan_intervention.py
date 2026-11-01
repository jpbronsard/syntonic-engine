# Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
# Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
"""Causal test of the acceleration reading: change the gesture's timing, keep its path, and see whether the error
follows. Written 31 Oct 2026 after the held-out confirmation (accel_reading.py, A1 to A3 held) and before any run of
a re-timed plan; go from JP ("Go pour le test causal").

Reading under test. The error the plan causes is the loop's response to the plan's second time derivative (the
motion's acceleration and, through the load at each pose, the load's): the engine's correction integrates the error
over its window, which removes a constant lag, so what remains scales like an acceleration. Re-timing the same path at
a tempo omega (gesture seconds per second) multiplies every second derivative by omega squared. So the part of the
tip error locked to the plan should scale by omega squared, while the part that is not locked to it (noise passed by
the loop) should not move when the window is held fixed. Rivals: an error driven by velocity scales by omega; an
error locked to the pose does not scale.

Plans. Same poses in the same order and the same load at a given pose (the recorded load is kept; inertial loads are
not rescaled, a limit of this bench). The plan is re-timed by a phase map: gesture time phi(t), tempo omega = dphi/dt.
The scored segment is phi in [0, 30 s] in every arm; a 3 s pre-roll runs at the arm's first tempo.
  B    the published plan (omega = 1)
  U70  uniform tempo 0.7 (42.9 s)            U50  uniform tempo 0.5 (60 s)
  EG   tempo from 0.5 to 1 with the rate linear in time (1/period affine: the e-geodesic of the tempo, a reading of
       Dual Geodesy Theorem 1 transposed to the beat, not a statement of the corpus; 40 s)
  PL   tempo from 0.5 to 1 with the period linear in time (the other straight line between the same tempos; 43.3 s)
The faster tempos are excluded: at 1.5 the plan already exceeds the slew limit of RJ4 (checked on the plan alone).

Engines. 'syntonic': the published bench engine (1 kHz, 3 ms frame, adaptive window). 'frozen': the same engine with
its window held, per joint, at the mean window it used over the five baseline runs of the excerpt; the loop is then
fixed, which isolates the plan's effect. Baton arm RJ1-RJ6, steady noise, fresh seeds 11 to 15 (1 to 10 were used
before), excerpts A and B: 2 engines x 5 arms x 5 seeds x 2 excerpts = 100 runs.

Statistics, per excerpt, engine and arm. Tip error vectors resampled at gesture times phi = 0, 0.01, ..., 29.99 s.
Plan-locked part: ebar(phi), the mean over the five seeds. Part not locked to the plan:
    N = (5/4) mean over seeds and phi of |e - ebar|^2.
Scaling coefficient of an arm I against the baseline B of the same engine, with m(phi) = omega_I(phi)^2 ebar_B(phi):
    k = sum ebar_I . m / (sum |m|^2 - (N_B / 5) sum omega_I^4)
(the subtraction removes the noise left in ebar_B; k = 1 when the plan-locked error scales by omega squared; for the
uniform tempos the rivals give k = 1/omega, 2 and 1.43, for velocity, and 1/omega^2, 4 and 2.04, for the pose).
Gate: the baseline arm of the syntonic engine reproduces omega_phase_rates.run (seed 11, both excerpts) exactly
(tip error to 1e-9 mm and every sampled angle), or nothing is judged.

Predictions:
K1  The causal line. With the window frozen, the plan-locked tip error under a slower tempo is the baseline's scaled
    by omega squared: 0.8 <= k <= 1.25 for U50 and U70, in both excerpts (4 cells).
K2  Along the tempo paths. With the window frozen, the same holds when the tempo changes along the run, following the
    local tempo: 0.8 <= k <= 1.25 for EG and PL, in both excerpts (4 cells). (The velocity rival gives about 1.3 to
    1.4 here, near the band's edge; K1 is the sharper line.)
K3  The law's trade. With the window frozen, the part not locked to the plan does not move with the tempo:
    0.9 <= N(U50) / N(B) <= 1.1 in both excerpts; with the adaptive window it falls, because a slower gesture drifts
    less and the law widens the window: N(U50) / N(B) <= 0.9 in both excerpts.
Described, not judged: k and the waveform correlation in every cell (the adaptive engine's k is expected above 1, the
widened window letting more of the acceleration through); P = mean |ebar|^2 - N / 5 and N; tip error by arm, engine
and seed; window used and law by joint; the joints' acceleration gains c_j (least squares of e_j on [1, qddot_j], as
in accel_reading.py) and their ratio to the baseline; EG against PL. A failed line is reported as such. Before this
file was committed it was run only on the baseline plan of seed 1 (exploratory/plan_baseline_explo.py: the identity
gate, the phase maps, and at baseline the gains c_j, the window used and the acceleration part's share of the tip mean
square) and as a smoke test of the gate; no re-timed plan was run.

Run (from the repository root, about 30 minutes):   python bench/cadence_2026-10-31/plan_intervention.py"""
import os, sys, json, time, hashlib
import numpy as np
from scipy.interpolate import CubicSpline
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from cadence_bench import CadenceEngine, Z, TAU_S, W0, ENC_SD, ENC_LSB, TS, DEG, tip
from lag_cancellation import jacobians, canc

DT = 1e-3; LIM = 5 * DEG; PRE = 3.0; SPAN = 30.0; J = 6
SEEDS = [11, 12, 13, 14, 15]
ARMS = ('B', 'U70', 'U50', 'EG', 'PL')


def tempo(arm):
    """Phase map of an arm: run length T (s) and phi, omega, omega_dot as functions of t_rel = t - t0 >= -PRE.
    phi is gesture time relative to t0; the scored segment is phi in [0, SPAN] in every arm."""
    if arm in ('B', 'U70', 'U50'):
        w = {'B': 1.0, 'U70': 0.7, 'U50': 0.5}[arm]
        return SPAN / w, (lambda t: w * t), (lambda t: np.full_like(t, w)), (lambda t: np.zeros_like(t))
    w0, w1 = 0.5, 1.0
    if arm == 'EG':                                                  # rate linear in time: 1/period affine
        T = 2.0 * SPAN / (w0 + w1); a = (w1 - w0) / T
        phi = lambda t: np.where(t < 0, w0 * t, w0 * t + 0.5 * a * t * t)
        om = lambda t: np.where(t < 0, w0, w0 + a * t)
        omd = lambda t: np.where(t < 0, 0.0, a)
        return T, phi, om, omd
    p0, p1 = 1.0 / w0, 1.0 / w1                                      # 'PL': period linear in time
    T = SPAN * (p1 - p0) / np.log(p1 / p0); b = (p1 - p0) / T
    p = lambda t: p0 + b * np.maximum(t, 0.0)
    phi = lambda t: np.where(t < 0, w0 * t, (1.0 / b) * np.log(p(t) / p0))
    om = lambda t: 1.0 / p(t)
    omd = lambda t: np.where(t < 0, 0.0, -b / p(t) ** 2)
    return T, phi, om, omd


def run(name, seed, engine, arm, tau_fix=None):
    """One steady-noise run of the published bench (the arithmetic of omega_phase_rates.run, baton arm RJ1-RJ6) on the
    plan re-timed by the arm's phase map: same poses in the same order, same load at a given pose. engine: 'syntonic'
    (as published: 3 ms frame, adaptive window), 'frozen' (the same engine with its window held at tau_fix per joint),
    'pid'. Samples at 100 Hz over the scored segment."""
    t0 = float(Z[f'{name}_t0']); T, phi, om, omd = tempo(arm)
    tt = np.arange(t0 - PRE, t0 + T, DT)
    ph = tt if arm == 'B' else t0 + phi(tt - t0)
    S = CubicSpline(Z[f'{name}_t30'], Z[f'{name}_q30'][:, :J], axis=0)
    q = S(ph); tau = CubicSpline(Z[f'{name}_t30'], Z[f'{name}_tau30'][:, :J], axis=0)(ph)
    sd = np.full(len(tt), ENC_SD)
    rng = np.random.default_rng(seed)
    tauS, w0 = TAU_S[:J], W0[:J]; k_inv = (1.0 * DEG) / tauS; hb = 0.5 * 0.3 * DEG
    th_g = q[0].copy(); cmd = q[0].copy()
    if engine in ('syntonic', 'frozen'):
        allowed = np.zeros((J, 5), bool); allowed[:, 2] = True
        eng = CadenceEngine(q[0], DT, w0, allowed, np.zeros(J, bool))
        if engine == 'frozen':
            tf = np.asarray(tau_fix, float).copy()
            eng.floor = tf.copy(); eng.tau_max = tf.copy(); eng.tau = tf.copy(); eng.tau_target = tf.copy()
    else:
        kp, ki, kd = (Z[f'{name}_pid'][:J, i] for i in range(3))
        I = np.zeros(J); ef = np.zeros(J); efp = np.zeros(J)
    keep, ks, win = [], [], []
    for k in range(len(tt)):
        load = tau[k]
        th_out = th_g - load * k_inv - hb * np.tanh(load / (0.02 * tauS))
        a = np.round((th_out + rng.normal(0.0, 1.0, J) * sd[k]) / ENC_LSB) * ENC_LSB
        b = np.round((th_out + rng.normal(0.0, 1.0, J) * sd[k]) / ENC_LSB) * ENC_LSB
        if engine != 'pid':
            cmd = eng.step(q[k], a, b)
        else:
            e = q[k] - 0.5 * (a + b)
            ef += (DT / 0.005) * (e - ef)
            raw = kp * e + ki * I + kd * (ef - efp) / DT
            I += np.where(np.abs(raw) >= LIM, 0.0, e * DT); efp = ef.copy()
            if k % 3 == 0: cmd = q[k] + np.clip(raw, -LIM, LIM)
        w_eff = w0 * np.maximum(0.05, 1.0 - np.abs(load) / tauS)
        th_g = th_g + DT * np.clip((cmd - th_g) / TS, -w_eff, w_eff)
        if k % 10 == 0 and tt[k] >= t0 - 1e-9:
            keep.append((th_out[:6].copy(), q[k, :6].copy())); ks.append(k)
            win.append(np.vstack([eng.tau[:6], eng.law[:6]]) if engine != 'pid' else np.full((2, 6), np.nan))
    th = np.array([x for x, _ in keep]); qp = np.array([y for _, y in keep]); ks = np.array(ks)
    e = np.array([tip(x) - tip(y) for x, y in keep])
    trel = tt[ks] - t0; w, wd = om(trel), omd(trel)
    if arm == 'B': w, wd = np.ones_like(trel), np.zeros_like(trel)
    qdd = S(ph[ks], 2) * w[:, None] ** 2 + S(ph[ks], 1) * wd[:, None]
    qd = S(ph[ks], 1) * w[:, None]
    win = np.array(win)
    return dict(phi=ph[ks] - t0, omega=w, tip_mm=float(1e3 * np.sqrt(np.mean(np.sum(e ** 2, 1)))), tip_ms=float(np.mean(np.sum(e ** 2, 1))),
                th=th, q=qp, qdd=qdd, qd=qd, tau_used=win[:, 0], law=win[:, 1], etip=e, T=T)


def fit(r, Js):
    """Per joint, least squares of e_j on [1, qdd_j]; gain c_j, R2_j; the acceleration part at the tip."""
    e = r['th'] - r['q']; qdd = r['qdd']; n = len(e); c = np.zeros(6); R2 = np.zeros(6); eacc = np.zeros_like(e)
    for j in range(6):
        X = np.column_stack([np.ones(n), qdd[:, j]]); co, *_ = np.linalg.lstsq(X, e[:, j], rcond=None)
        c[j] = co[1]; R2[j] = 1.0 - np.var(e[:, j] - X @ co) / np.var(e[:, j]); eacc[:, j] = c[j] * qdd[:, j]
    ix = list(range(n))
    return dict(c=c, R2=R2, A2=float(np.mean(np.sum(np.einsum('nij,nj->ni', Js, eacc) ** 2, 1))),
                C_err=canc(Js, e, ix), C_acc_part=canc(Js, eacc, ix), C_plan_acc=canc(Js, qdd, ix))


def A2_pred(c, r, Js):
    """Tip mean square of the acceleration part predicted from gains c and the run's planned accelerations."""
    return float(np.mean(np.sum(np.einsum('nij,nj->ni', Js, c[None, :] * r['qdd']) ** 2, 1)))


GRID = np.arange(3000) * 0.01


def on_grid(r, key='etip'):
    x = r[key]
    return np.column_stack([np.interp(GRID, r['phi'], x[:, i]) for i in range(x.shape[1])]) if x.ndim == 2 \
        else np.interp(GRID, r['phi'], x)


def cell_stats(G):
    G = np.asarray(G); eb = G.mean(0)
    N = 1.25 * float(np.mean(np.sum((G - eb[None]) ** 2, 2)))
    return eb, N, float(np.mean(np.sum(eb ** 2, 1))) - N / 5.0


if __name__ == '__main__':
    sha = hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest()
    from omega_phase_rates import run as run_ref
    smoke = '--smoke' in sys.argv
    t_start = time.time(); out = dict(script_sha256=sha, gate=[], tau_fix={}, runs={}, cells={})
    gate_ok = True
    for name in ('A', 'B'):
        grids, oms = {}, {}
        for engine in ('syntonic', 'frozen'):
            if engine == 'frozen':
                out['tau_fix'][name] = np.mean([out['runs'][f'{name} syntonic B {s_}']['window_used_mean'] for s_ in SEEDS], 0).tolist()
            for arm in ARMS:
                for seed in SEEDS:
                    r = run(name, seed, engine, arm, tau_fix=out['tau_fix'].get(name))
                    if engine == 'syntonic' and arm == 'B' and seed == SEEDS[0]:
                        ref = run_ref(name, 'steady', seed, 'syntonic', 6)
                        d = dict(excerpt=name, tip_diff=abs(r['tip_mm'] - ref['tip_mm']), th_diff=float(np.max(np.abs(r['th'] - ref['th']))))
                        out['gate'].append(d); gate_ok &= d['tip_diff'] <= 1e-9 and d['th_diff'] == 0.0
                        print('gate', d, flush=True)
                        if smoke: sys.exit(0)
                    Js = jacobians(r['q']); f = fit(r, Js)
                    out['runs'][f'{name} {engine} {arm} {seed}'] = dict(
                        tip_mm=r['tip_mm'], c=f['c'].tolist(), R2=f['R2'].tolist(), A2_over_tipMS=f['A2'] / r['tip_ms'],
                        C_err=f['C_err'], C_acc_part=f['C_acc_part'], C_plan_acc=f['C_plan_acc'],
                        window_used_mean=r['tau_used'].mean(0).tolist(), law_median=np.nanmedian(r['law'], 0).tolist())
                    grids.setdefault((engine, arm), []).append(on_grid(r)); oms[arm] = on_grid(r, 'omega')
                    print(name, engine, arm, seed, 'tip %.4f mm' % r['tip_mm'], 'R2', np.round(np.median(f['R2']), 2),
                          round(time.time() - t_start), 's', flush=True)
        for engine in ('syntonic', 'frozen'):
            ebB, NB, PB = cell_stats(grids[(engine, 'B')])
            for arm in ARMS:
                eb, N, P = cell_stats(grids[(engine, arm)]); w2 = oms[arm] ** 2
                m = w2[:, None] * ebB
                k = float(np.sum(eb * m) / (np.sum(m ** 2) - (NB / 5.0) * np.sum(w2 ** 2)))
                corr = float(np.corrcoef(eb.ravel(), m.ravel())[0, 1])
                tips = [out['runs'][f'{name} {engine} {arm} {s_}']['tip_mm'] for s_ in SEEDS]
                cr = np.median([np.array(out['runs'][f'{name} {engine} {arm} {s_}']['c']) / np.array(out['runs'][f'{name} {engine} B {s_}']['c'])
                                for s_ in SEEDS], 0)
                out['cells'][f'{name} {engine} {arm}'] = dict(k=k, corr=corr, N=N, P=P, N_over_NB=N / NB, tip_mm_median=float(np.median(tips)),
                                                             c_ratio_median_by_joint=cr.tolist())
                print(name, engine, arm, 'k %.3f corr %.3f N/NB %.3f P %.3g N %.3g tip %.4f' % (k, corr, N / NB, P, N, np.median(tips)),
                      'c ratio', np.round(cr, 2), flush=True)
    cl = out['cells']; inb = lambda x: 0.8 <= x <= 1.25
    out['gate_ok'] = bool(gate_ok)
    out['verdict'] = dict(
        K1=bool(gate_ok and all(inb(cl[f'{n} frozen {a}']['k']) for n in 'AB' for a in ('U50', 'U70'))),
        K1_k={f'{n} {a}': cl[f'{n} frozen {a}']['k'] for n in 'AB' for a in ('U50', 'U70')},
        K2=bool(gate_ok and all(inb(cl[f'{n} frozen {a}']['k']) for n in 'AB' for a in ('EG', 'PL'))),
        K2_k={f'{n} {a}': cl[f'{n} frozen {a}']['k'] for n in 'AB' for a in ('EG', 'PL')},
        K3=bool(gate_ok and all(0.9 <= cl[f'{n} frozen U50']['N_over_NB'] <= 1.1 and cl[f'{n} syntonic U50']['N_over_NB'] <= 0.9 for n in 'AB')),
        K3_ratios={f'{n} {e}': cl[f'{n} {e} U50']['N_over_NB'] for n in 'AB' for e in ('frozen', 'syntonic')})
    out['minutes'] = round((time.time() - t_start) / 60, 1)
    json.dump(out, open(os.path.join(HERE, 'results', 'plan_intervention.json'), 'w'), indent=1)
    print('verdict', out['verdict'])
