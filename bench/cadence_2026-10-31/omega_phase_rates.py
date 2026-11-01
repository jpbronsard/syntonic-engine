# Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
# Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
"""Theorem Omega on the joint bench: do the six joints' phase rates move together, and does their agreement go with
the cancellation of joint errors at the baton tip? Written 31 Oct 2026, before any run of this file and before
looking at any per-joint window series; go from JP.

Background. Syntony Principle V5.1, section 18 (Theorem Omega): inside a declared synchronization block, the realized
phase rates omega_a = 1 / T*_a agree, and a connected Laplacian coupling contracts their disagreement; the agreement
is a hypothesis of the theorem, not a conclusion. The six joints of the baton arm are a natural block: their errors
combine at the tip. The README records an open fact of the published bench: the syntonic engine's joint errors are
about 28 % larger than the tuned PID's under steady noise, yet its tip error is within 0 to 7 %, so its joint errors
largely cancel at the tip, and why is not known. Reading tested here (a hypothesis): the cancellation goes with the
agreement of the joints' phase rates.

Setup: the published bench (examples/servo_bench.py: 1 kHz loop, 3 ms frame, the bench engine with identification
'raw', and the PID tuned joint by joint), run through the plant loop of cadence_bench.py; baton arm (RJ1-RJ6).
Phase rate of joint j: omega_j(t) = 1 / tau_law_j(t), tau_law the law's window before it is clipped to the floor
(the T* of Theorem Omega; the clipped window would agree trivially whenever the joints sit on the same floor).
Sampled at 100 Hz over the scored 30 s. Joint error e_j = actual - planned angle; tip Jacobian J(q) by central
differences of arm_fk.tip at the planned pose. Cancellation ratio over a window W of samples:
    C_W = sum_t |J e|^2 / sum_t sum_j |J_j e_j|^2     (1 if the joint contributions to the tip were uncorrelated,
                                                       below 1 when they cancel, above 1 when they add).
Dispersion of the phase rates in a window: D_W = mean over t of the standard deviation across the six joints of
log omega_j(t). Windows of 1 s (30 per run). Should a joint's law not yet be valid at a sample, the sample is filled from its
nearest valid ones, and the count is recorded.
Runs: excerpts A and B; steady noise, seeds 1 to 5 (10 runs, judged); noise step x10 from 10 to 20 s, seeds 1 to 5
(described only, since the step drives every joint's noise at once).

Gate: with 12 joints and the published seeds (5 steady, 7 step), the syntonic and PID tip errors reproduce
examples/servo_bench.py in its four conditions to 1e-9 mm. If the gate fails, nothing is judged.

Predictions (steady runs):
W1  Co-movement (the weak one; the shared gesture alone may give it). The mean pairwise correlation across the six
    joints of log omega_j (100 Hz, demeaned per joint) exceeds the same statistic computed with joint j's series
    shifted circularly by 5 j seconds (j = 0..5), in at least 9 of the 10 runs.
W2  The reading. Within a run, the Spearman correlation across its 30 windows between D_W and C_W is positive (more
    dispersion of the phase rates, less cancellation): median over the 10 runs at least 0.2, and positive in at least
    8 of the 10.
Described, not judged: the pooled cancellation ratio C over the run for the syntonic engine and for the PID; W2 with
the window's mean tip speed partialled out; the same statistics under the noise step.
A failed line is reported as such. Before this file was committed it was run only as a smoke test (2 s of excerpt A,
gate path only; no statistic printed).

Run (from the repository root, about 6 minutes):   python bench/cadence_2026-10-31/omega_phase_rates.py
"""
import os, sys, json, time, hashlib
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.signal import butter, sosfiltfilt
from scipy.stats import spearmanr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from cadence_bench import CadenceEngine, Z, TAU_S, W0, ENC_SD, ENC_LSB, TS, DEG, tip

DT = 1e-3; LIM = 5 * DEG; SEEDS = [1, 2, 3, 4, 5]


def run(name, cond, seed, engine, J, t_end=30.0):
    """One run of the published bench (servo_bench.py arithmetic). Returns the tip error and, at 100 Hz over the
    scored window, actual and planned angles of RJ1-RJ6 and, for the syntonic engine, tau_law of RJ1-RJ6."""
    t0 = float(Z[f'{name}_t0']); tt = np.arange(t0 - 3.0, t0 + t_end, DT)
    q = CubicSpline(Z[f'{name}_t30'], Z[f'{name}_q30'][:, :J], axis=0)(tt)
    tau = CubicSpline(Z[f'{name}_t30'], Z[f'{name}_tau30'][:, :J], axis=0)(tt)
    sd = np.full(len(tt), ENC_SD)
    if cond == 'step': sd[(tt >= t0 + 10) & (tt < t0 + 20)] *= np.sqrt(10.0)
    rng = np.random.default_rng(seed)
    tauS, w0 = TAU_S[:J], W0[:J]; k_inv = (1.0 * DEG) / tauS; hb = 0.5 * 0.3 * DEG
    th_g = q[0].copy(); cmd = q[0].copy()
    if engine == 'syntonic':
        allowed = np.zeros((J, 5), bool); allowed[:, 2] = True                     # fixed 3 ms frame
        eng = CadenceEngine(q[0], DT, w0, allowed, np.zeros(J, bool))
    else:
        kp, ki, kd = (Z[f'{name}_pid'][:J, i] for i in range(3))
        I = np.zeros(J); ef = np.zeros(J); efp = np.zeros(J)
    keep, law = [], []
    for k in range(len(tt)):
        load = tau[k]
        th_out = th_g - load * k_inv - hb * np.tanh(load / (0.02 * tauS))
        a = np.round((th_out + rng.normal(0.0, 1.0, J) * sd[k]) / ENC_LSB) * ENC_LSB
        b = np.round((th_out + rng.normal(0.0, 1.0, J) * sd[k]) / ENC_LSB) * ENC_LSB
        if engine == 'syntonic':
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
            keep.append((th_out[:6].copy(), q[k, :6].copy()))
            law.append(eng.law[:6].copy() if engine == 'syntonic' else np.full(6, np.nan))
    th = np.array([x for x, _ in keep]); qp = np.array([y for _, y in keep])
    e = np.array([tip(x) - tip(y) for x, y in keep])
    return dict(tip_mm=float(1e3 * np.sqrt(np.mean(np.sum(e ** 2, 1)))), th=th, q=qp, law=np.array(law), etip=e)


def jacobian(q6, h=1e-6):
    Jm = np.zeros((3, 6))
    for j in range(6):
        d = np.zeros(6); d[j] = h
        Jm[:, j] = (tip(q6 + d) - tip(q6 - d)) / (2 * h)
    return Jm


def cancellation(r, idx):
    """C over the samples idx: |J e|^2 summed, over the sum of the per-joint contributions |J_j e_j|^2."""
    num = den = 0.0
    for t in idx:
        Jm = jacobian(r['q'][t]); ej = r['th'][t] - r['q'][t]
        contrib = Jm * ej[None, :]
        num += float(np.sum(contrib.sum(1) ** 2)); den += float(np.sum(contrib ** 2))
    return num / den


def phase_stats(r):
    law = r['law'].copy(); n_nan = int(np.isnan(law).sum())       # a joint whose law is not yet valid: filled from
    for j in range(6):                                            # its nearest valid sample (count recorded)
        v = law[:, j]; ok = ~np.isnan(v)
        if ok.any() and not ok.all(): law[:, j] = np.interp(np.arange(len(v)), np.flatnonzero(ok), v[ok])
    x = -np.log(law)                                             # log omega_j, (n, 6)
    xd = x - x.mean(0)
    corr = np.corrcoef(xd.T); iu = np.triu_indices(6, 1); comove = float(np.mean(corr[iu]))
    n = len(x); xs = np.stack([np.roll(xd[:, j], 500 * j) for j in range(6)], 1)
    cs = np.corrcoef(xs.T); comove_shift = float(np.mean(cs[iu]))
    W = 100; nw = n // W
    D = np.array([np.mean(np.std(x[w * W:(w + 1) * W], axis=1)) for w in range(nw)])
    C = np.array([cancellation(r, range(w * W, (w + 1) * W)) for w in range(nw)])
    spd = np.array([np.mean(np.linalg.norm(np.diff(np.array([tip(qq) for qq in r['q'][w * W:(w + 1) * W]]), axis=0), axis=1)) * 100
                    for w in range(nw)])
    rho = float(spearmanr(D, C).correlation)
    rD = spearmanr(D, spd).correlation; rC = spearmanr(C, spd).correlation           # partial Spearman on speed
    rho_partial = float((rho - rD * rC) / np.sqrt(max(1e-12, (1 - rD ** 2) * (1 - rC ** 2))))
    return dict(comove=comove, comove_shift=comove_shift, rho_DC=rho, rho_DC_partial_speed=rho_partial, n_nan_law=n_nan,
                C_run=cancellation(r, range(n)), D=D.tolist(), C=C.tolist())


if __name__ == '__main__':
    sha = hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest()
    if '--smoke' in sys.argv:
        r = run('A', 'steady', 5, 'syntonic', 6, t_end=2.0); print('smoke ok', r['th'].shape, r['law'].shape); sys.exit(0)
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(HERE)), 'examples'))
    import servo_bench as SB
    t_start = time.time(); out = dict(script_sha256=sha, gate=[], runs={})
    g = 0.0
    for name in ('A', 'B'):
        for cond, sbc, seed in (('steady', 'steady noise', 5), ('step', 'noise step', 7)):
            for eng in ('syntonic', 'pid'):
                ref = SB.run(eng, name, sbc, seed)[0]; mine = run(name, cond, seed, eng, 12)['tip_mm']
                g = max(g, abs(ref - mine)); out['gate'].append(dict(cell=f'{name} {cond} {eng}', ref=ref, this=mine))
    out['gate_max_diff_mm'] = g; out['gate_ok'] = g <= 1e-9
    print('gate', g, round(time.time() - t_start), 's', flush=True)
    if not out['gate_ok']:
        json.dump(out, open(os.path.join(HERE, 'results', 'omega_phase_rates.json'), 'w'), indent=1); sys.exit(1)
    for name in ('A', 'B'):
        for cond in ('steady', 'step'):
            for seed in SEEDS:
                rs = run(name, cond, seed, 'syntonic', 6); rp = run(name, cond, seed, 'pid', 6)
                st = phase_stats(rs)
                st['C_run_pid'] = cancellation(rp, range(len(rp['q']))); st['tip_mm'] = rs['tip_mm']; st['tip_mm_pid'] = rp['tip_mm']
                out['runs'][f'{name} {cond} {seed}'] = st
                print(name, cond, seed, {k: round(v, 3) for k, v in st.items() if not isinstance(v, list)},
                      round(time.time() - t_start), 's', flush=True)
    steady = [v for k, v in out['runs'].items() if ' steady ' in k]
    w1 = sum(v['comove'] > v['comove_shift'] for v in steady)
    rhos = [v['rho_DC'] for v in steady]
    verdict = dict(W1=bool(w1 >= 9), W1_runs=int(w1), W2=bool(np.median(rhos) >= 0.2 and sum(r > 0 for r in rhos) >= 8),
                   W2_median_rho=float(np.median(rhos)), W2_positive_runs=int(sum(r > 0 for r in rhos)))
    out['verdict'] = verdict; out['minutes'] = round((time.time() - t_start) / 60, 1)
    json.dump(out, open(os.path.join(HERE, 'results', 'omega_phase_rates.json'), 'w'), indent=1)
    print('verdict', verdict)
    for c in ('steady', 'step'):
        vv = [v for k, v in out['runs'].items() if f' {c} ' in k]
        print(c, 'C syntonic median', round(float(np.median([v['C_run'] for v in vv])), 3), 'C PID median',
              round(float(np.median([v['C_run_pid'] for v in vv])), 3), 'rho median', round(float(np.median([v['rho_DC'] for v in vv])), 3),
              'rho partial median', round(float(np.median([v['rho_DC_partial_speed'] for v in vv])), 3))
