# Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
# Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
"""Gate for the 0.2.9 engine options: the published SyntonicEngine with lead='servo+delay' (and travel='e') must
reproduce, to the last bit, the arms these benches ran with code of their own (written 2 Nov 2026 for 0.2.9; a check,
not a test).

  lead only      against BASE of transport_confirm.py (the F15D setting of eg732_hold.py: 2 kHz loop, fixed 1.5 ms
                 frames, lead (servo_tau + D) x planned velocity, the published travel)
  lead + 'e'     against E64 of transport_confirm.py (1/tau relaxing with T = 64 ms, the longest lag)

Same harness as window_transport.simulate (baton arm, identification 'raw'); seeds and cells: A steady 21, B step 25.
Run: python engine_options_gate.py <repo_root>
"""
import os, sys, json
import numpy as np
from scipy.interpolate import CubicSpline

ROOT = os.path.abspath(sys.argv[1]); HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
_argv = sys.argv; sys.argv = [_argv[0], ROOT, 'import']
import eg732_hold as H                                                # bench constants and arm geometry
from syntonic_engine import SyntonicEngine
sys.argv = _argv


def simulate(name, cond, seed, arms, t_end=30.0):
    """window_transport.simulate with the published engine in place of TransportEngine."""
    dt = 5e-4; joints = list(range(6)); J = 6; A = len(arms)
    t0 = float(H.Z[f'{name}_t0']); tt = np.arange(t0 - 3.0, t0 + t_end, dt)
    S = CubicSpline(H.Z[f'{name}_t30'], H.Z[f'{name}_q30'][:, joints], axis=0); qJ = S(tt); qdJ = S(tt, 1)
    tauJ = CubicSpline(H.Z[f'{name}_t30'], H.Z[f'{name}_tau30'][:, joints], axis=0)(tt)
    sd = np.full(len(tt), H.ENC_SD)
    if cond == 'step': sd[(tt >= t0 + 10) & (tt < t0 + 20)] *= np.sqrt(10.0)
    rng = np.random.default_rng(seed)
    tauS = np.tile(H.TAU_S[joints], A); w0 = np.tile(H.W0[joints], A)
    k_inv = (1.0 * H.DEG) / tauS; hb = 0.5 * 0.3 * H.DEG
    engs = [SyntonicEngine(qJ[0], dt, 1.5e-3, H.W0[joints], identification='raw', lead='servo+delay', travel=tr)
            for tr in arms]
    th_g = np.tile(qJ[0], A); every = int(round(0.01 / dt)); keep = []
    for k in range(len(tt)):
        load = np.tile(tauJ[k], A); q = qJ[k]; qd = qdJ[k]
        th_out = th_g - load * k_inv - hb * np.tanh(load / (0.02 * tauS))
        na = np.tile(rng.normal(0.0, 1.0, J), A); nb = np.tile(rng.normal(0.0, 1.0, J), A)
        ra = np.round((th_out + na * sd[k]) / H.ENC_LSB) * H.ENC_LSB
        rb = np.round((th_out + nb * sd[k]) / H.ENC_LSB) * H.ENC_LSB
        cmd = np.concatenate([e.step(q, ra[i * J:(i + 1) * J], rb[i * J:(i + 1) * J], qd_plan=qd) for i, e in enumerate(engs)])
        w_eff = w0 * np.maximum(0.05, 1.0 - np.abs(load) / tauS)
        th_g = th_g + dt * np.clip((cmd - th_g) / H.TS, -w_eff, w_eff)
        if tt[k] >= t0 - 1e-9 and k % every == 0: keep.append((th_out.copy(), np.tile(q, A)))
    out = {}
    for i, lab in enumerate(arms):
        s = slice(i * J, (i + 1) * J)
        e2 = np.array([np.sum((H.tip(a[s]) - H.tip(b[s])) ** 2) for a, b in keep])
        out[lab] = float(1e3 * np.sqrt(np.mean(e2)))
    return out


if __name__ == '__main__':
    ref = {}
    for fn in ('runs_A_steady.jsonl', 'runs_B_step.jsonl'):
        for line in open(os.path.join(HERE, 'results_confirm', fn)):
            x = json.loads(line); ref[x['key']] = x
    res = {}
    for key in ('A steady 21', 'B step 25'):
        name, cond, seed = key.split()
        r = simulate(name, cond, int(seed), ['rate', 'e'])
        res[key] = dict(lead=r['rate'], lead_ref=ref[key]['BASE']['tip_mm'], lead_e=r['e'], lead_e_ref=ref[key]['E64']['tip_mm'])
        print(key, res[key], flush=True)
    worst = max(max(abs(v['lead'] - v['lead_ref']), abs(v['lead_e'] - v['lead_e_ref'])) for v in res.values())
    print('max |engine option - bench arm| =', worst, 'mm;', 'ok' if worst == 0.0 else 'NOT exact')
    json.dump(dict(runs=res, max_diff_mm=worst), open(os.path.join(HERE, 'results_engine_options_gate.json'), 'w'), indent=1)
    sys.exit(0 if worst == 0.0 else 1)
