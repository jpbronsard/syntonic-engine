# Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
# Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
"""The PID given the same lead and the same frames as the syntonic engine, re-tuned for them. Written 2 Nov 2026, before
any run of this file except the smoke test; request from JP ("FAIT LE !") after the assistant wrote that the 30 % gain
of the best configuration (eg732_hold.py, F15D) was the engine against itself, not against the PID.

The comparison so far: the bench's PID (examples/servo_bench.py) runs at the published setting (1 kHz loop, 3 ms frames,
no lead) with gains tuned there; the engine's best configuration runs at a 2 kHz loop, fixed 1.5 ms frames and the
lead (servo_tau + D) x planned velocity. Here both get the same loop, frames and lead, and the PID is re-tuned for them.

SETTING (baton arm RJ1-RJ6, the servo model of examples/servo_bench.py; loop 2 kHz, fixed 1.5 ms frames)
  Lead for every arm that has one: (servo_tau + D) x planned velocity, D = loop_delay(dt, frame) = frame/2 + dt/2 +
  0.1 ms (26.1 ms in all), the same coefficient for the PID and the engine. The PID's command at a frame is plan + its
  clipped correction + the lead; the PID is unchanged otherwise (derivative on the error low-passed at 5 ms, integral
  frozen while the raw command exceeds the clip, clip at servo_tau x w0, the authority limit of 0.2.8).
  Arms:
    PIDD    PID with the lead, gains re-tuned here.                                     (judged)
    SYND    engine 0.2.9, lead='servo+delay', identification 'raw', published travel     (judged; = F15D)
    SYNDE   SYND with travel='e' (= E64 of transport_confirm.py)                         (described)
    PID15   PID without the lead, same loop and frames, gains re-tuned here              (described)
    SYN15   engine without the lead, same loop and frames (= F15)                        (described)
    PIDDp   PID with the lead and the published gains (Z['<excerpt>_pid'])               (described)
  Tuning of PIDD and PID15, the method of the published bench: joint by joint, the (Kp, Ki, Kd) of the published grid
  of 280 (Kp in {0, 0.5, 1, 2, 3, 4, 6, 8}, Ki in {0, 80, 160, 320, 640, 1280, 2560} /s, Kd in {0, 0.5, 1, 2, 5} ms)
  with the lowest joint error (rms of plan - joint angle over the 30 s), steady noise, per excerpt, mean square over
  tuning seeds 61, 62, 63. Declared now: the PID is tuned on the same two excerpts it is tested on (other noise draws),
  which favours it; the published PID was tuned on the rest of the movement, which this package does not contain.
  Edge rule, declared now: if any joint's choice in an arm sits on the top edge of an axis (Kp 8, Ki 2560, Kd 5 ms),
  that axis is extended once (Kp + {12, 16}; Ki + {5120}; Kd + {10 ms}) and the arm re-tuned on the extended grid,
  before any test run; recorded in tuned_gains.json.
  Test: excerpts A and B, steady noise and the x10 noise step (10 to 20 s), seeds 41 to 45 (never used): 20 paired
  runs, every arm sharing every encoder read.
  Gates: G1, at the published setting (12 joints, 1 kHz, 3 ms, no lead, published gains) this harness reproduces
  examples/servo_bench.py, PID and engine, on A steady seed 5 and B noise step seed 7, difference 0; G2, at the fast
  setting, SYND and SYN15 on A steady seed 1 reproduce F15D of eg732_hold.py and F15 of eg732_lead.py, difference 0.

PREDICTIONS (judged; tip = rms baton-tip error over the 30 s; pooled = rms over runs of the per-run tip)
R1  Noise step, the law's ground: pooled SYND <= 0.85 x pooled PIDD over the 10 step runs, SYND below in >= 9 of 10.
R2  Steady noise, parity: pooled SYND <= 1.10 x pooled PIDD over the 10 steady runs.
R3  Overall: pooled SYND <= pooled PIDD over the 20 runs, SYND below in >= 12 of 20.
R4  The lead helps the engine more than the PID: 1 - pooled SYND / pooled SYN15 > 1 - pooled PIDD / pooled PID15.
Described, not judged: per cell; SYNDE, PIDDp; joint errors (the PID's own target); jitter; the tuned gains.

Run: python pid_lead_fair.py <repo_root> smoke | gate | tune | run <excerpt> <steady|step> | summarize
"""
import os, sys, json, time, hashlib
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.signal import butter, sosfiltfilt

ROOT = os.path.abspath(sys.argv[1]); MODE = sys.argv[2]
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(ROOT, 'python')); sys.path.insert(0, os.path.join(ROOT, 'examples'))
from syntonic_engine import SyntonicEngine
from syntonic_engine.engine import loop_delay
from arm_fk import tip
OUT = os.path.join(HERE, 'results_pid_lead'); os.makedirs(OUT, exist_ok=True)
SHA = hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest()

DEG = np.pi / 180
Z = np.load(os.path.join(ROOT, 'examples', 'data', 'maestro_excerpts.npz'))
TAU_S, W0 = Z['stall_torque'], Z['no_load_speed']
ENC_SD, ENC_LSB = float(Z['enc_sd_rad']), float(Z['enc_lsb_rad'])
TS, HB = 0.025, 0.5 * 0.3 * DEG
K_INV = (1.0 * DEG) / TAU_S
KP = [0, 0.5, 1, 2, 3, 4, 6, 8]; KI = [0, 80, 160, 320, 640, 1280, 2560]; KD = [0, 0.0005, 0.001, 0.002, 0.005]
EXT = dict(kp=[12, 16], ki=[5120], kd=[0.01])
TUNE_SEEDS, TEST_SEEDS = [61, 62, 63], [41, 42, 43, 44, 45]
DT_F, FR_F = 5e-4, 1.5e-3
BATON = list(range(6))


def grid(kp=KP, ki=KI, kd=KD):
    return [(a, b, c) for a in kp for b in ki for c in kd]


def simulate(name, cond, seed, dt, frame, joints, pid, syn, tip_score=True):
    """pid: list of (label, joint index per column, kp, ki, kd, lead coefficient per column); syn: list of (label, options).
    Every column reads the same encoder noise draws as the published bench (one normal per joint and tick)."""
    J = len(joints); t0 = float(Z[f'{name}_t0']); tt = np.arange(t0 - 3.0, t0 + 30.0, dt)
    S = CubicSpline(Z[f'{name}_t30'], Z[f'{name}_q30'][:, joints], axis=0); qJ = S(tt); qdJ = S(tt, 1)
    tauJ = CubicSpline(Z[f'{name}_t30'], Z[f'{name}_tau30'][:, joints], axis=0)(tt)
    sd = np.full(len(tt), ENC_SD)
    if cond == 'step': sd[(tt >= t0 + 10) & (tt < t0 + 20)] *= np.sqrt(10.0)
    rng = np.random.default_rng(seed)
    tS, w0J, kinv = TAU_S[joints], W0[joints], K_INV[joints]
    fe = max(1, int(round(frame / dt)))
    # PID block
    jc = np.concatenate([p[1] for p in pid]) if pid else np.zeros(0, int)
    kp = np.concatenate([p[2] for p in pid]) if pid else np.zeros(0); ki = np.concatenate([p[3] for p in pid]) if pid else np.zeros(0)
    kd = np.concatenate([p[4] for p in pid]) if pid else np.zeros(0); ld = np.concatenate([p[5] for p in pid]) if pid else np.zeros(0)
    C = len(jc); LIM = TS * w0J[jc]
    th_p = qJ[0][jc].copy(); cmd_p = qJ[0][jc].copy(); I = np.zeros(C); ef = np.zeros(C); efp = np.zeros(C)
    # engines
    engs = [SyntonicEngine(qJ[0], dt, frame, w0J, identification='raw', **o) for _, o in syn]
    th_s = [qJ[0].copy() for _ in syn]
    every = int(round(0.01 / dt)); keep_p, keep_s, keep_q = [], [[] for _ in syn], []
    err_p = np.zeros(C); err_s = [np.zeros(J) for _ in syn]; m = 0
    for k in range(len(tt)):
        load = tauJ[k]; q = qJ[k]; qd = qdJ[k]
        na = rng.normal(0.0, 1.0, J); nb = rng.normal(0.0, 1.0, J)
        w_eff = w0J * np.maximum(0.05, 1.0 - np.abs(load) / tS)
        if C:
            lp = load[jc]
            out_p = th_p - lp * kinv[jc] - HB * np.tanh(lp / (0.02 * tS[jc]))
            a = np.round((out_p + na[jc] * sd[k]) / ENC_LSB) * ENC_LSB
            b = np.round((out_p + nb[jc] * sd[k]) / ENC_LSB) * ENC_LSB
            e = q[jc] - 0.5 * (a + b)
            ef += (dt / 0.005) * (e - ef)
            raw = kp * e + ki * I + kd * (ef - efp) / dt
            I += np.where(np.abs(raw) >= LIM, 0.0, e * dt); efp = ef.copy()
            corr = np.clip(raw, -LIM, LIM)
            if k % fe == 0: cmd_p = q[jc] + corr + ld * qd[jc]
            th_p = th_p + dt * np.clip((cmd_p - th_p) / TS, -w_eff[jc], w_eff[jc])
        outs = []
        for i, (lab, o) in enumerate(syn):
            out = th_s[i] - load * kinv - HB * np.tanh(load / (0.02 * tS)); outs.append(out)
            a = np.round((out + na * sd[k]) / ENC_LSB) * ENC_LSB
            b = np.round((out + nb * sd[k]) / ENC_LSB) * ENC_LSB
            cmd = engs[i].step(q, a, b, qd_plan=qd if o.get('lead') else None)
            th_s[i] = th_s[i] + dt * np.clip((cmd - th_s[i]) / TS, -w_eff, w_eff)
        if tt[k] >= t0:
            if C: err_p += (q[jc] - out_p) ** 2
            for i in range(len(syn)): err_s[i] += (q - outs[i]) ** 2
            m += 1
        if tip_score and k % every == 0 and tt[k] >= t0 - 1e-9:
            if C: keep_p.append(out_p.copy())
            for i in range(len(syn)): keep_s[i].append(outs[i][:6].copy())
            keep_q.append(q[:6].copy())
    res = {'pid_joint_deg': np.sqrt(err_p / m) / DEG}
    if tip_score:
        sos = butter(4, 5.0, 'highpass', fs=100.0, output='sos')
        def score(rows):
            e = np.array([tip(a) - tip(b) for a, b in zip(rows, keep_q)])
            hp = sosfiltfilt(sos, e, axis=0)
            return float(1e3 * np.sqrt(np.mean(np.sum(e ** 2, 1)))), float(1e3 * np.sqrt(np.mean(np.sum(hp ** 2, 1))))
        off = 0
        for p in pid:
            n = len(p[1]); cols = slice(off, off + n); off += n
            tp, jt = score([r[cols][:6] for r in keep_p])
            res[p[0]] = dict(tip_mm=tp, jitter_mm=jt, joint_deg=float(np.mean(res['pid_joint_deg'][cols][:6])))
        for i, (lab, _) in enumerate(syn):
            tp, jt = score(keep_s[i]); res[lab] = dict(tip_mm=tp, jitter_mm=jt, joint_deg=float(np.mean(np.sqrt(err_s[i] / m)[:6] / DEG)))
    return res


def pid_arm(label, gains, lead):
    """gains: (J, 3) per joint, baton arm."""
    g = np.asarray(gains, float); J = len(g)
    return (label, np.arange(J), g[:, 0], g[:, 1], g[:, 2], np.full(J, lead))


LEAD = TS + loop_delay(DT_F, FR_F)
SYN_ARMS = [('SYND', dict(lead='servo+delay')), ('SYNDE', dict(lead='servo+delay', travel='e')), ('SYN15', {})]


def tune_arm(lead, G):
    """per excerpt, per joint: the grid point with the lowest joint error, mean square over the tuning seeds."""
    G = np.array(G); n = len(G); out = {}
    for name in ('A', 'B'):
        cols = (f'T', np.repeat(np.arange(6), n), np.tile(G[:, 0], 6), np.tile(G[:, 1], 6), np.tile(G[:, 2], 6), np.full(6 * n, lead))
        J2 = np.mean([simulate(name, 'steady', s, DT_F, FR_F, BATON, [cols], [], tip_score=False)['pid_joint_deg'] ** 2 for s in TUNE_SEEDS], 0)
        J2 = J2.reshape(6, n); best = np.argmin(J2, 1)
        out[name] = dict(gains=G[best].tolist(), joint_deg=np.sqrt(J2[np.arange(6), best]).tolist())
    return out


def on_edge(t, ax):
    i = {'kp': 0, 'ki': 1, 'kd': 2}[ax]; top = {'kp': max(KP), 'ki': max(KI), 'kd': max(KD)}[ax]
    return any(g[i] >= top for name in ('A', 'B') for g in t[name]['gains'])


if __name__ == '__main__':
    if MODE == 'smoke':
        t = time.time(); g = np.array(Z['A_pid'][:6])
        r = simulate('A', 'step', 1, DT_F, FR_F, BATON, [pid_arm('PIDD', g, LEAD)], SYN_ARMS[:1])
        print('smoke ok', {k: v for k, v in r.items() if k != 'pid_joint_deg'}, round(time.time() - t, 1), 's'); sys.exit(0)
    if MODE == 'gate':
        sys.argv = [sys.argv[0]]; import servo_bench as SB
        G = {}
        for name, cond, seed in (('A', 'steady', 5), ('B', 'step', 7)):
            r = simulate(name, cond, seed, 1e-3, 3e-3, list(range(12)), [pid_arm('PID', Z[f'{name}_pid'], 0.0)], [('SYN', {})])
            c = 'steady noise' if cond == 'steady' else 'noise step'
            ref_p = SB.run('pid', name, c, seed)[0]; ref_s = SB.run('syntonic', name, c, seed)[0]
            G[f'G1 {name} {cond} {seed}'] = dict(pid=r['PID']['tip_mm'], pid_ref=ref_p, syn=r['SYN']['tip_mm'], syn_ref=ref_s,
                                                ok=bool(r['PID']['tip_mm'] == ref_p and r['SYN']['tip_mm'] == ref_s))
        r = simulate('A', 'steady', 1, DT_F, FR_F, BATON, [], SYN_ARMS[:1] + SYN_ARMS[2:])
        ref = lambda fn, arm: [json.loads(l) for l in open(os.path.join(HERE, fn))][0][arm]['tip_mm']
        f15d, f15 = ref('results_hold/runs_A_steady.jsonl', 'F15D'), ref('results/runs_A_steady.jsonl', 'F15')
        G['G2 A steady 1'] = dict(SYND=r['SYND']['tip_mm'], F15D=f15d, SYN15=r['SYN15']['tip_mm'], F15=f15,
                                  ok=bool(r['SYND']['tip_mm'] == f15d and r['SYN15']['tip_mm'] == f15))
        G['script_sha256'] = SHA; json.dump(G, open(os.path.join(OUT, 'gates.json'), 'w'), indent=1); print(json.dumps(G, indent=1)); sys.exit(0)
    if MODE == 'tune':
        res = {}
        for arm, lead in (('PIDD', LEAD), ('PID15', 0.0)):
            t = time.time(); tr = tune_arm(lead, grid()); ext = {}
            axes = [ax for ax in ('kp', 'ki', 'kd') if on_edge(tr, ax)]
            if axes:                                            # the declared edge rule, applied once
                ext = {ax: EXT[ax] for ax in axes}
                tr = tune_arm(lead, grid(KP + ext.get('kp', []), KI + ext.get('ki', []), KD + ext.get('kd', [])))
            res[arm] = dict(tuned=tr, extended=ext); print(arm, json.dumps(res[arm]), round(time.time() - t), 's', flush=True)
        res['script_sha256'] = SHA; json.dump(res, open(os.path.join(OUT, 'tuned_gains.json'), 'w'), indent=1); sys.exit(0)
    if MODE == 'run':
        name, cond = sys.argv[3], sys.argv[4]; T = json.load(open(os.path.join(OUT, 'tuned_gains.json')))
        pids = [pid_arm('PIDD', T['PIDD']['tuned'][name]['gains'], LEAD), pid_arm('PID15', T['PID15']['tuned'][name]['gains'], 0.0),
                pid_arm('PIDDp', Z[f'{name}_pid'][:6], LEAD)]
        for seed in TEST_SEEDS:
            t = time.time(); r = simulate(name, cond, seed, DT_F, FR_F, BATON, pids, SYN_ARMS); r.pop('pid_joint_deg')
            with open(os.path.join(OUT, f'runs_{name}_{cond}.jsonl'), 'a') as f: f.write(json.dumps(dict(key=f'{name} {cond} {seed}', sha=SHA, **r)) + '\n')
            print(name, cond, seed, {a: round(v['tip_mm'], 3) for a, v in r.items()}, round(time.time() - t), 's', flush=True)
        sys.exit(0)
    if MODE == 'summarize':
        runs = {}
        for fn in sorted(os.listdir(OUT)):
            if fn.startswith('runs_'):
                for l in open(os.path.join(OUT, fn)): x = json.loads(l); runs[x.pop('key')] = x
        assert len(runs) == 20, len(runs)
        G = json.load(open(os.path.join(OUT, 'gates.json'))); gates_ok = all(v['ok'] for k, v in G.items() if k.startswith('G'))
        arms = ['PIDD', 'SYND', 'SYNDE', 'PID15', 'SYN15', 'PIDDp']
        P = lambda ks, a: float(np.sqrt(np.mean([runs[k][a]['tip_mm'] ** 2 for k in ks])))
        allk = sorted(runs); st = [k for k in allk if ' steady ' in k]; sp = [k for k in allk if ' step ' in k]
        win = lambda ks, x, y: int(sum(runs[k][x]['tip_mm'] < runs[k][y]['tip_mm'] for k in ks))
        pooled = {a: dict(all=P(allk, a), steady=P(st, a), step=P(sp, a)) for a in arms}
        V = dict(R1=pooled['SYND']['step'] <= 0.85 * pooled['PIDD']['step'] and win(sp, 'SYND', 'PIDD') >= 9,
                 R2=pooled['SYND']['steady'] <= 1.10 * pooled['PIDD']['steady'],
                 R3=pooled['SYND']['all'] <= pooled['PIDD']['all'] and win(allk, 'SYND', 'PIDD') >= 12,
                 R4=(1 - pooled['SYND']['all'] / pooled['SYN15']['all']) > (1 - pooled['PIDD']['all'] / pooled['PID15']['all']))
        V = {k: bool(v and gates_ok) for k, v in V.items()}
        cells = {c: {a: P([k for k in allk if k.startswith(c)], a) for a in arms} for c in ('A steady', 'A step', 'B steady', 'B step')}
        S = dict(script_sha256=SHA, gates_ok=gates_ok, pooled_tip_mm=pooled, cells=cells, verdict=V,
                 ratios=dict(SYND_PIDD={g: pooled['SYND'][g] / pooled['PIDD'][g] for g in ('all', 'steady', 'step')},
                             SYNDE_PIDD={g: pooled['SYNDE'][g] / pooled['PIDD'][g] for g in ('all', 'steady', 'step')}),
                 wins=dict(SYND_vs_PIDD_step=win(sp, 'SYND', 'PIDD'), SYND_vs_PIDD_steady=win(st, 'SYND', 'PIDD'), SYND_vs_PIDD_all=win(allk, 'SYND', 'PIDD'),
                           SYNDE_vs_PIDD_all=win(allk, 'SYNDE', 'PIDD')),
                 lead_gain=dict(engine=1 - pooled['SYND']['all'] / pooled['SYN15']['all'], pid=1 - pooled['PIDD']['all'] / pooled['PID15']['all']),
                 joint_deg={a: float(np.sqrt(np.mean([runs[k][a]['joint_deg'] ** 2 for k in allk]))) for a in arms},
                 jitter_mm={a: float(np.sqrt(np.mean([runs[k][a]['jitter_mm'] ** 2 for k in allk]))) for a in arms})
        json.dump(S, open(os.path.join(OUT, 'summary.json'), 'w'), indent=1); print(json.dumps(S, indent=1))
