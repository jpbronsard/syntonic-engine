# Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
# Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
"""Pre-registered confirmation: the syntonic loop law against tuned PIDs on a realistic MAESTRO 1 servo bench, servo
frames 333 Hz, in every circumstance of a declared battery. Written 2 Nov 2026 before any run of this file except the
smoke test; goal set by JP: "trouve comment la loi peut a la fois etre a parite ou depasser un PID regle et egalement
egaler ou depasser un Kalman regle en toute circonstances. Servo 333 hz max. Test realiste".

FROZEN MODULES (sha256, checked at start):
  bench/realistic_bench.py  1fe422ea21f32633d88c88e274565f237fb6e40282a41d9f874043aa07f461e8  (written by a sub-agent:
      baton arm RJ1-RJ6 of MAESTRO 1, 1 kHz loop, 3 ms frames (333 Hz) held, first-order servo 25 ms with load-derated
      speed limit, 1 deg compliance at stall, 0.3 deg gear play, AS5048 read twice per tick, command quantization
      0.025 deg; reproduces examples/servo_bench.py bit for bit at the published plant: gate G1 below)
  bench/rivals.py           e4bf3df3479c7d86390d111e24b3971d3c3145204a5409f922e9d66132c68691  (the published PID with the
      lead and the frame-mean proportional term, and its tuning on the published 280 grid + edge rule)
  syntonic_loop_law.py      a7a234ccefb2c9bb963b87b9a3268493644127bfd43d44c387e1411b3910b2b3  (the law; its docstring
      states the rule; developed on seeds 91-93, never on 71-75)

COMMISSIONING (the hardware as installed; the PIDs are tuned here, the law measures its facts here):
  PWM servo latency: the servo knows a pulse width only at its falling edge, 1.5 ms after the frame start at mid-travel,
  plus 1 ms of internal processing: dead_time = 2.5 ms (sub-agent's figure). Other plant options: bench defaults.
  Rivals, tuned at commissioning with rivals.tune_pid (joint by joint, per excerpt, joint error, steady noise, tuning
  seeds 61-63, edge rule):
    PIDLM  PID + lead (servo_tau + D) x planned velocity + frame-mean P   (the primary rival, same input processing)
    PIDL   PID + lead                                                     (the conventional tuned PID)
  Law: SyntonicLoopLaw(dead=2.5 ms): its facts are dt, frame, servo_tau = 25 ms, w0 (data sheet) and the measured
  2.5 ms delay; the same lead coefficient as the PIDs (27.1 ms); nothing tuned.

CIRCUMSTANCES IN OPERATION (no re-tuning, no new facts; each = commissioning + the change):
  C1 nominal          C2 noise_step (read variance x10, 10-20 s)     C3 bursts (x30, 0.5 s every 2.5 s, 5-25 s)
  C4 latency_long     (dead time 3.5 ms: the far end of PWM travel)  C5 latency_short (1.5 ms: the near end)
  C6 servo_slow       (true time constant 35 ms)                     C7 servo_fast (18 ms)
  C8 load_x1.2        (loads x1.2: a heavier baton; x1.5 exceeds the stall torque and was dropped as infeasible)
  C9 quant_coarse     (command quantization 0.09 deg)
  Excerpts A and B, seeds 71-75 (reserved, never used): 10 paired runs per circumstance. Score: baton-tip rms (mm),
  pooled = rms over the 10 runs.

PREDICTIONS (judged; ratio = pooled LAW / pooled rival)
S1  Parity in every circumstance, by the corpus's criterion (regret under 10 %): ratio to PIDLM <= 1.10 in C1-C9.
S2  Better where the law re-measures what changed (noise C2, C3; servo time constant C6, C7; load C8): ratio to
    PIDLM <= 1.00 in each of them.
S3  Better than the conventional tuned PID in every circumstance: ratio to PIDL <= 1.00 in C1-C9.
Not claimed: better in C1 (the PIDs' tuning point), C4, C5 (the delay fact is wrong for both), C9.
Described, not judged: per excerpt, wins per circumstance, the law's gains.
Gates: G1, at the published plant (12 joints), the bench reproduces examples/servo_bench.py (PID A steady seed 5
0.6761419910844138 mm, engine 0.6839103668504339 mm; PID B step seed 7 1.3144217897870654, engine 1.0973289271118845),
difference 0; G0, module hashes as above.

Run: python confirm_servo.py smoke | gate | tune | run <C1..C9> | summarize
"""
import os, sys, json, time, hashlib
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, 'bench')); sys.path.insert(0, HERE)
import realistic_bench as RB, rivals as RV
from syntonic_loop_law import SyntonicLoopLaw

SHA = hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest()
HASHES = {'bench/realistic_bench.py': '1fe422ea21f32633d88c88e274565f237fb6e40282a41d9f874043aa07f461e8',
          'bench/rivals.py': 'e4bf3df3479c7d86390d111e24b3971d3c3145204a5409f922e9d66132c68691',
          'syntonic_loop_law.py': 'a7a234ccefb2c9bb963b87b9a3268493644127bfd43d44c387e1411b3910b2b3'}
for f, h in HASHES.items():
    got = hashlib.sha256(open(os.path.join(HERE, f), 'rb').read()).hexdigest()
    assert got == h, f'G0 failed: {f} {got}'
OUT = os.path.join(HERE, 'results_servo'); os.makedirs(OUT, exist_ok=True)
COMM = dict(dead_time=0.0025)
CIRC = {'C1': {}, 'C2': dict(noise='step'), 'C3': dict(noise='bursts'), 'C4': dict(dead_time=0.0035),
        'C5': dict(dead_time=0.0015), 'C6': dict(servo_tau=0.035), 'C7': dict(servo_tau=0.018), 'C8': dict(load_scale=1.2),
        'C9': dict(cmd_quant_deg=0.09)}
SEEDS = [71, 72, 73, 74, 75]
REF = {('A', 'steady', 5): (0.6761419910844138, 0.6839103668504339), ('B', 'step', 7): (1.3144217897870654, 1.0973289271118845)}


def controllers(T):
    return [('PIDLM', {n: RV.pid_factory(np.array(T['leadmean']['gains'][n]), lead=True, p_mean=True) for n in 'AB'}),
            ('PIDL', {n: RV.pid_factory(np.array(T['lead']['gains'][n]), lead=True, p_mean=False) for n in 'AB'}),
            ('LAW', lambda q0, w0: SyntonicLoopLaw(q0, w0, dt=RB.DT, frame=RB.FRAME, servo_tau=RB.TS_NOM, dead=0.0025))]


if __name__ == '__main__':
    MODE = sys.argv[1]
    if MODE == 'smoke':
        T = json.load(open(os.path.join(HERE, 'dev_tune_pwm.json')))
        t = time.time(); r = RB.simulate('A', dict(COMM), 91, controllers(T))
        print('smoke ok', {k: round(v['tip_mm'], 3) for k, v in r.items()}, round(time.time() - t), 's'); sys.exit(0)
    if MODE == 'gate':
        G = {}
        for (name, cond, seed), (rp, rs) in REF.items():
            r = RB.simulate(name, f'published_{cond}', seed, [('PID', RV.pid_factory(RV.published_gains(name), lead=False)),
                                                              ('SYN', RV.syntonic_factory(lead=False))], n_joints=12)
            G[f'G1 {name} {cond} {seed}'] = dict(pid=r['PID']['tip_mm'], syn=r['SYN']['tip_mm'], ok=bool(r['PID']['tip_mm'] == rp and r['SYN']['tip_mm'] == rs))
        G['script_sha256'] = SHA; json.dump(G, open(os.path.join(OUT, 'gates.json'), 'w'), indent=1); print(json.dumps(G, indent=1)); sys.exit(0)
    if MODE == 'tune':
        out = {}
        for arm, pm in (('leadmean', True), ('lead', False)):
            r = RV.tune_pid(lead=True, p_mean=pm, circumstance=dict(COMM), log=None); r.pop('_grids', None); out[arm] = r
            print(arm, r['gains'], r['extended'], flush=True)
        out['script_sha256'] = SHA; json.dump(out, open(os.path.join(OUT, 'tuned.json'), 'w'), indent=1); sys.exit(0)
    if MODE == 'run':
        c = sys.argv[2]; T = json.load(open(os.path.join(OUT, 'tuned.json'))); P = dict(COMM, **CIRC[c])
        for name in ('A', 'B'):
            for s in SEEDS:
                t = time.time(); r = RB.simulate(name, P, s, controllers(T), allow_reserved=True)
                row = dict(key=f'{c} {name} {s}', sha=SHA, **{k: dict(tip_mm=v['tip_mm'], jitter_mm=v['jitter_mm'], joint_deg=v['joint_deg']) for k, v in r.items()})
                with open(os.path.join(OUT, f'runs_{c}.jsonl'), 'a') as f: f.write(json.dumps(row) + '\n')
                print(c, name, s, {k: round(v['tip_mm'], 3) for k, v in r.items()}, round(time.time() - t), 's', flush=True)
        sys.exit(0)
    if MODE == 'summarize':
        G = json.load(open(os.path.join(OUT, 'gates.json'))); gok = all(v['ok'] for k, v in G.items() if k.startswith('G1'))
        rows = {}
        for c in CIRC:
            fn = os.path.join(OUT, f'runs_{c}.jsonl')
            rows[c] = [json.loads(l) for l in open(fn)] if os.path.exists(fn) else []
        P = lambda rs, a: float(np.sqrt(np.mean([r[a]['tip_mm'] ** 2 for r in rs])))
        tab = {c: {a: P(rs, a) for a in ('PIDLM', 'PIDL', 'LAW')} for c, rs in rows.items() if len(rs) == 10}
        ratio = {c: dict(LAW_PIDLM=v['LAW'] / v['PIDLM'], LAW_PIDL=v['LAW'] / v['PIDL'],
                         wins_vs_PIDLM=int(sum(r['LAW']['tip_mm'] < r['PIDLM']['tip_mm'] for r in rows[c]))) for c, v in tab.items()}
        complete = len(tab) == len(CIRC)
        V = dict(S1=complete and all(r['LAW_PIDLM'] <= 1.10 for r in ratio.values()),
                 S2=complete and all(ratio[c]['LAW_PIDLM'] <= 1.00 for c in ('C2', 'C3', 'C6', 'C7', 'C8')),
                 S3=complete and all(r['LAW_PIDL'] <= 1.00 for r in ratio.values()))
        V = {k: bool(v and gok) for k, v in V.items()}
        S = dict(script_sha256=SHA, gates_ok=gok, complete=complete, pooled_tip_mm=tab, ratios=ratio, verdict=V)
        json.dump(S, open(os.path.join(OUT, 'summary.json'), 'w'), indent=1); print(json.dumps(S, indent=1))
