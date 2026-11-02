"""Development (seeds 91, 92 only): the frozen-candidate law on the realistic bench, commissioned at the PWM latency
(servo internal delay 2.5 ms, the agent's mid-travel figure): PID + lead + frame-mean P re-tuned there (the agent's
tune_pid, 280 grid + edge rule, seeds 61-63); the law given dead = 2.5 ms as its delay fact."""
import sys, os, json, time
sys.path.insert(0, '/home/claude/work/goal/servo'); sys.path.insert(0, '/home/claude/work/goal/confirm')
import numpy as np
import realistic_bench as RB, rivals as RV
from syntonic_loop_law import SyntonicLoopLaw
COMM = dict(dead_time=0.0025)
if sys.argv[1] == 'tune':
    out = {}
    for arm, pm in (('leadmean', True), ('lead', False)):
        r = RV.tune_pid(lead=True, p_mean=pm, circumstance=COMM, log=None); r.pop('_grids', None); out[arm] = r
        print(arm, r['gains'], r['extended'], round(r['wall_s']), 's', flush=True)
    json.dump(out, open('dev_tune_pwm.json', 'w'), indent=1); sys.exit(0)
T = json.load(open('dev_tune_pwm.json'))
name, seed, circ = sys.argv[1], int(sys.argv[2]), json.loads(sys.argv[3])
C = dict(COMM, **circ)
ctr = [('PIDLM', {n: RV.pid_factory(np.array(T['leadmean']['gains'][n]), lead=True, p_mean=True) for n in 'AB'}),
       ('PIDL', {n: RV.pid_factory(np.array(T['lead']['gains'][n]), lead=True, p_mean=False) for n in 'AB'}),
       ('LAW', lambda q0, w0: SyntonicLoopLaw(q0, w0, dead=0.0025))]
t = time.time(); r = RB.simulate(name, C, seed, ctr)
print(name, seed, circ, {k: round(v['tip_mm'], 3) for k, v in r.items()}, round(time.time() - t), 's', flush=True)
with open('dev_pwm.jsonl', 'a') as f: f.write(json.dumps(dict(key=f'{name} {seed} {json.dumps(circ)}', **{k: {kk: (vv if not hasattr(vv, 'tolist') else vv.tolist()) for kk, vv in v.items()} for k, v in r.items()})) + '\n')
