"""Exploratory, after the verdicts (not pre-registered): the Kalman controller WITH the same setpoint feed-forward and motor
lead as SYN1 (KFFF), gains tuned on the training seeds on amendment 4's grid, then flown on the test seeds."""
import os, sys, json, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
R = os.path.abspath(os.path.join(HERE, '..', '..')); sys.argv = ['x', R, 'import']
import flight_bench as F
F.KAS = [8.0, 16.0, 32.0, 64.0]; F.KPS = [80.0 / F.B, 160.0 / F.B, 320.0 / F.B, 640.0 / F.B, 1280.0 / F.B]; F.KDS = [0.0, 5e-4, 1e-3, 2e-3]; F.KIM = [0.0, 2.0]
F.GRID = [(ka, kp, kd, km) for ka in F.KAS for kp in F.KPS for kd in F.KDS for km in F.KIM]
G = F.G
def records_a2():
    if not F._REC:
        for name, seg in (('idle', G.SEGS[0]), ('low', G.SEGS[1]), ('half', G.SEGS[2])):
            tilt, _, g1k = G.load(*seg[:4]); F._REC[name] = (tilt - tilt.mean(), g1k)
    return F._REC
F.records = records_a2
orig = F.arm_spec
def spec(name, fx, gains, tuned):
    if name == 'KFFF': return 1, 1.0, False, True, gains
    return orig(name, fx, gains, tuned)
F.arm_spec = spec
OUT = os.path.join(HERE, 'results_a2'); fx = json.load(open(os.path.join(OUT, 'facts.json'))); tuned = json.load(open(os.path.join(OUT, 'tuned_estimators.json')))
cols = [('KFFF', g) for g in F.GRID]
J = np.mean([np.square(F.fly(s, cols, fx, tuned)['J_att']) for s in F.TRAIN], 0) ** 0.5; i = int(np.argmin(J)); g = F.GRID[i]
print('KFFF tuned gains', g, 'J_train', round(float(J[i]), 3), flush=True)
rows = []
for s in F.TEST:
    r = F.fly(s, [('KFFF', g)], fx, tuned); rows.append(dict(seed=s, J_att=r['J_att'][0], J_est=r['J_est'][0], jit=r['jit'][0], reg=r['att_by_reg']))
P = lambda v: float(np.sqrt(np.mean(np.square(v))))
out = dict(gains=list(g), J_train=float(J[i]), pooled_att=P([r['J_att'] for r in rows]), pooled_est=P([r['J_est'] for r in rows]),
           jit=P([r['jit'] for r in rows]), by_reg={k: P([r['reg'][k][0] for r in rows]) for k in ('low', 'half')}, rows=rows)
json.dump(out, open(os.path.join(OUT, 'explore_kfff.json'), 'w'), indent=1)
syn = [json.loads(l) for l in open(os.path.join(OUT, 'test_runs.jsonl'))]
for name in ('SYN1', 'SYNB'):
    i2 = syn[0]['arms'].index(name); a = {r['seed']: r['J_att'][i2] for r in syn}
    print(name, 'vs KFFF wins', sum(a[r['seed']] < r['J_att'] for r in rows), '/ 10')
print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in out.items() if k != 'rows'})
