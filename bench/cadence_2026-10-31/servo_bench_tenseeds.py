# Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
# Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
"""The published bench (examples/servo_bench.py) over ten paired noise draws per condition, with the 0.2.8 correction
limit (servo time constant x no-load speed per joint, the same clip for the PID) and, for reference, the 0.2.7 limit
(5 degrees for both). Written 1 Nov 2026 for the 0.2.8 README; a description of the new default, not a test.
Per condition, seeds 101 to 110 (the same draws for both engines and both limits); per draw the ratio
rho = tip rms (syntonic) / tip rms (PID); per condition the geometric mean of rho and its 95 % interval from the t
distribution on log rho (n = 10), as in bench/scripts/flow_v10_seeds.py (whose harness and seeds differ).
Run (from the repository root, about 25 minutes):   python bench/cadence_2026-10-31/servo_bench_tenseeds.py"""
import os, sys, json, time
import numpy as np
from scipy.stats import t as tdist
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..', '..', 'examples'))
import servo_bench as SB

SEEDS = list(range(101, 111))
if __name__ == '__main__':
    t_start = time.time(); out = {}
    for lim, lab in (('authority', '0.2.8'), (np.radians(5.0), '0.2.7')):
        for name in ('A', 'B'):
            for cond in ('steady noise', 'noise step'):
                rows = []
                for s in SEEDS:
                    p = SB.run('pid', name, cond, s, lim); y = SB.run('syntonic', name, cond, s, lim)
                    rows.append(dict(seed=s, tip_pid=p[0], tip_syn=y[0], joint_pid=float(p[2].mean()), joint_syn=float(y[2].mean())))
                lr = np.log([r['tip_syn'] / r['tip_pid'] for r in rows]); m, se = lr.mean(), lr.std(ddof=1) / np.sqrt(len(lr))
                h = tdist.ppf(0.975, len(lr) - 1) * se
                jr = np.exp(np.mean(np.log([r['joint_syn'] / r['joint_pid'] for r in rows])))
                out[f'{lab} {name} {cond}'] = dict(rho=float(np.exp(m)), ci=[float(np.exp(m - h)), float(np.exp(m + h))], joint_ratio=float(jr),
                                                   tip_pid_median=float(np.median([r['tip_pid'] for r in rows])),
                                                   tip_syn_median=float(np.median([r['tip_syn'] for r in rows])), rows=rows)
                v = out[f'{lab} {name} {cond}']
                print(lab, name, cond, 'rho %.3f [%.3f, %.3f]' % (v['rho'], *v['ci']), 'joint x%.3f' % jr,
                      'tip PID %.3f syn %.3f' % (v['tip_pid_median'], v['tip_syn_median']), round(time.time() - t_start), 's', flush=True)
    json.dump(out, open(os.path.join(HERE, 'results', 'servo_bench_tenseeds.json'), 'w'), indent=1)
    print('done', round((time.time() - t_start) / 60, 1), 'min')
