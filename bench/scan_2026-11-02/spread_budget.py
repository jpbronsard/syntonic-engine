# Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
# Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
"""What the spreads of 0.3.2 cost on the realistic MAESTRO 1 bench (described, not pre-registered), as
bench/port_2026-11-02/spread_effect.py measured 250 and 499 for 0.3.1.

The spreads come from cost.json (cost_model.py), two joints, servo delay 2.5 ms (878 robust candidates):
  365  the smallest spread that keeps the recomputation within 0.45 ms of every tick, counting no pruning
       (guaranteed; no spread below 500 keeps it within 0.25 ms)
  188  and 344: the spreads for 0.45 and 0.25 ms if pruning were counted on at the largest share seen on the
       two-joint log (51.5 %); not guaranteed, not implemented, for information
The law with spread = m: the same gains as spread = 0, taking effect m ticks later (SpreadLoopLaw of
rust/tools/make_law_replay_log.py, which the C++ and Rust laws reproduce with and without the speed-ups). The bench is
that of spread_effect.py: baton arm, 6 joints, nominal and encoder noise x10 from 10 to 20 s, excerpts A and B,
development seeds 91-93, the columns paired (the same noise). Writes spread_budget.json. About 8 minutes.
   python bench/scan_2026-11-02/spread_budget.py
"""
import os, sys, json, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..'))
os.environ.setdefault('SYNTONIC_REPO', REPO)
sys.path.insert(0, os.path.join(REPO, 'python'))
sys.path.insert(0, os.path.join(REPO, 'rust', 'tools'))
sys.path.insert(0, os.path.join(REPO, 'bench', 'goal_2026-11-02', 'servo', 'bench'))
import realistic_bench as RB  # noqa: E402
from make_law_replay_log import SpreadLoopLaw  # noqa: E402
from syntonic_engine import SyntonicLoopLaw  # noqa: E402

DEAD = 0.0025
COST = json.load(open(os.path.join(HERE, 'cost.json')))
SP = [c for c in COST['cases'] if c['dead_ms'] == 2.5][0]['spreads']
GUAR = SP['0.45 ms']['after']
CRED = (SP['0.45 ms']['after_if_pruning_counted'], SP['0.25 ms']['after_if_pruning_counted'])
SPREADS = (0, GUAR) + CRED


def factory(m):
    kw = dict(dt=RB.DT, frame=RB.FRAME, servo_tau=RB.TS_NOM, dead=DEAD)
    if m == 0:
        return lambda q0, w0: SyntonicLoopLaw(q0, w0, **kw)
    return lambda q0, w0: SpreadLoopLaw(q0, w0, spread=m, **kw)


if __name__ == '__main__':
    print('spreads', SPREADS, flush=True)
    rows = []
    for circ, mod in (('nominal', {}), ('noise_step', dict(noise='step'))):
        for name in 'AB':
            for seed in (91, 92, 93):
                t = time.time()
                r = RB.simulate(name, dict(dead_time=DEAD, **mod), seed, [(f'S{m}', factory(m)) for m in SPREADS])
                row = dict(circ=circ, ex=name, seed=seed, **{k: r[k]['tip_mm'] for k in r})
                rows.append(row); print(json.dumps(row), f'{time.time() - t:.0f} s', flush=True)
    summary = {}
    for circ in ('nominal', 'noise_step'):
        R = [r for r in rows if r['circ'] == circ]
        for m in SPREADS[1:]:
            rat = np.array([r[f'S{m}'] / r['S0'] for r in R])
            summary[f'{circ} spread {m}'] = dict(mean=float(rat.mean()), lo=float(rat.min()), hi=float(rat.max()),
                                                 worse=int(np.sum(rat > 1)), n=len(rat))
            print(f'{circ}: spread {m}: tip error ratio to spread 0, mean {rat.mean():.4f}, '
                  f'range {rat.min():.4f} to {rat.max():.4f}, worse in {np.sum(rat > 1)} of {len(rat)}')
    json.dump(dict(spreads=SPREADS, rows=rows, summary=summary), open(os.path.join(HERE, 'spread_budget.json'), 'w'),
              indent=1)
