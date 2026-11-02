# Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
# Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
"""What the C++/Rust option spread costs on the realistic MAESTRO 1 bench (0.3.1; described, not pre-registered).

With spread = m the loop law's new gains take effect m ticks after the statistics they come from, instead of in the same
tick. This runs the law with m = 0 (the Python law), 250 and 499 side by side on the realistic bench of
bench/goal_2026-11-02/servo (baton arm, 6 joints, servo delay 2.5 ms as in the confirmation), nominal and with the
encoder noise x10 from 10 to 20 s, excerpts A and B, development seeds 91-93 (not the confirmation's seeds 71-75).
Paired: the three columns see the same noise. The reference of the option is SpreadLoopLaw in
rust/tools/make_law_replay_log.py, which the C++ and Rust laws reproduce (cpp/tests/replay_laws.cpp).
Writes spread_effect.json. About 5 minutes.   python bench/port_2026-11-02/spread_effect.py
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
SPREADS = (0, 250, 499)


def factory(m):
    kw = dict(dt=RB.DT, frame=RB.FRAME, servo_tau=RB.TS_NOM, dead=DEAD)
    if m == 0:
        return lambda q0, w0: SyntonicLoopLaw(q0, w0, **kw)
    return lambda q0, w0: SpreadLoopLaw(q0, w0, spread=m, **kw)


if __name__ == '__main__':
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
    json.dump(dict(rows=rows, summary=summary), open(os.path.join(HERE, 'spread_effect.json'), 'w'), indent=1)
