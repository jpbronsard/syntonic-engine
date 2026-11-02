# Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
# Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
"""Exploratory, 2 Nov 2026, after the verdicts of pid_lead_fair.py (not pre-registered): does the PID's advantage come from
being tuned on the excerpt it is tested on? PIDD with the gains tuned on the other excerpt (A's gains on B, B's on A),
same test seeds and cells. Run: python explore_pid_cross.py <repo_root> <excerpt> <steady|step>"""
import os, sys, json, importlib.util
import numpy as np
ROOT, NAME, COND = os.path.abspath(sys.argv[1]), sys.argv[2], sys.argv[3]
HERE = os.path.dirname(os.path.abspath(__file__))
sys.argv = [sys.argv[0], ROOT, 'import']
spec = importlib.util.spec_from_file_location('P', os.path.join(HERE, 'pid_lead_fair.py')); P = importlib.util.module_from_spec(spec); spec.loader.exec_module(P)
T = json.load(open(os.path.join(P.OUT, 'tuned_gains.json')))
other = 'B' if NAME == 'A' else 'A'
arms = [P.pid_arm('PIDD_cross', T['PIDD']['tuned'][other]['gains'], P.LEAD)]
for seed in P.TEST_SEEDS:
    r = P.simulate(NAME, COND, seed, P.DT_F, P.FR_F, P.BATON, arms, []); r.pop('pid_joint_deg')
    with open(os.path.join(P.OUT, 'explore_cross.jsonl'), 'a') as f: f.write(json.dumps(dict(key=f'{NAME} {COND} {seed}', **r)) + '\n')
    print(NAME, COND, seed, round(r['PIDD_cross']['tip_mm'], 3), flush=True)
