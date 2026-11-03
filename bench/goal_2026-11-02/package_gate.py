# Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
# Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
"""Gate for the two laws added to the Python package in 0.3.0 (written 2 Nov 2026 after both confirmations; a check, not
a prediction). Run inside the frozen harnesses of this folder, the package's classes must reproduce the recorded
confirmation runs exactly (difference 0.0):

  servo   syntonic_engine.SyntonicLoopLaw in place of the frozen syntonic_loop_law.py, with the two tuned PIDs, as
          confirm_servo.py builds them: servo/results_servo/runs_C1.jsonl, key 'C1 A 71' (all three columns)
  flight  syntonic_engine.RiccatiAttitudeLaw in place of ricv2.RicV, its facts read from the fact files:
          flight_realistic/results_confirm3/runs_C1.jsonl and runs_C3.jsonl, seeds 500 and 501, columns LAW_R (q_b at
          T_rec/2, the judged law) and LAW_Rqb4 (q_b at T_rec/4, described)
  all     in addition the LAW column of all 90 servo runs (the law alone in the harness) and LAW_R, LAW_Rqb4 in all
          60 flights of confirmation 3 (about 25 minutes on one core)

  --scan  (0.3.2) the servo runs with syntonic_engine.ScanLoopLaw in place of SyntonicLoopLaw: the gains found as the
          C++ and Rust laws find them, with the two speed-ups of 0.3.2 (pruning, the quadratic form). They must
          reproduce the same recorded runs exactly; prints the share of the points evaluated and the fallbacks.

The flight harness names absolute paths of the working session: run replay_setup.sh first (see NOTE.md). The servo
harness needs scipy. Seeds 500-509 are reserved for confirmation 3, which is complete; this check sets
GOAL_ALLOW_CONFIRM3=1 for itself only.
Run: python3 package_gate.py [servo | flight | all] [--scan]
"""
import os, sys, json, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, os.path.join(REPO, 'python'))
from syntonic_engine import SyntonicLoopLaw, RiccatiAttitudeLaw, ScanLoopLaw  # noqa: E402

SCAN = '--scan' in sys.argv
ARGS = [a for a in sys.argv[1:] if a != '--scan']
MODE = ARGS[0] if ARGS else 'default'
OUT = {}
LAWS = []                                      # with --scan: the ScanLoopLaw instances, for their statistics


def loop_law(q0, w0, **kw):
    if not SCAN:
        return SyntonicLoopLaw(q0, w0, **kw)
    law = ScanLoopLaw(q0, w0, **kw); LAWS.append(law)
    return law


def report(name, got, ref):
    d = float(got) - float(ref); OUT[name] = d
    print(f'{name:28s} {got!r:24s} recorded {ref!r:24s}', 'equal' if got == ref else f'DIFF {d:.3e}', flush=True)


def servo(all_runs):
    os.environ.setdefault('SYNTONIC_REPO', REPO)
    sys.path.insert(0, os.path.join(HERE, 'servo', 'bench'))
    import realistic_bench as RB, rivals as RV
    S = os.path.join(HERE, 'servo')
    T = json.load(open(os.path.join(S, 'results_servo', 'tuned.json')))
    law = ('LAW', lambda q0, w0: loop_law(q0, w0, dt=RB.DT, frame=RB.FRAME, servo_tau=RB.TS_NOM, dead=0.0025))
    ctr = [('PIDLM', {n: RV.pid_factory(np.array(T['leadmean']['gains'][n]), lead=True, p_mean=True) for n in 'AB'}),
           ('PIDL', {n: RV.pid_factory(np.array(T['lead']['gains'][n]), lead=True, p_mean=False) for n in 'AB'}), law]
    rows = {}
    for c in [f'C{i}' for i in range(1, 10)]:
        for line in open(os.path.join(S, 'results_servo', f'runs_{c}.jsonl')):
            d = json.loads(line); rows[d['key']] = d
    COMM = dict(dead_time=0.0025)
    CIRC = {'C1': {}, 'C2': dict(noise='step'), 'C3': dict(noise='bursts'), 'C4': dict(dead_time=0.0035),
            'C5': dict(dead_time=0.0015), 'C6': dict(servo_tau=0.035), 'C7': dict(servo_tau=0.018),
            'C8': dict(load_scale=1.2), 'C9': dict(cmd_quant_deg=0.09)}               # as in confirm_servo.py
    t0 = time.time(); r = RB.simulate('A', dict(COMM), 71, ctr, allow_reserved=True); ref = rows['C1 A 71']
    for k in ('PIDLM', 'PIDL', 'LAW'):
        for m in ('tip_mm', 'jitter_mm', 'joint_deg'):
            report(f'servo C1 A 71 {k} {m}', r[k][m], ref[k][m])
    print(f'  {time.time() - t0:.0f} s', flush=True)
    if all_runs:
        for c, mod in CIRC.items():
            for name in 'AB':
                for s in (71, 72, 73, 74, 75):
                    r = RB.simulate(name, dict(COMM, **mod), s, [law], allow_reserved=True)
                    report(f'servo {c} {name} {s} LAW tip_mm', r['LAW']['tip_mm'], rows[f'{c} {name} {s}']['LAW']['tip_mm'])


def flight(all_runs):
    snap = '/home/claude/work/goal/confirm/flight_snap'
    if not os.path.isdir(snap):
        sys.exit(f'{snap} not found: run replay_setup.sh first (see NOTE.md)')
    os.environ['GOAL_ALLOW_CONFIRM3'] = '1'
    F3 = os.path.join(HERE, 'flight_realistic')
    sys.path.insert(0, F3)
    import bench_bb as B
    B.install()
    F2 = json.load(open(os.path.join(HERE, 'flight_looped', 'flight_snap', 'results', 'facts2.json')))
    FX = json.load(open(os.path.join(REPO, 'bench', 'flight_2026-11-02', 'results_a2', 'facts.json')))
    levels = {r: dict(s2=F2['curves']['raw'][r]['s2_rec'], lam=FX['lam'][r], T_rec=F2['curves']['raw'][r]['T_rec'])
              for r in ('idle', 'low', 'half')}
    assert F2 == B.V.F2 and FX['lam'] == B.C.FX['lam'], 'the fact files differ from those of the harness'
    law = lambda qb: RiccatiAttitudeLaw(levels, q_b=qb, T_pre=B.C.T_PRE, dt=B.C.DT, dt_accel=B.C.DTA)
    rows = {}
    for c in ('C1', 'C2', 'C3', 'C4', 'C5', 'C6'):
        for line in open(os.path.join(F3, 'results_confirm3', f'runs_{c}.jsonl')):
            d = json.loads(line); rows[(c, d['seed'])] = dict(zip(d['cols'], d['J']))
    flights = [(c, s) for c in ('C1', 'C3') for s in (500, 501)]
    if all_runs:
        flights = [(c, s) for c in ('C1', 'C2', 'C3', 'C4', 'C5', 'C6') for s in range(500, 510)]
    for c, s in flights:
        r = B.run4(s, c, [law(F2['q_b_half']), law(F2['q_b'])])
        for name, j in zip(('LAW_R', 'LAW_Rqb4'), r['J']):
            report(f'flight {c} {s} {name}', float(np.atleast_1d(j)[0]), rows[(c, s)][name])


if __name__ == '__main__':
    t = time.time()
    if MODE in ('default', 'servo', 'all'):
        servo(MODE == 'all')
    if MODE in ('default', 'flight', 'all'):
        flight(MODE == 'all')
    if SCAN and LAWS:
        st = [x for law in LAWS for x in law.stats]
        print(f'ScanLoopLaw: {len(LAWS)} law(s), {len(st)} recomputations, points evaluated '
              f'{np.mean([x["points"] / x["total"] for x in st]):.1%} (largest {max(x["points"] / x["total"] for x in st):.1%}), '
              f'fallbacks {sum(law.fallbacks for law in LAWS)}, smallest gap best-next J seen {min(x["margin"] for x in st):.2e}')
    worst = max(abs(v) for v in OUT.values())
    print(f'{len(OUT)} values, largest difference {worst:.3e}:', 'PASS (exact)' if worst == 0.0 else 'FAIL',
          f'({time.time() - t:.0f} s)')
    sys.exit(0 if worst == 0.0 else 1)
