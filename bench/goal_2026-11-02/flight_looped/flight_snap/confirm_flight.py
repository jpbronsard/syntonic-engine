# Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
# Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
"""Pre-registered confirmation: syntonic attitude estimators against tuned Kalman filters on the flight bench, in every
circumstance of a declared battery. Written 2 Nov 2026 before any run of this file except the smoke test; goal set by JP:
"trouve comment la loi peut a la fois etre a parite ou depasser un PID regle et egalement egaler ou depasser un Kalman
regle en toute circonstances. Servo 333 hz max. Test realiste". The servo half is confirmed (confirm_servo.py, S1-S3 pass).

BENCH (frozen, sha256 checked at start): the flight bench of the repository (bench/flight_2026-11-02, amendment 2 noise:
the logged S6 idle/low and S5 half-throttle records of the X500, gyro at 1 kHz with its logged bias, tilt at 250 Hz
mean-removed, looped; one axis, plant at rest so the estimate is the estimation error; 5 s idle pre-roll, arming
calibration = pre-roll gyro mean), driven by a snapshot of the goal/flight harness (common.py, variants*.py, written by
a sub-agent, exploratory; check_harness: C3 noise = the bench's build_noise, run() = the bench's est_only, difference 0).

INFORMATION LEVELS (what each estimator may know before flight)
  I   the corpus bench's: throttle level (a known input), per-throttle noise facts from the records, the idle angle and
      rate random walks, the arming calibration.
  II  I + the per-throttle gyro bias table logged at commissioning (record mean minus idle mean), given to every
      estimator of level II.
RIVALS, tuned against truth (an oracle no real filter has) on training seeds 100-102, schedule C3, before this file:
  KF_rival       the bench's two-state Kalman, fixed R, Q_G, Q_B (results/rival.json: R = 1121 rad^2, Q_G = 0, Q_B =
                 3.12e-7, bias reset at arming; flat ridge, no edge)                                          level I
  KF_sched       the same Kalman with R scheduled on the throttle level, R = m rho s2_rec / DTA (law_flight.KFsched,
                 results/sched_tuned.json: m = 2048, rho_half = 0.354, Q_G = 0, Q_B = 0 (both at the physical limit),
                 J_train 0.0892 against 0.2388 for KF_rival)                                                  level I
  KF_rival_btab  KF_rival + the bias table (tuned without it)                                                level II
  KF_sched_btab  KFsched tuned with the bias table applied (results/sched_tuned.json: m = 131072, Q_G = 0, Q_B =
                 1.56e-7; rho_half ran to 1.05e6, i.e. the half-throttle tilt ignored and the gyro integrated with the
                 in-sample exact table; the surface is flat there, J_train 0.07713 -> 0.07708 over 8 refinements: kept at
                 the last value, the one declared exception to G2)                                             level II
LAWS, nothing tuned against truth; structures chosen on dev seeds 300-309 (declared post hoc, see the notes):
  LAW_I   b4_nm_ic2 (variants3): two-state syntonic law, steady-state Kalman-Bucy gains k2 = sqrt(q_b/s2),
          k1 = sqrt(q_theta/s2 + 2 k2) from the live noise reading s2 with per-throttle memory, q_theta = idle ARW,
          q_b = idle RRW bound at T_rec/4, information-weighted start-up cap, published engine constants    level I
  LAW_II  law_flight.LevelLaw(prior='btab', q_b at T_rec/2): Riccati recursion of (theta, b_idle, b_low, b_half) fed by
          facts only, live per-level noise reading, bias table with its own variance as the level prior      level II
CIRCUMSTANCES (common.CIRC, 60 s flights after the pre-roll): C1 steady low; C2 steady half; C3 half/low every 15 s
  (the training schedule); C4 vibration step (tilt noise x sqrt(10), 20-40 s, low); C5 half/low every 5 s; C6 gyro bias
  step +0.02 deg/s at 30 s, low. Seeds 400-409 (reserved, never used). Score J = rms estimation error (deg) after arming,
  pooled = rms over the 10 seeds; ratio = pooled law / pooled rival.

PREDICTIONS (judged)
 Level I
  F1  parity with the corpus's tuned Kalman in every circumstance: LAW_I / KF_rival <= 1.10 in C1-C6        (dev max 1.03)
  F2  better where dev was clear: LAW_I / KF_rival <= 1.00 in C1, C5, C6                            (dev 0.85, 0.71, 0.47)
  F3  parity with the throttle-scheduled Kalman in C1, C4, C6: LAW_I / KF_sched <= 1.10             (dev 0.97, 0.74, 0.79)
  F4  NOT parity with the throttle-scheduled Kalman in the switching flights: LAW_I / KF_sched > 1.10 in C3 and C5
      (dev 2.01, 1.86: the scheduled filter carries the low-throttle precision through the half-throttle phases; no
      fact-based law tried on dev closed that gap). F4 passing means the goal is NOT met against this rival.
 Level II
  F5  parity with the tuned Kalman given the same bias table: LAW_II / KF_rival_btab <= 1.10 in C1-C6        (dev max 1.08)
  F6  better: LAW_II / KF_rival_btab <= 1.00 in C3, C5, C6                                         (dev 0.66, 0.49, 0.61)
  F7  NOT parity with the scheduled Kalman given the table: LAW_II / KF_sched_btab > 1.10 in C2, C3, C5
      (dev 1.39, 1.59, 1.35). Passing means the goal is NOT met there.
  F8  parity with it in C4 and C6: LAW_II / KF_sched_btab <= 1.10                                  (dev 1.09, 0.64)
Described, not judged: every other ratio (C2 for F3/F4: dev 1.21), wins per seed, KF_sched vs KF_rival.
All dev figures above are against the commissioned rivals (flight_dev/tuned_rivals_dev.jsonl, seeds 300-309).
Gates: G0 module and data hashes; G1 the snapshot reproduces dev values exactly (seed 300: LAW_I and KF_rival as in
goal/flight/results/r3_runs.jsonl, LAW_II as in flight_dev/lvl_runs.jsonl, C1 and C2); G2 the commissioning file is the
one hashed here and no tuned parameter sits on a grid edge (0 is the physical limit, allowed).

Run: python3 confirm_flight.py smoke | gate | run <C1..C6> | summarize
"""
import os, sys, json, time, hashlib
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
SHA = hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest()
REPO = '/home/claude/repo/syntonic-engine'
HASHES = {
    'common.py': 'cfc69c2247f79ed65f75af3806631433f7d9b1af4dc172a3f70533d89dd0ae2a',
    'variants.py': 'e3d097aeb0189078215a25b7ef02e4543fd5e3a77ac8ab28a37351d4ebb29d8b',
    'variants2.py': 'b32a19e04646dec8613145562e642e6b07fd878c228fc6b77f466d447a039eaf',
    'variants3.py': '42b46494480c68186c0950e78034333d9a91f5787b1a79d6b0d12d6620b9910a',
    'law_flight.py': 'e0e9fbaa2039acd91ceae041f0c67c249932dc718df6acf131964812abbbeb78',
    'results/facts2.json': '3390b30ffd6ada79222004d2a4a18f6c4288f8e8a9a71c0db66058467fbbf3ee',
    'results/rival.json': '4d5f95ccd4c69b05ef65d0d65aba00a5327f730f24bcbe2fd15870fcfe64a79f',
    'results/sched_tuned.json': '124d1f8aa7fe5d5fd26b577f9170b6eb68703d758eb6c8da32793e6e72277d24',
    '/home/claude/repo/syntonic-engine/bench/flight_2026-11-02/flight_bench.py': '39e240d07833571ce0a978d29571fe9c635e335755580547161f9fb7f8c918e1',
    '/home/claude/repo/syntonic-engine/bench/flight_2026-11-02/results_a2/facts.json': '4f89486baacdc521330302fdd415b3dcae3b108f68e00fa4b50e938197232c9a',
    '/home/claude/repo/syntonic-engine/bench/flight_2026-11-02/results_a2/tuned_estimators.json': '08185ebcba19c61add4c90bea35f5186cca2051e2c4b752ceb250d1e52ab634e',
    '/home/claude/repo/syntonic-engine/bench/software_2026-10-26/step6g_real_noise.py': 'f3a9ee5b8d967981d2730688ea445329ac0a2c0f0c346c189ae674675d40b476',
    '/home/claude/repo/syntonic-engine/bench/software_2026-10-26/data/sess6.csv.gz': '643780cbb3fc55a5f867052ce621af2194904970cae4f71c47f8b7b503069038',
    '/home/claude/repo/syntonic-engine/bench/software_2026-10-26/data/syntonic_00000000_000000_00005.csv.gz': 'c06bf2822db1921b4558da1bf212f2833f59fc7c1fb7b0b743eebf0ad97d8438',
}
OUT = os.path.join(HERE, 'results_confirm'); os.makedirs(OUT, exist_ok=True)
SEEDS = list(range(400, 410)); CIRCS = ['C1', 'C2', 'C3', 'C4', 'C5', 'C6']


def g0():
    for f, h in HASHES.items():
        p = f if os.path.isabs(f) else os.path.join(HERE, f)
        got = hashlib.sha256(open(p, 'rb').read()).hexdigest()
        assert got == h, f'G0 failed: {f} {got}'
    return len(HASHES)


def build():
    import common as C, variants as V, variants3 as V3, law_flight as L
    ST = json.load(open(os.path.join(C.OUT, 'sched_tuned.json'))); RV = L.RV
    ests, cols = [], []
    ests.append(V.KF(np.array([RV['R']]), np.array([RV['QG']]), np.array([RV['QB']]), 'reset')); cols.append('KF_rival')
    ests.append(L.kf_sched_from(ST['KF_sched'], btab=False)); cols.append('KF_sched')
    ests.append(L.KFbtab(np.array([RV['R']]), np.array([RV['QG']]), np.array([RV['QB']]), 'reset')); cols.append('KF_rival_btab')
    ests.append(L.kf_sched_from(ST['KF_sched_btab'], btab=True)); cols.append('KF_sched_btab')
    e3, c3, _, _ = V3.build(['b4_nm_ic2']); assert c3 == ['b4_nm_ic2']; ests += e3; cols.append('LAW_I')
    ests.append(L.LevelLaw(prior='btab', qb=V.QB2)); cols.append('LAW_II')
    return ests, cols


def fly(seed, circ):
    import common as C
    ests, cols = build(); r = C.run(seed, circ, ests)
    return cols, np.concatenate(r['J']).tolist()


def ratios(P, a, b): return {c: P[c][a] / P[c][b] for c in P}


if __name__ == '__main__':
    MODE = sys.argv[1]
    if MODE == 'smoke':
        os.environ.pop('GOAL_ALLOW_CONFIRM', None)
        cols, J = fly(300, 'C1'); print('smoke C1 seed 300 (dev seed)', dict(zip(cols, np.round(J, 4))))
    elif MODE == 'gate':
        n = g0(); res = {'G0_files': n}
        r3 = {(d['circ'], d['seed']): d for d in map(json.loads, open('/home/claude/work/goal/flight/results/r3_runs.jsonl'))}
        lv = {(d['circ'], d['seed']): d for d in map(json.loads, open('/home/claude/work/goal/confirm/flight_dev/lvl_runs.jsonl'))}
        ok = True
        for circ in ('C1', 'C2'):
            cols, J = fly(300, circ); J = dict(zip(cols, J)); a, b = r3[(circ, 300)], lv[(circ, 300)]
            ref = {'LAW_I': a['J'][a['cols'].index('b4_nm_ic2')], 'KF_rival': a['J'][a['cols'].index('KF_rival')],
                   'LAW_II': b['J'][b['cols'].index('lvl_btab_qb2')]}
            d = {k: J[k] - v for k, v in ref.items()}; res[f'G1_{circ}'] = d; ok &= all(x == 0.0 for x in d.values())
        import common as C
        ST = json.load(open(os.path.join(C.OUT, 'sched_tuned.json')))
        e = {k: v['edges'] for k, v in ST.items()}; res['G2_edges'] = e
        ok &= all(x in ('interior', 'at 0') or (k == 'KF_sched_btab' and q == 'rho_h') for k, v in e.items() for q, x in v.items())
        res['ok'] = bool(ok); json.dump(res, open(os.path.join(OUT, 'gates.json'), 'w'), indent=1); print(res)
    elif MODE == 'run':
        g0(); assert json.load(open(os.path.join(OUT, 'gates.json')))['ok'], 'gates not passed'
        os.environ['GOAL_ALLOW_CONFIRM'] = '1'; circ = sys.argv[2]; path = os.path.join(OUT, f'runs_{circ}.jsonl')
        done = set(json.loads(l)['seed'] for l in open(path)) if os.path.exists(path) else set()
        for s in SEEDS:
            if s in done: continue
            t0 = time.time(); cols, J = fly(s, circ)
            with open(path, 'a') as f: f.write(json.dumps(dict(circ=circ, seed=s, sha=SHA, cols=cols, J=J)) + '\n')
            print(circ, s, ' '.join(f'{c} {j:.4f}' for c, j in zip(cols, J)), f'{time.time() - t0:.0f} s', flush=True)
    elif MODE == 'summarize':
        P, W, complete = {}, {}, True
        for circ in CIRCS:
            p = os.path.join(OUT, f'runs_{circ}.jsonl')
            rows = sorted([json.loads(l) for l in open(p)], key=lambda r: r['seed']) if os.path.exists(p) else []
            if [r['seed'] for r in rows] != SEEDS or any(r['sha'] != SHA for r in rows): complete = False; continue
            cols = rows[0]['cols']; J = np.array([r['J'] for r in rows]); pooled = np.sqrt(np.mean(J ** 2, 0))
            P[circ] = dict(zip(cols, pooled.tolist()))
            W[circ] = {f'{a}<{b}': int((J[:, cols.index(a)] < J[:, cols.index(b)]).sum())
                       for a, b in (('LAW_I', 'KF_rival'), ('LAW_I', 'KF_sched'), ('LAW_II', 'KF_rival_btab'), ('LAW_II', 'KF_sched_btab'))}
        R = {'I_rival': ratios(P, 'LAW_I', 'KF_rival'), 'I_sched': ratios(P, 'LAW_I', 'KF_sched'),
             'II_rival': ratios(P, 'LAW_II', 'KF_rival_btab'), 'II_sched': ratios(P, 'LAW_II', 'KF_sched_btab'),
             'sched_vs_rival': ratios(P, 'KF_sched', 'KF_rival')}
        V = {}
        if complete:
            V['F1'] = all(R['I_rival'][c] <= 1.10 for c in CIRCS); V['F2'] = all(R['I_rival'][c] <= 1.00 for c in ('C1', 'C5', 'C6'))
            V['F3'] = all(R['I_sched'][c] <= 1.10 for c in ('C1', 'C4', 'C6')); V['F4'] = all(R['I_sched'][c] > 1.10 for c in ('C3', 'C5'))
            V['F5'] = all(R['II_rival'][c] <= 1.10 for c in CIRCS); V['F6'] = all(R['II_rival'][c] <= 1.00 for c in ('C3', 'C5', 'C6'))
            V['F7'] = all(R['II_sched'][c] > 1.10 for c in ('C2', 'C3', 'C5')); V['F8'] = all(R['II_sched'][c] <= 1.10 for c in ('C4', 'C6'))
        S = dict(script_sha256=SHA, complete=complete, pooled_deg=P, ratios=R, wins=W, verdict=V)
        json.dump(S, open(os.path.join(OUT, 'summary.json'), 'w'), indent=1); print(json.dumps(S, indent=1))
