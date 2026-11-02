# Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
# Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
"""Pre-registered confirmation 3: the Riccati form of the syntonic law against tuned Kalman filters on the realistic flight
bench, every circumstance. Written 2 Nov 2026 before any run of this file except the smoke test, after confirmation 2
(confirm_flight.py, looped bench, seeds 400-409) FAILED its parity predictions; that result stands. Goal set by JP:
"trouve comment la loi peut a la fois etre a parite ou depasser un PID regle et egalement egaler ou depasser un Kalman
regle en toute circonstances. Servo 333 hz max. Test realiste" (servo half confirmed: confirm_servo.py, S1-S3 pass).

BENCH = amendments 3-4 (bench_bb.py, frozen): the repository's flight bench (logged X500 records, one axis, plant at rest,
5 s idle pre-roll, arming calibration) with (3) block-bootstrap replay instead of looping (validated: looping made 60 s
tilt means 5-8x too precise at half throttle and the gyro integral 8.6x too small; the bootstrap is within 1.0-1.5x of the
records' white / random-walk extrapolation) and (4) an unknown initial attitude theta0 ~ N(0, 0.5 deg) per flight.
RIVALS, tuned against truth (an oracle) on training seeds 100-102, schedule C3, before this file (rivals_bb.py, frozen):
  KF_rival3  the bench's two-state Kalman, fixed R, Q_G, Q_B (R ran up a flat ridge: J_train 0.426247 constant to 6 digits)
  KF_sched3  R scheduled on the throttle level, R = m rho s2_rec / DTA (m and Q_B ran down a flat ridge, J 0.151557)
  KF_full3   R and Q_G scheduled (converged to the same filter as KF_sched3: Q_G = 0)
  The ridge runaways are the declared exception to G2: the training J is constant along them.
LAW_R (ricv2.RicV, one column, nothing tuned against truth; structure chosen on dev seeds 300-309, declared post hoc:
  the sequence was ricv v1 -> v2 options sa and rm -> q_b bound, see r3/dev_r4*.jsonl):
  state (theta, b_idle, b_low, b_half); the active level's bias drives the gyro integration and random-walks at
  q_b = 3 sigma_A^2(T)/T at T = T_rec/2 of the idle record (9.23e-11 rad^2/s^3, the tighter of the two published bounds);
  q_theta(level) = lam(level) (the record's angle random walk at that throttle); z = tilt, R = max(s2_rec(level),
  live reading) / DTA (published Hadamard pipeline with per-level memory initialised at the fact); stationary alignment
  in the pre-roll (theta from the tilt alone, no gyro integration: the vehicle is on the ground, its gyro being
  calibrated); arming: biases 0, bias covariance (lam_idle / T_pre) ones + diag(lam(level) / T_rec(level)).
DESCRIBED (not judged): LAW_I (b4_nm_ic2, the steady-state law of confirmation 2), LAW_Rqb4 (LAW_R with q_b at T_rec/4).
CIRCUMSTANCES (60 s): C1 steady low; C2 steady half; C3 half/low every 15 s; C4 vibration step (tilt noise x sqrt(10),
20-40 s, low); C5 half/low every 5 s; C6 gyro bias step +0.02 deg/s at 30 s, low. Seeds 500-509 (reserved, never used).
Score: rms estimation error (deg) after arming; pooled = rms over the 10 seeds; ratio = pooled LAW_R / pooled rival.

PREDICTIONS (judged)
 R1  the goal: parity or better with every tuned Kalman in every circumstance: LAW_R / KF <= 1.10 for KF in {KF_rival3,
     KF_sched3, KF_full3} and C1-C6                                    (dev max 0.891, all bootstrap upper bounds <= 1.05)
 R2  better: LAW_R / KF_rival3 <= 1.00 in C1, C2, C3, C5, C6 and LAW_R / KF_sched3 <= 1.00 in C1, C2, C4, C5, C6
     (dev 0.85, 0.60, 0.25, 0.23, 0.56 and 0.55, 0.76, 0.36, 0.82, 0.89; not claimed: C4 vs KF_rival3 0.89 [0.76, 1.05],
     C3 vs KF_sched3 0.88 [0.74, 1.04])
 R3  the memory is what matters: the steady-state law LAW_I is NOT at parity with KF_sched3 in C3 and C5:
     LAW_I / KF_sched3 > 1.10 in both                                                                   (dev 2.17, 2.24)
Gates: G0 hashes; G1 seed 300, C1 and C4: every column reproduces r3/dev_r4c_runs.jsonl exactly (LAW_I and the rivals
also r3/dev_r4b_runs.jsonl); G2 rival edges interior, at 0, or on the declared flat ridges.
Run: python3 confirm_flight3.py smoke | gate | run <C1..C6> | summarize
"""
import os, sys, json, time, hashlib
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
SHA = hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest()
HASHES = {
    '/home/claude/work/goal/r3/bench_bb.py': '9d36615a6b231aeb2cf77b9d975b9c3d4ab0baca98644847d0290a4fbc4eca09',
    '/home/claude/work/goal/r3/rivals_bb.py': '3e68f31510fb83a2e27438e612fc5bee251b233b848ec774d4c1d1570699f88d',
    '/home/claude/work/goal/r3/ricv2.py': '3a32772c8b4366ddce94339118fced400154953e07b41b9f427492ecba594588',
    '/home/claude/work/goal/r3/results/rivals_r4_rival_sched.json': 'f8441e6c8053b0b87a2836114b76327755e52328708be4ab27949e376a8c0cf3',
    '/home/claude/work/goal/r3/results/rivals_r4_full.json': '38e1fa12b3dd48b28ffd3602d904e996b380c3cf7587a533d42c941ec4dee530',
    '/home/claude/work/goal/r3/dev_r4b_runs.jsonl': '8716e1bb3ce66cbe502004a8d8a343e5e52b7b45422aef864fa6073ac8488fde',
    '/home/claude/work/goal/r3/dev_r4c_runs.jsonl': 'e6dfb7b72ac7e3c2c0f567c0cc8acf64f2423df95167d261d7c24b8ecb2ff76b',
    '/home/claude/work/goal/confirm/flight_snap/common.py': 'cfc69c2247f79ed65f75af3806631433f7d9b1af4dc172a3f70533d89dd0ae2a',
    '/home/claude/work/goal/confirm/flight_snap/variants.py': 'e3d097aeb0189078215a25b7ef02e4543fd5e3a77ac8ab28a37351d4ebb29d8b',
    '/home/claude/work/goal/confirm/flight_snap/variants2.py': 'b32a19e04646dec8613145562e642e6b07fd878c228fc6b77f466d447a039eaf',
    '/home/claude/work/goal/confirm/flight_snap/variants3.py': '42b46494480c68186c0950e78034333d9a91f5787b1a79d6b0d12d6620b9910a',
    '/home/claude/work/goal/confirm/flight_snap/law_flight.py': 'e0e9fbaa2039acd91ceae041f0c67c249932dc718df6acf131964812abbbeb78',
    '/home/claude/work/goal/confirm/flight_snap/results/facts2.json': '3390b30ffd6ada79222004d2a4a18f6c4288f8e8a9a71c0db66058467fbbf3ee',
    '/home/claude/work/goal/confirm/flight_snap/results/rival.json': '4d5f95ccd4c69b05ef65d0d65aba00a5327f730f24bcbe2fd15870fcfe64a79f',
    '/home/claude/repo/syntonic-engine/bench/flight_2026-11-02/flight_bench.py': '39e240d07833571ce0a978d29571fe9c635e335755580547161f9fb7f8c918e1',
    '/home/claude/repo/syntonic-engine/bench/flight_2026-11-02/results_a2/facts.json': '4f89486baacdc521330302fdd415b3dcae3b108f68e00fa4b50e938197232c9a',
    '/home/claude/repo/syntonic-engine/bench/flight_2026-11-02/results_a2/tuned_estimators.json': '08185ebcba19c61add4c90bea35f5186cca2051e2c4b752ceb250d1e52ab634e',
    '/home/claude/repo/syntonic-engine/bench/software_2026-10-26/step6g_real_noise.py': 'f3a9ee5b8d967981d2730688ea445329ac0a2c0f0c346c189ae674675d40b476',
    '/home/claude/repo/syntonic-engine/bench/software_2026-10-26/data/sess6.csv.gz': '643780cbb3fc55a5f867052ce621af2194904970cae4f71c47f8b7b503069038',
    '/home/claude/repo/syntonic-engine/bench/software_2026-10-26/data/syntonic_00000000_000000_00005.csv.gz': 'c06bf2822db1921b4558da1bf212f2833f59fc7c1fb7b0b743eebf0ad97d8438',
}
OUT = os.path.join(HERE, 'results_confirm3'); os.makedirs(OUT, exist_ok=True)
SEEDS = list(range(500, 510)); CIRCS = ['C1', 'C2', 'C3', 'C4', 'C5', 'C6']; KFS = ['KF_rival3', 'KF_sched3', 'KF_full3']


def g0():
    for f, h in HASHES.items():
        got = hashlib.sha256(open(f, 'rb').read()).hexdigest(); assert got == h, f'G0 failed: {f} {got}'
    return len(HASHES)


def build():
    import rivals_bb as RB, bench_bb as B, ricv2 as R, variants3 as V3
    V = B.V; rp = {}
    for f in ('rivals_r4_rival_sched.json', 'rivals_r4_full.json'):
        for k, v in json.load(open(os.path.join(HERE, 'results', f))).items(): rp[k] = v['params']
    ests, cols = [], []
    for kind in ('rival', 'sched', 'full'):
        ests.append(RB.make(kind, np.array([rp[kind]]))); cols.append(f'KF_{kind}3')
    ests.append(R.RicV([dict(k=3, qt='level', prior='vb', qb=V.QB2, sa=True, rm='max')])); cols.append('LAW_R')
    e3, c3, _, _ = V3.build(['b4_nm_ic2']); ests += e3; cols.append('LAW_I')
    ests.append(R.RicV([dict(k=3, qt='level', prior='vb', qb=V.QB4, sa=True, rm='max')])); cols.append('LAW_Rqb4')
    return ests, cols, B


def fly(seed, circ):
    ests, cols, B = build(); r = B.run4(seed, circ, ests)
    return cols, np.concatenate([np.atleast_1d(j) for j in r['J']]).tolist()


if __name__ == '__main__':
    MODE = sys.argv[1]
    if MODE == 'smoke':
        cols, J = fly(300, 'C1'); print('smoke C1 seed 300 (dev)', dict(zip(cols, np.round(J, 4))))
    elif MODE == 'gate':
        res = {'G0_files': g0()}; ok = True
        dc = {(d['circ'], d['seed']): d for d in map(json.loads, open(os.path.join(HERE, 'dev_r4c_runs.jsonl')))}
        db = {(d['circ'], d['seed']): d for d in map(json.loads, open(os.path.join(HERE, 'dev_r4b_runs.jsonl')))}
        for circ in ('C1', 'C4'):
            cols, J = fly(300, circ); J = dict(zip(cols, J)); c, b = dc[(circ, 300)], db[(circ, 300)]
            ref = {'KF_rival3': c['J'][c['cols'].index('KF_rival3')], 'KF_full3': c['J'][c['cols'].index('KF_full3')],
                   'LAW_R': c['J'][c['cols'].index('L3_vb_sa_max_qb2')], 'LAW_Rqb4': c['J'][c['cols'].index('L3_vb_sa_max')],
                   'LAW_I': b['J'][b['cols'].index('LAW_I')]}
            d = {k: J[k] - v for k, v in ref.items()}; d['sched_minus_full'] = J['KF_sched3'] - J['KF_full3']
            res[f'G1_{circ}'] = d; ok &= all(x == 0.0 for x in d.values())
        E = {}
        for f in ('rivals_r4_rival_sched.json', 'rivals_r4_full.json'):
            for k, v in json.load(open(os.path.join(HERE, 'results', f))).items():
                hist = [h['J'] for h in v['history']]; flat = (max(hist[-4:]) - min(hist[-4:])) / min(hist[-4:]) < 1e-5
                E[k] = dict(edges=v['edges'], flat_last4=bool(flat)); ok &= all(e in ('interior', 'at 0') for e in v['edges']) or flat
        res['G2'] = E; res['ok'] = bool(ok); json.dump(res, open(os.path.join(OUT, 'gates.json'), 'w'), indent=1); print(res)
    elif MODE == 'run':
        g0(); assert json.load(open(os.path.join(OUT, 'gates.json')))['ok'], 'gates not passed'
        os.environ['GOAL_ALLOW_CONFIRM3'] = '1'; circ = sys.argv[2]; path = os.path.join(OUT, f'runs_{circ}.jsonl')
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
            cols = rows[0]['cols']; J = np.array([r['J'] for r in rows]); P[circ] = dict(zip(cols, np.sqrt(np.mean(J ** 2, 0)).tolist()))
            W[circ] = {f'{a}<{b}': int((J[:, cols.index(a)] < J[:, cols.index(b)]).sum()) for a in ('LAW_R', 'LAW_I') for b in KFS}
        R = {f'{a}/{b}': {c: P[c][a] / P[c][b] for c in P} for a in ('LAW_R', 'LAW_I', 'LAW_Rqb4') for b in KFS}
        V = {}
        if complete:
            V['R1'] = all(R[f'LAW_R/{k}'][c] <= 1.10 for k in KFS for c in CIRCS)
            V['R2'] = (all(R['LAW_R/KF_rival3'][c] <= 1.00 for c in ('C1', 'C2', 'C3', 'C5', 'C6'))
                       and all(R['LAW_R/KF_sched3'][c] <= 1.00 for c in ('C1', 'C2', 'C4', 'C5', 'C6')))
            V['R3'] = all(R['LAW_I/KF_sched3'][c] > 1.10 for c in ('C3', 'C5'))
        S = dict(script_sha256=SHA, complete=complete, pooled_deg=P, ratios=R, wins=W, verdict=V)
        json.dump(S, open(os.path.join(OUT, 'summary.json'), 'w'), indent=1); print(json.dumps(dict(verdict=V, complete=complete), indent=1))
