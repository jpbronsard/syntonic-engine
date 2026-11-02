# Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
# Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
"""AMENDMENT 2 to flight_bench.py (sha256 39e240d0...) and its amendment 1 (4876f5d0...), written 2 Nov 2026 after the
widened calibration and before any gain tuning or any test-seed run.

What the widened calibration showed: the tuned Kalman went to the far edges again (R = var(tilt) x 2^20, Q_G =
lambda_hat x 2^-10): it ignores the accelerometer and integrates the gyro, and the law's best factor is 4 with a flat
curve from 2.8 to 8. The cause is the noise model, not the estimators: flight_bench.py removed each record's mean from the
gyro, and the half-throttle record (8 s) is looped. A mean-free looped record integrates to a bounded, periodic angle:
the simulated gyro never drifts, so integrating it forever is optimal. A real gyro drifts. The logs show it: the gyro's
mean (its bias) is -0.037 deg/s at idle, -0.036 at low throttle and -0.026 at half throttle (-0.019 and -0.032 in the two
halves of the record): the bias moves with throttle, by about 0.011 deg/s at half throttle, after an arming calibration
made at idle.

Change: the gyro records keep their own means (each throttle level's logged bias); the tilt records keep their means
removed (as step 6g: truth carries the accelerometer's offset). The arming calibration (pre-roll mean, idle) is unchanged,
so after arming the gyro carries the logged bias shift of the throttle level it flies at. Everything else as flight_bench.py
with amendment 1's widened grids: predictions, arms, gain grid, seeds and scores unchanged. Recalibration, tuning and test
follow from here.

Run: python flight_bench_a2.py <repo_root> calibrate | tune <arm> | run [seeds] | summarize
"""
import os, sys, json, time, hashlib
import numpy as np
REPO, MODE = os.path.abspath(sys.argv[1]), sys.argv[2]; _args = list(sys.argv)
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
sys.argv = [_args[0], REPO, 'import']
import flight_bench as F
sys.argv = _args
G = F.G
SHA = hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest()
F.CAL = np.array([0.25, 0.35, 0.5, 0.7, 1.0, 1.4, 2.0, 2.8, 4.0, 5.6, 8.0])
OUT = os.path.join(HERE, 'results_a2'); os.makedirs(OUT, exist_ok=True); F.OUT = OUT


def records_a2():
    if not F._REC:
        for name, seg in (('idle', G.SEGS[0]), ('low', G.SEGS[1]), ('half', G.SEGS[2])):
            tilt, _, g1k = G.load(*seg[:4]); F._REC[name] = (tilt - tilt.mean(), g1k)      # gyro keeps its logged bias
    return F._REC
F.records = records_a2


def lj(name):
    return json.load(open(os.path.join(OUT, name)))


if __name__ == '__main__':
    if MODE == 'calibrate':
        fx = F.facts(); fx['amendment2_sha256'] = SHA; json.dump(fx, open(os.path.join(OUT, 'facts.json'), 'w'), indent=1)
        Rg = fx['var_tilt'] * 2.0 ** np.arange(-8, 21); Qg = fx['lam_hat'] * 2.0 ** np.arange(-10, 9)
        kR = np.repeat(Rg, len(Qg)); kQ = np.tile(Qg, len(Rg)); nk = len(kR)
        Jk = np.mean([F.est_only(s, fx, np.ones(nk, int), kR=kR, kQ=kQ) ** 2 for s in F.TRAIN], 0) ** 0.5
        i = int(np.argmin(Jk))
        Jc = np.mean([F.est_only(s, fx, np.full(len(F.CAL), 2), cfac=F.CAL) ** 2 for s in F.TRAIN], 0) ** 0.5
        j = int(np.argmin(Jc))
        tuned = dict(kf_R=float(kR[i]), kf_QG=float(kQ[i]), kf_R_mult=float(kR[i] / fx['var_tilt']), kf_QG_mult=float(kQ[i] / fx['lam_hat']),
                     kf_J_est=float(Jk[i]), kf_on_edge=bool(kR[i] in (Rg[0], Rg[-1]) or kQ[i] in (Qg[0], Qg[-1])),
                     c_bench=float(F.CAL[j]), c_on_edge=bool(j in (0, len(F.CAL) - 1)),
                     cal_curve={str(c): float(v) for c, v in zip(F.CAL, Jc)}, amendment2_sha256=SHA, script_sha256=F.SHA)
        json.dump(tuned, open(os.path.join(OUT, 'tuned_estimators.json'), 'w'), indent=1); print(json.dumps(fx)); print(json.dumps(tuned)); sys.exit(0)
    fx = lj('facts.json'); tuned = lj('tuned_estimators.json')
    if MODE == 'tune':
        arm = sys.argv[3]; cols = [(arm, g) for g in F.GRID]
        J = np.mean([np.square(F.fly(s, cols, fx, tuned)['J_att']) for s in F.TRAIN], 0) ** 0.5
        i = int(np.argmin(J))
        res = dict(arm=arm, gains=list(F.GRID[i]), J_train=float(J[i]), on_edge=dict(
            Ka=F.GRID[i][0] in (F.KAS[0], F.KAS[-1]), Kp=F.GRID[i][1] in (F.KPS[0], F.KPS[-1]), Kd=F.GRID[i][2] in (F.KDS[0], F.KDS[-1])),
            top5=[(list(F.GRID[j]), float(J[j])) for j in np.argsort(J)[:5]], amendment2_sha256=SHA)
        json.dump(res, open(os.path.join(OUT, f'gains_{arm}.json'), 'w'), indent=1); print(json.dumps(res)); sys.exit(0)
    if MODE == 'run':
        gP = tuple(lj('gains_PID.json')['gains']); gK = tuple(lj('gains_KF.json')['gains']); gS = tuple(lj('gains_SYN1.json')['gains'])
        cols = [('PID', gP), ('KF', gK), ('SYN1', gS), ('SYND', gS), ('SYNB', gS), ('PIDFF', gP), ('SYN1noFF', gS)]
        seeds = [int(s) for s in sys.argv[3:]] if len(sys.argv) > 3 else F.TEST
        for s in seeds:
            t0 = time.time(); r = F.fly(s, cols, fx, tuned)
            row = dict(seed=s, sha=SHA, arms=[c[0] for c in cols], **r)
            with open(os.path.join(OUT, 'test_runs.jsonl'), 'a') as f: f.write(json.dumps(row) + '\n')
            print(s, {c[0]: (round(r['J_att'][i], 3), round(r['J_est'][i], 3)) for i, c in enumerate(cols)}, round(time.time() - t0), 's', flush=True)
        sys.exit(0)
    if MODE == 'summarize':
        rows = [json.loads(l) for l in open(os.path.join(OUT, 'test_runs.jsonl'))]; assert len(rows) == 10, len(rows)
        arms = rows[0]['arms']
        A = {a: np.array([r['J_att'][i] for r in rows]) for i, a in enumerate(arms)}
        E = {a: np.array([r['J_est'][i] for r in rows]) for i, a in enumerate(arms)}
        P = lambda v: float(np.sqrt(np.mean(np.square(v))))
        PA = {a: P(v) for a, v in A.items()}; PE = {a: P(v) for a, v in E.items()}
        V = dict(F1=PE['SYN1'] <= 0.5 * PE['PID'], F2=PE['SYN1'] <= PE['KF'],
                 F3=PA['SYN1'] <= 0.8 * PA['PID'] and int(np.sum(A['SYN1'] < A['PID'])) >= 8,
                 F4=PA['SYN1'] <= PA['KF'] and int(np.sum(A['SYN1'] < A['KF'])) >= 7,
                 K1=PE['SYND'] >= 0.95 * PE['SYN1'],
                 K2=tuned['c_bench'] <= 0.7 and PE['SYNB'] <= 0.9 * PE['SYN1'],
                 JPq=0.8 <= tuned['c_bench'] <= 1.25)
        reg = {a: {r: P([row['att_by_reg'][r][i] for row in rows]) for r in ('low', 'half')} for i, a in enumerate(arms)}
        S = dict(amendment2_sha256=SHA, pooled_att=PA, pooled_est=PE,
                 wins=dict(SYN1_vs_PID=int(np.sum(A['SYN1'] < A['PID'])), SYN1_vs_KF=int(np.sum(A['SYN1'] < A['KF'])),
                           SYNB_vs_KF=int(np.sum(A['SYNB'] < A['KF'])), SYNB_vs_SYN1=int(np.sum(E['SYNB'] < E['SYN1']))),
                 by_regime=reg, jitter={a: P([r['jit'][i] for r in rows]) for i, a in enumerate(arms)},
                 verdict={k: bool(v) for k, v in V.items()}, tuned=tuned,
                 gains={a: lj(f'gains_{a}.json')['gains'] for a in ('PID', 'KF', 'SYN1')})
        json.dump(S, open(os.path.join(OUT, 'summary.json'), 'w'), indent=1); print(json.dumps(S, indent=1))
