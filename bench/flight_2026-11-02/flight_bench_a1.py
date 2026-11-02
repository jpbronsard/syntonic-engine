# Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
# Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
"""AMENDMENT 1 to flight_bench.py (sha256 39e240d0..., frozen), written 2 Nov 2026 after its calibration step and before
any gain tuning or any test-seed run.

What the calibration showed (training seeds 100-102, estimation error only): the tuned Kalman landed on the edges of
its pre-registered grid (R = var(tilt) x 256, the top; Q_G = lambda_hat x 0.25, the bottom), and the law's calibrated
factor landed on the top of its grid (c = 2, the error still falling: 0.32 deg at c = 1, 0.24 at c = 2). An optimum on
a grid edge is not an optimum: the Kalman would be under-tuned, which would flatter the syntonic controller, and the
calibrated factor would be truncated. The grids are therefore widened, and nothing else changes:
  Kalman: R in var(tilt) x 2^(-8..20), Q_G in lambda_hat x 2^(-10..8).
  Law calibration: c in {0.25, 0.35, 0.5, 0.7, 1, 1.4, 2, 2.8, 4, 5.6, 8}.
Predictions, arms, gain grid, seeds and scores: unchanged (read them in flight_bench.py). The calibration already seen
makes K2 (c_bench <= 0.7) and JPq (c_bench in [0.8, 1.25]) fail on the training seeds; they stay as written and are
judged on the widened calibration.

Run: python flight_bench_a1.py <repo_root> calibrate   (then flight_bench.py tune / run / summarize, unchanged; they read
results/tuned_estimators.json written here)
"""
import os, sys, json, time, hashlib
import numpy as np
REPO, MODE = os.path.abspath(sys.argv[1]), sys.argv[2]; _args = list(sys.argv)
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
sys.argv = [_args[0], REPO, 'import']
import flight_bench as F
sys.argv = _args
SHA = hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest()
F.CAL = np.array([0.25, 0.35, 0.5, 0.7, 1.0, 1.4, 2.0, 2.8, 4.0, 5.6, 8.0])
OUT = F.OUT

if __name__ == '__main__':
    fx = F.load_json('facts.json')
    if MODE == 'calibrate':
        Rg = fx['var_tilt'] * 2.0 ** np.arange(-8, 21); Qg = fx['lam_hat'] * 2.0 ** np.arange(-10, 9)
        kR = np.repeat(Rg, len(Qg)); kQ = np.tile(Qg, len(Rg)); nk = len(kR)
        Jk = np.mean([F.est_only(s, fx, np.ones(nk, int), kR=kR, kQ=kQ) ** 2 for s in F.TRAIN], 0) ** 0.5
        i = int(np.argmin(Jk))
        Jc = np.mean([F.est_only(s, fx, np.full(len(F.CAL), 2), cfac=F.CAL) ** 2 for s in F.TRAIN], 0) ** 0.5
        j = int(np.argmin(Jc))
        tuned = dict(kf_R=float(kR[i]), kf_QG=float(kQ[i]), kf_R_mult=float(kR[i] / fx['var_tilt']), kf_QG_mult=float(kQ[i] / fx['lam_hat']),
                     kf_J_est=float(Jk[i]), kf_on_edge=bool(kR[i] in (Rg[0], Rg[-1]) or kQ[i] in (Qg[0], Qg[-1])),
                     c_bench=float(F.CAL[j]), c_on_edge=bool(j in (0, len(F.CAL) - 1)),
                     cal_curve={str(c): float(v) for c, v in zip(F.CAL, Jc)}, amendment_sha256=SHA, script_sha256=F.SHA)
        json.dump(tuned, open(os.path.join(OUT, 'tuned_estimators.json'), 'w'), indent=1); print(json.dumps(tuned)); sys.exit(0)
