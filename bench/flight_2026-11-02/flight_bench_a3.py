# Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
# Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
"""AMENDMENT 3 to flight_bench.py (39e240d0...) after amendments 1 (4876f5d0...) and 2 (c89b386c...), written 2 Nov 2026
after the gain tuning of amendment 2 and before any test-seed run.

What the tuning showed (training seeds 100-102): all three controllers took the top of the pre-registered gain grid
(Ka = 16 /s, rate Kp = 80 / B; PID also the top Kd). Gains on the edge of their grid are not tuned; the comparison would
be between under-tuned controllers. The grid is widened, and nothing else changes:
  Ka in {8, 16, 32, 64} /s; rate Kp in {40, 80, 160, 320} / B; rate Kd in {0, 2e-4, 5e-4, 1e-3}; rate Ki in {0, 2 Kp}.
Everything else is amendment 2 (noise model with the logged gyro bias, widened estimator grids, its calibration:
Kalman R = var(tilt) x 8192 and Q_G at the bottom edge, about zero; law factor c_bench = 2.8). Predictions unchanged.

Run: python flight_bench_a3.py <repo_root> tune <arm> | run [seeds] | summarize   (results in results_a2/, gains
overwritten by this amendment's tuning)
"""
import os, sys, runpy
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
_args = list(sys.argv); sys.argv = [_args[0], _args[1], 'import']
import flight_bench as F
F.KAS = [8.0, 16.0, 32.0, 64.0]; F.KPS = [40.0 / F.B, 80.0 / F.B, 160.0 / F.B, 320.0 / F.B]; F.KDS = [0.0, 2e-4, 5e-4, 1e-3]; F.KIM = [0.0, 2.0]
F.GRID = [(ka, kp, kd, km) for ka in F.KAS for kp in F.KPS for kd in F.KDS for km in F.KIM]
sys.argv = _args
runpy.run_path(os.path.join(HERE, 'flight_bench_a2.py'), run_name='__main__')
