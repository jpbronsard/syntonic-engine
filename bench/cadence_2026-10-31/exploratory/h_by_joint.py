"""Exploratory (31 Oct 2026, after JP's question; no prediction): the smoothness H the engine identifies per joint on
the published bench (1 kHz, 3 ms frame, baton arm), steady noise, excerpts A and B, seeds 1 to 5. Do the six joints
share one H, i.e. one class in the 4D field?"""
import sys, numpy as np
sys.path.insert(0, sys.argv[1])
import cadence_bench as CB, omega_phase_rates as OM
REC = []
_step = CB.CadenceEngine.step
def step(self, *a):
    out = _step(self, *a); REC.append((self.H[:6].copy(), self.law[:6].copy())); return out
CB.CadenceEngine.step = step
rows = []
for name in ('A', 'B'):
    for seed in (1, 2, 3, 4, 5):
        REC.clear(); OM.run(name, 'steady', seed, 'syntonic', 6)
        H = np.array([h for h, _ in REC])[3000::10]                     # scored 30 s, 100 Hz
        med = np.median(H, 0); q1, q3 = np.percentile(H, [25, 75], axis=0)
        spread_t = np.median(np.std(H, axis=1))                        # across joints, at each instant
        rows.append((name, seed, med, q1, q3, spread_t, np.mean(H <= 0.051, 0), np.mean(H >= 0.999, 0)))
        print(name, seed, 'H median by joint', np.round(med, 2), 'IQR', np.round(q3 - q1, 2),
              'cross-joint sd (median over time) %.3f' % spread_t, 'at clip low', np.round(np.mean(H <= 0.051, 0), 2),
              'at clip high', np.round(np.mean(H >= 0.999, 0), 2), flush=True)
allmed = np.array([r[2] for r in rows])
print('median over runs, by joint:', np.round(np.median(allmed, 0), 3))
print('A:', np.round(np.median(allmed[:5], 0), 3), ' B:', np.round(np.median(allmed[5:], 0), 3))
