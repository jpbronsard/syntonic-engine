"""Check: build_noise_c(seed, 'C3') == F.build_noise(seed), and run() == F.est_only() on the same columns (seed 100,
a training seed). Prints max differences and timings. No score is used for any choice."""
import time, numpy as np
import common as C
F = C.F
s = 100
gn, tn, reg = F.build_noise(s); gn2, tn2, reg2 = C.build_noise_c(s, 'C3')
print('noise identical:', np.array_equal(gn, gn2), np.array_equal(tn, tn2), np.array_equal(reg, reg2), len(gn), len(tn))
tu = C.TUNED_A2
kind = np.array([1, 2, 2]); cf = np.array([1.0, 1.0, 2.8]); kR = np.full(3, tu['kf_R']); kQ = np.full(3, tu['kf_QG'])
t0 = time.time(); J0 = F.est_only(s, C.FX, kind, cfac=cf, kR=kR, kQ=kQ); t1 = time.time()
r = C.run(s, 'C3', [C.BenchEst(kind, cfac=cf, kR=kR, kQ=kQ)]); t2 = time.time()
print('est_only', J0, round(t1 - t0, 1), 's'); print('run     ', r['J'][0], round(t2 - t1, 1), 's')
print('max abs diff', float(np.max(np.abs(J0 - r['J'][0]))))
# timing with many columns
n = 2000; t0 = time.time()
C.run(s, 'C1', [C.BenchEst(np.ones(n, int), kR=np.full(n, 1.0), kQ=np.full(n, 1e-9), qB=np.full(n, 1e-11))])
print(f'{n} KF columns: {time.time() - t0:.1f} s per seed')
