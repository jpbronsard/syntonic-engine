"""Checks before any dev-seed run: (1) Law1 with the published rule == F.Est kind 2 (seed 100, C3, published case);
(2) Law2SS with the bias loop forced off and the published rule ~ the law; (3) the windows each rule gives at the records'
own s2 (no truth); (4) timing of one flight with all variants."""
import time, math, numpy as np
import common as C
import variants as V
F = C.F
# (1)
pub = lambda s2, reg: math.sqrt(max(s2, 1e-300) / V.LAM[reg])
r = C.run(100, 'C3', [C.BenchEst(np.array([2]), cfac=np.array([1.0])), V.Law1([pub])])
print('(1) Est law', r['J'][0], ' Law1 pub', r['J'][1], ' diff', float(abs(r['J'][0][0] - r['J'][1][0])))
# (3) windows at the records' s2
ests, cols, D, tabs = V.build()
print('(3) windows (s) at s2 = the record\'s own reading; two-state: 1/k1 and 1/sqrt(k2)')
for reg in V.REGS:
    s2 = V.F2['curves']['raw'][reg]['s2_rec']
    row = f'  {reg:5s} s2 {s2:.2e}: pub {pub(s2, reg):6.1f}'
    for nm in ('a_fp_raw', 'a_min_raw', 'a_fp_loop', 'a_min_loop'):
        row += f' | {nm} {float(V.lookup(tabs[nm], reg, s2)):6.1f}'
    k1, k2 = V.gains_closed(s2, V.QT, V.QB4); row += f' | b_ss {1 / k1:5.1f}/{1 / math.sqrt(k2):5.1f}'
    k1, k2 = V.gains_closed(s2, V.QT, V.QB2); row += f' | b_qb2 {1 / k1:5.1f}/{1 / math.sqrt(k2):5.1f}'
    for nm in ('c_raw', 'c_loop'):
        K1, K2 = tabs[nm][reg]; k1 = float(np.exp(np.interp(np.log(s2), np.log(V.S2G), np.log(K1)))); k2 = float(np.exp(np.interp(np.log(s2), np.log(V.S2G), np.log(K2))))
        row += f' | {nm} {1 / k1:5.1f}/{1 / math.sqrt(k2):5.1f}'
    print(row)
# (4) timing
t0 = time.time(); r = C.run(301, 'C1', ests); print(f'(4) one flight, {len(cols)} columns in {len(ests)} objects: {time.time() - t0:.1f} s')
print('columns:', cols)
