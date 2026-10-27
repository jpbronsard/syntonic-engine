"""EXPLORATORY: inspect the online identification of SYN-1S and SYN-noI at t = 100 s, one seed of A_steady."""
import numpy as np, sys
import step6_quad_attitude_estimator as S
S.SEG_S = 120.0
d = S.generate('A_steady', [600]); W, A1 = d['W'], d['A1']
for label, withI in (('SYN-1S', True), ('SYN-noI', False)):
    m = S.Mahony(A1[0], withI=withI); idr = S.SynID(1, 'one')
    for k in range(100000):
        t = idr.update(A1[k] - m.Th); m.step(W[k], A1[k], t)
    l = (idr.scales * S.DT)
    good = idr.nblk[idr.Aidx] >= 2
    A = np.min(idr.DA[idr.Aidx[good]] * l[idr.Aidx[good], None], axis=0)
    floor = A[None, :] / l[:, None]; D = idr.DA - floor
    print(f'\n{label} at 100 s: tau {idr.tau[0]:.3f}  target {idr.tau_t[0]:.3f}  H {idr.H[0]:.2f}  A {A[0]:.3e}'
          f'  lambda {idr.lam[0]:.3e}  bias_hat {m.b[0]/S.D2R:.4f} deg/s')
    print('   l(s)     DA          floor      D-floor    clears 2x  nblk')
    for i in range(len(l)):
        print(f'  {l[i]:7.3f}  {idr.DA[i,0]:10.3e}  {floor[i,0]:10.3e}  {D[i,0]:10.3e}  {str(D[i,0] > 2*floor[i,0]):6s}  {idr.nblk[i]}')
