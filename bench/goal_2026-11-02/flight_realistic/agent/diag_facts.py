"""diagnostics only (understanding the records under block bootstrap); no constant from here enters the law"""
import sys, numpy as np
sys.path.insert(0, '/home/claude/work/goal/r3')
import bench_bb as B
B.install(); C, V = B.C, B.V
F2 = V.F2; TG = np.array(F2['TG'])
for r in V.REGS:
    cv = F2['curves']['raw'][r]
    TL = cv['T_lim']; i = int(np.argmin(np.abs(TG - TL)))
    print(f"{r}: T_rec {cv['T_rec']:.2f} T_lim {TL:.2f} s2_rec {cv['s2_rec']:.3e} r_top {cv['r_top']:.3e} q_top {cv['q_top']:.3e} "
          f"b_res {cv['b_res']:.3e} LAM {V.LAM[r]:.3e}")
    for j in range(0, len(TG), 4):
        print(f"   T {TG[j]:8.3f}  r {cv['r'][j]:.3e}  lam {cv['lam'][j]:.3e}  lamA {cv['lamA'][j]:.3e}  sz {cv['sz'][j]:.3e}  meas {cv['measured'][j]}")
print('QT', V.QT, 'QB4', V.QB4, 'QB2', V.QB2, 'F.Q_B', C.F.Q_B, 'var_tilt', C.FX['var_tilt'], 'lam_hat', C.FX['lam_hat'])
# empirical long-run ARW and tilt density under the block bootstrap (understanding only)
R = C.F.records()
for r in V.REGS:
    tilt, g = R[r]; L = int(round(B.TLIM[r] / C.DTA)); Lg = L * C.NSUB
    bm = g.mean(); offs = np.arange(0, len(tilt) - L + 1)
    cs = np.concatenate([[0.0], np.cumsum(g - bm)]) * C.DT
    I = cs[offs * C.NSUB + Lg] - cs[offs * C.NSUB]
    ct = np.concatenate([[0.0], np.cumsum(tilt)])
    Mt = (ct[offs + L] - ct[offs]) / L
    print(f"{r}: BB long-run gyro ARW var(block int)/L = {np.mean(I**2)/B.TLIM[r]:.3e} (mean of block int {I.mean():.2e}); "
          f"tilt density var(block mean)*L = {np.mean(Mt**2)*B.TLIM[r]:.3e}; record bias {bm:.4e} rad/s; tilt sample var {np.mean(tilt**2):.3e}")
