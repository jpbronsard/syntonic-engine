"""diagnostics only: the expected value of the replayed gyro and tilt under the block bootstrap (edge under-weighting)"""
import sys, numpy as np
sys.path.insert(0, '/home/claude/work/goal/r3')
import bench_bb as B
B.install(); C, V = B.C, B.V
R = C.F.records(); D2R = C.D2R
gi_mean = R['idle'][1].mean()
for r in V.REGS:
    tilt, g = R[r]; n = len(tilt); L = int(round(B.TLIM[r] / C.DTA))
    # coverage weight of tilt sample i: number of offsets o in [0, n-L] with o <= i < o+L
    i = np.arange(n); w = np.minimum(i, n - L) - np.maximum(0, i - L + 1) + 1; w = w / w.sum()
    wg = np.repeat(w, C.NSUB) / C.NSUB
    gm_bb = (wg * g).sum(); tm_bb = (w * tilt).sum()
    print(f"{r}: n {n} L {L}: gyro record mean {g.mean():.4e}, BB-weighted mean {gm_bb:.4e} (diff {gm_bb-g.mean():.3e} rad/s = {(gm_bb-g.mean())/D2R:.4f} deg/s); "
          f"tilt BB-weighted mean {tm_bb:.4e} rad = {tm_bb/D2R:.4f} deg; BB gyro mean - idle record mean {gm_bb-gi_mean:.3e}")
# idle: the arming calibration is the mean over the 5 s pre-roll drawn by BB blocks: expected = idle BB-weighted mean
