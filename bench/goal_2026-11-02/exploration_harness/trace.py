"""Diagnostic trace on one dev seed: estimation error (deg) of the rival and of selected variants, and their internal
windows / gains / bias estimate, around a time window. Run: python trace.py <circ> <seed> <t0> <t1> [names...]"""
import sys, numpy as np
import common as C
import variants as V
circ, seed, t0, t1 = sys.argv[1], int(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])
names = sys.argv[5:] or ['KF_rival', 'b_ss', 'b_ss_qb2', 'LAW_c1']
ests, cols, D, _ = V.build(names)
gn, tn, reg = C.build_noise_c(seed, circ); N = len(gn)
lam_of = {r: C.FX['lam'][r] for r in ('idle', 'low', 'half')}
gs = 0.0; bc = 0.0; rows = []
for k in range(N):
    g0 = gn[k] - bc; a0 = tn[k // 4] if k % 4 == 0 else None
    outs = [e.step(g0, a0, lam_of[reg[k]], k * C.DT, reg[k], k) for e in ests]
    if k < C.KPRE:
        gs += gn[k]
        if k == C.KPRE - 1:
            bc = gs / C.KPRE
            for e in ests:
                if hasattr(e, 'arm'): e.arm(bc)
        continue
    tf = (k - C.KPRE) * C.DT
    if t0 <= tf < t1 and (k - C.KPRE) % 100 == 0:
        err = np.concatenate([o for o in outs]) / C.D2R
        inner = []
        for e in ests:
            if isinstance(e, V.Law2SS): inner += [f'1/k1={1 / x:5.1f} b={y / C.D2R:+.4f}' for x, y in zip(e.k1, e.b)]
            elif isinstance(e, V.Law1): inner += [f'tau={x:5.1f}' for x in e.tau]
            elif isinstance(e, V.Law2TV): inner += [f'k1/DTA={x:.2e} b={y / C.D2R:+.4f}' for x, y in zip(e.p11 / (e.p11 + (e.R if e.R is not None else 1)), e.b)]
            elif isinstance(e, V.KF): inner += [f'b={y / C.D2R:+.4f}' for y in e.e.b]
        rows.append((tf, reg[k], err, inner))
print(circ, seed, 'columns', cols)
for tf, rg, err, inner in rows:
    print(f'{tf:6.1f} {rg:4s} ' + ' '.join(f'{x:+7.3f}' for x in err) + '   ' + ' | '.join(inner))
