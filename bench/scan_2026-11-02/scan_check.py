# Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
# Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
"""The two speed-ups of 0.3.2 against the recorded choices of the Python law (a check, not a prediction).

Reads the loop-law replay logs of rust/tools/make_law_replay_log.py (spread 0). At every gain recomputation of the
log it runs the scan of syntonic_engine.loop_scan (the arithmetic of the C++ and Rust laws) on the statistics the
Python law used, warm-started from the gains in use, in the four combinations of the speed-ups, and checks that the
candidate chosen is the one the Python law chose (the log's gain index). It reports:
  - the share of the (candidate, frequency) points evaluated with pruning (mean and largest over the recomputations),
  - whether every quadratic-form choice was certified (no fallback),
  - the smallest relative gap between the best and the next J (scan without pruning, so the true runner-up), with
    the quadratic form and from L,
  - the largest relative difference between |S|^2 from the quadratic form and from L over all robust candidates.
Run from the folder that holds the logs:
  python bench/scan_2026-11-02/scan_check.py replay_loop_two.bin replay_loop_one.bin replay_loop_bench.bin
"""
import os, sys, json
import numpy as np
REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
sys.path.insert(0, os.path.join(REPO, 'python'))
from syntonic_engine.loop_law import LoopLaw, loop_delay  # noqa: E402
from syntonic_engine.loop_scan import LoopScan  # noqa: E402

MODES = (('plain', False, False), ('prune', True, False), ('quad', False, True), ('prune+quad', True, True))


def check(path):
    v = np.fromfile(path, '<f8'); h = v[:16]; N, T = int(v[1]), int(v[2])
    if int(h[8]) != 0:
        return None
    rows = v[16 + 2 * N:].reshape(T, 15 * N + 1); out = rows[:, 4 * N + 1:]
    dt, every = h[3], int(h[9]); D = loop_delay(dt, h[4]) + h[6]
    LL = LoopLaw(h[5], D); sc = LoopScan(h[5], D, LL)
    gi = out[:, 10 * N:11 * N].astype(int)
    s2d = 0.0
    for i, kis in sc.rows:
        a, b = sc.s2(i, kis, True), sc.s2(i, kis, False); s2d = max(s2d, float(np.max(np.abs(a / b - 1))))
    res = {m: dict(bad=0, share=[], fallback=0, margin=np.inf) for m, _, _ in MODES}
    nrec = 0
    for k in range(T):
        if k % every or gi[k, 0] < 0:
            continue
        nrec += 1
        A, gam, Oc, r = out[k, 5 * N:6 * N], out[k, 6 * N:7 * N], out[k, 7 * N:8 * N], out[k, 4 * N:5 * N] * dt
        warm = [int(c) for c in gi[k - 1] if c >= 0] if k > 0 else []
        for m, prune, quad in MODES:
            c, st = sc.choose(A, gam, Oc, r, warm, prune, quad)
            R = res[m]
            if c is None:
                R['fallback'] += 1
            elif np.any(c != gi[k]):
                R['bad'] += 1
            R['share'].append(st['points'] / st['total'])
            if not prune:
                R['margin'] = min(R['margin'], st['margin'])
    rep = dict(log=os.path.basename(path), joints=N, robust=len(sc.okc), recomputations=nrec, s2_rel_diff=s2d)
    for m, prune, quad in MODES:
        R = res[m]; s = np.array(R['share'])
        rep[m] = dict(differ=R['bad'], fallbacks=R['fallback'], share_mean=float(s.mean()), share_max=float(s.max()))
        if not prune:
            rep[m]['min_margin'] = R['margin']
    return rep


if __name__ == '__main__':
    reps = []
    for p in sys.argv[1:]:
        r = check(p)
        if r is None:
            print(f'{p}: a spread log, skipped'); continue
        reps.append(r)
        print(f"{r['log']}: {r['joints']} joint(s), {r['robust']} robust candidates, {r['recomputations']} recomputations; "
              f"|S|^2 quadratic form vs from L: largest relative difference {r['s2_rel_diff']:.2e}")
        for m, _, _ in MODES:
            x = r[m]
            print(f"    {m:11s} choice differs from Python's at {x['differ']} recomputation(s), fallbacks {x['fallbacks']}, "
                  f"points evaluated {x['share_mean']:.1%} (largest {x['share_max']:.1%})"
                  + (f", smallest gap best-next J {x['min_margin']:.2e}" if 'min_margin' in x else ''))
    out = os.environ.get('SCAN_CHECK_JSON')
    if out:
        json.dump(reps, open(out, 'w'), indent=1)
