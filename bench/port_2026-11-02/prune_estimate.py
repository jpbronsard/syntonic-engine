# Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
# Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
"""How much an exact pruning would save in the loop law's gain recomputation (0.3.1; a proposal, not in the law).

J of a candidate is a sum of positive terms plus r x MTr. Summed in a fixed order, its partial sum never decreases, so
once the partial sum plus r x MTr exceeds the smallest J found so far, the candidate cannot be chosen and the rest of
its sum can be skipped without changing the chosen gains. With the gains in use evaluated first as the bound, this
counts, at every recomputation of a replay log of rust/tools/make_law_replay_log.py, the share of the (candidate,
frequency) points that would still be evaluated (frequencies in ascending order; a candidate stops when every joint has
passed its bound).
Run from the folder that holds the logs:   python bench/port_2026-11-02/prune_estimate.py replay_loop_bench.bin replay_loop_two.bin
"""
import os, sys
import numpy as np
REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
sys.path.insert(0, os.path.join(REPO, 'python'))
from syntonic_engine.loop_law import LoopLaw, loop_delay, K_gamma  # noqa: E402


def share(path):
    v = np.fromfile(path, '<f8'); h = v[:16]; N, T = int(v[1]), int(v[2])
    rows = v[16 + 2 * N:].reshape(T, 15 * N + 1); out = rows[:, 4 * N + 1:]
    dt, every = h[3], int(h[9]); LL = LoopLaw(h[5], loop_delay(dt, h[4]) + h[6]); ok = list(np.flatnonzero(LL.ok))
    gi = out[:, 10 * N:11 * N]; nw = len(LL.W); total = done = 0
    for k in range(1, T):
        if k % every or gi[k - 1, 0] < 0:          # a recomputation with gains already in use (the bound)
            continue
        A, gam, Oc, r = out[k, 5 * N:6 * N], out[k, 6 * N:7 * N], out[k, 7 * N:8 * N], out[k, 4 * N:5 * N] * dt
        Phi = (A / K_gamma(gam))[:, None] * (LL.W ** 2 + Oc[:, None] ** 2) ** (-(1 + gam[:, None]) / 2)
        cum = np.cumsum(Phi[:, None, :] * LL.MS[None, ok, :], axis=2) + (r[:, None] * LL.MTr[None, ok])[:, :, None]
        B = np.array([cum[j, ok.index(int(gi[k - 1, j])), -1] for j in range(N)]); done += nw
        for p in range(len(ok)):
            need = 0
            for j in range(N):
                ex = np.flatnonzero(cum[j, p] > B[j]); need = max(need, ex[0] + 1 if len(ex) else nw)
            done += need; total += nw
            B = np.minimum(B, cum[:, p, -1])
    return N, done / total


if __name__ == '__main__':
    for path in sys.argv[1:]:
        N, s = share(path)
        print(f'{path}: {N} joint(s); exact pruning would evaluate {s:.1%} of the (candidate, frequency) points')
