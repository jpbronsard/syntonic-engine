"""For the colleague, after pre-registration: run named variants (any round) together with the rival KF_rival on
given seeds and circumstances, one pass per (circumstance, seed); prints pooled J_est (deg) and the ratio to the rival.
Seeds 400-499 need GOAL_ALLOW_CONFIRM=1 (a guard against accidental use during development).
  python3 confirm.py --names b4_nm_ic2,b2_nm_ic2 --seeds 400-409 [--circs C1,C2,C3,C4,C5,C6] [--out results/confirm_runs.jsonl]"""
import sys, os, json, time, argparse, numpy as np
import common as C
import variants as V, variants2 as V2, variants3 as V3
from evaluate import PHASES, CIRCS


def registry():
    r1 = V.build()[1]; r2 = [c for c in V2.build()[1] if c not in r1]; r3 = [c for c in V3.build()[1] if c not in r1]
    return r1, r2, r3


def build(names):
    r1, r2, r3 = registry(); names = list(dict.fromkeys(['KF_rival'] + names))
    ests, cols, D = [], [], {}
    for M, reg in ((V, r1), (V2, r2), (V3, r3)):
        sub = [n for n in names if n in reg]
        if sub:
            e, c, d, _ = M.build(sub); ests += e; cols += c; D.update(d)
    missing = [n for n in names if n not in cols]
    assert not missing, f'unknown variant names: {missing}'
    return ests, cols, D


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--names', required=True); ap.add_argument('--seeds', required=True)
    ap.add_argument('--circs', default=','.join(CIRCS)); ap.add_argument('--out', default=os.path.join(C.OUT, 'confirm_runs.jsonl'))
    a = ap.parse_args()
    lo, _, hi = a.seeds.partition('-'); seeds = list(range(int(lo), int(hi) + 1)) if hi else [int(x) for x in a.seeds.split(',')]
    names = a.names.split(','); res = {}
    for circ in a.circs.split(','):
        Js = []
        for s in seeds:
            t0 = time.time(); ests, cols, D = build(names); r = C.run(s, circ, ests, phases=PHASES[circ])
            J = np.concatenate(r['J']); Js.append(J)
            with open(a.out, 'a') as f: f.write(json.dumps(dict(circ=circ, seed=s, cols=cols, J=J.tolist())) + '\n')
            print(f'{circ} seed {s} {time.time() - t0:.0f} s', flush=True)
        P = C.pooled(np.array(Js)); iK = cols.index('KF_rival'); res[circ] = {c: (float(P[i]), float(P[i] / P[iK])) for i, c in enumerate(cols)}
    for c in cols:
        print(f'{c:18s} ' + ' '.join(f'{k} {res[k][c][0]:.3f} ({res[k][c][1]:.2f})' for k in res))
