"""Dev seeds 300-309 on the amendment-3 bench: rivals (as tuned) + LAW_I + Riccati-law columns. Appends to dev_bb_runs.jsonl.
Run: python3 dev_bb.py C1,C2,C3 [rivals: rival,full,sched]"""
import os, sys, json, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import rivals_bb as RB
import bench_bb as B
import ricv as R
C, V = B.C, B.V
import variants3 as V3
QB4, QB2 = V.QB4, V.QB2
VAR = [('L3_lvl_mag', dict(k=3, qt='level', prior='mag', qb=QB4)), ('L3_lvl_vb', dict(k=3, qt='level', prior='vb', qb=QB4)),
       ('L3_lvl_diff', dict(k=3, qt='level', prior='diff', qb=QB4)), ('L3_idle_mag', dict(k=3, qt='idle', prior='mag', qb=QB4)),
       ('L3_lvl_mag_qb2', dict(k=3, qt='level', prior='mag', qb=QB2)), ('L1_lvl', dict(k=1, qt='level', prior='cal', qb=QB4)),
       ('L1_lvl_sw', dict(k=1, qt='level', prior='cal', qb=QB4, sw=True)), ('L1_idle', dict(k=1, qt='idle', prior='cal', qb=QB4))]


def rival_params():
    d = {}
    for f in os.listdir(os.path.join(HERE, 'results')):
        if f.startswith('rivals_r4_') and f.endswith('.json'):
            for k, v in json.load(open(os.path.join(HERE, 'results', f))).items(): d[k] = v['params']
    return d


def build(kinds, laws=True):
    rp = rival_params(); ests, cols = [], []
    for kind in kinds:
        ests.append(RB.make(kind, np.array([rp[kind]]))); cols.append('KF_' + kind + '3')
    if laws:
        e3, c3, _, _ = V3.build(['b4_nm_ic2']); ests += e3; cols.append('LAW_I')
        ests.append(R.RicV([v for _, v in VAR])); cols += [n for n, _ in VAR]
    return ests, cols


if __name__ == '__main__':
    circs = sys.argv[1].split(','); kinds = sys.argv[2].split(',') if len(sys.argv) > 2 else ['rival', 'full']
    laws = len(sys.argv) <= 3 or sys.argv[3] != 'nolaws'
    out = os.path.join(HERE, 'dev_r4_runs.jsonl' if laws else 'dev_r4_rivals.jsonl')
    for circ in circs:
        for s in range(300, 310):
            t0 = time.time(); ests, cols = build(kinds, laws); r = B.run4(s, circ, ests)
            J = np.concatenate([np.atleast_1d(j) for j in r['J']]).tolist()
            with open(out, 'a') as f: f.write(json.dumps(dict(circ=circ, seed=s, cols=cols, J=J)) + '\n')
            print(circ, s, f'{time.time() - t0:.0f} s', flush=True)
