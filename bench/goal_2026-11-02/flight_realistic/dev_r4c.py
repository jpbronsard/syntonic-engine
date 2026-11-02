"""Dev seeds 300-309 on the amendment-3 bench: rivals (as tuned) + LAW_I + Riccati-law columns. Appends to dev_bb_runs.jsonl.
Run: python3 dev_bb.py C1,C2,C3 [rivals: rival,full,sched]"""
import os, sys, json, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import rivals_bb as RB
import bench_bb as B
import ricv2 as R
C, V = B.C, B.V
import variants3 as V3
QB4, QB2 = V.QB4, V.QB2
b = lambda **k: dict(dict(k=3, qt='level', prior='vb', qb=QB4), **k)
VAR = [('L3_vb_sa_max', b(sa=True, rm='max')), ('L3_vb_sa_max_qb2', b(sa=True, rm='max', qb=QB2)),
       ('L3_vb_sa_max_idle_qb2', b(sa=True, rm='max', qt='idle', qb=QB2)), ('L1_sa_max_qb2', b(k=1, prior='cal', sa=True, rm='max', qb=QB2))]


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
    out = os.path.join(HERE, 'dev_r4c_runs.jsonl' if laws else 'dev_r4_rivals.jsonl')
    for circ in circs:
        for s in range(300, 310):
            t0 = time.time(); ests, cols = build(kinds, laws); r = B.run4(s, circ, ests)
            J = np.concatenate([np.atleast_1d(j) for j in r['J']]).tolist()
            with open(out, 'a') as f: f.write(json.dumps(dict(circ=circ, seed=s, cols=cols, J=J)) + '\n')
            print(circ, s, f'{time.time() - t0:.0f} s', flush=True)
