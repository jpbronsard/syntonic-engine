import json, sys, numpy as np
fn, refs = sys.argv[1], sys.argv[2].split(',')
rows = [json.loads(l) for l in open(fn)]
ric = {(r['circ'], r['seed']): r for r in (json.loads(l) for l in open('ric_runs.jsonl'))}
cols = rows[0]['cols']; circs = ['C1','C2','C3','C4','C5','C6']; P = {}
for c in circs:
    rr = sorted([r for r in rows if r['circ'] == c], key=lambda r: r['seed'])
    if len(rr) < 10: print(c, 'incomplete', len(rr)); continue
    J = np.array([r['J'] for r in rr]); P[c] = dict(zip(cols, np.sqrt(np.mean(J ** 2, 0))))
    k = [ric[(c, r['seed'])] for r in rr]; JK = np.array([x['J'][x['cols'].index('KF_rival')] for x in k]); P[c]['KF_rival'] = float(np.sqrt(np.mean(JK ** 2)))
    P[c]['W'] = {n: int((J[:, i] < J[:, cols.index(refs[0])]).sum()) for i, n in enumerate(cols)} if refs[0] in cols else {}
cs = [c for c in circs if c in P]
for ref in refs:
    print('ratio to', ref, ' '.join('%7s' % c for c in cs))
    for n in cols:
        r = [P[c][n] / P[c][ref] for c in cs]; print('  %-20s' % n, ' '.join('%7.3f' % x for x in r), '  max %.3f' % max(r))
print('pooled', {c: {n: round(float(v), 4) for n, v in P[c].items() if n != 'W'} for c in cs})
