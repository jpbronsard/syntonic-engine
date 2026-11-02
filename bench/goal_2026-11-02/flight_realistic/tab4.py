import json, sys, numpy as np
fn = sys.argv[1] if len(sys.argv) > 1 else 'dev_r4_runs.jsonl'
rows = [json.loads(l) for l in open(fn)]; cols = rows[0]['cols']; circs = ['C1','C2','C3','C4','C5','C6']
P, Wn = {}, {}
for c in circs:
    rr = sorted([r for r in rows if r['circ'] == c], key=lambda r: r['seed']); assert len(rr) == 10, (c, len(rr))
    J = np.array([r['J'] for r in rr]); P[c] = np.sqrt(np.mean(J ** 2, 0)); Wn[c] = J
print('pooled deg')
for i, n in enumerate(cols): print('  %-15s' % n, ' '.join('%.4f' % P[c][i] for c in circs))
iR, iF = cols.index('KF_rival3'), cols.index('KF_full3')
print('ratio to the best tuned Kalman per circumstance; max; wins vs it')
for i, n in enumerate(cols):
    ib = {c: iR if P[c][iR] < P[c][iF] else iF for c in circs}
    r = [P[c][i] / P[c][ib[c]] for c in circs]; w = [int((Wn[c][:, i] < Wn[c][:, ib[c]]).sum()) for c in circs]
    print('  %-15s' % n, ' '.join('%6.3f' % x for x in r), '  max %.3f' % max(r), ' wins', w)
