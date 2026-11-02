"""Combined table from results/<tag>_runs.jsonl files: pooled J_est (deg) per variant and circumstance, ratio to the
rival (KF_rival, same seeds), counts of circumstances with ratio <= 1.00 and <= 1.05, and a paired bootstrap over the
10 dev seeds (2000 resamples) for the pooled ratio: [2.5 %, 97.5 %].
Run: python report_table.py dev r2 [--ci name1,name2,...]"""
import sys, json, os, numpy as np
import common as C
CIRCS = ['C1', 'C2', 'C3', 'C4', 'C5', 'C6']
tags = [a for a in sys.argv[1:] if not a.startswith('--')]
ci_names = []
for a in sys.argv[1:]:
    if a.startswith('--ci='): ci_names = a[5:].split(',')
data = {}   # name -> circ -> per-seed J (10,)
for tag in tags:
    rows = [json.loads(l) for l in open(os.path.join(C.OUT, f'{tag}_runs.jsonl'))]
    for r in rows:
        for c, j in zip(r['cols'], r['J']):
            data.setdefault(c, {}).setdefault(r['circ'], {})[r['seed']] = j
rng = np.random.default_rng(12345); B = 2000; idx = rng.integers(0, 10, (B, 10))
out = {}
print(f'{"variant":18s} ' + ' '.join(f'{c:>14s}' for c in CIRCS) + '  <=1.00 <=1.05  max')
for name in data:
    if not all(c in data[name] and len(data[name][c]) == 10 for c in CIRCS): continue
    row = {}; rat = []
    for c in CIRCS:
        J = np.array([data[name][c][s] for s in C.DEV]); K = np.array([data['KF_rival'][c][s] for s in C.DEV])
        p, pk = np.sqrt(np.mean(J ** 2)), np.sqrt(np.mean(K ** 2)); r = p / pk; rat.append(r)
        bs = np.sqrt(np.mean(J[idx] ** 2, 1)) / np.sqrt(np.mean(K[idx] ** 2, 1))
        row[c] = dict(J=float(p), ratio=float(r), ci=[float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))],
                      wins=int(np.sum(J < K)))
    out[name] = dict(circ=row, n100=int(sum(x <= 1.0 for x in rat)), n105=int(sum(x <= 1.05 for x in rat)), max_ratio=float(max(rat)))
    print(f'{name:18s} ' + ' '.join(f'{row[c]["J"]:6.3f} ({row[c]["ratio"]:5.2f})' for c in CIRCS) + f'  {out[name]["n100"]:5d} {out[name]["n105"]:6d}  {max(rat):5.2f}')
for name in ci_names:
    if name in out:
        print(f'  {name:18s} CI95 of ratio: ' + ' '.join(f'{c} [{out[name]["circ"][c]["ci"][0]:.2f},{out[name]["circ"][c]["ci"][1]:.2f}] w{out[name]["circ"][c]["wins"]}' for c in CIRCS))
C.dump('table_' + '_'.join(tags) + '.json', out)
