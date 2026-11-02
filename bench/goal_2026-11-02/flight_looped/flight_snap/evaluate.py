"""Evaluate every variant on the six circumstances, dev seeds 300-309 (never 400-499), one pass per (circumstance, seed)
with all columns at once. Appends rows to results/<tag>_runs.jsonl (resumable); 'summarize' writes results/<tag>_summary.json
and prints J_est pooled (rms over seeds of per-seed RMS, deg) and the ratio to the rival KF_rival.
Run: python evaluate.py run [tag] [C1 C2 ...]   |   python evaluate.py summarize [tag]"""
import sys, os, json, time, numpy as np
import common as C
import variants as V
import variants2 as V2
import variants3 as V3
PHASES = {'C1': [(0, 20), (20, 40), (40, 60)], 'C2': [(0, 20), (20, 40), (40, 60)], 'C3': [(0, 15), (15, 30), (30, 45), (45, 60)],
          'C4': [(0, 20), (20, 40), (40, 60)], 'C5': [(5 * i, 5 * i + 5) for i in range(12)], 'C6': [(0, 30), (30, 60)]}
CIRCS = ['C1', 'C2', 'C3', 'C4', 'C5', 'C6']


def main():
    mode = sys.argv[1]; tag = sys.argv[2] if len(sys.argv) > 2 else 'dev'
    path = os.path.join(C.OUT, f'{tag}_runs.jsonl')
    if mode == 'run':
        circs = sys.argv[3:] or CIRCS
        done = set()
        if os.path.exists(path):
            for l in open(path): d = json.loads(l); done.add((d['circ'], d['seed']))
        names = None; M = V3 if tag.startswith('r3') else (V2 if tag.startswith('r2') else V)
        if tag == 'ref': M = type('M', (), {'build': staticmethod(V3.build_ref)})
        if tag not in ('dev', 'r2', 'r3', 'ref'):
            names = json.load(open(os.path.join(C.OUT, f'{tag}_names.json')))
        for circ in circs:
            for s in C.DEV:
                if (circ, s) in done: continue
                t0 = time.time(); ests, cols, D, _ = M.build(names)
                r = C.run(s, circ, ests, phases=PHASES[circ])
                J = np.concatenate(r['J']); ph = np.concatenate([np.array(p).T for p in r['phase']], axis=0)   # (cols, phases)
                row = dict(circ=circ, seed=s, cols=cols, J=J.tolist(), phase=ph.tolist(), phases=PHASES[circ])
                with open(path, 'a') as f: f.write(json.dumps(row) + '\n')
                print(f'{circ} seed {s} done in {time.time() - t0:.0f} s', flush=True)
        return
    if mode == 'summarize':
        rows = [json.loads(l) for l in open(path)]
        cols = rows[0]['cols']; M = V3 if tag.startswith('r3') else (V2 if tag.startswith('r2') else V)
        _, _, D, _ = (V3.build_ref() if tag == 'ref' else M.build(cols))
        S = dict(cols=cols, desc=D, circ={})
        for circ in CIRCS:
            rr = sorted([r for r in rows if r['circ'] == circ], key=lambda r: r['seed'])
            if not rr: continue
            assert [r['seed'] for r in rr] == C.DEV, (circ, [r['seed'] for r in rr])
            J = np.array([r['J'] for r in rr]); P = C.pooled(J); iK = cols.index('KF_rival')
            ph = np.sqrt(np.mean(np.square(np.array([r['phase'] for r in rr])), 0))
            wins = (J < J[:, [iK]]).sum(0)
            S['circ'][circ] = dict(pooled={c: float(P[i]) for i, c in enumerate(cols)}, ratio={c: float(P[i] / P[iK]) for i, c in enumerate(cols)},
                                   wins_vs_rival={c: int(wins[i]) for i, c in enumerate(cols)},
                                   phase_pooled={c: ph[i].tolist() for i, c in enumerate(cols)}, phases=rr[0]['phases'])
        C.dump(f'{tag}_summary.json', S)
        circs = [c for c in CIRCS if c in S['circ']]
        print(f'{"variant":16s} ' + ' '.join(f'{c:>14s}' for c in circs) + '   n<=1.00 n<=1.05  max ratio')
        for c in cols:
            rat = [S['circ'][k]['ratio'][c] for k in circs]
            print(f'{c:16s} ' + ' '.join(f'{S["circ"][k]["pooled"][c]:6.3f} ({S["circ"][k]["ratio"][c]:5.2f})' for k in circs)
                  + f'   {sum(x <= 1.0 for x in rat):4d} {sum(x <= 1.05 for x in rat):7d}   {max(rat):6.2f}')


if __name__ == '__main__':
    main()
