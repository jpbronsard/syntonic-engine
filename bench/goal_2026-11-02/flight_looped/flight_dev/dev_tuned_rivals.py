"""dev seeds 300-309: the commissioned scheduled rivals (results/sched_tuned.json), to set the predictions"""
import sys, os, json
sys.path.insert(0, '/home/claude/work/goal/confirm/flight_snap')
import numpy as np, common as C, law_flight as L
ST = json.load(open(os.path.join(C.OUT, 'sched_tuned.json')))
out = '/home/claude/work/goal/confirm/flight_dev/tuned_rivals_dev.jsonl'
for circ in sys.argv[1].split(','):
    for s in range(300, 310):
        r = C.run(s, circ, [L.kf_sched_from(ST['KF_sched'], False), L.kf_sched_from(ST['KF_sched_btab'], True)])
        with open(out, 'a') as f: f.write(json.dumps(dict(circ=circ, seed=s, cols=['KF_sched_t', 'KF_sched_btab_t'], J=np.concatenate(r['J']).tolist())) + '\n')
