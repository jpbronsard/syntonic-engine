"""STEP 6g, EXPLORATORY (after the pre-registered run, sha 54f0c7fc, failed H6g.1 to H6g.3). Not a judged test.
Diagnosis to test: the failure is the real gyro's uncalibrated DC bias meeting the law's slow bias integral, plus the engine's
start-up ramp (0.5 s start, 0.4 s/s rate limit) on records of 8 to 40 s, rather than the vibration itself.
Two changes, declared: (1) gyro bias calibrated at rest = the mean rate over the record is subtracted (optimistic for held frames:
their own slow motion is removed too; this is what an arming calibration gives a real autopilot); (2) no rate limit and the window
starts at the law's first reading (TAU0 unchanged otherwise). Same seeds, same scoring, same comparators as the pre-registered run."""
import os, sys, json, time, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import step6g_real_noise as G

G.LOG.clear()
G.TAU_RATE = 1e9
_load = G.load
def load_cal(name, f, t0, t1, gyro_lp=False):
    tilt, w, g1k = _load(name, f, t0, t1, gyro_lp)
    return tilt, w - w.mean(), g1k - g1k.mean()
G.load = load_cal

if __name__ == '__main__':
    import hashlib
    G.P('step6g_exploratory.py sha256', hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest()[:12], '| bias calibrated, no rate limit | seeds', G.SEEDS)
    _, _, g_idle = _load(*G.SEGS[0][:4]); lam_hat, slope = G.allan_lambda(g_idle - g_idle.mean())
    G.P(f'gyro fact unchanged: lambda_hat = {lam_hat:.3e} rad^2/s')
    res = {}
    G.P(f'{"segment":10s} {"rec":>5s} {"bias":>6s} {"class":>6s} {"law tau":>8s} {"orc tau":>8s} {"regret":>7s} {"RMS law":>7s} {"RMS orc":>7s} {"RMS KF":>7s} {"law/KF":>6s} {"x2":>5s} {"/2":>5s}')
    for name, f, a, b, role in G.SEGS:
        tilt, w_raw, _ = _load(name, f, a, b); bias = w_raw.mean() / G.D2R
        r = G.run_segment(name, f, a, b, lam_hat); res[name] = r
        law = np.array(r['law_rms'][0]); orc = np.array(r['orc_rms']); kf = np.array(r['kf_rms']); L = np.array(r['law_rms'])
        regret = np.median(law / orc - 1); ltau = np.median(r['law_tau_med']); otau = np.median(r['orc_tau']); lim = '*' if otau > 0.25 * r['record_s'] else ' '
        G.P(f'{name:10s} {r["record_s"]:5.1f} {bias:+6.2f} {r["class_ratio"]:6.2f} {ltau:8.2f} {otau:8.2f}{lim} {regret*100:6.1f}% {law.mean():7.3f} {orc.mean():7.3f} {kf.mean():7.3f} {law.mean()/kf.mean():6.3f} {L[1].mean()/L[0].mean():5.3f} {L[2].mean()/L[0].mean():5.3f}')
    G.P('  bias = mean real gyro rate over the record, deg/s (removed in this run); * oracle at the record limit')
    json.dump(res, open(os.path.join(G.OUT, 'step6g_exploratory_rows.json'), 'w'), indent=1)
    open(os.path.join(G.OUT, 'step6g_exploratory.txt'), 'w').write('\n'.join(G.LOG) + '\n')
