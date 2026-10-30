"""STEP 6g, EXPLORATORY 2 (after the pre-registered FAIL, sha 54f0c7fc, and exploratory 1). Not a judged test.
Diagnosis: all filters start on the first accel sample (5 deg of noise at idle, 50 deg at half throttle); a 30 s window forgets it in about
a minute; the records are 8 to 40 s, scored after 10 s. The start transient dominates every long-window score. Here every filter (law,
oracle, Kalman) is initialised at the true angle, with the bias calibrated and no rate limit as in exploratory 1. This measures the
steady-state quality of the window on real noise, the question 6g asked; it does not measure start-up, which the cradle will."""
import os, sys, json, numpy as np, hashlib
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import step6g_real_noise as G
G.LOG.clear(); G.TAU_RATE = 1e9
_load = G.load
def load_cal(name, f, t0, t1, gyro_lp=False):
    tilt, w, g1k = _load(name, f, t0, t1, gyro_lp); return tilt, w - w.mean(), g1k - g1k.mean()
G.load = load_cal
# initialise at the truth: patch the generator of the first accel sample used as a0 by replacing A[0] with TRU[0] inside run_segment
_run = G.run_segment
def run_init_truth(name, f, t0, t1, lam_hat, gyro_lp=False):
    tilt, w_real, _ = G.load(name, f, t0, t1, gyro_lp); N = len(tilt); TH = G.truth(N, G.SEEDS); tru0 = TH[1] + tilt.mean()
    orig_truth = G.truth
    # make A[0] equal to the truth by shifting tilt[0] to its mean for the first sample only (a0 = TH[1] + tilt[0])
    tilt_mod = tilt.copy(); tilt_mod[0] = tilt.mean()
    G.load = lambda *a, **k: (tilt_mod, w_real, None)
    try: return _run(name, f, t0, t1, lam_hat, gyro_lp)
    finally: G.load = load_cal
if __name__ == '__main__':
    G.P('step6g_exploratory2.py sha256', hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest()[:12], '| init at truth, bias calibrated, no rate limit | seeds', G.SEEDS)
    _, _, g_idle = _load(*G.SEGS[0][:4]); lam_hat, _ = G.allan_lambda(g_idle - g_idle.mean()); res = {}
    G.P(f'{"segment":10s} {"rec":>5s} {"class":>6s} {"law tau":>8s} {"orc tau":>8s} {"regret":>7s} {"RMS law":>7s} {"RMS orc":>7s} {"RMS KF":>7s} {"law/KF":>6s} {"x2":>5s} {"/2":>5s}  KF tuning (R, QG)')
    for name, f, a, b, role in G.SEGS:
        r = run_init_truth(name, f, a, b, lam_hat); res[name] = r
        law = np.array(r['law_rms'][0]); orc = np.array(r['orc_rms']); kf = np.array(r['kf_rms']); L = np.array(r['law_rms'])
        regret = np.median(law / orc - 1); ltau = np.median(r['law_tau_med']); otau = np.median(r['orc_tau']); lim = '*' if otau > 0.25 * r['record_s'] else ' '
        G.P(f'{name:10s} {r["record_s"]:5.1f} {r["class_ratio"]:6.2f} {ltau:8.2f} {otau:8.2f}{lim} {regret*100:6.1f}% {law.mean():7.3f} {orc.mean():7.3f} {kf.mean():7.3f} {law.mean()/kf.mean():6.3f} {L[1].mean()/L[0].mean():5.3f} {L[2].mean()/L[0].mean():5.3f}  {np.median(r["kf_R_mult"]):.3g}, {np.median(r["kf_QG_mult"]):.3g}')
    G.P('  * oracle at the record limit (window > 25 % of record)')
    json.dump(res, open(os.path.join(G.OUT, 'step6g_exploratory2_rows.json'), 'w'), indent=1)
    open(os.path.join(G.OUT, 'step6g_exploratory2.txt'), 'w').write('\n'.join(G.LOG) + '\n')
