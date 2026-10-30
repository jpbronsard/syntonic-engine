"""STEP 6g, EXPLORATORY 3. The law's WINDOW scored as a fixed window (the filter's own start-up cannot be scored on 8 to 40 s records).
For each segment: s2 read by the engine's estimator over the whole record (min over 32/64/128 ms of D_H l / 6 on z_raw), lambda_hat
from S6-idle; law window tau_law = sqrt(s2 / lambda_hat), uncapped and capped at 30 s. Fixed-window Mahony initialised at the true
angle, bias calibrated; RMS after the 10 s / 30 % skip; 10 seeds; regret against the oracle on the same grid, refined. Not judged."""
import os, sys, json, numpy as np, hashlib
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import step6g_real_noise as G
G.LOG.clear()
def s2_engine(z):
    out = []
    for s in G.HSCALES:
        mb = z[: len(z) // s * s].reshape(-1, s).mean(1); dd = mb[2:] - 2 * mb[1:-1] + mb[:-2]; out.append(np.mean(dd ** 2) * (s * G.DT) / 6.0)
    return min(out), out
def fixed_rms(tilt, w, taus, seeds):
    N = len(tilt); ns = len(seeds); TH = G.truth(N, seeds); W = np.diff(TH, axis=0) / G.DT + w[:, None]; A = TH[1:] + tilt[:, None]; TRU = TH[1:] + tilt.mean()
    skip = int(round(min(10.0, 0.3 * N * G.DT) / G.DT)); nt = len(taus)
    m = G.Mahony(np.tile(TRU[0], nt)); tau = np.repeat(taus, ns); ss = 0.0; n = 0
    for k in range(N):
        e = m.step(np.tile(W[k], nt), np.tile(A[k], nt), tau) - np.tile(TRU[k], nt)
        if k >= skip: ss = ss + e ** 2; n += 1
    return (np.sqrt(ss / n) / G.D2R).reshape(nt, ns)
if __name__ == '__main__':
    G.P('step6g_exploratory3.py sha256', hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest()[:12], '| fixed windows from the true angle')
    _, _, g_idle = G.load(*G.SEGS[0][:4]); lam_hat, _ = G.allan_lambda(g_idle - g_idle.mean())
    G.P(f'{"segment":10s} {"s2 (rad^2 s)":>12s} {"tau law":>8s} {"capped":>7s} {"oracle":>7s} {"RMS law":>8s} {"RMS cap":>8s} {"RMS orc":>8s} {"regret law":>10s} {"regret cap":>10s}   RMS at 1, 3, 10 s')
    res = {}
    for name, f, a, b, role in G.SEGS[:3]:
        tilt, w, _ = G.load(name, f, a, b); w = w - w.mean()
        Th = np.cumsum(w) * G.DT; s2, by = s2_engine(tilt - Th)       # z_raw with the synthetic truth cancelling exactly
        tl = np.sqrt(s2 / lam_hat); tc = min(tl, 30.0)
        grid = np.unique(np.concatenate([G.TAU_GRID, [tl, tc, 1.0, 3.0, 10.0]]))
        R = fixed_rms(tilt, w, grid, G.SEEDS); med = np.median(R, 1); lg = np.log(grid)
        o_tau, o_rms = G.refine(lg, med); i_l = np.argmin(np.abs(grid - tl)); i_c = np.argmin(np.abs(grid - tc))
        r_l = med[i_l] / o_rms - 1; r_c = med[i_c] / o_rms - 1
        at = lambda v: med[np.argmin(np.abs(grid - v))]
        G.P(f'{name:10s} {s2:12.3e} {tl:8.1f} {tc:7.1f} {o_tau:7.1f} {med[i_l]:8.3f} {med[i_c]:8.3f} {o_rms:8.3f} {r_l*100:9.1f}% {r_c*100:9.1f}%   {at(1):.3f}, {at(3):.3f}, {at(10):.3f}')
        res[name] = dict(s2=s2, s2_by_scale=by, tau_law=tl, tau_cap=tc, oracle_tau=o_tau, oracle_rms=o_rms, grid=grid.tolist(), rms_med=med.tolist())
    G.P('  regret = RMS(window) / RMS(oracle) - 1, seed medians; record lengths 35, 40, 8 s so windows beyond ~10 s are only partly resolved')
    json.dump(res, open(os.path.join(G.OUT, 'step6g_exploratory3_rows.json'), 'w'), indent=1)
    open(os.path.join(G.OUT, 'step6g_exploratory3.txt'), 'w').write('\n'.join(G.LOG) + '\n')
