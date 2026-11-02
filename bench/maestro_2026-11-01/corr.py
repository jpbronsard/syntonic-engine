import numpy as np, sys
from scipy.stats import spearmanr
from scipy.ndimage import uniform_filter1d
for n in 'AB':
    d = np.load(f'out{n}.npz'); T=d['T']; W=d['win']; qdd=d['qdd']; tp=d['tip_plan']
    acc_tip = np.linalg.norm(np.gradient(np.gradient(tp,0.01,axis=0),0.01,axis=0),axis=1)
    sm = lambda x: uniform_filter1d(x, 100, axis=0)   # 1 s, like the engine's averaging
    m = T>1.0
    print(n, 'tip acc m/s2 median/p95', np.round(np.median(acc_tip),2), np.round(np.percentile(acc_tip,95),2))
    for j in range(6):
        r1 = spearmanr(W[m,j], sm(np.abs(qdd[:,j]))[m]).correlation
        r2 = spearmanr(W[m,j], sm(acc_tip)[m]).correlation
        atfloor = np.mean(W[m,j] < d['floor']*1.02)
        print(f'  RJ{j+1} tau ms [{W[m,j].min()*1e3:.2f},{np.median(W[m,j])*1e3:.2f},{W[m,j].max()*1e3:.2f}]  at floor {atfloor:.0%}  rho(tau, |own acc| 1s) {r1:+.2f}  rho(tau, tip acc 1s) {r2:+.2f}')
