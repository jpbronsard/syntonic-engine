"""EXPLORATORY: second differences of block means (Hadamard) on z built from the RAW gyro (open loop, no bias
feedback). A constant or slowly drifting gyro bias makes a ramp in z; second differences cancel it. One seed per
scenario-segment regime, batch; constants C2(H) and the noise factor computed numerically on synthetic series."""
import numpy as np
import step6_quad_attitude_estimator as S
dt = S.DT; rng = np.random.default_rng(1)
def had(x, s):
    m = x[: len(x) // s * s].reshape(-1, s).mean(1); return np.mean((m[2:] - 2 * m[1:-1] + m[:-2]) ** 2)
# calibrate: white unit-variance noise -> factor; Wiener (lambda = 1/s) -> C2(1/2)
wn = rng.standard_normal(2_000_000); wi = np.cumsum(rng.standard_normal(2_000_000)) * np.sqrt(dt)
for s in (64, 1024): print(f'calib s={s}: white had*l/(var*dt) = {had(wn, s)*s*dt/dt:.3f}   wiener had/(l) = {had(wi, s)/(s*dt):.4f}')
print('expected: white 6 (3 x the first-difference 2); Wiener C2(1/2) = 1 (first differences: 2/3)')
for scn, seg in (('A_steady', 0), ('B_throttle', 1), ('B_throttle', 2), ('C_gyrobias', 0)):
    S.SEG_S = 120.0
    d = S.generate(scn, [600]); n = d['nseg']; W, A1 = d['W'][:, 0], d['A1'][:, 0]
    sl = slice(seg * n, (seg + 1) * n)
    z = A1[sl] - (A1[sl][0] + np.cumsum(W[sl]) * dt)
    sc = [2 ** i for i in range(5, 15)]; L = np.array(sc) * dt
    Dh = np.array([had(z, s) for s in sc])
    A3 = np.min(Dh[:3] * L[:3])                         # noise term = 6 Var_eff dt / l = 3 A / l
    sig = Dh - A3 / L; ok = sig > 2 * A3 / L
    p = np.polyfit(np.log(L[ok]), np.log(sig[ok]), 1); H = p[0] / 2; lam = np.exp(p[1])   # C2(1/2) = 1
    s2 = A3 / 6
    tau = (s2 / lam) ** 0.5
    V = S.SCEN[scn]['V'][seg]
    print(f'{scn} seg {seg}: H {H:.2f} lambda/Q_G {lam/S.Q_G:.2f} s2/S0 {s2/S.S0(V):.2f} -> tau {tau:.2f} s'
          f' (sqrt(S0/Q_G) = {np.sqrt(S.S0(V)/S.Q_G):.2f}); scales used {L[ok].round(2).tolist()}')
