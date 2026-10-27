"""EXPLORATORY diagnosis after the step-6 verdict (not pre-registered): why does the one-stream identification
put tau 20x too short? Batch block-mean structure function of z on one seed of A_steady."""
import numpy as np, sys
sys.argv = ['x']
import step6_quad_attitude_estimator as S
d = S.generate('A_steady', [600])
W, A1, TH = d['W'][:, 0], d['A1'][:, 0], d['TH'][:, 0]
dt = S.DT
# z with the raw gyro (no bias feedback): slow part = gyro angle random walk only
Th = A1[0] + np.cumsum(W) * dt
z = A1 - Th
scales = [2 ** i for i in range(4, 16)]
print('true: A = 2 S0 =', 2 * S.S0(1.0), ' lambda C(1/2) = ', S.Q_G * 2 / 3)
print(' l(s)      D(l)         D*l         A/l (true)   signal part (true)')
for s in scales:
    m = z[: len(z) // s * s].reshape(-1, s).mean(1)
    D = np.mean(np.diff(m) ** 2); l = s * dt
    print(f'{l:8.3f}  {D:11.3e}  {D*l:11.3e}  {2*S.S0(1.0)/l:11.3e}  {S.Q_G*2/3*l:11.3e}')
