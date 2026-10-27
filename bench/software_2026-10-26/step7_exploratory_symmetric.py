"""EXPLORATORY (after step 7a): is T_sym = Q^-1 # Sigma the best SYMMETRIC window, and does T_K satisfy T Q T^T = Sigma?"""
import numpy as np
from scipy.optimize import minimize
import step7_axis_coupling as S7
rng = np.random.default_rng(11); rows = []
for _ in range(200):
    A = rng.standard_normal((2, 2)); B = rng.standard_normal((2, 2))
    Q = A @ A.T + 0.1 * np.eye(2); Sg = B @ B.T + 0.1 * np.eye(2)
    Td, Ts, Tk = S7.windows(Q, Sg)
    bal_k = np.linalg.norm(Tk @ Q @ Tk.T - Sg) / np.linalg.norm(Sg); bal_s = np.linalg.norm(Ts @ Q @ Ts - Sg) / np.linalg.norm(Sg)
    def f(p):
        L = np.array([[np.exp(p[0]), 0], [p[1], np.exp(p[2])]]); T = L @ L.T
        try: return S7.trP(T, Q, Sg)[0]
        except Exception: return 1e30
    L0 = np.linalg.cholesky(Ts); p0 = [np.log(L0[0, 0]), L0[1, 0], np.log(L0[1, 1])]
    r = minimize(f, p0, method='Nelder-Mead', options=dict(xatol=1e-10, fatol=1e-14, maxiter=4000))
    js, jk = S7.trP(Ts, Q, Sg)[0], S7.trP(Tk, Q, Sg)[0]
    rows.append((js / jk - 1, r.fun / jk - 1, bal_k, bal_s))
R = np.array(rows)
print('balance residual T Q T^T = Sigma: T_K max %.1e, T_sym max %.1e' % (R[:, 2].max(), R[:, 3].max()))
print('sym/K - 1: median %.3e max %.3e' % (np.median(R[:, 0]), R[:, 0].max()))
print('best symmetric/K - 1: median %.3e max %.3e' % (np.median(R[:, 1]), R[:, 1].max()))
print('T_sym vs best symmetric: (sym - bestsym)/K, median %.3e max %.3e' % (np.median(R[:, 0] - R[:, 1]), (R[:, 0] - R[:, 1]).max()))
