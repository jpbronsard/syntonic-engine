"""Conforming estimator of the syntonic window tau2 on a log-domain level series.

Implements EXPA_PREREG.md §1.3 with the amendments produced by the offline
falsifier (§1.3-bis): the original per-frame floor max(rF^2-nF,0) and the
stride W = tau are provably biased/bistable when per-frame SNR <~ 1
(SNR = W^2/(3 tau^2)).  Amended design, all constants class-tagged:

  level  x = (xA + xB)/2,  disagreement d = (xA - xB)/sqrt(2)
  n2  = EMA_{beta_s}(d^2),  beta_s = exp(-1/(2 tau))          [per-step anchor]
  stride W = clip(round(W_MULT * tau), 10, 5000)              [W_MULT = 3,
        STRUCTURAL: per-frame SNR = W^2/(3 tau^2) >= 3]
  per-frame paired, unbiased innovation estimate (noise reference d̄ is the
  SAME frame's disagreement mean, independent of the level channel):
        q_j = (rF^2 - dbar^2) / (C_FRAME * W),   C_FRAME = 2/3
        [STRUCTURAL: exact consecutive-block-mean factor of a random walk]
  q   = max( EMA_{beta_f}(q_j), q_floor ),  beta_f = exp(-1/N_MEM), N_MEM = 8
        [SAFETY-class averaging depth: memory = 8 strides = 24 tau]
  tau <- clip( sqrt(n2/(2 q)), rate-capped at |dlog tau| <= ln 2 per frame )
        [RATE_CAP: guard-rail class, kills the floor/cap bistability]

H8 Adam-style bias correction (running normaliser c) on every EMA.
Framework-free (numpy); vectorised over n streams sharing one W
(median of tau across streams).
"""
from __future__ import annotations
import numpy as np

PSI1_HALF = np.pi ** 2 / 2.0  # trigamma(1/2): Var(ln chi^2_1)


class ConformingTau:
    W_MULT = 3.0          # STRUCTURAL: stride/horizon ratio (per-frame SNR >= 3)
    C_FRAME = 2.0 / 3.0   # STRUCTURAL: RW consecutive block-mean factor
    N_MEM = 8.0           # SAFETY: innovation memory, in frames
    RATE_CAP = np.log(2.) # SAFETY: max |dlog tau| per frame (guard-rail class)

    def __init__(self, n: int = 1, tau0: float = 100.0,
                 w_bounds=(10, 5000), tau_bounds=(2.0, 2.0e4),
                 q_floor: float = 1e-12):
        self.n = int(n)
        self.tau = np.full(self.n, float(tau0))
        self.wmin, self.wmax = w_bounds
        self.tmin, self.tmax = tau_bounds
        self.q_floor = float(q_floor)
        self.W = self._stride()
        self.m_n2 = np.zeros(self.n); self.c_n2 = np.zeros(self.n)
        self.m_q = np.zeros(self.n);  self.c_q = np.zeros(self.n)
        self.m_nF = np.zeros(self.n); self.c_nF = np.zeros(self.n)  # telemetry
        self.x_sum = np.zeros(self.n); self.d_sum = np.zeros(self.n); self.k = 0
        self.prev_xbar = None
        self.t = 0
        self.frames = 0
        self.last_frame = None

    def _stride(self) -> int:
        return int(np.clip(round(self.W_MULT * float(np.median(self.tau))),
                           self.wmin, self.wmax))

    @property
    def n2_hat(self):
        return self.m_n2 / np.maximum(self.c_n2, 1e-300)

    @property
    def q_hat(self):
        return np.maximum(self.m_q / np.maximum(self.c_q, 1e-300), self.q_floor)

    def burned_in(self):
        """Boolean mask: t >= 5*tau (H8 burn-in, prereg §1.3.5)."""
        return self.t >= 5.0 * self.tau

    def update(self, xA, xB):
        xA = np.atleast_1d(np.asarray(xA, dtype=np.float64))
        xB = np.atleast_1d(np.asarray(xB, dtype=np.float64))
        d = (xA - xB) / np.sqrt(2.0)
        x = 0.5 * (xA + xB)

        bs = np.exp(-1.0 / (2.0 * self.tau))
        self.m_n2 = bs * self.m_n2 + (1.0 - bs) * d * d
        self.c_n2 = bs * self.c_n2 + (1.0 - bs)

        self.x_sum += x
        self.d_sum += d
        self.k += 1
        self.t += 1
        if self.k >= self.W:
            self._close_frame()
        return self.tau

    def _close_frame(self):
        W = self.k
        xbar = self.x_sum / W
        dbar = self.d_sum / W
        if self.prev_xbar is not None:
            rF = xbar - self.prev_xbar
            bf = np.exp(-1.0 / self.N_MEM)
            # unbiased per-frame innovation (paired noise reference, same frame)
            qj = (rF * rF - dbar * dbar) / (self.C_FRAME * W)
            self.m_q = bf * self.m_q + (1.0 - bf) * qj
            self.c_q = bf * self.c_q + (1.0 - bf)
            self.m_nF = bf * self.m_nF + (1.0 - bf) * dbar * dbar
            self.c_nF = bf * self.c_nF + (1.0 - bf)
            tau_raw = np.sqrt(self.n2_hat / (2.0 * self.q_hat))
            step = np.clip(np.log(np.maximum(tau_raw, self.tmin) / self.tau),
                           -self.RATE_CAP, self.RATE_CAP)
            self.tau = np.clip(self.tau * np.exp(step), self.tmin, self.tmax)
            self.frames += 1
            self.last_frame = dict(
                t=self.t, W=W, rF=rF.copy(),
                nF=(self.m_nF / np.maximum(self.c_nF, 1e-300)).copy(),
                q=self.q_hat.copy(), n2=self.n2_hat.copy(), tau=self.tau.copy(),
            )
        self.prev_xbar = xbar
        self.x_sum[:] = 0.0
        self.d_sum[:] = 0.0
        self.k = 0
        self.W = self._stride()


def log_level_pair(vA, vB, eps: float = 1e-12):
    """Log-domain transform of a pair of half-batch variance measurements."""
    return np.log(np.asarray(vA) + eps), np.log(np.asarray(vB) + eps)
