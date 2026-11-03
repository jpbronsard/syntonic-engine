# Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
# Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
"""The loop law's gain recomputation as the C++ and Rust laws do it, with the two speed-ups of 0.3.2, and ScanLoopLaw,
the Python law that uses it. This is the Python reference of cpp/SyntonicLoopLaw.h and rust/src/loop_law.rs.

SyntonicLoopLaw (loop_law.py, the frozen law) chooses every 500 ticks, per joint, the robust candidate (Kp, Ki) with the
smallest predicted loop error J = sum over w of Phi(w) |S(w)|^2 dW / pi, plus r MTr; numpy computes all the J by a
matrix product. The C++ and Rust laws find the same candidate by a scan: candidate after candidate in index order,
each J summed over the 1500 frequencies in ascending order, the first smallest kept. 0.3.2 makes that scan cheaper in
two ways, without changing the chosen gains:

  1. Pruning (prune=True). Every term of the sum is positive, so a partial sum never decreases, in floating point as
     in exact arithmetic. The candidates in use (one per joint) are evaluated first. A candidate's sum is checked every
     50 frequencies: once, for every joint, its partial sum plus r MTr exceeds the smaller of the J of the gains in
     use and the best J found so far, times (1 + DELTA), it cannot be chosen and the rest of its sum is skipped.
  2. The quadratic form (quad=True). With L = (Kp + Ki / jW) G,
         |1 + L|^2 = 1 + 2 Kp Re G + Kp^2 |G|^2 + 2 Ki Im G / W + Ki^2 |G|^2 / W^2 = P + Ki (Q + Ki R),
     P = 1 + Kp (2 Re G + Kp |G|^2) per Kp row and frequency, Q = 2 Im G / W and R = |G|^2 / W^2 per frequency. Then
     |S|^2 = 1 / (P + Ki (Q + Ki R)) costs four operations and a division instead of ten and a division. The two
     forms differ by rounding only: for a robust candidate (|1 + L| >= 1/2, and |G| <= 1) the five terms add up to at
     most 59 |1 + L|^2, so |S|^2 changes by at most about 1e-13 relative, and J, a sum of 1500 positive terms, by at
     most about 5e-13. A choice is certified when the runner-up's J exceeds the best's by more than DELTA = 1e-9 in
     relative terms (the margin is kept by the pruning, which never drops a candidate within DELTA of the best): then
     the exact J orders them the same way. Otherwise the choice is made again exactly (a fallback, counted).

DELTA is a declared numerical constant, three orders of magnitude above the rounding bound and three below the
smallest gap seen between the best and the next J on the replays of 0.3.1 (3.5e-6). The arithmetic below is that of
the C++ and Rust laws (the spectrum with dW folded in, the same expressions, the sums in the same order); the fallback
of ScanLoopLaw is LoopLaw.gains, the frozen law's own matrix product, so ScanLoopLaw chooses the gains of
SyntonicLoopLaw at every recomputation (checked by bench/goal_2026-11-02/package_gate.py --scan and
bench/port_2026-11-02/scan_check.py). Written 2 Nov 2026 (0.3.2).
"""
import numpy as np
from .loop_law import SyntonicLoopLaw, K_gamma, loop_delay

DELTA = 1e-9      # a choice is certified when the runner-up's J exceeds the best's by more than this fraction
FB = 50           # frequencies per block: the pruning test runs at the end of each block


class LoopScan:
    """The robust candidates of a LoopLaw and the per-frequency constants of the C++/Rust recomputation."""

    def __init__(self, servo_tau, D, LL):
        self.LL = LL
        W = LL.W; nw = len(W); self.nw, self.nb = nw, nw // FB
        dW = np.empty(nw); dW[0] = W[1] - W[0]; dW[-1] = W[-1] - W[-2]; dW[1:-1] = (W[2:] - W[:-2]) / 2.0
        self.dWf = dW / (2 * np.pi) * 2; self.W2 = W * W; self.iW = 1.0 / W
        cr = np.cos(W * D); ci = -np.sin(W * D); di = W * servo_tau; m = 1.0 + di * di
        self.Gr = (cr + ci * di) / m; self.Gi = (ci - cr * di) / m                      # exp(-s D) / (1 + s ts)
        self.Gr2 = 2.0 * self.Gr; self.g2 = self.Gr * self.Gr + self.Gi * self.Gi       # the quadratic form
        self.Q = 2.0 * (self.Gi * self.iW); self.R = (self.iW * self.iW) * self.g2
        nki = 60
        self.kp = LL.KP.reshape(-1, nki)[:, 0].copy(); self.ki = LL.KI[:nki].copy()
        self.okc = np.flatnonzero(LL.ok)
        self.rows = [(i, self.okc[self.okc // nki == i] % nki) for i in range(len(self.kp))]
        self.rows = [(i, k) for i, k in self.rows if len(k)]
        self.index = {(kp, ki): c for c, (kp, ki) in enumerate(zip(LL.KP, LL.KI))}
        self.MTr = np.zeros(len(LL.KP))                                                 # sum of |T|^2 dW / pi, in order
        for i, kis in self.rows:
            kiw = self.ki[kis][:, None] * self.iW[None, :]
            Lr = self.kp[i] * self.Gr + kiw * self.Gi; Li = self.kp[i] * self.Gi - kiw * self.Gr
            l2 = Lr * Lr + Li * Li; s2 = 1.0 / ((1.0 + Lr) * (1.0 + Lr) + Li * Li)
            self.MTr[i * nki + kis] = np.cumsum(l2 * s2 * self.dWf, axis=1)[:, -1]

    def s2(self, i, kis, quad):
        """|S|^2 of the candidates (row i, Ki numbers kis) at every frequency, by the quadratic form or from L"""
        kp, ki = self.kp[i], self.ki[kis][:, None]
        if quad:
            P = 1.0 + kp * (self.Gr2 + kp * self.g2)
            return 1.0 / (P[None, :] + ki * (self.Q[None, :] + ki * self.R[None, :]))
        kiw = ki * self.iW[None, :]
        Lr = kp * self.Gr + kiw * self.Gi; Li = kp * self.Gi - kiw * self.Gr
        return 1.0 / ((1.0 + Lr) * (1.0 + Lr) + Li * Li)

    def choose(self, A, gam, Oc, r, warm=(), prune=True, quad=True, delta=DELTA):
        """The scan of the C++/Rust laws. warm: the candidates in use. Returns (candidate per joint or None when the
        choice must be made again exactly, statistics)."""
        n = len(A); nki = len(self.ki); nw, nb = self.nw, self.nb
        PhiW = (A / K_gamma(gam))[:, None] * (self.W2[None, :] + (Oc * Oc)[:, None]) ** (-(1 + gam) / 2)[:, None] \
            * self.dWf[None, :]
        st = dict(points=0, total=len(self.okc) * nw, rows=0, certified=True, margin=np.inf, guard=False)
        if np.isnan(PhiW).any() or np.isnan(r).any():                                  # pruning needs positive terms
            st['guard'] = True; st['certified'] = False
            return None, st
        b1 = np.full(n, np.inf); b2 = np.full(n, np.inf); c1 = np.zeros(n, int)
        warmJ = {}
        if prune:
            for c in warm:
                if c in warmJ or not self.LL.ok[c]:
                    continue
                s2 = self.s2(c // nki, np.array([c % nki]), quad)[0]
                warmJ[c] = np.cumsum(PhiW * s2[None, :], axis=1)[:, -1] + r * self.MTr[c]
                st['points'] += nw; st['rows'] += quad
        bw = np.min(np.array(list(warmJ.values())), axis=0) if warmJ else np.full(n, np.inf)
        for i, kis in self.rows:
            cs = i * nki + kis
            s2 = self.s2(i, kis, quad); st['rows'] += quad
            cum = np.cumsum(PhiW[:, None, :] * s2[None, :, :], axis=2)                 # (joint, candidate, frequency)
            rM = r[:, None] * self.MTr[None, cs]
            for t, c in enumerate(cs):
                if c in warmJ:
                    J = warmJ[c]
                else:
                    if prune:
                        thr = np.minimum(bw, b1) * (1.0 + delta)
                        out = (cum[:, t, FB - 1::FB] + rM[:, t][:, None]) > thr[:, None]
                        first = np.where(out.any(1), out.argmax(1), nb)
                        stop = int(first.max())
                        if stop < nb:
                            st['points'] += (stop + 1) * FB
                            continue
                    st['points'] += nw
                    J = cum[:, t, -1] + rM[:, t]
                for j in range(n):
                    v = J[j]
                    if v < b1[j] or (v != v and b1[j] == b1[j]):
                        b2[j] = b1[j]; b1[j] = v; c1[j] = c
                    elif v < b2[j]:
                        b2[j] = v
        with np.errstate(invalid='ignore', divide='ignore'):
            st['margin'] = float(np.min((b2 - b1) / np.abs(b1)))
        if quad and not np.all(b2 > b1 * (1.0 + delta)):
            st['certified'] = False
            return None, st
        return c1, st


_SCANS = {}


class _ScanGains:
    """stands for the law's LoopLaw: gains() by the scan, warm-started from the gains in use; the rest is LoopLaw's"""

    def __init__(self, law, LL):
        self._law, self._LL = law, LL

    def __getattr__(self, name):
        return getattr(self._LL, name)

    def gains(self, A, gam, Oc, r):
        law = self._law; sc = law.scan
        warm = [sc.index[(kp, ki)] for kp, ki in zip(law.Kp, law.Ki) if (kp, ki) in sc.index]
        c, st = sc.choose(A, gam, Oc, r, warm, law.prune, law.quad, law.delta)
        law.stats.append(st)
        if c is None:                      # not certified, or a NaN: the frozen law's own computation
            law.fallbacks += 1
            return self._LL.gains(A, gam, Oc, r)
        return self._LL.KP[c], self._LL.KI[c]


class ScanLoopLaw(SyntonicLoopLaw):
    """SyntonicLoopLaw with the gains chosen by the scan of the C++ and Rust laws (0.3.2): the same gains at every
    recomputation, found with pruning (prune) and the quadratic form (quad). stats: one dict per recomputation (points
    evaluated, of total; rows of P computed; certified; the observed margin between the best and the next J);
    fallbacks: recomputations made again by LoopLaw.gains. delta: the declared margin DELTA (a larger value costs
    time and never changes the gains; tests use 1e3 to make every recomputation fall back). Slower than
    SyntonicLoopLaw in Python: it is the reference of the C++ and Rust code, not a speed-up of numpy."""

    def __init__(self, q0, w0, dt=1e-3, frame=3e-3, servo_tau=0.025, dead=0.0, lead=None, prune=True, quad=True,
                 delta=DELTA, **kw):
        if not delta >= 1e-10:
            raise ValueError('delta must be at least 1e-10, far above the rounding bound of the quadratic form')
        super().__init__(q0, w0, dt=dt, frame=frame, servo_tau=servo_tau, dead=dead, lead=lead, **kw)
        D = loop_delay(dt, frame) + dead; key = (servo_tau, round(D, 12))
        if key not in _SCANS:
            _SCANS[key] = LoopScan(servo_tau, D, self.LL)
        self.scan = _SCANS[key]; self.prune, self.quad, self.delta = bool(prune), bool(quad), float(delta)
        self.stats = []; self.fallbacks = 0
        self.LL = _ScanGains(self, self.LL)
