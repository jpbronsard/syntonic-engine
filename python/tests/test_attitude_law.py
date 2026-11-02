# Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
# Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
"""Unit tests of the Riccati attitude law (0.3.0).   python -m unittest discover -s python/tests
Its exact agreement with the confirmation's column LAW_R is checked by bench/goal_2026-11-02/package_gate.py, which
needs the frozen flight harness; these tests need numpy only. The synthetic world below is the law's own model (white
tilt noise, white gyro noise, one constant gyro bias per level) with facts of the size of the X500's."""
import os, sys, unittest
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from syntonic_engine import RiccatiAttitudeLaw

DT, DTA, NSUB, T_PRE, QB = 1e-3, 4e-3, 4, 5.0, 9.2e-11
LAM = {'idle': 9e-9, 'low': 3.3e-8, 'half': 1.2e-6}                   # gyro angle random walk (rad^2/s)
S2 = {'idle': 1.9e-5, 'low': 2.2e-5, 'half': 3.5e-3}                  # tilt noise density (rad^2 s)
TREC = {'low': 40.0, 'half': 8.0}
LEVELS = {r: dict(s2=S2[r], lam=LAM[r], **({'T_rec': TREC[r]} if r in TREC else {})) for r in LAM}
SD_G = {r: np.sqrt(LAM[r] / DT) for r in LAM}                          # per gyro sample (rad/s)
SD_A = {r: np.sqrt(S2[r] / DTA) for r in S2}                           # per accelerometer sample (rad)
BIAS = {'idle': 6e-4, 'low': 6e-4 + 2.7e-5, 'half': 6e-4 + 2.0e-4}     # gyro bias per level (rad/s)


def law():
    return RiccatiAttitudeLaw(LEVELS, q_b=QB, T_pre=T_PRE, dt=DT, dt_accel=DTA)


def fly(est, segs=(('low', 15.0), ('half', 8.0), ('low', 7.0)), seed=0, theta0=0.01, tilt_scale=None):
    """plant at rest at theta0; 5 s on the ground at idle, arming calibration = the pre-roll gyro mean, then the
    segments. tilt_scale = (level, t0, t1, factor) multiplies the tilt noise in flight time [t0, t1). Returns the
    estimation errors after arming, those of a 1 s average of the tilt, and the calibration."""
    rng = np.random.default_rng(seed)
    plan = [('idle', T_PRE)] + list(segs); kpre = int(round(T_PRE / DT)); gsum = cal = 0.0; k = 0
    errs, avg, buf = [], [], []
    for level, length in plan:
        for _ in range(int(round(length / DT))):
            g_raw = BIAS[level] + rng.normal(0.0, SD_G[level]); a = None
            if k % NSUB == 0:
                sd = SD_A[level]
                if tilt_scale and level == tilt_scale[0] and tilt_scale[1] <= k * DT - T_PRE < tilt_scale[2]:
                    sd = sd * tilt_scale[3]
                a = theta0 + rng.normal(0.0, sd); buf = (buf + [a])[-250:]
            th = est.step(g_raw - cal, a, level=level)
            if k < kpre:
                gsum += g_raw
                if k == kpre - 1:
                    cal = gsum / kpre; est.arm(cal)
            else:
                errs.append(th - theta0); avg.append(np.mean(buf) - theta0)
            k += 1
    return np.array(errs), np.array(avg), cal


class Facts(unittest.TestCase):
    def test_missing_facts(self):
        with self.assertRaises(ValueError):
            RiccatiAttitudeLaw({'idle': dict(s2=1e-5, lam=1e-8), 'low': dict(s2=1e-5, lam=1e-8)}, q_b=QB, T_pre=T_PRE)
        with self.assertRaises(ValueError):
            RiccatiAttitudeLaw({'idle': dict(s2=1e-5)}, q_b=QB, T_pre=T_PRE)
        with self.assertRaises(ValueError):
            RiccatiAttitudeLaw(LEVELS, q_b=QB, T_pre=T_PRE, ground='hover')
        with self.assertRaises(ValueError):
            law().step(0.0, None)                                                  # the level is needed

    def test_harness_signature(self):
        a, b = law(), law()
        self.assertEqual(a.n, 1)
        for k in range(40):
            acc = 0.01 if k % 4 == 0 else None
            self.assertEqual(a.step(0.002, acc, 1e-9, k * DT, 'idle', k), b.step(0.002, acc, level='idle'))
        a.arm(0.002); b.arm()
        np.testing.assert_array_equal(a.P, b.P)


class Recursion(unittest.TestCase):
    def test_stationary_alignment_before_arming(self):
        est = law(); rng = np.random.default_rng(1)
        for k in range(5000):                                                      # a large uncalibrated gyro rate
            est.step(0.05, 0.01 + rng.normal(0, SD_A['idle']) if k % NSUB == 0 else None, level='idle')
        self.assertLess(abs(est.theta - 0.01), 0.01)                               # the tilt alone, not 0.25 rad of gyro
        np.testing.assert_array_equal(est.x[1:], 0.0); np.testing.assert_array_equal(est.P[1:, :], 0.0)
        self.assertAlmostEqual(est.thr, 0.05 * 5.0, places=9)                      # the gyro feeds the noise reading only

    def test_arming_covariance(self):
        est = law()
        for k in range(400):
            est.step(0.0, 0.0 if k % NSUB == 0 else None, level='idle')
        est.arm()
        want = np.full((3, 3), LAM['idle'] / T_PRE) + np.diag([0.0, LAM['low'] / TREC['low'], LAM['half'] / TREC['half']])
        np.testing.assert_allclose(est.P[1:, 1:], want, rtol=1e-15)
        np.testing.assert_array_equal(est.P[0, 1:], 0.0); np.testing.assert_array_equal(est.P[1:, 0], 0.0)

    def test_prediction_is_F_P_Ft_plus_Q(self):
        est = law(); est.arm(); rng = np.random.default_rng(2)
        M = rng.normal(size=(4, 4)); est.P = M @ M.T * 1e-6; est.x = rng.normal(size=4) * 1e-3
        P0, x0 = est.P.copy(), est.x.copy(); g = 0.003
        est.step(g, None, level='half')                                            # the bias of 'half' is x[3]
        F = np.eye(4); F[0, 3] = -DT
        Q = np.zeros((4, 4)); Q[0, 0] = LAM['half'] * DT; Q[3, 3] = QB * DT
        np.testing.assert_allclose(est.P, F @ P0 @ F.T + Q, rtol=1e-12, atol=1e-24)
        np.testing.assert_allclose(est.x, x0 + np.array([(g - x0[3]) * DT, 0, 0, 0]), rtol=1e-15)

    def test_update_is_the_kalman_update(self):
        est = law()
        for k in range(2000):
            est.step(0.0, 0.003 if k % NSUB == 0 else None, level='idle')
        est.arm(); est.P[1:, 0] = est.P[0, 1:] = [1e-9, 2e-9, -1e-9]               # some cross terms
        x0, P0 = est.x.copy(), est.P.copy()
        est.step(0.0, 0.004, level='low')                                          # an accelerometer tick at 'low'
        F = np.eye(4); F[0, 2] = -DT
        Pp = F @ P0 @ F.T; Pp[0, 0] += LAM['low'] * DT; Pp[2, 2] += QB * DT
        xp = x0.copy(); xp[0] -= x0[2] * DT
        R = max(S2['low'], est.noise_reading('low')) / DTA
        K = Pp[:, 0] / (Pp[0, 0] + R)
        np.testing.assert_allclose(est.x, xp + K * (0.004 - xp[0]), rtol=1e-12)
        np.testing.assert_allclose(est.P, Pp - np.outer(K, Pp[0, :]), rtol=1e-10, atol=1e-30)


class Flight(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.est = law(); cls.err, cls.avg, cls.cal = fly(cls.est)

    def test_better_than_a_tilt_average(self):
        rms = np.sqrt(np.mean(self.err ** 2)); rms_avg = np.sqrt(np.mean(self.avg ** 2))
        self.assertLess(rms, 0.2 * rms_avg, (rms, rms_avg))                        # 0.03 to 0.08 on seeds 0-7

    def test_consistent_with_its_own_covariance(self):
        est = self.est
        z = [(est.theta - 0.01) / np.sqrt(est.P[0, 0])]
        z += [(est.bias(r) - (BIAS[r] - self.cal)) / np.sqrt(est.P[1 + j, 1 + j]) for j, r in enumerate(LEVELS)]
        self.assertTrue(np.all(np.abs(z) < 4.0), z)                                # within 1.8 on seeds 0-7

    def test_covariance_stays_a_covariance(self):
        P = self.est.P
        np.testing.assert_allclose(P, P.T, rtol=1e-9, atol=1e-30)
        self.assertTrue(np.all(np.linalg.eigvalsh(0.5 * (P + P.T)) > 0))

    def test_live_reading_raises_R_and_keeps_one_memory_per_level(self):
        est = law(); fly(est, segs=(('low', 10.0),))
        self.assertLess(est.noise_reading('low'), 2.0 * S2['low'])                 # white noise at the fact's level
        est = law(); fly(est, segs=(('low', 10.0),), tilt_scale=('low', 4.0, 10.0, 3.0))
        loud = est.noise_reading('low')
        self.assertGreater(loud, 3.0 * S2['low'])                                  # nine times the variance
        half0 = est.noise_reading('half'); rng = np.random.default_rng(5)
        for k in range(5000):                                                      # 5 s at half throttle
            est.step(rng.normal(0, SD_G['half']), 0.01 + rng.normal(0, SD_A['half']) if k % NSUB == 0 else None,
                     level='half')
        self.assertEqual(est.noise_reading('low'), loud)                           # low's memory left alone
        self.assertNotEqual(est.noise_reading('half'), half0)

    def test_regression(self):
        # values of the 0.3.0 package on this harness, seed 0
        self.assertAlmostEqual(float(np.sum(self.err ** 2)) / REF['sum_err2'], 1.0, places=9)
        self.assertAlmostEqual(self.est.theta / REF['theta'], 1.0, places=9)
        self.assertAlmostEqual(self.est.bias('low') / REF['b_low'], 1.0, places=9)
        self.assertAlmostEqual(float(self.est.P[0, 0]) / REF['p00'], 1.0, places=9)


REF = dict(sum_err2=0.15556646333379984, theta=0.013047539565128654, b_low=-2.33821250312898e-05,
           p00=2.9211834774105486e-06)

if __name__ == '__main__':
    unittest.main()
