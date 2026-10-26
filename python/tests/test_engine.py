# Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
# Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
"""Unit tests of the Python engine.   python -m unittest discover -s python/tests"""
import os, sys, unittest
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from syntonic_engine import SyntonicEngine, kappa_H, C_H, loop_delay, window_floor


class Formulas(unittest.TestCase):
    def test_kappa(self):
        self.assertAlmostEqual(float(kappa_H(0.5)), 1.0, places=12)            # classical square-root law
        self.assertAlmostEqual(float(kappa_H(1.0)), 4.0 ** (-1 / 3), places=12)  # (2 Gamma(3))^(-1/3)

    def test_allan_constant(self):
        self.assertAlmostEqual(float(C_H(0.5)), 2.0 / 3.0, places=12)       # Allan variance of a random walk
        self.assertAlmostEqual(float(C_H(1.0)), 1.0, places=12)             # (16 - 4) / (3 * 4)

    def test_floor(self):
        self.assertAlmostEqual(loop_delay(0.001, 0.003), 0.0021, places=15)
        self.assertAlmostEqual(window_floor(0.001, 0.003), 0.001 / 0.45, places=15)   # TH41-A wins: 2.22 ms
        self.assertAlmostEqual(window_floor(0.001, 0.008), 0.0046, places=15)         # slow frame: the delay wins


class Engine(unittest.TestCase):
    def run_engine(self, seconds=3.0, noise=1e-3, seed=0, identification='allan'):
        rng = np.random.default_rng(seed); dt = 0.001; n = 3
        t = np.arange(0.0, seconds, dt)
        q = np.stack([0.3 * np.sin(2 * np.pi * f * t) for f in (0.5, 1.0, 2.0)], 1)
        eng = SyntonicEngine(q[0], dt, 0.003, w0=np.full(n, 10.0), identification=identification)
        th = q[0].copy(); cmds = []
        for k in range(len(t)):
            a = th + rng.normal(0, noise, n); b = th + rng.normal(0, noise, n)
            cmd = eng.step(q[k], a, b)
            th = th + dt * np.clip((cmd - th) / 0.025, -10.0, 10.0)            # a simple servo
            cmds.append(cmd.copy())
        return eng, np.array(cmds), q, th

    def test_window_stays_in_bounds(self):
        eng, cmds, _, _ = self.run_engine()
        self.assertTrue(np.all(np.isfinite(cmds)))
        self.assertTrue(np.all(eng.tau >= eng.floor - 1e-15) and np.all(eng.tau <= 0.5))
        self.assertTrue(np.all((eng.H >= 0.05) & (eng.H <= 1.0)))

    def test_tracks_the_plan(self):
        eng, _, q, th = self.run_engine()
        self.assertLess(np.max(np.abs(th - q[-1])), 0.05)                      # within 0.05 rad at the end

    def test_command_held_between_frames(self):
        eng = SyntonicEngine(np.zeros(2), 0.001, 0.003, w0=np.full(2, 10.0))
        c1 = eng.step(np.array([0.1, 0.2]), np.zeros(2), np.zeros(2), frame=True).copy()
        c2 = eng.step(np.array([0.5, 0.5]), np.zeros(2), np.zeros(2), frame=False).copy()
        np.testing.assert_array_equal(c1, c2)
        c3 = eng.step(np.array([0.5, 0.5]), np.zeros(2), np.zeros(2), frame=True)
        self.assertFalse(np.array_equal(c2, c3))

    def test_raw_identification_is_the_0_1_x_engine(self):
        # reference values computed with syntonic-engine 0.1.1 before the 0.2.0 change (same harness, seed 0)
        eng, cmds, _, _ = self.run_engine(identification='raw')
        self.assertAlmostEqual(float(cmds.sum()), 187.66031755510681, places=9)
        np.testing.assert_allclose(eng.tau, [0.0035389, 0.00222222, 0.00222222], rtol=1e-4)
        np.testing.assert_allclose(eng.H, [1.0, 1.0, 0.99764172], rtol=1e-4)

    def test_bad_identification_name(self):
        with self.assertRaises(ValueError):
            SyntonicEngine(np.zeros(1), 0.001, 0.003, w0=10.0, identification='other')

    def test_deterministic(self):
        _, a, _, _ = self.run_engine(seconds=1.5, seed=3)
        _, b, _, _ = self.run_engine(seconds=1.5, seed=3)
        np.testing.assert_array_equal(a, b)


class Identification(unittest.TestCase):
    """The fit itself, fed a residual with a known structure function and known read noise: fractional Brownian
    motion of index H (lambda_H = 1 in rad^2 s^-2H) plus two white reads. Checks the case the software bench of
    26 Oct 2026 found broken in 0.1.x: read noise large enough that the signal is below the noise floor at
    every lag up to 64 ms (the U' optimum then sits at tens of ms), where the raw fit has no usable lag."""

    @staticmethod
    def fbm(n, H, rng, dt):
        k = np.arange(n + 1, dtype=float)
        g = 0.5 * (np.abs(k + 1) ** (2 * H) - 2 * np.abs(k) ** (2 * H) + np.abs(k - 1) ** (2 * H))
        r = np.concatenate([g, g[-2:0:-1]]); lam = np.maximum(np.fft.fft(r).real, 0.0); M = r.size
        W = rng.standard_normal(M) + 1j * rng.standard_normal(M)
        return dt ** H * np.concatenate([[0.0], np.cumsum((np.fft.ifft(np.sqrt(lam) * W) * np.sqrt(M)).real[:n])[:-1]])

    def identify(self, H, tau_target, identification, n=2 ** 15, seed=1):
        """Feeds the identification alone (no servo) and returns the time-medians of H and of the window it would
        set, over the second half of the run."""
        dt = 0.001; rng = np.random.default_rng(seed)
        kH = float(kappa_H(H)); s2 = (tau_target / kH) ** (2 * H + 1)        # noise PSD placing the U' optimum
        sig = np.sqrt(2 * s2 / dt)
        x = self.fbm(n, H, rng, dt); a = x + sig * rng.standard_normal(n); b = x + sig * rng.standard_normal(n)
        eng = SyntonicEngine(np.zeros(1), dt, 0.003, w0=1e9, identification=identification, tau_max=1.0)
        Hs, taus = [], []
        for k in range(n):
            d = (a[k] - b[k]) / np.sqrt(2.0); eng.n2 += eng.a_w * (d * d - eng.n2)
            r = np.array([0.5 * (a[k] + b[k])])
            fit = eng._identify_raw(r) if identification == 'raw' else eng._identify_allan(r)
            if fit is not None and fit[2][0] and k > n // 2:
                Hn, lam, _ = fit; s2_hat = 0.5 * eng.n2[0] * dt
                taus.append(float(kappa_H(Hn[0]) * (s2_hat / lam[0]) ** (1.0 / (2 * Hn[0] + 1)))); Hs.append(float(Hn[0]))
            eng.k += 1
        return float(np.median(Hs)), float(np.median(taus))

    def test_allan_finds_the_window_under_large_read_noise(self):
        for H in (0.5, 0.85):
            Hn, tau_hat = self.identify(H, 0.05, 'allan')                    # optimum at 50 ms; signal under the floor up to ~s
            self.assertTrue(0.05 / 1.5 < tau_hat < 0.05 * 1.5, (H, tau_hat))
            self.assertGreater(Hn, 0.3)                                      # not collapsed to the clip floor

    def test_raw_fails_there_and_works_with_small_noise(self):
        Hn, tau_hat = self.identify(0.5, 0.05, 'raw')                         # the 0.1.x limit: H collapses, window off
        self.assertLess(Hn, 0.1)
        self.assertFalse(0.05 / 1.5 < tau_hat < 0.05 * 1.5, tau_hat)
        Hn, tau_hat = self.identify(0.5, 0.005, 'raw')                        # fine encoder regime: raw is right
        self.assertTrue(0.005 / 1.5 < tau_hat < 0.005 * 1.5, tau_hat)
        self.assertGreater(Hn, 0.3)


if __name__ == '__main__':
    unittest.main()
