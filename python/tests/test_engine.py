# Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
# Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
"""Unit tests of the Python engine.   python -m unittest discover -s python/tests"""
import os, sys, unittest
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from syntonic_engine import SyntonicEngine, kappa_H, loop_delay, window_floor


class Formulas(unittest.TestCase):
    def test_kappa(self):
        self.assertAlmostEqual(float(kappa_H(0.5)), 1.0, places=12)            # classical square-root law
        self.assertAlmostEqual(float(kappa_H(1.0)), 4.0 ** (-1 / 3), places=12)  # (2 Gamma(3))^(-1/3)

    def test_floor(self):
        self.assertAlmostEqual(loop_delay(0.001, 0.003), 0.0021, places=15)
        self.assertAlmostEqual(window_floor(0.001, 0.003), 0.001 / 0.45, places=15)   # TH41-A wins: 2.22 ms
        self.assertAlmostEqual(window_floor(0.001, 0.008), 0.0046, places=15)         # slow frame: the delay wins


class Engine(unittest.TestCase):
    def run_engine(self, seconds=3.0, noise=1e-3, seed=0):
        rng = np.random.default_rng(seed); dt = 0.001; n = 3
        t = np.arange(0.0, seconds, dt)
        q = np.stack([0.3 * np.sin(2 * np.pi * f * t) for f in (0.5, 1.0, 2.0)], 1)
        eng = SyntonicEngine(q[0], dt, 0.003, w0=np.full(n, 10.0))
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

    def test_deterministic(self):
        _, a, _, _ = self.run_engine(seconds=1.5, seed=3)
        _, b, _, _ = self.run_engine(seconds=1.5, seed=3)
        np.testing.assert_array_equal(a, b)


if __name__ == '__main__':
    unittest.main()
