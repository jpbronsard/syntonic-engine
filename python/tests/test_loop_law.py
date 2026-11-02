# Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
# Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
"""Unit tests of the syntonic loop law (0.3.0).   python -m unittest discover -s python/tests"""
import os, sys, hashlib, importlib.util, unittest
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from syntonic_engine import SyntonicLoopLaw

FROZEN = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'bench', 'goal_2026-11-02', 'servo',
                      'syntonic_loop_law.py')
FROZEN_SHA = 'a7a234ccefb2c9bb963b87b9a3268493644127bfd43d44c387e1411b3910b2b3'


def closed_loop(law_cls, seconds=4.0, seed=0, correct=True):
    """three joints, 1 kHz loop, 3 ms frames, a first-order 25 ms servo with a speed limit and 2 ms of dead time, an
    output disturbance (slow random walk + a 0.01 rad step at 1.5 s), two noisy reads per tick. correct=False sends the
    plan and the lead only. Returns the law and the joint error per tick."""
    rng = np.random.default_rng(seed); dt, frame, ts, dead_ticks, n = 1e-3, 3e-3, 0.025, 2, 3
    w0 = np.full(n, 6.0); t = np.arange(0.0, seconds, dt); f = np.array([0.4, 0.7, 1.1])
    q = 0.3 * np.sin(2 * np.pi * f * t[:, None]); qd = 0.6 * np.pi * f * np.cos(2 * np.pi * f * t[:, None])
    dist = np.cumsum(rng.normal(0, 2e-4, (len(t), n)), 0) + np.where(t[:, None] >= 1.5, 0.01, 0.0)
    law = law_cls(q[0], w0, dt=dt, frame=frame, servo_tau=ts, dead=dead_ticks * dt)
    th = q[0].copy(); latched = q[0].copy(); pipe = [q[0].copy()] * (dead_ticks + 1); err = []
    for k in range(len(t)):
        y = th + dist[k]
        a = y + rng.normal(0, 2e-4, n); b = y + rng.normal(0, 2e-4, n)
        send = k % 3 == 0
        cmd = law.step(q[k], qd[k], a, b, send)
        if send:
            latched = (cmd if correct else q[k] + law.L * qd[k]).copy()
        pipe.append(latched.copy()); u = pipe.pop(0)
        th = th + dt * np.clip((u - th) / ts, -w0, w0)
        err.append(q[k] - y)
    return law, np.array(err)


class LoopLaw(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.law, cls.err = closed_loop(SyntonicLoopLaw)

    def test_corrects_what_the_lead_cannot(self):
        _, err_ff = closed_loop(SyntonicLoopLaw, correct=False)
        rms = np.sqrt(np.mean(self.err[-2000:] ** 2, 0)); rms_ff = np.sqrt(np.mean(err_ff[-2000:] ** 2, 0))
        self.assertTrue(np.all(rms < 0.5 * rms_ff), (rms, rms_ff))                 # 2.6 to 8.7 times smaller here

    def test_gains_are_computed_and_robust(self):
        law = self.law
        self.assertTrue(law.started)
        self.assertTrue(np.all(np.isfinite(law.Kp)) and np.all(np.isfinite(law.Ki)))
        chosen = [np.flatnonzero((law.LL.KP == kp) & (law.LL.KI == ki)) for kp, ki in zip(law.Kp, law.Ki)]
        self.assertTrue(all(len(c) == 1 and law.LL.ok[c[0]] for c in chosen))   # on the grid, phase margin and |S|
        self.assertTrue(np.all((law.gam >= 1.0) & (law.gam <= 3.0)))
        np.testing.assert_allclose(law.R, 2e-8, rtol=0.2)                          # read-pair noise: sd^2 / 2
        self.assertAlmostEqual(law.L, 0.025 + 0.0015 + 0.0005 + 1e-4, places=15)  # servo_tau + loop delay

    def test_regression(self):
        # values of the 0.3.0 package (identical to the frozen law, next test) on this harness, seed 0
        self.assertAlmostEqual(float(np.sum(self.err ** 2)), 0.05607977970139743, places=12)
        np.testing.assert_allclose(self.law.Kp, [4.48708032, 4.48708032, 0.83098111], rtol=1e-8)
        np.testing.assert_allclose(self.law.Ki, [93.35905456, 93.35905456, 5.0], rtol=1e-8)

    @unittest.skipUnless(os.path.exists(FROZEN), 'the frozen law of bench/goal_2026-11-02 is not here')
    def test_is_the_frozen_law(self):
        with open(FROZEN, 'rb') as fh:
            self.assertEqual(hashlib.sha256(fh.read()).hexdigest(), FROZEN_SHA)
        spec = importlib.util.spec_from_file_location('frozen_syntonic_loop_law', FROZEN)
        mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
        law_f, err_f = closed_loop(mod.SyntonicLoopLaw)
        np.testing.assert_array_equal(err_f, self.err)                              # bit for bit
        np.testing.assert_array_equal(law_f.Kp, self.law.Kp); np.testing.assert_array_equal(law_f.Ki, self.law.Ki)


if __name__ == '__main__':
    unittest.main()
