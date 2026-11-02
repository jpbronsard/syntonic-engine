# Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
# Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
"""Unit tests of the scan of the C++ and Rust loop law with the two speed-ups of 0.3.2 (loop_scan.py).
   python -m unittest discover -s python/tests"""
import os, sys, functools, unittest
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from syntonic_engine import SyntonicLoopLaw, ScanLoopLaw
from syntonic_engine.loop_law import LoopLaw, loop_delay
from syntonic_engine.loop_scan import LoopScan
from test_loop_law import closed_loop


class Scan(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.law, cls.err = closed_loop(SyntonicLoopLaw)

    def same(self, **kw):
        law, err = closed_loop(functools.partial(ScanLoopLaw, **kw))
        np.testing.assert_array_equal(err, self.err)                              # the same commands, bit for bit
        np.testing.assert_array_equal(law.Kp, self.law.Kp); np.testing.assert_array_equal(law.Ki, self.law.Ki)
        self.assertGreater(len(law.stats), 3)
        return law

    def test_same_gains_in_every_mode(self):
        for prune in (False, True):
            for quad in (False, True):
                law = self.same(prune=prune, quad=quad)
                self.assertEqual(law.fallbacks, 0)
                share = np.mean([s['points'] / s['total'] for s in law.stats])
                if prune:
                    self.assertLess(share, 0.9)                                     # 0.66 here
                else:
                    self.assertEqual(share, 1.0)

    def test_fallback_is_the_frozen_computation(self):
        law = self.same(delta=1e3)                                                  # no choice is certified
        self.assertEqual(law.fallbacks, len(law.stats))

    def test_delta_is_checked(self):
        with self.assertRaises(ValueError):
            ScanLoopLaw(np.zeros(1), 6.0, delta=1e-12)

    def test_quadratic_form_is_rounding(self):
        D = loop_delay(1e-3, 3e-3) + 0.0025; sc = LoopScan(0.025, D, LoopLaw(0.025, D))
        worst = max(float(np.max(np.abs(sc.s2(i, k, True) / sc.s2(i, k, False) - 1))) for i, k in sc.rows)
        self.assertLess(worst, 1e-13)                                               # 1.3e-15 here


if __name__ == '__main__':
    unittest.main()
