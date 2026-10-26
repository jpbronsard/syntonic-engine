# Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
# Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
"""Syntonic engine: joint correction with a time window from Theorem U' of the Syntony Principle (nothing tuned)."""
from .engine import SyntonicEngine, kappa_H, C_H, loop_delay, window_floor

__version__ = '0.2.0'
__all__ = ['SyntonicEngine', 'kappa_H', 'C_H', 'loop_delay', 'window_floor']
