# Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
# Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
"""Syntonic engine: one formula, three facts about the machine, no tuning. The time window of a correction or an
estimate from Theorem U' of the Syntony Principle."""
from .engine import SyntonicEngine, kappa_H, C_H, loop_delay, window_floor

__version__ = '0.2.2'
__all__ = ['SyntonicEngine', 'kappa_H', 'C_H', 'loop_delay', 'window_floor']
