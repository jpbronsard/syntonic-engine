# Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
# Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
"""Syntonic engine: one formula, three facts about the machine, no tuning. The time window of a correction or an
estimate from Theorem U' of the Syntony Principle. Since 0.3.0 also two laws whose gains are computed from measured
facts: SyntonicLoopLaw (a servo joint's PI gains) and RiccatiAttitudeLaw (a one-axis attitude estimator); since 0.3.1
also in C++ (cpp/) and Rust (rust/), checked against these classes. Since 0.3.2 ScanLoopLaw: the loop law with its gains
found as the C++ and Rust laws find them, with two speed-ups that leave them unchanged (the reference of those ports)."""
from .engine import SyntonicEngine, kappa_H, C_H, loop_delay, window_floor
from .loop_law import SyntonicLoopLaw
from .attitude_law import RiccatiAttitudeLaw
from .loop_scan import ScanLoopLaw

__version__ = '0.3.1'
__all__ = ['SyntonicEngine', 'kappa_H', 'C_H', 'loop_delay', 'window_floor', 'SyntonicLoopLaw', 'RiccatiAttitudeLaw',
           'ScanLoopLaw']
