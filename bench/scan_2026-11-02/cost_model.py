# Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
# Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
"""The gain recomputation's cost on a Teensy 4.1 (600 MHz) for two joints, before and after the speed-ups of 0.3.2, and
the spread each per-tick budget needs (a model, not a measurement; examples/teensy/LoopLawTiming.ino measures).

Cycles per point from m7/cycles.json (LLVM's Cortex-M7 model on arm-none-eabi-g++'s code): 0.3.1's candidate loop, the
quadratic form's, and the loop that fills a Kp row's P. Assumed: a spectrum point (newlib's pow) 500 cycles. Shares of
the points evaluated with pruning from scan_check.json (the realistic bench's two-joint log, RJ5 and RJ6).
  worst case  no point pruned: what a per-tick budget must allow for, since pruning cannot be counted on
  typical     the mean and the largest share evaluated on the two-joint log
The spread for a budget is the law's own rule (spread_for in cpp/SyntonicLoopLaw.h and rust/src/loop_law.rs): the
smallest m whose share of the work, counted as if nothing were pruned, plus one work item, fits the budget, with the
time of a candidate from the model. Writes cost.json.   python bench/scan_2026-11-02/cost_model.py
"""
import os, sys, json, math
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'python'))
from syntonic_engine.loop_law import LoopLaw, loop_delay  # noqa: E402

HZ, NW, N, POW = 600e6, 1500, 2, 500.0
CYC = json.load(open(os.path.join(HERE, 'm7', 'cycles.json')))
SHARE = {r['log']: r['prune+quad'] for r in json.load(open(os.path.join(HERE, 'scan_check.json')))}['replay_loop_bench2.bin']
BUDGETS = (0.00045, 0.00025)
EVERY = 500


def spread_for(U, unit_s, item, budget):
    """the law's rule: 0 if all fits in one tick, else the smallest m < every with (ceil(U / m) + item) units <= budget"""
    if U * unit_s <= budget:
        return 0
    for m in range(1, EVERY):
        if (math.ceil(U / m) + item) * unit_s <= budget:
            return m
    return None


def case(dead):
    D = loop_delay(1e-3, 3e-3) + dead; LL = LoopLaw(0.025, D)
    nok = int(LL.ok.sum()); nrows = int(LL.ok.reshape(41, 60).any(1).sum())
    spec = NW * N * POW
    before = nok * NW * CYC['before'] + spec
    worst = (nok + N) * NW * CYC['quad'] + (nrows + N) * NW * CYC['row'] + spec
    typ = {k: SHARE[k] * nok * NW * CYC['quad'] + (nrows + N) * NW * CYC['row'] + spec for k in ('share_mean', 'share_max')}
    # the law's units: 0.3.1 (spectrum 50 frequencies x 8 a joint, candidates) and 0.3.2 (rows of P, warm candidates)
    U_old, item_old = NW * 8 * N + nok * NW, 50 * 8
    U_new, item_new = NW * 8 * N + nok * NW + N * NW + (nrows + N) * 900, 80
    out = dict(dead_ms=dead * 1e3, robust=nok, rows=nrows, before_s=before / HZ, worst_s=worst / HZ,
               typical_mean_s=typ['share_mean'] / HZ, typical_max_s=typ['share_max'] / HZ, spreads={})
    for B in BUDGETS:
        m_old = spread_for(U_old, NW * CYC['before'] / HZ / NW, item_old, B)
        m_new = spread_for(U_new, NW * CYC['quad'] / HZ / NW, item_new, B)
        tick_new = None if not m_new else (math.ceil(U_new / m_new) + item_new) * CYC['quad'] / HZ
        # not guaranteed: the same rule with the work scaled by the largest share seen (pruning counted on)
        m_cred = spread_for(SHARE['share_max'] * U_new, NW * CYC['quad'] / HZ / NW, item_new, B)
        out['spreads'][f'{B * 1e3:.2f} ms'] = dict(before=m_old, after=m_new, after_tick_bound_ms=None if tick_new is None
                                                   else tick_new * 1e3, after_if_pruning_counted=m_cred)
    return out


if __name__ == '__main__':
    res = dict(cycles=CYC, assumed_pow_cycles=POW, shares=SHARE, cases=[case(0.0025), case(0.0)])
    for c in res['cases']:
        print(f"dead {c['dead_ms']:.1f} ms, {c['robust']} robust candidates in {c['rows']} Kp rows, two joints:")
        print(f"  recomputation: 0.3.1 {c['before_s']:.3f} s; 0.3.2 worst case {c['worst_s']:.3f} s, "
              f"typical {c['typical_mean_s']:.3f} s (largest {c['typical_max_s']:.3f} s)")
        for b, s in c['spreads'].items():
            print(f"  budget {b} per tick: spread 0.3.1 {s['before']}, 0.3.2 {s['after']}"
                  + (f" (at most {s['after_tick_bound_ms']:.3f} ms a tick)" if s['after'] else ' (none below every)')
                  + f"; if pruning were counted on at its largest share, {s['after_if_pruning_counted']}")
    json.dump(res, open(os.path.join(HERE, 'cost.json'), 'w'), indent=1)
