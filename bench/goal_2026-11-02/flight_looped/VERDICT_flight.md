# Flight confirmation, 2 Nov 2026: syntonic attitude estimators vs tuned Kalman filters

## Pre-registration and gates

- **Script:** confirm_flight.py, sha256 6eac13e14a5a0755bd1a6b773baf5bcf8888a4d95e0aa316819a1af2dcbcbaea.
- **Frozen with it:**
  - law_flight.py (e0e9fbaa…);
  - the snapshot of the goal/flight harness (common.py, variants.py, variants2.py, variants3.py);
  - facts2.json, rival.json and sched_tuned.json;
  - the repository bench files and the two logged records (14 hashes).
- **Gates** (G0 hashes, G1 dev values reproduced exactly, G2 commissioning edges) all pass. The one declared exception is KF_sched_btab's rho_half on a flat surface.
- **Seeds:** 400-409, never used before. 60 flights.

## Setup

- **Bench:** the repository's flight bench, amendment-2 noise. The records are the logged X500 S6 idle/low and S5 half-throttle records, looped. It models one axis with the plant at rest, a 5 s idle pre-roll, and the arming calibration as the gyro mean.
- **Circumstances:**
  - C1: steady low throttle.
  - C2: steady half throttle.
  - C3: half/low every 15 s.
  - C4: vibration step, tilt noise ×√10 from 20 to 40 s.
  - C5: half/low every 5 s.
  - C6: gyro bias step +0.02°/s at 30 s.
- **Rivals:** all are tuned against the truth, which is an oracle no real filter has. Tuning used training seeds 100-102 on the C3 schedule.
  - KF_rival: the corpus bench's two-state Kalman, fixed R and Q.
  - KF_sched: the same filter with R scheduled on the throttle level. Tuned: m = 2048, rho_half = 0.354, Q_G = Q_B = 0.
  - KF_rival_btab and KF_sched_btab: the same two filters given the per-throttle bias table.
- **Laws:** nothing is tuned against truth.
  - LAW_I = b4_nm_ic2: two-state law with the live noise reading and per-throttle memory.
  - LAW_II = the level law: Riccati recursion of (θ, b_idle, b_low, b_half) fed by facts and the bias table.

## Results (pooled rms estimation error, degrees)

| | C1 | C2 | C3 | C4 | C5 | C6 |
|---|---|---|---|---|---|---|
| KF_rival | 0.0699 | 0.2902 | 0.2794 | 0.0756 | 0.3401 | 0.2083 |
| KF_sched | 0.0609 | 0.2066 | 0.0994 | 0.1000 | 0.1029 | 0.1280 |
| KF_rival_btab | 0.0654 | 0.2537 | 0.2695 | 0.0694 | 0.3425 | 0.1885 |
| KF_sched_btab | 0.0649 | 0.1424 | 0.1007 | 0.0697 | 0.0846 | 0.1809 |
| LAW_I | 0.0609 | 0.2560 | 0.1959 | 0.0868 | 0.2169 | 0.0996 |
| LAW_II | 0.0696 | 0.1915 | 0.1370 | 0.0871 | 0.1259 | 0.1149 |

### Ratios

| | C1 | C2 | C3 | C4 | C5 | C6 |
|---|---|---|---|---|---|---|
| LAW_I / KF_rival | 0.872 | 0.882 | 0.701 | **1.149** | 0.638 | 0.478 |
| LAW_I / KF_sched | 1.001 | 1.239 | **1.972** | 0.868 | **2.108** | 0.778 |
| LAW_II / KF_rival_btab | 1.065 | 0.755 | 0.509 | **1.256** | 0.368 | 0.609 |
| LAW_II / KF_sched_btab | 1.073 | 1.344 | 1.361 | 1.251 | 1.488 | 0.635 |
| KF_sched / KF_rival | 0.872 | 0.712 | 0.356 | 1.323 | 0.303 | 0.614 |

Law wins out of 10 seeds:

- LAW_I vs KF_rival: C1 7, C2 8, C3 10, C4 2, C5 10, C6 10.
- LAW_I vs KF_sched: C1 4, C2 1, C3 1, C4 9, C5 1, C6 10.

## Verdicts

- **F1 FAILS.** Parity with the corpus's tuned Kalman in C1-C6. C4 is 1.149, above 1.10. The law is better in the other five circumstances, by 13 to 52 %.
- **F2 PASSES.** Better than KF_rival in C1, C5 and C6.
- **F3 PASSES.** Parity with the throttle-scheduled Kalman in C1, C4 and C6.
- **F4 PASSES, as predicted.** The law does not reach the scheduled Kalman in the switching flights: 1.97 in C3 and 2.11 in C5. **The goal is not met against this rival.**
- **F5 FAILS.** With the bias table given to both, C4 is 1.256.
- **F6 PASSES.** Better than KF_rival_btab in C3, C5 and C6.
- **F7 PASSES, as predicted.** The law is not at parity with KF_sched_btab in C2, C3 and C5. **The goal is not met there.**
- **F8 FAILS.** C4 is 1.251.

## What it means

**The flight half of the goal is not met.**

- **Against the corpus's tuned Kalman:** the law wins five circumstances and loses the vibration step. In development the vibration step was 1.03, with a bootstrap interval up to 1.17. That dev figure was the best of many variants and was optimistic.
- **Against a Kalman that also schedules R on the throttle level and is tuned on truth:** the law loses by 2× in switching flights.

**Common cause.** The steady-state law has no memory of the information it has accumulated. The scheduled Kalman carries the precision gathered at low throttle (or before the vibration) through the noisy phases.

**What I tried.** I tried three fact-based Riccati forms: fact R per throttle, a switch-aware bias variance, and a four-state per-level bias memory. None closed the gap on dev without breaking another circumstance.

**Bench caveats.**

- The looped records reward long memories: the truth-tuned filters choose Q_G = 0, and KF_sched also chooses Q_B = 0.
- KF_sched_btab, tuned with the in-sample exact bias table, ignores the half-throttle tilt altogether. A real filter would drift doing that.
- These caveats favour the truth-tuned filters, but the verdicts stand as pre-registered.
