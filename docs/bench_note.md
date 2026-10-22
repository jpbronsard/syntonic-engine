# MAESTRO joint bench: loop rate, window floor, and parity with PID (runs v6 to v10, 21 Oct 2026)

This note records a day of tests on the MAESTRO 1 joint bench. Each run had its predictions written down before it started, and the note keeps the ones that failed next to the ones that held.

The scripts are in `bench/scripts/`: `flow_v6_rate.py` runs the loop at 2 kHz (world and controller every 0.5 ms), `flow_v7_ctrl.py` runs the controller at 500 or 333 Hz with the world at 1 kHz, `flow_v8_floor.py` and `flow_v9_frame.py` test the window floor across loop rates and servo frames, and `flow_v10_seeds.py` repeats the comparison with PID over ten noise draws. Every result file is in `bench/results/`.

## In short

- At 1 kHz, the syntonic engine (Theorem U', nothing tuned) matches a PID tuned on the same movement at the baton tip when the noise is steady. The tip-error ratio is 1.006 on one excerpt and 1.063 on the other, both inside the ±10 % margin declared beforehand. When the noise jumps, the syntonic engine does better at the tip, by 29 % and 11 %. These figures come from ten paired noise draws per case (v10). Parity is what we expected: a PID tuned for fixed conditions sits close to the optimum for those conditions, and the syntonic engine reaches that optimum by derivation instead.
- Scored on joint angles, each joint on its own, the syntonic engine's error is about 28 % larger than the tuned PID's under steady noise. The PID's gains were chosen to minimise exactly that error. The syntonic engine's joint errors largely cancel at the tip (see "Tip against joints" below).
- The window needs a lower bound equal to the loop delay, D = half the servo PWM frame + half a controller tick + 0.1 ms for the AS5048. We wrote this rule down before run v8, and run v9 confirmed it across servo frames of 2, 3, 5 and 8 ms, within 0.4 % of the best floor found by sweeping. Because it comes from the hardware, the loop stays free of tuned parameters.
- A floor below the loop delay costs accuracy, and with slow servo frames it costs stability too. With an 8 ms frame, a 2 ms floor makes the loop ring and the tip error climbs to 13.5 mm, against 1.96 mm with the right floor.
- For MAESTRO (333 Hz frame, 1 kHz loop), D comes to 2.1 ms. The v4 video engine runs with a 2.22 ms floor, the TH41-A limit at 1 kHz, which sits just above D.
- Still open: the rule missed by 7 % with the loop at 333 Hz. The likely reason is that in the simulation the servo frame lines up exactly with the loop tick at that rate, so the correction is always fresh and the half-tick term overcounts. On the Teensy the PWM is not synchronised with the loop. This can be checked by shifting the frame off the tick in the simulation. The remaining doubt about parity comes from having only two musical excerpts; the noise draws barely move the result.

## Setup, fixed before each run

- The bench is the one from v5c (19 Oct). Two excerpts, A at 0:06 to 0:36 and B at 3:37 to 4:07, each run once with steady noise and once with the noise variance raised tenfold from 10 to 20 s. The scores are the rms baton-tip error and its jitter above 5 Hz, in mm.
- In v6 and v7 only the loop rate changes. The servo PWM frame stays at 333 Hz (every 3 ms) and carries the plan plus the latest correction. The noise per encoder read does not change, since the AS5048A refreshes at 10.2 to 12.4 kHz, and every time constant is kept in seconds. The noise term already carries dt, as the corpus states it: sigma^2 is the per-read variance times dt (in the code, s2 = 0.5 n2 dt).
- The PID with interpolation and the best fixed window are re-tuned at each rate on the v5 grids. At 500 and 333 Hz the derivative gain hit the top of its grid, so both were re-tuned on a wider grid; the PID improved by at most 2 %.
- The syntonic engine runs the v5 code unchanged, with nothing tuned.
- Each script first checks that it reproduces v5c with the original settings (1 kHz loop, 3 ms frame, 2 ms floor) before changing anything.

## v10: parity with PID over ten paired noise draws

Loop at 1 kHz, 3 ms frame, syntonic floor 2.22 ms, all declared before the run. Each case gets ten fresh noise draws, and both engines see exactly the same noise. The ratio rho is the syntonic tip rms divided by the PID's, reported as a geometric mean with its 95 % interval (t distribution, n = 10). "On par" was defined in advance as an interval inside [0.90, 1.10].

| Case | PID (mm) | Syntonic (mm) | rho [95 %] | Jitter ratio [95 %] | Syntonic better |
|---|---|---|---|---|---|
| A steady | 0.727 | 0.732 | 1.006 [1.000, 1.013] | 0.858 [0.850, 0.867] | 2/10 |
| A noise x10 | 1.327 | 0.942 | 0.709 [0.697, 0.722] | 0.563 [0.555, 0.571] | 10/10 |
| B steady | 1.194 | 1.270 | 1.063 [1.060, 1.066] | 1.074 [1.066, 1.083] | 0/10 |
| B noise x10 | 1.636 | 1.455 | 0.889 [0.882, 0.897] | 0.744 [0.731, 0.758] | 10/10 |

- E1, A steady inside the margin with an interval containing 1: held, at the edge.
- E2, B steady with the PID slightly ahead (1.05 to 1.08) and inside the margin: held at 1.063.
- E3, syntonic better under the noise step on both excerpts: held.
- From one draw to the next the result moves by about ±1 %, so the differences are systematic. The 6 % on B is real. Whether "on par" holds more widely depends on the excerpts, not on the noise.

## Tip against joints

This came up while building the open example. Every run above is scored on the baton tip, as fixed before the first run. Scoring each joint's angle error instead, averaged over RJ1 to RJ6 (single draws, syntonic floor 2.22 ms, `examples/servo_bench.py`), gives a different picture:

| Case | Tip error PID / syntonic (mm) | Joint error PID / syntonic (deg) |
|---|---|---|
| A steady | 0.732 / 0.734 (x1.00) | 0.0309 / 0.0395 (x1.28) |
| A noise x10 | 1.330 / 0.922 (x0.69) | 0.0538 / 0.0475 (x0.88) |
| B steady | 1.193 / 1.275 (x1.07) | 0.0489 / 0.0627 (x1.28) |
| B noise x10 | 1.631 / 1.450 (x0.89) | 0.0674 / 0.0701 (x1.04) |

- Joint by joint, the syntonic engine's error is about 28 % larger than the tuned PID's under steady noise. That is the error the PID's gains were chosen to minimise.
- Every joint moves the tip. The median lever, in mm per degree, is 5.2 for RJ1, 13.3 for RJ2, 13.0 for RJ3, 8.9 for RJ4, 12.0 for RJ5 and 8.2 for RJ6. Weighting the joint errors by these levers, as if they were independent, predicts a syntonic tip error about 25 % above the PID's. The measured tip error is only 0 to 7 % above. So the syntonic engine's joint errors largely cancel at the tip, while the PID's do not: on A, the independent estimate for the PID is 0.759 mm and the measured value 0.732 mm. We do not yet know why the syntonic errors cancel.
- What can be said in public: on par with the tuned PID at the baton tip under steady noise, better at the tip when the noise jumps, and on joint angles the tuned PID is tighter under steady noise.

## Loop rate (v6 and v7): tip error rms, mm

PID against the syntonic engine with a 2 ms floor (each cell reads PID / syntonic):

| Loop | A steady | A noise x10 | B steady | B noise x10 |
|---|---|---|---|---|
| 2 kHz | 0.67 / 0.71 | 1.13 / 0.86 | 1.21 / 1.30 | 1.52 / 1.46 |
| 1 kHz | 0.73 / 0.73 | 1.33 / 0.92 | 1.19 / 1.28 | 1.63 / 1.45 |
| 500 Hz | 0.84 / 0.86 | 1.58 / 1.15 | 1.27 / 1.38 | 1.79 / 1.61 |
| 333 Hz | 0.85 / 0.91 | 1.59 / 1.29 | 1.31 / 1.43 | 1.84 / 1.76 |

With the floor set to twice the tick instead, the syntonic engine does worse than with the 2 ms floor at every rate except 1 kHz, where the two are the same: at 2 kHz (1 ms floor) 0.91, 1.04, 1.84, 1.76; at 500 Hz (4 ms) 0.94, 1.20, 1.49, 1.68; at 333 Hz (6 ms) 1.18, 1.47, 1.79, 1.99.

Predictions for v6:
- Q1, every engine improves a little: partly failed.
- Q2, the syntonic engine gains most on B steady once the floor drops: failed (1.28 went to 1.84 mm).
- Q3, the window shrinks by about 0.79 away from its bounds: held.

Predictions for v7:
- S1, the PID gets a little worse: held.
- S2, the syntonic engine with the 4 ms floor loses more than the PID on steady noise: held.
- S3, the syntonic engine with the 2 ms floor stays close to its 1 kHz result: half held. The window grew by 1.22 to 1.29 as predicted, but the error rose by 8 to 25 %, about as much as the PID's.

## v8: floor set by the loop delay, across loop rates

Declared before the run: D = 1.85 ms at 2 kHz, 2.10 ms at 1 kHz, 2.60 ms at 500 Hz and 3.10 ms at 333 Hz, compared with a sweep of floors at 1, 1.5, 2, 2.5, 3 and 4 ms.

| Loop | D | Syntonic with floor D (A, A x10, B, B x10) | Best floor in the sweep | Largest gap of D |
|---|---|---|---|---|
| 2 kHz | 1.85 ms | 0.72, 0.86, 1.31, 1.47 | 2.5 ms | 3.1 % |
| 1 kHz | 2.10 ms | 0.73, 0.92, 1.28, 1.45 | 2 ms | 0.2 % |
| 500 Hz | 2.60 ms | 0.87, 1.16, 1.39, 1.62 | 1 ms (flat from 1 to 2 ms) | 1.6 % |
| 333 Hz | 3.10 ms | 0.95, 1.31, 1.48, 1.78 | 1 ms | 7.0 % |

- F1, D within 5 % of the best everywhere: held at 2 kHz, 1 kHz and 500 Hz, failed at 333 Hz.
- F2, the best floor rises as the rate falls: failed. It falls instead.
- F3, the 1 kHz run reproduces v5c: held.

After seeing these results, a second rule (never correct more than the full error within one servo frame) fitted all four rates within 1 %. Run v9 was designed to tell it apart from D and from a fixed 2 ms floor.

## v9: floor across servo frames

Loop and plant at 1 kHz, declared before the run. The servo frame is set to 2, 3, 5 and 8 ms. The 2 ms frame is faster than the BLS3355 allows (330 Hz), so it only serves as a simulation probe. The sweep covers floors of 1, 1.5, 2, 3, 4, 5, 6, 8 and 10 ms, scored by the root mean square of the tip rms over the four conditions.

| Frame | Best floor in the sweep | Delay rule D (gap) | Per-frame rule (gap) | Fixed 2 ms (gap) |
|---|---|---|---|---|
| 2 ms | 1 ms (1.005) | 1.60 ms (+0.3 %) | 1.44 ms (+0.2 %) | +0.7 % |
| 3 ms | 2 ms (1.131) | 2.10 ms (+0.1 %) | 2.47 ms (+0.3 %) | 0 % |
| 5 ms | 3 ms (1.436) | 3.10 ms (0.0 %) | 4.48 ms (+3.9 %) | +7.5 % |
| 8 ms | 5 ms (1.962) | 4.60 ms (-0.4 %) | 7.49 ms (+8.8 %) | +589 % (13.5 mm, the loop rings) |

- G1, the best floor rises with the frame: held (1, 2, 3 and 5 ms).
- G2, the per-frame rule stays within 2 % at every frame: failed. It overshoots, by 3.9 % at 5 ms and 8.8 % at 8 ms.
- G3, a fixed 2 ms floor falls more than 2 % behind at 5 and 8 ms: held. It is 7.5 % behind at 5 ms and unstable at 8 ms.
- The delay rule D, written down before v8, lands within 0.4 % of the best floor at every frame and within 1.4 % in every single condition.

## What it means

- The PID gets steadily better as the loop gets faster. The syntonic engine does best at 1 kHz (or at 2 kHz with the right floor) and gets worse below 1 kHz.
- Whenever the noise changes, the syntonic engine keeps its advantage at the tip, at every loop rate, and the advantage is largest at 1 kHz. With steady noise it matches the tuned PID at the tip.
- The syntonic engine follows the noise per unit time, just as its law says, and that gives the best estimate. But a closed loop adds a delay that the estimator does not model, since U'3 counts only tracking and filtering. A window shorter than that delay loses accuracy and, with slow servo frames, stability. This is a limit of the control loop, derived from the hardware and declared alongside U'. It does not refute U'.

## For the build

- Keep the loop at 1 kHz and the servo frame at 333 Hz. Use a floor of D = 2.1 ms (the half-tick term applies, because the Teensy's PWM is asynchronous), or the 2.22 ms TH41-A limit used in the v4 video engine. On the real arm, measure the full loop delay once, since the servo's internal latency is not in the model, and use the larger of the two values.
- If the servo frame or the loop rate changes, recompute D, and never let the window go below it.
- A 1 kHz loop fits every bus as it is: reading the shoulder encoders over I2C at 400 kHz takes about 36 % of the tick, and sending the angles over classic CAN about 25 % of the bus.

## Open code

The engine used in these runs is published with this note under the Apache 2.0 licence, in three languages that give the same numbers:

- Python (`python/syntonic_engine`), extracted from the bench engine and checked to be identical to it: zero difference in joint angles and windows over 33 s on all twelve joints.
- C++ (`cpp/SyntonicEngine.h`), a single header for Arduino and Teensy, which agrees with the Python engine to 1.5e-14 rad on a 33 s, 12-joint replay (`cpp/tests/replay.cpp`).
- Rust (`rust/`), a `no_std` crate without allocation, which agrees with the Python engine to 5.4e-14 rad on the same replay.

`examples/servo_bench.py` reproduces the tip and joint scores of the "Tip against joints" table without the MAESTRO 1 simulation model. Neither the C++ header nor the Rust crate has run on a microcontroller yet.
