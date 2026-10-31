# Cadence bench: the servo frame chosen live by the law (31 Oct 2026)

Script: `cadence_bench.py` (sha256 e2c687827b2ff4f708266afc14564d737dd2780ae597debaa66d03a995521181), pushed to GitHub
with its predictions before any scored run (commit 66892fc). Results: `results/summary.json`. Written and run by the
assistant in a cloud container (8.8 minutes); not yet reproduced by the author.

## Question

In the MAESTRO 1 videos the shoulders' window sits on its floor during continuous motion. On the floor the law decides
nothing: the engine integrates over a fixed window set by the hardware. The law itself says why: when its window falls
below the floor D, the loop is too slow for the motion. Can the machine follow that instruction live, taking at each
servo frame the slowest frame period whose floor stays under the law's window, and so stay in the range where the law
works?

## Setup

The self-contained joint bench of `examples/servo_bench.py` (same servo model, excerpts A and B of Beethoven's Fifth,
encoder model, baton-tip kinematics), baton arm (RJ1 to RJ6), loop and world at 2 kHz, the bench engine
(identification 'raw') with its own floor rule unchanged: floor(1, 2, 3, 5, 8 ms) = 1.11, 1.35, 1.85, 2.85, 4.35 ms.
Arms: fixed frames of 3 ms (today's hardware), 2 ms and 1 ms; AD, adaptive over {1, 2, 3} ms (may only speed up from
3 ms); ADS, adaptive over {1, 2, 3, 5, 8} ms (may also slow down). No margin, no hysteresis. Steady encoder noise and a
noise step (variance x10 from 10 to 20 s), five paired seeds each: 20 runs per arm, every arm on the same reads.
Frames faster than 3 ms are a simulation probe: the bench's servos accept about 330 Hz.

## Results

| Arm | Pooled tip error (mm) | Against 3 ms | Mean frame rate (Hz) | Tip jitter (mm) | Joint error (deg) |
|---|---|---|---|---|---|
| 3 ms fixed | 1.129 | 1 | 333 | 0.782 | 0.0510 |
| 2 ms fixed | 0.979 | 0.867 | 500 | 0.683 | 0.0436 |
| 1 ms fixed | 0.839 | 0.743 | 1000 | 0.580 | 0.0378 |
| AD (1 to 3 ms) | 0.898 | 0.796 | 541 | 0.655 | 0.0403 |
| ADS (1 to 8 ms) | 0.953 | 0.844 | 462 | 0.723 | 0.0425 |

By cell (seed-pooled tip error, mm; A steady, A step, B steady, B step): 3 ms 0.716, 0.869, 1.306, 1.457; 1 ms 0.516,
0.683, 0.945, 1.090; AD 0.566, 0.752, 0.972, 1.182; ADS 0.605, 0.856, 0.984, 1.252. AD's mean frame rate by cell: 485,
447, 652, 578 Hz. Paired, run by run: AD against 3 ms 0.74 to 0.87; AD against 1 ms 1.03 to 1.11; ADS against AD 1.01
to 1.16.

Time on the floor (the law's window below the floor of the current frame, baton-arm joint-time): 3 ms fixed 36 %, 26 %,
66 %, 50 % in the four cells; 1 ms fixed and AD 12 %, 9 %, 30 %, 22 %. AD spends 27 % of the time at 1 ms, 18 % at
2 ms, 56 % at 3 ms, and changes frame 0.4 times per joint per second. ADS spends 28 % of the time at 8 ms and changes
2.3 times per second.

## Predictions

- P0, gate (reproduce `examples/servo_bench.py` at 1 kHz, 3 ms, 12 joints, to 1e-9 mm): held, difference exactly
  zero in all four cells.
- P1, the wall (3 ms: 15 to 70 % of the time on the floor under steady noise, less under the step): held (36 % and
  66 %; 26 % and 50 % under the step).
- P2, a faster frame helps (1 ms at or below 0.95 of 3 ms): held, 0.743.
- P3, main (AD within 2 % of the best fixed arm, at most 75 % of its frame rate): **failed**. The frame-rate half
  holds (541 Hz, 54 % of the 1 ms arm) but AD's tip error is 7.1 % above the 1 ms arm's, and above it in every run
  (by 2.8 to 11 %).
- P4, slowing down costs (ADS more than 5 % above AD): held, 6.1 %.
- P5, speeding up only is safe (AD at or below 3 ms, within 1 % in every cell): held; AD is 20 % below 3 ms
  pooled and below it in every run.

## What it means

- Following the law's instruction pays: frames chosen live cut the tip error by 20 % against today's 3 ms frame, at
  a mean rate of 541 Hz instead of 1000. That is 80 % of the gain of a fixed 1 ms frame for 54 % of its frames. But
  the pre-registered claim, parity with the fastest frame, fails.
- The likely reason, a reading after the verdict and not tested: the servo frame does two things, it bounds the
  window (the floor) and it holds the plan between commands. The rule covers only the first. P4, which held, points
  the same way: slowing down when the window allows it still costs, because the plan is held longer. A rule that
  also accounts for the plan's motion within a frame would have to be derived and pre-registered before a new run.
- Even at 1 ms the law's window is under its floor 9 to 30 % of the time: on this arm, with these encoders, part of
  the motion asks for a faster loop than any frame tested.
- Frames under 3 ms need a faster actuator link than the bench's servos (about 330 Hz). On today's hardware the
  rule could only choose 3 ms or slower; that arm was not run, and P4 suggests it would trade accuracy for fewer
  frames.

## Run it

From the repository root (Python 3.9 or later, numpy, scipy): `python bench/cadence_2026-10-31/cadence_bench.py`
(about 10 minutes on two cores; it writes `results/summary.json`).
