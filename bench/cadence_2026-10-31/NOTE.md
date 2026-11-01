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

## Diagnostic: window against command hold (exploratory, after the verdict)

Script `diag_hold_window.py` (sha256 1b82a0acc90422c2591edb7763fe0ab916e4b0eedd9e3d6ad67a615dfb26f1fd), pushed before
its run (commit 66626da); results `results/diag_hold_window.json`. No prediction judged. A servo frame bounds the window
(its floor) and holds the command, plan and correction, between two sends. Each column here gets its own floor, send
period and correction refresh, on the same 20 paired runs. Gate: the 3 ms and 1 ms columns reproduce the cadence arms
run by run, difference exactly zero.

| Column | Command sent every | Floor | Pooled tip error (mm) | Cells (A steady, A step, B steady, B step) |
|---|---|---|---|---|
| A | 3 ms | 1.85 ms | 1.129 | 0.716, 0.869, 1.306, 1.457 |
| D | 1 ms | 1.85 ms | 0.851 | 0.538, 0.700, 0.956, 1.097 |
| E | 1 ms (correction refreshed every 3 ms) | 1.85 ms | 0.988 | 0.667, 0.828, 1.088, 1.262 |
| B | 1 ms | 1.11 ms | 0.839 | 0.516, 0.683, 0.945, 1.090 |
| C | 3 ms | 1.11 ms (a probe, not realizable) | 1.309 | 0.899, 1.020, 1.524, 1.637 |

- Of the 0.290 mm gained from 3 ms to 1 ms, sending the command every 1 ms gives 0.278 mm (96 %, A to D), about half
  from the plan (A to E, 0.140 mm) and half from the correction (E to D, 0.138 mm). Lowering the floor then gives
  0.012 mm (4 %, D to B), though D's window sits on its floor 44 % of the time. Run by run: D against A 0.73 to 0.82;
  B against D 0.95 to 1.00.
- With the command held 3 ms, a floor under the loop delay makes things worse: C is 16 % above A (1.09 to 1.30 run by
  run), as the delay rule of runs v8 and v9 says.
- On this arm, the window on its floor costs little (about 1.4 %); what limits is the time the command is held. That
  is why P3 failed: when the law's window allows 3 ms, the plan and the correction are still held 3 ms.

## Timing control: does the law know when to speed up? (pre-registered T1)

Script `timing_control.py` (sha256 d17bf7ec729e06f92fd0287b096170eebd12a0a44851b9f58790c0d666770157), pushed before
its run with its prediction (commit a0cbc91), design agreed with the author; results `results/timing_control.json`.
SH replays AD's own schedule from the same run (sends and frames, hence floors) with the scored 30 s shifted by 15 s:
same frames, same shares, same changes, no longer lined up with the motion. JF gives each joint one fixed frame, the
period on the 0.5 ms grid closest to its mean period in AD. Gates: the cadence arms reproduce `summary.json` and AD
replayed unshifted reproduces AD, run by run, both with a difference of exactly zero.

| Arm | Pooled tip error (mm) | Mean frame rate (Hz) | Cells (A steady, A step, B steady, B step) | Time on the floor |
|---|---|---|---|---|
| AD, the law's timing | 0.898 | 541 | 0.566, 0.752, 0.972, 1.182 | 18 % |
| JF, one frame per joint | 0.923 | 546 | 0.610, 0.769, 1.008, 1.196 | 24 % |
| SH, AD's frames shifted 15 s | 1.021 | 541 | 0.630, 0.822, 1.135, 1.344 | 35 % |

- T1 (AD at least 5 % below SH, and below it in at least 18 of 20 runs): **held**. AD is 12 % below SH (0.880),
  below it in 20 of 20 runs (0.85 to 0.92 run by run). With the same frames, the law's timing is worth 12 %: it
  puts its fast frames where the motion needs them.
- Described: AD is 2.7 % below JF and below it in all 20 runs (0.92 to 0.99); JF is 10 % below SH. Against a fixed
  frame per joint, chosen from the law's own average, the timing adds less: most of what SH loses is fast frames at
  the wrong moments. JF is not law-free: its periods come from AD (in excerpt A mostly 2, 3, 2.5, 1.5, 2.5 and 2 ms for
  RJ1 to RJ6).
- Reading: the law's window works as a measure of how much cadence the motion needs at each instant. As a cadence
  signal it is right in time (T1); the size of what it adds depends on the reference (12 % against the same frames
  mistimed, 2.7 % against the best fixed frame per joint, 7 % against the fixed-frame frontier at the same budget).

## Cadence along the e-geodesic (pre-registered E1 to E3)

Script `egeodesic_cadence.py` (sha256 d58a4ae9838dd914fa87e87319381419063034bf081568b9b12ca01b07127d49), pushed before
its run with its predictions (commit 152f6ee), idea and go from the author; results `results/egeodesic_cadence.json`.
EG: the frame follows a copy of the law's window carried on the e-geodesic of Dual Geodesy (1/tau stepped toward the
law's raw window at the engine's own bounded speed, without the floor), mapped by the engine's floor rule onto the
loop's finest grid (1, 1.5, 2, 2.5, 3 ms). GR: the same grid from the raw law. AD: the staircase {1, 2, 3} ms. EGsh:
EG's own schedule shifted 15 s. Gates (F3, F1 and AD reproduce `summary.json`; EG replayed unshifted reproduces EG):
difference exactly zero.

| Arm | Pooled tip error (mm) | Mean frame rate (Hz) | Gain against the fixed frames at that rate | Switches per joint per second | Mean switch |
|---|---|---|---|---|---|
| AD, staircase | 0.898 | 541 | 7.1 % | 0.44 | 1.0 ms |
| GR, fine grid, raw law | 0.911 | 504 | 6.7 % | 0.91 | 0.5 ms |
| EG, fine grid, e-geodesic | 0.912 | 504 | 6.7 % | 0.90 | 0.5 ms |
| EGsh, EG shifted 15 s | 1.028 | 504 | | | |

Cells (A steady, A step, B steady, B step): EG 0.574, 0.757, 0.982, 1.208; EGsh 0.638, 0.824, 1.147, 1.352. The script's
docstring quotes AD's gain as 7.2 %, from rounded values; computed exactly here it is 7.1 %.

- E1 (EG at or below 0.95 of EGsh, below it in at least 18 of 20 runs): **held**. EG is 11.3 % below EGsh (0.887),
  in 20 of 20 runs (0.85 to 0.93): the timing holds on the path.
- E2 (gain at least 5 %): **held**, 6.7 %.
- E3 (within 1 point of the staircase's gain): **held**, 6.7 % against 7.1 %. The path neither adds nor costs.
- Described: the path replaces the 1 ms steps by 0.5 ms steps, twice as many, for the same total travel (about
  0.45 ms of frame change per joint per second in both). EG and GR are within 1 % of each other in every run (0.997
  to 1.010): the e-geodesic's speed bound almost never binds, because the law's raw window already moves slowly (its
  statistics average over about a second). EG uses 7 % fewer frames than AD and is 0.5 to 2.9 % less accurate run by
  run, on the same efficiency line.
- Reading: the continuous path is the same cadence as the staircase, made smooth, at the same efficiency and with the
  same knowledge of when. On this bench the geometry gives the form of the cadence, not more accuracy.

## Theorem Omega on the joint bench (pre-registered W1, W2)

Script `omega_phase_rates.py` (sha256 134102de76e3a15889f7bb2401fadee37c54563bbce2b94c1d9cce9ae0c11257), pushed before
its run and before any per-joint window series was looked at (commit 315302c); results `results/omega_phase_rates.json`.
The published bench (1 kHz, 3 ms frame, the engine with identification 'raw', and the tuned PID), baton arm. Phase rate
of a joint: 1 / tau_law, the law's window before the floor. Cancellation ratio C over a window: |J e|^2 summed over the
sum of the per-joint contributions |J_j e_j|^2 (1 if the contributions to the tip were uncorrelated, below 1 when
they cancel). Gate (12 joints, published seeds, syntonic and PID): difference exactly zero.

- W1, co-movement (mean pairwise correlation of the joints' log phase rates above the same statistic with the joints
  shifted against each other, in at least 9 of 10 steady runs): **held**, 10 of 10 (0.85 to 0.87 in excerpt A, 0.97
  in B, against -0.12 to -0.19 shifted). The weak prediction: the shared gesture alone may give it.
- W2, the reading (more dispersion of the phase rates, less cancellation: median Spearman over the 1 s windows at
  least 0.2, positive in at least 8 of 10 runs): **failed**. Median -0.10, positive in 4 of 10. The sign splits by
  excerpt: -0.19 to -0.48 in A (more dispersion went with more cancellation), -0.02 to 0.24 in B. With the window's
  tip speed partialled out (described): median -0.12. The agreement of the phase rates does not account for the
  cancellation.
- Described: the cancellation is real and specific to the syntonic engine. Over a run, C is 0.58 to 0.63 (A) and 0.70
  to 0.72 (B) for the syntonic engine, 0.94 to 0.98 for the PID, under steady noise (step: 0.63 to 0.77 against 0.98
  to 1.02). The PID's joint errors reach the tip almost as independent errors; the syntonic engine's cancel by about
  a third. Why stays open.

## Why the joint errors cancel: a velocity lag? (pre-registered L1 to L3)

Script `lag_cancellation.py` (sha256 6c8c9ebfd4337c671dc2eadf81503060a72e6e467822f3df7c214ce3876094f6), pushed before its
run (commit 9fb34ca); results `results/lag_cancellation.json`. Same runs as the Omega test; gate (every run reproduces
the Omega test's tip errors): difference exactly zero.

- L1 (the error's cancellation follows the gesture's, window by window: median Spearman at least 0.5, PID at least
  0.2 lower): **failed**. Syntonic median 0.32, PID 0.16.
- L2 (a velocity lag explains at least half of the joint error): **failed**, and clearly: median R2 0.003 (PID 0.001).
  The error is not a lag along the planned velocity.
- L3 (the lag part reproduces the cancellation within 0.05): **failed**, median gap 0.24.
- Described: the gesture's own cancellation of joint velocities at the tip is 0.37 in both excerpts.

## Exploratory passes after the verdict (no prediction; scripts in `exploratory/`)

- Implementation check (one run): the samples align with the plan to 1e-11 rad. The error correlates with the planned
  velocity 140 to 150 ms ahead or behind (|r| 0.42 to 0.69, the sign reversing between the two sides), the pattern
  of an acceleration term, not of a lag; 18 to 60 % of the error's variance lies above 5 Hz.
- `acc_explo.py` (sha256 76bfa20020de9d53..., 4 steady runs): regressed on the planned acceleration, the syntonic engine's joint
  error has R2 0.26 to 0.69 by joint (the PID's 0.03 to 0.34). That acceleration part cancels strongly at the tip
  (ratio 0.29 to 0.30), close to the gesture's own cancellation of joint accelerations (0.36). Reading, to be
  pre-registered before it counts: the syntonic engine's error is largely the gesture's acceleration seen through
  the servo's dynamics, which inherits the gesture's kinematic cancellation; the PID adds errors of each joint's own.
  If it holds, how the gesture is planned (its accelerations) is the lever on the tip error.
- `h_by_joint.py` (sha256 4f0b06ed978b7607..., 10 steady runs): the smoothness H the engine identifies is close across the six
  joints, 0.95 to 0.99 by joint (excerpt A: 0.97, 0.97, 0.97, 0.99, 0.95, 0.97; B: 0.98, 0.98, 0.99, 0.99, 0.97, 0.98),
  the spread across joints 0.10 in A and 0.013 in B, and at the engine's clip of 1.0 for 3 to 35 % of the time. One
  class for the block in the sense of the 4D field, likely inherited from the spline-planned gesture; at the edge
  where a first-order window is weakest.

## The acceleration reading on held-out seeds (pre-registered A1 to A3)

Script `accel_reading.py` (sha256 a1113086a5886c43ced14d95aab99779263ad16b9b39210d02076966c65232f6), pushed before its
run (commit f2f1e77), after the exploratory pass on seeds 1 and 2; results `results/accel_reading.json`. Seeds 6 to 10,
excerpts A and B, steady noise: 10 runs.

- A1 (the planned acceleration explains the syntonic engine's joint error, median R2 at least 0.35, the PID's at most
  half of it): **held**. Median R2 0.44 (by joint 0.35 to 0.70 in A, 0.26 to 0.58 in B); PID 0.14.
- A2 (that part cancels at the tip like the gesture's own accelerations, within 0.10): **held**. Median gap 0.066:
  0.29 for the error's acceleration part against 0.36 for the gesture's accelerations.
- A3 (the syntonic engine's cancellation at least 0.2 below the PID's in every run): **held**, smallest gap 0.26
  (0.58 to 0.71 against 0.94 to 0.98).
- Scope: the held-out seeds change the noise draws, not the gesture; the two excerpts are those of the exploratory
  pass. This confirms the reading over noise, not over new gestures, and it is a correlation. The causal test is to
  change the gesture's accelerations and see the error follow: planning the gesture differently (the open test on
  the e-geodesic).
- What it means: the syntonic engine's joint error is largely the planned acceleration seen through the servo's
  dynamics, and that part inherits the gesture's own kinematic cancellation at the tip; the PID's error is mostly its
  joints' own. That is why the engine matches the tuned PID at the tip with joint errors about 28 % larger.

## The causal test: re-timing the gesture (pre-registered K1 to K3)

Script `plan_intervention.py` (sha256 38d4b60bb66a4075b077fbb888f5059c3b1c934477ecce04a3621038064249b6), pushed before
its run (383f18e); results `results/plan_intervention.json`. The same path re-timed: tempo 1 (B), uniform tempos 0.7
(U70) and 0.5 (U50), and two tempo paths from 0.5 to 1, rate linear in time (EG, the e-geodesic of the tempo, a
reading) and period linear in time (PL). Two engines: the published adaptive window, and the same engine with its
window frozen at its baseline mean. Seeds 11 to 15, excerpts A and B, 100 runs, 39 minutes; the gate reproduced
omega_phase_rates.run exactly in both excerpts. k = 1 when the tip error locked to the plan scales by the tempo squared.

- K1 (window frozen, slower tempos, 0.8 <= k <= 1.25): **failed**. k = 0.71 (U70) and 0.67 (U50) in A, 0.33 and 0.31
  in B. The failure is on the side opposite to the rivals (velocity 1.43 and 2, pose 2.04 and 4): the tip error locked
  to the plan falls faster than the tempo squared.
- K2 (window frozen, the tempo paths): **failed**. k = 0.94 (EG) and 0.92 (PL) in A, inside the band; 0.55 and 0.54
  in B, outside.
- K3 (the law's trade): **held**. Window frozen, the part not locked to the plan does not move: N(U50) / N(B) = 1.01
  (A) and 1.04 (B). Adaptive window, it falls: 0.69 and 0.77.
- Described. With the window frozen the joints' acceleration gains do not move with the tempo (ratio to baseline, by
  joint, median over seeds: 0.99 to 1.07 at U70 and 1.01 to 1.16 at U50 in A; 1.00 to 1.10 and 1.00 to 1.21 in B).
  With the adaptive window they grow as the window widens (1.29 to 2.36 at U50 in A, 1.18 to 1.95 in B). Share of the
  baseline tip mean square locked to the plan, P / (P + N): 0.76 (A) and 0.90 (B) frozen, 0.65 and 0.84 adaptive;
  the acceleration part alone was 0.42 and 0.24 frozen at baseline on seed 1 (exploratory), so in B most of the tip
  error locked to the plan is not the acceleration part. Tip error, median mm, arms B, U70, U50, EG, PL:
  A adaptive 0.741, 0.506, 0.502, 0.651, 0.607; A frozen 0.817, 0.533, 0.448, 0.668, 0.632;
  B adaptive 1.274, 0.625, 0.529, 0.670, 0.621; B frozen 1.359, 0.634, 0.512, 0.688, 0.638.
  The adaptive window beats the frozen one in every arm except U50, in both excerpts.
- A reading to test next, not a result: at the published tempo a joint nears its slew limit (RJ4's plan alone reaches
  0.84 of it in A and 0.92 in B, before the correction is added), and the error of those moments vanishes once the
  gesture slows. That would give k below 1 at both slow tempos, a stronger effect in B, the faster excerpt, and K2
  failing in B, where the tempo paths end at tempo 1. The check is the share of time at the slew limit by arm and its
  share of the tip error.
- What it means. The joint error follows the acceleration with a gain that belongs to the loop (described, the
  causal part at the joints), and the law trades noise for drift as Theorem U' says (K3). The tip error locked to the
  plan is not only the acceleration: at the published tempo of excerpt B it is mostly a part that disappears at tempo
  0.7. The prediction of a pure tempo-squared scaling at the tip is refuted.

## Exploratory passes after the causal test (readings, not tests)

Seeds 1 to 5 (all used before) and noise-free runs, frozen window as in `results/plan_intervention.json`; nothing here
is judged, and a test of what follows takes fresh seeds.

- The slew-limit reading is **refuted** (`exploratory/slew_explo.py`): no joint reaches its slew limit at any tempo in
  either excerpt, and removing the limit leaves k unchanged (0.72 and 0.67 in A, 0.33 and 0.31 in B, with or without).
- Without encoder noise the bench gives the plan-locked error itself and the same k (0.71 and 0.67; 0.33 and 0.31):
  the effect is deterministic (`exploratory/loop_dynamics_explo.py`).
- What it is: the engine's correction limit. At tempo 1 the correction of RJ4 and RJ6 sits at its 5 degree limit 0.6
  to 0.7 % of the time in A and 3.7 % in B, never at tempo 0.7 or 0.5. With the limit at 30 degrees
  (`exploratory/limit_explo.py`, noise-free), the tip error at tempo 1 falls from 0.71 to 0.64 mm in A and from 1.29 to
  0.81 mm in B, and k becomes 0.87 and 0.81 in A, 0.82 and 0.76 in B.
- What stays below 1 is the loop itself. The correction (u' = e / tau) and the servo (time constant TS = 25 ms) form a
  second-order loop, TS tau e'' + tau e' + e = TS tau qddot: static gain TS tau (7.5e-5 s^2 for a 3 ms window, the
  order of the measured gains), natural frequency 13 to 21 Hz and damping 0.5 sqrt(tau / TS) = 0.15 to 0.24 for the
  frozen windows, lightly damped. That linear loop alone gives k 0.90 to 0.95 in A and 0.87 to 0.92 in B, and
  reproduces 89 to 96 % of the noise-free tip error at the slow tempos; it has no frame hold and no quantization.
- Readings for the next act, not results. (1) The correction must carry the servo's lag, about TS times the joint
  speed, which exceeds 5 degrees on RJ4 and RJ6 at the published tempo. A limit stated in the units of the authority
  fact, TS times the no-load speed (12.5 degrees for RJ4, 25 for RJ6), would not saturate. The main statement of 26
  October lists the limiter among the design constants identical on every machine, so changing it touches the
  portability claim; whether to restate the limiter from the authority fact is JP's decision. (2) The law's windows
  (2 to 6 ms) are 4 to 11 times shorter than the servo's time constant, which leaves the loop lightly damped. The law
  decides when and the stability stays with the machine, as the main statement says; the damping is the machine's
  side of that split, now measured.
