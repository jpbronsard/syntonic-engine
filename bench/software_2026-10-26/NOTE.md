# Software bench of the syntonic engine against adaptive filters

Bench note, version 1.5. 26 and 27 October 2026. SyntonicAI Recherche.
Scripts, tables and per-seed rows: `bench/software_2026-10-26/` of the syntonic-engine repository (the
bank checkpoints, 110 MB each, are recomputed by the scripts).
Version 1.1 adds section 8: the engine rows re-run with the identification exactly as released in
syntonic-engine 0.2.0 (commits 17a1536 and a9cc957 on GitHub), which fits the six finest scales that clear the
noise floor where the bench scripts fitted all of them. No conclusion changes; the numbers that move are listed.
Version 1.2 (27 Oct) adds sections 9 to 11, the quadrotor substrate in software: step 5 (the rate loop, where the
law has nothing to decide), step 6 (the attitude estimator, where the law holds and the engine's identification
fails) and step 6b (the attitude estimator with the gyro drift as a fact measured once at rest, where it passes).
Version 1.3 (27 Oct) adds section 12, step 7: cross-axis coupling on logged data of the X500 (does the attitude
estimator, or the rate path, need the matrix law?), with the TH44 collapse certificate on three axes.
Version 1.4 (27 Oct) adds section 13, step 6c: a pre-registered ablation of the attitude-estimator law (the
exponent, a fixed horizon, a x256 rise of the vibration). Four of five hypotheses hold, the exponent test at the edge
of its band; the predicted window fails at x256.
Version 1.5 (27 Oct) adds section 14, step 6d: the mechanism of the step-6c failure, and whether the generalized law
U' predicts a change of exponent. Five of seven hypotheses hold; the two mechanism hypotheses fail.

## 1. Purpose and verdict

Question: does the syntonic engine (jpbronsard/syntonic-engine, Python v0.1.1), which sets its window from
Theorem U' with no tuned parameter, hold an advantage over classical adaptive estimators? Six pre-registered
steps on synthetic fractional Brownian targets, each with its hypotheses and decision rule written before the
run and the script's SHA-256 printed at the top of its results.

Verdict, in six lines.

1. Theorem U' holds exactly in this class: measured exponents 0.635, 0.501, 0.416, 0.371, 0.346 against
   predicted 0.625, 0.500, 0.417, 0.370, 0.345 for H = 0.30, 0.50, 0.70, 0.85, 0.95 (step 0).
2. The published engine's identification (raw-read structure function at 2 to 64 ms) fails when the read
   noise is large enough that the signal is below the noise floor at those lags; an Allan-type (block-mean)
   structure function fixes it: the engine then lands within 7 % of the swept optimum in 35 of 35 cells
   (step 0b, exploratory). This is the v0.2 change.
3. In stationary regimes the engine is at parity with likelihood-selected Kalman banks (1 to 7 % worse than
   the best fixed window; the banks 0 to 2 %). A random-walk Kalman bank finds the optimal window at every H
   without reading H. On smooth targets (H >= 0.85) a constant-velocity Kalman beats every first-order window
   by 25 % (step 1).
4. After step changes of regime, the engine beats plain Bayesian banks with forgetting (replicated), mainly
   because it tolerates sensor-noise increases that make a naive bank run to fast models; against a bank with
   a textbook change detector the advantage does not survive: the engine keeps noise increases and loses noise
   decreases (steps 2 and 3). A parameter-free META noise estimator (Identification Duality, Theorem 2) does
   not replace a detector for steps and makes the x16 decrease worse (step 3c).
5. Under slow log-random-walk drift the engine is at parity; under fast drift it is behind the banks, and the
   META level is worse still (step 4). All four drift hypotheses failed.
6. What remains of the engineering branch: a single filter with no per-application tuning, at parity with
   banks of 61 to 158 filters wherever its first-order class fits, graceful under noise increases, blind to its
   own class limits. Its cost advantage is plausible and unmeasured. The branch is closed.

## 2. Common setup

Targets: fractional Brownian motion of index H (Davies-Harte exact synthesis, checked to 2 % on 30 paths),
structure function lambda_H s^(2H), lambda_H = 1, dt = 1 ms, 2^16 samples. Two independent white Gaussian
reads per tick. Noise levels placed so that the U' optimum sits at 7 log-spaced windows between 5 and 160 ms.
Yardstick: the best fixed exponential window on the same stream (sweep of 40 windows, parabolic refinement),
or a per-segment / per-5-s-block version of it where the regime changes. Scores on t > 5 s. Ten seeds per cell
(0..9 for steps 0 to 1, 100..109 for steps 2 to 4, 200..209 for the replication).

Engine: the identification arithmetic of engine.py applied as a smoother (window on the pair mean), engine
defaults except tau_max raised from 0.5 to 1.0 s so that no 160 ms optimum is masked, and the floor set to
dt/0.45 (no servo). Step 0b onwards: block-mean (Allan) structure function at 2 ms to 4.096 s, noise floor
n2 dt / l subtracted analytically, C(H) = (2^(2H+2) - 4) / ((2H+1)(2H+2)).

Competitors: Bayesian multiple-model Kalman banks (Magill 1965), random-walk model (61 values of q/r) and
constant-velocity model (97 values of q dt^3/r), posterior by cumulative innovation likelihood, with
exponential forgetting (memory 1 s and 4 s) from step 2, and with a change detector from step 3 (noise input:
50 ms over 1 s ratio of the squared read difference beyond 2 or 1/2 for 100 ms; NIS input > 3 for 100 ms;
action: noise re-initialisation and posterior reset; refractory 1 s). Union of both banks as the model-class
adaptive competitor (named IMM in the scripts). Fixed window tuned on the pre-jump regime as the practitioner
baseline.

## 3. The steps

Step 0, pre-registered, sha 7d1d419d. H1 (U' exponent within 0.05): PASS at all five H. H2 to H4 (engine
within x1.5 of the optimum, regret <= 10 %, H within 0.1): FAIL above the two lowest noise levels at every H.
Cause, quantified: with raw reads, the lag at which the signal clears the noise floor grows from 0.1 s (level
0) to 100 s (level 6, H = 1/2); the engine's lags stop at 64 ms. The published engine is correct in the regime
it was validated in (servo encoders, where the signal exceeds the floor by orders of magnitude) and silent
about this limit. v0.2 should state it and carry the fix below.

Step 0b, exploratory, sha 936fc8b9. Allan identification. 35 of 35 cells: tau_eng / tau_opt in [0.93, 1.07],
regret <= 6.5 %, H within 0.1 except one cell (H = 0.70, level 6: 0.59). Negative regret at H = 0.95 (the
moving window slightly beats any fixed one).

Step 1, pre-registered, sha fc13474b. Method amendment before the run: innovation-covariance matching (first
draft) was noise-dominated in a sanity check at H = 1/2 (q under-estimated by 37 %, +44 % MSE) and replaced by
the Bayesian bank; no result existed when the change was made. C1 (parity at H = 1/2): PASS. C2 (wrong model
costs >= 20 %): FAIL, the random-walk bank is within 2 % of the optimum at every H. C3 (engine within 10 % of
the better single-model bank): FAIL at H = 0.95 only, the constant-velocity bank beating the best exponential
window by 25 % and the engine by 35 %, as named in the pre-registered risk. C4 (union finds the class): PASS.
Q5 (H = 0.70): the engine does not beat the union (3 to 5 % behind).

Step 2, pre-registered, sha 397c22b5, replicated on fresh seeds, sha 4e34029f (seeds and output folder as
environment variables, nothing else). Jumps: read noise x4 (J1) and x16 (J2), signal amplitude x16 (J3), H
switch 0.5 -> 0.85 -> 0.5 (J4), each at 22 s and back at 44 s. Decision rule: a method wins a cell when its
adaptation excess (5 s MSE over its own steady state) is lower by more than 10 % and its recovery is not
longer; the claim "re-adapts faster without a tuned memory" is supported only with a majority of the 14 cells
against both bank memories. Result: 8/14 and 13/14 (replication 9/14 and 12/14): supported. Pattern: all eight
noise-jump cells to the engine in both passes (after x16 the engine loses 3 to 9 % for 5 s and recovers in
0.6 s; the banks lose x6 to x7 and recover in 2 to 5 s); the amplitude jumps (J3) to the 1 s bank (engine
recovery up to 2 s); the H switch split. T2 (fixed window loses >= 20 % after a x16 jump): PASS, 78 to 154 %.
T1 (all adaptive methods re-converge within 10 %): marginal FAIL for the 1 s bank at one cell (10.5 %, 16.4 %
in the replication), never for the engine.

Step 3, pre-registered, sha 9086d5a3 (table reprinted unchanged as 609a7928 after a JSON cast fix). Hardened
banks. D1 (detector fires within 1 s after both noise jumps in >= 90 % of J1/J2 streams, < 3 false alarms):
PASS, 100 % and 0.00. T3h: 5/14 and 9/14: NOT supported. The engine keeps the four noise-increase cells (its
excess 1 to 6 % against the hardened banks' 11 to 23 %) and loses the four noise-decrease cells (its 1 s noise
average decays slowly; the detector re-initialises the bank at once). J3 and J4 as before.

Step 3b, drafted and withdrawn before running: a thresholded detector inside the engine would add parameters,
which contradicts the engine's definition. Recorded for the trail.

Step 3c, pre-registered, sha f214d0f0. The engine's noise estimate supplied by a second syntonic engine on
log d^2 (known scatter pi^2/2, offset -1.2704), the META level of the Identification Duality, with no new
constant. M1 (steady state not worse): PASS. Rule outcome: MATCHES (6 to 6 against the 1 s bank, 9 to 3 against
the 4 s one). But the written prediction failed where it mattered: on the x16 noise decrease the META engine
is worse than the plain one (excess x1.4 to x2.0 against +22 %, recovery 2.2 s against 0.7 s). Mechanism: the
meta engine reads the step in log d^2 as a large innovation, shortens its window to 25 to 45 ms and keeps it
short for seconds because its block averages stay contaminated; the noise estimate then carries +/- 40 % scatter
and the main window shakes with it. Theorem 2 assumes a log random walk; a step is outside its class.

Step 4, pre-registered, sha e5adce2e. Both log sigma^2 and log lambda_H as independent log random walks,
e-folding 10 s (slow) and 2 s (fast), weak mean reversion, H in {0.50, 0.85}; regret against the best fixed
window per 5 s block. R1 (detector quiet under drift): FAIL in fast drift (12 firings per stream). R2a (META
helps in fast drift): FAIL (+74 % against +57 % at H = 0.5). R2b (parity with the better bank memory within
10 %): FAIL by a hair at slow drift, H = 0.85 (engine +7 %, banks -4 % through the constant-velocity class), and
clearly in fast drift (+57 % and +41 % against +34 % and +24 %). R2c (beats the worse memory in 3 of 4 cells):
FAIL, 0 of 4; in slow drift the bank's memory does not matter (+8 % and +9 %).

## 4. What stands, what fell

Stands: the U' exponent, exactly, over five memory classes and 1.5 decades of window; a parameter-free engine
that lands on its class optimum across H once its identification works at the right scales; graceful behaviour
under sensor-noise increases; parity with far heavier adaptive banks wherever the first-order class fits.

Fell: "a wrong-model Kalman is off the optimum" (a likelihood-selected bank is not); "reading H is an
advantage" (it is one route among others); "the law alone re-adapts faster" (not against a detector); "the META
level replaces a detector" (it does not, for steps); "drift is the law's home ground" (parity at best).

Limits now quantified: first-order windows lose 25 % to a second-order filter on smooth targets (H >= 0.85);
the engine cannot notice when its class is wrong, because it computes instead of comparing; the engine's noise
average is slow to follow a decrease; the long-scale block averages make it slow after signal-amplitude jumps.

## 5. Recommendations for the engine

v0.2: Allan-type block-mean structure function with analytic floor subtraction and C(H) (step 0b code), scales
from 2 ms to a few seconds, per-stream block counters. README: state the validated domain of v0.1.1 (signal
above the read-noise floor within 64 ms, the servo case) and the extension in v0.2. CHANGELOG entry naming the
limit found in step 0. The "zero configuration" wording stays as an engineering property; the design constants
(scales, SNR gate 2, averaging 1 s, warm-up 1 s, rate limit 0.4 s/s, floor, ceiling) are listed as such.

Cost: the parsimony argument (one filter against 61 to 158 plus a detector) is plausible and unmeasured. Before
it is written anywhere: CPU time per tick, memory and allocations, for the engine, a lean bank (about ten models
per class with the same detector) and the full bank, on the laptop and on a Teensy 4.1.

## 6. Reproduction and audit status

All scripts are deterministic (seeded numpy, no threading dependence). In order, from `syntony_bench/`:
`step0_uprime_engine.py`, `step0b_allan_identification.py`, `step1_vs_kalman.py` (about 30 min),
`step2_regime_jumps.py` (about 12 min; replication with `STEP2_SEED0=200 STEP2_OUT=step2_rep`),
`step3_hardened_bank.py`, `step3c_meta_noise.py`, `step4_drift.py`. Each prints its SHA-256 first; the tables
should reproduce to the printed digits. Dependencies: numpy, scipy, and the engine's `engine.py` arithmetic
(copied verbatim into the smoother classes; the package itself is not imported).

Audit status: the code was written and run by the assistant inside the working session, with pre-registration
in the docstrings and hashes printed before each run; it has not been reproduced or audited independently.
The author's own re-run on his machine is the first audit and should be recorded here when done.

## 7. What follows

The engineering branch is closed with this note. The programme's claim moves to the one place where the law is
not a re-description of an optimum someone designed in: the same kernel, the same three hardware facts (loop
timing and delay; actuator authority; access to the read noise), on substrates of different physics, with the
same pass criterion as step 0b (window within a factor 1.5 of the swept optimum, regret under 10 %, exponent
against U'): ARM-1 (servo, real encoder noise and delay), the quadrotor on a gimbal (unstable plant, smooth
gyro signals where the first-order class may fail and the matrix law may be needed), and HUMAN-1 (the substrate
nobody designed, with the synthetic positive and negative controls of the estimator validated first). A failure
on any substrate is published as such and marks where the law stops.

## 8. Version 1.1: the engine rows with the released 0.2.0 identification

The engine shipped in syntonic-engine 0.2.0 differs from the bench smoother of steps 0b to 4 in one respect,
chosen after the servo bench of the package showed the all-scales fit losing a few percent in the servo regime:
the log-log fit of the block-mean structure function uses only the six finest scales that clear the noise
floor, so that the scaling it measures is the local one, at the operating scale (Theorem U' is stated at that
scale). Every engine row of the bench was recomputed with that identification against the stored bank, oracle
and fixed-window rows (`rerun_engine_v020.py`, log `rerun_engine_v020_results.txt`). The banks did not change.

Stationary (steps 0b and 1): window within 0.93 to 1.08 of the swept optimum in all 35 cells, regret at most
6.5 %; the ratio to the better single-model Kalman bank moves by at most 0.02 and keeps the same verdicts
(C3 still fails at H = 0.95: 1.27 to 1.37).

Step changes (step 2): 8 of 14 cells against the 1 s bank and 14 of 14 against the 4 s bank (version 1.0:
8 and 13); replication 10 and 12 (version 1.0: 9 and 12). T1 now passes for the engine in every segment. The
one visible gain is on the signal-amplitude jumps (J3): steady-state regret in the high-amplitude segment
falls from 4.4 % to 1.3 % at H = 0.5 and from 3.0 % to 0.2 % at H = 0.85, and the transient excess from 0.22
to 0.13; the 1 s bank still wins the first J3 jump. The finest-scale fit stops stale long-scale blocks from
holding the window back; it does not shorten any average, so nothing is paid in steady state.

Hardened banks (step 3): 5 of 14 and 10 of 14 (version 1.0: 5 and 9); the claim remains unsupported. The
META noise estimator (step 3c): 6 and 10, "matches" by the rule, and the x16 noise decrease remains its
failure (transient excess 1.5 and 2.1).

Drift (step 4): slow drift +6.4 % and +6.6 % (version 1.0: +7.1 %, +7.0 %); fast drift unchanged (+57 %,
+40 %); META worse in fast drift than in version 1.0 (+94 %, +56 %). All four drift hypotheses still fail.

Package state: Python, C++ and Rust carry both identifications; replays over 33 s and 12 joints agree with
Python to 7.8e-16 rad (C++) and 8.9e-16 rad (Rust) with block means, 1.5e-14 and 5.4e-14 rad with the raw
identification, unchanged. On the package's own servo bench the block-mean identification is a few percent
behind the raw one, mostly during the noise step (tip-error ratios 1.02 / 0.78 / 1.07 / 0.91 against
1.00 / 0.69 / 1.07 / 0.89); both figures are in the README.


## 9. Step 5 (27 Oct 2026): the rate-loop horizon of a quadrotor cascade

`step5_quad_rate_horizon.py` (SHA-256 454b9174, pre-registered; one amendment written into the docstring after a five-run smoke test and before the first full run: a second plant stiffness B = 30 000 deg/s^2 per unit was added, R2 was restated in the form of step 0b, and a steady-state score J1h was added for reporting). The plant is the one-axis body of `dRehmFlight_syntonic/sim_compare.c` at 2 kHz; the controller is that graft's cascade (angle loop fixed at 2 ms, rate loop with the horizon under test, authority 0.25). The horizon of the rate loop was set in turn by a fixed sweep (13 values, 2 ms to 0.5 s), by `Syntonic.h` as written (the raw proxy: sigma and lambda of the raw loop error), and by the 0.2.0 measurement layer (known gyro-noise variance, block-mean structure function with the floor subtracted) with H fixed at 1/2 or fitted, with and without a cap derived from the loop. Three gyro-noise levels, five seeds, two plants. Results in `results/step5_quad_rate_horizon_results.txt` and `results/step5_rows.json`.

- In every cell the fixed sweep puts the best horizon for the angle error at the floor (2 ms), and the error rises monotonically with the horizon until the loop goes unstable near 300 ms. The reason is structural: the body integrates the command, so rate-loop noise costs almost nothing in angle while any lag costs tracking. The steady-state score says the same. Only the command jitter J2 falls with the horizon.
- The raw proxy sits at the floor 83 to 100 % of the time (R1 failed only at the lowest noise on the soft plant, 83 %). The measurement layer sits there too at 1 and 3 deg/s and lifts to 5 to 9 ms at 10 deg/s, with +1 to +7 % angle error and 43 to 57 % less command jitter. By the pre-registered rule the measurement-layer claim for the quad rate loop fails (R3 missed: the jitter gain is under 30 % at 3 deg/s), and the loop cap is not supported (R4: the law never went long enough to need it; the fixed sweep shows where instability starts, so a cap at a third of the outer time constant remains harmless insurance, not a necessity). R2 is held trivially: no cell has an optimum off the floor to land on.
- U' changes nothing here (R5 held), and the fitted H on a closed-loop error is not a power-law exponent at all (medians 0.05 to 0.79 across cells): this signal is not in the law's class as a stationary process.

What it means. The law computes the window that best estimates the signal. In a quadrotor rate loop with an ideal actuator the best window for control is not the best window for estimation: the plant discards the rate noise, so the control optimum is the loop-delay floor whatever the noise. In flight, tau at the floor was therefore the right answer, not a symptom. The place where a rate-loop horizon matters is downstream of the loop: motor command noise, heat, frame resonances, actuator bandwidth, none of which this plant has. So the quadrotor substrate cannot test the horizon law by angle error alone; its criterion has to carry the cost of command noise (or a plant in which that noise costs angle), chosen before the sweep. That choice is a design weight, not a machine fact, and the master protocol has to name it as such. The soft X500-class plant (inner pole 4 /s under an outer pole of 7.5 /s) also shows that any inner lag is paid immediately, which is the measured form of the authority condition B g / Y >> outer pole.


## 10. Step 6 (27 Oct 2026): the attitude estimator of a quadrotor on a gimbal

`step6_quad_attitude_estimator.py` (SHA-256 035647c9, pre-registered in its docstring; a smoke test of two seeds printed only the oracle windows, no amendment, the full run carries the same hash; the smoke files in `results/` hold every row, unread before the full run). Judge `evaluate_step6.py` (its printed label was corrected once after the run, as for step 6b below: the file path and "script" SHA-256; no rule changed, the verdict regenerated identically). Results in `results/step6_results.txt` and `results/step6_rows.json`; exploratory diagnosis in `step6_diag1_structure_function.py` to `step6_diag3_hadamard.py`.

Step 5 placed the law out of the rate loop. The place where a window costs on both sides is the attitude estimator, the fusion of gyro and accelerometer. One axis at 1 kHz on a gimbal (no translation): an Ornstein-Uhlenbeck attitude (5 deg, 2 s); a gyro with angle random walk Q_G = (0.01 deg/s/sqrt Hz)^2 and a bias random walk; two IMUs whose accelerometer tilt carries common-mode vibration (a white low-frequency part scaled by a level V and a rotor line at 92 Hz V^(1/4)), independent sensor noise and an 80 Hz anti-alias filter. Three scenarios of three 120 s segments, ten seeds: steady (V = 1, 1, 1), throttle (V = 1, 16, 1/16), and an uncalibrated gyro bias of 0.3 deg/s. The values are representative, not measured; the gimbal bench will measure them.

The branchement, declared before the run: the two sensors reduce to one signal, z = accel tilt minus the open-loop integral of the gyro. Its slow part is the gyro drift (lambda) and its fast part the accel noise (sigma^2); the true attitude cancels. The window sets Kp = 1/tau of a Mahony filter, with the gyro-bias integral tied to it by the damping of the flown firmware (Ki = Kp^2 / (4 zeta^2), zeta = 1.826). Candidates: the oracle (fixed window swept, best per segment); the law with the true noise levels (SYN-K); the engine 0.2.0 identification on one stream, noise floor read at 16 to 128 ms (SYN-1S); the engine unchanged with its noise read from the two IMUs (SYN-2R); the law without the bias integral (witness); a two-state Kalman filter (angle, gyro bias) tuned once per stream (KF-OPT), with its noise fixed at nominal (KF-FIX, EKF2's fixed parameters) and switched with the true level (KF-CLAIR); a one-axis reduction of Madgwick; the flown `attitude_update` (0.5 s, no bias integral).

- The law holds (H6.1, 9 of 9 cells). With the true noise levels the window lands within 0.88 to 1.28 of the oracle, regret 0.2 to 8.1 %, at 1.00 to 1.01 of the tuned Kalman under steady vibration and within 1.2 % of the clairvoyant one when the vibration changes x16. This is the expected equivalence of the first-order window with the steady-state Kalman gain; what it establishes is the branchement.
- The engine's identification fails (H6.2, H6.3, H6.4, H6.5, H6.7). Measuring noise and drift itself, it sets the window 3 to 43 times too short (fitted H 0.05 to 0.07 in eight of nine cells), regret 33 to 456 %.
- The predicted failure happened (H6.6). Reading its noise from the difference of two IMUs, the engine sits at the 2.2 ms floor with errors 23 to 31 times the oracle: the difference cancels the common-mode vibration, which is precisely the noise of the accelerometer channel in attitude fusion. The same holds for any inter-IMU noise estimate used on that channel.
- What adaptation is worth here: the fixed-noise Kalman loses 36 % against the clairvoyant one when the vibration rises x16 and 12 % when it falls. The flown filter loses 14 % at nominal vibration and is about 14 times worse with the uncalibrated bias; the best one-axis Madgwick loses 38 %, dRehmFlight's beta = 0.04 about 2.6 times.

Exploratory diagnosis, after the verdict. The bias integral, tied to the window, removes the slow drift from z beyond about 13 tau, so the structure function bends down at long scales (8, 16, 32 s: 1.3e-6, 7.2e-7, 3.9e-7 rad^2 on one seed): the fit reads H = 0.05 and lambda a hundred times too large, the window shortens, the integral speeds up, and the bending grows. Without the bias correction the accumulated bias random walk makes a ramp in z (structure function rising as l^1.8 from 4 to 32 s), outside the class. And the drift is visible only beyond about 2.5 to 3.5 times the optimal window: lambda cannot be measured live from minutes of data, while the noise level is measured within 12 % (second differences of block means at 32 to 128 ms).

## 11. Step 6b (27 Oct 2026): the gyro drift as a fact, the accel noise measured live

`step6b_gyro_fact.py` (SHA-256 aa082341, pre-registered; smoke test printed only the rest measurement; no amendment). Judge `evaluate_step6b.py` (its printed label was corrected once after the run, "rows" to "script" SHA-256; no rule changed, the verdict regenerated identically). Results in `results/step6b_results.txt` and `results/step6b_rows.json`. Same generator, scenarios, seeds and scoring as step 6; the oracle, Kalman, Madgwick and flown-filter rows are read from step 6 unchanged, and the flown filter recomputed here matches them exactly (identity check, maximum difference 0).

SYN-6B: lambda is the gyro's angle random walk measured once from a 300 s static record (Allan variance at 0.128 to 1.024 s; lambda_hat / Q_G = 0.955 to 1.021, Allan-deviation slope -0.507), H = 1/2 by the sensor physics; the noise level is measured live from second differences of block means of z built on the raw gyro (no bias feedback), the minimum over 32, 64 and 128 ms; the law, the Mahony filter and its tied bias integral are those of step 6. Sensitivity: lambda multiplied by 2 and by 1/2.

| Cell | Oracle (window) | SYN-6B (window) | lambda x2 | lambda /2 | KF-OPT | KF-FIX | KF-CLAIR |
|---|---|---|---|---|---|---|---|
| steady, seg 0 | 10.8 (1.14 s) | 11.2 (0.87 s) | 12.0 (0.62 s) | 10.9 (1.23 s) | 10.8 | 10.9 | 10.9 |
| throttle, V x16 | 20.6 (3.72 s) | 21.5 (3.32 s) | 22.4 (2.34 s) | 23.0 (4.73 s) | 22.2 | 29.0 | 21.4 |
| throttle, V /16 | 7.1 (0.54 s) | 7.3 (0.40 s) | 7.9 (0.29 s) | 7.1 (0.57 s) | 10.6 | 7.9 | 7.0 |
| gyro bias, seg 0 | 11.2 (0.84 s) | 11.3 (0.87 s) | 12.0 (0.62 s) | 13.5 (1.23 s) | 10.8 | 10.9 | 10.9 |

RMS angle error in millidegrees on the steady parts, mean over ten seeds; windows are seed medians.

- H6b.1 holds: within 0.75 to 0.93 of the oracle window, regret 1.8 to 3.8 %, in the six cells of the steady and throttle scenarios. Per seed the regret exceeds 10 % in 3 of 60 seed-cells (10.1, 10.1 and 14.9 %).
- H6b.2 holds: across the throttle steps the window follows the noise with exponents 0.505 and 0.474 (0.5 predicted).
- H6b.3 holds: at 1.02 to 1.04 of the tuned Kalman under steady vibration.
- H6b.4 holds: 26 % better than the fixed-noise Kalman when the vibration rises x16 (10 of 10 seeds) and 7 % better when it falls; within 0.7 to 3.7 % of the clairvoyant Kalman. A single Kalman tuning per stream cannot serve both levels: SYN-6B beats it by 7 % at nominal and 31 % at low vibration.
- H6b.5 holds: with the uncalibrated gyro bias, regret 0.1 to 3.0 % once the bias is learned.
- H6b.6 fails: with lambda halved the criterion holds everywhere, with lambda doubled it fails in all six cells (window 0.54 to 0.65 of the oracle, regret 5.7 to 13 %). The live noise level is biased low by the minimum over three noisy scales, and an overestimated drift adds to it. The drift fact must be measured at rest, not read from a data sheet that may be off by two in that direction.
- Reported, not judged: after the x16 rise the transient RMS (first 40 s) is 23.7 against 31.9 for the fixed-noise Kalman and 22.5 for the clairvoyant one; at start with the uncalibrated bias the bias integral is slow, 86 against 12.6 for the Kalman, which carries an explicit bias state. The bias channel, with its gain tied to the window, is the weak point.

What it means. On the attitude estimator the window law holds and reaches parity with a tuned Kalman, and it keeps the optimum when the vibration changes, where a Kalman with fixed noise parameters loses up to a third. The engine's own identification does not carry to this substrate: the drift must enter as a fact measured once (the gyro's angle random walk), and only the noise is measured live. The portability hypothesis in its strict form, the same code identification included, is therefore refuted on this substrate in software; it survives as the same law with facts measured once per machine. Not tested here: the U' exponent measurement (H is fixed at 1/2 by the sensor), gyro noise that rises with throttle, and everything a real gimbal adds (lever arm, frame resonances, measured vibration spectra). The gimbal bench, with an encoder per axis as truth, is the next test.

## 12. Step 7 (27 Oct 2026): cross-axis coupling on logged data, and the TH44 certificate

Scripts `step7_axis_coupling.py` (step 7a, SHA-256 ff76fe25, pre-registered), `step7c_rate3_th44.py` (step 7c,
d1f43d53, pre-registered), `step7c_A1_exploratory.py` (36b249fb) and `step7c_A1_sensitivity.py` (02720e47),
the amendment written after step 7c, and `step7_exploratory_symmetric.py` (10ebc0bf). The two pre-registered
scripts print their SHA-256 at the top of their results; the three exploratory ones do not, and their hashes are
the ones listed here. Run all of them with `run_step7.py`, which unpacks `data/`. Results in `results/step7*`.
No smoke run; nothing was changed in a pre-registered script after it ran.

The question. Steps 6 and 6b run the law one axis at a time. The matrix law (Syntony Principle V5.1, Theorem
35.1: T Q T = kappa^2 Sigma, T = kappa (Q^-1 # Sigma), the geometric mean of Pusz and Woronowicz; kappa = 1
below) is proven for coupled axes but not in the engine. Before building it: are the axes of this machine coupled
enough to matter? An EKF estimates all axes jointly with a full covariance, so the same numbers say what that
joint estimate is worth against the law run axis by axis.

Method. No filter is run on the logs, which carry no truth. The noise and drift matrices are measured from
second differences of block means, and the steady-state error of each window is predicted from them (filter
xhat' = T^-1 (y - xhat) on a random walk measured with white noise, continuous Lyapunov equation, score tr P,
a mean-square quantity). The build rule was written in mean square; RMS equivalents are given beside it
(RMS = sqrt(1 + g) - 1 for a mean-square excess g).

Data. Three sessions of the X500's SD logger (Pixhawk 6X, accelerometer at 250 Hz, gyroscope at 1 kHz),
reduced in `data/` to the columns used, no position: session 6, motor test with 2 kg of ballast per leg,
frame still, 0.5 to 1.2 A; session 5, propulsion on a held frame, 15 to 23 s at 0.436 throttle and 9.8 A;
session 7, airborne from about 19 to 24.5 s, then a throttle cut and a crash whose cause lies outside the
estimator. The flight window is not steady: throttle 0.41 to 0.56, current 10 to 20 A. The telemetry of the
earlier flights (M5, M6, about 50 Hz) logs only per-axis noise and drift, with no cross-axis terms, so those
flights cannot be examined this way. The SD log carries no magnetometer values, so yaw is covered on the rate
side only.

Step 7a, the attitude estimator, roll and pitch (session 6). Sigma, the 2x2 accelerometer-tilt noise, and Q, the
2x2 gyro drift. Correlations: Sigma -0.05 and +0.09, Q +0.16 and -0.16 (idle, low throttle); commutator ratio
0.02 and 0.06. Per-axis windows lose 0.58 % (idle) and 0.82 % (low throttle) in mean square against a two-axis
Kalman filter (RMS 0.29 and 0.41 %). Rule D2 (5 % in mean square): the per-axis law is enough, no matrix-law
bench for the estimator on these data. Limits: low throttle only, a still frame, windows of 16 to 52 s
extrapolated beyond the scales measured (up to 1 s).

Step 7c, the rate path, roll, pitch and yaw (pre-registered). The 3x3 white-noise intensity of the gyro is
uncoupled in all three sessions: correlations at most 0.09 in magnitude, on the bench, at half throttle and in
flight. The pre-registered estimate of the rate's random-walk intensity Lambda failed in all three windows:
Lambda was not positive definite. This is attributed to the rotor vibration being narrowband rather than white,
so that subtracting the short-scale noise over-corrects: at 32 ms in flight (yaw) and, on the held frames, at
all three long scales (the roll term is negative at 32, 64 and 128 ms there); the attribution was not tested. As written, the
build rule and the TH44 certificate were therefore not computed. This is a failure of the estimator, recorded as
such.

Amendment A1, exploratory (written after that run): Lambda from the 128 ms scale alone, without noise
subtraction (the white-noise share of that scale in flight is 0.2 % on roll and pitch and 10 % on yaw). On the
held frames Lambda is vibration residue, not motion, so only the flight counts. Session 7, 19 to 24 s (36 second
differences at 128 ms), from `results/step7c_A1_results.txt` (90 % moving-block bootstrap, seed 78):

- Correlations of the rate innovations: roll-pitch +0.76 (0.38 to 0.91), pitch-yaw +0.41 (0.24 to 0.64),
  roll-yaw +0.11.
- Per-axis windows against a three-axis Kalman filter: 10.6 % worse in mean square (4.1 to 18.7 %), above the
  5 % line; RMS 5.2 % (2.0 to 8.9 %). The windows this implies are 1.9, 1.9 and 13.7 ms, about two gyro samples on
  roll and pitch (the engine's own logged rate windows were 36 to 63 ms), so the continuous-time comparison is at
  its edge there.
- The symmetric matrix window T_sym = Lambda^-1 # Sn comes within 0.32 % of that Kalman filter (RMS 0.16 %).
- TH44/VM-2 certificate (generator diag(1/tau) from the engine's logged windows, 42, 36 and 63 ms; R = Sigma =
  Lambda T_env, T_env 1, 4 and 16 ms): the coupling raises chi by 7.3, 7.6 and 9.1 % (intervals from 1.8 to
  20.9 %), and the collapse direction stays on yaw, share 0.997 (0.98 to 1.00). Yaw-pure, not mixed: yaw goes first
  because its window is the slowest. The absolute chi is not interpreted, T_env being unidentified on this machine.
- Shorter windows, from `results/step7c_A1_sensitivity.txt` (seed 79; it also recomputes 19 to 24 s with its own
  draws, whose intervals differ from the above in the second digit). 19 to 23 s and 20 to 23 s: per-axis loss 7.6
  and 8.0 % in mean square (3.0 to 11.8 and 3.5 to 10.4 %; RMS 3.7 and 3.9 %); roll-pitch +0.62 (-0.06 to 0.74) and
  +0.65 (-0.25 to 0.68), intervals that cross zero; pitch-yaw +0.52 and +0.63; certificate at T_env 4 ms raised by
  6.5 % (1.2 to 54.8 %) and 14.2 % (3.8 to 43.7 %); yaw share 0.994 and 0.983, with intervals reaching down to
  0.80 and 0.86, below the 0.90 single-axis rule. The yaw-pure reading is not secure in the shorter windows.

Reading: in the one flight logged at 1 kHz, the motion appeared coupled and the noise was not. A joint estimate
is predicted to gain about 4 to 5 % in RMS on the rate path there (8 to 11 % in mean square), and the symmetric
matrix law would capture nearly all of it. This rests on 3 to 5 s of an unsteady flight that ended in a crash, on
predictions rather than filters run against a truth, and on an estimator amended after the pre-registered one
failed: a hint, not a verdict.

A point of theory (7a, D3, and the exploratory script). The V5.1 window T_sym = Q^-1 # Sigma and the steady-state
Kalman window T_K = (P Sigma^-1)^-1, P = Sigma # Q, both satisfy T Q T^T = Sigma. T_sym is the best symmetric
window, numerically: on 200 random 2x2 pairs a Nelder-Mead search over symmetric positive-definite windows,
started at T_sym, finds nothing better than 2.4e-15. The unrestricted optimum T_K is not symmetric when Q and Sigma
do not commute. In mean square, on 1000 random pairs (7a) the gap has a median of 2.75 % and a maximum of 87 %,
and follows the commutator (correlation 0.96 in logarithms); on the 200 pairs of the exploratory script, median
2.64 % and maximum 94.9 %. On the measured matrices it is small: 0.004 to 0.03 % on the attitude estimator, 0.3 to
0.4 % on the flight rates.

What it means. On this machine, in the regimes logged, the attitude estimator does not need the matrix law, and
an EKF's joint estimate of roll and pitch is predicted to be worth under 0.5 % in RMS over the law axis by axis.
On the rate path in flight a coupling of the motion appeared that the matrix law would handle; it needs a longer
steady flight at 1 kHz to be measured with tight intervals, an innovation estimator written for narrowband
vibration before the data, and magnetometer values in the log for yaw attitude. If the coupling holds, step 7b is
a three-axis rate bench, per-axis law against matrix law, with the measured correlations.

## 13. Step 6c (27 Oct 2026): ablation of the attitude-estimator law

`step6c_ablation.py` (SHA-256 d7c8b23b) and its judge `evaluate_step6c.py` (905b80fd), both pre-registered. Their
contents were fixed before the smoke run (file times) and they were committed before any result was written (the
full run had started; the hashes it printed match the commit). The smoke run (two seeds, W only) printed only the
drift fact, the oracle windows and the runtime; at V = 256 its oracle window was 10.0 s against 16.0 s predicted,
and no hypothesis was amended. Results in `results/step6c_results.txt`, rows in `results/step6c_rows.json`, run log
in `results/step6c_run.log`.

The question. A critique of step 6b asks whether the square-root law contributes anything beyond the new noise
estimator and a well-chosen fixed horizon. SYN-6B has two separable parts: the live noise estimate s2 and the map from
s2 to the window. The ablation keeps the first and replaces the second by tau = sqrt(S0_nom / lambda_hat) (s2 /
S0_nom)^p, p in {0, 0.25, 0.375, 0.5, 0.625, 0.75, 1}, every variant calibrated at the true nominal level (which
favours the alternatives); p = 0.5 is the law, p = 0 its nominal window held fixed. BEST-FIX is the best single
fixed window over a whole stream, chosen with hindsight from the oracle sweep. Scenarios: B (V = 1, 16, 1/16, step 6)
and a new W (V = 1, 256, 1/16, segments 180 s, steady part the last 120 s), where the optimal window predicted by
sqrt(S0/Q_G) moves from 1.08 to 16 s. Same generator, filters, drift fact, constants and seeds as steps 6 and 6b.
Total error J = sum over segments of the seed-mean squared RMS (a mean-square quantity).

| Row | B: V 1 / 16 / 1/16 | W: V 1 / 256 / 1/16 |
|---|---|---|
| Oracle (best fixed window per segment) | 10.8 / 20.6 / 7.1 (1.14, 3.72, 0.54 s) | 10.4 / 49.2 / 7.0 (1.05, 8.05, 0.50 s) |
| Law, POW-0.5 | 11.2 / 21.5 / 7.3 (0.87, 3.32, 0.40 s) | 10.6 / 54.9 / 7.2 (0.89, 13.3, 0.40 s) |
| BEST-FIX (hindsight) | 11.9 / 22.9 / 10.3 (2.02 s) | 21.6 / 55.3 / 21.4 (5.50 s) |
| POW-0 (nominal window fixed) | 10.9 / 28.7 / 8.0 | 10.4 / 114.1 / 8.1 |
| KF-FIX (fixed noise, EKF2 style) | 10.9 / 29.0 / 7.9 | 10.4 / 116.0 / 8.0 |
| KF-CLAIR (told the true noise) | 10.9 / 21.4 / 7.0 | 10.4 / 52.1 / 7.0 |

RMS angle error in millidegrees on the steady parts, mean over ten seeds; windows are seed medians.

J relative to the law, by exponent p:

| p | 0 | 0.25 | 0.375 | 0.5 | 0.625 | 0.75 | 1 |
|---|---|---|---|---|---|---|---|
| B | 1.546 | 1.119 | 1.019 | 1 | 1.094 | 1.362 | 2.890 |
| W | 3.796 | 1.229 | 0.833 | 1 | 2.642 | 3.903 | 5.008 |

- H6c.1 holds: POW-0.5 reproduces step 6b's rows (largest difference 4e-16 deg).
- H6c.2 holds: the law beats the best single fixed window chosen with hindsight (J ratio 0.840 in B, 0.860 in W; lower
  RMS in every segment of both, but at V = 256 by 0.7 % only, in 6 seeds of 10) and its own nominal window held fixed
  (0.647 and 0.263). It does not beat the per-segment oracle, which is a fixed window chosen with hindsight for each
  segment.
- H6c.3 holds, at the edge: in W the J-minimising exponent is p* = 0.4005 (band [0.4, 0.6]), a pass within seed
  noise; the grid minimum is at 0.375, which beats 0.5 by 17 % in J; exponents 0.25 and 0.75 cost 1.23 and 3.90 times J. In B the square root is
  the grid minimum. The oracle's own window grows with exponent 0.387 from V = 1 to 256 and 0.455 from 1 to 1/16; the
  law's with 0.501 and 0.482.
- H6c.4 fails: in W the law's window is within the band at V = 1 (0.84 of the oracle, regret 1.5 %) and 1/16 (0.81,
  3.4 %), but at V = 256 it is 1.65 times the oracle window (13.3 against 8.0 s; limit 1.5) and the regret is 10.2 %
  (limit 10 %). Per seed at V = 256: oracle window 6.5 to 19.0 s, window ratio 0.69 to 2.01 (7 of 10 above 1.5),
  regret -21 to +62 % (5 of 10 above 10 %). The window half of the failure is robust; the regret half is within seed
  noise.
- H6c.5 holds: in W the law is within 1.8, 5.4 and 4.1 % of the clairvoyant Kalman and has 53 % less error than the
  fixed-noise Kalman at V = 256 (10 seeds of 10; 9 % less at 1/16, 2 % more at nominal). Reported: against KF-OPT (one noise setting for the
  whole stream) the law has 53 % and 67 % less error at V = 1 and 1/16 and 3 % more at 256.
- Reported, not judged: transients (first 60 s of each W segment). After the drop from V = 256 to 1/16 the law's RMS is
  27.2 against 14.5 for the fixed-noise Kalman, 10.9 for the clairvoyant one and 14.8 for its own nominal window held
  fixed: the window shortens at the rate limit.
- About 95 % of W's total J comes from its V = 256 segment, so the W verdicts rest on 120 s of slow bias excursions per
  seed; everything in B is robust to the seeds.

What it means. The square-root law contributes beyond the noise estimator and beyond any fixed horizon: the exponent
does the work, and the law stays near a Kalman filter told the true noise across a x256 rise of the vibration. It
does not predict the optimal window once that window is long: at V = 256 the optimum grows more slowly than the
square root, and the law, with kappa = 1 and no refit, overshoots it by a factor 1.65. The likely cause, not tested
here, is the gyro's bias random walk. The drift fact is measured at rest over 0.128 to 1.024 s, where the angle random
walk dominates; at windows of several seconds the bias walk adds drift, and it reaches the angle through the bias
integral, whose gain is tied to the window. A stationary calculation of this loop made during the review (not in the
repository, to be checked) suggests the integral path dominates: the Allan drift at the window scale is only 5 to 21 %
above the rest value, so measuring the drift out to tens of seconds would not close the gap alone. The next test must
address the bias channel (its gain, or the bias walk entered as a second fact) as well as the scale of the drift, or
let H be measured as U' allows, pre-registered against the same rows.

## 14. Step 6d (27 Oct 2026): the mechanism of the 6c failure, and the exponent of U'

`step6d_uprime_mechanism.py` (SHA-256 c6e49019) and its judge `evaluate_step6d.py` (8491bf7c), pre-registered and
pushed to GitHub (commit a1c22e9, 19:11 UTC) before the full run. The smoke run (U75, two seeds) printed only the
fractional-noise variance check, H and lambda_H measured at rest (0.70 against 0.75 on both seeds; 0.86 and 0.91 of
the true lambda_H) and whether each oracle optimum lay inside its grid; no hypothesis was amended. Results in
`results/step6d_results.txt`, rows in `results/step6d_rows.json`, run log in `results/step6d_run.log`. Exploratory
readings after the verdict: `step6d_exploratory.py` (01c5e644), output `results/step6d_exploratory.txt`.

The questions. (M) Step 6c attributed the overshoot at x256 to the gyro's bias random walk acting through the bias
integral tied to the window; arm M tests that. (U) Theorem U' gives the exponent 1/(2H+1) for a first-order window on
a drift of index H in (0, 1); step 6c's setting lies outside it (a bias random walk, two coupled timescales). Arm U
puts the attitude estimator inside U': the gyro's drift is fractional Brownian motion of index H = 0.6 or 0.75 with
lambda_H set so that the nominal window is step 6's 1.08 s, no bias, and the filter is Mahony without bias integral,
which is exactly a first-order window on z. Scenario W of step 6c throughout (V = 1, 256, 1/16; 180 s segments, steady
part the last 120 s), same accelerometer noise and attitude truth (same random streams), seeds 600..609. Criterion
(as H6b.1 and H6c.4): window within x1.5 of the swept oracle and regret <= 10 %, seed medians, per segment.

Arm M (6b's law: tau = sqrt(s2 / lambda_hat)). REPLAY reproduces step 6c (H6d.0 holds, 4e-15 deg). M0: the step-6c
data with the bias walk switched off, Mahony with the tied integral. M1: the step-6c data, the filter given the true
bias, no integral.

| Row | Oracle (window) | Law (window) | tau / oracle | regret |
|---|---|---|---|---|
| REPLAY, V 256 | 49.2 (8.05 s) | 54.9 (13.3 s) | 1.65 | 10.2 % |
| M0, V 256 | 39.9 (19.1 s) | 48.5 (13.3 s) | 0.71 | 27.4 % |
| M1, V 256 | 37.8 (20.5 s) | 41.7 (13.3 s) | 0.67 | 3.6 % |

RMS in millidegrees, steady part, mean over ten seeds; windows are seed medians. At V = 1 and 1/16 the law passes in
all three rows. Oracle exponent from V = 1 to 256: 0.387 (REPLAY), 0.530 (M0), 0.556 (M1).

- H6d.1 fails: in M0 the oracle exponent returns to 0.530 (inside [0.45, 0.55]), but the law's regret at V = 256 is
  27.4 %.
- H6d.2 fails: in M1 the law passes the window criterion in all three segments (regret 1.6, 3.6, 2.5 %), but the
  oracle exponent is 0.556, outside [0.45, 0.55].
- The pre-registered reading rule ("if H6d.1 fails, the cause lies elsewhere: the noise estimate or the vibration
  model") does not fit these rows: M1 uses the same noise estimate and the same vibration and passes the window
  criterion. Exploratory reading (step6d_exploratory.py, not pre-registered): at V = 256 the law in M0 is 11 % worse
  than a fixed window of its own length (7 seeds of 10), while in M1 it equals it (0.998). Removing the bias walk
  restores the long optimum (19 to 20 s instead of 8 s), so the walk is what shortened it in 6c; with the bias integral
  tied to the window, the law loses about 11 % when its window has to move fifteenfold, whether or not the walk is
  there. A plausible cause, not tested: the integral keeps a bias estimate learned with the short window and corrects
  it with a time constant of about 13 tau, several minutes at these windows.

Arm U (U' with the true inputs, ULAW-K; U' with H and lambda_H measured once at rest from 300 s and s2 live, ULAW; the
square-root law of 6b on the same gyro, SQRT).

| Row | H 0.6: windows (s) | H 0.6: regret | H 0.75: windows (s) | H 0.75: regret |
|---|---|---|---|---|
| Oracle | 1.06 / 13.3 / 0.53 | | 1.05 / 10.8 / 0.58 | |
| U' predicted | 1.08 / 12.5 / 0.51 | | 1.08 / 9.33 / 0.56 | |
| ULAW-K | as predicted | 0.1 / 5.9 / 0.0 % | as predicted | 0.2 / 6.0 / 0.0 % |
| ULAW | 0.92 / 11.3 / 0.44 | 1.2 / 5.4 / 2.0 % | 0.95 / 8.40 / 0.50 | 1.2 / 4.3 / 1.6 % |
| SQRT | 1.22 / 18.2 / 0.55 | 1.8 / 11.9 / 1.0 % | 2.08 / 29.6 / 0.94 | 16.4 / 26.2 / 28.3 % |

Windows at V = 1 / 256 / 1/16, seed medians; regret, seed medians.

- H6d.3 holds: ULAW-K passes in all segments (window 0.87 to 1.02 of the oracle).
- H6d.4 holds: the oracle exponent from V = 1 to 256 is 0.484 for H = 0.6 (U' 0.455) and 0.441 for H = 0.75 (U' 0.400),
  within 0.05.
- H6d.5 holds: ULAW passes in all segments (window 0.73 to 0.89 of the oracle, regret 1.2 to 5.4 %), with H measured at
  rest at 0.54 to 0.63 (true 0.6) and 0.70 to 0.79 (true 0.75), lambda_H at 0.82 to 1.09 of the truth.
- H6d.6 holds: the square-root law fails at V = 256 for H = 0.75 (window 2.77 times the oracle, regret 26.2 %), and in
  the other two segments as well; for H = 0.6 it fails at V = 256 (regret 11.9 %).
- Reported: every oracle exponent sits above U''s value from V = 1 to 256 (by 0.03 to 0.06, including H = 1/2 in M0 and
  M1) and below it from V = 1 to 1/16 (0.446, 0.435, 0.409, 0.346 against 0.500, 0.500, 0.455, 0.400). Exploratory,
  with M1 as the H = 1/2 reference of the same class: the change of exponent from H = 1/2 to 0.75 is -0.115 (1 to 256)
  and -0.089 (1 to 1/16) against -0.100 predicted; to 0.6, -0.072 and -0.026 against -0.045.

What it means. Inside its class, U' predicts the change of exponent on the attitude estimator: when the gyro's drift
is smoother than a random walk, the optimal window grows more slowly than the square root, as 1/(2H+1) says, and the
law with H and lambda_H measured once at rest and the noise measured live lands on the optimum where the square root
does not. The step-6c failure has two parts: the gyro's bias walk shortens the optimal window at long scales, and the
bias integral tied to the window costs the law when its window moves far. For the engine on this substrate: measure H
at the scale where the window sits (U'), and give the gyro bias its own channel instead of tying its gain to the
window. Limits: one simulated axis, representative parameters, one generator family (not an independent test);
fractional Brownian drift is a stand-in for real gyro low-frequency noise (angle random walk, bias instability, rate
random walk); arm U has no bias channel; M1 gives the filter the true bias, a diagnostic, not a filter one can build;
the constant offset of the oracle exponents (about +0.04 upward, -0.05 downward) is unexplained.

