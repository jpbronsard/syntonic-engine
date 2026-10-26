# Software bench of the syntonic engine against adaptive filters

Bench note, version 1.1. 26 October 2026. SyntonicAI Recherche.
Scripts, tables and per-seed rows: `bench/software_2026-10-26/` of the syntonic-engine repository (the
bank checkpoints, 110 MB each, are recomputed by the scripts).
Version 1.1 adds section 8: the engine rows re-run with the identification exactly as released in
syntonic-engine 0.2.0 (commits 17a1536 and a9cc957 on GitHub), which fits the six finest scales that clear the
noise floor where the bench scripts fitted all of them. No conclusion changes; the numbers that move are listed.

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

