# Syntonic Engine: one formula, three facts about the machine, no tuning

**One formula, three facts about the machine, no tuning.** The syntonic engine sets the time window of a correction or an estimate from statistics measured in the signal, by a law derived from first principles (Theorem U' of the Syntony Principle). Its only machine-specific inputs are three facts measured once: the loop timing and delay; the actuator's authority in a correction loop, or, in an estimator that fuses a fast drifting sensor with a noisy absolute one, the fast sensor's drift; and access to the read noise. The law decides when to correct, not how much: it sets a window, not a loop gain, and stability stays with the machine.

Jean-Pierre Bronsard · SyntonicAI Recherche · Montréal QC Canada  
ORCID: [0009-0008-6639-7553](https://orcid.org/0009-0008-6639-7553)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23549689.svg)](https://doi.org/10.5281/zenodo.23549689)

[![MAESTRO 1 conducting Beethoven's Fifth in a MuJoCo simulation, every arm joint driven by this engine](docs/media/maestro1_conducting.gif)](https://github.com/jpbronsard/syntonic-engine/releases/download/v0.1.0/MAESTRO1_Syntonic_Conductor_v4_web.mp4)

*MAESTRO 1 conducting Beethoven's Fifth in a MuJoCo simulation, every arm joint driven by this engine. Seven seconds from the full video: [download the full video](https://github.com/jpbronsard/syntonic-engine/releases/download/v0.1.0/MAESTRO1_Syntonic_Conductor_v4_web.mp4) (7:28, 1080p, 324 MB, CC BY-SA 4.0).*

> **27 October 2026: attitude estimation at Kalman parity, with no tuning. Pending the hardware test.**
>
> Placed in a plain Mahony attitude filter, the window law matches a two-state (angle, gyro bias) Kalman filter, the one-axis linear form of the fusion an EKF performs:
>
> - **In simulation** (bench step 6b, one axis, representative parameters) its RMS error is within 2 to 4 % of a Kalman filter tuned for the test under steady vibration, and within 4 % of one told the true noise level when the vibration changes. When the vibration rises sixteenfold it has 26 % less error than a Kalman filter with fixed noise settings, the way PX4's EKF2 is configured; at nominal vibration it is 3 % behind that filter.
> - **Ablation, in simulation** (step 6c, pre-registered; the same generator as step 6b, so not an independent test): with the same noise estimate, the square root beats any single fixed window held over the whole run, even the best one chosen with hindsight (total mean-square error 16 % and 14 % lower, the second margin within seed noise). Exponents far from one half do much worse: over the run with a 256-fold rise of the vibration, the nominal window held fixed costs 3.8 times the total mean-square error and an exponent of 1 five times. The best exponent there is 0.375 to 0.4, not one half (0.375 has 17 % less total mean-square error). Through that rise the law stays within 1.8 to 5.4 % of the Kalman filter told the true noise, and at 256-fold it has 53 % less error than the fixed-noise one.
> - **On logs of our X500 on the ground** (step 7a: motor test, frame still, low throttle), the noise and drift matrices measured on roll and pitch predict the law run axis by axis within 0.3 to 0.4 % (RMS) of a joint estimate of both axes, as an EKF makes.
> - **In the one flight logged at 1 kHz** (step 7c, amendment A1, exploratory, written after the pre-registered estimator failed; 3 to 5 s of data), the motion of roll and pitch appeared coupled. The measured matrices predict the law run axis by axis 4 to 5 % (RMS) behind a joint estimate of roll, pitch and yaw, and the matrix form of the law (Syntony Principle V5.1, Theorem 35.1) within 0.2 % of it.
>
> The formula has no tuning parameter. It measures the accelerometer noise live and needs one fact measured once at rest, the gyro's angle random walk. Per axis, parity is what the theory predicts, since the steady-state Kalman filter is the Gaussian case of the syntonic window; across coupled axes the matrix law is the best symmetric window, and a joint Kalman filter can do better (bench note, section 12). The benches show that the prediction holds where the law plugs in, and that when the noise changes, a window that follows the measured noise beats one that does not.
>
> **What this is not, yet.**
>
> - The law has not run in the estimator of a real machine. The step 7 figures are predictions from matrices measured on logs, not errors of filters run against a truth.
> - The comparison is with the Kalman core, not with PX4's EKF2 software, which also estimates heading, velocity and position (magnetometer, GPS, barometer). The law does not replace that part.
> - The simulation is one axis with representative parameters. The engine's own identification failed there (step 6): the law needs the gyro's drift measured once at rest, and a drift fact overestimated by two fails (step 6b). What runs is the same law with the engine's small constants, not the same code: the noise is read by a new estimator (second differences of block means, the minimum over 32, 64 and 128 ms), the smoothness H is fixed at 1/2 by the gyro's physics instead of being fitted, the window's ceiling is 30 s instead of the servo engine's 0.5 s, and the gyro-bias correction is tied to the window with the damping of the flown firmware (1.826). With an uncalibrated gyro bias the law is slow to learn it at start, where the Kalman filter, which carries a bias state, is fast.
> - The square root sets the window too long when the vibration rises 256-fold (step 6c): 13.3 s where the best window for that level was 8.0 s (seed median; 6.5 to 19 s per seed), 1.65 times against a limit of 1.5, and a regret of 10.2 % against a limit of 10 %. That pre-registered criterion failed; its window half fails in 7 seeds of 10. Step 6d tested the attribution: without the gyro's slow bias drift, or with the true bias given to the filter, the best window returns to 19 to 20 s, so the drift shortens it through the bias estimate; but both pre-registered mechanism hypotheses failed, and why the law itself overshoots is not settled. After the vibration drops back, the law is also slow to shorten its window.
> - The flight figures rest on 3 to 5 s of one unsteady flight (throttle 0.41 to 0.56) that ended in a crash; without its last second the roll-pitch coupling's interval crosses zero.
>
> The deciding test is the gimbal bench: the law and EKF2 fed the same IMU samples, an encoder on each axis as truth, rules written before the run. If parity fails there, the failure will be published here.
>
> Details: [Results on the attitude estimator](#results-on-the-attitude-estimator-of-a-quadrotor-simulation) and [Cross-axis coupling](#cross-axis-coupling-on-logged-data-of-the-x500-step-7) below, and the bench note, sections 10 to 16.

---

## What is claimed, and what is not

On synthetic benches (`bench/software_2026-10-26/`) the engine stays within a few percent of the best fixed window and of Bayesian Kalman banks fitted with a change detector wherever the signal belongs to its class (a first-order window on a signal with a power-law structure function), tolerates increases of the sensor noise better than they do, and stays behind them where a second-order model applies (smooth targets, H >= 0.85: 25 %) or when the regime drifts fast. It does not notice when its class is the wrong one, because it computes its window instead of comparing candidates. In a simulated quadrotor attitude estimator the law, with the gyro's angle random walk entered as a fact measured at rest (step 6b), is at parity with a tuned Kalman filter under steady vibration and better than a fixed-noise one when the vibration rises; on logs of a real quadrotor on the ground, the measured matrices predict that a joint estimate of roll and pitch, as an EKF makes, would gain under 0.5 % in RMS over the law run axis by axis (step 7). In simulation, on gyro drifts drawn from the laws the theory assumes, the optimal exponent from nominal to x256 vibration follows 1/(2H+1) across H for a fractional drift (step 6e; the exponent's gain in error is resolved at H = 0.9, not at 0.75), and the local form of the law, Theorem U'-Local, fed the drift's structure function measured once at rest, predicts the optimal window and its exponent through crossovers of the drift and a jump of the optimum, from predictions frozen before the evaluation (step 6f). The law has not yet run in the estimator of a real machine.

The hypothesis under test is portability: the same law reaching the optimum of its class on machines of different physics, each described by a few facts measured once and nothing else. In its strict form, the same code with its own identification unchanged, it has already failed once, in software: on the attitude estimator of a quadrotor the engine's identification of the drift ran away, and the law reached the optimum only with the gyro's angle random walk entered as a fact measured once at rest (bench steps 6 and 6b, below). The substrates are a servo joint (ARM-1), a quadrotor on a gimbal, and, for the law rather than the code, a human in a tracking task. The pass criterion on each is the one of the software bench: window within a factor 1.5 of the swept optimum, regret under 10 %, exponent against U'. A failure on any of them is published as such and marks where the law stops. Nothing of this has run on hardware yet.

What the law does not decide is the loop gain. The window law is invariant under any rescaling of the error (the time constant does not change when the error is multiplied by a constant), so the law and the gain are separable, and the gain is outside it. On a servo joint the gain lives in the servo's own loop and the engine only moves its setpoint; that is this code. On a bare plant such as a quadrotor a proportional gain remains between the filtered error and the actuators. The SyntonicQuad flew with no PID, no I and no D, but with that gain set by hand (0.25 on roll and pitch, 0.20 on yaw, 0.90 from angle to rate) and adjusted after flights. The meta-gain law written for it (May 2026) replaces the hand-set gain by a measured control effectiveness B = d(omega dot)/du, with g = G* / B and G* a loop-gain constant declared once; it exists as a note and a single-axis simulation and has not flown. So the authority fact means, on a servo, the data-sheet speed and reaction time, and on a bare plant, B measured, with G* counted among the design constants that must carry unchanged from one machine to the next. For the quadrotor substrate the actuation-side test is stated with the window test: G* fixed on the X500 must hold on the gimbal quadrotor without retouching, or portability fails on the allocation side and is published as such.

What is common across the substrates is the window law. Where the first-order element sits depends on the substrate: here it accumulates a correction on the setpoint of a servo that stabilises itself; on a quadrotor it belongs in the attitude estimator, as the crossover of the gyro and accelerometer fusion, and not in the rate loop, where the body integrates the command and the best window is the loop-delay floor (bench step 5). Same law, different coupling, and on each machine the facts it needs, measured once: on a servo joint, the data-sheet speed and reaction time and the encoder's read noise; on the attitude estimator, the gyro's angle random walk measured at rest, with only the accelerometer noise measured live.

I developed this engine for MAESTRO 1, a two-arm robot conductor I designed. It corrects each joint toward its planned motion, and the time window it uses for that correction comes from the Syntony Principle (Theorem U'). Nothing in it is tuned per machine; the design constants it carries (scales, averaging time, warm-up, rate limit, gate) are the same everywhere, and the portability hypothesis stands or falls with them.

Every joint controller has to decide how long to average its sensor before it corrects. Average too briefly and the joint chases sensor noise. Average too long and it falls behind the motion. A PID makes that decision through gains found by trial, one joint at a time. This engine measures what the decision depends on while the robot moves (the sensor noise, and how fast and how smoothly the motion changes) and computes the window directly:

```
tau* = kappa_H * (sigma^2 / lambda_H) ^ (1 / (2H + 1)),     kappa_H = (2H * Gamma(2H + 1)) ^ (-1 / (2H + 1))
```

Here sigma^2 is the encoder noise intensity (the variance of one read times the loop period), and lambda_H and H describe how strongly and how smoothly the motion changes. With H = 1/2 the formula reduces to the classical square-root law, tau* = sqrt(sigma^2 / lambda). The window never drops below the loop delay, because a correction cannot usefully act faster than the loop can respond.

## How it works

On every loop tick, for every joint:

1. Read the encoder twice. The difference between the two reads gives the noise, and their mean gives the measurement.
2. Compare the measurement with a nominal model of the servo, driven by the plan and by the commands actually sent. What is left over is what the correction has to cancel. Measuring it this way keeps the engine from reading its own corrections back as motion.
3. Track how that residual changes across time scales (its structure function), subtract what the read noise alone would produce, keep the finest scales that stand clear of the noise, and fit its smoothness H and intensity lambda_H. Since 0.2.0 this is done on the means of adjacent blocks of 2 ms to 4 s (the Allan-variance construction), which keeps the fit valid when the sensor is noisy; see *Two identifications* below.
4. Compute the window from the formula, keep it between the floor described below and 0.5 s, and let it move smoothly.
5. Integrate the error over that window and add the result to the planned position at each servo frame.

## What you need to provide

Three facts about your hardware. There is nothing to tune. In general terms: the loop timing and its delay (which set the window floor), the actuator's authority (for a servo, its no-load speed and reaction time from the data sheet; for a bare plant, its control effectiveness, measured) or, in an estimator, the drift of the fast sensor (for a gyro, its angle random walk measured once at rest; bench step 6b), and access to the read noise (two independent reads per tick, or a known noise statistic). For a servo joint they are:

| | MAESTRO 1 |
|---|---|
| loop period `dt` | 1 ms (1 kHz) |
| servo frame period `frame` | 3 ms (333 Hz PWM) |
| servo no-load speed `w0` and nominal reaction time `servo_tau` | from the data sheet; 25 ms if it gives none |

The window floor follows from these. It is the loop delay (half a frame, plus half a tick, plus 0.1 ms for the sensor) or the discretisation limit dt / 0.45, whichever is larger. For a 1 kHz loop driving 333 Hz servos that comes to 2.22 ms. The formula leaves out the servo's own internal latency, so on real hardware measure the full loop delay once and, if it is larger, pass it as `floor`.

Each encoder has to be read twice per tick, at least one internal sample apart so that the two reads carry independent noise. An AS5048 refreshes at about 11 kHz, so the reads should be about 0.1 ms apart, which fits easily in a 1 ms tick.

## Two identifications, and where each is valid

The 0.1.x engine measured the residual's structure function on the raw reads at lags of 2 to 64 ms. That is exact where the motion over 64 ms is far above the read noise, which is the case of the MAESTRO 1 joint bench (AS5048 encoders, conducting motion). It is not the case of a noisy sensor or a slow motion: a software bench on 26 October 2026 (fractional Brownian targets, H from 0.30 to 0.95, read noise placing the optimal window from 5 to 160 ms) found that as soon as the signal falls below the noise floor at 64 ms, the raw fit has no usable lag, H collapses to its clip and the window drifts by a factor of 2 to 3. The lag at which the signal clears the floor grows from 0.1 s to 100 s over that range, so no fixed raw-read lag set can cover it.

Version 0.2.0 therefore identifies on **block means** by default (`identification='allan'`): the structure function of adjacent block means of length l, l from 2 ms to 4.096 s in powers of two. For white read noise its noise part is exactly n2 dt / l (n2 the per-read variance measured from the two reads) and is subtracted analytically; the signal part is lambda_H C(H) l^(2H) with C(H) = (2^(2H+2) - 4) / ((2H+1)(2H+2)), which is 2/3 at H = 1/2 (the Allan variance of a random walk). The fit uses the six finest scales that clear the floor, so that the scaling it measures is the local one, at the operating scale. On the software bench this lands within 7 % of the best fixed window in all 35 cells (regret at most 6.5 %), where the raw identification failed in 21 of them.

The two identifications on the MAESTRO 1 joint bench (the table below, `python examples/servo_bench.py --allan` prints both):

| Condition | Tip error, PID / raw / allan | Tip jitter, PID / raw / allan |
|---|---|---|
| A, steady noise | 0.732 / 0.734 (x1.00) / 0.749 (x1.02) | 0.651 / 0.561 (x0.86) / 0.603 (x0.93) |
| A, noise step | 1.330 / 0.922 (x0.69) / 1.033 (x0.78) | 1.261 / 0.695 (x0.55) / 0.900 (x0.71) |
| B, steady noise | 1.193 / 1.275 (x1.07) / 1.274 (x1.07) | 0.713 / 0.769 (x1.08) / 0.775 (x1.09) |
| B, noise step | 1.631 / 1.450 (x0.89) / 1.487 (x0.91) | 1.284 / 0.958 (x0.75) / 1.056 (x0.82) |

In the servo regime the raw identification remains a few percent better, mostly during the noise step, because its averages react in one second at every lag while the block means of long scales need several blocks. The default is nevertheless `allan`, because it holds on any sensor; `identification='raw'` reproduces the 0.1.x results bit for bit (the example and the replay tool use it for that reason). All three implementations carry both identifications, selected by `identification` in Python, `Config::identification` in C++ (`syntonic::ALLAN` or `syntonic::RAW`) and in Rust (`Identification::Allan` or `Identification::Raw`).

The complete software bench (scripts with pre-registered hypotheses, SHA-256 of each script, tables, per-seed rows and the bench note) is in `bench/software_2026-10-26/`; its headline beyond the fix above: against Bayesian Kalman banks the engine is at parity in stationary and slowly drifting regimes, better after sensor-noise increases, behind after noise decreases and fast drift, and a constant-velocity Kalman beats any first-order window by 25 % on smooth targets (H >= 0.85). The engine does not notice when its class is the wrong one, because it computes its window instead of comparing candidates.

## Quick start in Python

```python
from syntonic_engine import SyntonicEngine           # pip install -e python   (needs only numpy)

# q0, no_load_speeds, read_encoders, q_plan, servo_frame_due and send stand for your own code

eng = SyntonicEngine(q0, dt=0.001, frame=0.003, w0=no_load_speeds)   # radians and seconds
while running:
    a, b = read_encoders(), read_encoders()          # two reads of every joint
    cmd = eng.step(q_plan(t), a, b)                  # planned angles in, command out
    if servo_frame_due(): send(cmd)
    # eng.tau holds each joint's window (s), eng.H the measured smoothness
```

## On a Teensy or another microcontroller (C++)

`cpp/SyntonicEngine.h` is a single header with the same arithmetic in double precision, which the Teensy 4.x handles in hardware.

```cpp
#include "SyntonicEngine.h"
syntonic::Engine<6> eng;  syntonic::Config cfg;  cfg.dt = 0.001;  cfg.frame = 0.003;
eng.begin(q0, w0, cfg);
// on each tick:  const double* cmd = eng.step(q_plan, read_a, read_b);
```

A sketch to start from is in `examples/teensy/SyntonicJoint.ino`. You fill in the encoder and servo calls for your board.

## In Rust

`rust/` is a `no_std` crate with no allocation (math from `libm`), so it suits embedded Rust as well as a PC:

```rust
use syntonic_engine::{Config, Engine};
let mut eng: Engine<6> = Engine::new(&q0, &w0, Config { dt: 0.001, frame: 0.003, ..Config::default() });
// on each tick:  let cmd = eng.step(&q_plan, &read_a, &read_b, None);
```

See `rust/README.md` for details.

## Results on the MAESTRO 1 joint bench (simulation)

These results were obtained with the raw identification (`identification='raw'`, the 0.1.x engine); the example reproduces them as printed. The bench simulates the twelve arm joints of MAESTRO 1 conducting the first movement of Beethoven's Fifth. Each joint has a servo model (333 Hz command frame, 25 ms internal loop, speed derating, compliance, gear play), AS5048 encoder noise, and loads from the inverse dynamics of the planned motion. The opponent is a PID with interpolation. Its gains were picked joint by joint from a grid of 280 (P, I, D) combinations, each tried on all twelve joints over the rest of the same movement with normal noise. Both run on two held-out 30 s excerpts, once with steady noise and once with the encoder noise variance raised tenfold for 10 s. The figures below are for the baton arm and come from `examples/servo_bench.py`. A ratio below 1 means the syntonic engine did better.

| Condition | Tip error (mm), PID / syntonic | Tip jitter (mm), PID / syntonic | Joint error (deg), PID / syntonic |
|---|---|---|---|
| A, steady noise | 0.732 / 0.734 (x1.00) | 0.651 / 0.561 (x0.86) | 0.0309 / 0.0395 (x1.28) |
| A, noise step | 1.330 / 0.922 (x0.69) | 1.261 / 0.695 (x0.55) | 0.0538 / 0.0475 (x0.88) |
| B, steady noise | 1.193 / 1.275 (x1.07) | 0.713 / 0.769 (x1.08) | 0.0489 / 0.0627 (x1.28) |
| B, noise step | 1.631 / 1.450 (x0.89) | 1.284 / 0.958 (x0.75) | 0.0674 / 0.0701 (x1.04) |

Repeated over ten paired noise draws per condition (see `bench/`), the tip ratios are 1.006 [1.000, 1.013], 0.709 [0.697, 0.722], 1.063 [1.060, 1.066] and 0.889 [0.882, 0.897], with 95 % intervals.

What this shows:
- The baton tip was the bench's main score, fixed before any run. There, the syntonic engine matches the tuned PID under steady noise (a tie on one excerpt, 7 % behind on the other) and does better when the noise jumps (11 % to 31 % less error). It does this with no tuning at all.
- Scored joint by joint, the syntonic engine's error is about 28 % larger than the PID's under steady noise. That is the error the PID's gains were chosen to minimise. The syntonic engine's larger joint errors largely cancel at the tip: if they were independent, its tip error would be about 25 % higher than the PID's, and it is only 0 to 7 % higher. We do not yet know why they cancel.
- Limits: this is one simulated servo model, and its 25 ms reaction time is an assumption that happens to match the engine's nominal model exactly. It covers two excerpts of one piece, and nothing has run on real hardware yet. Drift in the servo itself (slower response, worn gears) has not been tested.

The window floor was tested on its own (`bench/`, runs v8 and v9). Set to the loop delay, it lands within 0.4 % of the best floor found by sweeping, for servo frames of 2, 3, 5 and 8 ms. A window shorter than the loop delay costs accuracy, and with slow servo frames it makes the loop ring. One case is still open: with the loop itself at 333 Hz, the rule was 7 % off the best floor (see the bench note).

## Results on the attitude estimator of a quadrotor (simulation)

The software bench's steps 6 and 6b (`bench/software_2026-10-26/`, note sections 10 and 11) place the law in the roll and pitch estimator of a quadrotor on a gimbal: one axis at 1 kHz, a gyro with angle random walk and bias drift, two IMUs whose accelerometers carry common-mode vibration with a rotor line, and a vibration level that changes with throttle (x1, x16, x1/16). The window sets the gain of a Mahony filter. The yardsticks are the best fixed window found by sweeping and a two-state Kalman filter (angle, gyro bias) tuned once, with its noise fixed at nominal as in PX4's EKF2, or switched with the true level (clairvoyant). Hypotheses and decision rules were written in each script before it ran; the parameters are representative, not measured.

| Vibration | Best fixed window | Law, drift as a fact | Kalman, fixed noise | Kalman, clairvoyant |
|---|---|---|---|---|
| nominal | 10.8 (1.14 s) | 11.2 (0.87 s) | 10.9 | 10.9 |
| x16 | 20.6 (3.72 s) | 21.5 (3.32 s) | 29.0 | 21.4 |
| x1/16 | 7.1 (0.54 s) | 7.3 (0.40 s) | 7.9 | 7.0 |

RMS angle error in millidegrees, mean over ten seeds; windows are seed medians.

What this shows:
- Given the true noise levels, the law lands on the best fixed window and at parity with a tuned Kalman filter (step 6). That is the expected equivalence of a first-order window with the steady-state Kalman gain; what it establishes is where the law plugs in.
- The engine's own identification does not carry to this substrate (step 6): the gyro-bias integral, tied to the window, erases the drift it measures, and the window ran 3 to 43 times too short. Reading the noise from the difference of two IMUs fails too, as predicted: the difference cancels the common-mode vibration, which is the noise of the accelerometer channel.
- With the gyro's angle random walk measured once at rest and only the accelerometer noise measured live (step 6b), the window lands within 2 to 4 % of the best fixed one, follows the vibration with exponents 0.505 and 0.474 (0.5 predicted), stays within 4 % of the clairvoyant Kalman, and beats the fixed-noise Kalman by 26 % when the vibration rises x16 (ten seeds out of ten). With the drift fact overestimated by a factor 2 the window falls outside the pass band; underestimated by 2, it passes.
- What runs in step 6b is the same law (kappa = 1) with the engine's window floor rule, rate limit (0.4 s/s), warm-up, starting window and block averaging, not the same code. Three parts differ from the servo engine, all fixed before the run and none fitted to its results: the drift is a fact measured at rest instead of identified live; the noise is read by an estimator written after step 6 (second differences of block means at 32, 64 and 128 ms, the minimum of the three), with H fixed at 1/2 by the gyro's physics instead of fitted; and the window's ceiling is 30 s instead of 0.5 s, since the best attitude windows run from 0.4 to 3.7 s. The gyro-bias integral is tied to the window by the damping of the flown firmware (1.826, from Kp = 2, Ki = 0.3), a constant specific to the Mahony filter.
- An ablation (step 6c, pre-registered; the same generator and seeds as step 6b, so not an independent test of the horizons) keeps the live noise estimate and replaces the square root by other exponents, each calibrated at the true nominal level. Over each run (rises of the vibration x16 and x256) the square root beats the best single fixed window chosen with hindsight (total mean-square error 16 and 14 % lower; at x256 its margin is 0.7 %, 6 seeds of 10) and its own nominal window held fixed (x256: 54.9 against 114.1 millidegrees). Over the x256 run, exponents 0 and 1 cost 3.8 and 5.0 times the total mean-square error. The law stays within 1.8 to 5.4 % of the clairvoyant Kalman and has 53 % less error than the fixed-noise one at x256 (10 seeds of 10). But over the x256 run the best exponent is 0.375 to 0.4 (0.375 beats 0.5 by 17 % in mean square; the fitted optimum, 0.4005, clears the pre-registered band [0.4, 0.6] only within seed noise), the per-segment optimal window grows with exponent 0.39 from nominal to x256, and the law's window there, 13.3 s, is 1.65 times the optimal one, 8.0 s (limit 1.5; 7 seeds of 10 above it), with a regret of 10.2 % (limit 10 %; 5 seeds of 10 above): the pre-registered horizon criterion fails at x256 and holds at nominal and at x1/16. Attributed, not tested: the gyro's bias random walk acting through the bias integral tied to the window, which the drift fact measured at rest over 0.1 to 1 s does not capture. After the drop from x256 back to low vibration the law's window shortens at the rate limit, and its first 60 s are worse than both Kalman filters and than its own nominal window held fixed (27.2 against 14.5, 10.9 and 14.8 millidegrees).
- Step 6d (pre-registered; rules pushed to GitHub a few seconds before the run; the same generator, not an independent test) asked what caused that failure and whether the generalized law U' predicts a change of exponent. Five of seven hypotheses hold, the identity check included; the two mechanism hypotheses fail. Mechanism: without the gyro's bias walk (filter unchanged), or with the true bias given and no bias integral, the best window at x256 returns from 8.0 s to 19 to 20.5 s, which supports the attribution of the shortened optimum to the walk. But the law with the tied integral and no walk loses 27.4 % there (seed median), and with the true bias it passes the window criterion only just (window 0.674 of the optimum against a 0.667 bound) while the best exponent, 0.556, falls outside its band. The pre-registered reading pointed elsewhere and does not fit; why the law itself overshoots is not settled. U': with a synthetic fractional gyro drift (H = 0.6 and 0.75, a pure power law, no bias, a first-order window), the best window's exponent from nominal to x256 is 0.484 and 0.441 against 0.455 and 0.400 predicted, inside the ±0.05 set in advance, and U' with H and lambda_H measured once at rest passes the window criterion in every segment on seed medians (regret at most 5.4 %; mean regret at x256 22 and 23 %, 4 seeds of 10 above 10 %). The test has little power: the per-seed exponents spread with a standard deviation of 0.14, their 95 % intervals contain 0.5, at H = 0.6 the square-root law is not told apart from U' (mean RMS within 1.3 %), and at H = 0.75 the square-root control fails mostly because its drift fact is miscalibrated at nominal, not because of its exponent. Consistent with U', not a confirmation of its exponent.
- Step 6e (pre-registered; rules pushed to GitHub before the run; the same generator and filter, not an independent test) is the powered version: 100 paired seeds, H = 0.5, 0.6, 0.75 and 0.9, the optimum taken on the mean error over all seeds. Four of five hypotheses hold. The optimal exponent from nominal to x256 follows 1/(2H+1): 0.509, 0.465, 0.413 and 0.369 against 0.500, 0.455, 0.400 and 0.357, slope 0.98 [0.89, 1.07], and the square root is excluded at H = 0.75 and 0.9. U' with H and lambda_H measured once at rest stays within 4.1 % of the optimum. The exponent pays in error only at H = 0.9 (20 % less than the square root at x256); at H = 0.75 the gain is not resolved, and that pre-registered hypothesis fails. The downward exponents fall below the prediction; a reading after the verdict, whose mechanism step 6f tested, attributes about 60 % of that at H = 0.5 to 0.75 (point estimates) to the accelerometer noise of this generator, whose shape changes with the vibration level; with the shape accounted for, every downward exponent lies within its 95 % interval, and the rest of the offset is not resolved.
- Step 6f (pre-registered; rules, then predictions frozen from a calibration, pushed to GitHub before the evaluation; design by the author) tests the local form of the law, Theorem U'-Local (a note whose deposit is pending), against the square root and against U' read with a raw local index (exponent 1/(2H+1), H half the log-slope of the drift's structure function at the window's scale). All three are fed the same structure function measured once at rest; the noise level changes without changing the noise's shape. Seven of seven hypotheses hold on 200 held-out seeds, the validity check included. U'-Local's exponents match the measured ones within three standard errors in all 37 secants, through a crossover of the drift from H = 1/2 to 0.9 and to a constant random bias; its window is within 2.3 % of the optimum (regret at most 0.04 %). The raw reading overestimates the exponent at every secant of the two crossover drifts, by 0.015 to 0.053, and loses to U'-Local; the square root is far off (regret up to 140 % at the highest level of the bias arm). Where the drift carries an oscillation the optimum jumps, as U'-Local predicts, by a factor 4.23 between noise levels 2 % below and above the predicted one (4.25 predicted on a fine grid, 4.19 on the evaluation's); the square root and the raw fixed-point reading predict no jump. With step 6e's noise, whose shape changes with the level, the measured exponent follows the prediction that accounts for the shape (0.4703 against 0.4715), not the one that ignores it (0.4994). The drifts are drawn from the laws the theorem assumes, so this tests the mathematics, the calibration chain and the code, not a gyroscope; the filter has no bias integral, and the noise level is given rather than measured live.
- Limits: one simulated axis, no lever arm or frame resonance, gyro noise independent of throttle, H fixed at 1/2 by the sensor rather than measured (steps 6d to 6f measured it only on synthetic drifts), and a slow start when the gyro bias is not calibrated (the Kalman carries an explicit bias state). Nothing has run on hardware; the gimbal bench, with an encoder per axis as truth, is the next test.

## Cross-axis coupling on logged data of the X500 (step 7)

Steps 6 and 6b run the law one axis at a time. An EKF estimates all axes together with a full covariance, and the matrix form of the law (Syntony Principle V5.1, Theorem 35.1, with kappa = 1: T Q T = Sigma, T = Q^-1 # Sigma, the matrix geometric mean) is proven but not in this engine. Step 7 measures, on logged data of the X500 (Pixhawk 6X, gyroscope at 1 kHz), whether the axes are coupled enough to matter before building it. It runs no filter on the logs, which carry no truth: it measures the noise and drift matrices and predicts from them the steady-state error of each window. The figures below are RMS; the bench note also gives the mean-square values, in which the build rule was written (5 %). The logs are in `bench/software_2026-10-26/data/`, reduced to the columns the scripts use; `run_step7.py` reproduces every figure.

| | Data | Law axis by axis, behind the joint Kalman (predicted RMS) | Matrix law, behind the joint Kalman |
|---|---|---|---|
| Attitude, roll and pitch (7a, pre-registered) | motor test on the ground, frame still, 55 s idle and 85 s low throttle | 0.29 and 0.41 % | 0.002 and 0.016 % |
| Rate, roll, pitch and yaw (7c amendment A1, exploratory) | flight, 19 to 24 s | 5.2 % (2.0 to 8.9) | 0.16 % |
| Same, shorter windows | flight, 19 to 23 s and 20 to 23 s | 3.7 and 3.9 % (1.5 to 5.8) | 0.17 and 0.22 % |

In brackets, 90 % moving-block bootstrap intervals.

What this shows:
- The sensor noise is not coupled across axes: gyro white-noise correlations at most 0.09 in magnitude (bench, half throttle, flight), accelerometer-tilt noise 0.09, gyro drift on the still frame 0.16. On the attitude estimator the per-axis law is enough: under the build threshold set before the run (5 % in mean square; 0.6 and 0.8 % found).
- In the one flight the motion appeared coupled: the rate changes of roll and pitch correlated +0.76 (0.38 to 0.91), pitch and yaw +0.41; without the last second the roll-pitch interval crosses zero. There, a joint estimate is predicted 4 to 5 % better in RMS (7.6 to 10.6 % in mean square, above the build threshold), and the matrix law would capture nearly all of it.
- The TH44 collapse certificate on the three rate axes, the open diagnostic of the flight record, points to yaw alone over the full window (share 0.997, 0.98 to 1.00), not to a mixed direction: yaw goes first because its window is the slowest. The coupling raises the certificate by 7 to 9 % without turning its direction. In the shorter windows the share's interval reaches down to 0.80, so the reading is not secure there.
- A point of theory: the matrix law's window is the best symmetric window (numerically, on 200 random pairs), and the joint Kalman window is not symmetric when the noise and drift matrices do not commute. On random matrices the mean-square gap has a median near 3 % and reaches 87 to 95 %; on this machine's matrices it is 0.004 to 0.4 %.
- Limits: the figures are predictions from measured matrices, not errors of filters run on the logs. The pre-registered estimator of the rate's motion failed on all three logs (attributed to the narrowband rotor vibration), so the flight figures come from an amendment written afterwards. They rest on 3 to 5 s of one unsteady flight (throttle 0.41 to 0.56, 10 to 20 A) that ended in a crash, and the per-axis rate windows they imply are about 2 ms on roll and pitch, two gyro samples, where the continuous-time comparison is at its edge. On the attitude side: low throttle, and windows of 16 to 52 s extrapolated beyond the scales measured (up to 1 s). There are no magnetometer values in the log, so yaw attitude is not covered. A steady flight of a minute or more logged at 1 kHz, with an estimator for narrowband vibration written before the data, settles it; if the coupling holds, a three-axis bench of the matrix law follows.

## Checks

- The Python engine in `python/syntonic_engine` with `identification='raw'` reproduces the bench engine exactly: zero difference in joint angles and windows over 33 s on 12 joints; a unit test pins it to the 0.1.1 numbers.
- On the same 33 s, 12-joint replay, the C++ header agrees with the Python engine to 1.5e-14 rad (raw identification) and 7.8e-16 rad (block means); the Rust crate to 5.4e-14 rad and 8.9e-16 rad.
- The block-mean identification is unit-tested on a fractional Brownian target with read noise large enough to defeat the raw one (`python/tests`).
- `examples/arm_fk.py`, the baton-tip kinematics in plain numpy, agrees with MuJoCo to 7e-16 m over 200 random poses.
- Neither the C++ header nor the Rust crate has run on a microcontroller yet. The Rust crate builds as `no_std` but has not been cross-compiled for a Cortex-M board, and the Teensy sketch has only been compiled on a PC against stand-ins for the Arduino functions.

To run the checks yourself, from the repository root (Python 3.9 or later, a C++11 compiler, and a recent stable Rust):

```
pip install -e python scipy
python -m unittest discover -s python/tests          # unit tests of the Python engine
python examples/servo_bench.py                       # the results table above, about a minute; add --allan for both identifications
python rust/tools/make_replay_log.py                 # writes replay_log.bin (raw), replay_log_allan.bin and replay_w0.bin
g++ -O2 -std=c++11 -Icpp cpp/tests/replay.cpp -o replay_cpp
./replay_cpp replay_log.bin replay_w0.bin raw && ./replay_cpp replay_log_allan.bin replay_w0.bin allan
cd rust && cargo test && cargo run --release --example replay -- ../replay_log.bin ../replay_w0.bin raw \
        && cargo run --release --example replay -- ../replay_log_allan.bin ../replay_w0.bin allan
```

## What is in the repository

```
python/syntonic_engine/   the engine (numpy), with unit tests in python/tests/
cpp/SyntonicEngine.h      the engine as a single C++ header (both identifications), with a replay test in cpp/tests/
rust/                     the engine as a no_std Rust crate, with a replay check against Python
examples/                 servo_bench.py (the table above, numpy + scipy), arm_fk.py, data/, teensy/
bench/                    the MAESTRO 1 joint bench scripts and every result file (the scripts need the MAESTRO 1
                          simulation model, which is not included, but all results are), and the software bench
                          of 26 Oct 2026 against adaptive estimators, with the quadrotor steps 5 to 6f and
                          step 7 on logged X500 data (software_2026-10-26/, self-contained, logs in data/)
docs/bench_note.md        the bench note: setup, the predictions written before each run, results and failures
```

## Citing

If you use this work, please credit me (see NOTICE and CITATION.cff):

- this software: J.-P. Bronsard, *Syntonic Engine*, Zenodo, https://doi.org/10.5281/zenodo.23549689 (this DOI always points to the latest version; each release also has its own DOI on Zenodo)
- the theory: J.-P. Bronsard, *The Syntony Principle*, Zenodo, https://doi.org/10.5281/zenodo.17254395

## Licence

Apache License 2.0 (see LICENSE). Any redistribution or derivative work must keep the NOTICE file.

## Authorship

Code and documentation were drafted with language models under the author's direction, as solution architecture and integration of independent components. The author alone approves every statement and is accountable for it; no language model is an author.

Jean-Pierre Bronsard, SyntonicAI Recherche, 2026.
