# The syntonic laws against tuned PIDs and tuned Kalman filters, pre-registered (simulation), 2 Nov 2026

Request from the author: "trouve comment la loi peut à la fois être à parité ou dépasser un PID réglé et également égaler ou dépasser un Kalman réglé en toute circonstances. Servo 333 Hz max. Test réaliste." Written and run by the assistant on 2 Nov 2026, with two sub-agents (the realistic servo bench, and the exploration harness of the flight part), not yet reproduced by the author. Simulation only.

## In short

- **Servo (confirmation 1, held).** The syntonic loop law computes a joint's PI gains from statistics it measures live; nothing is tuned. On a realistic MAESTRO 1 bench with 333 Hz servo frames, it is at parity with or better than PIDs tuned at commissioning in all nine circumstances, on fresh seeds. S1, S2 and S3 held.
- **Flight, looped bench (confirmation 2, failed, on record).** The steady-state law lost parity with the corpus's tuned Kalman under a vibration step (1.149), and lost by about 2x to a Kalman scheduled on the throttle level in the switching flights, as predicted. F1, F5 and F8 failed.
- **Flight, realistic bench (confirmation 3, held).** Two amendments to the bench were declared after confirmation 2 failed: block-bootstrap replay of the logged noise instead of looping it, and an unknown initial attitude. On that bench the Riccati form of the law (LAW_R) has a lower error than all three truth-tuned Kalman filters in all six circumstances, on fresh seeds. R1, R2 and R3 held. The steady-state law of confirmation 2 still loses to the scheduled Kalman in switching flights (1.76 and 1.68), as predicted (R3): what wins is the memory the Riccati recursion carries.
- **In the package since 0.3.0 (Python only):** `syntonic_engine.SyntonicLoopLaw` and `syntonic_engine.RiccatiAttitudeLaw`. `package_gate.py` checks that they reproduce the confirmations exactly.

## The sequence (times UTC, file times of the working session)

1. **Exploration, not pre-registered.** The servo law was developed on seeds 91-93 of the published plant, then on seeds 91-92 of the realistic bench (`servo/dev_pwm.py`). A sub-agent's flight harness (`exploration_harness/`) ran on dev seeds 300-309.
2. **Confirmation 1, servo.** `servo/confirm_servo.py` was frozen at 14:15:44 and run from 14:16 to 14:44. S1-S3 held.
3. **Flight development.** This step is in `flight_looped/flight_dev/`. It added a stronger rival, the throttle-scheduled Kalman. `flight_looped/flight_snap/confirm_flight.py` was frozen at 16:02:39 and run from 16:03 to 16:12:52 on the looped bench. F1, F5 and F8 failed.
4. **After the failure.** The author said "poursuit". The work after that point (first file 16:24) has four parts:
   - amendment 3, the block-bootstrap replay;
   - amendment 4, the unknown initial attitude;
   - the rivals re-tuned on the amended bench;
   - the Riccati form developed on dev seeds 300-309: `ricv.py`, then `ricv2.py` with two options, then the bound on q_b.
5. **Confirmation 3, flight.** `flight_realistic/confirm_flight3.py` was frozen at 17:24:55 and run from 17:25 to 17:36 on fresh seeds 500-509. R1-R3 held.

Confirmation 2's failure stands. The realistic bench and the law that passed on it were both chosen after seeing that failure. Only the seeds of confirmation 3 were fresh.

## Servo: the syntonic loop law against tuned PIDs (confirmation 1)

**The law** (`servo/syntonic_loop_law.py`, sha256 a7a234cc...). A joint correction whose PI gains are computed, every 0.5 s, as the gains that minimise the loop error predicted from statistics measured live:

- the disturbance spectrum, from second differences of block means of the loop residual at lags of 2 to 256 ms;
- the read noise, from the two encoder reads.

The candidate gains are restricted to those that keep the loop robust: phase margin at least 30 degrees, and max |S| at most 2.

- **Facts.** dt, the frame, servo_tau and w0 from the data sheet, plus the servo's internal delay measured once at commissioning (2.5 ms).
- **Shared with the primary rival.** The lead (27.1 ms x planned velocity) and the frame-mean proportional term.
- **Design constants, declared.** Averaging 1 s, warm-up 1 s, noise gate 2, at most 3 fitted lags, update every 0.5 s.

**Bench** (`servo/bench/realistic_bench.py`, written by a sub-agent, sha256 1fe422ea...). It models the baton arm RJ1-RJ6 of MAESTRO 1:

- a 1 kHz loop, with 3 ms frames (333 Hz) held between sends;
- a first-order 25 ms servo with a load-derated speed limit;
- 1 degree of compliance at stall and 0.3 degree of gear play;
- the AS5048 encoder read twice per tick;
- command quantization of 0.025 degree;
- a PWM latency of 2.5 ms.

Gate G1 checks that, at the published plant, the bench reproduces `examples/servo_bench.py` bit for bit. It held exactly.

**Rivals** (`servo/bench/rivals.py`, sha256 e4bf3df3...). Both were tuned at commissioning on the published 280-point grid with the edge rule, on tuning seeds 61-63.

- **PIDLM** is the PID with the same lead and the same frame-mean proportional term as the law. It is the primary rival.
- **PIDL** is the PID with the lead only, the conventional tuned PID.

**Circumstances in operation.** Nothing is re-tuned and no new fact is given.

- C1: nominal.
- C2: noise step, read variance x10 from 10 to 20 s.
- C3: bursts, x30 for 0.5 s every 2.5 s.
- C4: latency 3.5 ms.
- C5: latency 1.5 ms.
- C6: servo 35 ms.
- C7: servo 18 ms.
- C8: load x1.2.
- C9: command quantization 0.09 degree.

**Runs and score.** Fresh seeds 71-75 x excerpts A and B, so 10 paired runs per circumstance. Score: baton-tip rms (mm), pooled over the 10 runs.

| Pooled tip rms (mm) | PIDLM | PIDL | LAW | LAW / PIDLM | LAW / PIDL | law better than PIDLM |
|---|---|---|---|---|---|---|
| C1 nominal | 0.491 | 0.565 | 0.488 | 0.994 | 0.864 | 5 of 10 |
| C2 noise step x10 | 0.839 | 0.922 | 0.653 | 0.778 | 0.708 | 10 of 10 |
| C3 bursts x30 | 0.916 | 1.008 | 0.813 | 0.887 | 0.806 | 10 of 10 |
| C4 latency 3.5 ms | 0.541 | 0.623 | 0.526 | 0.971 | 0.843 | 10 of 10 |
| C5 latency 1.5 ms | 0.453 | 0.522 | 0.457 | 1.009 | 0.876 | 5 of 10 |
| C6 servo 35 ms | 0.566 | 0.714 | 0.526 | 0.930 | 0.737 | 10 of 10 |
| C7 servo 18 ms | 0.568 | 0.642 | 0.551 | 0.970 | 0.858 | 10 of 10 |
| C8 load x1.2 | 1.029 | 1.081 | 0.976 | 0.948 | 0.902 | 10 of 10 |
| C9 quantization 0.09 deg | 0.507 | 0.583 | 0.500 | 0.985 | 0.858 | 8 of 10 |

| Prediction | Bar | Result | Held? |
|---|---|---|---|
| S1 parity with PIDLM in every circumstance | LAW / PIDLM <= 1.10 in C1-C9 | worst 1.009 (C5) | held |
| S2 better where the law re-measures what changed | LAW / PIDLM <= 1.00 in C2, C3, C6, C7, C8 | 0.778 to 0.970 | held |
| S3 better than the conventional tuned PID everywhere | LAW / PIDL <= 1.00 in C1-C9 | worst 0.902 (C8) | held |

Not claimed beforehand: better in C1 (the PIDs' tuning point), C4 and C5 (the delay fact is wrong for both), C9.

**Described, not judged.**

- **Per excerpt.** Excerpt A is the law's weaker one. LAW / PIDLM is 1.023 on A and 0.973 on B in C1, 1.035 and 0.990 in C5, 1.002 and 0.973 in C9.
- **Jitter.** The law has less jitter (tip rms above 5 Hz) than PIDLM in every circumstance. The ratios run from 0.651 (C2) to 0.926 (C6), and 0.873 at nominal.
- **Limits of the battery.** Load x1.5 exceeds the stall torque and was dropped as infeasible. A resonant second-order servo destabilises every controller on this bench. An unmodeled dead time breaks the law unless the delay is measured once at commissioning, as the README asks for the full loop delay; the PIDs receive the same fact through their tuning.

## Flight: the attitude laws against tuned Kalman filters

**What every flight bench here shares.** It is the repository's flight bench (`bench/flight_2026-11-02/`):

- one axis, with the plant at rest, so the estimate is the estimation error;
- the X500's logged noise: S6 idle and low, and S5 half throttle (the only half-throttle record, 7.96 s);
- the gyro at 1 kHz and the accelerometer tilt at 250 Hz;
- a 5 s idle pre-roll on the ground, with the arming calibration as the pre-roll gyro mean.

**Circumstances.** Each is 60 s after the pre-roll.

- C1: steady low throttle.
- C2: steady half throttle.
- C3: half/low every 15 s.
- C4: vibration step, tilt noise x sqrt(10) from 20 to 40 s, at low throttle.
- C5: half/low every 5 s.
- C6: gyro bias step +0.02 deg/s at 30 s, at low throttle.

**Rivals.** All the Kalman rivals are tuned against the true angle, an oracle no real filter has, on training seeds 100-102 and the C3 schedule. They are not re-tuned per circumstance. The laws are not tuned against truth.

### Confirmation 2, looped bench (failed)

**Pre-registration.** `flight_looped/flight_snap/confirm_flight.py`, sha256 6eac13e14a5a0755bd1a6b773baf5bcf8888a4d95e0aa316819a1af2dcbcbaea. 14 hashed files. Seeds 400-409.

**Noise.** Amendment 2 of the flight bench: the records are looped from a random offset.

**Rivals.**

- KF_rival: the corpus's two-state Kalman, with fixed R, Q_G and Q_B.
- KF_sched: the same filter with R scheduled on the throttle level.
- KF_rival_btab and KF_sched_btab: the same two filters given the per-throttle bias table.

**Laws.**

- LAW_I: b4_nm_ic2. Steady-state Kalman-Bucy gains from the live noise reading, with per-throttle memory.
- LAW_II: a Riccati recursion of (theta, b_idle, b_low, b_half) with the bias table.

| Pooled rms error (deg) | C1 | C2 | C3 | C4 | C5 | C6 |
|---|---|---|---|---|---|---|
| KF_rival | 0.0699 | 0.2902 | 0.2794 | 0.0756 | 0.3401 | 0.2083 |
| KF_sched | 0.0609 | 0.2066 | 0.0994 | 0.1000 | 0.1029 | 0.1280 |
| KF_rival_btab | 0.0654 | 0.2537 | 0.2695 | 0.0694 | 0.3425 | 0.1885 |
| KF_sched_btab | 0.0649 | 0.1424 | 0.1007 | 0.0697 | 0.0846 | 0.1809 |
| LAW_I | 0.0609 | 0.2560 | 0.1959 | 0.0868 | 0.2169 | 0.0996 |
| LAW_II | 0.0696 | 0.1915 | 0.1370 | 0.0871 | 0.1259 | 0.1149 |
| LAW_I / KF_rival | 0.872 | 0.882 | 0.701 | **1.149** | 0.638 | 0.478 |
| LAW_I / KF_sched | 1.001 | 1.239 | **1.972** | 0.868 | **2.108** | 0.778 |
| LAW_II / KF_rival_btab | 1.065 | 0.755 | 0.509 | **1.256** | 0.368 | 0.609 |
| LAW_II / KF_sched_btab | 1.073 | 1.344 | 1.361 | 1.251 | 1.488 | 0.635 |

**Verdicts.**

- F1 failed: C4 is 1.149, above 1.10.
- F2 held.
- F3 held.
- F4 held as predicted: the law is not at parity with KF_sched in C3 and C5, so the goal is not met against that rival.
- F5 failed: C4 is 1.256.
- F6 held.
- F7 held as predicted: the goal is not met against KF_sched_btab in C2, C3 and C5.
- F8 failed: C4 is 1.251.

**Common cause.** The steady-state law has no memory of the information it has accumulated. The scheduled Kalman carries the precision gathered at low throttle through the noisy phases.

Details: `flight_looped/VERDICT_flight.md`, `flight_looped/NOTES_flight_dev.md`.

### Amendments 3 and 4 (declared after confirmation 2 failed)

**Amendment 3: block-bootstrap replay** (`flight_realistic/bench_bb.py`, sha256 9d36615a...).

- **What was wrong.** The half-throttle record is 7.96 s long, so a looped 60 s flight replays it 7.5 times. The tilt noise then sums to exactly zero over every 7.96 s, and the gyro integral never random-walks beyond it. The truth-tuned filters exploited this: they set Q_G = 0 and Q_B = 0, and with the bias table they ignored the half-throttle tilt altogether.
- **The change.** Each throttle segment is filled with consecutive blocks drawn at random offsets of its record. The block length is T_lim = T_rec / 4, the longest scale at which the record's own noise facts are measured.
- **Validation on 60 s windows.**
  - At half throttle, looping made the tilt mean 5.7 times too precise against white-noise theory (8.4 times against the bootstrap), and the gyro integral 8.6 times too small against random-walk theory.
  - At low throttle, the looped tilt mean was 4.6 times too precise.
  - The bootstrap sits within 1.0 to 1.5 times of theory.

**Amendment 4: unknown initial attitude.**

- **The change.** The plant rests at theta0 ~ N(0, 0.5 deg), drawn per flight. Every estimator starts at 0 and is scored on the estimate minus theta0.
- **Why.** Without it, a filter that ignores the tilt and starts at the true angle won the low-throttle circumstances by never moving. No real flight controller knows its initial attitude.

**Rivals re-tuned on the amended bench** (`flight_realistic/rivals_bb.py`, then `results/rivals_r4_*.json`). The grids were widened and recentred until every parameter was interior, at its physical limit 0, or on a declared flat ridge.

- **KF_rival3.** Fixed R, Q_G and Q_B. R ran up a flat ridge; J_train was 0.426247, constant to 6 digits.
- **KF_sched3.** R = m rho s2_rec / DTA, scheduled on the throttle level. m and Q_B ran down a flat ridge; J_train 0.151557.
- **KF_full3.** R and Q_G both scheduled. It converged to the same filter as KF_sched3 (Q_G = 0).

### Confirmation 3, realistic bench (held)

**Pre-registration.** `flight_realistic/confirm_flight3.py`, sha256 181406c194111be47720d622e472d7f4e9cc8372dae9621a6004de25b71f0066. 20 hashed files. Fresh seeds 500-509.

**Gates.**

- G0: the hashes.
- G1: seed 300, C1 and C4. Every column reproduced the dev runs exactly.
- G2: the rival edges are interior, at 0, or on the declared flat ridges.

**LAW_R** is `ricv2.RicV` with dict(k=3, qt='level', prior='vb', qb = q_b at T_rec/2, sa=True, rm='max'). Nothing in it is tuned against truth.

- **State.** (theta, b_idle, b_low, b_half). The active level's bias drives the gyro integration and random-walks at q_b. q_b = 3 sigma_A^2(T) / T at T = T_rec / 2 of the idle record, 9.23e-11 rad^2/s^3, the tighter of the two published bounds.
- **Angle random walk.** q_theta(level) = lam(level), the record's angle random walk at that throttle.
- **Measurement.** z is the tilt, with R = max(s2_rec(level), live reading) / DTA. The live reading is the published Hadamard pipeline, with one memory per level initialised at the fact.
- **Stationary alignment in the pre-roll.** The vehicle is on the ground while its gyro is calibrated, so theta is estimated from the tilt alone and the gyro is not integrated.
- **Arming.** Biases 0. Bias covariance (lam_idle / T_pre) x ones, plus lam(level) / T_rec(level) on the low and half diagonal.
- **How the structure was chosen.** On dev seeds 300-309, in sequence: `ricv.py`, then `ricv2.py` with the options sa and rm, then the q_b bound (`dev_r4*.jsonl`). The choice was declared post hoc in the docstring.

**Described, not judged.**

- LAW_I: the steady-state law of confirmation 2.
- LAW_Rqb4: LAW_R with q_b at T_rec/4.

| Pooled rms error (deg) | C1 | C2 | C3 | C4 | C5 | C6 |
|---|---|---|---|---|---|---|
| KF_rival3 | 0.0702 | 0.5638 | 0.5209 | 0.0866 | 0.4429 | 0.2133 |
| KF_sched3 (= KF_full3) | 0.1102 | 0.5963 | 0.1728 | 0.1716 | 0.1519 | 0.1387 |
| **LAW_R** | **0.0636** | **0.4064** | **0.1387** | **0.0853** | **0.1258** | **0.1192** |
| LAW_I (described) | 0.0789 | 0.4638 | 0.3045 | 0.1347 | 0.2552 | 0.1111 |
| LAW_Rqb4 (described) | 0.0661 | 0.4108 | 0.1388 | 0.0942 | 0.1260 | 0.0977 |
| LAW_R / KF_rival3 | 0.906 | 0.721 | 0.266 | 0.985 | 0.284 | 0.559 |
| LAW_R / KF_sched3 | 0.577 | 0.681 | 0.803 | 0.497 | 0.829 | 0.859 |
| LAW_I / KF_sched3 | 0.716 | 0.778 | **1.762** | 0.785 | **1.681** | 0.801 |
| LAW_R better than KF_rival3 (of 10) | 7 | 7 | 10 | 6 | 10 | 10 |
| LAW_R better than KF_sched3 (of 10) | 9 | 10 | 8 | 10 | 9 | 8 |

| Prediction | Bar | Result | Held? |
|---|---|---|---|
| R1 the goal: parity or better with every tuned Kalman in every circumstance | LAW_R / KF <= 1.10 for the three KF, C1-C6 | worst 0.985 (C4, KF_rival3) | held |
| R2 better | LAW_R / KF_rival3 <= 1.00 in C1, C2, C3, C5, C6 and LAW_R / KF_sched3 <= 1.00 in C1, C2, C4, C5, C6 | 0.266 to 0.906 and 0.497 to 0.859 | held |
| R3 the memory is what matters | LAW_I / KF_sched3 > 1.10 in C3 and C5 | 1.762 and 1.681 | held |

Dev figures (seeds 300-309) written in the script before the run:

- R1: dev max 0.891, all bootstrap upper bounds at most 1.05.
- R2: 0.85, 0.60, 0.25, 0.23, 0.56 and 0.55, 0.76, 0.36, 0.82, 0.89.
- Not claimed beforehand: C4 against KF_rival3 (0.89 [0.76, 1.05]) and C3 against KF_sched3 (0.88 [0.74, 1.04]). On the test seeds both came out below 1 (0.985 and 0.803).
- R3: 2.17 and 2.24.

**Described.**

- LAW_I loses to the fixed-R KF_rival3 in C1 (1.124) and C4 (1.554).
- LAW_Rqb4, with the faster bias walk, does better than LAW_R in the bias step C6 (0.0977 against 0.1192) and worse in C4 (1.087 against KF_rival3).
- On this bench the scheduled filter, tuned on the switching schedule, is worse than the fixed one at steady throttle (C1, C2, C4).

## What it means, in plain terms

- **Servo.** On this bench the goal is met. The law matches or beats PIDs tuned at commissioning, with nothing tuned, and does better where conditions change after commissioning, because it re-measures them.
- **Flight.** The goal is met on the realistic bench, by a law with memory: the Riccati recursion of the model the facts describe. The steady-state law, which only follows the live noise reading, is not enough when the throttle switches.
- **What the looped bench showed.** It flattered the long memories of the truth-tuned filters. That is why the bench was amended, but the amendment came after the failure.

## Limits

- **Simulation only.** No hardware ran.
  - The servo bench carries the imperfections listed above, but the servo is exactly first order at the time constant the law assumes, except in C6 and C7.
  - The flight bench is one axis with a linear tilt model and the plant at rest. Nothing flies; it scores estimation only.
- **The flight bench and the law that passed on it were chosen after confirmation 2 failed.**
  - The two amendments are argued from the records' own facts and validated against them. But they were declared after the failure, and they changed the ranking.
  - The law's structure was chosen on dev seeds that had already been run on the amended bench.
  - Only the test seeds 500-509 were fresh.
- **The rivals are tuned against truth.** That is an oracle no real filter has, and here it favours the rivals. They were tuned once on the switching schedule and not re-tuned per circumstance. KF_sched3 and KF_full3 converged to the same filter, so "three Kalman filters" are two distinct ones.
- **The noise comes from three short records of one airframe** (35, 40 and 8 s); the half-throttle record is 8 s. The decisive test is a real log of 60 s or more at half throttle, then the gimbal bench with encoders as truth.
- **The laws are in the Python package only.** They are not yet in the Rust crate or the C++ header. RiccatiAttitudeLaw is an estimator, not a flight controller.

## In the package (0.3.0) and the gate

- **`syntonic_engine.SyntonicLoopLaw`** (`python/syntonic_engine/loop_law.py`). It is `servo/syntonic_loop_law.py` with a provenance paragraph added to its docstring; the code is unchanged.
- **`syntonic_engine.RiccatiAttitudeLaw`** (`python/syntonic_engine/attitude_law.py`). It is LAW_R rewritten as a self-contained single column:
  - its facts are constructor arguments: per level s2, lam and T_rec; q_b; T_pre; dt and dt_accel;
  - the noise pipeline's constants are those of the harness: block scales 8, 16 and 32 samples, EMA 1 s, a second difference used once 4 blocks have completed after a level change;
  - its `step(g, a, reg_lam, t, level, k)` and `arm(bias_cal)` are the harness's.
- **`package_gate.py` (a check, not a prediction).** Run inside the frozen harnesses, the package classes reproduce the recorded runs, with a difference of exactly 0.0.
  - By default, servo run 'C1 A 71' (all three columns, the law from the package; LAW 0.46094973372876746 mm), and LAW_R and LAW_Rqb4 for seeds 500 and 501 in C1 and C3.
  - With `all`, also the LAW column of all 90 servo runs and LAW_R and LAW_Rqb4 in all 60 flights of confirmation 3.
  - On 2 Nov 2026 both modes gave a largest difference of 0.0. For one flight (C5, seed 500), the 60,000 per-tick estimates were identical to those of `ricv2.RicV`.
- **Unit tests** in `python/tests/test_loop_law.py` and `test_attitude_law.py` (numpy only).

## Files

| Folder | What | sha256 of the frozen scripts |
|---|---|---|
| `servo/` | `confirm_servo.py` (pre-registration), `syntonic_loop_law.py` (the law), `bench/realistic_bench.py`, `bench/rivals.py` (sub-agent), `results_servo/` (gates, tuning, the 90 runs, summary, chain logs), `VERDICT_servo.md`, `dev_pwm.py` with its tuning `dev_tune_pwm.json` and runs `dev_pwm.jsonl` | af4996a0efaae1d8bcab93aeda59bc1ce8cd2b957dceb444034a173662f6f567, a7a234cc..., 1fe422ea..., e4bf3df3... |
| `flight_looped/` | confirmation 2 (failed): `flight_snap/` (`confirm_flight.py`, `law_flight.py`, the frozen snapshot of the exploration harness, `results/` facts and tuned rivals, `results_confirm/`), `flight_dev/` (the development before it), `VERDICT_flight.md`, `NOTES_flight_dev.md` | 6eac13e14a5a0755bd1a6b773baf5bcf8888a4d95e0aa316819a1af2dcbcbaea, e0e9fbaa... |
| `flight_realistic/` | confirmation 3 (held): `confirm_flight3.py`, `bench_bb.py` (amendments 3-4), `rivals_bb.py`, `ricv.py`, `ricv2.py` (the law), the dev scripts and their `dev_*_runs.jsonl`, `tab4.py`, `results/` (tuned rivals), `results_confirm3/`, `agent/` (two diagnostics of the sub-agent; no constant from them enters the law). The timing logs of the dev passes are not kept | 181406c194111be47720d622e472d7f4e9cc8372dae9621a6004de25b71f0066, 9d36615a..., 3e68f315..., 3a32772c... |
| `exploration_harness/` | the sub-agent's flight harness, exploratory, dev seeds 300-309 (origin of `common.py` and `variants*.py`), with its results | |
| `package_gate.py`, `replay_setup.sh` | the gate above, and the replay helper below | |

The files of the working session are copied byte for byte (128 files checked against their sources); none was edited. The full hashes are in the scripts' headers and `results_*/frozen.sha256`.

## How to replay

These runs used Python 3.13, numpy 2.5.3 and scipy 1.18.1. Other versions may differ in the last digits.

**Servo.** `confirm_servo.py` finds its files relative to its own folder, and the repository through the environment variable SYNTONIC_REPO (default: the working session's path):

```
cd bench/goal_2026-11-02/servo
SYNTONIC_REPO=../../.. python3 confirm_servo.py gate        # G0 and G1; rewrites results_servo/gates.json
SYNTONIC_REPO=../../.. python3 confirm_servo.py summarize   # the verdicts from the recorded runs
```

`run C1` appends to `results_servo/runs_C1.jsonl` (about 30 s per run, 10 runs per circumstance), so run it on a copy of the folder.

**Flight.** The frozen flight scripts name absolute paths of the working session: `/home/claude/repo/syntonic-engine` and `/home/claude/work/goal/...`. Editing them would break their hashes. `replay_setup.sh` creates those paths as symbolic links to this folder (or to a copy of it) and to the repository. It needs write access to `/home/claude` (as root, or in a container). It refuses to replace anything that is not a link, and `--remove` takes the links away.

```
sudo sh bench/goal_2026-11-02/replay_setup.sh                 # or, to rerun the run modes, a copy:
cp -r bench/goal_2026-11-02 /tmp/goal && sudo sh bench/goal_2026-11-02/replay_setup.sh /tmp/goal
python3 bench/goal_2026-11-02/package_gate.py                 # about 30 s; 'all' about 20 minutes
cd /home/claude/work/goal/r3 && python3 confirm_flight3.py gate        # G0 (20 hashes), G1, G2; about a minute
python3 confirm_flight3.py summarize                                    # the verdicts from the recorded runs
cd /home/claude/work/goal/confirm/flight_snap && python3 confirm_flight.py gate   # confirmation 2's gates
sudo sh bench/goal_2026-11-02/replay_setup.sh --remove
```

`gate` and `summarize` rewrite `gates.json` and `summary.json`, with the same content when the run is reproduced. The `run` modes skip seeds already recorded. To rerun the 60 flights of a confirmation, move its `runs_*.jsonl` aside in the copy; each flight takes about 20 s. These replay steps were checked on 2 Nov 2026 in an isolated mount namespace, with a copy of the repository outside `/home/claude`. The rewritten `gates.json` and `summary.json` of all three confirmations were identical to the recorded ones, and so was `runs_C1.jsonl` of confirmation 3 after its 10 flights were run again.
