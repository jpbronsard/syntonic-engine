# A syntonic flight controller against a PID and a PID with a tuned Kalman filter (one axis, the X500's logged noise), 2 Nov 2026

Request from the author ("monte-moi le meilleur contrôleur de vol syntonique avec toutes nos connaissances, k = 1, et compare-le à un PID dans les mêmes conditions et à un Kalman ajusté"); scope he chose: one axis, the X500's real noise, three full controllers; then "peut-être k dynamique ici". Written and run by the assistant on 2 Nov 2026, not yet reproduced by the author. Simulation only.

Pre-registration: the predictions are in the docstring of `flight_bench.py`, hashed and frozen (sha256 39e240d0...) before any run except a smoke test. Four amendments followed, each written and hashed before the step it changes and before any test-seed run (table below). The two `explore_*.py` scripts ran after the verdicts and are not judged. The scripts need `bench/software_2026-10-26` (the step 6g loader and the X500 logs in its `data/`).

## Bench

- **Plant (roll axis):** motor lag 30 ms, B = 4000 deg/s² per unit, damping 2 /s, command limit ±0.25, 1 kHz loop.
- **Disturbance:** Ornstein-Uhlenbeck torque (60 deg/s², 1 s) plus 4 gusts of ±300 deg/s² for 0.25 s.
- **Flight:** 5 s on the ground (arming gyro calibration at idle), then 60 s of flight. Throttle alternates half / low / half / low every 15 s. Square manoeuvres of ±10° every 2.5 s between 10 and 50 s, shaped by a prefilter at 10 rad/s.
- **Sensors:** the X500's logged noise. Gyro at 1 kHz; accelerometer tilt at 250 Hz. Records used: S6-idle, S6-low (40 s) and S5-half (8 s, looped; it is the only half-throttle record).
- **Seeds:** training 100-102, test 200-209 (never used in tuning).
- **Controllers:** all use the same cascade, angle P, then rate PID, then the motor. Each arm's gains are tuned on the training seeds over the same grid.

| Arm | Attitude estimator | Feed-forward |
|---|---|---|
| PID | Mahony, fixed 0.5 s window of the flown firmware, its bias rule | no |
| KF | two-state Kalman (angle, bias), R and Q_G tuned against truth on the training seeds | no |
| SYN1 | the law at kappa = 1: tau = sqrt(sigma²/lambda); sigma² read live (step 6g); lambda the drift fact per throttle level (6h); start-up schedule; bias loop off; 1/tau travel at 128 ms; floor dt/0.45 | setpoint rate and torque, with the motor lead |
| SYND | SYN1 with a dynamic kappa in flight: a bank of five windows c x tau, c in {0.5, 0.7, 1, 1.4, 2}, selected by a 10 s EMA of the accelerometer innovation | as SYN1 |
| SYNB | SYN1 with c calibrated once against truth on the training seeds, then fixed (the gimbal-bench recipe) | as SYN1 |
| PIDFF, SYN1noFF | diagnostics: the PID with SYN1's feed-forward, and SYN1 without it | |

## Amendments (all written before any test-seed run)

| # | sha256 (start) | What was seen | Change |
|---|---|---|---|
| 0 | 39e240d0 | (pre-registration, predictions F1-F4, K1, K2, JPq) | |
| 1 | 4876f5d0 | Kalman R, Q_G and the law's c on the edges of their grids | grids widened (R x 2^-8..20, Q_G x 2^-10..8, c up to 8) |
| 2 | c89b386c | Kalman back on the far edges: it ignored the accelerometer. Cause: the gyro records had their mean removed, and a mean-free looped record never drifts. | gyro records keep their logged bias (-0.037 / -0.036 / -0.026 deg/s idle / low / half); results in `results_a2/` |
| 3 | 0475868d | all three arms on the top edge of the gain grid | gain grid widened |
| 4 | 465c138f | KF on the top rate Kp, PID on the top Kd | gain grid widened again; all three arms then interior |

After amendment 2 the Kalman's Q_G stayed on the bottom edge, about zero: the Kalman wants no angle random walk at all. Widening further cannot change anything, so this was declared and left as it was.

Tuned values:

- **Kalman:** R = var(tilt) x 8192, Q_G about 0. Estimation error on training 0.239°.
- **Law factor:** c_bench = 2.8 (training curve: 0.25 → 0.628°, 0.5 → 0.453, 0.7 → 0.365, 1 → 0.288, 1.4 → 0.247, 2 → 0.232, 2.8 → 0.231, 4 → 0.235, 5.6 → 0.239, 8 → 0.242).
- **Gains (Ka, rate Kp, Kd, Ki multiplier):** PID (16, 160/B, 1e-3, 2); KF (32, 320/B, 5e-4, 0); SYN1 (16, 160/B, 5e-4, 2).

## Results (test seeds 200-209; pooled RMS over the seeds, in degrees)

| Arm | Attitude error | Estimation error | Attitude low / half throttle | Motor jitter |
|---|---|---|---|---|
| PID | 2.511 | 2.353 | 1.121 / 3.369 | 0.120 |
| PIDFF | 2.330 | 2.353 | 0.700 / 3.220 | 0.119 |
| KF | 0.603 | 0.231 | 0.600 / 0.606 | 0.154 |
| SYN1noFF | 0.991 | 0.273 | 0.996 / 0.987 | 0.144 |
| **SYN1 (kappa = 1)** | **0.476** | 0.273 | 0.449 / 0.501 | 0.141 |
| SYND (dynamic kappa) | 0.481 | 0.282 | 0.441 / 0.518 | 0.141 |
| SYNB (c = 2.8) | 0.455 | 0.233 | 0.441 / 0.467 | 0.141 |
| KFFF (exploratory: KF + the same feed-forward) | 0.453 | 0.231 | 0.442 / 0.464 | 0.141 |

Wins in attitude: SYN1 beats PID in 10/10 runs and KF in 10/10. SYNB beats KF in 10/10. In estimation SYNB beats SYN1 in 7/10. Against KFFF (exploratory), SYN1 wins 3/10 and SYNB 4/10.

## Verdicts

| Prediction | Bar | Result | Held? |
|---|---|---|---|
| F1 the law's estimator against the flown fixed window | <= 0.5 x PID | 0.273 against 2.353 (-88 %) | held |
| F2 the law's estimator against the tuned Kalman | <= KF | 0.273 against 0.231 (+18 %) | **failed** |
| F3 SYN1 flies better than PID | <= 0.8 x, >= 8/10 | -81 %, 10/10 | held |
| F4 SYN1 flies at least as well as PID + tuned Kalman | <= KF, >= 7/10 | -21 %, 10/10 | held |
| K1 (assistant) the dynamic kappa in flight hears nothing | SYND >= 0.95 x SYN1 | 0.282 against 0.273 | held |
| K2 (assistant) pipeline bias, c_bench <= 0.7 | | c_bench = 2.8, the opposite direction | **failed** |
| JPq (author) kappa tends to 1, c_bench in [0.8, 1.25] | | c_bench = 2.8 | **failed** |

## Exploratory (after the verdicts, not judged)

- **The 21 % over the PID + Kalman comes from the feed-forward, not from the law.** The Kalman given the same feed-forward and motor lead (KFFF) flies at 0.453°, slightly better than SYN1 (0.476°). Removing the feed-forward from SYN1 doubles its error (0.991°). Adding it to the PID barely helps (2.33°), because the PID's problem is its estimator.
- **What the law does without truth:** a Kalman set from the same facts, with no tuning against truth (R = per-sample tilt variance, Q_G = drift fact), estimates at 0.583° (Q at idle) and 0.589° (Q at half throttle). On the same test seeds the law at kappa = 1 gets 0.272°, 53 % better. The truth-tuned Kalman needed a 551-point search against the true angle. A drone in flight does not have that truth; a gimbal bench with encoders does.
- **The law at c = 2.8 matches the tuned Kalman:** 0.233° against 0.231°. Once its factor is calibrated against truth, the law reaches the same optimum.
- **Why the dynamic kappa hears nothing (rough estimate):** going from c = 1 to c = 2.8 changes the prior variance by about 0.03 deg². The accelerometer innovation carries the per-sample tilt noise, about 31 deg² at low throttle and about 2500 at half. That makes the signal 1e-3 to 1e-5 of the noise. A 10 s EMA at 250 Hz has a relative error of at least about 2 %. The signal is 20 to 2000 times below what the selector can resolve.
- **Why c = 2.8 rather than 1 (hypothesis, untested):** the tuned Kalman wants Q_G about 0, so the true angle random walk is far smaller than the drift fact lambda says. lambda at half throttle comes from one Allan point at 1 s on an 8 s looped record (7 differences), and it sits on top of the bias wander (-0.019 then -0.032 deg/s in the record's two halves). An inflated lambda shortens the window by sqrt(inflation): c² = 7.8 would mean lambda is overestimated about 8 times. The assistant's K2 had predicted the bias in the other direction (from step 6h); it was wrong.

## Status in plain terms

- **Best syntonic flight controller built here:** SYN1, the law at kappa = 1 with everything the benches taught. In attitude it flies 81 % better than the PID with the flown fixed window and 21 % better than PID + truth-tuned Kalman, in 10/10 runs both times.
- **Fair reading:**
  - Against the PID, the gain is the law's: the estimator goes from 2.35° to 0.27°.
  - Against the Kalman, the 21 % is the feed-forward. On equal feed-forward it is a tie, and the Kalman is slightly ahead (KFFF 0.453°).
  - As a pure estimator, the truth-tuned Kalman is 18 % better than the law at kappa = 1.
- **The law's real advantage:** it needs no truth. Set from the same facts, a Kalman is twice as bad (0.58° against 0.27°). Calibrated once against truth (c = 2.8), the law equals the tuned Kalman.
- **kappa:** kappa = 1 is not the optimum on this bench (JPq failed). The likely cause is an inflated drift fact at half throttle, not the formula; this is untested.
- **Dynamic kappa in flight:** useless here, as predicted; the signal is far below the noise.

## Caveats

- One axis and a linear tilt model.
- A simulated plant.
- The only half-throttle record is 8 s, and it is looped.
- Four amendments, all before the test seeds.
- The Kalman's Q_G sits on the bottom edge, about zero.
- KFFF and the facts-set Kalman are exploratory, run after the verdicts.

## Next, in suggested order

1. A long half-throttle log (60 s or more) to measure lambda properly. Prediction to pre-register: c_bench moves toward 1.
2. Gimbal bench with encoders: calibrate c once in the real world (the SYNB recipe).
3. Any public claim against a Kalman must use KFFF (equal feed-forward) as the rival.

## Files and how to rerun

| File | Role | sha256 |
|---|---|---|
| `flight_bench.py` | pre-registration: bench, controllers, predictions, modes smoke / facts / calibrate / tune / run / summarize | 39e240d07833571ce0a978d29571fe9c635e335755580547161f9fb7f8c918e1 |
| `flight_bench_a1.py` | amendment 1, widened estimator grids (calibrate) | 4876f5d0b08c4d3a290eed6dfecfe2a0e949208c8d678531260df06e09dbdc10 |
| `flight_bench_a2.py` | amendment 2, the gyro keeps its logged bias; all later steps write to `results_a2/` | c89b386c02be84260bef6a4685985cb545dc9216dc288e383fc65c867cbdc8a7 |
| `flight_bench_a3.py` | amendment 3, widened gain grid (tune) | 0475868d52b91fce045e7fd5ec5ca1bf78d65e8f96ec1f08028bdde138d1b5c0 |
| `flight_bench_a4.py` | amendment 4, gain grid widened again (tune, run, summarize) | 465c138fe6c354916ba07d8d5b6ae071a8c0339ed6af22ec76cac6bcfa2c5624 |
| `explore_kfff.py`, `explore_kf_facts.py` | exploratory, after the verdicts | (paths made relative for the repository) |

`results/` holds amendment 1's calibration (superseded by amendment 2; the first calibration of `flight_bench.py` was overwritten by it, and its edge values are quoted in amendment 1's docstring). `results_a2/` holds everything judged: `facts.json`, `tuned_estimators.json`, `gains_*.json`, `test_runs.jsonl` (the ten test seeds) and `summary.json` (verdicts), plus the two exploratory outputs.

From the repository root, in order (each `tune` and the `run` take tens of minutes):

```
python bench/flight_2026-11-02/flight_bench_a2.py . calibrate
python bench/flight_2026-11-02/flight_bench_a4.py . tune PID     # then KF, SYN1
python bench/flight_2026-11-02/flight_bench_a4.py . run
python bench/flight_2026-11-02/flight_bench_a4.py . summarize
python bench/flight_2026-11-02/explore_kfff.py
python bench/flight_2026-11-02/explore_kf_facts.py
```

`run` appends to `results_a2/test_runs.jsonl`; move the recorded file aside first.
