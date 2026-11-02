# MAESTRO 1 joint bench, 1-2 Nov 2026: the servo lead, the e-geodesic, the size of the effect, and kappa

Simulation only, on the model of `examples/servo_bench.py` and `bench/cadence_2026-10-31/cadence_bench.py` (baton arm RJ1-RJ6 of MAESTRO 1, an ideal first-order servo at 25 ms, AS5048 encoder noise, loads from the inverse dynamics of the planned motion). Engine 0.2.8 (`identification='raw'`, the authority limit). Written and run by the assistant; the author proposed the ideas and gave each go. Not yet reproduced by the author.

**Method.** Each pre-registered script states its predictions in its docstring. It was hashed and frozen (file made read-only) after a smoke test, before its gates and runs; the hashes are below. Every gate reproduced its reference with a difference of exactly zero. Exploratory scripts are marked as such; nothing they found is a result until a pre-registered script on fresh seeds confirms it.

**Run.** Each script takes the repository root as its first argument, for example `python bench/maestro_2026-11-01/eg732_hold.py . run A steady`. Results land next to the script, in the `results*/` folders kept here.

## Summary

| Test | Question | Held | Failed |
|---|---|---|---|
| eg732_lead | Servo lead and e-geodesic cadence at the servos' ceiling (667 Hz on the 2 kHz grid) | L1: lead -10 % at 667 Hz. A1: the fastest fixed frame is the most accurate. | L2: lead on the e-geodesic cadence -7.4 % (bar -8 %). E2, E2L: cadence efficiency 3.5 % and 1.1 % (bar 5 %). |
| eg732_hold | Lead with the loop delay, (servo_tau + D) x velocity | D3: -4 % at every frame | D1, D2, D4: the slow/fast frame gap does not move; the cadence falls to the slowest frame |
| window_vs_frozen | Adaptive window against the window frozen at its own mean | W1-W3: 15 % without lead, 48 % with it | none; but that rival is weak (exploratory below: 4.5 % against the best fixed window) |
| window_transport | Window travel in 1/tau (the e-geodesic) against travel in tau | P2a: -19 % at 1 s (20/20). P2b: -3 % at 0.1 s (18/20). | none |
| transport_confirm (fresh seeds, T = 64 ms declared) | Same, and against the published travel | C2: 1/tau beats tau, 20/20. C3. | C1: against the published travel -3.9 % pooled, but better in only 11/20 runs |
| confirm_dose (fresh seeds) | The law against the best fixed window under large and frequent noise changes | K1: -14.6 % at a x100 noise step (10/10). K2: -16.4 % with 1/tau travel (10/10). K4: -6.2 % in bursts (10/10). | K3: 1/tau travel against the published travel in bursts, -1.2 % (bar -5 %) |
| pid_lead_fair (fresh seeds) | The PID given the same lead and frames, re-tuned for them | none | R1-R4: the PID is 29 % better under steady noise, level under the noise step; the lead helps it more |
| kappa_dyn2 (fresh seeds) | Is kappa_H placed right? A dynamic kappa per joint and shared; the author's predictions | Q1: published bench, best static factor 1.0 (regret 0). Q2: best configuration, best factor 0.7 (regret 3.6 %). JP1: the shared dynamic factor ends at 0.99 and 0.96. JP2: one shared kappa does as well as six. | D1: the dynamic kappa costs about 1 %. D2: in the best configuration it does not find 0.7. |

**In plain terms.**

- **Best configuration on this bench.** Fixed frames at the servos' fastest rate, a 2 kHz loop, and the lead (servo_tau + D) x planned velocity: 0.587 mm pooled, against 0.834 mm for 333 Hz without the lead (-30 %). That gain is classical engineering (feed-forward, faster frames), not the law.
- **What the law itself adds.** About 4.5 % against the best fixed window in ordinary conditions, about 15 % when the noise changes a hundredfold, and about 6 % under frequent bursts, with nothing tuned. Same profile as against the PID: parity in steady conditions, better when conditions change.
- **kappa_H.** From Theorem U' it places the window exactly right on the published bench. With the lead, the optimum moves to about 0.7, at a 3.6 % cost.
- **Dual Geodesy.** First measurable difference from the geometry, confirmed twice: the window should travel in 1/tau rather than in tau. Against the published travel it is a robustness gain, not a uniform one.
- **Against a PID given the same lead and frames, re-tuned for them** (section 12): the PID is 29 % more accurate under steady noise and level under the noise step. Every prediction of that test failed. What remains for the law is that it needs no tuning.
- **0.2.9.** The lead and the 1/tau travel are now options of the engine, off by default (section 11).

## Files

Pre-registered, unchanged since their runs:

| Script | sha256 | Results |
|---|---|---|
| `eg732_lead.py` | 507791a9cdd09b3fea9a570a1fbf3927eecbdf10bae7ae905dd29b6bb3732d64 | `results/` |
| `eg732_hold.py` | 64134c25d23890690727d34a358872d21f0aa37f592211521070c3c5844f7007 | `results_hold/` |
| `window_vs_frozen.py` | 6c7fccfa7c4f4c4f74411a5db7789b8ad9320a5aa5deef08451c3c29c44ea29b | `results_window/` |
| `window_transport.py` | 101655fa7c2ac85f2c6d3902d361a448f92e2eca05e2dcd54162a125b01304d9 | `results_transport/` |
| `transport_confirm.py` | 88c7d505db7abb407b6be9e285aa926d4f9d98c1a8b8c611f81af7481f6445c7 | `results_confirm/` |
| `confirm_dose.py` | ac0df048daba0cce0c8d91ea4362b725bbabe2e3e71689ea71a2995b9dfcd5ac | `results_confirm_dose/` |
| `kappa_scale.py` | 235ebf89671482c678023c71bd24156e556f7cfd18fe2721c8f45d0e001fb2ee | frozen, never run (the author asked for a dynamic kappa next) |
| `kappa_dyn.py` | 9e90b61ce8fb3f247ffb460fbd7734fd22073f887bbff1d3baec881a6c56b065 | `results_kdyn/` (gates only; superseded before its runs, when the author's two statements came) |
| `kappa_dyn2.py` | 26e1a497cabcdfd14e9c86a1a569a5f891a25ea1629a2823f0234e12b61f435a | `results_kdyn2/` |
| `pid_lead_fair.py` | 2b0d71f27454b842a9110b2da7689bda6620ede9d8caa94351c36db12af19536 | `results_pid_lead/` |

Exploratory (not judged), with the hard-coded paths of the conversation workspace made relative for the repository:

- **Acceleration and rate grid:** `accel_tau.py`, `accel_tau_async.py`, `accel_tau_rate.py`, `corr.py`, `plot.py` (figure `maestro_accel_window.png`), `mc_cell.py` (`grid.jsonl`).
- **Servo lead:** `ff_test.py`, `pid_ff.py`.
- **Follow-ups to the pre-registered tests:**
  - `diag_seed2.py`;
  - `sweep_fixed.py` and `confirm_fixed.py`, which write into `results_window/`;
  - `explore_dose.py` and `explore_dose2.py` (`results_dose/`, `results_dose_sweep/`);
  - `explore_ema.py`, `explore_oracle.py`, `explore_clair.py` (`results_ema/`, `results_oracle/`, `results_clair/`).
- **The 0.2.9 gate:** `engine_options_gate.py` (`results_engine_options_gate.json`).
- **After pid_lead_fair:** `explore_pid_cross.py` (`results_pid_lead/explore_cross.jsonl`).

## 1. The window against the acceleration (exploratory)

Run to make a figure for the public. Published servo bench arithmetic, steady noise. The harness reproduces the published syntonic tip error exactly (A, seed 5, 0.684 mm). Choices made before looking: excerpt A, seed 5, shoulder RJ2.

The window against the joint's own acceleration (1 s average, matching the engine's averaging time), Spearman, after the 1 s warm-up:

| Joint | A: rho | A: time at floor | B: rho | B: time at floor |
|---|---|---|---|---|
| RJ1 | -0.76 | 46 % | -0.79 | 63 % |
| RJ2 | -0.58 | 4 % | -0.83 | 50 % |
| RJ3 | -0.77 | 20 % | -0.82 | 53 % |
| RJ4 | -0.34 | 93 % | -0.36 | 96 % |
| RJ5 | -0.69 | 9 % | -0.84 | 51 % |
| RJ6 | -0.75 | 48 % | -0.79 | 66 % |

The window narrows when the gesture accelerates, on every joint. The elbow RJ4 sits on the 2.22 ms floor almost all the time.

**Faster servo frames** (seed 5; only the command frame and the loop period change). The floor is max(loop delay, dt/0.45). At a 1 kHz loop it stays at 2.22 ms (the TH41-A term) whatever the servo frame, so the loop must go to 2 kHz for the floor to drop, to 1.11 ms. The 732 Hz runs send frames asynchronously (the `frame=` flag of `SyntonicEngine.step`).

| Case | Floor | RJ4 at floor (A / B) | Tip adaptive / frozen, A (mm) | Tip adaptive / frozen, B (mm) |
|---|---|---|---|---|
| 333 Hz servo, 1 kHz loop (published) | 2.22 ms | 93 / 96 % | 0.684 / 0.757 | 0.852 / 0.917 |
| 732 Hz servo, 1 kHz loop | 2.22 ms | 93 / 96 % | 0.613 / 0.707 | 0.759 / 0.854 |
| 732 Hz servo, 2 kHz loop | 1.11 ms | 44 / 67 % | 0.501 / 0.570 | 0.630 / 0.699 |
| 1 kHz servo, 2 kHz loop | 1.11 ms | 44 / 67 % | 0.477 / 0.555 | 0.592 / 0.676 |

"Frozen" is the same engine with its window held per joint at the mean of its own adaptive run. That is a weak rival (section 4). The PID was not rerun at these rates.

## 2. The servo lead (exploratory)

The command sent at each servo frame is plan + servo_tau x planned velocity + correction. servo_tau (25 ms) is the authority fact the engine already declares. The engine's internal servo and plan models see the same lead, so its residual stays exogenous. One seed per case. Both scripts reproduce the published figures without the lead.

| Case | PID | PID + lead | Engine | Engine + lead |
|---|---|---|---|---|
| A steady | 0.676 | 0.660 | 0.684 | 0.599 |
| B steady | 0.704 | 0.669 | 0.852 | 0.737 |
| A noise step | 1.299 | 1.287 | 0.881 | 0.841 |
| B noise step | 1.314 | 1.294 | 1.097 | 0.839 |

- **For the engine:** the lead lowers its tip error by 12 to 13 % under steady noise and by 5 to 24 % under the noise step.
- **For the PID:** the same lead helps it by 1 to 5 % only, but the PID was not re-tuned with the lead.
- **Model match:** the simulated servo is exactly first order at the 25 ms the lead assumes, so the lead is as good as it can be. On hardware the data-sheet time constant will be off.

## 3. The e-geodesic cadence at the servos' ceiling, with the lead (`eg732_lead.py`, pre-registered)

**Setting.**

- 2 kHz loop.
- The servos accept at most 732 Hz (the author). On the 0.5 ms grid the fastest frame under that is 1.5 ms (667 Hz).
- Arms: fixed frames F3, F2 and F15, and the e-geodesic cadence EG (31 Oct mechanism), each without and with the lead (suffix L).
- Excerpts A and B, steady noise and the x10 noise step, seeds 1 to 5: 20 paired runs.

**Gates.**

- G1: `examples/servo_bench.py` reproduced exactly.
- G2: `ff_test.py` reproduced.

| Line | Prediction | Result | Verdict |
|---|---|---|---|
| L1 | F15L <= 0.92 x F15, below in >= 18 of 20 | 0.897, 18 of 20 | held |
| L2 | EGL <= 0.92 x EG, below in >= 18 of 20 | 0.926, 15 of 20 | failed (narrowly) |
| E2 | gain(EG) >= 5 % against the fixed frontier | 3.5 % | failed |
| E2L | gain(EGL) >= 5 % against the lead frontier | 1.1 % | failed |
| A1 | F15L <= EGL (form, not accuracy) | 0.613 against 0.687, 20 of 20 | held |

Pooled tip error (mm), with the mean frame rate in brackets:

| Without lead | | With lead | |
|---|---|---|---|
| F3 (333 Hz) | 0.834 | F3L (333 Hz) | 0.695 |
| F2 (500 Hz) | 0.735 | F2L (500 Hz) | 0.632 |
| F15 (667 Hz) | 0.684 | F15L (667 Hz) | 0.613 |
| EG (444 Hz) | 0.742 | EGL (334 Hz) | 0.687 |

**What it says.** With the lead on, the law's window never falls below any floor, so the e-geodesic cadence picks the slowest frame almost all the time. Yet the fastest fixed frame with the lead is 11 % more accurate (20 of 20). Once the lead is on, the window no longer measures what fast frames buy.

**Exploratory** (`diag_seed2.py`): the run where the lead lost (A steady, seed 2) shows a transient at 8.6 s on RJ2, RJ3 and RJ5. It is not the correction limit. Guess, unverified: the load-derated speed limit.

## 4. The lead with the loop delay (`eg732_hold.py`, pre-registered)

**Reading tested.** What fast frames still buy is the hold of the frame. The lead becomes (servo_tau + D) x planned velocity, with D = loop_delay(dt, frame) = frame/2 + dt/2 + 0.1 ms, the engine's declared delay. No new constant.

**Gates.** G1 (`servo_bench.py`) and G3 (the four lead arms of `eg732_lead.py`, run by run): reproduced exactly.

| Line | Prediction | Result | Verdict |
|---|---|---|---|
| D1 | F3D <= 1.05 x F15D (the frame gap closes) | 1.137 (was 1.133) | failed |
| D2 | F3D <= 0.92 x F3L, below in >= 18 of 20 | 0.961, 18 of 20 | failed |
| D3 | F15D <= F15L | 0.587 against 0.613, 18 of 20 | held |
| D4 | EGD <= 1.03 x F15D at <= 500 Hz | 1.106 at 337 Hz | failed |

Pooled tip error (mm): F3D 0.668, F2D 0.610, F15D 0.587, EGD 0.650.

The reading failed. Leading the hold improves every frame by about 4 %, but the slow/fast gap does not move. The law's cadence still falls to the slowest frame. With any lead, on this bench, the e-geodesic cadence has no use; its 31 Oct timing result stands as it was. Best configuration: F15D, 0.587 mm, 30 % below F3 without the lead.

## 5. The adaptive window against frozen and fixed windows (`window_vs_frozen.py`, pre-registered, and exploratory follow-ups)

**Setting.** F15 (no lead) and F15D (lead), adaptive. F15z and F15Dz are the same engines with the window frozen per joint at the mean of their own adaptive run. 20 paired runs.

**Gates.** G4 (F15D of `eg732_hold.py`) and G5 (F15 of `eg732_lead.py`): reproduced exactly.

| Line | Prediction | Result | Verdict |
|---|---|---|---|
| W1 | gap(F15) >= 5 %, below in >= 18 of 20 | 15.0 %, 20 of 20 | held |
| W2 | gap(F15D) >= 5 %, below in >= 18 of 20 | 47.9 %, 20 of 20 | held |
| W3 | gap(F15D) >= 0.5 x gap(F15) | 47.9 % against 15.0 % | held |

Pooled tip error (mm): F15 0.684, F15z 0.805, F15D 0.587, F15Dz 1.129. Mean windows with the lead (ms), RJ1-RJ6: 10.8, 67.7, 32.4, 7.6, 15.4, 3.8.

**Exploratory: is the frozen rival fair?** With the lead the windows spread widely, so a window frozen at its mean is far too slow in the fast passages.

- **Sweep** (`sweep_fixed.py`, seed 1): one fixed window for all joints, swept over 1.11 to 13 ms. The best was 5 ms: 0.624 mm pooled against 0.596 for the adaptive window.
- **Held out** (`confirm_fixed.py`): 5 ms on seeds 2 to 5, all four cells, 16 runs. Adaptive 0.585 mm, fixed 0.613 mm: the adaptive window is 4.5 % better, below in 12 of 16.
  - A steady: -3.7 % (the fixed window is better).
  - A step: +2.1 %.
  - B steady: +8.8 %.
  - B step: +7.2 %.

**What it says.** Against the best single fixed window, the law's window is about 4.5 % better. The fixed window needed a search over seven values; the law needed none. The 7 to 14 % gaps of section 1 were against the weak rival.

## 6. Where the e-geodesic makes a difference: the window's travel (`window_transport.py`, pre-registered)

The request: "trouve-moi la place de l'e-géodésique ici et qui fait une différence". Dual Geodesy, Theorem 1 (the Register's repaired P2): which coordinate relaxes exponentially, theta = 1/tau (e) or tau?

**Setting.** F15D. The law's target is unchanged in every arm; only the window's travel toward it changes.

- BASE: the published travel (1/tau stepped, at most 0.4 s/s; rarely binds).
- E1 and E01: theta = 1/tau relaxes exponentially, time constant 1 s (the engine's ema) or 0.1 s.
- N1 and N01: tau itself relaxes, same time constants.

**Gate.** G6 (F15D of `eg732_hold.py`): reproduced exactly.

| Line | Prediction | Result | Verdict |
|---|---|---|---|
| P2a | E1 <= N1 pooled, below in >= 15 of 20 | 0.5915 against 0.7262 (-18.6 %), 20 of 20 | held |
| P2b | E01 <= N01 pooled, below in >= 15 of 20 | 0.5704 against 0.5900 (-3.3 %), 18 of 20 | held |
| Null clause | both ratios within 1 % | 0.814 and 0.967 | not met: the coordinate matters |

Pooled tip error (mm): BASE 0.5874, E1 0.5915, E01 0.5704. Travelling in 1/tau shortens the window fast and lengthens it slowly, and on this bench lag costs more than noise. This is prescriptive (which coordinate the window should relax in). The repaired P2 asks which coordinate a physical stack does relax in (from telemetered tau(t)); not answered here.

## 7. The same on fresh seeds, with a declared time constant (`transport_confirm.py`, pre-registered)

**Declared rule.** T = the longest lag of the identification, 64 ms. Stated openly in the script: the rule was chosen after seeing that 0.1 s helped and 1 s did not, so only the fresh-seed result counts.

**Setting.** Seeds 21 to 25, 20 paired runs. Arms: BASE, E64, N64, and E01 (replication).

**Gate.** G7 (`window_transport.py`'s BASE): reproduced exactly.

| Line | Prediction | Result | Verdict |
|---|---|---|---|
| C1 | E64 <= 0.98 x BASE pooled, below in >= 14 of 20 | 0.961 pooled, but 11 of 20 | failed (consistency) |
| C2 | E64 <= N64 pooled, below in >= 15 of 20 | 0.971, 20 of 20 | held |
| C3 | E01 <= 0.98 x BASE pooled | 0.964 (12 of 20) | held |

Pooled tip error (mm): BASE 0.580, E64 0.557, N64 0.574, E01 0.559. Worst run: BASE 0.709, E64 0.668.

**What it says.** The coordinate result is replicated: travelling in 1/tau beats travelling in tau, 20 of 20. Against the published travel, the e travel removes the bad episodes the published travel sometimes has on excerpt A, and costs about 1 % in ordinary runs. It is a robustness gain, not a uniform one.

## 8. A more marked effect: the size of the change in conditions (`explore_dose*.py`, then `confirm_dose.py`, pre-registered)

The goal set by the author: "comment obtenir un effet plus marqué ?".

**Noise profiles** on the encoder read variance:

- steady;
- step10: x10 from 10 to 20 s;
- step100: x100 from 10 to 20 s;
- bursts: x30 for 0.5 s every 2.5 s.

**Arms.**

- BASE: the published travel.
- E64 and N64: 1/tau and tau relaxing at 64 ms.
- FIXv: the window fixed at v ms on every joint.

**Exploratory** (seeds 21 to 23):

- Against the 5 ms window, the law's lead grows with the size of the change: +4.4 % steady, +8.5 % at x10, +32 % at x100.
- Against the best fixed window for each profile (13 ms at x100, 5 ms in bursts, chosen on seed 21): +13.5 % at x100, and -5.2 % in bursts, where E64 turns it into +4.2 %.

**Confirmation** (seeds 26 to 30, 10 paired runs per profile). Gate G8: the exploratory values reproduced exactly.

| Line | Prediction | Result | Verdict |
|---|---|---|---|
| K1 | step100: BASE <= 0.92 x FIX13, below in >= 9 of 10 | 0.854 (-14.6 %), 10 of 10 | held |
| K2 | step100: E64 <= 0.90 x FIX13, below in >= 9 of 10 | 0.836 (-16.4 %), 10 of 10 | held |
| K3 | bursts: E64 <= 0.95 x BASE, below in >= 9 of 10 | 0.988 (-1.2 %), 6 of 10 | failed |
| K4 | bursts: E64 <= FIX5, below in >= 7 of 10 | 0.938 (-6.2 %), 10 of 10 | held |

**What it says.** The lever is the size of the change in conditions. Against the best fixed window chosen for the same conditions, the law's lead goes from about 4.5 % to about 15 % when the noise variance changes a hundredfold, and is about 5 to 6 % under frequent bursts. The exploratory 9 % of the e travel in bursts did not replicate. The X500's logs show a tilt-noise variance ratio of about 77 between idle and half throttle: the quadrotor's attitude estimator sits in the regime where the law's lead is largest.

## 9. How much better can it get? (exploratory, `explore_ema.py`, `explore_oracle.py`, `explore_clair.py`)

**Averaging time of the law's statistics** (ema, published 1 s), against the fair fixed window per profile, 6 runs each. Faster statistics matter under bursts, where a 1 s average cannot see a 0.5 s burst; not for a single large step.

| Profile | ema 1 s | ema 0.25 s | ema 64 ms | 1/tau travel + ema 0.25 s |
|---|---|---|---|---|
| steady | +4.4 % | +3.2 % | +4.6 % | +7.4 % |
| step100 | +11.2 % | +10.4 % | +11.7 % | +11.3 % |
| bursts | -3.9 % | -0.4 % | +5.8 % | +9.4 % |

**An oracle noise estimate failed.** Giving the engine the true noise variance every tick, while the structure function stays averaged over 1 s, made it much worse. The noise and drift estimates must share one time scale: a design constraint, not a ceiling.

**Clairvoyant two-window schedule** (one window in low noise, one in high noise, switched exactly at the regime changes; seeds 22-23):

| Profile | Best fixed | Published law | 1/tau travel + ema 0.25 s | Clairvoyant | Share of the possible gain captured |
|---|---|---|---|---|---|
| step100 | 1.177 | 1.018 | 1.023 | 0.958 | published 72 %, travel + ema 70 % |
| bursts | 0.740 | 0.779 | 0.649 | 0.643 | published -40 %, travel + ema 94 % |

The 0.25 s averaging time was chosen by hand. It needs a rule from declared facts and a pre-registered confirmation before it counts. A markedly larger effect has to come from the conditions the machine meets, not from more tuning on this bench.

## 10. kappa: static map, a dynamic kappa per joint and shared (`kappa_dyn2.py`, pre-registered)

The author's statements: "kappa = 1 is canonical only"; "k tend vers 1 et k est identique sur tous segments du même ensemble (4D)".

**Setting.** The engine takes kappa_H from Theorem U'. On the baton arm H is about 0.97, so kappa_H is about 0.64. The law's target is multiplied by a factor c (kappa -> c x kappa_H).

- **Static:** c = 0.5, 0.7, 1, 1.4, 2.
- **DYNJ:** per joint, extremum seeking on that joint's encoder error. Constants declared, not fitted.
- **DYNS:** one c for the six joints; cost = the measured tip error.
- **Learning signal:** neither learns from kappa_obs, which the Register V1.3 shows is circular.
- **Settings:** S1, the published bench (1 kHz, 3 ms frames, no lead); S2, the best configuration (2 kHz, 1.5 ms frames, lead (servo_tau + D)).
- **Runs:** seeds 31 to 35, 20 paired runs per setting.

**Gates.** G9 (`servo_bench.py`) and G10 (F15D of `eg732_hold.py`): reproduced exactly.

| Line | Prediction | S1 | S2 | Verdict |
|---|---|---|---|---|
| Q1 / Q2 | best static c in {0.7, 1, 1.4}, regret of c = 1 under 10 % | c* = 1.0, regret 0 % | c* = 0.7, regret 3.6 % | held, both |
| D1 | shared dynamic kappa no worse than c = 1 | +1.1 % | +1.3 % | failed, both |
| D2 | shared dynamic kappa within 2 % of the best static | 0.894 against 0.902 bar | 0.603 against 0.585 bar | held S1, failed S2 |
| JP1 (author) | shared final c in [0.8, 1.25] | 0.988 | 0.956 | held, both |
| JP2 (author) | per-joint final c within x1.25, shared within 2 % of per-joint | spread 1.14, 0.894 vs 0.893 | spread 1.11, 0.603 vs 0.603 | held, both |

Pooled tip error (mm), static c = 0.5 / 0.7 / 1 / 1.4 / 2:

- S1: 0.953 / 0.911 / 0.884 / 0.893 / 0.977.
- S2: 0.598 / 0.574 / 0.595 / 0.642 / 0.734.

**What it says.**

- **Published bench:** kappa_H places the window exactly right; no calibration is needed.
- **With the lead:** the optimum moves to about 0.7. The lead removes the acceleration part of the residual, so the class picture behind kappa_H no longer fits exactly.
- **The dynamic kappa** settles near 1, and one shared factor does as well as six separate ones. This supports treating the arm as one class block. It is one simulated arm, not a change to the Canonical Edition's non-claim ("no common kappa imposed").
- **Read with care:** the learner starts at 1, its gradient is weak, and in S2 it did not reach 0.7.

## 11. The lead and the e travel as engine options (0.2.9)

The two pieces that held up are now options of `SyntonicEngine` (Python, Rust and C++), off by default. With both off the engine is that of 0.2.8 bit for bit: the replay logs of `rust/tools/make_replay_log.py` are byte-identical, and so is the output of `examples/servo_bench.py`.

- `lead='servo+delay'` (a negative number in Rust and C++): the command adds (servo_tau + loop_delay(dt, frame, sensor_delay)) x the planned velocity. The engine's servo and plan models see the same lead. Pass `qd_plan` to `step` (Rust and C++: `step_lead`). This is the F15D lead of section 4.
- `travel='e'`: theta = 1/tau relaxes exponentially toward 1/tau_target with `travel_T`, by default the longest lag (64 ms). The window takes the target at the first fit. This is the E64 arm of section 7.

`engine_options_gate.py` runs the published engine with these options in the harness of `window_transport.py`. It reproduces BASE (lead only) and E64 (lead + e travel) of `transport_confirm.py` with a difference of exactly 0.0 mm (A steady seed 21, B step seed 25).

## 12. The PID given the same lead and frames, re-tuned (`pid_lead_fair.py`, pre-registered)

Asked by the author ("FAIT LE !") after the assistant pointed out that the 30 % of sections 3-4 was the engine against itself. Script sha256 2b0d71f27454b842a9110b2da7689bda6620ede9d8caa94351c36db12af19536, frozen after a smoke test, before its gates, tuning and runs.

**Setting.**

- Baton arm, 2 kHz loop, fixed 1.5 ms frames.
- The same lead for both, (servo_tau + D) x planned velocity: 26.1 ms.
- The PID re-tuned the published way: joint by joint, on the published grid of 280 (Kp, Ki, Kd), lowest joint error, steady noise, per excerpt, tuning seeds 61-63.
- Declared beforehand: the PID was tuned on the same two excerpts it is tested on (other noise draws), which favours it.
- Edge rule declared beforehand; no joint landed on an edge.
- Test seeds 41-45: 20 paired runs.

**Gates.** G1 (`examples/servo_bench.py`, PID and engine, A steady seed 5 and B step seed 7) and G2 (F15D and F15 at the fast setting): reproduced exactly.

**Tuned gains.** With the lead the PID wants lower gains than the published ones. Kp 0.5 to 4, Ki 160 to 640 /s, Kd 0, against Kp 3 to 4 and Ki 640 to 1280 published. Without the lead, at the same frames: Kp 2 to 4, Ki 640 to 1280.

| Line | Prediction | Result | Verdict |
|---|---|---|---|
| R1 | noise step: SYND <= 0.85 x PIDD, below in >= 9 of 10 | 0.630 against 0.615 (x1.02), 4 of 10 | failed |
| R2 | steady noise: SYND <= 1.10 x PIDD | 0.491 against 0.350 (x1.40), 0 of 10 | failed |
| R3 | overall: SYND <= PIDD, below in >= 12 of 20 | 0.565 against 0.500 (x1.13), 4 of 20 | failed |
| R4 | the lead helps the engine more than the PID | engine -16.7 %, PID -29.3 % | failed |

Pooled tip error (mm):

| Arm | Steady | Noise step | All |
|---|---|---|---|
| PIDD (PID + lead, re-tuned) | 0.350 | 0.615 | 0.500 |
| SYND (engine + lead) | 0.491 | 0.630 | 0.565 |
| SYNDE (engine + lead + e travel) | 0.492 | 0.615 | 0.557 |
| PID15 (PID, no lead, re-tuned) | 0.472 | 0.882 | 0.707 |
| SYN15 (engine, no lead) | 0.580 | 0.764 | 0.678 |
| PIDDp (PID + lead, published gains) | 0.430 | 0.835 | 0.664 |

By cell (PIDD / SYND / SYNDE): A steady 0.327 / 0.441 / 0.440; A step 0.580 / 0.611 / 0.580; B steady 0.371 / 0.537 / 0.539; B step 0.648 / 0.648 / 0.648. Joint error (deg): PIDD 0.0205, SYND 0.0238. Tip jitter (mm): PIDD 0.463, SYND 0.456.

**Exploratory, after the verdicts** (`explore_pid_cross.py`): the PID with the gains tuned on the other excerpt (A's on B, B's on A) loses almost nothing: 0.364 steady, 0.623 step, 0.510 overall. The engine is below it in 0 of 10 steady runs and 4 of 10 step runs. The PID's lead does not come from being tuned on the test excerpts.

**What it says.**

- **Same lead, same frames, PID re-tuned:** the PID wins. Under steady noise its tip error is 29 % lower than the engine's (0 of 10 runs for the engine). Under the x10 noise step the two are level: the engine is 2 % behind, and the engine with the e travel exactly level.
- **The engine's advantage under the noise step at the published setting was largely the PID's high gains.** Those gains pass the extra noise through. With the lead, the PID's best gains drop, and the noise step costs it much less (0.882 to 0.615).
- **At the fast frames alone, without the lead,** the re-tuned PID also beats the engine under steady noise (0.472 against 0.580), and the engine still wins under the noise step (0.764 against 0.882).
- **What survives for the law:** it needed no tuning. The PID needed a 280-point search per joint and per configuration. Its published gains, carried to the new configuration, are 33 % worse than its re-tuned ones (0.664 against 0.500). The law's case on this bench is portability (no re-tuning when the machine changes), not accuracy against a re-tuned PID.
- **For any public claim:** at the fast setting with the lead, a re-tuned PID is more accurate than the syntonic engine on this simulated arm.

**Limits of everything above.**

- Simulation, with the servo exactly first order at the 25 ms the lead assumes. On hardware the lead will cancel less.
- 667 Hz stands in for the servos' 732 Hz.
- Two excerpts of one piece.
- The PID, given the same lead and frames and re-tuned, is more accurate than the engine (section 12).
