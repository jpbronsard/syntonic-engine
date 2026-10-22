# Syntonic Engine

**A joint correction loop for robots whose time window comes from the Syntony Principle, with nothing tuned.**

Jean-Pierre Bronsard · SyntonicAI Recherche · Montréal  
ORCID: [0009-0008-6639-7553](https://orcid.org/0009-0008-6639-7553)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23549690.svg)](https://doi.org/10.5281/zenodo.23549690)

---

I developed this engine for MAESTRO 1, a two-arm robot conductor I designed. It corrects each joint toward its planned motion, and the time window it uses for that correction comes from the Syntony Principle (Theorem U'). Nothing in it is tuned.

Every joint controller has to decide how long to average its sensor before it corrects. Average too briefly and the joint chases sensor noise. Average too long and it falls behind the motion. A PID makes that decision through gains found by trial, one joint at a time. This engine measures what the decision depends on while the robot moves (the sensor noise, and how fast and how smoothly the motion changes) and computes the window directly:

```
tau* = kappa_H * (sigma^2 / lambda_H) ^ (1 / (2H + 1)),     kappa_H = (2H * Gamma(2H + 1)) ^ (-1 / (2H + 1))
```

Here sigma^2 is the encoder noise intensity (the variance of one read times the loop period), and lambda_H and H describe how strongly and how smoothly the motion changes. With H = 1/2 the formula reduces to the classical square-root law, tau* = sqrt(sigma^2 / lambda). The window never drops below the loop delay, because a correction cannot usefully act faster than the loop can respond.

## How it works

On every loop tick, for every joint:

1. Read the encoder twice. The difference between the two reads gives the noise, and their mean gives the measurement.
2. Compare the measurement with a nominal model of the servo, driven by the plan and by the commands actually sent. What is left over is what the correction has to cancel. Measuring it this way keeps the engine from reading its own corrections back as motion.
3. Track how that residual changes over lags of 2 to 64 ms, keeping only the lags that stand clear of the noise, and fit its smoothness H and intensity lambda_H.
4. Compute the window from the formula, keep it between the floor described below and 0.5 s, and let it move smoothly.
5. Integrate the error over that window and add the result to the planned position at each servo frame.

## What you need to provide

Three facts about your hardware. There is nothing to tune.

| | MAESTRO 1 |
|---|---|
| loop period `dt` | 1 ms (1 kHz) |
| servo frame period `frame` | 3 ms (333 Hz PWM) |
| servo no-load speed `w0` and nominal reaction time `servo_tau` | from the data sheet; 25 ms if it gives none |

The window floor follows from these. It is the loop delay (half a frame, plus half a tick, plus 0.1 ms for the sensor) or the discretisation limit dt / 0.45, whichever is larger. For a 1 kHz loop driving 333 Hz servos that comes to 2.22 ms. The formula leaves out the servo's own internal latency, so on real hardware measure the full loop delay once and, if it is larger, pass it as `floor`.

Each encoder has to be read twice per tick, at least one internal sample apart so that the two reads carry independent noise. An AS5048 refreshes at about 11 kHz, so the reads should be about 0.1 ms apart, which fits easily in a 1 ms tick.

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

The bench simulates the twelve arm joints of MAESTRO 1 conducting the first movement of Beethoven's Fifth. Each joint has a servo model (333 Hz command frame, 25 ms internal loop, speed derating, compliance, gear play), AS5048 encoder noise, and loads from the inverse dynamics of the planned motion. The opponent is a PID with interpolation. Its gains were picked joint by joint from a grid of 280 (P, I, D) combinations, each tried on all twelve joints over the rest of the same movement with normal noise. Both run on two held-out 30 s excerpts, once with steady noise and once with the encoder noise variance raised tenfold for 10 s. The figures below are for the baton arm and come from `examples/servo_bench.py`. A ratio below 1 means the syntonic engine did better.

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

## Checks

- The Python engine in `python/syntonic_engine` reproduces the bench engine exactly: zero difference in joint angles and windows over 33 s on 12 joints.
- On the same 33 s, 12-joint replay, the C++ header agrees with the Python engine to 1.5e-14 rad and the Rust crate to 5.4e-14 rad.
- `examples/arm_fk.py`, the baton-tip kinematics in plain numpy, agrees with MuJoCo to 7e-16 m over 200 random poses.
- Neither the C++ header nor the Rust crate has run on a microcontroller yet. The Rust crate builds as `no_std` but has not been cross-compiled for a Cortex-M board, and the Teensy sketch has only been compiled on a PC against stand-ins for the Arduino functions.

To run the checks yourself, from the repository root (Python 3.9 or later, a C++11 compiler, and a recent stable Rust):

```
pip install -e python scipy
python -m unittest discover -s python/tests          # unit tests of the Python engine
python examples/servo_bench.py                       # the results table above, about a minute
python rust/tools/make_replay_log.py                 # writes replay_log.bin and replay_w0.bin
g++ -O2 -std=c++11 -Icpp cpp/tests/replay.cpp -o replay_cpp && ./replay_cpp replay_log.bin replay_w0.bin
cd rust && cargo test && cargo run --release --example replay -- ../replay_log.bin ../replay_w0.bin
```

## What is in the repository

```
python/syntonic_engine/   the engine (numpy), with unit tests in python/tests/
cpp/SyntonicEngine.h      the engine as a single C++ header, with a replay test in cpp/tests/
rust/                     the engine as a no_std Rust crate, with a replay check against Python
examples/                 servo_bench.py (the table above, numpy + scipy), arm_fk.py, data/, teensy/
bench/                    the full bench scripts and every result file; the scripts need the MAESTRO 1
                          simulation model, which is not included, but all results are
docs/bench_note.md        the bench note: setup, the predictions written before each run, results and failures
```

## Citing

If you use this work, please credit me (see NOTICE and CITATION.cff):

- this software: J.-P. Bronsard, *Syntonic Engine* v0.1.0, Zenodo, https://doi.org/10.5281/zenodo.23549690
- the theory: J.-P. Bronsard, *The Syntony Principle*, Zenodo, https://doi.org/10.5281/zenodo.17254395

## Licence

Apache License 2.0 (see LICENSE). Any redistribution or derivative work must keep the NOTICE file.

## Authorship

Code and documentation were drafted with language models under the author's direction, as solution architecture and integration of independent components. The author alone approves every statement and is accountable for it; no language model is an author.

Jean-Pierre Bronsard, SyntonicAI Recherche, 2026.
