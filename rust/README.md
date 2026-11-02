# syntonic-engine (Rust)

The syntonic engine as a Rust crate. It is `no_std`, allocates nothing, and takes its math from `libm`, so the same code runs on a PC and on a microcontroller. The arithmetic is the same as in the Python and C++ versions.

```rust
use syntonic_engine::{Config, Engine};

let cfg = Config { dt: 0.001, frame: 0.003, ..Config::default() };   // 1 kHz loop, 333 Hz servo frames
let mut eng: Engine<6> = Engine::new(&q0, &w0, cfg);                  // q0: start angles (rad), w0: no-load speeds (rad/s)

// on each tick, with two reads of every encoder:
let cmd = eng.step(&q_plan, &read_a, &read_b, None);                 // None: the engine counts servo frames itself
// send cmd to the servos on each frame; eng.tau holds each joint's window (s)
```

Two options, off by default since 0.2.9 (with both off the arithmetic is that of 0.2.8): `Config::lead`, a servo lead on the planned velocity (negative: servo_tau plus the loop delay, from the facts the engine already takes), used through `eng.step_lead(&q_plan, &qd_plan, &read_a, &read_b, None)` (`step` panics if the lead is on); and `Config::travel`, `Travel::E` to let the window travel in 1/tau (the e-geodesic) with `Config::travel_t` (default: the longest lag, 64 ms). `Config` gained three fields, so build it with `..Config::default()`.

`Engine<N, RING>` handles `N` joints. `RING` is the length of the residual history in loop ticks and must be larger than the longest lag in ticks: 65 is enough at 1 kHz, and the default of 129 also covers 2 kHz. For 12 joints the engine takes about 14 KB with the default `RING` (8 KB with `RING = 65`), so on a small board keep it in a `static` rather than on the stack.

## Checking it against Python

`tools/make_replay_log.py` runs the Python engine in the servo model of `examples/servo_bench.py` for 33 s and records every input and output. The replay example feeds the same inputs to the Rust engine and compares the commands:

```
python rust/tools/make_replay_log.py                 # from the repository root
cd rust && cargo run --release --example replay -- ../replay_log.bin ../replay_w0.bin raw
cd rust && cargo run --release --example replay -- ../replay_log_allan.bin ../replay_w0.bin allan
cd rust && cargo run --release --example replay -- ../replay_log_lead.bin ../replay_w0.bin raw-lead          # options on
cd rust && cargo run --release --example replay -- ../replay_log_lead_allan.bin ../replay_w0.bin allan-lead
```

On our machine the largest difference over the 33,000 ticks was 5.4e-14 rad with the raw identification and 8.9e-16 rad with the block-mean one (the default since 0.2.0), with identical final windows on all twelve joints in both cases; with the lead and the e travel on, 6.2e-14 and 1.3e-14 rad. `Config::identification` selects between `Identification::Allan` and `Identification::Raw`. `cargo test` runs the unit tests and the example in the documentation.

The crate builds as `no_std` on a PC. It has not yet been cross-compiled for or run on a Cortex-M board such as the Teensy 4.x (target `thumbv7em-none-eabihf`).

## The two laws of 0.3.0 (since 0.3.1)

`SyntonicLoopLaw` and `RiccatiAttitudeLaw` are in the crate with the facts and the arithmetic of the Python classes; see their module documentation and the Python docstrings.

```rust
use syntonic_engine::{AttitudeConfig, AttitudeLevel, LoopLawConfig, RiccatiAttitudeLaw, SyntonicLoopLaw};

// a servo joint's PI gains, recomputed every 0.5 s from the live statistics; dead = the servo's delay, measured once
let cfg = LoopLawConfig { dt: 0.001, frame: 0.003, servo_tau: 0.025, dead: 0.0025, ..LoopLawConfig::default() };
let mut law: Box<SyntonicLoopLaw<2>> = Box::new(SyntonicLoopLaw::new(&q0, &w0, cfg));   // about 114 KB: not on the stack
let cmd = law.step(&q_plan, &qd_plan, &read_a, &read_b, frame_now);                    // every tick

// a one-axis attitude estimator; per level of the known input that changes the noise (throttle), the facts measured once
let levels = [AttitudeLevel { s2: s2_idle, lam: lam_idle, t_rec: 0.0 },                // level 0: the ground level
              AttitudeLevel { s2: s2_low, lam: lam_low, t_rec: t_low },
              AttitudeLevel { s2: s2_half, lam: lam_half, t_rec: t_half }];
let mut att = RiccatiAttitudeLaw::new(&levels, q_b, 5.0, AttitudeConfig::default());
let theta = att.step(gyro, tilt, 1);   // every gyro sample; tilt: Option<f64>, None between accelerometer samples
// att.arm() at the end of the ground calibration; from then on pass the gyro with the calibration subtracted
```

The loop law's gain recomputation evaluates about 1.3 million complex responses. With `LoopLawConfig::spread = 0` (the default) it runs in one tick, as in Python: about 4 ms for two joints on our PC, and by LLVM's Cortex-M7 model about 0.2 to 0.3 s on a Teensy 4.1. `spread = m` divides it over the m ticks that follow; the gains are the same but take effect m ticks later, which raises the tip error on the realistic MAESTRO 1 bench by about 1.5 to 2 % at m = 250 and 4 to 5 % at m = 499 (`bench/port_2026-11-02/`).

`tools/make_law_replay_log.py` runs the two Python laws on deterministic inputs and records every tick; the replay example compares the Rust laws with them:

```
python rust/tools/make_law_replay_log.py              # from the repository root; add --flight-bench for the flight replay
cd rust && cargo run --release --example replay_laws -- ../replay_loop_bench.bin ../replay_loop_bench_spread.bin \
        ../replay_loop_two.bin ../replay_loop_one.bin ../replay_att_synth.bin
```

On our machine the loop law chose the same gains as Python at every tick of every replay (also with `spread = 400`), and its other outputs agreed to 2.4e-13 of their scale; the attitude law agreed exactly (difference 0.0), also on the realistic flight bench.
