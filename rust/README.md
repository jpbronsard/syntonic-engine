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

`Engine<N, RING>` handles `N` joints. `RING` is the length of the residual history in loop ticks and must be larger than the longest lag in ticks: 65 is enough at 1 kHz, and the default of 129 also covers 2 kHz. For 12 joints the engine takes about 14 KB with the default `RING` (8 KB with `RING = 65`), so on a small board keep it in a `static` rather than on the stack.

## Checking it against Python

`tools/make_replay_log.py` runs the Python engine in the servo model of `examples/servo_bench.py` for 33 s and records every input and output. The replay example feeds the same inputs to the Rust engine and compares the commands:

```
python rust/tools/make_replay_log.py                 # from the repository root
cd rust && cargo run --release --example replay -- ../replay_log.bin ../replay_w0.bin
```

On our machine the largest difference over the 33,000 ticks was 5.4e-14 rad, with identical final windows on all twelve joints. `cargo test` runs the unit tests and the example in the documentation.

The crate builds as `no_std` on a PC. It has not yet been cross-compiled for or run on a Cortex-M board such as the Teensy 4.x (target `thumbv7em-none-eabihf`).
