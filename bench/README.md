# MAESTRO 1 joint bench: scripts and results

These are the scripts behind every number in the main README and in `docs/bench_note.md`, together with their result files.

The scripts import the MAESTRO 1 simulation model (presto.py, sim.py, maestro_model.py and the CAD data), which is not part of this package. They are here so that every step can be read and checked. If you want to run something yourself, `examples/servo_bench.py` reproduces the bench's main table without the model.

| Script | What it tests | Results |
|---|---|---|
| flow_vs_pid.py, flow_v3.py | the plant, the plan, the PID with interpolation, the best fixed window, and the tuning on the rest of the movement | tuning_v5_pair_*.json |
| flow_v5_uprime.py | the syntonic engine (Theorem U') against the tuned baselines | summary_v5c.json |
| flow_v6_rate.py | the loop at 2 kHz, with world and controller every 0.5 ms | summary_v6_*.json, tuning_v6_*.json |
| flow_v7_ctrl.py | the controller at 1 kHz, 500 Hz and 333 Hz, with the world at 1 kHz | summary_v7_*.json, summary_v7w_*.json |
| flow_v8_floor.py | the window floor set by the loop delay, across loop rates | summary_v8_*.json |
| flow_v9_frame.py | the window floor across servo frames of 2, 3, 5 and 8 ms | summary_v9_frame.json |
| flow_v10_seeds.py | the syntonic engine against the PID over ten paired noise draws | summary_v10_seeds.json |

Each script says in its header what was fixed before it ran and which predictions it tested. `docs/bench_note.md` records which of those predictions held and which failed.

`cadence_2026-10-31/` is self-contained (it runs on the model of `examples/servo_bench.py`, no MAESTRO 1 simulation needed): the servo frame chosen live by the law, each joint taking the slowest frame whose floor stays under its window, against fixed frames of 1, 2 and 3 ms; then, to 1 Nov 2026, Theorem Omega on the joints' phase rates, why the joint errors cancel at the tip (a velocity lag, then the acceleration), the causal test by re-timing the gesture, the correction limit derived from the authority fact, and the published bench over ten paired draws with that limit. Its `NOTE.md` records the predictions and which held. Its scripts up to `authority_limit.py` ran with the 5 degree limit of 0.2.7; with the 0.2.8 default they no longer reproduce their recorded results, so check out commit 23ea17f to rerun them.

`maestro_2026-11-01/` is self-contained in the same way (it imports `cadence_2026-10-31/cadence_bench.py`): the servo lead and the e-geodesic cadence at the servos' ceiling, the lead with the loop delay, the adaptive window against frozen and fixed windows, the window's travel in 1/tau against tau (twice, the second on fresh seeds), the law's effect against the size of the change in conditions, its headroom, and kappa (static map, a dynamic kappa per joint and shared, the author's predictions). Its `NOTE.md` records the predictions, which held and which failed, the exploratory passes, and `engine_options_gate.py`, which checks that the 0.2.9 engine options reproduce the bench arms they come from.

`flight_2026-11-02/` (it imports `software_2026-10-26/step6g_real_noise.py` and its logs): a one-axis flight controller with the law's attitude estimator at kappa = 1, against a PID with the flown fixed window and a PID with a Kalman filter tuned against truth, on the X500's logged vibration noise; with its four declared amendments and two exploratory passes. Its `NOTE.md` has the verdicts.

`goal_2026-11-02/` holds the pre-registered confirmations of 2 Nov 2026, every frozen file byte for byte. First, the syntonic loop law against PIDs tuned at commissioning, on a realistic MAESTRO 1 servo bench with 333 Hz frames (held). Then the attitude laws against Kalman filters tuned against truth: on the looped flight bench (failed, kept on record), then on a realistic one, with block-bootstrap noise and an unknown initial attitude, declared after that failure (held). The exploration harness they came from is there too. Its `NOTE.md` gives the sequence, the verdicts, the limits and how to replay. The frozen flight scripts name absolute paths of the working session, so `replay_setup.sh` creates them as links. `package_gate.py` checks that the two laws added to the Python package in 0.3.0 reproduce the confirmations exactly.

`port_2026-11-02/` goes with the C++ and Rust ports of the two laws (0.3.1). `spread_effect.py` measures on the realistic servo bench what the option `spread` costs: the same gains, applied up to 499 ticks late (described, development seeds). `prune_estimate.py` measures how many points an exact pruning of the gain recomputation would skip (a proposal, not in the law). `m7/cycle_model.sh` estimates the recomputation's cost on a Teensy 4.1 with LLVM's Cortex-M7 model, on the code arm-none-eabi-g++ produces. `m7/run_qemu.sh` replays the laws' logs on an emulated Cortex-M7 with newlib.
