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
