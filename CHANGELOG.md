# Changelog

## 0.2.0 (26 Oct 2026)

- Identification on block means (the Allan-variance construction, scales 2 ms to 4.096 s, analytic noise floor n2 dt / l, constant C(H)), the new default `identification='allan'`. It fixes a limit of 0.1.x found by the software bench of 26 Oct 2026: with a noisy sensor or a slow motion, the signal falls below the read-noise floor at every lag up to 64 ms, the raw fit has no usable lag, H collapses and the window drifts by a factor of 2 to 3. With block means the engine lands within 7 % of the best fixed window in all 35 bench cells (fractional Brownian targets, H 0.30 to 0.95, optimal windows 5 to 160 ms).
- The same identification in the C++ header (`Config::identification`, `syntonic::ALLAN` / `syntonic::RAW`) and in the Rust crate (`Identification::Allan` / `Identification::Raw`), both default to block means. Replay against Python over 33 s and 12 joints: C++ 7.8e-16 rad, Rust 8.9e-16 rad with block means; 1.5e-14 and 5.4e-14 rad with the raw identification, unchanged. `make_replay_log.py` now writes both logs.
- `identification='raw'` keeps the 0.1.x arithmetic bit for bit (pinned by a unit test); `examples/servo_bench.py` and the raw replay use it, so the published tables and the earlier replay figures are unchanged. `servo_bench.py --allan` prints both identifications: in the servo regime the raw one stays a few percent better, mostly during the noise step (tip error ratios 1.02 / 0.78 / 1.07 / 0.91 against 1.00 / 0.69 / 1.07 / 0.89).
- New `C_H(H)` in the public API; new unit tests of the identification.
- README: the validated domain of each identification, and the findings of the software bench against Bayesian Kalman banks (parity in stationary and slow-drift regimes; better after noise increases; behind after noise decreases and fast drift; 25 % behind a constant-velocity Kalman on smooth targets).

## 0.1.1 (22 Oct 2026)

- Joint-error comparison stated correctly: under steady noise the syntonic engine's joint error is about 28 % larger than the tuned PID's (0.1.0 wrongly said the PID was 28 % tighter). The tip results are unchanged.
- The PID tuning grid is described precisely: 280 (P, I, D) combinations, each scored on all twelve joints.
- Teensy template: the second encoder read now waits for a new AS5048 sample (0.1 ms), so the two reads carry independent noise.
- Added Python unit tests (`python/tests`), a C++ replay test (`cpp/tests/replay.cpp`), and a README section with the commands that reproduce every check.
- Corrected the replay paths in `rust/README.md`, the memory sizes, the encoder read spacing, and the MuJoCo agreement figure in `examples/arm_fk.py`.
- A seven-second preview of MAESTRO 1 conducting at the top of the README (`docs/media`), linked to the full video attached to the v0.1.0 release.
- The README and CITATION.cff now use the Zenodo DOI for all versions, 10.5281/zenodo.23549689.

## 0.1.0 (22 Oct 2026)

First release: the engine in Python, C++ and Rust, the self-contained servo bench example, the full bench scripts with every result file, and the bench note. Archived on Zenodo as 10.5281/zenodo.23549690.
