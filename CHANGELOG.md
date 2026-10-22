# Changelog

## 0.1.1 (22 Oct 2026)

- Joint-error comparison stated correctly: under steady noise the syntonic engine's joint error is about 28 % larger than the tuned PID's (0.1.0 wrongly said the PID was 28 % tighter). The tip results are unchanged.
- The PID tuning grid is described precisely: 280 (P, I, D) combinations, each scored on all twelve joints.
- Teensy template: the second encoder read now waits for a new AS5048 sample (0.1 ms), so the two reads carry independent noise.
- Added Python unit tests (`python/tests`), a C++ replay test (`cpp/tests/replay.cpp`), and a README section with the commands that reproduce every check.
- Corrected the replay paths in `rust/README.md`, the memory sizes, the encoder read spacing, and the MuJoCo agreement figure in `examples/arm_fk.py`.
- A seven-second preview of MAESTRO 1 conducting at the top of the README (`docs/media`).
- The README and CITATION.cff now use the Zenodo DOI for all versions, 10.5281/zenodo.23549689.

## 0.1.0 (22 Oct 2026)

First release: the engine in Python, C++ and Rust, the self-contained servo bench example, the full bench scripts with every result file, and the bench note. Archived on Zenodo as 10.5281/zenodo.23549690.
