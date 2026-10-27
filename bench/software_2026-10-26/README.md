# Software bench of 26 October 2026

The syntonic engine against adaptive estimators on synthetic fractional Brownian targets (steps 0 to 4), the
quadrotor substrate in simulation (steps 5 to 6b) and cross-axis coupling on logged data of the X500 (step 7).
Each pre-registered script carries its hypotheses and decision rule in its docstring and prints its own SHA-256
at the top of its results; the hashes of the exploratory scripts of step 7 are listed in NOTE.md, section 12. `NOTE.md` is the bench note (version 1.3), `results/` holds every table as printed
(`*_results.txt`) and the per-seed rows (`*_rows.json`). The 110 MB bank checkpoints are not included; the
scripts recompute them.

Run from this folder, in order (numpy and scipy, and pandas for step 7; the engine's identification arithmetic is copied into the
smoother classes, the package is not imported):

```
python step0_uprime_engine.py            # Theorem U' and the 0.1.x identification      ~2 min
python step0b_allan_identification.py    # the block-mean identification (exploratory)  ~2 min
python step1_vs_kalman.py                # stationary, against Kalman banks            ~30 min
python step2_regime_jumps.py             # step changes of regime                      ~12 min
STEP2_SEED0=200 STEP2_OUT=step2_rep python step2_regime_jumps.py      # replication
python step3_hardened_bank.py            # banks with a change detector                ~12 min
python step3c_meta_noise.py              # the engine with a law-driven noise estimator ~1 min
python step4_drift.py                    # slow and fast drift                         ~10 min
python rerun_engine_v020.py              # engine rows with the released 0.2.0 fit      ~4 min
python step5_quad_rate_horizon.py        # the rate loop of a quadrotor cascade
python step6_quad_attitude_estimator.py  # the attitude estimator of a quadrotor         ~8 min
python evaluate_step6.py                 # its pre-registered verdict
python step6b_gyro_fact.py               # the same, gyro drift as a fact at rest       ~2 min
python evaluate_step6b.py                # its pre-registered verdict (needs step 6 rows)
python step6_diag1_structure_function.py # exploratory diagnosis of step 6 (also diag2, diag3)
python run_step7.py                      # step 7: coupling on logged X500 data (data/)  ~10 s
```

Steps 5 to 6b are the quadrotor substrate in software (NOTE.md, sections 9 to 11); step 6b reads the oracle and
competitor rows of step 6 and checks that its own recomputation of the flown filter matches them exactly.
Step 7 (NOTE.md, section 12) runs on three logged sessions of the X500, reduced in `data/` to the columns the
scripts use (no position); `run_step7.py` unpacks them and runs the five step-7 scripts unchanged.

Every script is deterministic; the tables should reproduce to the printed digits. The code was written and
run inside one working session with the author; the author's own re-run is the first independent check.
