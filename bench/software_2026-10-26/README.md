# Software bench of 26 October 2026

The syntonic engine against adaptive estimators on synthetic fractional Brownian targets: six pre-registered
steps, each script carrying its hypotheses and decision rule in its docstring and printing its own SHA-256 at
the top of its results. `NOTE.md` is the bench note (version 1.1), `results/` holds every table as printed
(`*_results.txt`) and the per-seed rows (`*_rows.json`). The 110 MB bank checkpoints are not included; the
scripts recompute them.

Run from this folder, in order (numpy and scipy; the engine's identification arithmetic is copied into the
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
```

Every script is deterministic; the tables should reproduce to the printed digits. The code was written and
run inside one working session with the author; the author's own re-run is the first independent check.
