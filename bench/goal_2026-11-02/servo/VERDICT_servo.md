# Servo confirmation, 2 Nov 2026: syntonic loop law vs tuned PIDs, 333 Hz frames, realistic bench

Pre-registered in confirm_servo.py (sha256 af4996a0efaae1d8bcab93aeda59bc1ce8cd2b957dceb444034a173662f6f567). The file was frozen read-only before any run except the smoke test. Its frozen modules are checked at start:

- realistic_bench.py: 1fe422ea…
- rivals.py: e4bf3df3…
- syntonic_loop_law.py: a7a234cc…

Gates: G0 (hashes) and G1 (the bench reproduces examples/servo_bench.py bit for bit at the published plant) both passed exactly. The commissioning tuning reproduced the development tuning file exactly.

## Setup

- **Plant:** baton arm RJ1-RJ6 of MAESTRO 1. 1 kHz loop with 3 ms frames held (333 Hz). First-order servo, 25 ms, with a load-derated speed limit. 1° compliance at stall, 0.3° gear play. AS5048 read twice per tick. Command quantization 0.025°.
- **Commissioning:** PWM latency is 2.5 ms (the pulse falling edge at mid-travel plus 1 ms of processing).
  - The PIDs are tuned at commissioning on the published 280 grid with the edge rule, tuning seeds 61-63.
  - The law is given the measured 2.5 ms delay and nothing else. Its other facts are dt, frame, servo_tau and w0 from the data sheet.
- **Rivals:**
  - PIDLM is the PID with the same lead (27.1 ms × planned velocity) and the same frame-mean proportional term as the law. This is the primary rival.
  - PIDL is the PID with the lead only, the conventional tuned PID.
- **Circumstances in operation:** no re-tuning and no new facts.
  - C1: nominal.
  - C2: noise step, read variance ×10 from 10 to 20 s.
  - C3: bursts, ×30 for 0.5 s every 2.5 s.
  - C4: latency 3.5 ms.
  - C5: latency 1.5 ms.
  - C6: servo 35 ms.
  - C7: servo 18 ms.
  - C8: load ×1.2.
  - C9: command quantization 0.09°.
- **Runs:** fresh seeds 71-75 × excerpts A and B, 10 paired runs per circumstance. Score: baton-tip rms (mm), pooled over the 10 runs.

## Results (pooled tip rms, mm)

| | PIDLM | PIDL | LAW | LAW/PIDLM | LAW/PIDL | law wins vs PIDLM |
|---|---|---|---|---|---|---|
| C1 nominal | 0.491 | 0.565 | 0.488 | 0.994 | 0.864 | 5/10 |
| C2 noise step ×10 | 0.839 | 0.922 | 0.653 | **0.778** | 0.708 | 10/10 |
| C3 bursts ×30 | 0.916 | 1.008 | 0.813 | **0.887** | 0.806 | 10/10 |
| C4 latency 3.5 ms | 0.541 | 0.623 | 0.526 | 0.971 | 0.843 | 10/10 |
| C5 latency 1.5 ms | 0.453 | 0.522 | 0.457 | 1.009 | 0.876 | 5/10 |
| C6 servo 35 ms | 0.566 | 0.714 | 0.526 | **0.930** | 0.737 | 10/10 |
| C7 servo 18 ms | 0.568 | 0.642 | 0.551 | **0.970** | 0.858 | 10/10 |
| C8 load ×1.2 | 1.029 | 1.081 | 0.976 | **0.948** | 0.902 | 10/10 |
| C9 quantization 0.09° | 0.507 | 0.583 | 0.500 | 0.985 | 0.858 | 8/10 |

## Verdicts

- **S1 PASSES.** Parity with PIDLM, ratio ≤ 1.10, in C1-C9. The worst case is C5 at 1.009.
- **S2 PASSES.** The law beats PIDLM (≤ 1.00) where it re-measures what changed: C2, C3, C6, C7 and C8.
- **S3 PASSES.** The law beats the conventional tuned PID (≤ 1.00) in C1-C9. The worst case is C8 at 0.902.

## Described, not judged

- **Per excerpt (LAW/PIDLM):** excerpt A is the law's weaker one. C1 A 1.023 / B 0.973; C5 A 1.035 / B 0.990; C9 A 1.002 / B 0.973.
- **Jitter (tip rms of the high-frequency part):** the law has less jitter than PIDLM in every circumstance. Ratios run from 0.651 (C2) to 0.926 (C6), and 0.873 at nominal.

## Limits, stated plainly

- This is a simulation, but a realistic one: it carries every imperfection listed under Setup.
- The two cases where the law is not better (C1 at the PIDs' tuning point, C5 where the delay fact is too long by 1 ms) are parity within 1 %.
- Outside the battery, and true for every controller:
  - load ×1.5 exceeds the stall torque and is infeasible;
  - a resonant second-order servo destabilises every controller on this bench.
- An unmodeled dead time breaks the law unless the delay is measured once at commissioning. This is the README's "measure the full loop delay once", and the PIDs receive the same fact through their tuning.
