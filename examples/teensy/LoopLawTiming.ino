// Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
// Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
//
// Measures on the board what the repository could only estimate (0.3.1): the cost of the syntonic loop law on a
// Teensy 4.x for two joints (as on the ARM-1 bench), with synthetic inputs, by the cycle counter. Prints over USB serial:
// the time of begin(), of an ordinary tick, of a gain-update tick with spread = 0 (the recomputation in one tick, as in
// Python), and the longest tick with spread = 499 (the recomputation spread over the 499 ticks that follow).
// Copy cpp/SyntonicEngine.h and cpp/SyntonicLoopLaw.h next to this sketch. For the PC's arithmetic, add
// -ffp-contract=off to the compiler flags (not needed for the timing). Not yet run on hardware.
#include "SyntonicLoopLaw.h"

const int N = 2;
static syntonic::SyntonicLoopLaw<N> law;   // about 95 KB: static, never on the stack
static uint32_t rs = 2463534242u;
static double noise() { rs ^= rs << 13; rs ^= rs >> 17; rs ^= rs << 5; return (rs / 4294967296.0 - 0.5) * 6.9e-4; }   // sd ~2e-4 rad

static void run(int spread) {
  syntonic::LoopLawConfig cfg;
  cfg.dt = 0.001; cfg.frame = 0.003; cfg.servo_tau = 0.025; cfg.dead = 0.0025; cfg.spread = spread;
  double q0[N] = {0, 0}, w0[N] = {6.0, 6.0};
  uint32_t c0 = ARM_DWT_CYCCNT;
  bool ok = law.begin(q0, w0, cfg);
  double t_begin = (ARM_DWT_CYCCNT - c0) / (F_CPU / 1e6);
  double th[N] = {0, 0}, dist[N] = {0, 0}, latched[N] = {0, 0}, pipe[3][N] = {{0}};
  double sum = 0, worst = 0, sum_up = 0; long n = 0, n_up = 0;
  for (long k = 0; k < 3000; k++) {
    double t = k * 1e-3, q[N], qd[N], a[N], b[N];
    for (int j = 0; j < N; j++) {
      double f = 0.4 + 0.3 * j;
      q[j] = 0.3 * sin(2 * M_PI * f * t); qd[j] = 0.6 * M_PI * f * cos(2 * M_PI * f * t);
      dist[j] += noise();
      double y = th[j] + dist[j]; a[j] = y + noise(); b[j] = y + noise();
    }
    bool send = k % 3 == 0;
    c0 = ARM_DWT_CYCCNT;
    const double* cmd = law.step(q, qd, a, b, send);
    double us = (ARM_DWT_CYCCNT - c0) / (F_CPU / 1e6);
    bool update = spread == 0 && law.started && (law.tick() - 1) % cfg.every == 0;
    if (update) { sum_up += us; n_up++; } else { sum += us; n++; if (us > worst) worst = us; }
    for (int j = 0; j < N; j++) {                                         // a 25 ms servo with 2 ms of delay
      if (send) latched[j] = cmd[j];
      double u = pipe[k % 3][j]; pipe[k % 3][j] = latched[j];
      th[j] += 1e-3 * syntonic::clampr((u - th[j]) / 0.025, -6.0, 6.0);
    }
  }
  Serial.printf("spread %d: begin %s, %.1f ms, %d robust candidates\n", spread, ok ? "ok" : "FAILED", t_begin / 1e3,
                law.robust_count());
  Serial.printf("  ordinary ticks: mean %.1f us, longest %.1f us\n", sum / n, worst);
  if (n_up) Serial.printf("  gain-update ticks: mean %.1f ms (%ld)\n", sum_up / n_up / 1e3, n_up);
}

void setup() {
  Serial.begin(115200);
  while (!Serial && millis() < 4000) {}
  ARM_DEMCR |= ARM_DEMCR_TRCENA;                 // cycle counter on
  ARM_DWT_CTRL |= ARM_DWT_CTRL_CYCCNTENA;
  Serial.printf("syntonic loop law, %d joints, F_CPU %lu MHz\n", N, (unsigned long)(F_CPU / 1000000));
  run(0);
  run(499);
}

void loop() {}
