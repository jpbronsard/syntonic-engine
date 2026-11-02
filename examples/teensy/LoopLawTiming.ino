// Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
// Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
//
// Measures on the board what the repository could only estimate: the cost of the syntonic loop law on a Teensy 4.x for
// two joints (as on the ARM-1 bench), with synthetic inputs, by the cycle counter. Prints over USB serial:
//   1. for each combination of the two speed-ups of 0.3.2 (prune, quad; 0.3.1 had neither), with spread = 0: the time
//      of begin(), of an ordinary tick, and of a gain-update tick (the whole recomputation in one tick, as in Python),
//      in cycles and ms, and the share of the (candidate, frequency) points evaluated;
//   2. the time of one candidate measured by time_candidate() (with and without quad), and the spread that
//      spread_for() gives for 0.45 ms and 0.25 ms per tick;
//   3. the longest tick with the spread for 0.45 ms, and whether a recomputation ran late.
// Copy cpp/SyntonicEngine.h and cpp/SyntonicLoopLaw.h next to this sketch. The first line below gives the law the PC's
// arithmetic (GCC's -ffp-contract=off: no fused multiply-adds), as on the replays; without it GCC's default fuses them,
// which changes results by rounding only and the timing slightly. Build as usual (Teensyduino, Teensy 4.1, 600 MHz).
// Not yet run on hardware.
#pragma GCC optimize("fp-contract=off")
#include "SyntonicLoopLaw.h"

const int N = 2;
static syntonic::SyntonicLoopLaw<N> law;   // about 150 KB: static, never on the stack
static uint32_t rs = 2463534242u;
static double noise() { rs ^= rs << 13; rs ^= rs >> 17; rs ^= rs << 5; return (rs / 4294967296.0 - 0.5) * 6.9e-4; }   // sd ~2e-4 rad
static uint32_t cycles() { return ARM_DWT_CYCCNT; }
static double ms(uint32_t c) { return c * 1e3 / F_CPU_ACTUAL; }

// 5 s of a closed loop: a 25 ms servo with 2 ms of delay, a random-walk disturbance, two noisy reads per tick
static void run(bool prune, bool quad, int spread) {
  syntonic::LoopLawConfig cfg;
  cfg.dt = 0.001; cfg.frame = 0.003; cfg.servo_tau = 0.025; cfg.dead = 0.0025;
  cfg.prune = prune; cfg.quad = quad; cfg.spread = spread;
  double q0[N] = {0, 0}, w0[N] = {6.0, 6.0};
  uint32_t c0 = ARM_DWT_CYCCNT;
  bool ok = law.begin(q0, w0, cfg);
  uint32_t c_begin = ARM_DWT_CYCCNT - c0;
  rs = 2463534242u;
  double th[N] = {0, 0}, dist[N] = {0, 0}, latched[N] = {0, 0}, pipe[3][N] = {{0}};
  double sum = 0, sum_up = 0, share = 0; uint32_t worst = 0, worst_up = 0; long n = 0, n_up = 0;
  for (long k = 0; k < 5000; k++) {
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
    uint32_t c = ARM_DWT_CYCCNT - c0;
    bool update = spread == 0 && law.started && (law.tick() - 1) % cfg.every == 0;
    if (update) {
      sum_up += c; n_up++; if (c > worst_up) worst_up = c;
      share += (double)law.points_evaluated() / law.points_total();
    } else { sum += c; n++; if (c > worst) worst = c; }
    for (int j = 0; j < N; j++) {
      if (send) latched[j] = cmd[j];
      double u = pipe[k % 3][j]; pipe[k % 3][j] = latched[j];
      th[j] += 1e-3 * syntonic::clampr((u - th[j]) / 0.025, -6.0, 6.0);
    }
  }
  Serial.printf("prune %d quad %d spread %d: begin %s, %.1f ms, %d robust candidates\n", prune, quad, spread,
                ok ? "ok" : "FAILED", ms(c_begin), law.robust_count());
  Serial.printf("  ordinary ticks: mean %.1f us, longest %.1f us\n", ms(sum / n) * 1e3, ms(worst) * 1e3);
  if (n_up)
    Serial.printf("  gain-update ticks (%ld): mean %.0f cycles = %.1f ms, longest %.1f ms; points evaluated %.1f %%\n",
                  n_up, sum_up / n_up, ms(sum_up / n_up), ms(worst_up), 100 * share / n_up);
  if (spread) Serial.printf("  late recomputations %ld, fallbacks %ld\n", law.late(), law.fallbacks());
  else Serial.printf("  fallbacks %ld\n", law.fallbacks());
}

void setup() {
  Serial.begin(115200);
  while (!Serial && millis() < 4000) {}
  ARM_DEMCR |= ARM_DEMCR_TRCENA;                 // cycle counter on
  ARM_DWT_CTRL |= ARM_DWT_CTRL_CYCCNTENA;
  Serial.printf("syntonic loop law, %d joints, F_CPU_ACTUAL %lu MHz\n", N, (unsigned long)(F_CPU_ACTUAL / 1000000));
  // 1. the recomputation in one tick, with each speed-up off and on
  for (int mode = 0; mode < 4; mode++) run(mode & 1, mode & 2, 0);
  // 2. one candidate on this board, and the spread for a budget per tick (a scheduling choice, opt-in)
  int m45 = -1;
  for (int quad = 0; quad < 2; quad++) {
    syntonic::LoopLawConfig cfg; cfg.dead = 0.0025; cfg.quad = quad;
    double q0[N] = {0, 0}, w0[N] = {6.0, 6.0};
    law.begin(q0, w0, cfg);
    double s = law.time_candidate(cycles, F_CPU_ACTUAL);
    int a = law.spread_for(s, 0.00045), b = law.spread_for(s, 0.00025);
    Serial.printf("quad %d: one candidate %.1f us (%.1f cycles a point); spread for 0.45 ms a tick %d, for 0.25 ms %d "
                  "(-1: none below every)\n", quad, s * 1e6, s * F_CPU_ACTUAL / 1500, a, b);
    if (quad) m45 = a;
  }
  // 3. the spread for 0.45 ms, both speed-ups on: the longest tick
  if (m45 > 0) run(true, true, m45);
}

void loop() {}
