// Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
// Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
//
// The time of the loop law's gain recomputation on this PC, before and after the speed-ups of 0.3.2 (built by
// pc_timing.sh: once against the 0.3.1 header, once against 0.3.2's with PRUNE and QUAD set). Replays a log of
// rust/tools/make_law_replay_log.py (spread 0) REPS times and times every gain-update tick; per recomputation keeps the
// shortest of the REPS runs (other jobs share the machine), then prints their mean and largest, and the mean ordinary
// tick. Usage: pc_timing LOG [REPS]
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include "SyntonicLoopLaw.h"

typedef std::chrono::steady_clock Clock;

template <int N>
static int run(const std::vector<double>& v, int reps) {
  static syntonic::SyntonicLoopLaw<N> law;
  const double* h = &v[0];
  const long T = (long)h[2];
  syntonic::LoopLawConfig cfg;
  cfg.dt = h[3]; cfg.frame = h[4]; cfg.servo_tau = h[5]; cfg.dead = h[6]; cfg.lead = h[7]; cfg.every = (int)h[9];
#ifdef PRUNE
  cfg.prune = PRUNE; cfg.quad = QUAD;
#endif
  const double *w0 = h + 16, *q0 = w0 + N, *rows = q0 + N;
  const int C = 15 * N + 1;
  std::vector<double> rec, tick(reps, 0.0);
  long nt = 0;
  for (int r = 0; r < reps; r++) {
    law.begin(q0, w0, cfg);
    size_t i = 0; nt = 0;
    for (long t = 0; t < T; t++) {
      const double* x = rows + t * C;
      Clock::time_point s = Clock::now();
      law.step(x, x + N, x + 2 * N, x + 3 * N, x[4 * N] != 0);
      double d = std::chrono::duration<double>(Clock::now() - s).count();
      if (law.started && (law.tick() - 1) % cfg.every == 0) {
        if (r == 0) rec.push_back(d); else { rec[i] = std::min(rec[i], d); }
        i++;
      } else { tick[r] += d; nt++; }
    }
  }
  double sum = 0, mx = 0;
  for (size_t i = 0; i < rec.size(); i++) { sum += rec[i]; mx = std::max(mx, rec[i]); }
  double tk = *std::min_element(tick.begin(), tick.end()) / nt;
  std::printf("%d joint(s), %d robust candidates, %zu recomputations: mean %.3f ms, largest %.3f ms; ordinary tick %.2f us\n",
              N, law.robust_count(), rec.size(), 1e3 * sum / rec.size(), 1e3 * mx, 1e6 * tk);
  return 0;
}

int main(int argc, char** argv) {
  if (argc < 2) { std::fprintf(stderr, "usage: %s LOG [REPS]\n", argv[0]); return 2; }
  int reps = argc > 2 ? std::atoi(argv[2]) : 5;
  std::FILE* f = std::fopen(argv[1], "rb");
  if (!f) return 2;
  std::vector<double> v; double buf[4096]; size_t m;
  while ((m = std::fread(buf, sizeof(double), 4096, f)) > 0) v.insert(v.end(), buf, buf + m);
  std::fclose(f);
  if (v[0] != 1 || v[8] != 0) { std::fprintf(stderr, "%s: not a loop-law log with spread 0\n", argv[1]); return 2; }
  int n = (int)v[1];
  return n == 1 ? run<1>(v, reps) : n == 2 ? run<2>(v, reps) : n == 6 ? run<6>(v, reps) : 2;
}
