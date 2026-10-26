// Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
// Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
//
// Replays a log written by rust/tools/make_replay_log.py through the C++ engine and compares every command with
// the Python engine's. From the repository root:
//   python rust/tools/make_replay_log.py
//   g++ -O2 -std=c++11 -Wall -Wextra -pedantic -Icpp cpp/tests/replay.cpp -o replay_cpp
//   ./replay_cpp replay_log.bin replay_w0.bin raw            (the 0.1.x identification)
//   ./replay_cpp replay_log_allan.bin replay_w0.bin allan    (the block-mean identification, default)
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>
#include "SyntonicEngine.h"

int main(int argc, char** argv) {
  const int N = 12, C = 4 * N;
  if (argc < 3) { std::fprintf(stderr, "usage: %s replay_log.bin replay_w0.bin [raw|allan]\n", argv[0]); return 2; }
  const bool raw = argc < 4 || std::string(argv[3]) == "raw";
  std::FILE* f = std::fopen(argv[1], "rb");
  std::FILE* g = std::fopen(argv[2], "rb");
  if (!f || !g) { std::fprintf(stderr, "cannot open the log files\n"); return 2; }
  std::vector<double> log;
  double row[C];
  while (std::fread(row, sizeof(double), C, f) == static_cast<size_t>(C)) log.insert(log.end(), row, row + C);
  double w0[N];
  if (std::fread(w0, sizeof(double), N, g) != static_cast<size_t>(N)) { std::fprintf(stderr, "bad w0 file\n"); return 2; }
  std::fclose(f); std::fclose(g);
  const size_t T = log.size() / C;
  if (T == 0) { std::fprintf(stderr, "empty log\n"); return 2; }

  static syntonic::Engine<N, 65> eng;   // static: the engine takes about 8 KB
  syntonic::Config cfg;                 // defaults: 1 kHz loop, 333 Hz frames, floor = window_floor(...)
  cfg.identification = raw ? syntonic::RAW : syntonic::ALLAN;
  eng.begin(&log[0], w0, cfg);
  double maxd = 0;
  size_t worst = 0;
  for (size_t t = 0; t < T; t++) {
    const double* r = &log[t * C];
    const double* cmd = eng.step(r, r + N, r + 2 * N);
    for (int j = 0; j < N; j++) {
      double d = std::fabs(cmd[j] - r[3 * N + j]);
      if (d > maxd) { maxd = d; worst = t; }
    }
  }
  std::printf("%zu ticks, %s identification, floor %.3f ms, max |cmd C++ - cmd Python| = %.3e rad (tick %zu)\n",
              T, raw ? "raw" : "allan", eng.floor_s() * 1e3, maxd, worst);
  return maxd < 1e-9 ? 0 : 1;
}
