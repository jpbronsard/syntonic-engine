// Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
// Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
//
// Replays a log written by rust/tools/make_law_replay_log.py through the C++ loop law or attitude law and compares
// every output with the Python law's. Prints, per output, the largest absolute difference and that difference
// relative to the output's largest magnitude over the run; for the loop law also the ticks where the chosen gains
// differ (none expected) and the time per tick. From the repository root:
//   python rust/tools/make_law_replay_log.py
//   g++ -O2 -std=c++11 -Wall -Wextra -pedantic -Icpp cpp/tests/replay_laws.cpp -o replay_laws_cpp
//   ./replay_laws_cpp replay_loop_bench.bin replay_loop_bench_spread.bin replay_loop_two.bin replay_loop_one.bin
//   ./replay_laws_cpp replay_att_synth.bin replay_att_bench.bin     (the second with --flight-bench)
// Exit status 0 when every log agrees: the same gains at every tick and every output within 1e-12 of its scale.
#include <chrono>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>
#include "SyntonicLoopLaw.h"
#include "RiccatiAttitudeLaw.h"

typedef std::chrono::steady_clock Clock;

struct Group {
  const char* name; double maxd, scale;
  Group(const char* n) : name(n), maxd(0), scale(0) {}
  void add(double got, double want) {
    double d = std::fabs(got - want); if (d > maxd) maxd = d;
    if (std::fabs(want) > scale) scale = std::fabs(want);
  }
  double rel() const { return scale > 0 ? maxd / scale : maxd; }
};

static bool report(std::vector<Group>& g) {
  bool ok = true;
  for (size_t i = 0; i < g.size(); i++) {
    std::printf("    %-8s max |C++ - Python| %.3e   relative to its scale %.3e\n", g[i].name, g[i].maxd, g[i].rel());
    if (!(g[i].rel() <= 1e-12)) ok = false;
  }
  return ok;
}

template <int N>
static bool loop_replay(const std::vector<double>& v, const char* path) {
  static syntonic::SyntonicLoopLaw<N> law;   // static: about 70 KB + 12 KB per joint
  const double* h = &v[0];
  const long T = (long)h[2];
  syntonic::LoopLawConfig cfg;
  cfg.dt = h[3]; cfg.frame = h[4]; cfg.servo_tau = h[5]; cfg.dead = h[6]; cfg.lead = h[7]; cfg.spread = (int)h[8];
  cfg.every = (int)h[9];
  const double *w0 = h + 16, *q0 = w0 + N, *rows = q0 + N;
  const int C = 15 * N + 1;
  if ((long)v.size() != 16 + 2 * N + T * C) { std::printf("%s: bad size\n", path); return false; }
  Clock::time_point t0 = Clock::now();
  if (!law.begin(q0, w0, cfg)) { std::printf("%s: begin() refused the facts\n", path); return false; }
  double t_begin = std::chrono::duration<double>(Clock::now() - t0).count();
  const char* names[10] = {"cmd", "Kp", "Ki", "I", "R", "A", "gam", "Oc", "th_m", "th_ff"};
  std::vector<Group> g; for (int i = 0; i < 10; i++) g.push_back(Group(names[i]));
  long idx_bad = 0, nrec = 0; double t_tick = 0, t_rec = 0, t_max = 0; long n_tick = 0;
  for (long t = 0; t < T; t++) {
    const double* r = rows + t * C;
    const double *q = r, *qd = r + N, *a = r + 2 * N, *b = r + 3 * N, *out = r + 4 * N + 1;
    bool send = r[4 * N] != 0;
    Clock::time_point s = Clock::now();
    const double* cmd = law.step(q, qd, a, b, send);
    double dts = std::chrono::duration<double>(Clock::now() - s).count();
    if (dts > t_max) t_max = dts;
    const bool rec = cfg.spread == 0 && law.started && (law.tick() - 1) % cfg.every == 0;   // gains recomputed now
    if (rec) { t_rec += dts; nrec++; } else { t_tick += dts; n_tick++; }
    for (int j = 0; j < N; j++) {
      g[0].add(cmd[j], out[j]); g[1].add(law.Kp[j], out[N + j]); g[2].add(law.Ki[j], out[2 * N + j]);
      g[3].add(law.I[j], out[3 * N + j]); g[4].add(law.R[j], out[4 * N + j]); g[5].add(law.A[j], out[5 * N + j]);
      g[6].add(law.gam[j], out[6 * N + j]); g[7].add(law.Oc[j], out[7 * N + j]);
      g[8].add(law.th_m[j], out[8 * N + j]); g[9].add(law.th_ff[j], out[9 * N + j]);
      if (law.gain_index(j) != (int)out[10 * N + j]) idx_bad++;
    }
  }
  std::printf("%s: loop law, %d joint(s), %ld ticks, dead %.1f ms, spread %d, %d robust candidates of %d\n", path, N, T,
              cfg.dead * 1e3, cfg.spread, law.robust_count(), syntonic::LOOP_NC);
  std::printf("    chosen gains differ from Python's on %ld joint-ticks\n", idx_bad);
  bool ok = report(g) && idx_bad == 0;
  std::printf("    time: begin %.1f ms; step %.2f us per ordinary tick", t_begin * 1e3, 1e6 * t_tick / (n_tick ? n_tick : 1));
  if (nrec) std::printf(", %.2f ms per gain-update tick (%ld)", 1e3 * t_rec / nrec, nrec);
  std::printf(", longest tick %.3f ms\n", t_max * 1e3);
  return ok;
}

template <int L>
static bool att_replay(const std::vector<double>& v, const char* path) {
  const int S = 1 + L;
  const double* h = &v[0];
  const long T = (long)h[2];
  syntonic::AttitudeConfig cfg;
  cfg.dt = h[3]; cfg.dt_accel = h[4]; cfg.ground = (int)h[7]; cfg.settle = (int)h[8]; cfg.ema = h[9];
  for (int i = 0; i < 3; i++) cfg.block_scales[i] = (int)h[10 + i];
  syntonic::AttitudeLevel lv[L];
  for (int r = 0; r < L; r++) { lv[r].s2 = h[16 + r]; lv[r].lam = h[16 + L + r]; lv[r].T_rec = h[16 + 2 * L + r]; }
  const double* rows = h + 16 + 3 * L;
  const int C = 4 + S + S * S + L + 1;
  if ((long)v.size() != 16 + 3 * L + T * C) { std::printf("%s: bad size\n", path); return false; }
  syntonic::RiccatiAttitudeLaw<L> law;
  if (!law.begin(lv, h[5], h[6], cfg)) { std::printf("%s: begin() refused the facts\n", path); return false; }
  std::vector<Group> g;
  g.push_back(Group("theta")); g.push_back(Group("biases")); g.push_back(Group("P")); g.push_back(Group("reading"));
  g.push_back(Group("thr"));
  long arms = 0;
  for (long t = 0; t < T; t++) {
    const double* r = rows + t * C;
    const bool have_a = r[1] == r[1];
    double th = law.step(r[0], have_a, have_a ? r[1] : 0.0, (int)r[2]);
    if (r[3] != 0) { law.arm(); arms++; }
    const double* out = r + 4;
    g[0].add(th, out[0]);
    for (int i = 1; i < S; i++) g[1].add(law.x[i], out[i]);
    for (int i = 0; i < S; i++) for (int m = 0; m < S; m++) g[2].add(law.P[i][m], out[S + i * S + m]);
    for (int i = 0; i < L; i++) g[3].add(law.noise_reading(i), out[S + S * S + i]);
    g[4].add(law.thr, out[S + S * S + L]);
  }
  std::printf("%s: attitude law, %d levels, %ld gyro samples, armed %ld time(s)\n", path, L, T, arms);
  return report(g);
}

int main(int argc, char** argv) {
  if (argc < 2) { std::fprintf(stderr, "usage: %s replay_loop_*.bin replay_att_*.bin ...\n", argv[0]); return 2; }
  bool all = true;
  for (int i = 1; i < argc; i++) {
    std::FILE* f = std::fopen(argv[i], "rb");
    if (!f) { std::fprintf(stderr, "cannot open %s\n", argv[i]); return 2; }
    std::vector<double> v; double buf[4096]; size_t m;
    while ((m = std::fread(buf, sizeof(double), 4096, f)) > 0) v.insert(v.end(), buf, buf + m);
    std::fclose(f);
    if (v.size() < 16) { std::fprintf(stderr, "%s: too short\n", argv[i]); return 2; }
    const int kind = (int)v[0], n = (int)v[1];
    bool ok = false;
    if (kind == 1 && n == 1) ok = loop_replay<1>(v, argv[i]);
    else if (kind == 1 && n == 2) ok = loop_replay<2>(v, argv[i]);
    else if (kind == 1 && n == 6) ok = loop_replay<6>(v, argv[i]);
    else if (kind == 2 && n == 3) ok = att_replay<3>(v, argv[i]);
    else std::printf("%s: kind %d with %d joints or levels is not compiled in this test\n", argv[i], kind, n);
    std::printf("    %s\n", ok ? "agrees" : "DIFFERS");
    all = all && ok;
  }
  return all ? 0 : 1;
}
