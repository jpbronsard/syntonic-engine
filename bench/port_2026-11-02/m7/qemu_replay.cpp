// Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
// Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
//
// The replays of cpp/tests/replay_laws.cpp on an emulated Cortex-M7 (QEMU mps2-an500) with newlib, the C library of
// the Teensy toolchain: the loop law (two and six joints, spread 0 and 400) and the attitude law, read from the logs
// of rust/tools/make_law_replay_log.py row by row (the board has 4 MB of RAM). Built and run by run_qemu.sh. Since
// 0.3.2 the loop law runs with its defaults, the two speed-ups of the gain recomputation on, and the replay also takes
// the two-joint log of the realistic bench and reports the points evaluated and the fallbacks (none expected).
#include <stdio.h>
#include <math.h>
#include "SyntonicLoopLaw.h"
#include "RiccatiAttitudeLaw.h"
static double rel(double d, double s) { return s > 0 ? d / s : d; }
template <int N> static int loop_log(FILE* f, const double* h) {
  static syntonic::SyntonicLoopLaw<N> law; static double row[15 * 6 + 1];
  syntonic::LoopLawConfig c; c.dt = h[3]; c.frame = h[4]; c.servo_tau = h[5]; c.dead = h[6]; c.lead = h[7];
  c.spread = (int)h[8]; c.every = (int)h[9];
  double w0[N], q0[N]; fread(w0, 8, N, f); fread(q0, 8, N, f);
  if (!law.begin(q0, w0, c)) { printf("begin failed\n"); return 1; }
  long T = (long)h[2], bad = 0, pts = 0, tot = 0; double md[3] = {0, 0, 0}, sc[3] = {0, 0, 0}; bool was = false;
  for (long t = 0; t < T; t++) {
    if (fread(row, 8, 15 * N + 1, f) != (size_t)(15 * N + 1)) { printf("  short log\n"); return 1; }
    const double* o = row + 4 * N + 1;
    const double* cmd = law.step(row, row + N, row + 2 * N, row + 3 * N, row[4 * N] != 0);
    if ((c.spread == 0 && law.started && (law.tick() - 1) % c.every == 0) || (was && !law.pending())) {
      pts += law.points_evaluated(); tot += law.points_total();
    }
    was = law.pending();
    for (int j = 0; j < N; j++) {
      double v[3] = {cmd[j], law.A[j], law.gam[j]}, w[3] = {o[j], o[5 * N + j], o[6 * N + j]};
      for (int i = 0; i < 3; i++) { double d = fabs(v[i] - w[i]); if (d > md[i]) md[i] = d; if (fabs(w[i]) > sc[i]) sc[i] = fabs(w[i]); }
      if (law.gain_index(j) != (int)o[10 * N + j]) bad++;
    }
  }
  printf("  loop law N=%d spread %d: gains differ on %ld joint-ticks; cmd %.2e, A %.2e, gam %.2e (relative to scale)\n",
         N, c.spread, bad, rel(md[0], sc[0]), rel(md[1], sc[1]), rel(md[2], sc[2]));
  printf("  points evaluated %.1f %%, fallbacks %ld, late %ld\n", 100.0 * pts / (tot ? tot : 1), law.fallbacks(), law.late());
  return bad || law.fallbacks() || law.late() ? 1 : 0;
}
static int att_log(FILE* f, const double* h) {
  syntonic::AttitudeConfig c; c.dt = h[3]; c.dt_accel = h[4]; c.ground = (int)h[7]; c.settle = (int)h[8]; c.ema = h[9];
  double fx[9]; fread(fx, 8, 9, f);
  syntonic::AttitudeLevel lv[3]; for (int r = 0; r < 3; r++) { lv[r].s2 = fx[r]; lv[r].lam = fx[3 + r]; lv[r].T_rec = fx[6 + r]; }
  syntonic::RiccatiAttitudeLaw<3> law; law.begin(lv, h[5], h[6], c);
  long T = (long)h[2]; double row[28], md = 0, mp = 0, st = 0, sp = 0;
  for (long t = 0; t < T; t++) {
    fread(row, 8, 28, f);
    bool ha = row[1] == row[1];
    double th = law.step(row[0], ha, ha ? row[1] : 0.0, (int)row[2]);
    if (row[3] != 0) law.arm();
    if (fabs(th - row[4]) > md) md = fabs(th - row[4]);
    if (fabs(row[4]) > st) st = fabs(row[4]);
    for (int i = 0; i < 4; i++) for (int m = 0; m < 4; m++) {
      double d = fabs(law.P[i][m] - row[8 + 4 * i + m]); if (d > mp) mp = d;
      if (fabs(row[8 + 4 * i + m]) > sp) sp = fabs(row[8 + 4 * i + m]);
    }
  }
  printf("  attitude law: max |theta diff| %.3e rad (%.2e of scale), max |P diff| %.3e (%.2e of scale)\n", md, rel(md, st),
         mp, rel(mp, sp));
  return rel(md, st) <= 1e-12 && rel(mp, sp) <= 1e-12 ? 0 : 1;
}
int main() {
  const char* files[] = {"replay_loop_two.bin", "replay_loop_bench2.bin", "replay_loop_bench.bin", "replay_loop_bench_spread.bin",
                         "replay_att_synth.bin", "replay_att_bench.bin"};
  int bad = 0;
  for (int i = 0; i < 6; i++) {
    FILE* f = fopen(files[i], "rb"); if (!f) { printf("%s: not found, skipped\n", files[i]); continue; }
    double h[16]; fread(h, 8, 16, f); printf("%s\n", files[i]);
    if (h[0] == 1 && h[1] == 2) bad += loop_log<2>(f, h);
    else if (h[0] == 1 && h[1] == 6) bad += loop_log<6>(f, h);
    else bad += att_log(f, h);
    fclose(f);
  }
  printf("%s\n", bad ? "DIFFERS" : "agrees");
  return bad;
}
