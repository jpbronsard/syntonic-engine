// Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
// Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
//
// The two laws of 0.3.0 in C++ (0.3.1): the loop law's facts, its gains on the robust grid in a closed loop, the
// option spread (the same gains, m ticks later); the attitude law's facts, its prediction and update against the
// textbook formulas, the arming covariance. Since 0.3.2 the two speed-ups of the gain recomputation (the same commands
// at every tick with them off and on, and with the fallback forced) and the scheduling helpers. Their agreement with
// Python is checked by replay_laws.cpp.
//   g++ -O2 -std=c++11 -Wall -Wextra -pedantic -Icpp cpp/tests/laws.cpp -o laws_cpp && ./laws_cpp
#include <cmath>
#include <chrono>
#include <cstdio>
#include <vector>
#include "SyntonicLoopLaw.h"
#include "RiccatiAttitudeLaw.h"

static const double PI = 3.141592653589793;
static unsigned long long rs = 88172645463325252ULL;
static double unif() { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return ((rs >> 11) + 0.5) / 9007199254740992.0; }
static double gauss() { return std::sqrt(-2.0 * std::log(unif())) * std::cos(2 * PI * unif()); }

// two joints, 1 kHz, 3 ms frames, a 25 ms servo with 2 ms of delay, a random-walk disturbance, two noisy reads
template <int N>
static void closed_loop(syntonic::SyntonicLoopLaw<N>& law, long ticks, int* first_index = 0, long* first_tick = 0,
                        std::vector<double>* trace = 0) {
  double th[N], dist[N], pipe[3][N], latched[N], q[N], qd[N], a[N], b[N];
  for (int j = 0; j < N; j++) { th[j] = latched[j] = dist[j] = 0; for (int h = 0; h < 3; h++) pipe[h][j] = 0; }
  if (first_tick) *first_tick = -1;
  for (long k = 0; k < ticks; k++) {
    double t = k * 1e-3;
    for (int j = 0; j < N; j++) {
      double f = 0.4 + 0.3 * j;
      q[j] = 0.3 * std::sin(2 * PI * f * t); qd[j] = 0.6 * PI * f * std::cos(2 * PI * f * t);
      dist[j] += 2e-4 * gauss() + (k == 1500 ? 0.01 : 0.0);
      double y = th[j] + dist[j]; a[j] = y + 2e-4 * gauss(); b[j] = y + 2e-4 * gauss();
    }
    bool send = k % 3 == 0;
    const double* cmd = law.step(q, qd, a, b, send);
    if (first_tick && *first_tick < 0 && law.gain_index(0) >= 0) { *first_tick = k; *first_index = law.gain_index(0); }
    if (trace) for (int j = 0; j < N; j++) trace->push_back(cmd[j]);
    for (int j = 0; j < N; j++) {
      if (send) latched[j] = cmd[j];
      double u = pipe[k % 3][j]; pipe[k % 3][j] = latched[j];          // 2 ticks of delay
      th[j] += 1e-3 * syntonic::clampr((u - th[j]) / 0.025, -6.0, 6.0);
    }
  }
}

int main() {
  int bad = 0;
  // ---- loop law: facts
  {
    static syntonic::SyntonicLoopLaw<2, 4> small;
    double q0[2] = {0, 0}, w0[2] = {6, 6};
    syntonic::LoopLawConfig c; c.dead = 0.003;                              // 3 ticks + 2 > PIPE = 4
    if (small.begin(q0, w0, c)) bad++;
    syntonic::LoopLawConfig s; s.spread = 500;                              // spread must be < every
    if (small.begin(q0, w0, s)) bad++;
    syntonic::LoopLawConfig dl; dl.delta = 1e-12;                           // below the rounding bound
    if (small.begin(q0, w0, dl)) bad++;
    syntonic::LoopLawConfig d;
    if (!small.begin(q0, w0, d) || std::fabs(small.L - (0.025 + 0.0015 + 0.0005 + 1e-4)) > 1e-15) bad++;
    if (std::fabs(small.limit(1) - 0.025 * 6) > 1e-15) bad++;
  }
  // ---- loop law: closed loop, the gains on the robust grid
  static syntonic::SyntonicLoopLaw<2> law, spread;
  syntonic::LoopLawConfig c; c.dead = 0.002;
  double q0[2] = {0, 0}, w0[2] = {6, 6};
  law.begin(q0, w0, c);
  rs = 88172645463325252ULL; int i0 = -1; long k0 = -1;
  closed_loop(law, 4000, &i0, &k0);
  if (!law.started) bad++;
  for (int j = 0; j < 2; j++) {
    int g = law.gain_index(j);
    if (g < 0 || !law.robust(g) || law.Kp[j] != law.grid_kp(g) || law.Ki[j] != law.grid_ki(g)) bad++;
    if (!(law.gam[j] >= 1 && law.gam[j] <= 3) || std::fabs(law.R[j] / 2e-8 - 1) > 0.2) bad++;   // sd^2 / 2
  }
  // ---- spread: the same first gains, exactly m ticks later
  const int ms[3] = {1, 37, 499};
  for (int i = 0; i < 3; i++) {
    syntonic::LoopLawConfig cs = c; cs.spread = ms[i];
    spread.begin(q0, w0, cs);
    rs = 88172645463325252ULL; int i1 = -1; long k1 = -1;
    closed_loop(spread, 4000, &i1, &k1);
    if (i1 != i0 || k1 != k0 + ms[i]) { bad++; std::printf("spread %d: first gains %d at tick %ld, not %d at %ld\n", ms[i], i1, k1, i0, k0 + ms[i]); }
  }
  // ---- 0.3.2: the speed-ups give the same commands at every tick: off and on, and with every choice falling back
  std::vector<double> ref;
  {
    syntonic::LoopLawConfig cr = c; cr.prune = false; cr.quad = false;    // the scan of 0.3.1
    law.begin(q0, w0, cr); rs = 88172645463325252ULL; closed_loop(law, 4000, 0, 0, &ref);
  }
  long recs = 0;
  for (int mode = 1; mode < 5; mode++) {
    syntonic::LoopLawConfig cm = c; cm.prune = mode & 1 || mode == 4; cm.quad = mode & 2 || mode == 4;
    if (mode == 4) cm.delta = 1e3;
    std::vector<double> tr;
    law.begin(q0, w0, cm); rs = 88172645463325252ULL; closed_loop(law, 4000, 0, 0, &tr);
    if (tr != ref) { bad++; std::printf("prune %d quad %d delta %g: the commands differ\n", cm.prune, cm.quad, cm.delta); }
    if (mode == 3 && (law.fallbacks() != 0 || !(law.points_evaluated() < law.points_total()))) bad++;
    if (mode == 3) recs = (4000 - k0 + c.every - 1) / c.every;            // recomputations from tick k0 on
    if (mode == 4 && law.fallbacks() != recs) { bad++; std::printf("forced: %ld fallbacks, not %ld\n", law.fallbacks(), recs); }
  }
  // ---- 0.3.2: scheduling
  {
    syntonic::LoopLawConfig cd = c; law.begin(q0, w0, cd);
    const double cand = 1500 * 70 / 600e6;                               // 70 cycles per point at 600 MHz
    const int m = law.spread_for(cand, 0.00045);
    if (!(m > 0 && m < 500) || law.spread_for(cand, 10.0) != 0 || law.spread_for(cand, 1e-5) != -1) bad++;
    if (!law.set_spread(m) || law.spread() != m || law.set_spread(500) || law.spread() != m) bad++;
    typedef std::chrono::steady_clock SC;
    double t = law.time_candidate([]() { return (unsigned long long)SC::now().time_since_epoch().count(); },
                                  (double)SC::period::den / SC::period::num);
    if (!(t > 0 && t < 0.1)) bad++;
    std::printf("loop law, scheduling: one candidate %.1f us on this machine; at 70 cycles per point on a Teensy 4.1, "
                "spread %d for 0.45 ms per tick\n", t * 1e6, m);
  }
  std::printf("loop law: %s (%d robust candidates, first gains Kp %.4g Ki %.4g at tick %ld)\n", bad ? "FAIL" : "ok",
              law.robust_count(), law.grid_kp(i0), law.grid_ki(i0), k0);

  // ---- attitude law
  int bad_a = 0;
  const double DT = 1e-3, DTA = 4e-3, QB = 9.2e-11, TP = 5.0;
  syntonic::AttitudeLevel lv[3] = {{1.9e-5, 9e-9, 0}, {2.2e-5, 3.3e-8, 40.0}, {3.5e-3, 1.2e-6, 8.0}};
  {
    syntonic::RiccatiAttitudeLaw<3> x;
    syntonic::AttitudeLevel miss[3] = {{1.9e-5, 9e-9, 0}, {2.2e-5, 3.3e-8, 0}, {3.5e-3, 1.2e-6, 8.0}};   // no T_rec
    if (x.begin(miss, QB, TP)) bad_a++;
    syntonic::AttitudeConfig g; g.ground = 3;
    if (x.begin(lv, QB, TP, g)) bad_a++;
    if (!x.begin(lv, QB, TP) || !std::isnan(x.step(0.0, false, 0.0, 3))) bad_a++;   // NaN for a level not among them
  }
  syntonic::RiccatiAttitudeLaw<3> est; est.begin(lv, QB, TP);
  for (int k = 0; k < 400; k++) est.step(0.0, k % 4 == 0, 0.0, 0);
  est.arm();
  for (int i = 0; i < 3; i++) for (int m = 0; m < 3; m++) {                                     // arming covariance
    double want = lv[0].lam / TP + (i == m && i > 0 ? lv[i].lam / lv[i].T_rec : 0.0);
    if (std::fabs(est.P[1 + i][1 + m] - want) > 1e-15 * want) bad_a++;
  }
  for (int i = 1; i < 4; i++) if (est.P[0][i] != 0 || est.P[i][0] != 0) bad_a++;
  // prediction: P <- F P F^T + Q with F = I - dt e_0 e_3^T (the bias of level 2 is x[3])
  double M[4][4], P0[4][4], x0[4];
  for (int i = 0; i < 4; i++) for (int m = 0; m < 4; m++) M[i][m] = 1e-3 * gauss();
  for (int i = 0; i < 4; i++) for (int m = 0; m < 4; m++) {
    double s = 0; for (int r = 0; r < 4; r++) s += M[i][r] * M[m][r]; est.P[i][m] = P0[i][m] = s;
  }
  for (int i = 0; i < 4; i++) est.x[i] = x0[i] = 1e-3 * gauss();
  est.step(0.003, false, 0.0, 2);
  double F[4][4] = {{1, 0, 0, -DT}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}}, FP[4][4], want[4][4];
  for (int i = 0; i < 4; i++) for (int m = 0; m < 4; m++) { FP[i][m] = 0; for (int r = 0; r < 4; r++) FP[i][m] += F[i][r] * P0[r][m]; }
  for (int i = 0; i < 4; i++) for (int m = 0; m < 4; m++) { want[i][m] = 0; for (int r = 0; r < 4; r++) want[i][m] += FP[i][r] * F[m][r]; }
  want[0][0] += lv[2].lam * DT; want[3][3] += QB * DT;
  for (int i = 0; i < 4; i++) for (int m = 0; m < 4; m++) if (std::fabs(est.P[i][m] - want[i][m]) > 1e-12 * std::fabs(want[i][m]) + 1e-24) bad_a++;
  if (std::fabs(est.x[0] - (x0[0] + (0.003 - x0[3]) * DT)) > 1e-15) bad_a++;
  // update: the Kalman update of the angle with R = max(s2, live reading) / dt_accel
  double Pp[4][4], xp[4];
  for (int i = 0; i < 4; i++) { xp[i] = est.x[i]; for (int m = 0; m < 4; m++) Pp[i][m] = est.P[i][m]; }
  syntonic::RiccatiAttitudeLaw<3> upd = est;
  upd.armed = false;                                          // no prediction: the update alone
  upd.step(0.0, true, 0.004, 1);
  double live = upd.noise_reading(1), R = (live > lv[1].s2 ? live : lv[1].s2) / DTA, S = Pp[0][0] + R;
  for (int i = 0; i < 4; i++) {
    double K = Pp[i][0] / S;
    if (std::fabs(upd.x[i] - (xp[i] + K * (0.004 - xp[0]))) > 1e-12 * std::fabs(xp[i]) + 1e-18) bad_a++;
    for (int m = 0; m < 4; m++) if (std::fabs(upd.P[i][m] - (Pp[i][m] - K * Pp[0][m])) > 1e-10 * std::fabs(Pp[i][m]) + 1e-30) bad_a++;
  }
  std::printf("attitude law: %s (facts, arming covariance, prediction, update)\n", bad_a ? "FAIL" : "ok");
  return bad || bad_a ? 1 : 0;
}
