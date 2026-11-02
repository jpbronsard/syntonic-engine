// Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
// Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
//
// The two options of 0.2.9, off by default: the servo lead (command = plan + lead x planned velocity + correction)
// and the e travel of the window (1/tau relaxes exponentially toward 1/target, 64 ms, the target taken at the first fit).
//   g++ -O2 -std=c++11 -Wall -Wextra -pedantic -Icpp cpp/tests/options.cpp -o options_cpp && ./options_cpp
#include <cmath>
#include <cstdio>
#include "SyntonicEngine.h"

int main() {
  int bad = 0;
  syntonic::Config def;
  if (def.lead != 0 || def.travel != syntonic::TRAVEL_RATE) bad++;                     // off by default

  // a joint held exactly on its plan: zero correction, command = plan + (servo_tau + loop delay) x velocity
  syntonic::Config cl; cl.lead = -1;
  syntonic::Engine<1> eng; double q0[1] = {0}, w0[1] = {10}, q[1] = {0.2}, qd[1] = {0.5};
  eng.begin(q0, w0, cl);
  double coef = 0.025 + syntonic::loop_delay(0.001, 0.003, 1e-4);
  if (std::fabs(eng.step_lead(q, qd, q, q)[0] - (0.2 + coef * 0.5)) > 1e-15) bad++;

  // e travel: after the first fit, 1/tau_new = 1/tau_old + (dt / 64 ms) (1/target - 1/tau_old)
  syntonic::Config ce; ce.travel = syntonic::TRAVEL_E; ce.identification = syntonic::RAW;
  syntonic::Engine<1> e2; e2.begin(q0, w0, ce);
  const double ae = 0.001 / 0.064;
  int checked = 0; bool first_ok = false;
  for (int k = 0; k < 3000; k++) {
    double x[1] = {0.3 * std::sin(k * 0.004)}, n = (k % 2 == 0) ? 1e-4 : -1e-4;
    double a[1] = {x[0] + n}, b[1] = {x[0] - n};
    double old = e2.tau[0];
    e2.step(x, a, b);
    if (k == 1001) first_ok = old == 0.5 && e2.tau[0] == e2.target(0);         // warm-up 1000 ticks, first fit at 1001
    if (k > 1001) {
      double want = 1.0 / old + ae * (1.0 / e2.target(0) - 1.0 / old);
      if (std::fabs(1.0 / e2.tau[0] - want) > 1e-9 * want) bad++;
      checked++;
    }
  }
  if (!first_ok || checked < 1000) bad++;
  std::printf("options (lead, e travel): %s\n", bad ? "FAIL" : "ok");
  return bad ? 1 : 0;
}
