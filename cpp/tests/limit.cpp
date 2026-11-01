// Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
// Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
//
// The correction limit (0.2.8): by default per joint servo_tau x w0 (the authority fact); a fixed limit when set >= 0.
// A blocked joint (plan at 1 rad, readings at 0) must stop its correction at that limit.
//   g++ -O2 -std=c++11 -Wall -Wextra -pedantic -Icpp cpp/tests/limit.cpp -o limit_cpp && ./limit_cpp
#include <cmath>
#include <cstdio>
#include "SyntonicEngine.h"

static int blocked(double limit, const double w0[2], const double expect[2]) {
  syntonic::Config cfg; cfg.limit = limit;
  syntonic::Engine<2> eng; double q0[2] = {0, 0}, q[2] = {1, 1}, z[2] = {0, 0};
  eng.begin(q0, w0, cfg);
  for (int k = 0; k < 2000; k++) eng.step(q, z, z);
  int bad = 0;
  for (int j = 0; j < 2; j++) if (std::fabs(std::fabs(eng.u[j]) - expect[j]) > 1e-12) bad++;
  return bad;
}

int main() {
  const double w0[2] = {5.0, 17.5};
  const double auth[2] = {0.025 * 5.0, 0.025 * 17.5}, fixed[2] = {0.1, 0.1};
  int bad = blocked(-1, w0, auth) + blocked(0.1, w0, fixed);
  std::printf("correction limit: %s\n", bad ? "FAIL" : "ok (authority and fixed)");
  return bad ? 1 : 0;
}
