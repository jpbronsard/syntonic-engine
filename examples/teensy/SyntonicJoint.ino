// Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
// Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
//
// Template: the syntonic engine on a Teensy 4.x for N joints, 1 kHz loop, 333 Hz servo frames.
// Fill in readEncoder(), plannedAngle() and writeServo() for your board. Not yet run on hardware.
#include "SyntonicEngine.h"   // copy cpp/SyntonicEngine.h next to this sketch

const int N = 2;                                  // e.g. ARM-1: J5 and J6
syntonic::Engine<N> eng;
IntervalTimer loopTimer;
volatile bool tick = false;

double readEncoder(int j) { (void)j; return 0.0; }                       // TODO: joint angle in rad (AS5048A over SPI)
double plannedAngle(int j, double t) { (void)j; (void)t; return 0.0; }   // TODO: the planned trajectory
void writeServo(int j, double angle) { (void)j; (void)angle; }           // TODO: angle (rad) to PWM pulse width

void onTick() { tick = true; }

void setup() {
  syntonic::Config cfg;
  cfg.dt = 0.001;              // loop period
  cfg.frame = 0.003;           // servo PWM frame (333 Hz)
  cfg.servo_tau = 0.025;       // nominal servo reaction time (data sheet, else 25 ms)
  double q0[N], w0[N];
  for (int j = 0; j < N; j++) { q0[j] = readEncoder(j); w0[j] = 10.47; }   // w0: no-load speed (rad/s), data sheet
  eng.begin(q0, w0, cfg);
  loopTimer.begin(onTick, 1000);   // 1 kHz
}

void loop() {
  if (!tick) return;
  tick = false;
  static uint32_t k = 0;
  double t = k * 0.001, q[N], a[N], b[N];
  for (int j = 0; j < N; j++) { q[j] = plannedAngle(j, t); a[j] = readEncoder(j); }
  delayMicroseconds(100);                              // the AS5048 refreshes every ~90 us: wait for a new sample
  for (int j = 0; j < N; j++) b[j] = readEncoder(j);   // second read, independent noise
  const double* cmd = eng.step(q, a, b);
  if (k % 3 == 0) for (int j = 0; j < N; j++) writeServo(j, cmd[j]);
  k++;
}
