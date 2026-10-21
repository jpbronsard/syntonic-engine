// Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
// Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
//
// Syntonic engine, header-only C++ (Arduino / Teensy 4.x and desktop). Same arithmetic as the Python reference
// (python/syntonic_engine/engine.py) and the Rust crate (rust/): per joint, the correction window follows Theorem U' of the Syntony Principle,
//   tau* = kappa_H (sigma^2 / lambda_H)^(1 / (2H + 1)),   kappa_H = (2H Gamma(2H + 1))^(-1 / (2H + 1)),
// bounded below by the loop delay. Nothing is tuned.
//
// Usage (12 joints, 1 kHz loop, 333 Hz servo frames):
//   syntonic::Engine<12> eng;  syntonic::Config cfg;  cfg.dt = 0.001; cfg.frame = 0.003;
//   eng.begin(q0, w0, cfg);                       // q0: start angles (rad), w0: no-load speeds (rad/s)
//   each tick:  const double* cmd = eng.step(q_plan, read_a, read_b);   // two reads of every encoder
//               on a servo frame, send cmd[j] to servo j
#pragma once
#include <math.h>

namespace syntonic {

typedef double real;   // Teensy 4.x has a double-precision FPU; float also works with small differences

static inline real clampr(real x, real lo, real hi) { return x < lo ? lo : (x > hi ? hi : x); }

static inline real kappa_H(real H) { return pow(2.0 * H * tgamma(2.0 * H + 1.0), -1.0 / (2.0 * H + 1.0)); }

// read-to-action delay: half the servo frame (zero-order hold) + half a loop tick + the sensor's delay
static inline real loop_delay(real dt, real frame, real sensor_delay = 1e-4) { return frame / 2 + dt / 2 + sensor_delay; }

// window lower bound: the larger of the loop delay and the discretisation limit dt / tau <= 0.45 (TH41-A)
static inline real window_floor(real dt, real frame, real sensor_delay = 1e-4) {
  real d = loop_delay(dt, frame, sensor_delay), t = dt / 0.45; return d > t ? d : t;
}

struct Config {
  real dt = 0.001;            // loop period (s)
  real frame = 0.003;         // servo command frame period (s)
  real servo_tau = 0.025;     // nominal servo reaction time (s); 25 ms when the data sheet gives none
  real floor = -1;            // window lower bound (s); negative = window_floor(dt, frame, sensor_delay)
  real sensor_delay = 1e-4;   // encoder delay (s), AS5048 about 0.1 ms
  real tau_max = 0.5;         // window upper bound (s)
  real tau_rate = 0.4;        // largest change of the window (s per s)
  real snr_gate = 2.0;        // lag used only this many noise floors above the floor
  real ema = 1.0;             // averaging time of the estimates (s)
  real warmup = 1.0;          // data collected before the window is first set (s)
  real limit = 0.08726646259971647;   // correction limit (rad), 5 degrees
  real lags[8] = {0.002, 0.004, 0.008, 0.016, 0.032, 0.064, 0, 0};   // structure-function lags (s)
  int nlags = 6;
};

// N joints; MAXLAG = longest lag in loop ticks (64 at 1 kHz, 128 at 2 kHz)
template <int N, int MAXLAG = 128>
class Engine {
 public:
  real tau[N], H[N], u[N], cmd[N], n2[N];   // window (s), smoothness, correction (rad), command (rad), noise var

  void begin(const real q0[N], const real w0[N], const Config& c) {
    cfg_ = c; dt_ = c.dt;
    frame_every_ = (int)lround(c.frame / c.dt); if (frame_every_ < 1) frame_every_ = 1;
    floor_ = c.floor >= 0 ? c.floor : window_floor(c.dt, c.frame, c.sensor_delay);
    nl_ = 0;                                                    // lags in ticks: at least 1, unique, ascending
    for (int i = 0; i < c.nlags; i++) {
      long t = lround(c.lags[i] / c.dt); if (t < 1) t = 1; if (t > MAXLAG) t = MAXLAG;
      bool dup = false; for (int m = 0; m < nl_; m++) if (lag_[m] == t) dup = true;
      if (!dup) lag_[nl_++] = (int)t;
    }
    for (int i = 1; i < nl_; i++) for (int m = i; m > 0 && lag_[m] < lag_[m - 1]; m--) { int s = lag_[m]; lag_[m] = lag_[m - 1]; lag_[m - 1] = s; }
    maxlag_ = lag_[nl_ - 1];
    for (int i = 0; i < nl_; i++) logl_[i] = log(lag_[i] * c.dt);
    aw_ = c.dt / c.ema; warm_ = lround(c.warmup / c.dt); k_ = 0; head_ = 0;
    for (int j = 0; j < N; j++) {
      w0_[j] = w0[j]; cmd[j] = q0[j]; thm_[j] = q0[j]; thff_[j] = q0[j]; u[j] = 0;
      tau[j] = 0.5; H[j] = 0.5; tgt_[j] = 0.5; n2[j] = 0;
      for (int i = 0; i < nl_; i++) Dk_[i][j] = 0;
      for (int h = 0; h <= MAXLAG; h++) hist_[h][j] = 0;
    }
  }

  // one tick: q = planned angles, a and b = two independent reads of every encoder (rad);
  // frame = 1 when a servo frame is sent now, 0 when not, -1 to count frames internally. Returns the command.
  const real* step(const real q[N], const real a[N], const real b[N], int frame = -1) {
    const real dt = dt_;
    head_ = (head_ + 1) % (MAXLAG + 1);
    bool send = frame < 0 ? (k_ % frame_every_ == 0) : (frame != 0);
    for (int j = 0; j < N; j++) {
      real meas = 0.5 * (a[j] + b[j]), e = q[j] - meas;
      real d = (a[j] - b[j]) / sqrt(2.0);
      n2[j] += aw_ * (d * d - n2[j]);
      real r = (q[j] - thff_[j]) + (meas - thm_[j]);
      hist_[head_][j] = r;
      if (k_ >= maxlag_)
        for (int i = 0; i < nl_; i++) {
          real dd = r - hist_[(head_ - lag_[i] + MAXLAG + 1) % (MAXLAG + 1)][j];
          Dk_[i][j] += aw_ * (dd * dd - Dk_[i][j]);
        }
      if (k_ > warm_) {
        real D[8], logD[8]; bool okg[8], oku[8], ok[8]; int cg = 0;
        for (int i = 0; i < nl_; i++) {
          D[i] = Dk_[i][j] - n2[j]; okg[i] = D[i] > cfg_.snr_gate * n2[j]; oku[i] = D[i] > 0; cg += okg[i];
        }
        real Sw = 0, Sx = 0, Sy = 0, Sxx = 0, Sxy = 0; int nok = 0;
        for (int i = 0; i < nl_; i++) {
          ok[i] = cg >= 2 ? okg[i] : oku[i]; nok += ok[i];
          logD[i] = ok[i] ? log(D[i] > 1e-30 ? D[i] : 1e-30) : 0.0;
          real W = ok[i] ? 1.0 : 0.0, L = logl_[i];
          Sw += W; Sx += W * L; Sy += W * logD[i]; Sxx += W * L * L; Sxy += W * L * logD[i];
        }
        real den = Sw * Sxx - Sx * Sx;
        real slope = den > 1e-12 ? (Sw * Sxy - Sx * Sy) / (den > 1e-12 ? den : 1e-12) : 1.0;
        real Hn = clampr(0.5 * slope, 0.05, 1.0);
        real s2 = 0.5 * n2[j] * dt, sl = 0;
        for (int i = 0; i < nl_; i++) if (ok[i]) sl += logD[i] - 2 * Hn * logl_[i];
        real lam = exp(sl / (Sw > 1 ? Sw : 1));
        real ts = kappa_H(Hn) * pow((s2 > 1e-30 ? s2 : 1e-30) / (lam > 1e-30 ? lam : 1e-30), 1.0 / (2 * Hn + 1));
        ts = clampr(ts, floor_, cfg_.tau_max);
        if (nok >= 2) { tgt_[j] = ts; H[j] = Hn; }
        real rr = 1.0 / tau[j], drm = cfg_.tau_rate * dt / (tau[j] * tau[j]);
        rr += clampr(1.0 / tgt_[j] - rr, -drm, drm); tau[j] = 1.0 / rr;
      }
      real g = 1.0 - exp(-dt / tau[j]);
      u[j] = clampr(u[j] + g * e, -cfg_.limit, cfg_.limit);
      if (send) cmd[j] = q[j] + u[j];
      thm_[j] += dt * clampr((cmd[j] - thm_[j]) / cfg_.servo_tau, -w0_[j], w0_[j]);
      thff_[j] += dt * clampr((q[j] - thff_[j]) / cfg_.servo_tau, -w0_[j], w0_[j]);
    }
    k_++;
    return cmd;
  }

  real floor_s() const { return floor_; }

 private:
  Config cfg_;
  real dt_, floor_, aw_, w0_[N], thm_[N], thff_[N], tgt_[N], logl_[8], Dk_[8][N], hist_[MAXLAG + 1][N];
  int lag_[8], nl_, maxlag_, frame_every_, head_;
  long k_, warm_;
};

}  // namespace syntonic
