// Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
// Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
//
// Syntonic engine, header-only C++ (Arduino / Teensy 4.x and desktop). Same arithmetic as the Python reference
// (python/syntonic_engine/engine.py) and the Rust crate (rust/): per joint, the correction window follows Theorem U' of the Syntony Principle,
//   tau* = kappa_H (sigma^2 / lambda_H)^(1 / (2H + 1)),   kappa_H = (2H Gamma(2H + 1))^(-1 / (2H + 1)),
// bounded below by the loop delay. Nothing is tuned.
//
// Two identifications of (H, lambda_H), selected by Config::identification (see the Python module docstring):
//   ALLAN (default since 0.2.0): structure function of the residual on the means of adjacent blocks, 2 ms to 4.096 s,
//         the read-noise floor n2 dt / l subtracted analytically, the fit on the finest scales that clear the floor;
//         valid whether the read noise is small (servo encoders) or large (noisy sensors).
//   RAW   (the 0.1.x engine, bit for bit): structure function of the raw residual at lags of 2 to 64 ms; valid where
//         the motion over 64 ms is far above the read noise, as on the MAESTRO 1 joint bench.
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

// block-mean (Allan) constant: E[(m2 - m1)^2] = lambda_H C(H) l^(2H) for adjacent block means of length l; C(1/2) = 2/3
static inline real C_H(real H) { return (pow(2.0, 2.0 * H + 2.0) - 4.0) / ((2.0 * H + 1.0) * (2.0 * H + 2.0)); }

// read-to-action delay: half the servo frame (zero-order hold) + half a loop tick + the sensor's delay
static inline real loop_delay(real dt, real frame, real sensor_delay = 1e-4) { return frame / 2 + dt / 2 + sensor_delay; }

// window lower bound: the larger of the loop delay and the discretisation limit dt / tau <= 0.45 (TH41-A)
static inline real window_floor(real dt, real frame, real sensor_delay = 1e-4) {
  real d = loop_delay(dt, frame, sensor_delay), t = dt / 0.45; return d > t ? d : t;
}

enum Identification { ALLAN = 0, RAW = 1 };

struct Config {
  real dt = 0.001;            // loop period (s)
  real frame = 0.003;         // servo command frame period (s)
  real servo_tau = 0.025;     // nominal servo reaction time (s); 25 ms when the data sheet gives none
  real floor = -1;            // window lower bound (s); negative = window_floor(dt, frame, sensor_delay)
  real sensor_delay = 1e-4;   // encoder delay (s), AS5048 about 0.1 ms
  real tau_max = 0.5;         // window upper bound (s)
  real tau_rate = 0.4;        // largest change of the window (s per s)
  real snr_gate = 2.0;        // a scale is used only this many noise floors above the floor
  real ema = 1.0;             // averaging time of the estimates (s)
  real warmup = 1.0;          // data collected before the window is first set (s)
  real limit = 0.08726646259971647;   // correction limit (rad), 5 degrees
  Identification identification = ALLAN;
  real lags[8] = {0.002, 0.004, 0.008, 0.016, 0.032, 0.064, 0, 0};   // RAW: structure-function lags (s)
  int nlags = 6;                                                     // RAW: lags used; ALLAN: finest scales fitted
  real scales[16] = {0.002, 0.004, 0.008, 0.016, 0.032, 0.064, 0.128, 0.256, 0.512, 1.024, 2.048, 4.096, 0, 0, 0, 0};
  int nscales = 12;                                                  // ALLAN: block lengths (s)
};

// N joints; MAXLAG = longest RAW lag in loop ticks (64 at 1 kHz, 128 at 2 kHz); unused by ALLAN but still allocated
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
    nsc_ = 0;                                                   // block lengths in ticks: at least 1, unique, ascending
    for (int i = 0; i < c.nscales && i < 16; i++) {
      long t = lround(c.scales[i] / c.dt); if (t < 1) t = 1;
      bool dup = false; for (int m = 0; m < nsc_; m++) if (sc_[m] == t) dup = true;
      if (!dup) sc_[nsc_++] = t;
    }
    for (int i = 1; i < nsc_; i++) for (int m = i; m > 0 && sc_[m] < sc_[m - 1]; m--) { long s = sc_[m]; sc_[m] = sc_[m - 1]; sc_[m - 1] = s; }
    for (int i = 0; i < nsc_; i++) {
      logs_[i] = log(sc_[i] * c.dt);
      real b = c.dt * sc_[i] / (c.ema > 8 * c.dt * sc_[i] ? c.ema : 8 * c.dt * sc_[i]);   // at least 8 blocks
      ablk_[i] = b < 1.0 ? b : 1.0;
      cnt_[i] = 0; hasprev_[i] = false;
    }
    finest_ = c.nlags;
    aw_ = c.dt / c.ema; warm_ = lround(c.warmup / c.dt); k_ = 0; head_ = 0;
    for (int j = 0; j < N; j++) {
      w0_[j] = w0[j]; cmd[j] = q0[j]; thm_[j] = q0[j]; thff_[j] = q0[j]; u[j] = 0;
      tau[j] = 0.5; H[j] = 0.5; tgt_[j] = 0.5; n2[j] = 0;
      for (int i = 0; i < nl_; i++) Dk_[i][j] = 0;
      for (int h = 0; h <= MAXLAG; h++) hist_[h][j] = 0;
      for (int i = 0; i < nsc_; i++) { acc_[i][j] = 0; prev_[i][j] = 0; DA_[i][j] = 0; nblk_[i][j] = 0; }
    }
  }

  // one tick: q = planned angles, a and b = two independent reads of every encoder (rad);
  // frame = 1 when a servo frame is sent now, 0 when not, -1 to count frames internally. Returns the command.
  const real* step(const real q[N], const real a[N], const real b[N], int frame = -1) {
    const real dt = dt_;
    const bool raw = cfg_.identification == RAW;
    if (raw) head_ = (head_ + 1) % (MAXLAG + 1);
    bool send = frame < 0 ? (k_ % frame_every_ == 0) : (frame != 0);
    bool done[16];
    if (!raw) for (int i = 0; i < nsc_; i++) { cnt_[i]++; done[i] = cnt_[i] == sc_[i]; }
    for (int j = 0; j < N; j++) {
      real meas = 0.5 * (a[j] + b[j]), e = q[j] - meas;
      real d = (a[j] - b[j]) / sqrt(2.0);
      n2[j] += aw_ * (d * d - n2[j]);
      real r = (q[j] - thff_[j]) + (meas - thm_[j]);
      real Hn = 0, lam = 0; int nok = 0; bool fitted = false;
      if (raw) {
        hist_[head_][j] = r;
        if (k_ >= maxlag_)
          for (int i = 0; i < nl_; i++) {
            real dd = r - hist_[(head_ - lag_[i] + MAXLAG + 1) % (MAXLAG + 1)][j];
            Dk_[i][j] += aw_ * (dd * dd - Dk_[i][j]);
          }
        if (k_ > warm_) {
          real D[8], fl[8];
          for (int i = 0; i < nl_; i++) { D[i] = Dk_[i][j] - n2[j]; fl[i] = n2[j]; }   // noise floor removed
          fit_(D, fl, nl_, logl_, 0, Hn, lam, nok); fitted = true;
        }
      } else {
        for (int i = 0; i < nsc_; i++) {
          acc_[i][j] += r;
          if (done[i]) {
            real mean = acc_[i][j] / sc_[i];
            if (hasprev_[i]) {
              real dd = mean - prev_[i][j];
              real aw = nblk_[i][j] > 0 ? ablk_[i] : 1.0;
              DA_[i][j] += aw * (dd * dd - DA_[i][j]); nblk_[i][j]++;
            }
            prev_[i][j] = mean; acc_[i][j] = 0;
          }
        }
        if (k_ > warm_) {
          real D[16], fl[16];
          for (int i = 0; i < nsc_; i++) {
            fl[i] = n2[j] * dt / (sc_[i] * dt);                                   // n2 dt / l, white read noise
            D[i] = nblk_[i][j] >= 2 ? DA_[i][j] - fl[i] : -1.0;
          }
          fit_(D, fl, nsc_, logs_, finest_, Hn, lam, nok); lam /= C_H(Hn); fitted = true;
        }
      }
      if (fitted) {
        real s2 = 0.5 * n2[j] * dt;
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
    if (!raw) for (int i = 0; i < nsc_; i++) if (done[i]) { cnt_[i] = 0; hasprev_[i] = true; }
    k_++;
    return cmd;
  }

  real floor_s() const { return floor_; }

 private:
  // log-log fit of the noise-corrected structure function D over the scales that clear the floor by snr_gate
  // (fallback: all positive D); with finest > 0, only the finest that many of them. Same arithmetic as engine.py.
  void fit_(const real* D, const real* fl, int m, const real* L, int finest, real& Hn, real& lam, int& nok) const {
    bool okg[16], oku[16], ok[16]; real logD[16]; int cg = 0;
    for (int i = 0; i < m; i++) { okg[i] = D[i] > cfg_.snr_gate * fl[i]; oku[i] = D[i] > 0; cg += okg[i]; }
    real Sw = 0, Sx = 0, Sy = 0, Sxx = 0, Sxy = 0; nok = 0; int used = 0;
    for (int i = 0; i < m; i++) {
      ok[i] = cg >= 2 ? okg[i] : oku[i];
      if (ok[i]) { used++; if (finest > 0 && used > finest) ok[i] = false; }
      nok += ok[i];
      logD[i] = ok[i] ? log(D[i] > 1e-30 ? D[i] : 1e-30) : 0.0;
      real W = ok[i] ? 1.0 : 0.0;
      Sw += W; Sx += W * L[i]; Sy += W * logD[i]; Sxx += W * L[i] * L[i]; Sxy += W * L[i] * logD[i];
    }
    real den = Sw * Sxx - Sx * Sx;
    real slope = den > 1e-12 ? (Sw * Sxy - Sx * Sy) / (den > 1e-12 ? den : 1e-12) : 1.0;
    Hn = clampr(0.5 * slope, 0.05, 1.0);
    real sl = 0;
    for (int i = 0; i < m; i++) if (ok[i]) sl += logD[i] - 2 * Hn * L[i];
    lam = exp(sl / (Sw > 1 ? Sw : 1));
  }

  Config cfg_;
  real dt_, floor_, aw_, w0_[N], thm_[N], thff_[N], tgt_[N], logl_[8], Dk_[8][N], hist_[MAXLAG + 1][N];
  real logs_[16], ablk_[16], acc_[16][N], prev_[16][N], DA_[16][N];
  long sc_[16], cnt_[16]; int nblk_[16][N]; bool hasprev_[16];
  int lag_[8], nl_, maxlag_, frame_every_, head_, nsc_, finest_;
  long k_, warm_;
};

}  // namespace syntonic
