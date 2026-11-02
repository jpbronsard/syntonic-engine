// Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
// Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
//
// The Riccati attitude law: a one-axis attitude estimator whose gains are the Riccati gains of the model that the
// machine's measured facts describe. Nothing in it is tuned against truth. Header-only C++ (Arduino / Teensy 4.x and
// desktop), the same arithmetic as the Python reference (python/syntonic_engine/attitude_law.py, class
// RiccatiAttitudeLaw) and the Rust crate (rust/src/attitude_law.rs). Ported in 0.3.1; see the Python module docstring
// for the law, its confirmation and its limits.
//
// FACTS (measured once). Per level r of a known input that changes the noise (on a multicopter: the throttle level),
// from a record at rest: s2 (the accelerometer tilt's white-noise density, rad^2 s), lam (the gyro's angle random walk,
// rad^2/s), T_rec (the record's length, s; not used for the ground level). Also q_b (the gyro bias random walk,
// rad^2/s^3), T_pre (the length of the arming calibration on the ground, s), dt and dt_accel (the gyro and accelerometer
// sample periods, s). Levels are numbered 0 .. L-1 in the order given; the ground level (default 0) is the one the
// vehicle is armed at, on the ground, at rest.
//
// MODEL. State x = (theta, b_0, ..., b_{L-1}): the angle and one gyro bias per level. The active level's bias drives
// the integration and random-walks at q_b; the other biases are held. Every gyro sample: theta += (g - b_active) dt and
// P <- F P F^T + Q. Every accelerometer sample a (the tilt): the Kalman update with R = max(s2(level), s2_live) / dt_accel,
// where s2_live is the live Hadamard noise reading (block means of a - integral of g at 8, 16, 32 samples, one memory per
// level, a second difference used only once 4 blocks of its scale have completed since the last change of level).
// Before arming: stationary alignment (theta from the tilt alone). arm(): the arming covariance of the biases.
// After arming, g is the gyro rate with the arming calibration subtracted, as the flight controller does.
//
// Plain arithmetic only (no exp, log or pow), so it agrees with Python exactly (difference 0.0 on the replays of
// rust/tools/make_law_replay_log.py, also on an emulated Cortex-M7, bench/port_2026-11-02/m7/run_qemu.sh). With GCC,
// build with -ffp-contract=off to keep that; GCC's default for the gnu++ dialects fuses multiply-adds, which moves the
// angle by about 1e-17 rad. About 0.1 us per gyro sample on our PC.
//
// Usage (three levels: idle on the ground, low, half):
//   syntonic::AttitudeLevel lv[3] = {{s2_idle, lam_idle, 0}, {s2_low, lam_low, T_low}, {s2_half, lam_half, T_half}};
//   syntonic::RiccatiAttitudeLaw<3> att;  att.begin(lv, q_b, 5.0);    // false if a fact is missing
//   every gyro sample:  double theta = att.step(g, have_accel, tilt, level);   // level = 0, 1 or 2
//   at the end of the ground calibration:  att.arm();
#pragma once
#include <math.h>
#include "SyntonicEngine.h"   // real

namespace syntonic {

struct AttitudeLevel {
  real s2;      // the accelerometer tilt's white-noise density at this level (rad^2 s), the record's fact s2_rec
  real lam;     // the gyro's angle random walk at this level (rad^2/s)
  real T_rec;   // the length of this level's record (s); not used for the ground level
};

struct AttitudeConfig {
  real dt = 1e-3;                    // gyro sample period (s)
  real dt_accel = 4e-3;              // accelerometer sample period (s)
  int ground = 0;                    // the level the vehicle is armed at
  int block_scales[3] = {8, 16, 32};  // noise reading: block lengths (accelerometer samples)
  real ema = 1.0;                    // noise reading: averaging time (s)
  int settle = 4;                    // noise reading: blocks of a scale to complete after a change of level
};

// L levels; the state has 1 + L components
template <int L>
class RiccatiAttitudeLaw {
 public:
  static const int S = 1 + L;
  real x[S];        // (theta, b_0, ..., b_{L-1}): angle (rad) and gyro bias per level (rad/s, after the calibration)
  real P[S][S];     // covariance
  bool armed;
  real thr;         // integral of the gyro (rad), for the noise reading

  // levels: the facts per level; q_b, T_pre: see above. Returns false, and leaves the law unusable, if a fact is
  // missing or out of range (a T_rec that is not positive on a level other than the ground one, a ground level that is
  // not among the levels).
  bool begin(const AttitudeLevel levels[L], real q_b, real T_pre, const AttitudeConfig& c = AttitudeConfig()) {
    ready_ = false;
    if (L < 1 || c.ground < 0 || c.ground >= L || !(c.dt > 0) || !(c.dt_accel > 0)) return false;
    for (int r = 0; r < L; r++) {
      if (levels[r].s2 != levels[r].s2 || levels[r].lam != levels[r].lam) return false;
      if (r != c.ground && !(levels[r].T_rec > 0)) return false;
      lv_[r] = levels[r];
    }
    for (int i = 0; i < 3; i++) if (c.block_scales[i] < 1) return false;
    cfg_ = c; q_b_ = q_b; T_pre_ = T_pre; dt_ = c.dt; dta_ = c.dt_accel; ground_ = c.ground;
    for (int i = 0; i < 3; i++) {
      H_[i] = c.block_scales[i]; hd_[i] = H_[i] * dta_;
      real e8 = 8 * dta_ * H_[i], b = dta_ * H_[i] / (c.ema > e8 ? c.ema : e8);
      ablk_[i] = b < 1.0 ? b : 1.0;
      acc_[i] = 0; cnt_[i] = 0; nblk_[i] = 0; since_[i] = 0; p1_[i] = p2_[i] = 0;
    }
    for (int r = 0; r < L; r++) for (int i = 0; i < 3; i++) DH_[r][i] = 6.0 * lv_[r].s2 / (H_[i] * dta_);
    cur_ = -1; thr = 0.0;
    for (int i = 0; i < S; i++) { x[i] = 0; for (int m = 0; m < S; m++) P[i][m] = 0; }
    P[0][0] = lv_[ground_].s2 / dta_;
    armed = false; ready_ = true;
    return true;
  }

  real theta() const { return x[0]; }                      // the angle estimate (rad)
  real bias(int level) const { return x[1 + level]; }     // the gyro-bias estimate of a level (rad/s)
  bool ready() const { return ready_; }                    // begin() succeeded

  // the live noise reading of a level (rad^2 s), as the next update would floor it with the fact
  real noise_reading(int level) const { return min3_(DH_[level]); }

  // end of the ground calibration (the caller subtracts the calibration from g from then on)
  void arm() {
    armed = true;
    for (int i = 1; i < S; i++) { x[i] = 0.0; P[0][i] = 0.0; P[i][0] = 0.0; }
    for (int i = 0; i < L; i++) for (int m = 0; m < L; m++) P[1 + i][1 + m] = lv_[ground_].lam / T_pre_;
    for (int r = 0; r < L; r++) if (r != ground_) P[1 + r][1 + r] += lv_[r].lam / lv_[r].T_rec;
  }

  // one gyro sample g (rad/s); have_a = true when an accelerometer sample arrives, a = the tilt (rad); level = the
  // current level (0 .. L-1). Returns the angle estimate (rad); NaN, with nothing changed, if the level is not one of
  // the L levels.
  real step(real g, bool have_a, real a, int level) {
    if (level < 0 || level >= L) return NAN;
    const int i = 1 + level; const real dt = dt_;
    thr += g * dt;
    if (armed) {
      x[0] += (g - x[i]) * dt;
      real c[S]; for (int r = 0; r < S; r++) c[r] = P[r][i];
      const real Pii = P[i][i];
      for (int m = 0; m < S; m++) P[0][m] -= dt * P[i][m];
      for (int r = 0; r < S; r++) P[r][0] -= dt * c[r];
      P[0][0] += dt * dt * Pii;
      P[0][0] += lv_[level].lam * dt;
      P[i][i] += q_b_ * dt;
    }
    if (have_a) {
      const real live = reading_(a - thr, level);
      const real R = (live > lv_[level].s2 ? live : lv_[level].s2) / dta_;
      const real Sv = P[0][0] + R, e = a - x[0];
      real Kg[S], row0[S];
      for (int r = 0; r < S; r++) { Kg[r] = P[r][0] / Sv; row0[r] = P[0][r]; }
      for (int r = 0; r < S; r++) x[r] += Kg[r] * e;
      for (int r = 0; r < S; r++) for (int m = 0; m < S; m++) P[r][m] = P[r][m] - Kg[r] * row0[m];
    }
    return x[0];
  }

 private:
  real min3_(const real* DH) const {
    real m = DH[0] * hd_[0] / 6.0;
    for (int i = 1; i < 3; i++) { real v = DH[i] * hd_[i] / 6.0; if (v < m) m = v; }
    return m;
  }

  real reading_(real z, int level) {
    if (level != cur_) { cur_ = level; for (int i = 0; i < 3; i++) since_[i] = 0; }
    real* DH = DH_[level];
    for (int i = 0; i < 3; i++) {
      acc_[i] += z; cnt_[i]++;
      if (cnt_[i] == H_[i]) {
        real mean = acc_[i] / H_[i]; since_[i]++;
        if (nblk_[i] >= 2 && since_[i] >= cfg_.settle) {
          real dd = mean - 2 * p1_[i] + p2_[i];
          DH[i] += ablk_[i] * (dd * dd - DH[i]);
        }
        p2_[i] = p1_[i]; p1_[i] = mean; acc_[i] = 0.0; cnt_[i] = 0; nblk_[i]++;
      }
    }
    return min3_(DH);
  }

  AttitudeConfig cfg_;
  AttitudeLevel lv_[L];
  bool ready_;
  int ground_, cur_;
  real q_b_, T_pre_, dt_, dta_;
  long H_[3], cnt_[3], nblk_[3], since_[3];
  real hd_[3], ablk_[3], acc_[3], p1_[3], p2_[3];
  real DH_[L][3];
};

}  // namespace syntonic
