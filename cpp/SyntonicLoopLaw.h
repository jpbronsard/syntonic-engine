// Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
// Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
//
// The syntonic loop law: a joint correction whose PI gains are not tuned but computed, every 0.5 s, as the gains that
// minimise the loop error predicted from statistics measured live. Header-only C++ (Arduino / Teensy 4.x and desktop),
// the same arithmetic as the Python reference (python/syntonic_engine/loop_law.py, class SyntonicLoopLaw) and the Rust
// crate (rust/src/loop_law.rs). Ported in 0.3.1; see the Python module docstring for the law and its confirmation.
//
// FACTS (measured once): loop period dt, servo frame period, servo time constant servo_tau (data sheet; 25 ms), no-load
// speed w0 (data sheet), and the servo's internal delay `dead` (measured once at commissioning). Design constants,
// declared: averaging time 1 s, warm-up 1 s, noise gate 2, block lags 2 to 256 ticks, at most 3 fitted lags, gain update
// every 500 ticks, phase margin >= 30 deg and max |S| <= 2, the candidate grid of (Kp, Ki) (41 x 60) and the frequency
// grid (1500 points, 0.1 to 1e4 rad/s).
//
// EACH TICK, per joint: two encoder reads give the joint angle and the read-pair noise R; the nominal servo models give
// the error the loop must cancel; its block-mean structure function gives S(u) = A u^gamma and a corner Oc. Every 500
// ticks the gains are chosen among the robust candidates as those minimising the predicted loop error J. PI: the integral
// accumulates Ki e dt every tick (frozen at the authority limit servo_tau x w0); on a frame tick the proportional term uses
// the mean error since the last frame; command = q + L x qd + clip(Kp e_mean + I, +-limit).
//
// THE GAIN RECOMPUTATION. J is evaluated for every robust candidate (about 900 to 1250 of the 2460) over the 1500
// frequencies: about 1.3 million complex evaluations per recomputation for all joints together. The frequency responses
// are not stored (that would take 30 MB); they are recomputed each time from 60 KB of per-frequency constants, so the
// law takes about 90 KB + 12 KB per joint (114 KB for two joints) and must live in a static, not on the stack.
//   LoopLawConfig::spread = 0 (default): the gains are recomputed and applied in the tick where the Python law does it.
//     This is the Python law exactly. On our PC it takes 4 to 8 ms for two joints; on a Teensy 4.1 at 600 MHz about
//     0.2 to 0.3 s by LLVM's Cortex-M7 model (bench/port_2026-11-02/m7/cycle_model.sh; not measured on the board), far
//     more than one tick of a 1 kHz loop.
//   spread = m (1 <= m < every): the statistics are taken in the same tick, the work is divided evenly over the m ticks
//     that follow (blocks of 50 frequencies), and the new gains take effect m ticks later than in the Python law. The
//     gains are the same; only when they take effect changes, and that costs: on the realistic MAESTRO 1 bench, gains
//     250 ticks late raise the tip error by about 1.5 to 2 %, 499 ticks late by about 4 to 5 %
//     (bench/port_2026-11-02/spread_effect.py). The Python reference of the option is rust/tools/make_law_replay_log.py.
//
// Numerics: double precision throughout. On the replays of rust/tools/make_law_replay_log.py the robust candidates and
// the chosen gains agree with Python exactly at every tick, and the continuous values to 2.4e-13 of their scale; also on
// an emulated Cortex-M7 with newlib (bench/port_2026-11-02/m7/run_qemu.sh). With GCC, -ffp-contract=off keeps the PC's
// arithmetic; GCC's default for the gnu++ dialects fuses multiply-adds, with differences of rounding size only.
//
// Usage (two joints, 1 kHz loop, 333 Hz servo frames):
//   static syntonic::SyntonicLoopLaw<2> law;  syntonic::LoopLawConfig cfg;  cfg.dead = 0.0025;  cfg.spread = 499;
//   law.begin(q0, w0, cfg);                     // q0: start angles (rad), w0: no-load speeds (rad/s); false if a fact is wrong
//   each tick:  const double* cmd = law.step(q_plan, qd_plan, read_a, read_b, frame_now);
//               on a frame tick, send cmd[j] to servo j
#pragma once
#include <math.h>
#include "SyntonicEngine.h"

namespace syntonic {

static const int LOOP_NW = 1500;                  // frequency grid, logspace(-1, 4, 1500) rad/s
static const int LOOP_NKP = 41;                   // Kp: 0 and geomspace(0.05, 12, 40)
static const int LOOP_NKI = 60;                   // Ki: geomspace(5, 5000, 60)
static const int LOOP_NC = LOOP_NKP * LOOP_NKI;   // candidate c = iKp x 60 + iKi
static const int LOOP_NLAG = 8;                   // block lags 2, 4, ..., 256 ticks

// block-mean constant K(gamma) of the second difference of adjacent block means (K(3) = 11/20 exact), linear between
// the tabulated points as numpy.interp
static inline real loop_K_gamma(real g) {
  static const real X[5] = {1.0, 1.5, 2.0, 2.5, 3.0};
  static const real Y[5] = {0.99999999987, 0.74370490905, 0.59538447860, 0.52771124824, 0.55};
  if (!(g > X[0])) return Y[0];
  if (!(g < X[4])) return Y[4];
  int j = 0;
  while (g >= X[j + 1]) j++;
  if (g == X[j]) return Y[j];
  real slope = (Y[j + 1] - Y[j]) / (X[j + 1] - X[j]);
  return slope * (g - X[j]) + Y[j];
}

// numpy.geomspace for positive ends: 10^(linspace(log10 a, log10 b, n)) with the ends set exactly
static inline real loop_geom(real a, real b, int n, int i) {
  if (i == 0) return a;
  if (i == n - 1) return b;
  real la = log10(a), lb = log10(b), st = (lb - la) / (n - 1);
  return pow(10.0, i * st + la);
}

struct LoopLawConfig {
  real dt = 0.001;          // loop period (s)
  real frame = 0.003;       // servo command frame period (s)
  real servo_tau = 0.025;   // servo time constant (s), data sheet; 25 ms when it gives none
  real dead = 0;            // the servo's internal delay (s), measured once at commissioning
  real lead = -1;           // lead coefficient L (s); negative = servo_tau + loop_delay(dt, frame) (Python lead=None)
  real ema = 1.0;           // averaging time of the estimates (s)
  real warmup = 1.0;        // data collected before the first fit (s)
  real gate = 2.0;          // a lag is used only this many white-noise floors above its floor
  int nfit = 3;             // at most this many fitted lags
  int every = 500;          // gain update period (ticks)
  real pm_min = 30.0;       // phase margin (deg), at least
  real ms_max = 2.0;        // max |S|, at most
  int spread = 0;           // 0: gains recomputed and applied in one tick (Python); m: spread over m ticks, applied m later
};

// N joints; PIPE = length of the nominal servo models' delay line in ticks, at least dead / dt + 2
template <int N, int PIPE = 16>
class SyntonicLoopLaw {
 public:
  real Kp[N], Ki[N];          // gains in use
  real I[N];                  // integral term (rad)
  real R[N];                  // read-pair noise variance (rad^2): EMA_1s((a - b)^2 / 2) / 2
  real A[N], gam[N], Oc[N];   // fitted structure function S(u) = A u^gam and corner (rad/s)
  real th_m[N], th_ff[N];     // nominal servo models: driven by the commands sent, by the plan feed-forward alone
  real cmd[N];                // command (rad), held between frames
  real L;                     // lead coefficient (s)
  bool started;               // true once a fit has used two lags or more on some joint

  // q0: start angles (rad), w0: no-load speeds (rad/s). Builds the frequency responses of the candidates (on a Teensy
  // 4.1 about a second). Returns false, and leaves the law unusable, if a fact is out of range.
  bool begin(const real q0[N], const real w0[N], const LoopLawConfig& c) {
    ready_ = false;
    if (!(c.dt > 0) || !(c.servo_tau > 0) || c.dead < 0 || c.every < 1 || c.nfit < 1 || c.spread < 0 ||
        (c.spread > 0 && c.spread >= c.every))
      return false;
    real ndf = c.dead / c.dt;
    nd_ = (int)floor(ndf + 1e-9); fd_ = ndf - nd_;
    if (nd_ + 2 > PIPE) return false;
    cfg_ = c; dt_ = c.dt; ts_ = c.servo_tau;
    D_ = loop_delay(c.dt, c.frame) + c.dead;
    L = c.lead < 0 ? c.servo_tau + loop_delay(c.dt, c.frame) : c.lead;
    aw_ = c.dt / c.ema; warm_ = (long)nearbyint(c.warmup / c.dt);
    for (int i = 0; i < LOOP_NLAG; i++) {
      lag_[i] = 2 << i; u_[i] = lag_[i] * c.dt; lu_[i] = log(u_[i]);
      real e8 = 8 * c.dt * lag_[i], b = c.dt * lag_[i] / (c.ema > e8 ? c.ema : e8);
      ab_[i] = b < 1.0 ? b : 1.0;
      cnt_[i] = 0; nblk_[i] = 0; nb_[i] = 0;
    }
    for (int j = 0; j < N; j++) {
      w0_[j] = w0[j]; lim_[j] = c.servo_tau * w0[j];
      th_m[j] = th_ff[j] = cmd[j] = cff_[j] = q0[j];
      n2_[j] = 0; R[j] = 1e-14; Kp[j] = 0; Ki[j] = 50.0; I[j] = 0; esum_[j] = 0;
      gam[j] = 3.0; A[j] = 0; Oc[j] = 1.0; best_[j] = -1;
      for (int i = 0; i < LOOP_NLAG; i++) acc_[i][j] = m1_[i][j] = m2_[i][j] = S2_[i][j] = 0;
      for (int h = 0; h < PIPE; h++) pc_[h][j] = pf_[h][j] = q0[j];
    }
    head_ = 0; ecnt_ = 0; k_ = 0; started = false; job_ = false;
    grid_();
    ready_ = true;
    return true;
  }

  // one tick: q, qd = planned angles (rad) and velocities (rad/s) now; a, b = two reads of every encoder (rad);
  // send = true on a servo frame tick. Returns the command (rad), which changes on frame ticks only.
  const real* step(const real q[N], const real qd[N], const real a[N], const real b[N], bool send) {
    const real dt = dt_;
    real y[N], z[N];
    for (int j = 0; j < N; j++) {
      y[j] = 0.5 * (a[j] + b[j]);
      real dn = (a[j] - b[j]) / sqrt(2.0);
      n2_[j] += aw_ * (dn * dn - n2_[j]);
      real h = 0.5 * n2_[j]; R[j] = h > 1e-14 ? h : 1e-14;
      z[j] = (q[j] - th_ff[j]) + (th_m[j] - y[j]);
    }
    identify_(z);
    if (started && k_ % cfg_.every == 0) {
      for (int j = 0; j < N; j++) { sA_[j] = A[j]; sg_[j] = gam[j]; sO_[j] = Oc[j]; r_[j] = R[j] * dt; }
      start_job_();
      if (cfg_.spread == 0) { while (work_one_()) {} apply_(); }
      else { job_ = true; job_k_ = k_; }
    } else if (job_) {
      long s = k_ - job_k_;                                      // slice s of spread, s = 1 .. spread
      long target = (work_total_ * s + cfg_.spread - 1) / cfg_.spread;
      while (work_done_ < target && work_one_()) {}
      if (s >= cfg_.spread) { while (work_one_()) {} apply_(); job_ = false; }
    }
    real e[N];
    for (int j = 0; j < N; j++) { e[j] = q[j] - y[j]; esum_[j] += e[j]; }
    ecnt_++;
    for (int j = 0; j < N; j++) I[j] += fabs(Kp[j] * e[j] + I[j]) >= lim_[j] ? 0.0 : Ki[j] * e[j] * dt;
    if (send) {
      for (int j = 0; j < N; j++) {
        real eb = esum_[j] / ecnt_; esum_[j] = 0.0;
        cmd[j] = q[j] + L * qd[j] + clampr(Kp[j] * eb + I[j], -lim_[j], lim_[j]); cff_[j] = q[j] + L * qd[j];
      }
      ecnt_ = 0;
    }
    head_ = (head_ + 1) % PIPE;
    const int h0 = (head_ - nd_ + PIPE) % PIPE, h1 = (head_ - nd_ - 1 + PIPE) % PIPE;
    for (int j = 0; j < N; j++) {
      pc_[head_][j] = cmd[j]; pf_[head_][j] = cff_[j];
      real cm = fd_ == 0 ? pc_[h0][j] : (1 - fd_) * pc_[h0][j] + fd_ * pc_[h1][j];
      real cf = fd_ == 0 ? pf_[h0][j] : (1 - fd_) * pf_[h0][j] + fd_ * pf_[h1][j];
      th_m[j] = th_m[j] + dt * clampr((cm - th_m[j]) / ts_, -w0_[j], w0_[j]);
      th_ff[j] = th_ff[j] + dt * clampr((cf - th_ff[j]) / ts_, -w0_[j], w0_[j]);
    }
    k_++;
    return cmd;
  }

  bool ready() const { return ready_; }               // begin() succeeded
  long tick() const { return k_; }                    // ticks done
  int robust_count() const { return nok_; }           // candidates passing the robustness tests
  bool robust(int c) const { return ok_[c]; }         // candidate c passes them
  real grid_kp(int c) const { return kp_[c / LOOP_NKI]; }
  real grid_ki(int c) const { return ki_[c % LOOP_NKI]; }
  int gain_index(int j) const { return best_[j]; }    // candidate in use on joint j, -1 before the first choice
  real limit(int j) const { return lim_[j]; }         // authority limit servo_tau x w0 (rad)
  real delay() const { return D_; }                   // loop delay of the model, loop_delay(dt, frame) + dead (s)
  bool pending() const { return job_; }               // spread: a recomputation is under way

 private:
  // frequency responses of the PI loop G = exp(-s D) / (servo_tau s + 1) over the candidate grid; robustness tests
  void grid_() {
    const real pi = 3.141592653589793, step = (4.0 - (-1.0)) / (LOOP_NW - 1), twopi = 2 * pi;
    real W[LOOP_NW];
    for (int w = 0; w < LOOP_NW; w++) W[w] = pow(10.0, w == LOOP_NW - 1 ? 4.0 : w * step + (-1.0));
    for (int w = 0; w < LOOP_NW; w++) {
      real dW = w == 0 ? W[1] - W[0] : (w == LOOP_NW - 1 ? W[w] - W[w - 1] : (W[w + 1] - W[w - 1]) / 2.0);
      dWf_[w] = dW / twopi * 2; W2_[w] = W[w] * W[w]; iW_[w] = 1.0 / W[w];
      real cr = cos(W[w] * D_), ci = -sin(W[w] * D_), dr = 1.0, di = W[w] * ts_;   // exp(-s D) / (1 + s ts)
      real m = dr * dr + di * di;
      Gr_[w] = (cr * dr + ci * di) / m; Gi_[w] = (ci * dr - cr * di) / m;
    }
    kp_[0] = 0.0;
    for (int i = 1; i < LOOP_NKP; i++) kp_[i] = loop_geom(0.05, 12.0, LOOP_NKP - 1, i - 1);
    for (int i = 0; i < LOOP_NKI; i++) ki_[i] = loop_geom(5.0, 5000.0, LOOP_NKI, i);
    const real r2d = 180.0 / pi, ms2 = cfg_.ms_max * cfg_.ms_max;
    nok_ = 0;
    for (int c = 0; c < LOOP_NC; c++) {
      const real kp = kp_[c / LOOP_NKI], ki = ki_[c % LOOP_NKI];
      int cross = -1; real s2max = 0, mtr = 0, L0 = 0, Lrc = 0, Lic = 0;
      for (int w = 0; w < LOOP_NW; w++) {
        real kiw = ki * iW_[w];
        real Lr = kp * Gr_[w] + kiw * Gi_[w], Li = kp * Gi_[w] - kiw * Gr_[w];
        real l2 = Lr * Lr + Li * Li, s2 = 1.0 / ((1.0 + Lr) * (1.0 + Lr) + Li * Li);
        if (w == 0) { L0 = l2; Lrc = Lr; Lic = Li; }
        if (cross < 0 && l2 < 1.0) { cross = w; Lrc = Lr; Lic = Li; }
        if (s2 > s2max) s2max = s2;
        mtr += l2 * s2 * dWf_[w];                              // |T|^2 = |L|^2 |S|^2
      }
      real pm = 180 + atan2(Lic, Lrc) * r2d;                     // at the first frequency where |L| < 1 (else the first)
      if (pm > 180) pm -= 360;
      ok_[c] = pm >= cfg_.pm_min && s2max <= ms2 && L0 > 1.0;
      if (ok_[c]) { okc_[nok_] = (short)c; MTr_[nok_] = mtr; nok_++; }
    }
    work_total_ = (long)LOOP_NW * PHI_COST * N + (long)nok_ * LOOP_NW;
  }

  // block-mean (Hadamard) identification of the error the loop must cancel: second differences of adjacent block
  // means at lags 2 .. 256 ticks, the white floor 6 R dt / u removed, S(u) = A u^gamma on the finest growing lags
  void identify_(const real z[N]) {
    for (int i = 0; i < LOOP_NLAG; i++) {
      cnt_[i]++;
      for (int j = 0; j < N; j++) acc_[i][j] += z[j];
      if (cnt_[i] == lag_[i]) {
        const bool upd = nblk_[i] >= 2;
        const real a = nb_[i] > 0 ? ab_[i] : 1.0;
        for (int j = 0; j < N; j++) {
          real mean = acc_[i][j] / lag_[i];
          if (upd) { real dd = mean - 2 * m1_[i][j] + m2_[i][j]; S2_[i][j] += a * (dd * dd - S2_[i][j]); }
          m2_[i][j] = m1_[i][j]; m1_[i][j] = mean; acc_[i][j] = 0.0;
        }
        if (upd) nb_[i]++;
        nblk_[i]++; cnt_[i] = 0;
      }
    }
    if (k_ <= warm_) return;
    bool any = false;
    for (int j = 0; j < N; j++) {
      real sig[LOOP_NLAG], ls[LOOP_NLAG]; bool valid[LOOP_NLAG], ok[LOOP_NLAG];
      int used = 0;
      for (int i = 0; i < LOOP_NLAG; i++) {
        real fl = 6.0 * R[j] * dt_ / u_[i];
        sig[i] = S2_[i][j] - fl;
        valid[i] = sig[i] > cfg_.gate * fl && nb_[i] >= 2;
        ls[i] = log(sig[i] > 1e-30 ? sig[i] : 1e-30);
        ok[i] = valid[i] && (i == 0 || ls[i] > ls[i - 1]);
        if (ok[i]) { used++; if (used > cfg_.nfit) ok[i] = false; }
      }
      real tw[LOOP_NLAG], tx[LOOP_NLAG], ty[LOOP_NLAG], txx[LOOP_NLAG], txy[LOOP_NLAG];
      for (int i = 0; i < LOOP_NLAG; i++) {
        real W = ok[i] ? 1.0 : 0.0;
        tw[i] = W; tx[i] = W * lu_[i]; ty[i] = W * ls[i]; txx[i] = W * lu_[i] * lu_[i]; txy[i] = W * lu_[i] * ls[i];
      }
      real Sw = sum8_(tw), Sx = sum8_(tx), Sy = sum8_(ty), Sxx = sum8_(txx), Sxy = sum8_(txy);
      real den = Sw * Sxx - Sx * Sx;
      real g = den > 1e-12 ? (Sw * Sxy - Sx * Sy) / (den > 1e-12 ? den : 1e-12) : 3.0;
      g = clampr(g, 1.0, 3.0);
      real lA = (Sy - g * Sx) / (Sw > 1 ? Sw : 1);
      if (Sw >= 2) { gam[j] = g; A[j] = exp(lA); any = true; }
      real ug = u_[LOOP_NLAG - 1];
      int i0 = 0;
      while (i0 < LOOP_NLAG && !valid[i0]) i0++;
      if (i0 < LOOP_NLAG) {
        while (i0 + 1 < LOOP_NLAG && valid[i0 + 1] && ls[i0 + 1] > ls[i0]) i0++;
        ug = u_[i0];
      }
      Oc[j] = 1.0 / ug;
    }
    if (any) started = true;
  }

  // the sum over the 8 lags in numpy's order: pairwise for one joint (a contiguous column), sequential otherwise
  static real sum8_(const real* t) {
    if (N == 1) return ((t[0] + t[1]) + (t[2] + t[3])) + ((t[4] + t[5]) + (t[6] + t[7]));
    real s = t[0];
    for (int i = 1; i < LOOP_NLAG; i++) s += t[i];
    return s;
  }

  // the recomputation as a list of work items of LOOP_FB frequencies each: first the spectrum (Phi of every joint),
  // then, candidate after candidate, its J of every joint, summed over the 1500 frequencies in order. The order of the
  // arithmetic does not depend on how the items are divided between ticks.
  static const int LOOP_FB = 50, LOOP_NB = LOOP_NW / LOOP_FB;   // 30 blocks of 50 frequencies
  static const int PHI_COST = 8;   // a spectrum point per joint, in (candidate, frequency) units; load balance only
  void start_job_() {
    for (int j = 0; j < N; j++) {
      C_[j] = sA_[j] / loop_K_gamma(sg_[j]); ex_[j] = -(1 + sg_[j]) / 2; o2_[j] = sO_[j] * sO_[j];
      bestJ_[j] = INFINITY; bnew_[j] = 0;
    }
    item_ = 0; work_done_ = 0;
  }

  bool work_one_() {
    if (item_ < LOOP_NB) {
      for (int w = item_ * LOOP_FB; w < (item_ + 1) * LOOP_FB; w++)
        for (int j = 0; j < N; j++) PhiW_[j][w] = C_[j] * pow(W2_[w] + o2_[j], ex_[j]) * dWf_[w];
      work_done_ += (long)LOOP_FB * PHI_COST * N;
    } else if (item_ < LOOP_NB * (1 + nok_)) {
      const int q = item_ - LOOP_NB, p = q / LOOP_NB, blk = q % LOOP_NB, c = okc_[p];
      const real kp = kp_[c / LOOP_NKI], ki = ki_[c % LOOP_NKI];
      real J[N];
      for (int j = 0; j < N; j++) J[j] = blk == 0 ? 0.0 : Jacc_[j];
      for (int w = blk * LOOP_FB; w < (blk + 1) * LOOP_FB; w++) {
        real kiw = ki * iW_[w];
        real Lr = kp * Gr_[w] + kiw * Gi_[w], Li = kp * Gi_[w] - kiw * Gr_[w];
        real s2 = 1.0 / ((1.0 + Lr) * (1.0 + Lr) + Li * Li);   // |S|^2
        for (int j = 0; j < N; j++) J[j] += PhiW_[j][w] * s2;
      }
      for (int j = 0; j < N; j++) Jacc_[j] = J[j];
      if (blk == LOOP_NB - 1) {
        for (int j = 0; j < N; j++) {
          real Jj = J[j] + r_[j] * MTr_[p];
          if (Jj < bestJ_[j] || (Jj != Jj && bestJ_[j] == bestJ_[j])) { bestJ_[j] = Jj; bnew_[j] = c; }   // first minimum
        }
      }
      work_done_ += LOOP_FB;
    } else {
      return false;
    }
    item_++;
    return true;
  }

  void apply_() {
    for (int j = 0; j < N; j++) { best_[j] = bnew_[j]; Kp[j] = grid_kp(bnew_[j]); Ki[j] = grid_ki(bnew_[j]); }
  }

  LoopLawConfig cfg_;
  bool ready_, job_;
  real dt_, ts_, D_, aw_, fd_;
  long warm_, k_, job_k_, work_total_, work_done_;
  int nd_, head_, ecnt_, nok_, item_;
  real w0_[N], lim_[N], n2_[N], esum_[N], cff_[N];
  real pc_[PIPE][N], pf_[PIPE][N];
  // identification
  int lag_[LOOP_NLAG]; long cnt_[LOOP_NLAG], nblk_[LOOP_NLAG], nb_[LOOP_NLAG];
  real u_[LOOP_NLAG], lu_[LOOP_NLAG], ab_[LOOP_NLAG];
  real acc_[LOOP_NLAG][N], m1_[LOOP_NLAG][N], m2_[LOOP_NLAG][N], S2_[LOOP_NLAG][N];
  // grid
  real W2_[LOOP_NW], iW_[LOOP_NW], Gr_[LOOP_NW], Gi_[LOOP_NW], dWf_[LOOP_NW];
  real kp_[LOOP_NKP], ki_[LOOP_NKI], MTr_[LOOP_NC];
  bool ok_[LOOP_NC]; short okc_[LOOP_NC];
  // recomputation
  real sA_[N], sg_[N], sO_[N], r_[N], C_[N], ex_[N], o2_[N], bestJ_[N], Jacc_[N];
  int best_[N], bnew_[N];
  real PhiW_[N][LOOP_NW];
};

}  // namespace syntonic
