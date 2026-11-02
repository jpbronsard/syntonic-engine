// Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
// Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.

//! The Riccati attitude law: a one-axis attitude estimator whose gains are the Riccati gains of the model that the
//! machine's measured facts describe. Nothing in it is tuned against truth. Same arithmetic as the Python reference
//! (`python/syntonic_engine/attitude_law.py`, class `RiccatiAttitudeLaw`) and the C++ header
//! (`cpp/RiccatiAttitudeLaw.h`). Ported in 0.3.1; see the Python module docstring for the law, its confirmation and
//! its limits.
//!
//! Facts (measured once). Per level r of a known input that changes the noise (on a multicopter: the throttle level),
//! from a record at rest: `s2` (the accelerometer tilt's white-noise density, rad^2 s), `lam` (the gyro's angle random
//! walk, rad^2/s), `t_rec` (the record's length, s; not used for the ground level). Also `q_b` (the gyro bias random
//! walk, rad^2/s^3), `t_pre` (the length of the arming calibration on the ground, s), `dt` and `dt_accel` (the gyro and
//! accelerometer sample periods, s). Levels are numbered 0 .. L-1 in the order given; the ground level (default 0) is
//! the one the vehicle is armed at, on the ground, at rest.
//!
//! Model. State x = (theta, b_0, ..., b_{L-1}): the angle and one gyro bias per level. The active level's bias drives
//! the integration and random-walks at q_b; the other biases are held. Every gyro sample: theta += (g - b_active) dt and
//! P <- F P F^T + Q. Every accelerometer sample a (the tilt): the Kalman update with R = max(s2(level), s2_live) /
//! dt_accel, where s2_live is the live Hadamard noise reading (block means of a - integral of g at 8, 16, 32 samples,
//! one memory per level, a second difference used only once 4 blocks of its scale have completed since the last change
//! of level). Before arming: stationary alignment (theta from the tilt alone). `arm()`: the arming covariance of the
//! biases. After arming, g is the gyro rate with the arming calibration subtracted, as the flight controller does.
//!
//! Plain arithmetic only (no exp, log or pow): it agrees with Python exactly (difference 0.0 on the replays of
//! `tools/make_law_replay_log.py`).

// index loops mirror the arithmetic of the Python class, element by element
#![allow(clippy::needless_range_loop, clippy::neg_cmp_op_on_partial_ord)]

/// The facts of one level.
#[derive(Clone, Copy, Debug)]
pub struct AttitudeLevel {
    /// The accelerometer tilt's white-noise density at this level (rad^2 s), the record's fact s2_rec.
    pub s2: f64,
    /// The gyro's angle random walk at this level (rad^2/s).
    pub lam: f64,
    /// The length of this level's record (s); not used for the ground level.
    pub t_rec: f64,
}

/// Sample periods and declared design constants of the attitude law.
#[derive(Clone, Copy, Debug)]
pub struct AttitudeConfig {
    /// Gyro sample period (s).
    pub dt: f64,
    /// Accelerometer sample period (s).
    pub dt_accel: f64,
    /// The level the vehicle is armed at.
    pub ground: usize,
    /// Noise reading: block lengths (accelerometer samples).
    pub block_scales: [u64; 3],
    /// Noise reading: averaging time (s).
    pub ema: f64,
    /// Noise reading: blocks of a scale to complete after a change of level.
    pub settle: u64,
}

impl Default for AttitudeConfig {
    fn default() -> Self {
        AttitudeConfig { dt: 1e-3, dt_accel: 4e-3, ground: 0, block_scales: [8, 16, 32], ema: 1.0, settle: 4 }
    }
}

/// The attitude law with `L` levels. The covariance is kept in blocks (angle, angle-bias, bias-angle, bias-bias);
/// `p(r, c)` reads it as the full (1 + L) x (1 + L) matrix of the Python class.
pub struct RiccatiAttitudeLaw<const L: usize> {
    /// The angle estimate (rad).
    pub theta: f64,
    /// The gyro-bias estimate per level (rad/s), relative to the arming calibration.
    pub b: [f64; L],
    /// True after `arm()`.
    pub armed: bool,
    /// Integral of the gyro (rad), for the noise reading.
    pub thr: f64,
    p00: f64,
    p0b: [f64; L],
    pb0: [f64; L],
    pbb: [[f64; L]; L],
    lv: [AttitudeLevel; L],
    ground: usize,
    q_b: f64,
    t_pre: f64,
    dt: f64,
    dta: f64,
    settle: u64,
    h: [u64; 3],
    hd: [f64; 3],
    ablk: [f64; 3],
    acc: [f64; 3],
    cnt: [u64; 3],
    nblk: [u64; 3],
    since: [u64; 3],
    p1: [f64; 3],
    p2: [f64; 3],
    dh: [[f64; 3]; L],
    cur: Option<usize>,
}

impl<const L: usize> RiccatiAttitudeLaw<L> {
    /// `levels`: the facts per level; `q_b`, `t_pre`: see the module docs. Panics if a fact is missing or out of range
    /// (as Python raises): a `t_rec` that is not positive on a level other than the ground one, a ground level that is
    /// not among the levels.
    pub fn new(levels: &[AttitudeLevel; L], q_b: f64, t_pre: f64, cfg: AttitudeConfig) -> Self {
        assert!(L >= 1, "levels: at least one level is needed");
        assert!(cfg.ground < L, "the ground level is not among the levels");
        for (r, f) in levels.iter().enumerate() {
            assert!(!f.s2.is_nan() && !f.lam.is_nan(), "a level's fact s2 or lam is missing");
            assert!(r == cfg.ground || f.t_rec > 0.0, "a level's fact t_rec is missing");
        }
        let dta = cfg.dt_accel;
        let (mut hd, mut ablk) = ([0.0; 3], [0.0; 3]);
        for i in 0..3 {
            let s = cfg.block_scales[i] as f64;
            hd[i] = s * dta;
            let e8 = 8.0 * dta * s;
            let b = dta * s / if cfg.ema > e8 { cfg.ema } else { e8 };
            ablk[i] = if b < 1.0 { b } else { 1.0 };
        }
        let mut dh = [[0.0; 3]; L];
        for r in 0..L {
            for i in 0..3 { dh[r][i] = 6.0 * levels[r].s2 / (cfg.block_scales[i] as f64 * dta); }
        }
        RiccatiAttitudeLaw {
            theta: 0.0, b: [0.0; L], armed: false, thr: 0.0,
            p00: levels[cfg.ground].s2 / dta, p0b: [0.0; L], pb0: [0.0; L], pbb: [[0.0; L]; L],
            lv: *levels, ground: cfg.ground, q_b, t_pre, dt: cfg.dt, dta, settle: cfg.settle,
            h: cfg.block_scales, hd, ablk, acc: [0.0; 3], cnt: [0; 3], nblk: [0; 3], since: [0; 3],
            p1: [0.0; 3], p2: [0.0; 3], dh, cur: None,
        }
    }

    /// The gyro-bias estimate of a level (rad/s), relative to the arming calibration.
    pub fn bias(&self, level: usize) -> f64 { self.b[level] }

    /// Component `i` of the state (theta, b_0, ..., b_{L-1}).
    pub fn x(&self, i: usize) -> f64 { if i == 0 { self.theta } else { self.b[i - 1] } }

    /// Element (r, c) of the covariance of (theta, b_0, ..., b_{L-1}).
    pub fn p(&self, r: usize, c: usize) -> f64 {
        match (r, c) {
            (0, 0) => self.p00,
            (0, c) => self.p0b[c - 1],
            (r, 0) => self.pb0[r - 1],
            (r, c) => self.pbb[r - 1][c - 1],
        }
    }

    /// The live noise reading of a level (rad^2 s), as the next update would floor it with the fact.
    pub fn noise_reading(&self, level: usize) -> f64 { self.min3(&self.dh[level]) }

    fn min3(&self, dh: &[f64; 3]) -> f64 {
        let mut m = dh[0] * self.hd[0] / 6.0;
        for i in 1..3 { let v = dh[i] * self.hd[i] / 6.0; if v < m { m = v; } }
        m
    }

    /// End of the ground calibration (the caller subtracts the calibration from g from then on).
    pub fn arm(&mut self) {
        self.armed = true;
        self.b = [0.0; L];
        self.p0b = [0.0; L];
        self.pb0 = [0.0; L];
        self.pbb = [[self.lv[self.ground].lam / self.t_pre; L]; L];
        for r in 0..L {
            if r != self.ground { self.pbb[r][r] += self.lv[r].lam / self.lv[r].t_rec; }
        }
    }

    /// One gyro sample `g` (rad/s); `a` = the tilt (rad) when an accelerometer sample arrives, else `None`; `level` =
    /// the current level (0 .. L-1; panics otherwise). Returns the angle estimate (rad).
    pub fn step(&mut self, g: f64, a: Option<f64>, level: usize) -> f64 {
        assert!(level < L, "the level is not among the levels");
        let (l, dt) = (level, self.dt);
        self.thr += g * dt;
        if self.armed {
            self.theta += (g - self.b[l]) * dt;
            let c0 = self.p0b[l];                                    // column of the active bias, before the update
            let pii = self.pbb[l][l];
            self.p00 -= dt * self.pb0[l];                            // row 0 -= dt x row i
            for m in 0..L { self.p0b[m] -= dt * self.pbb[l][m]; }
            self.p00 -= dt * c0;                                     // column 0 -= dt x (old) column i
            for r in 0..L { self.pb0[r] -= dt * self.pbb[r][l]; }
            self.p00 += dt * dt * pii;
            self.p00 += self.lv[l].lam * dt;
            self.pbb[l][l] += self.q_b * dt;
        }
        if let Some(a) = a {
            let live = self.reading(a - self.thr, l);
            let r = (if live > self.lv[l].s2 { live } else { self.lv[l].s2 }) / self.dta;
            let s = self.p00 + r;
            let e = a - self.theta;
            let k0 = self.p00 / s;
            let mut kb = [0.0; L];
            for (q, kq) in kb.iter_mut().enumerate() { *kq = self.pb0[q] / s; }
            let (r00, r0b) = (self.p00, self.p0b);                   // row 0 before the update
            self.theta += k0 * e;
            for q in 0..L { self.b[q] += kb[q] * e; }
            self.p00 -= k0 * r00;
            for m in 0..L { self.p0b[m] -= k0 * r0b[m]; }
            for q in 0..L {
                self.pb0[q] -= kb[q] * r00;
                for m in 0..L { self.pbb[q][m] -= kb[q] * r0b[m]; }
            }
        }
        self.theta
    }

    fn reading(&mut self, z: f64, level: usize) -> f64 {
        if self.cur != Some(level) { self.cur = Some(level); self.since = [0; 3]; }
        for i in 0..3 {
            self.acc[i] += z;
            self.cnt[i] += 1;
            if self.cnt[i] == self.h[i] {
                let mean = self.acc[i] / self.h[i] as f64;
                self.since[i] += 1;
                if self.nblk[i] >= 2 && self.since[i] >= self.settle {
                    let dd = mean - 2.0 * self.p1[i] + self.p2[i];
                    self.dh[level][i] += self.ablk[i] * (dd * dd - self.dh[level][i]);
                }
                self.p2[i] = self.p1[i];
                self.p1[i] = mean;
                self.acc[i] = 0.0;
                self.cnt[i] = 0;
                self.nblk[i] += 1;
            }
        }
        self.min3(&self.dh[level])
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    const DT: f64 = 1e-3;
    const QB: f64 = 9.2e-11;
    const TP: f64 = 5.0;
    const LV: [AttitudeLevel; 3] = [
        AttitudeLevel { s2: 1.9e-5, lam: 9e-9, t_rec: 0.0 },
        AttitudeLevel { s2: 2.2e-5, lam: 3.3e-8, t_rec: 40.0 },
        AttitudeLevel { s2: 3.5e-3, lam: 1.2e-6, t_rec: 8.0 },
    ];

    fn law() -> RiccatiAttitudeLaw<3> { RiccatiAttitudeLaw::new(&LV, QB, TP, AttitudeConfig::default()) }

    #[test]
    #[should_panic]
    fn missing_t_rec() {
        let mut lv = LV;
        lv[1].t_rec = 0.0;
        RiccatiAttitudeLaw::new(&lv, QB, TP, AttitudeConfig::default());
    }

    #[test]
    #[should_panic]
    fn ground_among_the_levels() {
        RiccatiAttitudeLaw::new(&LV, QB, TP, AttitudeConfig { ground: 3, ..AttitudeConfig::default() });
    }

    #[test]
    fn stationary_alignment_then_arming_covariance() {
        let mut est = law();
        for k in 0..400 { est.step(0.05, if k % 4 == 0 { Some(0.01) } else { None }, 0); }
        assert!((est.theta - 0.01).abs() < 1e-3);                                      // the tilt alone, not the gyro
        assert!((est.thr - 0.05 * 0.4).abs() < 1e-12);
        est.arm();
        for i in 0..3 {
            for m in 0..3 {
                let want = LV[0].lam / TP + if i == m && i > 0 { LV[i].lam / LV[i].t_rec } else { 0.0 };
                assert!((est.p(1 + i, 1 + m) - want).abs() <= 1e-15 * want);
            }
            assert_eq!(est.p(0, 1 + i), 0.0);
            assert_eq!(est.p(1 + i, 0), 0.0);
        }
    }

    #[test]
    fn prediction_is_f_p_ft_plus_q_and_update_is_kalman() {
        let mut est = law();
        est.arm();
        est.p00 = 2e-6; est.p0b = [1e-9, 2e-9, -1e-9]; est.pb0 = est.p0b;
        est.b = [1e-4, -2e-4, 3e-4];
        let full = |e: &RiccatiAttitudeLaw<3>| { let mut p = [[0.0; 4]; 4]; for r in 0..4 { for c in 0..4 { p[r][c] = e.p(r, c); } } p };
        let p0 = full(&est);
        let (th0, b2) = (est.theta, est.b[2]);
        est.step(0.003, None, 2);                                                       // the bias of level 2 is x[3]
        let mut fp = [[0.0; 4]; 4];
        for c in 0..4 { for r in 0..4 { fp[r][c] = p0[r][c] - if r == 0 { DT * p0[3][c] } else { 0.0 }; } }
        let mut want = [[0.0; 4]; 4];
        for r in 0..4 { for c in 0..4 { want[r][c] = fp[r][c] - if c == 0 { DT * fp[r][3] } else { 0.0 }; } }
        want[0][0] += LV[2].lam * DT;
        want[3][3] += QB * DT;
        let p1 = full(&est);
        for r in 0..4 { for c in 0..4 { assert!((p1[r][c] - want[r][c]).abs() <= 1e-12 * want[r][c].abs() + 1e-24); } }
        assert!((est.theta - (th0 + (0.003 - b2) * DT)).abs() < 1e-18);
        // the update alone, at level 1
        est.armed = false;
        let (pp, xp) = (full(&est), [est.theta, est.b[0], est.b[1], est.b[2]]);
        est.step(0.0, Some(0.004), 1);
        let live = est.noise_reading(1);
        let r = (if live > LV[1].s2 { live } else { LV[1].s2 }) / 4e-3;
        let s = pp[0][0] + r;
        for i in 0..4 {
            let k = pp[i][0] / s;
            assert!((est.x(i) - (xp[i] + k * (0.004 - xp[0]))).abs() <= 1e-12 * xp[i].abs() + 1e-18);
            for m in 0..4 { assert!((est.p(i, m) - (pp[i][m] - k * pp[0][m])).abs() <= 1e-10 * pp[i][m].abs() + 1e-30); }
        }
    }
}
