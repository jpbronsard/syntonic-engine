// Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
// Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.

//! Syntonic engine: one formula, three facts about the machine, no tuning. A correction loop whose time window
//! comes from the Syntony Principle (Theorem U'). Same arithmetic as the Python reference (`python/syntonic_engine`) and the C++ header
//! (`cpp/SyntonicEngine.h`).
//!
//! For each joint, the window is
//!
//! ```text
//! tau* = kappa_H (sigma^2 / lambda_H)^(1 / (2H + 1)),   kappa_H = (2H Gamma(2H + 1))^(-1 / (2H + 1))
//! ```
//!
//! where sigma^2 is the encoder noise per unit time and lambda_H, H describe how strongly and how smoothly the
//! motion changes, all measured while the robot moves. The window never drops below the loop delay.
//!
//! Two identifications of (H, lambda_H), selected by `Config::identification`: `Allan` (default since 0.2.0),
//! the structure function of the residual on the means of adjacent blocks of 2 ms to 4.096 s with the read-noise
//! floor n2 dt / l subtracted analytically, valid whether the read noise is small or large; and `Raw` (the 0.1.x
//! engine, bit for bit), the structure function of the raw residual at lags of 2 to 64 ms, valid where the motion
//! over 64 ms is far above the read noise, as on the MAESTRO 1 joint bench.
//!
//! Two options, off by default (since 0.2.9; with both off the engine is that of 0.2.8, bit for bit):
//! `Config::lead`, a servo lead on the planned velocity (feed-forward, use `Engine::step_lead`), and
//! `Config::travel`, how the window travels toward the law's target (`Travel::E`: 1/tau relaxes exponentially,
//! the e-geodesic of Dual Geodesy). See the Python module docstring for what each did on the MAESTRO 1 bench.
//!
//! Since 0.3.1 the crate also holds the two laws of 0.3.0, ported from Python with the same arithmetic:
//! [`SyntonicLoopLaw`], a servo joint's PI gains recomputed every 0.5 s from the statistics measured live, and
//! [`RiccatiAttitudeLaw`], a one-axis attitude estimator with one gyro bias per throttle level. See their modules.
//!
//! The crate is `no_std` and allocation-free (math from `libm`), so it runs on microcontrollers as well as on a PC.
//!
//! ```
//! use syntonic_engine::{Config, Engine};
//! let q0 = [0.0; 2];
//! let w0 = [10.47; 2];                                  // no-load servo speed (rad/s), from the data sheet
//! let cfg = Config { dt: 0.001, frame: 0.003, ..Config::default() };
//! let mut eng: Engine<2> = Engine::new(&q0, &w0, cfg);
//! // each tick: two reads of every encoder, then
//! let cmd = eng.step(&[0.1, 0.2], &[0.0, 0.0], &[0.0, 0.0], None);
//! assert_eq!(cmd.len(), 2);
//! ```
#![cfg_attr(not(test), no_std)]

use libm::{exp, log, pow, round, sqrt, tgamma};

pub mod attitude_law;
pub mod loop_law;
pub use attitude_law::{AttitudeConfig, AttitudeLevel, RiccatiAttitudeLaw};
pub use loop_law::{loop_k_gamma, LoopLawConfig, SyntonicLoopLaw, LOOP_NC, LOOP_NKI, LOOP_NKP, LOOP_NW};

/// Kernel prefactor of Theorem U' for an exponential window on a signal of structure function lambda_H u^(2H).
pub fn kappa_h(h: f64) -> f64 {
    pow(2.0 * h * tgamma(2.0 * h + 1.0), -1.0 / (2.0 * h + 1.0))
}

/// Delay from an encoder read to the servo acting on it: half the servo frame (zero-order hold), half a loop tick
/// (age of the latest correction when an asynchronous frame takes it), and the sensor's own delay.
pub fn loop_delay(dt: f64, frame: f64, sensor_delay: f64) -> f64 {
    frame / 2.0 + dt / 2.0 + sensor_delay
}

/// Lower bound of the window: the larger of the loop delay and the discretisation limit dt / tau <= 0.45 (TH41-A).
/// 2.22 ms for a 1 kHz loop driving 333 Hz servos.
pub fn window_floor(dt: f64, frame: f64, sensor_delay: f64) -> f64 {
    let d = loop_delay(dt, frame, sensor_delay);
    let t = dt / 0.45;
    if d > t { d } else { t }
}

/// Block-mean (Allan) constant: E[(m2 - m1)^2] = lambda_H C(H) l^(2H) for adjacent block means of length l on a
/// signal of structure function lambda_H u^(2H). C(1/2) = 2/3.
pub fn c_h(h: f64) -> f64 {
    (pow(2.0, 2.0 * h + 2.0) - 4.0) / ((2.0 * h + 1.0) * (2.0 * h + 2.0))
}

#[inline]
fn clamp(x: f64, lo: f64, hi: f64) -> f64 {
    if x < lo { lo } else if x > hi { hi } else { x }
}

/// Which structure function identifies (H, lambda_H).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Identification {
    /// Block means, 2 ms to 4.096 s, analytic noise floor (default since 0.2.0).
    Allan,
    /// Raw residual at lags of 2 to 64 ms (the 0.1.x engine).
    Raw,
}

/// How the window travels toward the law's target (since 0.2.9).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Travel {
    /// 1/tau stepped toward 1/tau_target, at most `tau_rate` seconds per second (the published travel, default).
    Rate,
    /// theta = 1/tau relaxes exponentially toward 1/tau_target with time constant `travel_t` (the e-geodesic);
    /// the window takes the target at the first fit; `tau_rate` is unused.
    E,
}

/// Hardware facts and fixed design constants. Only `dt`, `frame` and `servo_tau` describe your hardware;
/// the rest are the declared defaults of the engine and are not meant to be tuned.
#[derive(Clone, Copy, Debug)]
pub struct Config {
    /// Loop period (s).
    pub dt: f64,
    /// Servo command frame period (s), e.g. 0.003 for 333 Hz PWM.
    pub frame: f64,
    /// Nominal servo reaction time (s); 0.025 when the data sheet gives none.
    pub servo_tau: f64,
    /// Window lower bound (s); `None` = `window_floor(dt, frame, sensor_delay)`.
    pub floor: Option<f64>,
    /// Encoder delay (s); AS5048 about 0.1 ms.
    pub sensor_delay: f64,
    /// Window upper bound (s).
    pub tau_max: f64,
    /// Largest change of the window (seconds per second).
    pub tau_rate: f64,
    /// A lag is used only where its structure function stands this many noise floors above the floor.
    pub snr_gate: f64,
    /// Averaging time of the noise and structure-function estimates (s).
    pub ema: f64,
    /// Data collected before the window is first set (s).
    pub warmup: f64,
    /// Correction limit (rad). Negative (the default since 0.2.8): per joint `servo_tau * w0`, the largest lag the
    /// servo itself can produce, from the data sheet's values; zero or positive: a fixed limit (5 degrees until 0.2.7).
    pub limit: f64,
    /// Which identification (see the crate docs).
    pub identification: Identification,
    /// `Raw`: structure-function lags (s); the first `nlags` are used.
    pub lags: [f64; 8],
    /// `Raw`: number of lags used. `Allan`: number of finest scales fitted among those that clear the floor.
    pub nlags: usize,
    /// `Allan`: block lengths (s); the first `nscales` are used.
    pub scales: [f64; 16],
    /// `Allan`: number of block lengths used.
    pub nscales: usize,
    /// Servo lead (s), off by default (0, since 0.2.9). Negative: `servo_tau + loop_delay(dt, frame, sensor_delay)`,
    /// the servo's lag plus the read-to-action delay; positive: a fixed lead. When on, call `step_lead` with the
    /// planned velocity; the engine's servo and plan models see the same lead.
    pub lead: f64,
    /// How the window travels toward the law's target (since 0.2.9).
    pub travel: Travel,
    /// `Travel::E`: relaxation time constant (s); `None` = the longest of the `lags` used, 64 ms by default.
    pub travel_t: Option<f64>,
}

impl Default for Config {
    fn default() -> Self {
        Config {
            dt: 0.001,
            frame: 0.003,
            servo_tau: 0.025,
            floor: None,
            sensor_delay: 1e-4,
            tau_max: 0.5,
            tau_rate: 0.4,
            snr_gate: 2.0,
            ema: 1.0,
            warmup: 1.0,
            limit: -1.0,
            identification: Identification::Allan,
            lags: [0.002, 0.004, 0.008, 0.016, 0.032, 0.064, 0.0, 0.0],
            nlags: 6,
            scales: [0.002, 0.004, 0.008, 0.016, 0.032, 0.064, 0.128, 0.256, 0.512, 1.024, 2.048, 4.096,
                     0.0, 0.0, 0.0, 0.0],
            nscales: 12,
            lead: 0.0,
            travel: Travel::Rate,
            travel_t: None,
        }
    }
}

/// The engine for `N` joints. `RING` is the length of the residual history in loop ticks and must be larger than
/// the longest lag in ticks (65 is enough at 1 kHz with the default lags, 129 at 2 kHz).
pub struct Engine<const N: usize, const RING: usize = 129> {
    /// Window of each joint (s).
    pub tau: [f64; N],
    /// Measured smoothness of each joint's residual.
    pub h: [f64; N],
    /// Correction of each joint (rad).
    pub u: [f64; N],
    /// Command of each joint (rad), held between servo frames.
    pub cmd: [f64; N],
    /// Per-read encoder noise variance of each joint (rad^2).
    pub n2: [f64; N],
    cfg: Config,
    floor: f64,
    aw: f64,
    w0: [f64; N],
    lim: [f64; N],
    thm: [f64; N],
    thff: [f64; N],
    tgt: [f64; N],
    logl: [f64; 8],
    lag: [usize; 8],
    nl: usize,
    maxlag: usize,
    frame_every: u64,
    dk: [[f64; N]; 8],
    hist: [[f64; N]; RING],
    head: usize,
    k: u64,
    warm: u64,
    // block-mean identification
    sc: [usize; 16],
    nsc: usize,
    logs: [f64; 16],
    ablk: [f64; 16],
    cnt: [usize; 16],
    hasprev: [bool; 16],
    acc: [[f64; N]; 16],
    prev: [[f64; N]; 16],
    da: [[f64; N]; 16],
    nblk: [[u32; N]; 16],
    finest: usize,
    // options (0.2.9)
    lead: f64,
    a_e: f64,
    started: bool,
}

impl<const N: usize, const RING: usize> Engine<N, RING> {
    /// `q0`: initial joint positions (rad); `w0`: no-load servo speeds (rad/s), from the data sheet.
    pub fn new(q0: &[f64; N], w0: &[f64; N], cfg: Config) -> Self {
        let dt = cfg.dt;
        let fe = round(cfg.frame / dt);
        let frame_every = if fe < 1.0 { 1 } else { fe as u64 };
        let floor = match cfg.floor {
            Some(f) => f,
            None => window_floor(dt, cfg.frame, cfg.sensor_delay),
        };
        // lags in ticks: at least one, unique, ascending
        let mut lag = [0usize; 8];
        let mut nl = 0usize;
        for i in 0..cfg.nlags.min(8) {
            let mut t = round(cfg.lags[i] / dt);
            if t < 1.0 { t = 1.0; }
            let t = t as usize;
            assert!(t < RING, "RING must exceed the longest lag in ticks");
            if !lag[..nl].contains(&t) { lag[nl] = t; nl += 1; }
        }
        assert!(nl >= 2, "at least two distinct lags are needed");
        lag[..nl].sort_unstable();
        let mut logl = [0.0; 8];
        for i in 0..nl { logl[i] = log(lag[i] as f64 * dt); }
        // block lengths in ticks: at least one, unique, ascending
        let mut sc = [0usize; 16];
        let mut nsc = 0usize;
        for i in 0..cfg.nscales.min(16) {
            let mut t = round(cfg.scales[i] / dt);
            if t < 1.0 { t = 1.0; }
            let t = t as usize;
            if !sc[..nsc].contains(&t) { sc[nsc] = t; nsc += 1; }
        }
        assert!(nsc >= 2, "at least two distinct scales are needed");
        sc[..nsc].sort_unstable();
        let mut logs = [0.0; 16];
        let mut ablk = [0.0; 16];
        for i in 0..nsc {
            logs[i] = log(sc[i] as f64 * dt);
            let eight = 8.0 * dt * sc[i] as f64;
            let b = dt * sc[i] as f64 / if cfg.ema > eight { cfg.ema } else { eight };   // at least 8 blocks
            ablk[i] = if b < 1.0 { b } else { 1.0 };
        }
        let mut lim = [0.0; N];
        for j in 0..N { lim[j] = if cfg.limit >= 0.0 { cfg.limit } else { cfg.servo_tau * w0[j] }; }
        let lead = if cfg.lead < 0.0 { cfg.servo_tau + loop_delay(dt, cfg.frame, cfg.sensor_delay) } else { cfg.lead };
        let travel_t = match cfg.travel_t {
            Some(t) => t,
            None => cfg.lags[..cfg.nlags.min(8)].iter().fold(0.0f64, |m, &x| if x > m { x } else { m }),
        };
        assert!(travel_t > 0.0, "travel_t must be positive");
        Engine {
            tau: [0.5; N], h: [0.5; N], u: [0.0; N], cmd: *q0, n2: [0.0; N],
            cfg, floor, aw: dt / cfg.ema, w0: *w0, lim, thm: *q0, thff: *q0, tgt: [0.5; N],
            logl, lag, nl, maxlag: lag[nl - 1], frame_every,
            dk: [[0.0; N]; 8], hist: [[0.0; N]; RING], head: 0, k: 0,
            warm: round(cfg.warmup / dt) as u64,
            sc, nsc, logs, ablk, cnt: [0; 16], hasprev: [false; 16],
            acc: [[0.0; N]; 16], prev: [[0.0; N]; 16], da: [[0.0; N]; 16], nblk: [[0; N]; 16],
            finest: cfg.nlags,
            lead, a_e: dt / travel_t, started: false,
        }
    }

    /// Log-log fit of the noise-corrected structure function `d` over the scales that clear the floor `fl` by
    /// `snr_gate` (fallback: all positive); with `finest > 0`, only the finest that many of them. Returns
    /// (H, lambda before the C(H) correction, number of scales used). Same arithmetic as engine.py.
    fn fit(&self, d: &[f64], fl: &[f64], l: &[f64], finest: usize) -> (f64, f64, usize) {
        let m = d.len();
        let mut okg = [false; 16];
        let mut oku = [false; 16];
        let mut cg = 0;
        for i in 0..m {
            okg[i] = d[i] > self.cfg.snr_gate * fl[i];
            oku[i] = d[i] > 0.0;
            if okg[i] { cg += 1; }
        }
        let (mut sw, mut sx, mut sy, mut sxx, mut sxy) = (0.0, 0.0, 0.0, 0.0, 0.0);
        let mut ok = [false; 16];
        let mut logd = [0.0; 16];
        let mut nok = 0;
        let mut used = 0;
        for i in 0..m {
            ok[i] = if cg >= 2 { okg[i] } else { oku[i] };
            if ok[i] { used += 1; if finest > 0 && used > finest { ok[i] = false; } }
            if ok[i] { nok += 1; }
            logd[i] = if ok[i] { log(if d[i] > 1e-30 { d[i] } else { 1e-30 }) } else { 0.0 };
            let w = if ok[i] { 1.0 } else { 0.0 };
            sw += w; sx += w * l[i]; sy += w * logd[i]; sxx += w * l[i] * l[i]; sxy += w * l[i] * logd[i];
        }
        let den = sw * sxx - sx * sx;
        let slope = if den > 1e-12 { (sw * sxy - sx * sy) / den } else { 1.0 };
        let hn = clamp(0.5 * slope, 0.05, 1.0);
        let mut sl = 0.0;
        for i in 0..m { if ok[i] { sl += logd[i] - 2.0 * hn * l[i]; } }
        let lam = exp(sl / if sw > 1.0 { sw } else { 1.0 });
        (hn, lam, nok)
    }

    /// The window lower bound in use (s).
    pub fn floor(&self) -> f64 { self.floor }

    /// One loop tick. `q`: planned angles now (rad); `a`, `b`: two independent reads of every encoder (rad);
    /// `frame`: `Some(true)` when a servo frame is sent now, `Some(false)` when not, `None` to count frames
    /// internally (one every frame period, starting with the first tick). Returns the command to send (rad).
    /// Panics if the lead is on: use `step_lead`.
    pub fn step(&mut self, q: &[f64; N], a: &[f64; N], b: &[f64; N], frame: Option<bool>) -> &[f64; N] {
        assert!(self.lead == 0.0, "the lead needs step_lead with the planned velocity");
        self.tick(q, None, a, b, frame)
    }

    /// One loop tick with the servo lead (since 0.2.9): as `step`, with `qd` the planned velocities now (rad/s).
    pub fn step_lead(&mut self, q: &[f64; N], qd: &[f64; N], a: &[f64; N], b: &[f64; N], frame: Option<bool>) -> &[f64; N] {
        self.tick(q, Some(qd), a, b, frame)
    }

    #[allow(clippy::manual_is_multiple_of)]
    fn tick(&mut self, q: &[f64; N], qd: Option<&[f64; N]>, a: &[f64; N], b: &[f64; N], frame: Option<bool>) -> &[f64; N] {
        let dt = self.cfg.dt;
        let first = !self.started;
        let raw = self.cfg.identification == Identification::Raw;
        if raw { self.head = (self.head + 1) % RING; }
        let send = match frame { Some(f) => f, None => self.k % self.frame_every == 0 };
        let mut done = [false; 16];
        if !raw {
            for (i, dn) in done.iter_mut().enumerate().take(self.nsc) { self.cnt[i] += 1; *dn = self.cnt[i] == self.sc[i]; }
        }
        for j in 0..N {
            let meas = 0.5 * (a[j] + b[j]);
            let e = q[j] - meas;
            let d = (a[j] - b[j]) / sqrt(2.0);
            self.n2[j] += self.aw * (d * d - self.n2[j]);                       // per-read noise variance
            let r = (q[j] - self.thff[j]) + (meas - self.thm[j]);                // residual, exogenous channels
            let mut fitted: Option<(f64, f64, usize)> = None;
            if raw {
                self.hist[self.head][j] = r;
                if self.k >= self.maxlag as u64 {
                    for i in 0..self.nl {
                        let dd = r - self.hist[(self.head + RING - self.lag[i]) % RING][j];
                        self.dk[i][j] += self.aw * (dd * dd - self.dk[i][j]);
                    }
                }
                if self.k > self.warm {
                    let mut dv = [0.0; 8];
                    let mut fl = [0.0; 8];
                    for i in 0..self.nl { dv[i] = self.dk[i][j] - self.n2[j]; fl[i] = self.n2[j]; }   // noise floor removed
                    let (logl, nl) = (self.logl, self.nl);
                    fitted = Some(self.fit(&dv[..nl], &fl[..nl], &logl[..nl], 0));
                }
            } else {
                for (i, &dn) in done.iter().enumerate().take(self.nsc) {
                    self.acc[i][j] += r;
                    if dn {
                        let mean = self.acc[i][j] / self.sc[i] as f64;
                        if self.hasprev[i] {
                            let dd = mean - self.prev[i][j];
                            let aw = if self.nblk[i][j] > 0 { self.ablk[i] } else { 1.0 };
                            self.da[i][j] += aw * (dd * dd - self.da[i][j]);
                            self.nblk[i][j] += 1;
                        }
                        self.prev[i][j] = mean;
                        self.acc[i][j] = 0.0;
                    }
                }
                if self.k > self.warm {
                    let mut dv = [0.0; 16];
                    let mut fl = [0.0; 16];
                    for i in 0..self.nsc {
                        fl[i] = self.n2[j] * dt / (self.sc[i] as f64 * dt);                  // n2 dt / l, white read noise
                        dv[i] = if self.nblk[i][j] >= 2 { self.da[i][j] - fl[i] } else { -1.0 };
                    }
                    let (logs, nsc, finest) = (self.logs, self.nsc, self.finest);
                    let (hn, lam, nok) = self.fit(&dv[..nsc], &fl[..nsc], &logs[..nsc], finest);
                    fitted = Some((hn, lam / c_h(hn), nok));
                }
            }
            if let Some((hn, lam, nok)) = fitted {
                let s2 = 0.5 * self.n2[j] * dt;                                   // noise PSD of the pair mean
                let ts = kappa_h(hn) * pow((if s2 > 1e-30 { s2 } else { 1e-30 }) / (if lam > 1e-30 { lam } else { 1e-30 }),
                                           1.0 / (2.0 * hn + 1.0));
                let ts = clamp(ts, self.floor, self.cfg.tau_max);
                if nok >= 2 { self.tgt[j] = ts; self.h[j] = hn; }
                let rr = 1.0 / self.tau[j];
                if self.cfg.travel == Travel::E {                                  // theta = 1/tau relaxes (e-geodesic)
                    self.tau[j] = if first { self.tgt[j] } else { 1.0 / (rr + self.a_e * (1.0 / self.tgt[j] - rr)) };
                } else {
                    let drm = self.cfg.tau_rate * dt / (self.tau[j] * self.tau[j]);
                    self.tau[j] = 1.0 / (rr + clamp(1.0 / self.tgt[j] - rr, -drm, drm));
                }
            }
            let g = 1.0 - exp(-dt / self.tau[j]);
            self.u[j] = clamp(self.u[j] + g * e, -self.lim[j], self.lim[j]);
            let st = self.cfg.servo_tau;
            match qd {
                Some(v) if self.lead != 0.0 => {                                   // servo lead on the planned velocity
                    let ld = self.lead * v[j];
                    if send { self.cmd[j] = q[j] + ld + self.u[j]; }
                    self.thm[j] += dt * clamp((self.cmd[j] - self.thm[j]) / st, -self.w0[j], self.w0[j]);
                    self.thff[j] += dt * clamp((q[j] + ld - self.thff[j]) / st, -self.w0[j], self.w0[j]);
                }
                _ => {
                    if send { self.cmd[j] = q[j] + self.u[j]; }
                    self.thm[j] += dt * clamp((self.cmd[j] - self.thm[j]) / st, -self.w0[j], self.w0[j]);
                    self.thff[j] += dt * clamp((q[j] - self.thff[j]) / st, -self.w0[j], self.w0[j]);
                }
            }
        }
        if self.k > self.warm { self.started = true; }
        if !raw {
            for (i, &dn) in done.iter().enumerate().take(self.nsc) { if dn { self.cnt[i] = 0; self.hasprev[i] = true; } }
        }
        self.k += 1;
        &self.cmd
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn floor_is_loop_delay_or_th41a() {
        assert!((window_floor(0.001, 0.003, 1e-4) - 0.001 / 0.45).abs() < 1e-15);   // 2.22 ms
        assert!((loop_delay(0.001, 0.003, 1e-4) - 0.0021).abs() < 1e-15);
        assert!((window_floor(0.001, 0.008, 1e-4) - 0.0046).abs() < 1e-15);         // slow frame: the delay wins
    }

    #[test]
    fn allan_constant() {
        assert!((c_h(0.5) - 2.0 / 3.0).abs() < 1e-12);                                // Allan variance of a random walk
        assert!((c_h(1.0) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn limit_from_the_authority_fact() {
        let w0 = [5.0, 17.5];                                                          // rad/s
        let mut eng: Engine<2> = Engine::new(&[0.0; 2], &w0, Config::default());
        for _ in 0..2000 { eng.step(&[1.0; 2], &[0.0; 2], &[0.0; 2], None); }          // a blocked joint
        for j in 0..2 { assert!((eng.u[j].abs() - 0.025 * w0[j]).abs() < 1e-12); }     // servo_tau x w0
        let mut fixed: Engine<2> = Engine::new(&[0.0; 2], &w0, Config { limit: 0.1, ..Config::default() });
        for _ in 0..2000 { fixed.step(&[1.0; 2], &[0.0; 2], &[0.0; 2], None); }
        for j in 0..2 { assert!((fixed.u[j].abs() - 0.1).abs() < 1e-12); }
    }

    #[test]
    fn options_off_by_default() {
        let c = Config::default();
        assert_eq!(c.lead, 0.0);
        assert_eq!(c.travel, Travel::Rate);
        assert!(c.travel_t.is_none());
    }

    #[test]
    fn lead_adds_the_planned_velocity() {
        // a joint held exactly on its plan: the correction stays zero and the command is plan + lead x velocity
        let cfg = Config { lead: -1.0, ..Config::default() };
        let coef = 0.025 + loop_delay(0.001, 0.003, 1e-4);                             // servo_tau + read-to-action delay
        let mut eng: Engine<1> = Engine::new(&[0.0], &[10.0], cfg);
        let cmd = eng.step_lead(&[0.2], &[0.5], &[0.2], &[0.2], None)[0];
        assert!((cmd - (0.2 + coef * 0.5)).abs() < 1e-15);
        let mut fixed: Engine<1> = Engine::new(&[0.0], &[10.0], Config { lead: 0.01, ..Config::default() });
        assert!((fixed.step_lead(&[0.2], &[0.5], &[0.2], &[0.2], None)[0] - 0.205).abs() < 1e-15);
    }

    #[test]
    #[should_panic]
    fn lead_needs_step_lead() {
        let mut eng: Engine<1> = Engine::new(&[0.0], &[10.0], Config { lead: -1.0, ..Config::default() });
        eng.step(&[0.0], &[0.0], &[0.0], None);
    }

    #[test]
    fn e_travel_relaxes_in_one_over_tau() {
        // the first fit sets the window to the target; afterwards 1/tau relaxes toward 1/target with 64 ms
        let cfg = Config { travel: Travel::E, identification: Identification::Raw, ..Config::default() };
        let mut eng: Engine<1> = Engine::new(&[0.0], &[10.0], cfg);
        let a_e = 0.001 / 0.064;
        let mut checked = 0;
        for k in 0..3000u32 {
            let x = 0.3 * libm::sin(k as f64 * 0.004);
            let n = if k % 2 == 0 { 1e-4 } else { -1e-4 };
            let (was_started, old) = (eng.started, eng.tau[0]);
            eng.step(&[x], &[x + n], &[x - n], None);
            if eng.started && !was_started { assert_eq!(eng.tau[0], eng.tgt[0]); }
            if was_started {
                let want = 1.0 / old + a_e * (1.0 / eng.tgt[0] - 1.0 / old);
                assert!((1.0 / eng.tau[0] - want).abs() <= 1e-9 * want);
                checked += 1;
            }
        }
        assert!(checked > 1000);
    }

    #[test]
    fn kappa_matches_closed_forms() {
        assert!((kappa_h(0.5) - 1.0).abs() < 1e-12);                                  // canonical square-root law
        assert!((kappa_h(1.0) - pow(4.0, -1.0 / 3.0)).abs() < 1e-12);                 // H = 1: (2 Gamma(3))^(-1/3) = 0.63
    }
}
