// Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
// Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.

//! Syntonic engine: a joint correction loop whose time window comes from the Syntony Principle (Theorem U'),
//! with nothing tuned. Same arithmetic as the Python reference (`python/syntonic_engine`) and the C++ header
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

#[inline]
fn clamp(x: f64, lo: f64, hi: f64) -> f64 {
    if x < lo { lo } else if x > hi { hi } else { x }
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
    /// Correction limit (rad).
    pub limit: f64,
    /// Structure-function lags (s); the first `nlags` are used.
    pub lags: [f64; 8],
    /// Number of lags used.
    pub nlags: usize,
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
            limit: 5.0 * core::f64::consts::PI / 180.0,
            lags: [0.002, 0.004, 0.008, 0.016, 0.032, 0.064, 0.0, 0.0],
            nlags: 6,
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
        Engine {
            tau: [0.5; N], h: [0.5; N], u: [0.0; N], cmd: *q0, n2: [0.0; N],
            cfg, floor, aw: dt / cfg.ema, w0: *w0, thm: *q0, thff: *q0, tgt: [0.5; N],
            logl, lag, nl, maxlag: lag[nl - 1], frame_every,
            dk: [[0.0; N]; 8], hist: [[0.0; N]; RING], head: 0, k: 0,
            warm: round(cfg.warmup / dt) as u64,
        }
    }

    /// The window lower bound in use (s).
    pub fn floor(&self) -> f64 { self.floor }

    /// One loop tick. `q`: planned angles now (rad); `a`, `b`: two independent reads of every encoder (rad);
    /// `frame`: `Some(true)` when a servo frame is sent now, `Some(false)` when not, `None` to count frames
    /// internally (one every frame period, starting with the first tick). Returns the command to send (rad).
    #[allow(clippy::manual_is_multiple_of)]
    pub fn step(&mut self, q: &[f64; N], a: &[f64; N], b: &[f64; N], frame: Option<bool>) -> &[f64; N] {
        let dt = self.cfg.dt;
        self.head = (self.head + 1) % RING;
        let send = match frame { Some(f) => f, None => self.k % self.frame_every == 0 };
        for j in 0..N {
            let meas = 0.5 * (a[j] + b[j]);
            let e = q[j] - meas;
            let d = (a[j] - b[j]) / sqrt(2.0);
            self.n2[j] += self.aw * (d * d - self.n2[j]);                       // per-read noise variance
            let r = (q[j] - self.thff[j]) + (meas - self.thm[j]);                // residual, exogenous channels
            self.hist[self.head][j] = r;
            if self.k >= self.maxlag as u64 {
                for i in 0..self.nl {
                    let dd = r - self.hist[(self.head + RING - self.lag[i]) % RING][j];
                    self.dk[i][j] += self.aw * (dd * dd - self.dk[i][j]);
                }
            }
            if self.k > self.warm {
                let mut dv = [0.0; 8];
                let mut okg = [false; 8];
                let mut oku = [false; 8];
                let mut cg = 0;
                for i in 0..self.nl {
                    dv[i] = self.dk[i][j] - self.n2[j];                           // noise floor removed
                    okg[i] = dv[i] > self.cfg.snr_gate * self.n2[j];
                    oku[i] = dv[i] > 0.0;
                    if okg[i] { cg += 1; }
                }
                let (mut sw, mut sx, mut sy, mut sxx, mut sxy) = (0.0, 0.0, 0.0, 0.0, 0.0);
                let mut ok = [false; 8];
                let mut logd = [0.0; 8];
                let mut nok = 0;
                for i in 0..self.nl {
                    ok[i] = if cg >= 2 { okg[i] } else { oku[i] };
                    if ok[i] { nok += 1; }
                    logd[i] = if ok[i] { log(if dv[i] > 1e-30 { dv[i] } else { 1e-30 }) } else { 0.0 };
                    let w = if ok[i] { 1.0 } else { 0.0 };
                    let l = self.logl[i];
                    sw += w; sx += w * l; sy += w * logd[i]; sxx += w * l * l; sxy += w * l * logd[i];
                }
                let den = sw * sxx - sx * sx;
                let slope = if den > 1e-12 { (sw * sxy - sx * sy) / den } else { 1.0 };
                let hn = clamp(0.5 * slope, 0.05, 1.0);                           // smoothness H
                let s2 = 0.5 * self.n2[j] * dt;                                   // noise PSD of the pair mean
                let mut sl = 0.0;
                for i in 0..self.nl { if ok[i] { sl += logd[i] - 2.0 * hn * self.logl[i]; } }
                let lam = exp(sl / if sw > 1.0 { sw } else { 1.0 });
                let ts = kappa_h(hn) * pow((if s2 > 1e-30 { s2 } else { 1e-30 }) / (if lam > 1e-30 { lam } else { 1e-30 }),
                                           1.0 / (2.0 * hn + 1.0));
                let ts = clamp(ts, self.floor, self.cfg.tau_max);
                if nok >= 2 { self.tgt[j] = ts; self.h[j] = hn; }
                let rr = 1.0 / self.tau[j];
                let drm = self.cfg.tau_rate * dt / (self.tau[j] * self.tau[j]);
                self.tau[j] = 1.0 / (rr + clamp(1.0 / self.tgt[j] - rr, -drm, drm));
            }
            let g = 1.0 - exp(-dt / self.tau[j]);
            self.u[j] = clamp(self.u[j] + g * e, -self.cfg.limit, self.cfg.limit);
            if send { self.cmd[j] = q[j] + self.u[j]; }
            let st = self.cfg.servo_tau;
            self.thm[j] += dt * clamp((self.cmd[j] - self.thm[j]) / st, -self.w0[j], self.w0[j]);
            self.thff[j] += dt * clamp((q[j] - self.thff[j]) / st, -self.w0[j], self.w0[j]);
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
    fn kappa_matches_closed_forms() {
        assert!((kappa_h(0.5) - 1.0).abs() < 1e-12);                                  // canonical square-root law
        assert!((kappa_h(1.0) - pow(4.0, -1.0 / 3.0)).abs() < 1e-12);                 // H = 1: (2 Gamma(3))^(-1/3) = 0.63
    }
}
