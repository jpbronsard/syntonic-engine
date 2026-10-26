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
//! Two identifications of (H, lambda_H), selected by `Config::identification`: `Allan` (default since 0.2.0),
//! the structure function of the residual on the means of adjacent blocks of 2 ms to 4.096 s with the read-noise
//! floor n2 dt / l subtracted analytically, valid whether the read noise is small or large; and `Raw` (the 0.1.x
//! engine, bit for bit), the structure function of the raw residual at lags of 2 to 64 ms, valid where the motion
//! over 64 ms is far above the read noise, as on the MAESTRO 1 joint bench.
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
            identification: Identification::Allan,
            lags: [0.002, 0.004, 0.008, 0.016, 0.032, 0.064, 0.0, 0.0],
            nlags: 6,
            scales: [0.002, 0.004, 0.008, 0.016, 0.032, 0.064, 0.128, 0.256, 0.512, 1.024, 2.048, 4.096,
                     0.0, 0.0, 0.0, 0.0],
            nscales: 12,
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
        Engine {
            tau: [0.5; N], h: [0.5; N], u: [0.0; N], cmd: *q0, n2: [0.0; N],
            cfg, floor, aw: dt / cfg.ema, w0: *w0, thm: *q0, thff: *q0, tgt: [0.5; N],
            logl, lag, nl, maxlag: lag[nl - 1], frame_every,
            dk: [[0.0; N]; 8], hist: [[0.0; N]; RING], head: 0, k: 0,
            warm: round(cfg.warmup / dt) as u64,
            sc, nsc, logs, ablk, cnt: [0; 16], hasprev: [false; 16],
            acc: [[0.0; N]; 16], prev: [[0.0; N]; 16], da: [[0.0; N]; 16], nblk: [[0; N]; 16],
            finest: cfg.nlags,
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
    #[allow(clippy::manual_is_multiple_of)]
    pub fn step(&mut self, q: &[f64; N], a: &[f64; N], b: &[f64; N], frame: Option<bool>) -> &[f64; N] {
        let dt = self.cfg.dt;
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
    fn kappa_matches_closed_forms() {
        assert!((kappa_h(0.5) - 1.0).abs() < 1e-12);                                  // canonical square-root law
        assert!((kappa_h(1.0) - pow(4.0, -1.0 / 3.0)).abs() < 1e-12);                 // H = 1: (2 Gamma(3))^(-1/3) = 0.63
    }
}
