// Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
// Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.

//! The syntonic loop law: a joint correction whose PI gains are not tuned but computed, every 0.5 s, as the gains
//! that minimise the loop error predicted from statistics measured live. Same arithmetic as the Python reference
//! (`python/syntonic_engine/loop_law.py`, class `SyntonicLoopLaw`) and the C++ header (`cpp/SyntonicLoopLaw.h`).
//! Ported in 0.3.1; see the Python module docstring for the law and its confirmation.
//!
//! Facts (measured once): loop period `dt`, servo frame period, servo time constant `servo_tau` (data sheet; 25 ms),
//! no-load speed `w0` (data sheet), and the servo's internal delay `dead` (measured once at commissioning). Design
//! constants, declared: averaging time 1 s, warm-up 1 s, noise gate 2, block lags 2 to 256 ticks, at most 3 fitted
//! lags, gain update every 500 ticks, phase margin >= 30 deg and max |S| <= 2, the candidate grid of (Kp, Ki) (41 x 60)
//! and the frequency grid (1500 points, 0.1 to 1e4 rad/s).
//!
//! Each tick, per joint: two encoder reads give the joint angle and the read-pair noise R; the nominal servo models
//! give the error the loop must cancel; its block-mean structure function gives S(u) = A u^gamma and a corner Oc. Every
//! 500 ticks the gains are chosen among the robust candidates as those minimising the predicted loop error J. PI: the
//! integral accumulates Ki e dt every tick (frozen at the authority limit servo_tau x w0); on a frame tick the
//! proportional term uses the mean error since the last frame; command = q + L x qd + clip(Kp e_mean + I, +-limit).
//!
//! The gain recomputation evaluates J for the robust candidates (about 900 to 1250 of the 2460) over the 1500
//! frequencies, up to about 1.3 million complex evaluations for all joints together. The frequency responses are
//! recomputed each time from 84 KB of per-frequency constants rather than stored (30 MB), so the law takes about
//! 125 KB + 13 KB per joint. Since 0.3.2 the scan is cheaper in two ways that leave the chosen gains unchanged (both on
//! by default; their Python reference is `python/syntonic_engine/loop_scan.py`):
//! - `LoopLawConfig::prune`: every term of J is positive, so a partial sum never decreases. The candidates in use are
//!   evaluated first; every 50 frequencies a candidate is dropped once, for every joint, its partial J exceeds the
//!   smaller of their J and the best J so far, times 1 + `delta` (1e-9): it cannot be chosen. On the replays this
//!   evaluated 39 to 64 % of the points on average, but up to 96 % in one recomputation, so a per-tick budget cannot
//!   count on it.
//! - `LoopLawConfig::quad`: |1 + L|^2 = P + Ki (Q + Ki R), with P per Kp row and frequency, Q and R per frequency: four
//!   operations and a division per point instead of ten and a division (70 instead of 97 cycles by LLVM's Cortex-M7
//!   model). It changes J by rounding only (at most about 5e-13 relative for a robust candidate); the choice is
//!   certified when the next J exceeds the best by more than `delta`, otherwise the scan is done again from L, as in
//!   0.3.1 (a fallback, counted by `fallbacks()`; it never ran on the replays).
//!
//! `LoopLawConfig::spread` = 0 (default) recomputes and applies the gains in the tick where the Python law does: for two
//! joints on the realistic MAESTRO 1 bench about 0.8 ms on our PC in C++ (4.2 ms in 0.3.1), and by LLVM's Cortex-M7
//! model about 0.07 s on a Teensy 4.1 (0.16 s if nothing were pruned; 0.21 s in 0.3.1). `spread` = m divides the work,
//! counted as if nothing were pruned, over the m ticks that follow and applies the same gains m ticks later, which on
//! that bench raises the tip error by about 1.5 to 2 % at m = 250, 2.5 to 3.3 % at m = 365 and 4 to 5 % at m = 499
//! (`bench/port_2026-11-02/`, `bench/scan_2026-11-02/`). `time_candidate()` measures the time of one candidate on the
//! board and `spread_for()` gives the smallest spread whose share fits a declared time per tick (by the model, 365 for
//! 0.45 ms with two joints and 878 candidates): a scheduling choice, opt-in.
//!
//! On the replays of `tools/make_law_replay_log.py` the robust candidates and the chosen gains agree with Python
//! exactly at every tick, with the speed-ups off and on, and the continuous values to 2.4e-13 of their scale.

// index loops mirror the arithmetic of the Python class, element by element
#![allow(clippy::needless_range_loop, clippy::neg_cmp_op_on_partial_ord)]

use crate::loop_delay;
use libm::{atan2, cos, exp, fabs, floor, log, log10, pow, rint, sin, sqrt};

/// Points of the frequency grid, logspace(-1, 4, 1500) rad/s.
pub const LOOP_NW: usize = 1500;
/// Kp values of the candidate grid: 0 and geomspace(0.05, 12, 40).
pub const LOOP_NKP: usize = 41;
/// Ki values of the candidate grid: geomspace(5, 5000, 60).
pub const LOOP_NKI: usize = 60;
/// Candidates; candidate c has Kp number c / 60 and Ki number c % 60.
pub const LOOP_NC: usize = LOOP_NKP * LOOP_NKI;
const NLAG: usize = 8;
const FB: usize = 50;        // frequencies per work item
const NB: usize = LOOP_NW / FB;   // 30 blocks
const SPEC_FB: usize = 10;   // frequencies of one joint's spectrum per work item
// work units: one (candidate, frequency) point; the weights are for load balance only
const PHI_COST: u64 = 8;     // a spectrum point (a pow), in units
const ROW_UNITS: u64 = 30;   // FB frequencies of a row's P (42 against 70 cycles a point by the Cortex-M7 model)
const ITEM_MAX: u64 = if SPEC_FB as u64 * PHI_COST > FB as u64 + ROW_UNITS { SPEC_FB as u64 * PHI_COST } else { FB as u64 + ROW_UNITS };
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Phase { Spec, Warm, Scan, Done }

/// Block-mean constant K(gamma) of the second difference of adjacent block means (K(3) = 11/20 exact), linear between
/// the tabulated points as `numpy.interp`.
pub fn loop_k_gamma(g: f64) -> f64 {
    const X: [f64; 5] = [1.0, 1.5, 2.0, 2.5, 3.0];
    const Y: [f64; 5] = [0.99999999987, 0.74370490905, 0.59538447860, 0.52771124824, 0.55];
    if !(g > X[0]) { return Y[0]; }
    if !(g < X[4]) { return Y[4]; }
    let mut j = 0;
    while g >= X[j + 1] { j += 1; }
    if g == X[j] { return Y[j]; }
    let slope = (Y[j + 1] - Y[j]) / (X[j + 1] - X[j]);
    slope * (g - X[j]) + Y[j]
}

// numpy.geomspace for positive ends: 10^(linspace(log10 a, log10 b, n)) with the ends set exactly
fn geom(a: f64, b: f64, n: usize, i: usize) -> f64 {
    if i == 0 { return a; }
    if i == n - 1 { return b; }
    let (la, lb) = (log10(a), log10(b));
    let st = (lb - la) / (n - 1) as f64;
    pow(10.0, i as f64 * st + la)
}

#[inline]
fn clamp(x: f64, lo: f64, hi: f64) -> f64 {
    if x < lo { lo } else if x > hi { hi } else { x }
}

/// One block of frequencies `w0..w1` of one candidate, `j[i] += phiw[i][w] |S(w)|^2`, with |S|^2 from
/// L = (kp + ki / jW) G, as in 0.3.1.
#[inline]
#[allow(clippy::too_many_arguments)]
pub fn loop_block_l<const N: usize>(kp: f64, ki: f64, iw: &[f64; LOOP_NW], gr: &[f64; LOOP_NW], gi: &[f64; LOOP_NW],
                                    phiw: &[[f64; LOOP_NW]; N], w0: usize, w1: usize, j: &mut [f64; N]) {
    let mut jl = *j;
    for w in w0..w1 {
        let kiw = ki * iw[w];
        let lr = kp * gr[w] + kiw * gi[w];
        let li = kp * gi[w] - kiw * gr[w];
        let s2 = 1.0 / ((1.0 + lr) * (1.0 + lr) + li * li);   // |S|^2
        for i in 0..N { jl[i] += phiw[i][w] * s2; }
    }
    *j = jl;
}

/// The same by the quadratic form (0.3.2): |1 + L|^2 = P + ki (Q + ki R), P = 1 + kp (2 Re G + kp |G|^2) of the
/// candidate's Kp row.
#[inline]
#[allow(clippy::too_many_arguments)]
pub fn loop_block_quad<const N: usize>(ki: f64, p: &[f64; LOOP_NW], q: &[f64; LOOP_NW], r: &[f64; LOOP_NW],
                                       phiw: &[[f64; LOOP_NW]; N], w0: usize, w1: usize, j: &mut [f64; N]) {
    let mut jl = *j;
    for w in w0..w1 {
        let s2 = 1.0 / (p[w] + ki * (q[w] + ki * r[w]));        // |S|^2
        for i in 0..N { jl[i] += phiw[i][w] * s2; }
    }
    *j = jl;
}

/// Facts and declared design constants of the loop law. Only `dt`, `frame`, `servo_tau` and `dead` describe your
/// hardware; the rest are the law's declared constants and are not meant to be tuned.
#[derive(Clone, Copy, Debug)]
pub struct LoopLawConfig {
    /// Loop period (s).
    pub dt: f64,
    /// Servo command frame period (s).
    pub frame: f64,
    /// Servo time constant (s), from the data sheet; 0.025 when it gives none.
    pub servo_tau: f64,
    /// The servo's internal delay (s), measured once at commissioning.
    pub dead: f64,
    /// Lead coefficient L (s); `None` = `servo_tau + loop_delay(dt, frame, 1e-4)` (Python `lead=None`).
    pub lead: Option<f64>,
    /// Averaging time of the estimates (s).
    pub ema: f64,
    /// Data collected before the first fit (s).
    pub warmup: f64,
    /// A lag is used only this many white-noise floors above its floor.
    pub gate: f64,
    /// At most this many fitted lags.
    pub nfit: usize,
    /// Gain update period (ticks).
    pub every: u64,
    /// Phase margin (deg), at least.
    pub pm_min: f64,
    /// Max |S|, at most.
    pub ms_max: f64,
    /// 0: the gains are recomputed and applied in one tick, as in Python. m (1 <= m < every): the statistics are taken
    /// in the same tick, the work is spread evenly over the m ticks that follow, and the same gains take effect m ticks
    /// later (Rust and C++ only; its Python reference is `rust/tools/make_law_replay_log.py`).
    pub spread: u64,
    /// 0.3.2: skip the rest of a candidate's sum once it cannot be chosen (the same gains).
    pub prune: bool,
    /// 0.3.2: |1 + L|^2 as a quadratic in Ki (the same gains: certified, else from L again).
    pub quad: bool,
    /// 0.3.2: relative margin of the pruning bound and of `quad`'s certificate, declared; at least 1e-10 (far above the
    /// rounding bound 5e-13). A larger value costs time and never changes the gains.
    pub delta: f64,
}

impl Default for LoopLawConfig {
    fn default() -> Self {
        LoopLawConfig {
            dt: 0.001, frame: 0.003, servo_tau: 0.025, dead: 0.0, lead: None, ema: 1.0, warmup: 1.0, gate: 2.0,
            nfit: 3, every: 500, pm_min: 30.0, ms_max: 2.0, spread: 0, prune: true, quad: true, delta: 1e-9,
        }
    }
}

/// The loop law for `N` joints. `PIPE` is the length of the nominal servo models' delay line in ticks and must be at
/// least `dead / dt + 2`. It takes about 125 KB + 13 KB per joint: on a small board keep it in a `static`.
pub struct SyntonicLoopLaw<const N: usize, const PIPE: usize = 16> {
    /// Proportional gain in use per joint (Python `Kp`).
    pub kp: [f64; N],
    /// Integral gain in use per joint (Python `Ki`, 1/s).
    pub ki: [f64; N],
    /// Integral term per joint (rad; Python `I`).
    pub integral: [f64; N],
    /// Read-pair noise variance per joint (rad^2; Python `R`): EMA_1s((a - b)^2 / 2) / 2.
    pub noise: [f64; N],
    /// Fitted structure function S(u) = amp u^gam per joint (Python `A`, `gam`).
    pub amp: [f64; N],
    /// Exponent of the fitted structure function.
    pub gam: [f64; N],
    /// Corner of the disturbance spectrum per joint (rad/s; Python `Oc`).
    pub oc: [f64; N],
    /// Nominal servo model driven by the commands sent (rad).
    pub th_m: [f64; N],
    /// Nominal servo model driven by the plan feed-forward alone (rad).
    pub th_ff: [f64; N],
    /// Command per joint (rad), held between frames.
    pub cmd: [f64; N],
    /// Lead coefficient (s; Python `L`).
    pub lead: f64,
    /// True once a fit has used two lags or more on some joint.
    pub started: bool,
    cfg: LoopLawConfig,
    d: f64,
    aw: f64,
    fd: f64,
    nd: usize,
    warm: u64,
    k: u64,
    job: bool,
    job_k: u64,
    work_total: u64,
    work_done: u64,
    item: usize,
    phase: Phase,
    head: usize,
    ecnt: u64,
    nok: usize,
    w0: [f64; N],
    lim: [f64; N],
    n2: [f64; N],
    esum: [f64; N],
    cff: [f64; N],
    pc: [[f64; N]; PIPE],
    pf: [[f64; N]; PIPE],
    // identification
    lag: [u64; NLAG],
    cnt: [u64; NLAG],
    nblk: [u64; NLAG],
    nb: [u64; NLAG],
    u: [f64; NLAG],
    lu: [f64; NLAG],
    ab: [f64; NLAG],
    acc: [[f64; N]; NLAG],
    m1: [[f64; N]; NLAG],
    m2: [[f64; N]; NLAG],
    s2: [[f64; N]; NLAG],
    // grid
    w2: [f64; LOOP_NW],
    iw: [f64; LOOP_NW],
    gr: [f64; LOOP_NW],
    gi: [f64; LOOP_NW],
    dwf: [f64; LOOP_NW],
    kpv: [f64; LOOP_NKP],
    kiv: [f64; LOOP_NKI],
    mtr: [f64; LOOP_NC],
    ok: [bool; LOOP_NC],
    okc: [u16; LOOP_NC],
    // the quadratic form: per frequency, and P of the Kp row `prow`
    qf: [f64; LOOP_NW],
    rf: [f64; LOOP_NW],
    pf_row: [f64; LOOP_NW],
    nrows: usize,
    // recomputation
    s_amp: [f64; N],
    s_gam: [f64; N],
    s_oc: [f64; N],
    s_r: [f64; N],
    c_j: [f64; N],
    ex: [f64; N],
    o2: [f64; N],
    best: [i32; N],
    phiw: [[f64; LOOP_NW]; N],
    prow: usize,
    pvalid: usize,
    p: usize,
    blk: usize,
    wk: usize,
    nwarm: usize,
    nan: bool,
    quad_now: bool,
    prune_now: bool,
    warmc: [usize; N],
    warmp: [usize; N],
    c1: [usize; N],
    out: [bool; N],
    warm_j: [[f64; N]; N],
    bw: [f64; N],
    b1: [f64; N],
    b2: [f64; N],
    thr: [f64; N],
    jp: [f64; N],
    evaluated: u64,
    fallbacks: u64,
    late: u64,
    margin: f64,
    sink: f64,
}

impl<const N: usize, const PIPE: usize> SyntonicLoopLaw<N, PIPE> {
    /// `q0`: start angles (rad); `w0`: no-load speeds (rad/s), from the data sheet. Builds the frequency responses of
    /// the candidates. Panics if a fact is out of range (as Python raises), in particular if `PIPE < dead / dt + 2` or
    /// `spread >= every`.
    pub fn new(q0: &[f64; N], w0: &[f64; N], cfg: LoopLawConfig) -> Self {
        assert!(cfg.dt > 0.0 && cfg.servo_tau > 0.0 && cfg.dead >= 0.0, "dt, servo_tau must be positive, dead not negative");
        assert!(cfg.every >= 1 && cfg.nfit >= 1, "every and nfit must be at least 1");
        assert!(cfg.spread == 0 || cfg.spread < cfg.every, "spread must be smaller than every");
        assert!(cfg.delta >= 1e-10, "delta must be at least 1e-10");
        let dt = cfg.dt;
        let ndf = cfg.dead / dt;
        let nd = floor(ndf + 1e-9) as usize;
        let fd = ndf - nd as f64;
        assert!(nd + 2 <= PIPE, "PIPE must be at least dead / dt + 2");
        let mut lag = [0u64; NLAG];
        let (mut u, mut lu, mut ab) = ([0.0; NLAG], [0.0; NLAG], [0.0; NLAG]);
        for i in 0..NLAG {
            lag[i] = 2 << i;
            u[i] = lag[i] as f64 * dt;
            lu[i] = log(u[i]);
            let e8 = 8.0 * dt * lag[i] as f64;
            let b = dt * lag[i] as f64 / if cfg.ema > e8 { cfg.ema } else { e8 };
            ab[i] = if b < 1.0 { b } else { 1.0 };
        }
        let mut lim = [0.0; N];
        for j in 0..N { lim[j] = cfg.servo_tau * w0[j]; }
        let mut law = SyntonicLoopLaw {
            kp: [0.0; N], ki: [50.0; N], integral: [0.0; N], noise: [1e-14; N], amp: [0.0; N], gam: [3.0; N],
            oc: [1.0; N], th_m: *q0, th_ff: *q0, cmd: *q0,
            lead: match cfg.lead { Some(l) => l, None => cfg.servo_tau + loop_delay(dt, cfg.frame, 1e-4) },
            started: false, cfg,
            d: loop_delay(dt, cfg.frame, 1e-4) + cfg.dead,
            aw: dt / cfg.ema, fd, nd, warm: rint(cfg.warmup / dt) as u64, k: 0,
            job: false, job_k: 0, work_total: 0, work_done: 0, item: 0, phase: Phase::Done, head: 0, ecnt: 0, nok: 0,
            w0: *w0, lim, n2: [0.0; N], esum: [0.0; N], cff: *q0, pc: [*q0; PIPE], pf: [*q0; PIPE],
            lag, cnt: [0; NLAG], nblk: [0; NLAG], nb: [0; NLAG], u, lu, ab,
            acc: [[0.0; N]; NLAG], m1: [[0.0; N]; NLAG], m2: [[0.0; N]; NLAG], s2: [[0.0; N]; NLAG],
            w2: [0.0; LOOP_NW], iw: [0.0; LOOP_NW], gr: [0.0; LOOP_NW], gi: [0.0; LOOP_NW], dwf: [0.0; LOOP_NW],
            kpv: [0.0; LOOP_NKP], kiv: [0.0; LOOP_NKI], mtr: [0.0; LOOP_NC], ok: [false; LOOP_NC], okc: [0; LOOP_NC],
            qf: [0.0; LOOP_NW], rf: [0.0; LOOP_NW], pf_row: [0.0; LOOP_NW], nrows: 0,
            s_amp: [0.0; N], s_gam: [0.0; N], s_oc: [0.0; N], s_r: [0.0; N], c_j: [0.0; N], ex: [0.0; N],
            o2: [0.0; N], best: [-1; N], phiw: [[0.0; LOOP_NW]; N],
            prow: usize::MAX, pvalid: 0, p: 0, blk: 0, wk: 0, nwarm: 0, nan: false, quad_now: false, prune_now: false,
            warmc: [0; N], warmp: [0; N], c1: [0; N], out: [false; N], warm_j: [[0.0; N]; N], bw: [0.0; N],
            b1: [0.0; N], b2: [0.0; N], thr: [0.0; N], jp: [0.0; N], evaluated: 0, fallbacks: 0, late: 0,
            margin: f64::INFINITY, sink: 0.0,
        };
        law.grid();
        law
    }

    // frequency responses of the PI loop G = exp(-s D) / (servo_tau s + 1) over the candidate grid; robustness tests
    fn grid(&mut self) {
        let pi = core::f64::consts::PI;
        let step = (4.0 - (-1.0)) / (LOOP_NW - 1) as f64;
        let twopi = 2.0 * pi;
        let ts = self.cfg.servo_tau;
        let w_at = |w: usize| -> f64 { pow(10.0, if w == LOOP_NW - 1 { 4.0 } else { w as f64 * step + (-1.0) }) };
        for w in 0..LOOP_NW {
            let ww = w_at(w);
            let dw = if w == 0 { w_at(1) - w_at(0) } else if w == LOOP_NW - 1 { ww - w_at(w - 1) }
                     else { (w_at(w + 1) - w_at(w - 1)) / 2.0 };
            self.dwf[w] = dw / twopi * 2.0;
            self.w2[w] = ww * ww;
            self.iw[w] = 1.0 / ww;
            let (cr, ci, dr, di) = (cos(ww * self.d), -sin(ww * self.d), 1.0, ww * ts);   // exp(-s D) / (1 + s ts)
            let m = dr * dr + di * di;
            self.gr[w] = (cr * dr + ci * di) / m;
            self.gi[w] = (ci * dr - cr * di) / m;
            self.qf[w] = 2.0 * (self.gi[w] * self.iw[w]);                   // 0.3.2: the quadratic form
            self.rf[w] = (self.iw[w] * self.iw[w]) * (self.gr[w] * self.gr[w] + self.gi[w] * self.gi[w]);
        }
        self.kpv[0] = 0.0;
        for i in 1..LOOP_NKP { self.kpv[i] = geom(0.05, 12.0, LOOP_NKP - 1, i - 1); }
        for i in 0..LOOP_NKI { self.kiv[i] = geom(5.0, 5000.0, LOOP_NKI, i); }
        let r2d = 180.0 / pi;
        let ms2 = self.cfg.ms_max * self.cfg.ms_max;
        self.nok = 0;
        for c in 0..LOOP_NC {
            let (kp, ki) = (self.kpv[c / LOOP_NKI], self.kiv[c % LOOP_NKI]);
            let mut cross = false;
            let (mut s2max, mut mtr, mut l0, mut lrc, mut lic) = (0.0f64, 0.0f64, 0.0f64, 0.0f64, 0.0f64);
            for w in 0..LOOP_NW {
                let kiw = ki * self.iw[w];
                let lr = kp * self.gr[w] + kiw * self.gi[w];
                let li = kp * self.gi[w] - kiw * self.gr[w];
                let l2 = lr * lr + li * li;
                let s2 = 1.0 / ((1.0 + lr) * (1.0 + lr) + li * li);
                if w == 0 { l0 = l2; lrc = lr; lic = li; }
                if !cross && l2 < 1.0 { cross = true; lrc = lr; lic = li; }
                if s2 > s2max { s2max = s2; }
                mtr += l2 * s2 * self.dwf[w];                          // |T|^2 = |L|^2 |S|^2
            }
            let mut pm = 180.0 + atan2(lic, lrc) * r2d;                // at the first frequency where |L| < 1
            if pm > 180.0 { pm -= 360.0; }
            self.ok[c] = pm >= self.cfg.pm_min && s2max <= ms2 && l0 > 1.0;
            if self.ok[c] { self.okc[self.nok] = c as u16; self.mtr[self.nok] = mtr; self.nok += 1; }
        }
        self.nrows = (0..LOOP_NKP).filter(|&i| (0..LOOP_NKI).any(|k| self.ok[i * LOOP_NKI + k])).count();
    }

    /// One tick. `q`, `qd`: planned angles (rad) and velocities (rad/s) now; `a`, `b`: two reads of every encoder (rad);
    /// `send`: true on a servo frame tick. Returns the command (rad), which changes on frame ticks only.
    #[allow(clippy::manual_is_multiple_of)]
    pub fn step(&mut self, q: &[f64; N], qd: &[f64; N], a: &[f64; N], b: &[f64; N], send: bool) -> &[f64; N] {
        let dt = self.cfg.dt;
        let mut y = [0.0; N];
        let mut z = [0.0; N];
        for j in 0..N {
            y[j] = 0.5 * (a[j] + b[j]);
            let dn = (a[j] - b[j]) / sqrt(2.0);
            self.n2[j] += self.aw * (dn * dn - self.n2[j]);
            let h = 0.5 * self.n2[j];
            self.noise[j] = if h > 1e-14 { h } else { 1e-14 };
            z[j] = (q[j] - self.th_ff[j]) + (self.th_m[j] - y[j]);
        }
        self.identify(&z);
        if self.started && self.k % self.cfg.every == 0 {
            for j in 0..N {
                self.s_amp[j] = self.amp[j]; self.s_gam[j] = self.gam[j]; self.s_oc[j] = self.oc[j];
                self.s_r[j] = self.noise[j] * dt;
            }
            self.start_job();
            if self.cfg.spread == 0 {
                while self.work_one() {}
                self.apply();
            } else {
                self.job = true;
                self.job_k = self.k;
            }
        } else if self.job {
            let s = self.k - self.job_k;                                 // slice s of spread, s = 1 .. spread
            let target = (self.work_total * s).div_ceil(self.cfg.spread);
            while self.work_done < target && self.work_one() {}
            if s >= self.cfg.spread {                                     // the last slice: whatever is left
                if self.phase != Phase::Done {
                    self.late += 1;
                    while self.work_one() {}
                }
                self.apply();
                self.job = false;
            }
        }
        let mut e = [0.0; N];
        for j in 0..N { e[j] = q[j] - y[j]; self.esum[j] += e[j]; }
        self.ecnt += 1;
        for j in 0..N {
            self.integral[j] += if fabs(self.kp[j] * e[j] + self.integral[j]) >= self.lim[j] { 0.0 } else { self.ki[j] * e[j] * dt };
        }
        if send {
            for j in 0..N {
                let eb = self.esum[j] / self.ecnt as f64;
                self.esum[j] = 0.0;
                self.cmd[j] = q[j] + self.lead * qd[j] + clamp(self.kp[j] * eb + self.integral[j], -self.lim[j], self.lim[j]);
                self.cff[j] = q[j] + self.lead * qd[j];
            }
            self.ecnt = 0;
        }
        self.head = (self.head + 1) % PIPE;
        let h0 = (self.head + PIPE - self.nd) % PIPE;
        let h1 = (self.head + PIPE - self.nd - 1) % PIPE;
        let (fd, ts) = (self.fd, self.cfg.servo_tau);
        for j in 0..N {
            self.pc[self.head][j] = self.cmd[j];
            self.pf[self.head][j] = self.cff[j];
            let cm = if fd == 0.0 { self.pc[h0][j] } else { (1.0 - fd) * self.pc[h0][j] + fd * self.pc[h1][j] };
            let cf = if fd == 0.0 { self.pf[h0][j] } else { (1.0 - fd) * self.pf[h0][j] + fd * self.pf[h1][j] };
            self.th_m[j] += dt * clamp((cm - self.th_m[j]) / ts, -self.w0[j], self.w0[j]);
            self.th_ff[j] += dt * clamp((cf - self.th_ff[j]) / ts, -self.w0[j], self.w0[j]);
        }
        self.k += 1;
        &self.cmd
    }

    /// Ticks done.
    pub fn tick(&self) -> u64 { self.k }
    /// Number of candidates passing the robustness tests.
    pub fn robust_count(&self) -> usize { self.nok }
    /// Candidate `c` passes the robustness tests.
    pub fn robust(&self, c: usize) -> bool { self.ok[c] }
    /// Kp of candidate `c`.
    pub fn grid_kp(&self, c: usize) -> f64 { self.kpv[c / LOOP_NKI] }
    /// Ki of candidate `c`.
    pub fn grid_ki(&self, c: usize) -> f64 { self.kiv[c % LOOP_NKI] }
    /// Candidate in use on joint `j`, `None` before the first choice.
    pub fn gain_index(&self, j: usize) -> Option<usize> { if self.best[j] < 0 { None } else { Some(self.best[j] as usize) } }
    /// Authority limit `servo_tau * w0` of joint `j` (rad).
    pub fn limit(&self, j: usize) -> f64 { self.lim[j] }
    /// Loop delay of the model, `loop_delay(dt, frame, 1e-4) + dead` (s).
    pub fn delay(&self) -> f64 { self.d }
    /// With `spread` > 0: a recomputation is under way.
    pub fn pending(&self) -> bool { self.job }
    /// The spread in use.
    pub fn spread(&self) -> u64 { self.cfg.spread }
    /// 0.3.2: (candidate, frequency) points the last recomputation evaluated.
    pub fn points_evaluated(&self) -> u64 { self.evaluated }
    /// The same without the speed-ups (as in 0.3.1).
    pub fn points_total(&self) -> u64 { (self.nok * LOOP_NW) as u64 }
    /// Recomputations done again from L (not certified, or a NaN).
    pub fn fallbacks(&self) -> u64 { self.fallbacks }
    /// With `spread` > 0: recomputations not finished before their last slice.
    pub fn late(&self) -> u64 { self.late }
    /// The smallest relative gap seen between the best J and the next among the candidates evaluated to the end
    /// (with `prune`, the next is often dropped early: then this is larger than the true gap).
    pub fn tie_margin(&self) -> f64 { self.margin }

    /// 0.3.2, scheduling. The time one robust candidate's 1500 points take on this board, in seconds, by the caller's
    /// clock: `now()` returns a 32-bit count that increases `hz` times a second and may wrap (on a Cortex-M7, the
    /// cycle counter DWT CYCCNT). With the options of `new()` (quad or not). Call it when no recomputation is under way;
    /// the shortest of `reps` runs, or `None`.
    pub fn time_candidate<F: FnMut() -> u32>(&mut self, mut now: F, hz: f64, reps: usize) -> Option<f64> {
        if self.job || self.nok == 0 || !(hz > 0.0) { return None; }
        for j in 0..N {
            for w in 0..LOOP_NW { self.phiw[j][w] = 1e-6 * pow(self.w2[w] + 1.0, -2.0) * self.dwf[w]; }   // gamma 3
        }
        let c = self.okc[self.nok / 2] as usize;
        self.quad_now = self.cfg.quad;
        self.prow = c / LOOP_NKI;
        self.pvalid = 0;
        if self.quad_now { for b in 0..NB { self.row(self.prow, b); self.pvalid += 1; } }   // P ready: the candidate alone
        let mut best: Option<f64> = None;
        for _ in 0..reps.max(1) {
            self.jp = [0.0; N];
            let t0 = now();
            for b in 0..NB { self.block(c, b); }
            let t1 = now();
            let sec = t1.wrapping_sub(t0) as f64 / hz;
            best = Some(match best { Some(x) if x <= sec => x, _ => sec });
        }
        self.sink = self.jp[0];
        self.prow = usize::MAX;
        self.work_done = 0;
        self.evaluated = 0;
        best
    }

    /// The smallest spread whose share of the recomputation takes at most `budget` seconds per tick (on top of the
    /// ordinary tick), given the time of one candidate (`time_candidate()`). The work is counted as if nothing were
    /// pruned, so the budget holds in every recomputation. `Some(0)` if the whole recomputation fits in one tick,
    /// `None` if no spread below `every` fits. A scheduling choice from a measured cost and a declared budget.
    pub fn spread_for(&self, seconds_per_candidate: f64, budget: f64) -> Option<u64> {
        if !(seconds_per_candidate > 0.0) || !(budget > 0.0) { return None; }
        let u = self.work_bound(N);
        let unit = seconds_per_candidate / LOOP_NW as f64;
        if u as f64 * unit <= budget { return Some(0); }
        (1..self.cfg.every).find(|&m| (u.div_ceil(m) + ITEM_MAX) as f64 * unit <= budget)
    }

    /// Changes the spread between recomputations; false (and nothing changed) during one or if `m` is out of range.
    pub fn set_spread(&mut self, m: u64) -> bool {
        if self.job || (m > 0 && m >= self.cfg.every) { return false; }
        self.cfg.spread = m;
        true
    }

    // block-mean (Hadamard) identification of the error the loop must cancel: second differences of adjacent block
    // means at lags 2 .. 256 ticks, the white floor 6 R dt / u removed, S(u) = A u^gamma on the finest growing lags
    fn identify(&mut self, z: &[f64; N]) {
        for i in 0..NLAG {
            self.cnt[i] += 1;
            for j in 0..N { self.acc[i][j] += z[j]; }
            if self.cnt[i] == self.lag[i] {
                let upd = self.nblk[i] >= 2;
                let a = if self.nb[i] > 0 { self.ab[i] } else { 1.0 };
                for j in 0..N {
                    let mean = self.acc[i][j] / self.lag[i] as f64;
                    if upd {
                        let dd = mean - 2.0 * self.m1[i][j] + self.m2[i][j];
                        self.s2[i][j] += a * (dd * dd - self.s2[i][j]);
                    }
                    self.m2[i][j] = self.m1[i][j];
                    self.m1[i][j] = mean;
                    self.acc[i][j] = 0.0;
                }
                if upd { self.nb[i] += 1; }
                self.nblk[i] += 1;
                self.cnt[i] = 0;
            }
        }
        if self.k <= self.warm { return; }
        let dt = self.cfg.dt;
        let mut any = false;
        for j in 0..N {
            let (mut valid, mut ok) = ([false; NLAG], [false; NLAG]);
            let mut ls = [0.0; NLAG];
            let mut used = 0;
            for i in 0..NLAG {
                let fl = 6.0 * self.noise[j] * dt / self.u[i];
                let sig = self.s2[i][j] - fl;
                valid[i] = sig > self.cfg.gate * fl && self.nb[i] >= 2;
                ls[i] = log(if sig > 1e-30 { sig } else { 1e-30 });
                ok[i] = valid[i] && (i == 0 || ls[i] > ls[i - 1]);
                if ok[i] { used += 1; if used > self.cfg.nfit { ok[i] = false; } }
            }
            let (mut tw, mut tx, mut ty, mut txx, mut txy) = ([0.0; NLAG], [0.0; NLAG], [0.0; NLAG], [0.0; NLAG], [0.0; NLAG]);
            for i in 0..NLAG {
                let w = if ok[i] { 1.0 } else { 0.0 };
                tw[i] = w; tx[i] = w * self.lu[i]; ty[i] = w * ls[i];
                txx[i] = w * self.lu[i] * self.lu[i]; txy[i] = w * self.lu[i] * ls[i];
            }
            let (sw, sx, sy, sxx, sxy) = (sum8::<N>(&tw), sum8::<N>(&tx), sum8::<N>(&ty), sum8::<N>(&txx), sum8::<N>(&txy));
            let den = sw * sxx - sx * sx;
            let g = if den > 1e-12 { (sw * sxy - sx * sy) / den } else { 3.0 };
            let g = clamp(g, 1.0, 3.0);
            let la = (sy - g * sx) / if sw > 1.0 { sw } else { 1.0 };
            if sw >= 2.0 { self.gam[j] = g; self.amp[j] = exp(la); any = true; }
            let mut ug = self.u[NLAG - 1];
            if let Some(mut i0) = valid.iter().position(|&v| v) {
                while i0 + 1 < NLAG && valid[i0 + 1] && ls[i0 + 1] > ls[i0] { i0 += 1; }
                ug = self.u[i0];
            }
            self.oc[j] = 1.0 / ug;
        }
        if any { self.started = true; }
    }

    // The recomputation as a list of work items: the spectrum (Phi of every joint, SPEC_FB frequencies of one joint per
    // item); with prune, the candidates in use (each once); then the robust candidates in index order. A candidate's J
    // is summed over the 1500 frequencies in order, FB per item; with quad, the item first fills those FB frequencies
    // of the Kp row's P if no candidate of the row has reached them yet. The order of the arithmetic does not depend on
    // how the items are divided between ticks.
    fn work_bound(&self, nwarm: usize) -> u64 {
        let mut u = LOOP_NW as u64 * PHI_COST * N as u64 + (self.nok * LOOP_NW) as u64;
        if self.cfg.prune { u += (nwarm * LOOP_NW) as u64; }
        if self.cfg.quad { u += (self.nrows + if self.cfg.prune { nwarm } else { 0 }) as u64 * NB as u64 * ROW_UNITS; }
        u
    }

    fn start_job(&mut self) {
        for j in 0..N {
            self.c_j[j] = self.s_amp[j] / loop_k_gamma(self.s_gam[j]);
            self.ex[j] = -(1.0 + self.s_gam[j]) / 2.0;
            self.o2[j] = self.s_oc[j] * self.s_oc[j];
        }
        self.nwarm = 0;
        if self.cfg.prune {
            for j in 0..N {                                              // the candidates in use, each once, if robust
                if self.best[j] < 0 { continue; }
                let c = self.best[j] as usize;
                if !self.ok[c] || self.warmc[..self.nwarm].contains(&c) { continue; }
                self.warmc[self.nwarm] = c;
                self.warmp[self.nwarm] = self.okc[..self.nok].binary_search(&(c as u16)).unwrap();
                self.nwarm += 1;
            }
        }
        self.work_total = self.work_bound(self.nwarm);
        self.work_done = 0;
        self.evaluated = 0;
        self.phase = Phase::Spec;
        self.item = 0;
        self.nan = false;
    }

    // a pass over the candidates, by the quadratic form or from L; a NaN in the spectrum or the noise turns both
    // speed-ups off (the partial sums are then not monotone), as 0.3.1
    fn begin_pass(&mut self, quad: bool) {
        self.quad_now = quad && !self.nan;
        self.prune_now = self.cfg.prune && !self.nan;
        self.b1 = [f64::INFINITY; N];
        self.b2 = [f64::INFINITY; N];
        self.c1 = [0; N];
        self.bw = [f64::INFINITY; N];
        self.p = 0;
        self.blk = 0;
        self.wk = 0;
        self.prow = usize::MAX;
        if self.prune_now && self.nwarm > 0 { self.phase = Phase::Warm; } else { self.phase = Phase::Scan; self.advance(); }
    }

    // to the scan's next candidate: those in use are not evaluated twice; after the last, the end of the pass
    fn advance(&mut self) {
        while self.p < self.nok {
            let k = if self.prune_now { self.warmp[..self.nwarm].iter().position(|&x| x == self.p) } else { None };
            match k {
                Some(k) => { let jw = self.warm_j[k]; self.complete(self.okc[self.p] as usize, &jw); self.p += 1; }
                None => return,
            }
        }
        self.end_pass();
    }

    // the choice: certified if the next J exceeds the best by more than delta (needed with quad only)
    fn end_pass(&mut self) {
        let mut cert = true;
        for j in 0..N {
            if !(self.b2[j] > self.b1[j] * (1.0 + self.cfg.delta)) { cert = false; }
            let m = (self.b2[j] - self.b1[j]) / fabs(self.b1[j]);
            if m < self.margin { self.margin = m; }
        }
        if self.quad_now && !cert {
            self.fallbacks += 1;
            self.begin_pass(false);
            return;
        }
        self.phase = Phase::Done;
    }

    fn complete(&mut self, c: usize, jv: &[f64; N]) {                    // the first smallest J, as numpy.argmin
        for j in 0..N {
            let v = jv[j];
            if v < self.b1[j] || (v.is_nan() && !self.b1[j].is_nan()) {
                self.b2[j] = self.b1[j];
                self.b1[j] = v;
                self.c1[j] = c;
            } else if v < self.b2[j] {
                self.b2[j] = v;
            }
        }
    }

    fn row(&mut self, i: usize, blk: usize) {                            // P of Kp row i, block blk
        let kp = self.kpv[i];
        for w in blk * FB..(blk + 1) * FB {
            let (gr, gi) = (self.gr[w], self.gi[w]);
            self.pf_row[w] = 1.0 + kp * (2.0 * gr + kp * (gr * gr + gi * gi));
        }
        self.work_done += ROW_UNITS;
    }

    fn block(&mut self, c: usize, blk: usize) {
        let (w0, w1) = (blk * FB, (blk + 1) * FB);
        if self.quad_now {
            if c / LOOP_NKI != self.prow { self.prow = c / LOOP_NKI; self.pvalid = 0; }
            if blk == self.pvalid { self.row(self.prow, blk); self.pvalid += 1; }   // the blocks are reached in order
        }
        if self.quad_now {
            loop_block_quad::<N>(self.kiv[c % LOOP_NKI], &self.pf_row, &self.qf, &self.rf, &self.phiw, w0, w1, &mut self.jp);
        } else {
            loop_block_l::<N>(self.kpv[c / LOOP_NKI], self.kiv[c % LOOP_NKI], &self.iw, &self.gr, &self.gi, &self.phiw,
                              w0, w1, &mut self.jp);
        }
        self.work_done += FB as u64;
        self.evaluated += FB as u64;
    }

    fn work_one(&mut self) -> bool {
        match self.phase {
            Phase::Spec => {
                let (blk, j) = (self.item / N, self.item % N);
                for w in blk * SPEC_FB..(blk + 1) * SPEC_FB {
                    let v = self.c_j[j] * pow(self.w2[w] + self.o2[j], self.ex[j]) * self.dwf[w];
                    self.phiw[j][w] = v;
                    if v.is_nan() { self.nan = true; }
                }
                self.work_done += SPEC_FB as u64 * PHI_COST;
                self.item += 1;
                if self.item == LOOP_NW / SPEC_FB * N {
                    if self.s_r.iter().any(|x| x.is_nan()) { self.nan = true; }
                    if self.nan { self.fallbacks += 1; }
                    self.begin_pass(self.cfg.quad);
                }
                true
            }
            Phase::Warm => {
                let c = self.warmc[self.wk];
                if self.blk == 0 { self.jp = [0.0; N]; }
                self.block(c, self.blk);
                self.blk += 1;
                if self.blk == NB {
                    for j in 0..N {
                        let v = self.jp[j] + self.s_r[j] * self.mtr[self.warmp[self.wk]];
                        self.warm_j[self.wk][j] = v;
                        if v < self.bw[j] { self.bw[j] = v; }
                    }
                    self.blk = 0;
                    self.wk += 1;
                    if self.wk == self.nwarm { self.phase = Phase::Scan; self.advance(); }
                }
                true
            }
            Phase::Scan => {
                let c = self.okc[self.p] as usize;
                if self.blk == 0 {
                    for j in 0..N {
                        self.jp[j] = 0.0;
                        self.out[j] = false;
                        self.thr[j] = (if self.bw[j] < self.b1[j] { self.bw[j] } else { self.b1[j] }) * (1.0 + self.cfg.delta);
                    }
                }
                self.block(c, self.blk);
                self.blk += 1;
                if self.prune_now {                                      // dropped once no joint can choose it
                    let mut all = true;
                    for j in 0..N {
                        if !self.out[j] {
                            if self.jp[j] + self.s_r[j] * self.mtr[self.p] > self.thr[j] { self.out[j] = true; } else { all = false; }
                        }
                    }
                    if all { self.p += 1; self.blk = 0; self.advance(); return true; }
                }
                if self.blk == NB {
                    let mut jv = [0.0; N];
                    for j in 0..N { jv[j] = self.jp[j] + self.s_r[j] * self.mtr[self.p]; }
                    self.complete(c, &jv);
                    self.p += 1;
                    self.blk = 0;
                    self.advance();
                }
                true
            }
            Phase::Done => false,
        }
    }

    fn apply(&mut self) {
        for j in 0..N {
            let c = self.c1[j];
            self.best[j] = c as i32;
            self.kp[j] = self.kpv[c / LOOP_NKI];
            self.ki[j] = self.kiv[c % LOOP_NKI];
        }
    }
}

// the sum over the 8 lags in numpy's order: pairwise for one joint (a contiguous column), sequential otherwise
#[inline]
fn sum8<const N: usize>(t: &[f64; NLAG]) -> f64 {
    if N == 1 {
        return ((t[0] + t[1]) + (t[2] + t[3])) + ((t[4] + t[5]) + (t[6] + t[7]));
    }
    let mut s = t[0];
    for x in t.iter().skip(1) { s += *x; }
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    struct Rng(u64);
    impl Rng {
        fn unif(&mut self) -> f64 {
            self.0 ^= self.0 << 13; self.0 ^= self.0 >> 7; self.0 ^= self.0 << 17;
            ((self.0 >> 11) as f64 + 0.5) / 9007199254740992.0
        }
        fn gauss(&mut self) -> f64 { (-2.0 * self.unif().ln()).sqrt() * (2.0 * core::f64::consts::PI * self.unif()).cos() }
    }

    // two joints, 1 kHz, 3 ms frames, a 25 ms servo with 2 ms of delay, a random-walk disturbance, two noisy reads;
    // returns the first tick at which joint 0 has gains and its candidate
    fn closed_loop(law: &mut SyntonicLoopLaw<2>, ticks: u64) -> (u64, usize) { closed_loop_traced(law, ticks, None) }

    fn closed_loop_traced(law: &mut SyntonicLoopLaw<2>, ticks: u64, mut trace: Option<&mut [[f64; 2]]>) -> (u64, usize) {
        let mut rng = Rng(88172645463325252);
        let (mut th, mut dist, mut latched, mut pipe) = ([0.0; 2], [0.0; 2], [0.0; 2], [[0.0; 2]; 3]);
        let mut first = None;
        for k in 0..ticks {
            let t = k as f64 * 1e-3;
            let (mut q, mut qd, mut a, mut b) = ([0.0; 2], [0.0; 2], [0.0; 2], [0.0; 2]);
            for j in 0..2 {
                let f = 0.4 + 0.3 * j as f64;
                q[j] = 0.3 * (2.0 * core::f64::consts::PI * f * t).sin();
                qd[j] = 0.6 * core::f64::consts::PI * f * (2.0 * core::f64::consts::PI * f * t).cos();
                dist[j] += 2e-4 * rng.gauss() + if k == 1500 { 0.01 } else { 0.0 };
                let y = th[j] + dist[j];
                a[j] = y + 2e-4 * rng.gauss();
                b[j] = y + 2e-4 * rng.gauss();
            }
            let send = k % 3 == 0;
            let cmd = *law.step(&q, &qd, &a, &b, send);
            if first.is_none() { if let Some(c) = law.gain_index(0) { first = Some((k, c)); } }
            if let Some(t) = trace.as_deref_mut() { t[k as usize] = cmd; }
            for j in 0..2 {
                if send { latched[j] = cmd[j]; }
                let u = pipe[(k % 3) as usize][j];
                pipe[(k % 3) as usize][j] = latched[j];                     // 2 ticks of delay
                th[j] += 1e-3 * clamp((u - th[j]) / 0.025, -6.0, 6.0);
            }
        }
        first.expect("no gains chosen")
    }

    #[test]
    fn facts() {
        let law: Box<SyntonicLoopLaw<2>> = Box::new(SyntonicLoopLaw::new(&[0.0; 2], &[6.0; 2], LoopLawConfig::default()));
        assert!((law.lead - (0.025 + 0.0015 + 0.0005 + 1e-4)).abs() < 1e-15);       // servo_tau + loop delay
        assert!((law.limit(1) - 0.025 * 6.0).abs() < 1e-15);                         // servo_tau x w0
        let fixed: Box<SyntonicLoopLaw<1>> = Box::new(SyntonicLoopLaw::new(&[0.0], &[6.0], LoopLawConfig { lead: Some(0.02), ..LoopLawConfig::default() }));
        assert_eq!(fixed.lead, 0.02);
    }

    #[test]
    #[should_panic]
    fn delay_line_too_short() {
        let _law: Box<SyntonicLoopLaw<1, 4>> = Box::new(SyntonicLoopLaw::new(&[0.0], &[6.0], LoopLawConfig { dead: 0.003, ..LoopLawConfig::default() }));
    }

    #[test]
    #[should_panic]
    fn spread_below_every() {
        let _law: Box<SyntonicLoopLaw<1>> = Box::new(SyntonicLoopLaw::new(&[0.0], &[6.0], LoopLawConfig { spread: 500, ..LoopLawConfig::default() }));
    }

    #[test]
    fn gains_on_the_robust_grid_and_spread_delays_them_only() {
        let cfg = LoopLawConfig { dead: 0.002, ..LoopLawConfig::default() };
        let mut law: Box<SyntonicLoopLaw<2>> = Box::new(SyntonicLoopLaw::new(&[0.0; 2], &[6.0; 2], cfg));
        let (k0, c0) = closed_loop(&mut law, 4000);
        assert!(law.started);
        for j in 0..2 {
            let c = law.gain_index(j).unwrap();
            assert!(law.robust(c) && law.kp[j] == law.grid_kp(c) && law.ki[j] == law.grid_ki(c));
            assert!(law.gam[j] >= 1.0 && law.gam[j] <= 3.0);
            assert!((law.noise[j] / 2e-8 - 1.0).abs() < 0.2);                        // read-pair noise: sd^2 / 2
        }
        for m in [1u64, 37, 499] {
            let mut s: Box<SyntonicLoopLaw<2>> = Box::new(SyntonicLoopLaw::new(&[0.0; 2], &[6.0; 2], LoopLawConfig { spread: m, ..cfg }));
            assert_eq!(closed_loop(&mut s, 4000), (k0 + m, c0));                      // the same gains, m ticks later
        }
    }

    #[test]
    #[should_panic]
    fn delta_above_rounding() {
        let _law: Box<SyntonicLoopLaw<1>> = Box::new(SyntonicLoopLaw::new(&[0.0], &[6.0], LoopLawConfig { delta: 1e-12, ..LoopLawConfig::default() }));
    }

    #[test]
    fn speedups_give_the_same_commands() {
        let base = LoopLawConfig { dead: 0.002, ..LoopLawConfig::default() };
        let (mut want, mut got) = (vec![[0.0; 2]; 4000], vec![[0.0; 2]; 4000]);
        let mut law: Box<SyntonicLoopLaw<2>> = Box::new(SyntonicLoopLaw::new(&[0.0; 2], &[6.0; 2], LoopLawConfig { prune: false, quad: false, ..base }));
        let (k0, _) = closed_loop_traced(&mut law, 4000, Some(&mut want));                // the scan of 0.3.1
        let recs = (4000 - k0).div_ceil(base.every);
        for (prune, quad, delta) in [(true, false, 1e-9), (false, true, 1e-9), (true, true, 1e-9), (true, true, 1e3)] {
            let mut l: Box<SyntonicLoopLaw<2>> = Box::new(SyntonicLoopLaw::new(&[0.0; 2], &[6.0; 2], LoopLawConfig { prune, quad, delta, ..base }));
            closed_loop_traced(&mut l, 4000, Some(&mut got));
            assert!(got == want, "prune {} quad {} delta {}: the commands differ", prune, quad, delta);
            if delta > 1.0 { assert_eq!(l.fallbacks(), recs); }                                  // every choice fell back
            else { assert_eq!(l.fallbacks(), 0); }
            if prune && delta < 1.0 { assert!(l.points_evaluated() < l.points_total()); }
        }
    }

    #[test]
    fn scheduling() {
        let mut law: Box<SyntonicLoopLaw<2>> = Box::new(SyntonicLoopLaw::new(&[0.0; 2], &[6.0; 2], LoopLawConfig { dead: 0.002, ..LoopLawConfig::default() }));
        let cand = 1500.0 * 70.0 / 600e6;                                                     // 70 cycles a point at 600 MHz
        let m = law.spread_for(cand, 0.00045).unwrap();
        assert!(m > 0 && m < 500);
        assert_eq!(law.spread_for(cand, 10.0), Some(0));
        assert_eq!(law.spread_for(cand, 1e-5), None);
        assert!(law.set_spread(m) && law.spread() == m && !law.set_spread(500) && law.spread() == m);
        let mut t = 0u32;
        let sec = law.time_candidate(|| { t = t.wrapping_add(4_000_000_000); t }, 1e9, 3).unwrap();   // wraps
        assert!((sec - 4.0).abs() < 1e-12);
    }

    #[test]
    fn k_gamma_table() {
        assert_eq!(loop_k_gamma(3.0), 0.55);                                           // 11/20
        assert_eq!(loop_k_gamma(1.0), 0.99999999987);
        assert!((loop_k_gamma(2.75) - 0.5 * (0.52771124824 + 0.55)).abs() < 1e-15);
    }
}
