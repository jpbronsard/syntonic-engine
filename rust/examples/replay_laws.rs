// Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
// Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
//
// Replays a log written by tools/make_law_replay_log.py through the Rust loop law or attitude law and compares every
// output with the Python law's. Prints, per output, the largest absolute difference and that difference relative to the
// output's largest magnitude over the run; for the loop law also the ticks where the chosen gains differ (none
// expected) and the time per tick. From the rust/ folder, after `python tools/make_law_replay_log.py` in the root:
//   cargo run --release --example replay_laws -- ../replay_loop_bench.bin ../replay_loop_bench_spread.bin
//   cargo run --release --example replay_laws -- ../replay_loop_two.bin ../replay_loop_one.bin ../replay_att_synth.bin
// Exit status 0 when every log agrees: the same gains at every tick and every output within 1e-12 of its scale.
#![allow(clippy::needless_range_loop, clippy::neg_cmp_op_on_partial_ord, clippy::manual_is_multiple_of)]
use std::time::Instant;
use std::{env, fs, process};
use syntonic_engine::{AttitudeConfig, AttitudeLevel, LoopLawConfig, RiccatiAttitudeLaw, SyntonicLoopLaw, LOOP_NC};

fn read_f64(path: &str) -> Vec<f64> {
    fs::read(path).expect("cannot read file").chunks_exact(8).map(|c| f64::from_le_bytes(c.try_into().unwrap())).collect()
}

struct Group { name: &'static str, maxd: f64, scale: f64 }

impl Group {
    fn new(name: &'static str) -> Self { Group { name, maxd: 0.0, scale: 0.0 } }
    fn add(&mut self, got: f64, want: f64) {
        let d = (got - want).abs();
        if d > self.maxd { self.maxd = d; }
        if want.abs() > self.scale { self.scale = want.abs(); }
    }
    fn rel(&self) -> f64 { if self.scale > 0.0 { self.maxd / self.scale } else { self.maxd } }
}

fn report(g: &[Group]) -> bool {
    let mut ok = true;
    for x in g {
        println!("    {:<8} max |Rust - Python| {:.3e}   relative to its scale {:.3e}", x.name, x.maxd, x.rel());
        if !(x.rel() <= 1e-12) { ok = false; }
    }
    ok
}

fn arr<const N: usize>(s: &[f64]) -> [f64; N] { let mut a = [0.0; N]; a.copy_from_slice(&s[..N]); a }

fn loop_replay<const N: usize>(v: &[f64], path: &str) -> bool {
    let h = v;
    let t_len = h[2] as usize;
    let cfg = LoopLawConfig {
        dt: h[3], frame: h[4], servo_tau: h[5], dead: h[6], lead: if h[7] < 0.0 { None } else { Some(h[7]) },
        spread: h[8] as u64, every: h[9] as u64, ..LoopLawConfig::default()
    };
    let (w0, q0) = (arr::<N>(&v[16..]), arr::<N>(&v[16 + N..]));
    let rows = &v[16 + 2 * N..];
    let c = 15 * N + 1;
    if rows.len() != t_len * c { println!("{}: bad size", path); return false; }
    let t0 = Instant::now();
    let mut law: Box<SyntonicLoopLaw<N>> = Box::new(SyntonicLoopLaw::new(&q0, &w0, cfg));
    let t_begin = t0.elapsed().as_secs_f64();
    let names = ["cmd", "Kp", "Ki", "I", "R", "A", "gam", "Oc", "th_m", "th_ff"];
    let mut g: Vec<Group> = names.iter().map(|n| Group::new(n)).collect();
    let (mut idx_bad, mut nrec, mut n_tick) = (0usize, 0usize, 0usize);
    let (mut t_tick, mut t_rec, mut t_max) = (0.0f64, 0.0f64, 0.0f64);
    for row in rows.chunks_exact(c) {
        let (q, qd, a, b) = (arr::<N>(row), arr::<N>(&row[N..]), arr::<N>(&row[2 * N..]), arr::<N>(&row[3 * N..]));
        let send = row[4 * N] != 0.0;
        let out = &row[4 * N + 1..];
        let s = Instant::now();
        let cmd = *law.step(&q, &qd, &a, &b, send);
        let dts = s.elapsed().as_secs_f64();
        if dts > t_max { t_max = dts; }
        if cfg.spread == 0 && law.started && (law.tick() - 1) % cfg.every == 0 { t_rec += dts; nrec += 1; }
        else { t_tick += dts; n_tick += 1; }
        for j in 0..N {
            g[0].add(cmd[j], out[j]); g[1].add(law.kp[j], out[N + j]); g[2].add(law.ki[j], out[2 * N + j]);
            g[3].add(law.integral[j], out[3 * N + j]); g[4].add(law.noise[j], out[4 * N + j]);
            g[5].add(law.amp[j], out[5 * N + j]); g[6].add(law.gam[j], out[6 * N + j]); g[7].add(law.oc[j], out[7 * N + j]);
            g[8].add(law.th_m[j], out[8 * N + j]); g[9].add(law.th_ff[j], out[9 * N + j]);
            let gi = law.gain_index(j).map(|x| x as f64).unwrap_or(-1.0);
            if gi != out[10 * N + j] { idx_bad += 1; }
        }
    }
    println!("{}: loop law, {} joint(s), {} ticks, dead {:.1} ms, spread {}, {} robust candidates of {}", path, N, t_len,
             cfg.dead * 1e3, cfg.spread, law.robust_count(), LOOP_NC);
    println!("    chosen gains differ from Python's on {} joint-ticks", idx_bad);
    let ok = report(&g) && idx_bad == 0;
    print!("    time: new {:.1} ms; step {:.2} us per ordinary tick", t_begin * 1e3, 1e6 * t_tick / n_tick.max(1) as f64);
    if nrec > 0 { print!(", {:.2} ms per gain-update tick ({})", 1e3 * t_rec / nrec as f64, nrec); }
    println!(", longest tick {:.3} ms", t_max * 1e3);
    ok
}

fn att_replay<const L: usize>(v: &[f64], path: &str) -> bool {
    let h = v;
    let t_len = h[2] as usize;
    let cfg = AttitudeConfig {
        dt: h[3], dt_accel: h[4], ground: h[7] as usize, settle: h[8] as u64, ema: h[9],
        block_scales: [h[10] as u64, h[11] as u64, h[12] as u64],
    };
    let mut lv = [AttitudeLevel { s2: 0.0, lam: 0.0, t_rec: 0.0 }; L];
    for r in 0..L { lv[r] = AttitudeLevel { s2: h[16 + r], lam: h[16 + L + r], t_rec: h[16 + 2 * L + r] }; }
    let rows = &v[16 + 3 * L..];
    let s = 1 + L;
    let c = 4 + s + s * s + L + 1;
    if rows.len() != t_len * c { println!("{}: bad size", path); return false; }
    let mut law = RiccatiAttitudeLaw::<L>::new(&lv, h[5], h[6], cfg);
    let mut g: Vec<Group> = ["theta", "biases", "P", "reading", "thr"].iter().map(|n| Group::new(n)).collect();
    let mut arms = 0;
    for row in rows.chunks_exact(c) {
        let a = if row[1].is_nan() { None } else { Some(row[1]) };
        let th = law.step(row[0], a, row[2] as usize);
        if row[3] != 0.0 { law.arm(); arms += 1; }
        let out = &row[4..];
        g[0].add(th, out[0]);
        for i in 1..s { g[1].add(law.x(i), out[i]); }
        for i in 0..s { for m in 0..s { g[2].add(law.p(i, m), out[s + i * s + m]); } }
        for i in 0..L { g[3].add(law.noise_reading(i), out[s + s * s + i]); }
        g[4].add(law.thr, out[s + s * s + L]);
    }
    println!("{}: attitude law, {} levels, {} gyro samples, armed {} time(s)", path, L, t_len, arms);
    report(&g)
}

fn main() {
    let args: Vec<String> = env::args().skip(1).collect();
    if args.is_empty() { eprintln!("usage: replay_laws replay_loop_*.bin replay_att_*.bin ..."); process::exit(2); }
    let mut all = true;
    for path in &args {
        let v = read_f64(path);
        let (kind, n) = (v[0] as i32, v[1] as usize);
        let ok = match (kind, n) {
            (1, 1) => loop_replay::<1>(&v, path),
            (1, 2) => loop_replay::<2>(&v, path),
            (1, 6) => loop_replay::<6>(&v, path),
            (2, 3) => att_replay::<3>(&v, path),
            _ => { println!("{}: kind {} with {} joints or levels is not compiled in this example", path, kind, n); false }
        };
        println!("    {}", if ok { "agrees" } else { "DIFFERS" });
        all &= ok;
    }
    process::exit(if all { 0 } else { 1 });
}
