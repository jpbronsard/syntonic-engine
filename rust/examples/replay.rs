// Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
// Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
//
// Replays a log written by tools/make_replay_log.py through the Rust engine and compares every command with the
// Python engine's. From the rust/ folder:
//   cargo run --release --example replay -- ../replay_log.bin ../replay_w0.bin raw              (the 0.1.x identification)
//   cargo run --release --example replay -- ../replay_log_allan.bin ../replay_w0.bin allan      (the block-mean one, default)
//   cargo run --release --example replay -- ../replay_log_lead.bin ../replay_w0.bin raw-lead    (0.2.9: lead and e travel on)
//   cargo run --release --example replay -- ../replay_log_lead_allan.bin ../replay_w0.bin allan-lead
use std::{env, fs};
use syntonic_engine::{Config, Engine, Identification, Travel};

fn read_f64(path: &str) -> Vec<f64> {
    fs::read(path).expect("cannot read file").chunks_exact(8).map(|c| f64::from_le_bytes(c.try_into().unwrap())).collect()
}

fn main() {
    const N: usize = 12;
    let args: Vec<String> = env::args().collect();
    let log = read_f64(&args[1]);
    let w0v = read_f64(&args[2]);
    let mode = if args.len() < 4 { "raw" } else { args[3].as_str() };
    let raw = mode == "raw" || mode == "raw-lead";
    let lead = mode == "raw-lead" || mode == "allan-lead";
    let (cols, ia) = if lead { (5 * N, 2 * N) } else { (4 * N, N) };
    let rows: Vec<&[f64]> = log.chunks_exact(cols).collect();
    let mut q0 = [0.0; N]; q0.copy_from_slice(&rows[0][..N]);
    let mut w0 = [0.0; N]; w0.copy_from_slice(&w0v[..N]);
    let mut cfg = Config { identification: if raw { Identification::Raw } else { Identification::Allan }, ..Config::default() };
    if lead { cfg.lead = -1.0; cfg.travel = Travel::E; }
    let mut eng: Box<Engine<N, 65>> = Box::new(Engine::new(&q0, &w0, cfg));
    let (mut maxd, mut worst_tick) = (0.0f64, 0usize);
    for (t, row) in rows.iter().enumerate() {
        let (mut q, mut qd, mut a, mut b) = ([0.0; N], [0.0; N], [0.0; N], [0.0; N]);
        q.copy_from_slice(&row[..N]); a.copy_from_slice(&row[ia..ia + N]); b.copy_from_slice(&row[ia + N..ia + 2 * N]);
        let cmd = if lead { qd.copy_from_slice(&row[N..2 * N]); eng.step_lead(&q, &qd, &a, &b, None) } else { eng.step(&q, &a, &b, None) };
        for j in 0..N {
            let d = (cmd[j] - row[ia + 2 * N + j]).abs();
            if d > maxd { maxd = d; worst_tick = t; }
        }
    }
    println!("{} ticks, mode {}, floor {:.3} ms, max |cmd Rust - cmd Python| = {:.3e} rad (tick {}), final windows (ms): {:?}",
             rows.len(), mode, eng.floor() * 1e3, maxd, worst_tick,
             eng.tau.iter().map(|x| (x * 1e5).round() / 100.0).collect::<Vec<_>>());
}
