"""v8: the window floor set by the loop delay. Written 21 Oct 2026, before any run of this file.

Hypothesis (after v6 and v7, where 2 ms was the best floor at every loop rate): the lower bound U' needs is the delay of
the loop, from the encoder read to the servo acting on the correction, not the read period.
Declared from the hardware, nothing tuned: D = T_frame / 2 + dt_c / 2 + t_sensor, with T_frame = 3 ms (servo PWM frame,
zero-order hold), dt_c = controller tick (age of the latest correction at a frame, on average), t_sensor = 0.1 ms
(AS5048 system propagation delay, data sheet). D = 1.85 ms at 2 kHz, 2.10 at 1 kHz, 2.60 at 500 Hz, 3.10 at 333 Hz.

Runs: U' derived (v5 code, exactly as run in v6 at 2 kHz and v7 at 1 kHz and below) with floor D, and a floor sweep
{1, 1.5, 2, 2.5, 3, 4} ms to locate the best floor at each rate. Baselines: the PID and fixed window of v6 / v7w, unchanged.
Built-in check: the 2 ms row must reproduce v6 (2 kHz) and v7 (500, 333 Hz), and v5c at 1 kHz.

Predictions:
F1  Floor = D is within 5 % (tip rms) of the best floor of the sweep in every condition at every rate.
F2  The best floor of the sweep (lowest mean square over the four conditions) rises as the rate falls, as D does.
F3  At 1 kHz, floor = D (2.1 ms) stays within 3 % of v5c.
A failed line is reported as such."""
import os, sys, json, time, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
HZ = int(os.environ['BENCH_HZ'])
if HZ > 1000:
    import flow_v6_rate as M
    DTC = M.DT
    def sim(tt, q, tau, sd, seed, fl): return M.simulate_u(np.ones((1, 12)), tt, q, tau, sd, seed, fl)
else:
    import flow_v7_ctrl as M
    DTC = M.DTC
    def sim(tt, q, tau, sd, seed, fl): return M.simulate('uprime', dict(km=np.ones((1, 12))), tt, q, tau, sd, seed, fl, store=True)
F = M.F
D = 0.0015 + DTC / 2 + 0.0001
FLOORS = [('D', D)] + [(f'{x:g} ms', x * 1e-3) for x in (1, 1.5, 2, 2.5, 3, 4)]

if __name__ == '__main__':
    t_start = time.time(); out = dict(hz=HZ, D_ms=D * 1e3, results={})
    print(f'{HZ} Hz: tick {DTC * 1e3:g} ms, D = {D * 1e3:.2f} ms', flush=True)
    for name, t0 in (('A', 6.0), ('B', 217.5)):
        for cond in ('real', 'step'):
            tt, q, tau = F.plan(t0 - 3.0, t0 + 30.0)
            sd = np.full(len(tt), F.ENC_SD)
            if cond == 'step': sd[(tt >= t0 + 10) & (tt < t0 + 20)] *= np.sqrt(10.0)
            seed = 7 if cond == 'step' else 5
            res = {}
            for lab, fl in FLOORS:
                r = sim(tt, q, tau, sd, seed, fl); k = r['t'] >= t0 - 1e-9
                s_, e = F.score(r['th'][k], r['t'][k], t0); tau_ = r['tau'][k][:, :6]
                res[lab] = dict(floor_ms=fl * 1e3, rms=s_['tip_rms_mm'], jit=s_['jitter_mm'],
                                on_floor=float(np.mean(tau_ <= fl * 1.001)), tau_median_ms=float(1e3 * np.median(tau_)))
            out['results'][f'{name} {cond}'] = res
            print(f'RESULT {name} {cond}', {lab: (round(v['rms'], 2), round(v['jit'], 2)) for lab, v in res.items()}, round(time.time() - t_start), 's', flush=True)
            json.dump(out, open(f'{F.OUT}/summary_v8_{HZ}Hz.json', 'w'), indent=1)
    print('done', round(time.time() - t_start), 's', flush=True)
