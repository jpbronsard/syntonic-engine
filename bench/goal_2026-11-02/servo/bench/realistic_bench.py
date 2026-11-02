"""Realistic simulation bench for one arm of MAESTRO 1 (the baton arm, RJ1-RJ6).

Built on examples/servo_bench.py of syntonic-engine 0.2.9 (branch release-0.2.9): the same planned motion
(excerpts A and B of examples/data/maestro_excerpts.npz, Beethoven's Fifth), the same encoder, compliance, gear-play
and speed-derating models and the same scores, with the hardware the author fixed on 2 Nov 2026 and a set of
realism options. At the published plant options (see PUBLISHED_PLANT) and 12 joints it reproduces
examples/servo_bench.py bit for bit (run_baselines.py gate).

HARDWARE (author, 2 Nov 2026)
  Control loop 1 kHz (Teensy CAN node); servo command frames every 3 ms (333 Hz; the DSSERVO BLS servos accept
  50-330 Hz PWM), held between frames. RJ1, RJ2: DSSERVO BLS53180 (24 V); RJ3, RJ5: BLS3355; RJ4: BLS3375;
  RJ6: BLS3325. Their stall torque and no-load speed are those of the data file. Servo deadband set to 0 by the author.

PLANT, per joint (vectorised over columns; a column is one joint driven by one controller)
  load      planned load torque (inverse dynamics of the planned motion) x load_scale.
  command   latched on frame ticks only (every 3rd tick from tick 0), rounded to cmd_quant_deg, and seen by the servo
            dead_time later. A dead time that is not a whole number of ticks gives, over the tick, the time average
            of the delayed held command.
  servo     servo_order=1 (published): th += dt * clip((u - th) / servo_tau, -w_eff, w_eff), explicit Euler at the
              loop period exactly as published. Its 63 % step rise is 24.5 ms for servo_tau = 25 ms (discretisation);
              its low-frequency lag is exactly servo_tau.
            servo_order=2: th'' = wn^2 (u - th) - 2 zeta wn th', wn set so that the 63 % rise time of the continuous
              step response equals servo_tau (wn = 65.5 rad/s for 25 ms at zeta 0.6). Integrated exactly over each
              tick for a held input, then the position step is clipped at w_eff dt and the velocity at w_eff.
              NOTE: its low-frequency lag is 2 zeta / wn (18.3 ms at 25 ms, zeta 0.6), shorter than servo_tau.
            speed limit (published derating): w_eff = w0 x max(0.05, 1 - |load| / stall torque).
            deadband: the servo does not act while |u - th| <= deadband_deg / 2 (full width; hard deadband).
  joint     th - load x (1 deg / stall torque) - 0.15 deg x tanh(load / (0.02 x stall torque))
            (published: 1 deg of compliance at stall torque, 0.3 deg of gear play).
  encoder   (published, AS5048 on the joint output) two reads per tick; each read = joint + N(0, sd^2), rounded to
            0.0219 deg (14 bits); sd = 0.06 deg x sqrt(variance factor of the circumstance's noise windows).
            Noise draws: per tick, one standard normal per joint for read a, then one per joint for read b, shared
            by every column on that joint (paired comparison), as in pid_lead_fair.py.

  The plant options and their defaults are PLANT_DEFAULTS. Assumption stated: the default command quantization
  0.025 deg = about 0.25 us of PWM resolution at about 0.1 deg/us (not measured on the BLS servos).

CONTROLLER INTERFACE
  A controller is built by a factory, factory(q0, w0) -> controller: q0 = planned positions at the first tick (rad),
  w0 = data-sheet no-load speeds (rad/s), one per joint (RJ1..RJ6). The bench's other facts are module constants:
  DT = 1 ms, FRAME = 3 ms, TS_NOM = 25 ms (the servo time constant controllers may assume; not the true one).
  A controller never sees the circumstance. Every tick:
      command = controller.step(q_plan, qd_plan, read_a, read_b, send)
  q_plan, qd_plan: planned positions (rad) and velocities (rad/s) now; read_a, read_b: the two encoder reads (rad);
  send: True on frame ticks. The plant latches `command` (rad, shape (n,)) only when send is True. Angles in rad.
  Where a rival needs per-excerpt gains (the PIDs), pass {'A': factory_A, 'B': factory_B} instead of one factory.

SCORES (published definitions, examples/servo_bench.py), over the 30 s of the excerpt (3 s of run-in before it):
  tip_mm     rms distance between the actual and the planned baton tip (arm_fk.tip, joints RJ1-RJ6), at 100 Hz;
  jitter_mm  rms of that error above 5 Hz (4th-order Butterworth high-pass, forward-backward);
  joint_deg  rms of each joint's angle error (plan - joint output), averaged over RJ1-RJ6.

Usage:
  import realistic_bench as RB, rivals as RV
  T = RV.load_tuned('tuned_gains.json')
  res = RB.simulate('A', 'servo_slow', 91, [('PID+lead+meanP', RV.pid_factory(T['leadmean'], lead=True, p_mean=True)),
                                            ('mine', lambda q0, w0: MyController(q0, w0))])
  res['mine'] -> {'tip_mm': ..., 'jitter_mm': ..., 'joint_deg': ..., 'joint_rms_deg': array (6,)}

Seeds 41-45 and 71-75 are reserved for the author's confirmation runs: simulate() refuses them unless
allow_reserved=True.
"""
import os
import sys

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.linalg import expm
from scipy.optimize import brentq
from scipy.signal import butter, sosfiltfilt

REPO = os.path.abspath(os.environ.get('SYNTONIC_REPO', '/home/claude/repo/syntonic-engine'))
for _p in (os.path.join(REPO, 'python'), os.path.join(REPO, 'examples')):
    if _p not in sys.path:
        sys.path.insert(0, _p)
from arm_fk import tip  # noqa: E402  (the published tip kinematics, unchanged)

__all__ = ['simulate', 'simulate_columns', 'resolve', 'CIRCUMSTANCES', 'GATE_CIRCUMSTANCES',
           'PLANT_DEFAULTS', 'PUBLISHED_PLANT', 'NOISE_SCHEDULES', 'RESERVED_SEEDS', 'DT', 'FRAME', 'TS_NOM',
           'W0', 'STALL', 'JOINT_NAMES', 'Z', 'rise63_x', 'second_order_wn', 'step_response']

# ---------------------------------------------------------------------------------------------- published constants
DEG = np.pi / 180
Z = np.load(os.path.join(REPO, 'examples', 'data', 'maestro_excerpts.npz'))
JOINT_NAMES = tuple(str(j) for j in Z['joints'])           # RJ1..RJ6, LJ1..LJ6
STALL = Z['stall_torque']                                   # N m, data sheet
W0 = Z['no_load_speed']                                     # rad/s, data sheet
ENC_SD, ENC_LSB = float(Z['enc_sd_rad']), float(Z['enc_lsb_rad'])   # 0.06 deg per read, 0.0219 deg steps
K_INV = (1.0 * DEG) / STALL                                 # compliance: 1 deg at stall torque
HB = 0.5 * 0.3 * DEG                                        # half of the 0.3 deg gear play
TS_NOM = 0.025                                              # servo time constant the controllers assume (s)
DT, FRAME = 1e-3, 3e-3                                      # loop period, servo frame period (s)
RESERVED_SEEDS = frozenset(range(41, 46)) | frozenset(range(71, 76))

# ---------------------------------------------------------------------------------------------- options and battery
NOISE_SCHEDULES = {
    'steady': (),
    # (start, end, variance factor), seconds from the excerpt start (published noise step: x10, 10 to 20 s)
    'step': ((10.0, 20.0, 10.0),),
    # x30 for 0.5 s every 2.5 s from 5 to 25 s: onsets 5.0, 7.5, ..., 22.5 s (8 bursts, the last ends at 23 s)
    'bursts': tuple((5.0 + 2.5 * i, 5.5 + 2.5 * i, 30.0) for i in range(8)),
}

PLANT_DEFAULTS = dict(
    servo_order=1,        # 1: first order (published model); 2: second order
    zeta=0.6,             # damping ratio of the second-order servo
    servo_tau=0.025,      # true servo time constant (s), scalar or one per joint; 2nd order: its 63 % rise time
    dead_time=0.001,      # internal servo dead time (s), from the frame latch to the servo acting on it
    cmd_quant_deg=0.025,  # command quantization (deg), 0 = none; ASSUMPTION: ~0.25 us PWM step x ~0.1 deg/us
    deadband_deg=0.0,     # servo deadband, full width (deg); set to 0 on the robot by the author
    load_scale=1.0,       # multiplies the planned load torques
    noise='steady',       # a key of NOISE_SCHEDULES or a tuple of (start, end, variance factor) windows
)

# the published plant (examples/servo_bench.py): used by the gate only
PUBLISHED_PLANT = dict(PLANT_DEFAULTS, dead_time=0.0, cmd_quant_deg=0.0)

CIRCUMSTANCES = {
    'nominal': dict(),                                   # steady noise, nominal servo
    'noise_step': dict(noise='step'),                    # encoder read variance x10 from 10 to 20 s (published)
    'bursts': dict(noise='bursts'),                      # variance x30 for 0.5 s every 2.5 s from 5 to 25 s
    'servo_slow': dict(servo_tau=0.035),                 # true time constant 35 ms on every joint
    'servo_fast': dict(servo_tau=0.018),                 # 18 ms
    'servo_2nd': dict(servo_order=2, zeta=0.6),          # second order, same 63 % rise time
    'load_x1.5': dict(load_scale=1.5),                   # loads x1.5
    'quant_coarse': dict(cmd_quant_deg=0.09),            # 0.09 deg command steps (~1 us at ~0.09 deg/us)
    # Added: a PWM servo only knows the pulse width at the pulse's falling edge, 0.5-2.5 ms after the frame starts
    # (1.5 ms at mid-travel), before its own processing; 1 ms internal + 1.5 ms of pulse = 2.5 ms.
    'deadtime_long': dict(dead_time=0.0025),
}

GATE_CIRCUMSTANCES = {
    'published_steady': dict(PUBLISHED_PLANT),
    'published_step': dict(PUBLISHED_PLANT, noise='step'),
}


def resolve(circumstance):
    """Plant options of a circumstance: a name (CIRCUMSTANCES or GATE_CIRCUMSTANCES) or a dict of overrides of
    PLANT_DEFAULTS. Returns the full options dict, with the noise schedule as a tuple of windows."""
    if isinstance(circumstance, str):
        if circumstance in CIRCUMSTANCES:
            over = CIRCUMSTANCES[circumstance]
        elif circumstance in GATE_CIRCUMSTANCES:
            over = GATE_CIRCUMSTANCES[circumstance]
        else:
            raise KeyError(f'unknown circumstance {circumstance!r}; known: {sorted(CIRCUMSTANCES)}')
    else:
        over = dict(circumstance)
    bad = set(over) - set(PLANT_DEFAULTS)
    if bad:
        raise KeyError(f'unknown plant options {sorted(bad)}; known: {sorted(PLANT_DEFAULTS)}')
    P = dict(PLANT_DEFAULTS, **over)
    if isinstance(P['noise'], str):
        P['noise'] = NOISE_SCHEDULES[P['noise']]
    P['noise'] = tuple(tuple(float(x) for x in w) for w in P['noise'])
    if P['servo_order'] not in (1, 2):
        raise ValueError('servo_order must be 1 or 2')
    if not np.all(np.asarray(P['servo_tau'], float) > 0):
        raise ValueError('servo_tau must be positive')
    if P['dead_time'] < 0 or P['cmd_quant_deg'] < 0 or P['deadband_deg'] < 0 or P['load_scale'] < 0:
        raise ValueError('dead_time, cmd_quant_deg, deadband_deg and load_scale must be >= 0')
    if P['servo_order'] == 2 and not P['zeta'] > 0:
        raise ValueError('zeta must be positive')
    return P


# ---------------------------------------------------------------------------------------------- second-order servo
def _step2(x, zeta):
    """Unit-step response at normalised time x = wn t of wn^2 / (s^2 + 2 zeta wn s + wn^2)."""
    if zeta < 1.0:
        wd = np.sqrt(1.0 - zeta * zeta)
        return 1.0 - np.exp(-zeta * x) * (np.cos(wd * x) + zeta / wd * np.sin(wd * x))
    if zeta == 1.0:
        return 1.0 - np.exp(-x) * (1.0 + x)
    r = np.sqrt(zeta * zeta - 1.0); s1, s2 = -zeta + r, -zeta - r
    return 1.0 - (s2 * np.exp(s1 * x) - s1 * np.exp(s2 * x)) / (s2 - s1)


def rise63_x(zeta):
    """wn x (63 % rise time) of the second-order step response: the first x with response 1 - 1/e."""
    target = 1.0 - np.exp(-1.0)
    hi = np.pi / np.sqrt(1.0 - zeta * zeta) if zeta < 1.0 else 50.0 * zeta   # first peak (> 1) or far out
    return brentq(lambda x: _step2(x, zeta) - target, 1e-9, hi, xtol=1e-14)


def second_order_wn(servo_tau, zeta):
    """Natural frequency (rad/s) such that the 63 % rise time of the step response equals servo_tau."""
    return rise63_x(zeta) / np.asarray(servo_tau, float)


def _second_order_coeffs(wn, zeta, dt):
    """Exact one-tick transition of (th - u, th') for a held input u; returns (A00 - 1, A01, A10, A11) per joint."""
    out = np.zeros((4, len(wn)))
    for j, w in enumerate(wn):
        A = expm(np.array([[0.0, 1.0], [-w * w, -2.0 * zeta * w]]) * dt)
        out[:, j] = (A[0, 0] - 1.0, A[0, 1], A[1, 0], A[1, 1])
    return out


# ---------------------------------------------------------------------------------------------- plan and noise
def _plan(name, joints):
    """Planned positions, velocities and load torques on the loop grid, as pid_lead_fair.py builds them (the
    12-joint gate reproduces examples/servo_bench.py exactly with this construction)."""
    t0 = float(Z[f'{name}_t0'])
    tt = np.arange(t0 - 3.0, t0 + 30.0, DT)
    S = CubicSpline(Z[f'{name}_t30'], Z[f'{name}_q30'][:, joints], axis=0)
    q = S(tt); qd = S(tt, 1)
    tau = CubicSpline(Z[f'{name}_t30'], Z[f'{name}_tau30'][:, joints], axis=0)(tt)
    return t0, tt, q, qd, tau


def _noise_sd(tt, t0, windows):
    sd = np.full(len(tt), ENC_SD)
    for a, b, f in windows:
        sd[(tt >= t0 + a) & (tt < t0 + b)] *= np.sqrt(f)
    return sd


# ---------------------------------------------------------------------------------------------- seeds
def _check_seed(seed, allow_reserved):
    if int(seed) in RESERVED_SEEDS and not allow_reserved:
        raise ValueError(f'seed {seed} is reserved for the confirmation runs (41-45, 71-75); '
                         'pass allow_reserved=True only for those runs')


# ---------------------------------------------------------------------------------------------- simulation core
def simulate_columns(name, circumstance, seed, blocks, n_joints=6, tip_score=True, allow_reserved=False):
    """Low-level simulation: every block is a controller over its own columns, all on one shared plant model and
    one shared set of encoder noise draws.

    name          excerpt, 'A' or 'B'
    circumstance  a name of CIRCUMSTANCES / GATE_CIRCUMSTANCES, or a dict of plant options (see PLANT_DEFAULTS)
    seed          noise seed (numpy default_rng)
    blocks        list of (label, factory, cols): cols = joint index of each column (None: joints 0..n_joints-1);
                  factory(q0, w0) gets the block's columns' first planned positions and no-load speeds, or is a
                  dict {excerpt: factory}
    n_joints      joints simulated (6: the baton arm; 12: both arms, for the gate)
    tip_score     compute the tip scores for the blocks that drive a whole arm (cols = 0..n_joints-1)

    Returns {label: {'joint_rms_deg': per-column rms angle error (deg),
                     and for whole-arm blocks: 'joint_deg' (mean over RJ1-RJ6), 'tip_mm', 'jitter_mm'}}.
    """
    _check_seed(seed, allow_reserved)
    P = resolve(circumstance)
    J = int(n_joints)
    joints = list(range(J))
    t0, tt, q, qd, tau = _plan(name, joints)
    T = len(tt)
    sd = _noise_sd(tt, t0, P['noise'])
    tS, w0J, kinv = STALL[:J], W0[:J], K_INV[:J]
    # ---- controllers
    labels, ctrls, cols, slices, full = [], [], [], [], []
    off = 0
    for label, factory, cj in blocks:
        cj = np.arange(J) if cj is None else np.asarray(cj, dtype=int)
        if cj.ndim != 1 or len(cj) == 0 or cj.min() < 0 or cj.max() >= J:
            raise ValueError(f'block {label!r}: bad column map')
        fac = factory[name] if isinstance(factory, dict) else factory
        labels.append(label); ctrls.append(fac(q[0][cj].copy(), w0J[cj].copy())); cols.append(cj)
        slices.append(slice(off, off + len(cj))); off += len(cj)
        full.append(len(cj) == J and bool(np.all(cj == np.arange(J))))
    if len(set(labels)) != len(labels):
        raise ValueError('controller labels must be unique')
    jc = np.concatenate(cols); C = len(jc); nblk = len(ctrls)
    # ---- plant options
    ts_true = np.broadcast_to(np.asarray(P['servo_tau'], float), (J,)).astype(float)
    tsc = ts_true[jc]
    order = P['servo_order']
    if order == 2:
        co = _second_order_coeffs(second_order_wn(ts_true, P['zeta']), P['zeta'], DT)[:, jc]
        a00m1, a01, a10, a11 = co
    qc = P['cmd_quant_deg'] * DEG
    db_half = 0.5 * P['deadband_deg'] * DEG
    ls = float(P['load_scale'])
    nd = P['dead_time'] / DT
    n_int = int(np.floor(nd + 1e-9)); frac = nd - n_int
    if frac < 1e-9:
        frac = 0.0
    delay = n_int > 0 or frac > 0.0
    # ---- state
    th = q[0][jc].copy()                    # servo output (gear) angle
    v = np.zeros(C)                         # servo velocity (second order)
    latched = q[0][jc].copy()               # command held by the servo (after quantization)
    L = n_int + 2
    hist = np.tile(latched, (L, 1)) if delay else None
    rng = np.random.default_rng(seed)
    fe = max(1, int(round(FRAME / DT)))
    every = int(round(0.01 / DT))
    err2 = np.zeros(C); m = 0
    keep, keep_q = [], []
    for k in range(T):
        load = tau[k] if ls == 1.0 else tau[k] * ls
        c1 = load * kinv                                        # compliance
        c2 = HB * np.tanh(load / (0.02 * tS))                   # gear play taken up by the load
        w_eff = w0J * np.maximum(0.05, 1.0 - np.abs(load) / tS)
        out = th - c1[jc] - c2[jc]                              # joint output angle
        na = rng.normal(0.0, 1.0, J); nb = rng.normal(0.0, 1.0, J)
        ra = np.round((out + na[jc] * sd[k]) / ENC_LSB) * ENC_LSB
        rb = np.round((out + nb[jc] * sd[k]) / ENC_LSB) * ENC_LSB
        send = (k % fe == 0)
        qk = q[k]; qdk = qd[k]
        for i in range(nblk):
            s = slices[i]; cj = cols[i]
            cmd = ctrls[i].step(qk[cj], qdk[cj], ra[s], rb[s], send)
            if send:
                latched[s] = cmd
        if send and qc > 0.0:
            latched = np.round(latched / qc) * qc
        if delay:
            hist[k % L] = latched
            u = hist[(k - n_int) % L]
            if frac > 0.0:
                u = (1.0 - frac) * u + frac * hist[(k - n_int - 1) % L]
        else:
            u = latched
        e_s = u - th
        if db_half > 0.0:
            e_s = np.where(np.abs(e_s) <= db_half, 0.0, e_s)
        we = w_eff[jc]
        if order == 1:
            th = th + DT * np.clip(e_s / tsc, -we, we)
        else:
            y0 = -e_s                                           # th - u (deadbanded)
            dth = a00m1 * y0 + a01 * v
            v = np.clip(a10 * y0 + a11 * v, -we, we)
            th = th + np.clip(dth, -we * DT, we * DT)
        if tt[k] >= t0:
            err2 += (qk[jc] - out) ** 2; m += 1
        if tip_score and k % every == 0 and tt[k] >= t0 - 1e-9:
            keep.append(out.copy()); keep_q.append(qk[:6].copy())
    # ---- scores
    rms = np.sqrt(err2 / m) / DEG
    res = {}
    tq = None
    if tip_score and any(full) and J >= 6:
        tq = np.array([tip(b) for b in keep_q])
        sos = butter(4, 5.0, 'highpass', fs=100.0, output='sos')
    for i, lab in enumerate(labels):
        s = slices[i]
        r = {'joint_rms_deg': rms[s].copy()}
        if full[i] and J >= 6:
            r['joint_deg'] = float(np.mean(rms[s][:6]))
            if tq is not None:
                e = np.array([tip(row[s][:6]) for row in keep]) - tq
                hp = sosfiltfilt(sos, e, axis=0)
                r['tip_mm'] = float(1e3 * np.sqrt(np.mean(np.sum(e ** 2, 1))))
                r['jitter_mm'] = float(1e3 * np.sqrt(np.mean(np.sum(hp ** 2, 1))))
        res[lab] = r
    return res


def simulate(name, circumstance, seed, controllers, n_joints=6, allow_reserved=False):
    """Paired simulation of several controllers on one excerpt, circumstance and noise seed.

    controllers  list of (label, factory); factory(q0, w0) builds an object with
                 step(q_plan, qd_plan, read_a, read_b, send) -> command, over the n_joints joints (rad);
                 a dict {excerpt: factory} may stand for a factory (per-excerpt gains).
    Returns {label: {'tip_mm', 'jitter_mm', 'joint_deg', 'joint_rms_deg'}}: tip rms and tip jitter (mm, the published
    definitions, baton arm), mean over RJ1-RJ6 of the joints' rms angle error (deg), and the per-joint rms (deg).
    """
    if n_joints < 6:
        raise ValueError('the tip score needs the six baton-arm joints')
    return simulate_columns(name, circumstance, seed, [(lab, fac, None) for lab, fac in controllers],
                            n_joints=n_joints, tip_score=True, allow_reserved=allow_reserved)


# ---------------------------------------------------------------------------------------------- servo diagnostics
def step_response(circumstance='nominal', joint=0, amplitude_deg=1.0, duration=0.2):
    """Open-loop response of one servo of the plant to a commanded step at t = 0 (command latched at tick 0, no load),
    for checking the servo options. Returns (t, response / amplitude)."""
    P = resolve(circumstance)
    ts = float(np.broadcast_to(np.asarray(P['servo_tau'], float), (6,))[joint])
    w = float(W0[joint]); a = amplitude_deg * DEG
    n = int(round(duration / DT)); nd = P['dead_time'] / DT
    n_int = int(np.floor(nd + 1e-9)); frac = nd - n_int if nd - n_int > 1e-9 else 0.0
    if P['servo_order'] == 2:
        a00m1, a01, a10, a11 = _second_order_coeffs([second_order_wn(ts, P['zeta'])], P['zeta'], DT)[:, 0]
    th, v, out = 0.0, 0.0, []
    for k in range(n):
        out.append(th)
        u = a if k > n_int else (a * (1.0 - frac) if k == n_int else 0.0)   # the plant's delay line, step at tick 0
        e = u - th
        if P['servo_order'] == 1:
            th = th + DT * np.clip(e / ts, -w, w)
        else:
            dth = a00m1 * (-e) + a01 * v
            v = np.clip(a10 * (-e) + a11 * v, -w, w)
            th = th + np.clip(dth, -w * DT, w * DT)
    return np.arange(n) * DT, np.array(out) / a
