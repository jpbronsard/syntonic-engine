"""Shared harness, goal/flight exploration (2 Nov 2026, exploratory, not pre-registered).

Imports the flight bench read-only from the repository (the sys.argv trick of explore_kfff.py), installs amendment 2's
noise model (gyro records keep their logged bias, tilt records mean-removed), builds the six circumstances by a copy of
build_noise with a schedule and modifiers, and runs vectorised estimator columns with the logic of est_only() (plant at
rest: the estimate is the estimation error; arming calibration = pre-roll mean of the gyro at idle).
"""
import os, sys, json, math
import numpy as np

sys.dont_write_bytecode = True                       # never write into the repository
REPO = '/home/claude/repo/syntonic-engine'
BENCH = os.path.join(REPO, 'bench', 'flight_2026-11-02')
WORK = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(WORK, 'results'); os.makedirs(OUT, exist_ok=True)
sys.path.insert(0, BENCH)
_argv = list(sys.argv); sys.argv = ['x', REPO, 'import']
import flight_bench as F                              # noqa: E402
sys.argv = _argv
G = F.G
D2R, DT, DTA, NSUB, T_PRE, T_FLY = F.D2R, F.DT, F.DTA, F.NSUB, F.T_PRE, F.T_FLY
KPRE = int(round(T_PRE / DT))


def records_a2():
    if not F._REC:
        for name, seg in (('idle', G.SEGS[0]), ('low', G.SEGS[1]), ('half', G.SEGS[2])):
            tilt, _, g1k = G.load(*seg[:4]); F._REC[name] = (tilt - tilt.mean(), g1k)   # gyro keeps its logged bias
    return F._REC
F.records = records_a2
FX = json.load(open(os.path.join(BENCH, 'results_a2', 'facts.json')))
TUNED_A2 = json.load(open(os.path.join(BENCH, 'results_a2', 'tuned_estimators.json')))

TRAIN, TEST, DEV = [100, 101, 102], list(range(200, 210)), list(range(300, 310))
FORBIDDEN = set(range(400, 500))

# ---------------------------------------------------------------------------------------------------- circumstances
SQ10 = math.sqrt(10.0)
CIRC = {
    'C1': dict(segs=[('low', 0, 60)], desc='steady low throttle'),
    'C2': dict(segs=[('half', 0, 60)], desc='steady half throttle'),
    'C3': dict(segs=list(F.REGS), desc='flight bench schedule, half/low every 15 s'),
    'C4': dict(segs=[('low', 0, 60)], tilt_amp=(20.0, 40.0, SQ10), desc='vibration step: tilt noise x sqrt(10), 20-40 s, low'),
    'C5': dict(segs=[(('half', 'low')[i % 2], 5 * i, 5 * i + 5) for i in range(12)], desc='fast switching half/low every 5 s'),
    'C6': dict(segs=[('low', 0, 60)], gyro_step=(30.0, 0.02), desc='bias step +0.02 deg/s on the gyro from 30 s, low'),
}


def build_noise_c(seed, circ):
    """copy of F.build_noise with a schedule and modifiers. Returns gyro noise per 1 ms tick (rad/s), tilt noise per
    4 ms sample (rad), throttle label per tick. For C3 it is identical to F.build_noise (checked in check_harness)."""
    assert seed not in FORBIDDEN or os.environ.get('GOAL_ALLOW_CONFIRM') == '1', \
        'seeds 400-499 are reserved for the confirmation (set GOAL_ALLOW_CONFIRM=1 only for the pre-registered run)'
    spec = CIRC[circ]; R = F.records(); rng = np.random.default_rng(seed)
    segs = [('idle', 0.0, T_PRE)] + [(n, T_PRE + a, T_PRE + b) for n, a, b in spec['segs']]
    gn, tn, reg = [], [], []
    for name, a, b in segs:
        tilt, g = R[name]; ng = int(round((b - a) / DT)); na = ng // NSUB
        o = int(rng.integers(0, len(tilt))); idx_a = (o + np.arange(na)) % len(tilt)
        idx_g = (o * NSUB + np.arange(ng)) % len(g)
        gn.append(g[idx_g]); tn.append(tilt[idx_a]); reg += [name] * ng
    gn = np.concatenate(gn); tn = np.concatenate(tn).copy(); reg = np.array(reg)
    t_g = np.arange(len(gn)) * DT - T_PRE; t_a = np.arange(len(tn)) * DTA - T_PRE
    if 'tilt_amp' in spec:
        a, b, f = spec['tilt_amp']; tn[(t_a >= a) & (t_a < b)] *= f
    if 'gyro_step' in spec:
        a, db = spec['gyro_step']; gn = gn + np.where(t_g >= a, db * D2R, 0.0)
    return gn, tn, reg


# ---------------------------------------------------------------------------------------------------- harness
def run(seed, circ, ests, fx=FX, phases=None, trace_every=0):
    """ests: list of estimator objects with .n columns and .step(g, a, reg_lam, t_since, reg_name, k) -> th (rad),
    optional .arm(bias_cal). Same logic as F.est_only. Returns list of J_est arrays (deg), one per object; with phases
    (list of (t0, t1) in flight time) also the per-phase RMS; with trace_every the estimates every that many ticks."""
    gn, tn, reg = build_noise_c(seed, circ); N = len(gn)
    lam_of = {r: fx['lam'][r] for r in ('idle', 'low', 'half')}
    gsum = 0.0; gcnt = 0; bias_cal = 0.0
    ss = [np.zeros(e.n) for e in ests]; m = 0
    ph = phases or []; ps = [[np.zeros(e.n) for e in ests] for _ in ph]; pm = [0] * len(ph)
    tr = [[] for _ in ests]
    for k in range(N):
        g0 = gn[k] - bias_cal; a0 = tn[k // NSUB] if k % NSUB == 0 else None
        rl = lam_of[reg[k]]; tk = k * DT
        outs = [e.step(g0, a0, rl, tk, reg[k], k) for e in ests]
        if k < KPRE:
            gsum += gn[k]; gcnt += 1
            if k == KPRE - 1:
                bias_cal = gsum / gcnt
                for e in ests:
                    if hasattr(e, 'arm'): e.arm(bias_cal)
            continue
        tf = tk - T_PRE
        for i, th in enumerate(outs): ss[i] += (th / D2R) ** 2
        m += 1
        for j, (p0, p1) in enumerate(ph):
            if p0 <= tf < p1:
                for i, th in enumerate(outs): ps[j][i] += (th / D2R) ** 2
                pm[j] += 1
        if trace_every and (k - KPRE) % trace_every == 0:
            for i, th in enumerate(outs): tr[i].append(np.array(th, float) / D2R)
    J = [np.sqrt(s / m) for s in ss]
    res = dict(J=J)
    if ph: res['phase'] = [[np.sqrt(ps[j][i] / max(pm[j], 1)) for j in range(len(ph))] for i in range(len(ests))]
    if trace_every: res['trace'] = [np.array(t) for t in tr]
    return res


class BenchEst:
    """wraps F.Est (kinds 1 Kalman, 2 law) so it plugs into run(); per-column Q_B allowed (F.Q_B is read as a module
    global inside F.Est.step, so it is swapped in and out around each step)."""
    def __init__(self, kind, cfac=None, kR=None, kQ=None, qB=None):
        kind = np.asarray(kind); self.n = len(kind)
        self.e = F.Est(self.n, kind, kf_R=kR, kf_QG=kQ, cfac=cfac)
        self.qB = F.Q_B if qB is None else np.asarray(qB, float)
    def step(self, g, a, reg_lam, t_since, reg_name, k):
        q0 = F.Q_B; F.Q_B = self.qB
        try:
            gg = np.full(self.n, g); aa = None if a is None else np.full(self.n, a)
            return self.e.step(gg, aa, reg_lam, t_since)
        finally:
            F.Q_B = q0


def pooled(J):
    """J: array (seeds, cols) of per-seed RMS -> pooled RMS over seeds"""
    return np.sqrt(np.mean(np.square(np.asarray(J)), 0))


def dump(name, obj):
    p = os.path.join(OUT, name)
    json.dump(obj, open(p, 'w'), indent=1, default=lambda o: o.tolist() if hasattr(o, 'tolist') else str(o))
    return p
