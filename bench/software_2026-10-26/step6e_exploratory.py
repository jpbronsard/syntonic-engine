"""EXPLORATORY, written after the step-6e verdict (not pre-registered). Three readings of step 6e:
(a) the pairing: the slope of H6e.1 with the seed draw shared across the four arms (as judged) against independent draws
    per arm; the evaluation keeps the pairing, and re-estimates every optimum on every draw.
(b) the level: the population optimum against U''s window with the true H and lambda_H, per arm and segment.
(c) the noise shape: step 6e's accelerometer noise is not a white noise scaled by V (a 92 V^(1/4) Hz line of amplitude
    proportional to sqrt(V), a constant sensor floor, an 80 Hz Butterworth). The optimal windows predicted from the
    exact error of the discrete filter with that noise, against those with a white noise of intensity S0(V); the up and
    down exponents of each. Step 6f (pre-registered after this reading) tests (c) at constant V.
Reads results/step6e_rows.json; writes results/step6e_exploratory.txt."""
import json, math, os
import numpy as np
from scipy.optimize import minimize_scalar
from scipy.signal import butter, freqz

import step6_quad_attitude_estimator as S6
import step6d_uprime_mechanism as S6D

HERE = os.path.dirname(os.path.abspath(__file__))
D = json.load(open(os.path.join(HERE, 'results', 'step6e_rows.json')))
DT = S6.DT; S0 = lambda V: np.asarray(V, float) * S6.Q_G + S6.S_SENS
V = D['V']; LG = np.log(np.array(D['grid'])); ARMS = ['U50', 'U60', 'U75', 'U90']; B = 2000
out = []; P_ = lambda *a: out.append(' '.join(str(x) for x in a))


def popt(J):
    i = int(np.argmin(J))
    if i in (0, len(J) - 1):
        return float(np.exp(LG[i]))
    c = np.polyfit(LG[i - 1:i + 2], J[i - 1:i + 2], 2)
    return float(np.exp(-c[1] / (2 * c[0]))) if c[0] > 0 else float(np.exp(LG[i]))


O = {a: np.array(D[a]['orc']) ** 2 for a in ARMS}; pH = np.array([D[a]['p_H'] for a in ARMS]); ns = O['U50'].shape[2]


def ups(idx):
    e = []
    for a, ix in zip(ARMS, idx):
        t = [popt(O[a][g][:, ix].mean(1)) for g in (0, 1)]; e.append(math.log(t[1] / t[0]) / math.log(S0(V[1]) / S0(V[0])))
    return np.polyfit(pH, e, 1)[0]


rng = np.random.default_rng(2026)
paired = [ups([i] * 4) for i in (rng.integers(0, ns, ns) for _ in range(B))]
rng = np.random.default_rng(2027)
unpaired = [ups([rng.integers(0, ns, ns) for _ in ARMS]) for _ in range(B)]
P_('(a) slope of the upward exponents on p_H, 95 % percentile intervals over 2000 draws')
P_(f'    paired (one seed draw for the four arms): [{np.quantile(paired, 0.025):.3f}, {np.quantile(paired, 0.975):.3f}]')
P_(f'    unpaired (independent draws per arm):     [{np.quantile(unpaired, 0.025):.3f}, {np.quantile(unpaired, 0.975):.3f}]')

P_("\n(b) population optimum / U' window with the true H and lambda_H (V = 1, 256, 1/16)")
for a in ARMS:
    H = D[a]['H']; lam = D[a]['lamH']
    pred = [S6.kappa_H(H) * (S0(v) / lam) ** (1 / (2 * H + 1)) for v in V]
    t = [popt(O[a][g].mean(1)) for g in range(3)]
    P_(f'    {a}: {np.round(np.array(t) / np.array(pred), 3).tolist()}')

P_('\n(c) predicted exponents, white noise of intensity S0(V) against step 6e\'s noise (line + floor + Butterworth)')
bb, ab = butter(2, 80.0 / (0.5 / DT)); th = np.linspace(0, np.pi, 200001); _, Bw = freqz(bb, ab, worN=th); B2 = np.abs(Bw) ** 2
K2 = lambda tau, w: (DT / tau) ** 2 / np.abs(1 - (1 - DT / tau) * np.exp(-1j * w)) ** 2


def noise(tau, v, shape):
    if not shape:
        return S0(v) / (2 * tau - DT)
    white = S0(v) / DT * np.trapezoid(B2 * K2(tau, th), th) / np.pi
    w = 2 * np.pi * S6.LINE_F * v ** 0.25 * DT; _, b = freqz(bb, ab, worN=[w])
    return white + (S6.LINE_AMP ** 2 * v / 2) * abs(b[0]) ** 2 * K2(tau, np.array([w]))[0]


def track(tau, H, lam, M=400000):
    a = DT / tau; j = np.arange(1, M); return a * (1 - a) / (2 - a) * np.sum((1 - a) ** j * lam * (j * DT) ** (2 * H))


meas = json.load(open(os.path.join(HERE, 'results', 'step6e_rows.json')))
for a in ARMS:
    H = D[a]['H']; lam = D[a]['lamH']; row = []
    for shape in (False, True):
        t = [math.exp(minimize_scalar(lambda x: math.log(track(math.exp(x), H, lam) + noise(math.exp(x), v, shape)),
                                      bounds=(math.log(0.02), math.log(60)), method='bounded', options={'xatol': 1e-7}).x) for v in V]
        row += [math.log(t[1] / t[0]) / math.log(S0(V[1]) / S0(V[0])), math.log(t[2] / t[0]) / math.log(S0(V[2]) / S0(V[0]))]
    P_(f'    {a}: white up {row[0]:.3f} down {row[1]:.3f}; 6e noise up {row[2]:.3f} down {row[3]:.3f}')
txt = '\n'.join(out); print(txt)
open(os.path.join(HERE, 'results', 'step6e_exploratory.txt'), 'w').write(txt + '\n')
