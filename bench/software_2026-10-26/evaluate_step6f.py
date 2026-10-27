"""Judge of step 6f: the pre-registered rules of step6f_uprime_local.py applied to its frozen predictions and its
evaluation rows; also the pilot's seed-count rule. Written with step6f_uprime_local.py, before any run."""
import hashlib, json, math, os, sys
import numpy as np
from scipy.interpolate import CubicSpline

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import step6f_uprime_local as F

SCALED_JUDGED = ('C50', 'C90', 'X9', 'X10')
RNG_EV = 2026


# ------------------------------------------------------------------------------------------------ loading
def load_predictions(tag=F.TAG):
    P = json.load(open(os.path.join(F.RES, f'step6f_predictions{tag}.json')))
    Z = np.load(os.path.join(F.RES, f'step6f_pred_draws{tag}.npz')); draws = {}
    for k in Z.files:
        arm, key = k.split('__'); draws.setdefault(arm, {})[key] = Z[k]
    return P, draws


def load_rows(tag=F.TAG):
    Z = np.load(os.path.join(F.RES, f'step6f_rows{tag}.npz')); rows = {}
    for k in Z.files:
        if k == 'meta':
            continue
        arm, key = k.split('__'); rows.setdefault(arm, {})[key] = Z[k]
    return rows, json.loads(str(Z['meta']))


# ------------------------------------------------------------------------------------------------ measurement
def curves(r, arm, lev, idx=None):
    """population mean-square error curves on GRID at each level (seed means, or means over the draw idx)."""
    A, B, X = (r[k] if idx is None else r[k][idx] for k in ('A', 'B', 'X'))
    Am, Bm, Xm = A.mean(0), B.mean(0), X.mean(0); out = []
    for i, s2 in enumerate(lev):
        if F.ARMS[arm]['kind'] == 'scaled':
            c = s2 / F.Q; out.append(Am + c * Bm[0] + 2 * math.sqrt(c) * Xm[0])
        else:
            out.append(Am + Bm[i] + 2 * Xm[i])
    return out


def optima(Js):
    t, j = zip(*[F.argmin_log(J, F.GRID) for J in Js]); return np.array(t), np.array(j)


def regret(J, t):
    """sqrt(J(t) / J_min) - 1 on the population curve, both read from one cubic spline in log tau through the grid
    values (so that the regret is not negative by interpolation); nan outside the grid."""
    lg = np.log(F.GRID)
    if not (np.isfinite(t) and lg[0] <= math.log(t) <= lg[-1]):
        return float('nan')
    cs = CubicSpline(lg, J); jmin = float(cs(np.linspace(lg[0], lg[-1], 8001)).min())
    return math.sqrt(float(cs(math.log(t))) / jmin) - 1


def expo(arm, tau, lev):
    return F.secants(tau, lev) if F.ARMS[arm]['kind'] == 'scaled' else F.updown(tau, lev)


def rms(v):
    v = np.asarray(v, float); return float(np.sqrt(np.mean(v * v)))


def branch_preds(A, lev):
    """J arm: U'-Local's lower-branch and upper-branch minimisers at each level (TF below / above 1 s)."""
    lo = F.TF < 1.0; hi = ~lo; out = []
    for s2 in lev:
        J = A + s2 * F.beta_disc(F.TF)
        out.append((F.argmin_log(J[lo], F.TF[lo])[0], F.argmin_log(J[hi], F.TF[hi])[0]))
    return np.array(out)


# ------------------------------------------------------------------------------------------------ statistics per draw
def arm_stats(arm, r, P, Dr, idx=None, b=None):
    """measured exponents and the predictors' exponents; idx an evaluation draw (None: all seeds), b a calibration
    draw (None: the point prediction)."""
    lev = np.array(P[arm]['levels']); t, j = optima(curves(r, arm, lev, idx)); e = expo(arm, t, lev)
    pr = lambda k: (np.array(P[arm][k], float) if b is None else Dr[arm][k][b].astype(float))
    s = {'t': t, 'j': j, 'e': e}
    for k in ('LOCAL', 'LOCAL_G', 'AWARE', 'AWARE_G', 'RAWF', 'ENGINE', 'SQRT', 'LOCALC'):
        if k in P[arm] and (b is None or k in Dr[arm]):
            s['t_' + k] = pr(k); s['e_' + k] = expo(arm, s['t_' + k], lev)
    if F.ARMS[arm]['kind'] == 'scaled':
        s['e_RAWS'] = pr('RAWS_exp')
    return s


def boot(arm, r, P, Dr, mode='joint'):
    """B draws: 'joint' (evaluation and calibration draws together), 'ev' (evaluation only), 'cal' (calibration only)."""
    ns = len(r['A']); g = np.random.default_rng(RNG_EV + F.ARM_ID[arm]); out = []
    for b in range(F.B_BOOT):
        idx = g.integers(0, ns, ns)
        out.append(arm_stats(arm, r, P, Dr, idx=None if mode == 'cal' else idx, b=None if mode == 'ev' else b))
    return out


def ci(v):
    v = np.asarray(v, float); v = v[np.isfinite(v)]
    return (float(np.quantile(v, 0.025)), float(np.quantile(v, 0.975))) if len(v) else (float('nan'), float('nan'))


def sd(v):
    v = np.asarray(v, float); v = v[np.isfinite(v)]; return float(np.std(v, ddof=1)) if len(v) > 1 else float('nan')


def dropped(v):
    """number of bootstrap draws not finite (dropped by ci and sd), reported with each judged statistic."""
    return int(np.sum(~np.isfinite(np.asarray(v, float))))


# the statistics the hypotheses use; U'-Local's curves are read on the evaluation grid as the measured ones are (_G)
d_raws = lambda s: rms(s['e'] - s['e_RAWS']) - rms(s['e'] - s['e_LOCAL_G'])                # H6f.2
d_sqrt = lambda s: rms(s['e'] - 0.5) - rms(s['e'] - s['e_LOCAL_G'])                        # H6f.3
blind_gap = lambda s: s['e_LOCAL_G'][1] - s['e'][1]                                       # H6f.5, down exponent
aware_dev = lambda s: s['e'][1] - s['e_AWARE_G'][1]


# ------------------------------------------------------------------------------------------------ pilot rule
def n_required(se_ev, se_cal, gap, n0):
    """smallest N (evaluation seeds) with 3 sqrt(se_cal^2 + se_ev^2 n0 / N) <= gap, se_ev measured on n0 pilot seeds;
    inf if unreachable."""
    room = (gap / 3) ** 2 - se_cal ** 2
    return float('inf') if not (room > 0) else n0 * se_ev ** 2 / room


def pilot_report(rows, sha, runtime):
    P, Dr = load_predictions(); out = []; W = out.append
    W(f'step6f pilot  script sha256 {sha}  judge sha256 {F.sha_of(os.path.abspath(__file__))}')
    W(f'predictions sha256 {F.sha_of(os.path.join(F.RES, f"step6f_predictions{F.TAG}.json"))}  pilot seeds {F.PILOT_SEEDS[0]}..{F.PILOT_SEEDS[-1]}  runtime {runtime:.0f} s')
    W('Only standard errors, edge flags and checks are printed; no window, exponent, difference or verdict.')
    for arm in F.ARMS:
        s = arm_stats(arm, rows[arm], P, Dr); W(f'  {arm}: identity check {rows[arm]["ident"]:.1e} rad; any population optimum at a grid edge: {bool(np.any(~np.isfinite(s["t"])))}')
    req = {}
    for arm, stat, gapf in (('X9', d_raws, lambda p: rms(np.array(p['RAWS_exp']) - np.array(p['LOCAL_G_exp']))),
                            ('X10', d_raws, lambda p: rms(np.array(p['RAWS_exp']) - np.array(p['LOCAL_G_exp']))),
                            ('S50', blind_gap, lambda p: p['LOCAL_G_exp'][1] - p['AWARE_G_exp'][1])):
        ev = [stat(x) for x in boot(arm, rows[arm], P, Dr, 'ev')]; cal = [stat(x) for x in boot(arm, rows[arm], P, Dr, 'cal')]
        se_ev, se_cal, gap = sd(ev), sd(cal), gapf(P[arm]); n0 = len(rows[arm]['A'])
        assert F.SMOKE or (np.isfinite(gap) and np.isfinite(se_ev) and np.isfinite(se_cal)), f'{arm}: pilot statistic not finite'
        n = n_required(se_ev, se_cal, gap, n0); req[arm] = n
        W(f'  {arm}: predicted gap {gap:.4f} (from the frozen predictions); SE at {n0} pilot seeds: evaluation {se_ev:.4f},'
          f' calibration {se_cal:.4f}; seeds required {"unreachable" if not np.isfinite(n) else f"{n:.0f}"}')
    need = max([v for v in req.values() if np.isfinite(v)], default=F.N_MIN)
    N = int(min(F.N_MAX, max(F.N_MIN, 100 * math.ceil(need / 100))))
    nd = [a for a, v in req.items() if not (np.isfinite(v) and v <= F.N_MAX)]
    W(f'N_EVAL = {N}  (rule: the largest requirement, rounded up to a multiple of 100, within [{F.N_MIN}, {F.N_MAX}])')
    W('non-discriminant, reported and not judged: ' + (', '.join(f'{a} ({"H6f.2" if a != "S50" else "H6f.5"})' for a in nd) if nd else 'none'))
    return '\n'.join(out)


# ------------------------------------------------------------------------------------------------ judge
def main():
    P, Dr = load_predictions(); rows, meta = load_rows(); assert F.SMOKE or not meta['smoke']
    pp = os.path.join(F.RES, f'step6f_pilot{F.TAG}.txt'); assert F.SMOKE or os.path.exists(pp), 'the pilot file is missing'
    pilot = open(pp).read() if os.path.exists(pp) else ''
    if not F.SMOKE:
        assert f'N_EVAL = {meta["N"]} ' in pilot, 'N differs from the pilot rule'
        assert F.sha_of(os.path.join(F.RES, 'step6f_predictions.json')) in pilot, 'the predictions differ from those the pilot used'
    nondisc = pilot.split('non-discriminant, reported and not judged:')[1].strip() if 'non-discriminant' in pilot else 'none'
    out = []; W = out.append; V = {}
    W(f'evaluate_step6f.py sha256 {F.sha_of(os.path.abspath(__file__))}')
    W(f'step6f script sha256 (evaluation) {meta["sha256"]}; predictions by script sha256 {P["sha256"]}; N = {meta["N"]} seeds per arm')
    W(f'non-discriminant per the pilot: {nondisc}')
    S, BS = {}, {}
    for arm in F.ARMS:
        S[arm] = arm_stats(arm, rows[arm], P, Dr); BS[arm] = boot(arm, rows[arm], P, Dr, 'joint')
    # H6f.0
    W('\n== H6f.0 validity')
    ident = max(float(rows[a]['ident']) for a in F.ARMS); edges = {a: bool(np.any(~np.isfinite(S[a]['t']))) for a in F.ARMS}
    W(f'  identity of the filter decomposition with step 6\'s Mahony: max {ident:.1e} rad (<= 1e-9: {ident <= 1e-9})')
    W(f'  population optima at a grid edge: {edges}')
    for arm in F.ARMS:
        Bm = rows[arm]['B'].mean(0)
        if F.ARMS[arm]['kind'] == 'scaled':
            w = Bm[0] * (2 * F.GRID - F.DT) / F.Q; W(f'  {arm}: white-noise check B (2 tau - DT) / Q over the grid: {w.min():.4f} .. {w.max():.4f};'
                                                    ' the same unit path at every level, so the noise shape is identical across levels by construction')
        else:
            for i, v in enumerate(F.V_S):
                w = Bm[i] * (2 * F.GRID - F.DT) / F.S6.S0(v)
                W(f'  {arm} V = {v:g}: B (2 tau - DT) / S0(V) at tau = {F.GRID[0]:.3f}, 0.1, 1, 10 s: '
                  f'{w[0]:.3f}, {np.interp(0.1, F.GRID, w):.3f}, {np.interp(1.0, F.GRID, w):.3f}, {np.interp(10.0, F.GRID, w):.3f}')
    V['H6f.0'] = ident <= 1e-9 and not any(edges.values())
    # per arm tables
    for arm in F.ARMS:
        s = S[arm]; lev = np.array(P[arm]['levels']); W(f'\n== {arm}  levels sigma^2 / Q = {np.round(lev / F.Q, 4).tolist()}')
        Js = curves(rows[arm], arm, lev)
        W(f'  measured tau_opt (s) {np.round(s["t"], 4).tolist()}')
        for k in ('LOCAL', 'LOCAL_G', 'AWARE', 'AWARE_G', 'RAWF', 'ENGINE', 'SQRT', 'LOCALC'):
            if 't_' + k in s:
                reg = [regret(J, tp) for J, tp in zip(Js, s['t_' + k])]
                W(f'  {k:6s} tau {np.round(s["t_" + k], 4).tolist()}  ratio to tau_opt {np.round(s["t_" + k] / s["t"], 3).tolist()}'
                  f'  regret % {np.round(100 * np.array(reg), 2).tolist()}')
        if arm != 'J':
            W(f'  measured exponents {np.round(s["e"], 4).tolist()}' + ('  (up, down)' if F.ARMS[arm]['kind'] == 'shape' else ''))
            for k in ('LOCAL_G', 'LOCAL', 'AWARE_G', 'AWARE', 'RAWS', 'RAWF', 'ENGINE', 'SQRT', 'LOCALC'):
                if 'e_' + k in s:
                    dev = s['e'] - s['e_' + k]; se = [sd([x['e'][q] - x['e_' + k][q] for x in BS[arm] if 'e_' + k in x]) if k != 'LOCALC' else float('nan') for q in range(len(dev))]
                    W(f'  {k:6s} exponents {np.round(s["e_" + k], 4).tolist()}  measured - predicted {np.round(dev, 4).tolist()}  SE {np.round(se, 4).tolist()}')
    W('\n== Hypotheses')
    ok = []; small = []
    for arm in SCALED_JUDGED:
        s = S[arm]; dev = s['e'] - s['e_LOCAL_G']; Dq = [[x['e'][q] - x['e_LOCAL_G'][q] for x in BS[arm]] for q in range(len(dev))]
        se = np.array([sd(v) for v in Dq]); good = bool(np.all(np.abs(dev) <= 3 * se)); small.append(bool(np.all(np.abs(dev) < 0.005)))
        ok.append(good); W(f'H6f.1 {arm}: max |e - e_LOCAL_G| = {np.max(np.abs(dev)):.4f}; max |dev| / SE = {np.max(np.abs(dev) / se):.2f};'
                           f' all within 3 SE: {good}; draws dropped {max(dropped(v) for v in Dq)}')
    V['H6f.1'] = all(ok); ok = []
    W(f'   reading rule: every deviation under 0.005 in all four arms: {all(small)}')
    for arm in ('X9', 'X10'):
        bd = [d_raws(x) for x in BS[arm]]; d = d_raws(S[arm]); c = ci(bd); good = d > 0 and c[0] > 0
        tag = ' (non-discriminant: reported, not judged)' if arm in nondisc else ''
        W(f'H6f.2 {arm}: RMS(e - e_RAWS) - RMS(e - e_LOCAL_G) = {d:.4f} [{c[0]:.4f}, {c[1]:.4f}]; positive with interval above 0: {good};'
          f' draws dropped {dropped(bd)}{tag}')
        if arm not in nondisc:
            ok.append(good)
    V['H6f.2'] = all(ok) if ok else None; ok = []
    for arm in ('C90', 'X9', 'X10'):
        bd = [d_sqrt(x) for x in BS[arm]]; d = d_sqrt(S[arm]); c = ci(bd); good = d > 0 and c[0] > 0; ok.append(good)
        W(f'H6f.3 {arm}: RMS(e - 1/2) - RMS(e - e_LOCAL_G) = {d:.4f} [{c[0]:.4f}, {c[1]:.4f}]; interval above 0: {good}; draws dropped {dropped(bd)}')
    V['H6f.3'] = all(ok); ok = []
    for arm in SCALED_JUDGED:
        s = S[arm]; lev = np.array(P[arm]['levels']); Js = curves(rows[arm], arm, lev)
        ratio = s['t_LOCAL'] / s['t']; reg = np.array([regret(J, tp) for J, tp in zip(Js, s['t_LOCAL'])])
        good = bool(np.all((ratio >= 1 / 1.1) & (ratio <= 1.1)) and np.all(reg <= 0.01)); ok.append(good)
        W(f'H6f.4 {arm}: LOCAL window / tau_opt {ratio.min():.3f} .. {ratio.max():.3f}; regret max {100 * reg.max():.2f} %; within x1.1 and 1 %: {good}')
    V['H6f.4'] = all(ok)
    s = S['S50']; ad = [aware_dev(x) for x in BS['S50']]; bg = [blind_gap(x) for x in BS['S50']]
    a_dev = aware_dev(s); a_se = sd(ad); b_gap = blind_gap(s); b_se = sd(bg)
    good = abs(a_dev) <= 3 * a_se and b_gap > 3 * b_se
    tag = ' (non-discriminant: reported, not judged)' if 'S50' in nondisc else ''
    W(f'H6f.5 S50: down exponent {s["e"][1]:.4f}; minus AWARE {a_dev:.4f} (SE {a_se:.4f}); BLIND minus measured {b_gap:.4f} (SE {b_se:.4f});'
      f' AWARE within 3 SE and BLIND beyond 3 SE: {good}; draws dropped {max(dropped(ad), dropped(bg))}{tag}')
    V['H6f.5'] = None if 'S50' in nondisc else good
    s = S['S90']
    W(f'   reported: S90 down exponent {s["e"][1]:.4f}; AWARE {s["e_AWARE_G"][1]:.4f}; BLIND {s["e_LOCAL_G"][1]:.4f}')
    # H6f.6
    lev = np.array(P['J']['levels']); t = S['J']['t']; r = t[1:] / t[:-1]; bp = branch_preds(np.array(P['J']['A_pooled']), lev)
    big = np.where(r >= 2.0)[0]; one = len(big) == 1 and int(np.sum(r >= 1.5)) == 1
    loc = one and F.J_REL[big[0]] >= 0.92 and F.J_REL[big[0] + 1] <= 1.08
    on_branch = np.where(t > 1.0, bp[:, 1], bp[:, 0]); fit = t / on_branch; fit_ok = bool(np.all((fit >= 1 / 1.1) & (fit <= 1.1)))
    good = bool(one and loc and fit_ok)
    W(f'H6f.6 J: measured tau_opt {np.round(t, 3).tolist()}')
    W(f'   adjacent ratios {np.round(r, 3).tolist()}; one ratio >= 2 and no other >= 1.5: {one}; between relative levels 0.92 and 1.08: {loc};'
      f' each window within x1.1 of U\'-Local\'s branch: {fit_ok} ({fit.min():.3f} .. {fit.max():.3f})')
    for k in ('RAWF', 'SQRT', 'ENGINE'):
        tk = np.array(P['J'][k], float); rr = tk[1:] / tk[:-1]
        W(f'   reported: {k} largest adjacent ratio ' + (f'{np.nanmax(rr):.3f}' if np.any(np.isfinite(rr)) else 'undefined (no window)'))
    V['H6f.6'] = good
    W('\n== Verdict')
    if not V['H6f.0']:
        W('H6f.0 FAILS: the run is void; nothing else is judged (the lines above are reported only).')
        for k in V:
            W(f'{k}: ' + ('FAILS' if k == 'H6f.0' else 'void'))
    else:
        for k, v in V.items():
            W(f'{k}: {"not judged (non-discriminant)" if v is None else ("HOLDS" if v else "FAILS")}')
    txt = '\n'.join(out); print(txt); open(os.path.join(F.RES, f'step6f_results{F.TAG}.txt'), 'w').write(txt + '\n')


if __name__ == '__main__':
    main()
