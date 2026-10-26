"""Re-run of the ENGINE rows of the software bench with the exact identification of syntonic-engine 0.2.0:
block means with the fit restricted to the six finest scales that clear the floor (the bench scripts fitted all
of them). Banks, oracles and fixed windows are taken from the stored results of each step (unchanged). The
scripts themselves are not modified; this file records the engine rows for bench note v1.1.
"""
import json, os, sys, time, hashlib
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import step0_uprime_engine as S0
import step0b_allan_identification as S0b
import step2_regime_jumps as S2
import step3c_meta_noise as S3c
import step4_drift as S4
FINEST = 6


class Allan6(S0b.AllanSmoother):
    """step0b's smoother with the 0.2.0 restriction: the fit uses the finest FINEST scales that clear the floor."""

    def step(self, a, b):
        dt = self.dt
        meas = 0.5 * (a + b); d = (a - b) / np.sqrt(2.0)
        self.n2 += self.a_w * (d * d - self.n2)
        return self._rest(meas)

    def _rest(self, meas):
        dt = self.dt
        self.acc += meas[None, :]; self.cnt += 1
        done = self.cnt == self.scales
        if done.any():
            for i in np.flatnonzero(done):
                mean = self.acc[i] / self.scales[i]
                if not np.isnan(self.prev[i, 0]):
                    dd = mean - self.prev[i]
                    aw = self.a_blk[i] if self.nblk[i] > 0 else 1.0
                    self.DA[i] += aw * (dd * dd - self.DA[i]); self.nblk[i] += 1
                self.prev[i] = mean; self.acc[i] = 0.0; self.cnt[i] = 0
        if self.k > self.warm:
            floor = self.n2[None, :] * dt / (self.scales * dt)[:, None]
            D = self.DA - floor
            have = (self.nblk >= 2)[:, None] & np.ones((1, self.n), bool)
            ok_g = have & (D > self.snr_gate * floor); ok_u = have & (D > 0)
            use_g = ok_g.sum(0) >= 2
            ok = np.where(use_g[None], ok_g, ok_u); ok &= np.cumsum(ok, axis=0) <= FINEST; nok = ok.sum(0)
            logD = np.where(ok, np.log(np.maximum(D, 1e-30)), 0.0)
            W = ok.astype(float); L = self.logl[:, None]
            Sw = W.sum(0); Sx = (W * L).sum(0); Sy = (W * logD).sum(0)
            Sxx = (W * L ** 2).sum(0); Sxy = (W * L * logD).sum(0)
            den = Sw * Sxx - Sx ** 2
            slope = np.where(den > 1e-12, (Sw * Sxy - Sx * Sy) / np.maximum(den, 1e-12), 1.0)
            Hn = np.clip(0.5 * slope, 0.05, 1.0); valid = nok >= 2
            s2 = 0.5 * self.n2 * dt
            lam = np.exp((W * (logD - 2 * Hn[None] * L)).sum(0) / np.maximum(Sw, 1)) / S0b.C_H(Hn)
            ts = S0.kappa_H(Hn) * (np.maximum(s2, 1e-30) / np.maximum(lam, 1e-30)) ** (1.0 / (2 * Hn + 1))
            ts = np.clip(ts, self.floor, self.tau_max)
            self.tau_target = np.where(valid, ts, self.tau_target); self.H = np.where(valid, Hn, self.H)
            rr = 1.0 / self.tau; dr_max = self.tau_rate * dt / self.tau ** 2
            rr = rr + np.clip(1.0 / self.tau_target - rr, -dr_max, dr_max); self.tau = 1.0 / rr
        g = 1.0 - np.exp(-dt / self.tau)
        self.xhat = meas.copy() if self.xhat is None else self.xhat + g * (meas - self.xhat)
        self.k += 1
        return self.xhat


class FixedNoise6(Allan6):
    def __init__(self, n, dt, floor, tau_max, tau_rate, var_known):
        super().__init__(n, dt, floor, tau_max, tau_rate); self.n2_fixed = 2.0 * var_known

    def step(self, y):
        self.n2[:] = self.n2_fixed; return self._rest(y)


class Meta6(Allan6):
    def __init__(self, n, dt, floor, tau_max, tau_rate):
        super().__init__(n, dt, floor, tau_max, tau_rate)
        self.meta = FixedNoise6(n, dt, floor, tau_max, tau_rate, S3c.LOGCHI2_VAR); self.n2_ema = np.zeros(n)

    def step(self, a, b):
        d = (a - b) / np.sqrt(2.0); d2 = np.maximum(d * d, 1e-300)
        self.n2_ema += self.a_w * (d2 - self.n2_ema)
        ym = self.meta.step(np.log(d2)); n2_meta = np.exp(ym + (-S3c.LOGCHI2_MEAN))
        use_meta = self.meta.k > self.meta.warm
        meas = 0.5 * (a + b)
        self.n2 = np.where(use_meta, n2_meta, self.n2_ema)
        return self._rest(meas)


def eng_args():
    return (S0.DT, S0.ENGINE['floor'], S0.ENGINE['tau_max'], S0.ENGINE['tau_rate'])


def run_engine(cls, A, B, X, keep=None):
    n = X.shape[1]; eng = cls(n, *eng_args()); E = np.empty_like(X); TAU = np.empty_like(X); HE = np.empty_like(X)
    for k in range(X.shape[0]):
        xh = eng.step(A[k], B[k]); E[k] = (xh - X[k]) ** 2; TAU[k] = eng.tau; HE[k] = eng.H
    return E, TAU, HE


def step0b_and_1():
    print('\n=== step 0b / step 1 engine rows (0.2.0 identification) ===')
    rows1 = json.load(open(os.path.join(HERE, 'step1', 'step1_results.json')))['rows']
    DT, N = S0.DT, S0.N; t = np.arange(N) * DT; keep = t > S0.T_SKIP
    print('H     level: tau_eng/tau_opt  regret | step1 C3 ratio ENGINE/minKF (new) [old]')
    worst = dict(ratio=(1.0, 1.0), regret=0.0, c3=0.0)
    for H in S0.H_LIST:
        kH = float(S0.kappa_H(H)); s2l = S0.LAMBDA_H * (S0.TAU_TARGETS / kH) ** (2 * H + 1); sig = np.sqrt(2 * s2l / DT)
        X, A, B, cells = [], [], [], []
        for li, sg in enumerate(sig):
            for seed in S0.SEEDS:
                rp = np.random.default_rng(10_000 * seed + int(round(H * 100))); rn = np.random.default_rng(10_000 * seed + 100 * (li + 1) + int(round(H * 100)))
                x = np.sqrt(S0.LAMBDA_H) * S0.fbm(N, H, rp, DT); X.append(x); A.append(x + sg * rn.standard_normal(N)); B.append(x + sg * rn.standard_normal(N)); cells.append((li, seed))
        X = np.array(X).T; A = np.array(A).T; B = np.array(B).T
        E, TAU, HE = run_engine(Allan6, A, B, X)
        mse = E[keep].mean(0); te = np.median(TAU[keep], 0); he = np.median(HE[keep], 0)
        line = []
        for li in range(7):
            idx = [c for c, (l, s) in enumerate(cells) if l == li]
            r1 = [r for r in rows1 if r['H'] == H and r['level'] == li]
            ratios, regs, c3n, c3o = [], [], [], []
            for c in idx:
                li_, seed = cells[c]
                row = next(r for r in r1 if r['seed'] == seed)
                regs.append(mse[c] / row['mse_opt'] - 1); ratios.append(te[c] / row['tau_opt'])
                minkf = min(row['mse_KF_RW'], row['mse_KF_CV']); c3n.append(mse[c] / minkf); c3o.append(row['mse_ENGINE'] / minkf)
            q = float(np.exp(np.mean(np.log(ratios)))); rg = float(np.mean(regs))
            cn = float(np.exp(np.mean(np.log(c3n)))); co = float(np.exp(np.mean(np.log(c3o))))
            worst['regret'] = max(worst['regret'], rg); worst['ratio'] = (min(worst['ratio'][0], q), max(worst['ratio'][1], q)); worst['c3'] = max(worst['c3'], cn)
            line.append(f'{q:4.2f}/{rg:+.3f} | {cn:5.3f} [{co:5.3f}]')
        print(f'H {H:.2f}: ' + '   '.join(line))
    print('worst over 35 cells: tau ratio range %.2f..%.2f, regret %.3f, ENGINE/minKF %.3f' % (worst['ratio'][0], worst['ratio'][1], worst['regret'], worst['c3']))


def jump_measures(E, c, orc, steady_windows):
    steady = [float(np.mean(E[S2.seg_mask(*S2.STEADY[i]), c])) for i in range(3)]
    out = dict(steady_regret=[steady[i] / orc[c, i] - 1 for i in range(3)])
    for j, tj in enumerate((S2.T1, S2.T2)):
        tr5, rc = S2.transient_measures(E[:, c], tj, steady[j + 1])
        out[f'excess{j+1}'] = tr5 / orc[c, j + 1] - 1; out[f'adapt{j+1}'] = tr5 / steady[j + 1] - 1; out[f'recovery{j+1}'] = rc
    return out


def steps_2_3(seed0, tag, bank_json, bank_names):
    """Engine (and META for step 3c) rows against the stored bank rows of step 2/replication/3; 14-cell rule."""
    print(f'\n=== {tag}: engine rows (0.2.0) against stored banks {bank_names} ===')
    rows_b = json.load(open(bank_json))['rows']
    seeds = list(range(seed0, seed0 + 10))
    X, A, B, meta = [], [], [], []
    for kind, H in S2.CELLS:
        for seed in seeds:
            x, a, b, _ = S2.make_stream(kind, H, seed); X.append(x); A.append(a); B.append(b); meta.append((kind, H, seed))
    X = np.array(X).T; A = np.array(A).T; B = np.array(B).T; n = X.shape[1]
    MEAS = 0.5 * (A + B); orc = np.zeros((n, 3))
    for c in range(n):
        for i in range(3): orc[c, i] = S2.oracle_segment(MEAS[:, c], X[:, c], i)[1]
    E, _, _ = run_engine(Allan6, A, B, X)
    EM = run_engine(Meta6, A, B, X)[0] if 'META' in tag else None
    def M(kind, H, name, key, j=None, src=None):
        if src is None: rs = [r for r in rows_b if r['kind'] == kind and r['H'] == H and r['method'] == name]
        else: rs = [r for r in src if r['kind'] == kind and r['H'] == H]
        return float(np.mean([r[key][j] for r in rs])) if key == 'steady_regret' else float(np.mean([r[key] for r in rs]))
    new = []; newm = []
    for c, (kind, H, seed) in enumerate(meta):
        new.append(dict(kind=kind, H=H, seed=seed, **jump_measures(E, c, orc, None)))
        if EM is not None: newm.append(dict(kind=kind, H=H, seed=seed, **jump_measures(EM, c, orc, None)))
    for label, src in (('ENGINE', new), ('ENGINE_META', newm)):
        if not src: continue
        cells = []
        print(f'  {label}: cell results (adapt_excess, recovery) vs banks; T1 steady regret')
        t1 = True
        for kind, H in S2.CELLS:
            sr = [M(kind, H, None, 'steady_regret', i, src) for i in range(3)]
            if kind != 'J4': t1 &= all(s <= 0.10 for s in sr)
            old = [M(kind, H, 'ENGINE', 'steady_regret', i) for i in range(3)]
            for j in (1, 2):
                aE, rE = M(kind, H, None, f'adapt{j}', src=src), M(kind, H, None, f'recovery{j}', src=src); res = {}
                for bank in bank_names:
                    aB, rB = M(kind, H, bank, f'adapt{j}'), M(kind, H, bank, f'recovery{j}')
                    res[bank] = label if (aE < 0.9 * aB and rE <= rB) else (bank if (aB < 0.9 * aE and rB <= rE) else 'tie')
                cells.append((kind, H, j, res, round(aE, 3), round(rE, 2)))
            print(f'    {kind} {H}: steady {sr[0]:+.3f}/{sr[1]:+.3f}/{sr[2]:+.3f} (old engine {old[0]:+.3f}/{old[1]:+.3f}/{old[2]:+.3f})')
        nE = {b: sum(1 for c in cells if c[3][b] == label) for b in bank_names}
        nB = {b: sum(1 for c in cells if c[3][b] == b) for b in bank_names}
        print(f'  {label}: wins {nE}, bank wins {nB}, T1 {"PASS" if t1 else "FAIL"}')
        for c in cells: print('     ', c)


def step4():
    print('\n=== step 4 drift: engine rows (0.2.0) ===')
    rows = json.load(open(os.path.join(HERE, 'step4', 'step4_results.json')))['rows']
    X, A, B, meta = [], [], [], []
    for Te, H in S4.CELLS:
        for seed in S4.SEEDS:
            x, a, b = S4.make_stream(Te, H, seed); X.append(x); A.append(a); B.append(b); meta.append((Te, H, seed))
    X = np.array(X).T; A = np.array(A).T; B = np.array(B).T; n = X.shape[1]
    MEAS = 0.5 * (A + B); orc = np.zeros(n); scored = None
    for c in range(n): orc[c], _, scored = S4.local_oracle_and_best(MEAS[:, c], X[:, c])
    E, _, _ = run_engine(Allan6, A, B, X); EM = run_engine(Meta6, A, B, X)[0]
    tot = E[scored].sum(0); totm = EM[scored].sum(0)
    print('  T_e   H     ENGINE new [old]   META new [old]   HARD_1s  HARD_4s')
    for Te, H in S4.CELLS:
        idx = [c for c, (te, h, s) in enumerate(meta) if te == Te and h == H]
        rs = [r for r in rows if r['Te'] == Te and r['H'] == H]
        rn = np.mean([tot[c] / orc[c] - 1 for c in idx]); rm = np.mean([totm[c] / orc[c] - 1 for c in idx])
        print(f'  {Te:4.0f}  {H:.2f}   {rn:+.3f} [{np.mean([r["regret_ENGINE"] for r in rs]):+.3f}]     {rm:+.3f} [{np.mean([r["regret_ENGINE_META"] for r in rs]):+.3f}]    '
              f'{np.mean([r["regret_HARD_1s"] for r in rs]):+.3f}   {np.mean([r["regret_HARD_4s"] for r in rs]):+.3f}')


if __name__ == '__main__':
    print('rerun_engine_v020.py sha256', hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest()[:12], flush=True)
    t0 = time.time()
    step0b_and_1(); print(f'[{time.time() - t0:.0f} s]', flush=True)
    steps_2_3(100, 'step 2 (seeds 100-109)', os.path.join(HERE, 'step2', 'step2_results.json'), ['IMM_1s', 'IMM_4s']); print(f'[{time.time() - t0:.0f} s]', flush=True)
    steps_2_3(200, 'step 2 replication (seeds 200-209)', os.path.join(HERE, 'step2_rep', 'step2_results.json'), ['IMM_1s', 'IMM_4s']); print(f'[{time.time() - t0:.0f} s]', flush=True)
    steps_2_3(100, 'step 3 / 3c hardened banks, with META', os.path.join(HERE, 'step3', 'step3_results.json'), ['IMM_1s', 'IMM_4s']); print(f'[{time.time() - t0:.0f} s]', flush=True)
    step4(); print(f'[{time.time() - t0:.0f} s]', flush=True)
