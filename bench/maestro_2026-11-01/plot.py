import os, sys
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.ndimage import uniform_filter1d

d = np.load(sys.argv[1])          # written by accel_tau.py
T, W, WF, qdd = d['T'], d['win'] * 1e3, d['win_frozen'] * 1e3, d['qdd']
floor = float(d['floor']) * 1e3
J = 1                                    # RJ2, the shoulder shown in the video
acc = np.abs(qdd[:, J]); acc_s = uniform_filter1d(acc, 100)   # 1 s average, the engine's averaging time
m = T >= 1.0                             # after the engine's 1 s warm-up

SURF, INK, INK2, MUTED, GRID = '#1a1a19', '#f0efec', '#c3c2b7', '#8f8e87', '#2c2c2a'
BLUE, ORANGE = '#3987e5', '#d95926'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 15, 'text.color': INK, 'axes.labelcolor': INK2,
                     'xtick.color': MUTED, 'ytick.color': MUTED, 'axes.edgecolor': GRID})

fig, (a1, a2) = plt.subplots(2, 1, figsize=(12, 7.5), dpi=160, sharex=True,
                             gridspec_kw=dict(hspace=0.32, top=0.80, bottom=0.145, left=0.085, right=0.80))
fig.patch.set_facecolor(SURF)
for a in (a1, a2):
    a.set_facecolor(SURF); a.grid(axis='y', color=GRID, lw=1); a.set_axisbelow(True)
    for s in ('top', 'right'): a.spines[s].set_visible(False)
    a.spines['left'].set_visible(False); a.tick_params(length=0)

fig.text(0.085, 0.935, 'When the gesture speeds up, the window tightens', fontsize=24, weight='bold', color=INK)
fig.text(0.085, 0.885, 'MAESTRO 1 conducting Beethoven\'s Fifth (simulation), right shoulder RJ2, excerpt A, steady sensor noise',
         fontsize=13.5, color=INK2)

# 1. the gesture
a1.plot(T[m], acc[m], color=BLUE, lw=1, alpha=0.30)
a1.plot(T[m], acc_s[m], color=BLUE, lw=2.4)
a1.set_ylabel('acceleration (rad/s²)')
a1.set_ylim(0, np.percentile(acc[m], 99.5) * 1.05)
a1.text(1.01, 0.62, 'Shoulder\nacceleration', transform=a1.transAxes, color=INK, fontsize=15, weight='bold', va='center')
a1.text(1.01, 0.32, '1 s average (bold)\ninstant (faint)', transform=a1.transAxes, color=MUTED, fontsize=12, va='center')

# 2. the window
a2.plot(T[m], WF[m, J], color=MUTED, lw=2, ls=(0, (6, 4)))
a2.axhline(floor, color='#4a4a46', lw=1.5)
a2.plot(T[m], W[m, J], color=ORANGE, lw=2.4)
a2.set_ylabel('window τ* (ms)')
a2.set_ylim(0, max(W[m, J].max(), 1) * 1.08)
a2.set_xlabel('time in the excerpt (s)')
a2.set_xlim(0, 30)
a2.text(1.01, 0.70, 'Window set\nby the formula', transform=a2.transAxes, color=INK, fontsize=15, weight='bold', va='center')
yF = WF[m, J].mean()
a2.text(1.01, (yF) / a2.get_ylim()[1], 'Frozen at its\naverage', transform=a2.transAxes, color=MUTED, fontsize=13, va='center')
a2.text(1.01, floor / a2.get_ylim()[1] - 0.05, f'floor {floor:.2f} ms\n(loop delay)', transform=a2.transAxes,
        color=MUTED, fontsize=11.5, va='center')

fig.text(0.085, 0.015, 'Baton-tip error over the 30 s: formula 0.68 mm, frozen 0.76 mm (10 % lower; 7 to 12 % over 5 noise seeds × 2 excerpts).\n'
         'Simulation of the servo bench, not hardware. Code: github.com/jpbronsard/syntonic-engine', fontsize=11.5, color=MUTED)
fig.savefig(sys.argv[2] if len(sys.argv) > 2 else 'maestro_accel_window.png', facecolor=SURF)
print('saved')
