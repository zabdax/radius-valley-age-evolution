"""Regenerate the power-versus-effect figure with a monotonic x-axis.

Reads the authoritative 10-point curve and sorts by injected effect before
plotting, so neither curve doubles back on itself.
"""
import json, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from _paths import BUNDLE, REPO, SCI, AUDIT, FIGS, PROV, DATA  # portable paths

V2 = str(BUNDLE)
d = json.load(open(os.path.join(V2, 'results', 'science', 'power_vs_effect_summary_FINAL.json')))
pts = sorted(d['points'], key=lambda p: p['beta_inj'], reverse=True)  # 0 -> -1.20
x = np.array([-p['beta_inj'] for p in pts])
sep = np.array([p['arm_separation_nullSD'] for p in pts])
det = np.array([100*p['detect_rate'] for p in pts])
elo = np.array([100*p['detect_rate_lo_wilson'] for p in pts])
ehi = np.array([100*p['detect_rate_hi_wilson'] for p in pts])

fig, ax1 = plt.subplots(figsize=(7.4, 4.6))
c1, c2 = '#0072B2', '#D55E00'
l1, = ax1.plot(x, sep, '-o', color=c1, lw=2, ms=5, label='signed arm separation (left axis)')
ax1.set_xlabel(r'injected effect $|\beta_{\rm inj}|$ per SD of $v_{\rm tan}$')
ax1.set_ylabel('arm separation (null SD, signed)', color=c1)
ax1.tick_params(axis='y', labelcolor=c1)
ax1.grid(alpha=0.22)
ax2 = ax1.twinx()
l2, = ax2.plot(x, det, '-s', color=c2, lw=2, ms=4.5, label='detection rate (right axis)')
ax2.fill_between(x, elo, ehi, color=c2, alpha=0.15)
ax2.set_ylabel('detection rate [%] (Wilson 95%)', color=c2)
ax2.tick_params(axis='y', labelcolor=c2)
ax2.set_ylim(0, 100)
# Label only three cleanly separated reference points; the rest are listed in
# the caption. Offsets are chosen to avoid both curves and each other.
named = {
    'Kamulali_primary': ('Kamulali', (7, 14), 'left', 'bottom'),
    # Both Gaidos labels are moved into open space with thin leaders so neither
    # touches the blue curve or the orange Wilson band at their true markers.
    'Gaidos_rho050_vtan': (r'Gaidos $\rho=0.5$', (12, -44), 'left', 'top'),
    'Gaidos_raw_true_age': (r'Gaidos true age', (-8, -42), 'right', 'top'),
}
for xi, p in zip(x, pts):
    if p['point'] not in named:
        continue
    label, offset, ha, va = named[p['point']]
    ax1.annotate(label, (xi, p['arm_separation_nullSD']), xytext=offset,
                 textcoords='offset points', ha=ha, va=va, fontsize=8,
                 bbox=dict(facecolor='white', alpha=0.9, edgecolor='none', pad=1.2),
                 arrowprops=dict(arrowstyle='-', color='#666666', lw=0.6,
                                 shrinkA=2, shrinkB=4))
ax1.legend([l1, l2], ['signed arm separation', 'detection rate'],
           loc='center left', bbox_to_anchor=(0.02, 0.45), fontsize=8,
           framealpha=0.95)
fig.suptitle('Power versus injected effect: production cool-host gate (150 replicates per arm)', fontsize=10.5)
fig.tight_layout(rect=(0, 0.06, 1, 0.94))
out = os.path.join(V2, 'figures', 'fig7_power_vs_effect.png')
fig.savefig(out, dpi=200)
print('wrote', out)
