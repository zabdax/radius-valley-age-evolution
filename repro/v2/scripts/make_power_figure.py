import json, pathlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

SCI = pathlib.Path(r'REPO\manuscripts\q1-paper\v2-2026-11\results\science')
check = json.load(open(SCI / 'effect_size_derivation.json'))['bracket_points']
order = [
    ('null', 0.0),
    ('half_primary', check['half_primary']),
    ('Kamulali_primary', check['Kamulali_primary']),
    ('Gaidos_rho025_vtan', check['Gaidos_rho025_vtan']),
    ('Berger_FGK', check['Berger_FGK']),
    ('Gaidos_rho050_vtan', check['Gaidos_rho050_vtan']),
    ('Gaidos_raw_true_age', check['Gaidos_raw_true_age']),
]
x = [-v for _, v in order]
y = [-0.19 * (v / 0.14) for _, v in order]  # linear response scaled from archived -0.14 point
labels = [n.replace('_', ' ') for n, _ in order]

fig, ax = plt.subplots(figsize=(7, 4.4))
ax.plot(x, y, marker='o', color='#0072B2', lw=1.8)
ax.axhline(-1.0, color='#999999', ls=':', lw=1)
for xi, yi, lab in zip(x, y, labels):
    ax.annotate(lab, (xi, yi), xytext=(4, -12), textcoords='offset points', fontsize=7)
ax.set_xlabel('Injected effect magnitude |beta| per SD of v_tan')
ax.set_ylabel('Arm separation (null SD, signed)')
ax.set_title('Gate sensitivity vs effect size (linear-response estimate)')
ax.grid(alpha=0.25)
fig.tight_layout()
fig.savefig(SCI / 'power_vs_effect.png', dpi=180)
print('wrote', SCI / 'power_vs_effect.png')
