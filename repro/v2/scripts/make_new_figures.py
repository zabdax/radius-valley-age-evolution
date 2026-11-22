"""Publication-ready figures for the P2/P3 science upgrades.
Okabe-Ito palette; no em dashes in labels; explicit uncertainties.
"""
import json, pathlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from _paths import BUNDLE, REPO, SCI, AUDIT, FIGS, PROV, DATA  # portable paths

SCI = SCI
FIG = FIGS
FIG.mkdir(exist_ok=True)
BLUE, ORANGE, GREY = '#0072B2', '#D55E00', '#777777'

# ---------- Figure 5: TIME Table age coverage ----------
cm = json.load(open(SCI / 'time_table_crossmatch.json'))
matched = cm['matched_hosts']
rows = []
for m in matched:
    statuses = m.get('age_statuses', ['age_estimate'] * len(m['age_Gyr']))
    # One point per sample host: K2-148 A/B are two TIME Table rows for one host.
    for name, status, age, err, teff in zip(m['sample_hosts'], statuses, m['age_Gyr'],
                                            m['e_age_Gyr'], m['teff_K']):
        if status != 'age_estimate' or age is None or not np.isfinite(age):
            continue
        rows.append(dict(name=name, age=float(age), err=float(err) if err is not None else np.nan,
                         teff=float(teff), existing=bool(m['already_one_of_12_gyro'])))
df = pd.DataFrame(rows)
fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.6))
a = axes[0]
for existing, color, label in [(False, BLUE, '62 additional matches'), (True, ORANGE, '3 pipeline gyro hosts with finite ages')]:
    s = df[df.existing == existing]
    a.errorbar(s.teff, s.age, yerr=s.err, fmt='o', ms=3.2, elinewidth=0.6,
               color=color, label=label, alpha=0.8, capsize=0)
a.axvline(3800, color=GREY, ls='--', lw=1, label='Angus+19 validity floor')
a.set_xlabel(r'$T_{\rm eff}$ [K]'); a.set_ylabel('TIME Table age [Gyr]')
a.set_title('(a) Gyro ages for matched cool hosts')
a.legend(fontsize=6.5, loc='upper right'); a.grid(alpha=0.2)
b = axes[1]
frac = df.err / df.age
b.hist(frac[~df.existing], bins=np.arange(0, 0.9, 0.05), color=BLUE, alpha=0.75,
       label='new matches')
b.hist(frac[df.existing], bins=np.arange(0, 0.9, 0.05), color=ORANGE, alpha=0.8,
       label='pipeline gyro hosts')
b.axvline(0.20, color=GREY, ls='--', lw=1)
b.text(0.30, b.get_ylim()[1]*0.92, 'median 20%', fontsize=7, color=GREY)
b.set_xlabel('fractional age uncertainty'); b.set_ylabel('hosts')
b.set_title('(b) Age-precision distribution'); b.legend(fontsize=7); b.grid(alpha=0.2)
fig.tight_layout()
fig.savefig(FIG / 'fig5_time_table_ages.png', dpi=180)
print('wrote fig5', len(df), 'finite-age matches')

# ---------- Figure 6: host distribution and the mass boundary ----------
D = DATA
ps = pd.read_csv(D / 'planet_sample_M.csv')
hosts = ps.drop_duplicates('hostname')
teff = hosts.st_teff.dropna().values
mass = hosts.st_mass.dropna().values
fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.4))
a = axes[0]
a.hist(teff, bins=np.arange(2550, 4250, 100), color=BLUE, alpha=0.85)
a.axvline(3400, color=ORANGE, ls='--', lw=1.3)
a.axvline(4200, color=GREY, ls=':', lw=1.3)
a.text(3410, a.get_ylim()[1]*0.85, '3400 K', color=ORANGE, fontsize=7)
a.text(4120, a.get_ylim()[1]*0.6, '4200 K', color=GREY, fontsize=7)
a.set_xlabel(r'host $T_{\rm eff}$ [K]'); a.set_ylabel('hosts')
a.set_title('(a) Temperature (303 hosts)'); a.grid(alpha=0.2)
b = axes[1]
b.hist(mass, bins=np.arange(0.15, 0.85, 0.04), color=BLUE, alpha=0.85)
for x, lab in [(0.3, '0.3'), (0.4, '0.4')]:
    b.axvline(x, color=ORANGE, ls='--', lw=1.3)
    b.text(x+0.005, b.get_ylim()[1]*0.88, lab+r' $M_\odot$', color=ORANGE, fontsize=7)
b.set_xlabel(r'host mass [$M_\odot$]'); b.set_ylabel('hosts')
b.set_title(f'(b) Mass ({len(mass)} hosts with usable mass)'); b.grid(alpha=0.2)
fig.tight_layout()
fig.savefig(FIG / 'fig6_host_distribution.png', dpi=180)
print('wrote fig6', len(teff), len(mass))
