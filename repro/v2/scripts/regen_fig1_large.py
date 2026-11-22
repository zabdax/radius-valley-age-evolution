"""Regenerate Figure 1 (radius distributions) at publication-readable size.

Frozen inputs only, read in place. Weighted histograms split at the median
tangential velocity for the FGK Kepler-only control and the cool-star sample.
"""
import pathlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from _paths import BUNDLE, REPO, SCI, AUDIT, FIGS, PROV, DATA  # portable paths

D = DATA
OUT = FIGS / 'fig1_radius_valley.png'

plt.rcParams.update({
    'font.size': 13, 'axes.titlesize': 15, 'axes.labelsize': 14,
    'xtick.labelsize': 12, 'ytick.labelsize': 12, 'legend.fontsize': 12,
    'figure.titlesize': 16, 'axes.linewidth': 1.2,
})
BLUE, ORANGE = '#0072B2', '#D55E00'


def load(tag, kepler_only):
    ps = pd.read_csv(D / f'planet_sample_{tag}.csv')
    kin = pd.read_csv(D / f'hosts_kinematics_{tag}.csv')
    wt = pd.read_csv(D / f'completeness_weights_{tag}.csv').rename(columns={'w': 'wgt'})
    d = ps.merge(kin, on='hostname', how='inner').merge(wt, on='pl_name', how='inner')
    if kepler_only:
        d = d[d.hostname.str.startswith('Kepler')].copy()
    cap = d.wgt.quantile(0.95)
    d['ww'] = np.minimum(d.wgt.fillna(0), cap)
    return d


def weighted_hist(ax, d, color, label, cut):
    young = d[d.vtan <= d.vtan.median()]
    old = d[d.vtan > d.vtan.median()]
    bins = np.linspace(0.5, 4.0, 30)
    for df, lab, style in [(young, 'young (low $v_{\\rm tan}$)', '-'),
                           (old, 'old (high $v_{\\rm tan}$)', '--')]:
        ax.hist(df.pl_rade, bins=bins, weights=df.ww, density=True,
                histtype='step', linewidth=2.2, linestyle=style, color=color,
                label=lab)
    ax.axvline(cut, color='#555555', ls=':', lw=1.8)
    ax.set_xlabel(r'planet radius [$R_\oplus$]')
    ax.set_ylabel('normalized weighted density')
    ax.legend(loc='upper left', frameon=False)
    ax.grid(alpha=0.2)


fgk = load('FGK', True)
cool = load('M', False)
fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.6))
weighted_hist(axes[0], fgk, BLUE, 'FGK', 1.88)
axes[0].set_title(f'FGK control (N={len(fgk)}, ESS$\\simeq$1645)')
axes[0].text(1.90, axes[0].get_ylim()[1]*0.62, '1.88 $R_\\oplus$ cut',
             color='#555555', fontsize=10, rotation=90, va='bottom',
             bbox=dict(facecolor='white', alpha=0.78, edgecolor='none', pad=1))
weighted_hist(axes[1], cool, ORANGE, 'cool', 1.85)
axes[1].set_title(f'Cool hosts (N={len(cool)}, ESS$\\simeq$140)')
axes[1].text(1.87, axes[1].get_ylim()[1]*0.62, '1.85 $R_\\oplus$ cut',
             color='#555555', fontsize=10, rotation=90, va='bottom',
             bbox=dict(facecolor='white', alpha=0.78, edgecolor='none', pad=1))
fig.suptitle('Completeness-weighted radius distributions split at median $v_{\\rm tan}$',
             fontsize=16)
fig.tight_layout(rect=(0, 0, 1, 0.94))
fig.savefig(OUT, dpi=200)
print('wrote', OUT, 'size', OUT.stat().st_size)
