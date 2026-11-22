"""Regenerate the calibrated-gate figure from the production 1000-replicate CSVs.

Left:  FGK Kepler-only control, highrep/task2_gate_rep1000_replicates_raw.csv
Right: M/cool pooled, stratified_gate/M/strat_gate_M_vmag_*_r1000_raw.csv
All labels (n, arm separation, likelihood ratio) are computed from these CSVs.
"""
import json, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from _paths import BUNDLE, REPO, SCI, AUDIT, FIGS, PROV, DATA  # portable paths

FIG = str(FIGS / 'fig2_gate.png')
RES = str(REPO / 'results')
jt = json.load(open(os.path.join(RES, 'highrep', 'task2_calibration_kepler_only_trim0.95_rep1000.json')))
fgk = pd.read_csv(os.path.join(RES, 'highrep', 'task2_gate_rep1000_replicates_raw.csv'), keep_default_na=False)
b0 = fgk[fgk.arm == 'null'].beta_hat_per_iqr.values.astype(float)
b1 = fgk[fgk.arm == 'signal'].beta_hat_per_iqr.values.astype(float)
mn = pd.read_csv(os.path.join(RES, 'stratified_gate', 'M', 'strat_gate_M_vmag_null_r1000_raw.csv'))
ms = pd.read_csv(os.path.join(RES, 'stratified_gate', 'M', 'strat_gate_M_vmag_signal_r1000_raw.csv'))
m0 = mn.beta_pooled.dropna().values.astype(float)
m1 = ms.beta_pooled.dropna().values.astype(float)

real_f, real_m = jt['beta_real'], 0.04576344675022682

def stats_of(null, sig, real):
    sep = (np.median(sig) - np.median(null)) / np.std(null, ddof=1)
    k0, k1 = gaussian_kde(null), gaussian_kde(sig)
    lr = float(k1(real)[0] / k0(real)[0])
    pct = 100*np.mean(null <= real)
    return sep, lr, pct

sf, lf, pf = stats_of(b0, b1, real_f)
sm, lm, pm = stats_of(m0, m1, real_m)
print('FGK n=%d sep=%.3f LR=%.3f pct=%.1f' % (len(b0), sf, lf, pf))
print('M   n=%d sep=%.3f LR=%.3f pct=%.1f' % (len(m0), sm, lm, pm))

plt.rcParams.update({'font.size': 9})
fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
for ax, (null, sig, real, sep, lr, pct, title) in zip(axes, [
    (b0, b1, real_f, sf, lf, pf, 'FGK Kepler-only control'),
    (m0, m1, real_m, sm, lm, pm, 'Cool-host sample')]):
    bins = np.linspace(min(null.min(), sig.min()), max(null.max(), sig.max()), 45)
    ax.hist(sig, bins=bins, density=True, alpha=0.42, color='#D55E00', label='signal injected, n=1000')
    ax.hist(null, bins=bins, density=True, alpha=0.42, color='#7f7f7f', label='null injected, n=1000')
    ax.axvline(real, color='k', ls='--', lw=2, label='real coefficient')
    txt = ('real = %+.3f\nnull pct = %.1f\narm sep = %+.2f SD\nLR = %.2f' % (real, pct, sep, lr))
    ax.text(0.97, 0.95, txt, transform=ax.transAxes, ha='right', va='top', fontsize=8,
            bbox=dict(facecolor='white', alpha=0.88, edgecolor='gray', boxstyle='round,pad=0.35'))
    ax.set_xlabel(r'$\beta_{\rm age}$ per IQR of $v_{\rm tan}$')
    ax.set_ylabel('empirical density')
    ax.set_title(title)
    ax.legend(loc='upper left', fontsize=7.5, framealpha=0.9)
    ax.grid(alpha=0.2)
fig.suptitle('Calibrated statistical gate: production 1000-replicate ensembles', fontsize=11)
fig.tight_layout(rect=(0, 0, 1, 0.92))
fig.savefig(FIG, dpi=200)
print('wrote', FIG)
