import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import json
import os
from pathlib import Path

# Repo root (script lives in code/). NOTE: fig2_gate.png and fig4_power.png
# in figures/ are the authoritative high-replicate versions produced by
# regen_manuscript_figs.py from results/highrep/; this script regenerates the
# original 120/150-replicate versions of all four figures from results/
# (archived for the record).
base_dir = Path(__file__).resolve().parents[1]
out_dir = base_dir / "figures"
out_dir.mkdir(exist_ok=True)
data_dir = base_dir / "data"
res_dir = base_dir / "results"

plt.style.use('default')

# -----------------
# FIG 1: ESS Comparison
# -----------------
def fig1():
    schemes = ['Naive all-survey', 'Kepler-only\n(untrimmed)', 'Kepler-only\n(trimmed)']
    ess_vals = [487, 511, 1645]
    n_vals = [2848, 2135, 2135]
    rent = [17, 24, 77]
    
    fig, ax = plt.subplots(figsize=(8, 4.5))
    bars = ax.barh(schemes, ess_vals, color=['#aaaaaa', '#aaaaaa', '#1f77b4'],
                   height=0.6)
    
    # Add text labels at the right end of each bar (well clear of the y-axis labels)
    for bar, r in zip(bars, rent):
        width = bar.get_width()
        ax.text(width + 30, bar.get_y() + bar.get_height()/2,
                f'ESS = {width}  ({r}% retained)',
                ha='left', va='center', color='black', fontweight='bold')
                
    ax.set_xlabel('Effective Sample Size (ESS)')
    ax.set_title('Completeness Weighting Schemes')
    ax.set_xlim(0, max(ess_vals) * 1.30)
    ax.margins(y=0.15)
    plt.subplots_adjust(left=0.28)
    plt.tight_layout()
    plt.savefig(out_dir / 'fig1_ess.png', dpi=300)
    plt.close()

# -----------------
# FIG 2: Calibrated Gate Visualization
# -----------------
def fig2():
    with open(res_dir / 'task2_calibration_kepler_only_trim0.95.json') as f:
        d = json.load(f)
    
    beta_real = d['beta_real']
    null_med = d['null']['median']
    null_sd = (d['null']['q975'] - d['null']['q025']) / 3.92
    
    sig_med = d['signal']['median']
    sig_sd = (d['signal']['q975'] - d['signal']['q025']) / 3.92
    
    x = np.linspace(-0.5, 0.4, 500)
    y_null = np.exp(-0.5 * ((x - null_med) / null_sd)**2) / (null_sd * np.sqrt(2*np.pi))
    y_sig = np.exp(-0.5 * ((x - sig_med) / sig_sd)**2) / (sig_sd * np.sqrt(2*np.pi))
    
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.fill_between(x, y_sig, alpha=0.4, color='#d62728', label='Signal-Injected (SWEET-size effect)')
    ax.fill_between(x, y_null, alpha=0.4, color='#7f7f7f', label='Null-Injected (No effect)')
    
    ax.plot(x, y_sig, color='#d62728', lw=2)
    ax.plot(x, y_null, color='#7f7f7f', lw=2)
    
    ax.axvline(beta_real, color='k', linestyle='--', lw=2.5, label=rf'Real Data $\beta_{{age}}$ = {beta_real:.3f}')
    
    # Annotate percentiles in the upper-right where there is no curve mass
    ax.text(0.32, max(y_null)*0.55, 'Real coefficient:\n  43.3rd %ile of null\n  90.8th %ile of signal',
            ha='left', va='top', fontsize=10,
            bbox=dict(facecolor='white', alpha=0.85, edgecolor='gray', boxstyle='round,pad=0.4'))
    # Arrow from the annotation box down-left to the dashed line
    ax.annotate('', xy=(beta_real, max(y_null)*0.30), xytext=(0.30, max(y_null)*0.45),
                arrowprops=dict(arrowstyle='-', color='gray', lw=1))
            
    ax.set_xlabel(r'Age Coefficient $\beta_{{age}}$ (per IQR of $v_{{tan}}$)')
    ax.set_ylabel('Density')
    ax.set_title('Calibrated Gate: Genuinely Inconclusive Verdict')
    # Move legend to the lower-right where the curves have low density
    ax.legend(loc='lower right', framealpha=0.95)
    
    plt.tight_layout()
    plt.savefig(out_dir / 'fig2_gate.png', dpi=300)
    plt.close()

# -----------------
# FIG 3: Period Mixing 
# -----------------
def fig3():
    # Load actual data for true period distribution and SN fractions
    pl = pd.read_csv(data_dir / 'planet_sample_FGK.csv')
    k = pd.read_csv(data_dir / 'hosts_kinematics_FGK.csv')
    w = pd.read_csv(data_dir / 'completeness_weights_FGK.csv')
    
    # ensure string formats match for 'pl_name' matching
    df = pl.merge(k[['hostname', 'vtan']], on='hostname')
    df = df.merge(w, on='pl_name')
    df = df[df.hostname.str.startswith('Kepler')]
    
    med_vtan = df.vtan.median()
    df['kinematic_age'] = np.where(df.vtan >= med_vtan, 'Old (high v_tan)', 'Young (low v_tan)')
    df['is_SN'] = df.pl_rade >= 1.88
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Panel a: ECDFs of Period
    import matplotlib.ticker as ticker
    
    for age, color in [('Young (low v_tan)', '#1f77b4'), ('Old (high v_tan)', '#ff7f0e')]:
        subset = df[df.kinematic_age == age]
        p = np.sort(subset.pl_orbper.dropna())
        med_p = np.median(p)
        y = np.arange(1, len(p)+1) / len(p)
        ax1.plot(p, y, label=f'{age}', color=color, lw=2)
        
    ax1.set_xscale('log')
    ax1.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: str(int(x)) if x >= 1 else str(x)))
    ax1.set_xlabel('Orbital Period [days]')
    ax1.set_ylabel('Cumulative Probability')
    ax1.set_title('(a) Period Distribution Shift vs. Kinematic Age', pad=12)
    ax1.legend(loc='upper left')
    # Place median labels well-separated vertically so they cannot collide
    ax1.annotate('Young median\n10.05 d',
                 xy=(10.05, 0.50), xytext=(1.8, 0.10),
                 ha='center', color='#1f77b4', weight='bold', fontsize=10,
                 arrowprops=dict(arrowstyle='->', color='#1f77b4', lw=1))
    ax1.annotate('Old median\n10.69 d',
                 xy=(10.69, 0.50), xytext=(60, 0.10),
                 ha='center', color='#ff7f0e', weight='bold', fontsize=10,
                 arrowprops=dict(arrowstyle='->', color='#ff7f0e', lw=1))
    ax1.axvline(10.05, color='#1f77b4', linestyle=':', lw=1.2, alpha=0.7)
    ax1.axvline(10.69, color='#ff7f0e', linestyle=':', lw=1.2, alpha=0.7)
    
    # Panel b: SN Fraction vs Period
    bins = [0, 3, 10, 30, 100]
    bin_labels = ['<3', '3-10', '10-30', '30-100']
    
    fracs, errs = [], []
    for i in range(len(bins)-1):
        mask = (df.pl_orbper > bins[i]) & (df.pl_orbper <= bins[i+1])
        sub = df[mask]
        wt_tot = sub.w.sum()
        wt_sn = sub[sub.is_SN].w.sum()
        f = wt_sn / wt_tot if wt_tot > 0 else 0
        n_eff = (wt_tot**2) / (sub.w**2).sum() if len(sub) > 0 else 1
        se = np.sqrt(f * (1-f) / n_eff) if n_eff > 0 else 0
        fracs.append(f)
        errs.append(se)
        
    x = np.arange(len(bin_labels))
    ax2.bar(x, fracs, yerr=errs, capsize=5, color='#9467bd', alpha=0.8)
    ax2.set_xticks(x)
    ax2.set_xticklabels(bin_labels)
    ax2.set_xlabel('Orbital Period Bin [days]')
    ax2.set_ylabel('Completeness-Weighted SN Fraction')
    ax2.set_title(r'(b) Sub-Neptune Fraction vs. Orbital Period (fixed 1.88 R$_\oplus$ cut)', pad=12)
    ax2.set_ylim(0, 1.0)
    
    # Label the outer bins that manufacturing the spurious effect (0.141 and 0.837)
    for i, (fr, actual) in enumerate(zip(fracs, [0.141, None, None, 0.837])):
        if actual is not None:
             ax2.text(i, fr + 0.04, f'{actual:.3f}', ha='center', va='bottom', color='#4a2a6b', fontweight='bold')
        
    plt.tight_layout()
    plt.savefig(out_dir / 'fig3_period.png', dpi=300)
    plt.close()

# -----------------
# FIG 4: Power Grid
# -----------------
def fig4():
    with open(res_dir / 'power_analysis.json') as f:
        data = json.load(f)
        
    scales = {'0.5x': 0.5, '1x': 1, '2x': 2, '4x': 4, '8x': 8}
    df = pd.DataFrame(data)
    df['scale_num'] = df['scale'].map(scales)
    
    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    
    null_df = df[df.effect == 'null'].sort_values('scale_num')
    sweet_df = df[df.effect == 'SWEET'].sort_values('scale_num')
    strong_df = df[df.effect == 'strong'].sort_values('scale_num')
    
    ax.plot(null_df.scale_num, null_df.detect_rate * 100, marker='o', color='gray', linestyle=':', label='Null False Positive Rate', lw=2)
    ax.plot(sweet_df.scale_num, sweet_df.detect_rate * 100, marker='s', color='#d62728', label='SWEET-size (Literature) Effect', lw=2)
    ax.plot(strong_df.scale_num, strong_df.detect_rate * 100, marker='^', color='#2ca02c', label='Strong (2x Literature) Effect', lw=2)
    
    ax.axhline(80, color='k', linestyle='--', lw=1.5, alpha=0.5, label='Conventional 80% Power')
    ax.axhline(5, color='gray', linestyle='--', lw=1, alpha=0.5)
    
    # Mark 1x and 8x with thin vertical guide lines (not wide bands that hide data)
    ax.axvline(1, color='blue', linestyle=':', lw=1.5, alpha=0.7)
    ax.axvline(8, color='blue', linestyle=':', lw=1.5, alpha=0.7)
    
    # Place annotations HORIZONTALLY above the plot area, well clear of curves
    ax.annotate('Current Sample\n(~417 planets)',
                xy=(1, 95), xytext=(1, 95),
                ha='center', va='top', fontsize=10, color='blue', weight='bold',
                bbox=dict(facecolor='white', alpha=0.9, edgecolor='blue',
                          boxstyle='round,pad=0.3', linewidth=1))
    ax.annotate('Largest Scale Tested\n(~3,340 planets)',
                xy=(8, 95), xytext=(8, 95),
                ha='center', va='top', fontsize=10, color='blue', weight='bold',
                bbox=dict(facecolor='white', alpha=0.9, edgecolor='blue',
                          boxstyle='round,pad=0.3', linewidth=1))
    
    ax.set_xscale('log')
    ax.set_xticks([0.5, 1, 2, 4, 8])
    ax.set_xticklabels(['0.5x', '1x\n(Current)', '2x', '4x', '8x'])
    ax.set_xlabel('Sample Size Multiplier')
    ax.set_ylabel('Detection Rate (%)')
    ax.set_title('Power Grid vs. Sample Size for M-Dwarf Samples')
    
    ax.set_ylim(0, 125)
    ax.legend(loc='lower right', bbox_to_anchor=(1.0, 0.0), framealpha=0.95)
    
    plt.tight_layout()
    plt.savefig(out_dir / 'fig4_power.png', dpi=300)
    plt.close()

if __name__ == '__main__':
    fig1()
    fig2()
    fig3()
    fig4()
    print("All figures generated successfully.")