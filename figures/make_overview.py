"""make_overview.py — hero overview figure for the open-source release.

Panels (all from frozen, committed inputs — no synthetic values):
  (a) Calibrated FGK gate: 1,000 null- vs 1,000 signal-injected replicate
      beta_age/IQR distributions + observed coefficient.
  (b) M-dwarf power grid: per-cell detection rates (500 reps/cell) with
      Clopper-Pearson 95% intervals + 80% conventional-power reference.
  (c) Sample inventory: planet/host counts from the frozen CSVs.

Provenance:
  raw_data: data/planet_sample_FGK.csv, data/planet_sample_M.csv,
            results/highrep/task2_gate_rep1000_replicates_raw.csv,
            results/highrep/power_grid_rep500_replicates_raw.csv,
            results/highrep/task2_calibration_kepler_only_trim0.95_rep1000.json
  uncertainty: (a) empirical 2.5/97.5 percentiles; (b) Clopper-Pearson 95% CI
  style: Okabe-Ito palette, color + marker/linestyle redundancy, constrained layout
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from scipy import stats as st

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent

HIGHREP = ROOT / "results" / "highrep"
DATA = ROOT / "data"

# Okabe-Ito (colorblind-safe) on white
NULL = "#999999"
SIGNAL = "#D55E00"
SWEET = "#D55E00"
STRONG = "#009E73"
OBS = "black"
GRID_NULL = "#777777"

gate = pd.read_csv(HIGHREP / "task2_gate_rep1000_replicates_raw.csv",
                   keep_default_na=False)
t2 = json.loads((HIGHREP / "task2_calibration_kepler_only_trim0.95_rep1000.json").read_text())
beta_real = float(t2["beta_real"])
se_real = float(t2["se_real"])
b_null = gate[gate.arm == "null"].sort_values("rep").beta_hat_per_iqr.values.astype(float)
b_sig = gate[gate.arm == "signal"].sort_values("rep").beta_hat_per_iqr.values.astype(float)
pct_null = 100.0 * float(np.mean(b_null <= beta_real))
pct_sig = 100.0 * float(np.mean(b_sig <= beta_real))
null_med, null_lo, null_hi = float(np.median(b_null)), float(np.percentile(b_null, 2.5)), float(np.percentile(b_null, 97.5))
sig_med, sig_lo, sig_hi = float(np.median(b_sig)), float(np.percentile(b_sig, 2.5)), float(np.percentile(b_sig, 97.5))

grid = pd.read_csv(HIGHREP / "power_grid_rep500_replicates_raw.csv", keep_default_na=False)
SCALE = {"0.5x": 0.5, "1x": 1, "2x": 2, "4x": 4, "8x": 8}
rows = []
for s, sv in SCALE.items():
    for e in ("null", "SWEET", "strong"):
        c = grid[grid.cell == f"{s}|{e}"]
        n = len(c)
        dr = float(c.detect.mean())
        x = int(round(dr * n))
        lo = float(st.beta.ppf(0.025, x, n - x + 1)) if x > 0 else 0.0
        hi = float(st.beta.ppf(0.975, x + 1, n - x)) if x < n else 1.0
        rows.append({"scale": s, "num": sv, "effect": e, "rate": dr, "lo": lo, "hi": hi, "n": n,
                     "npl": float(c.n_planets.mean())})
powdf = pd.DataFrame(rows)

pl_fgk = pd.read_csv(DATA / "planet_sample_FGK.csv")
pl_m = pd.read_csv(DATA / "planet_sample_M.csv")
n_fgk, n_m = len(pl_fgk), len(pl_m)
n_h_fgk, n_h_m = pl_fgk.hostname.nunique(), pl_m.hostname.nunique()

plt.rcParams.update({
    "figure.dpi": 150, "savefig.dpi": 300,
    "font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10,
    "xtick.labelsize": 9, "ytick.labelsize": 9, "legend.fontsize": 8,
    "axes.spines.top": False, "axes.spines.right": False,
})
fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.6), layout="constrained",
                         gridspec_kw={"width_ratios": [1.15, 1.0, 0.85]})

# (a) gate
ax = axes[0]
hbins = np.linspace(-0.65, 0.45, 37)
ax.hist(b_sig, bins=hbins, density=True, alpha=0.45, color=SIGNAL,
        label="Signal-injected (n=1,000)", edgecolor="white", linewidth=0.4)
ax.hist(b_null, bins=hbins, density=True, alpha=0.5, color=NULL,
        label="Null-injected (n=1,000)", edgecolor="white", linewidth=0.4)
ax.axvline(beta_real, color=OBS, linestyle="--", lw=2,
           label=f"Observed $\\beta_{{age}}$ = {beta_real:.3f} ± {se_real:.3f}")
ax.text(0.97, 0.97, f"{pct_null:.1f}th %ile of null\n{pct_sig:.1f}th %ile of signal\nverdict: FAIL (inconclusive)",
        transform=ax.transAxes, ha="right", va="top", fontsize=8,
        bbox=dict(facecolor="white", alpha=0.9, edgecolor="#BBBBBB", boxstyle="round,pad=0.35"))
ax.set_xlabel(r"$\beta_{\mathrm{age}}$ (per IQR of $v_{\mathrm{tan}}$)")
ax.set_ylabel("Empirical density")
ax.set_title("(a) FGK calibrated gate — inconclusive", fontweight="bold")
ax.legend(loc="upper left", framealpha=0.95)

# (b) power grid with 95% CI
ax = axes[1]
for eff, color, marker, ls, lab in [
    ("null", GRID_NULL, "o", ":", "Null false-positive rate"),
    ("SWEET", SWEET, "s", "-", "Literature-sized effect"),
    ("strong", STRONG, "^", "-", "Strong effect (2×)"),
]:
    sel = powdf[powdf.effect == eff].sort_values("num")
    yerr = np.vstack([sel.rate - sel.lo, sel.hi - sel.rate])
    ax.errorbar(sel.num, sel.rate * 100, yerr=yerr * 100, marker=marker, color=color,
                linestyle=ls, lw=2, capsize=3, label=lab)
ax.axhline(80, color="black", linestyle="--", lw=1.2, alpha=0.6, label="Conventional 80% power")
ax.axvline(8, color="#0072B2", linestyle=":", lw=1.4, alpha=0.7)
ax.annotate("8× (~3,337 planets)\nSWEET 23.2%, strong 60.6%", xy=(7.7, 63), xytext=(1.1, 68),
            ha="left", fontsize=8, arrowprops=dict(arrowstyle="->", color="gray", lw=1))
ax.set_xscale("log")
ax.set_xticks([0.5, 1, 2, 4, 8])
ax.set_xticklabels(["0.5×", "1×", "2×", "4×", "8×"])
ax.set_ylim(0, 100)
ax.set_xlabel("Sample-size multiplier")
ax.set_ylabel("Detection rate (%) with 95% CI")
ax.set_title("(c)" if False else "(b) M-dwarf power grid", fontweight="bold")
ax.legend(loc="upper left", framealpha=0.95)

# (c) sample inventory (counts from frozen CSVs — no inference)
ax = axes[2]
cats = ["FGK control\nplanets", "FGK hosts", "M-dwarf\nplanets", "M hosts"]
vals = [n_fgk, n_h_fgk, n_m, n_h_m]
colors = ["#0072B2", "#56B4E9", "#E69F00", "#F0C060"]
hatches = ["", "///", "", "///"]
bars = ax.bar(cats, vals, color=colors, edgecolor="black", linewidth=0.8)
for b, h in zip(bars, hatches):
    b.set_hatch(h)
for x, v in zip(cats, vals):
    ax.text(x, v + max(vals) * 0.02, f"{v:,}", ha="center", va="bottom", fontsize=9, fontweight="bold")
ax.set_ylabel("Count (frozen CSVs)")
ax.set_title("(c) Samples behind the result", fontweight="bold")
ax.text(0.02, 0.96, f"FGK Kepler-only gate: N=2,135, ESS=1,645\nM regime: {n_m} planets / {n_h_m} hosts",
        transform=ax.transAxes, va="top", fontsize=8,
        bbox=dict(facecolor="white", alpha=0.9, edgecolor="#BBBBBB", boxstyle="round,pad=0.35"))

fig.suptitle("Radius valley age evolution around M dwarfs — calibrated gate + power + samples", fontsize=12, fontweight="bold")
out_png = HERE / "overview.png"
out_pdf = HERE / "overview.pdf"
fig.savefig(out_png)
fig.savefig(out_pdf)
print(f"wrote {out_png} | {out_pdf}")
print(f"gate: beta={beta_real:.4f}±{se_real:.4f} pct_null={pct_null:.1f} pct_sig={pct_sig:.1f}")
print(f"samples: FGK {n_fgk}/{n_h_fgk} M {n_m}/{n_h_m}")
