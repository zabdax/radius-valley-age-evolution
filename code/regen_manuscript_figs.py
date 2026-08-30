"""regen_manuscript_figs.py — Regenerate the two result figures from the
high-replicate ground-truth data (raw replicate CSVs / rep500 grid JSON).

fig2_gate.png   : empirical distributions of the 1,000 null- and 1,000
                  signal-injected gate replicates (raw CSV), with the real
                  coefficient and its empirical percentile placements.
fig4_power.png  : 15-cell power grid from the 500-replicate results
                  (identical plotting style to generate_figures.fig4).

Outputs: aastex_ms/fig2_gate.png, aastex_ms/fig4_power.png (manuscript
copies) and provenance copies in figures/ with _rep1000/_rep500 suffixes.
"""
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

BASE = Path(__file__).resolve().parents[1]   # repo root (script lives in code/)
HG = BASE / "results" / "highrep"            # authoritative high-replicate outputs
MS = BASE / "manuscript"
FIGS = BASE / "figures"
plt.style.use("default")

# ------------------------------------------------------------------
# FIG 2 — empirical gate distributions (n=1000 per arm)
# ------------------------------------------------------------------
gate = pd.read_csv(HG / "task2_gate_rep1000_replicates_raw.csv",
                   keep_default_na=False)
t2 = json.loads((HG / "task2_calibration_kepler_only_trim0.95_rep1000.json").read_text())
beta_real = t2["beta_real"]
b_null = gate[gate.arm == "null"].beta_hat_per_iqr.values.astype(float)
b_sig = gate[gate.arm == "signal"].beta_hat_per_iqr.values.astype(float)
pct_null = float(np.mean(b_null <= beta_real)) * 100
pct_sig = float(np.mean(b_sig <= beta_real)) * 100

fig, ax = plt.subplots(figsize=(8, 5))
bins = np.linspace(-0.65, 0.45, 37)
ax.hist(b_sig, bins=bins, density=True, alpha=0.40, color="#d62728",
        label="Signal-Injected, SWEET-size effect (n=1,000)")
ax.hist(b_null, bins=bins, density=True, alpha=0.45, color="#7f7f7f",
        label="Null-Injected, no effect (n=1,000)")
ax.axvline(beta_real, color="k", linestyle="--", lw=2.5,
           label=rf"Real Data $\beta_{{\rm age}}$ = {beta_real:.3f}")
ax.text(0.18, 2.55,
        f"Real coefficient:\n  {pct_null:.1f}th %ile of null\n  "
        f"{pct_sig:.1f}th %ile of signal",
        ha="left", va="top", fontsize=10,
        bbox=dict(facecolor="white", alpha=0.85, edgecolor="gray",
                  boxstyle="round,pad=0.4"))
ax.set_xlabel(r"Age Coefficient $\beta_{\rm age}$ (per IQR of $v_{\rm tan}$)")
ax.set_ylabel("Empirical Density")
ax.set_title("Calibrated Gate: Genuinely Inconclusive Verdict")
ax.legend(loc="upper left", framealpha=0.95)
plt.tight_layout()
plt.savefig(MS / "fig2_gate.png", dpi=300)
plt.savefig(FIGS / "fig2_gate_rep1000.png", dpi=300)
plt.close()
print(f"fig2 regenerated: null pct={pct_null:.1f}, signal pct={pct_sig:.1f}")

# ------------------------------------------------------------------
# FIG 4 — 500-replicate power grid (style identical to generate_figures)
# ------------------------------------------------------------------
with open(HG / "power_analysis_rep500.json") as f:
    data = json.load(f)
df = pd.DataFrame(data)
scales = {"0.5x": 0.5, "1x": 1, "2x": 2, "4x": 4, "8x": 8}
df["scale_num"] = df["scale"].map(scales)

fig, ax = plt.subplots(figsize=(8.5, 5.5))
null_df = df[df.effect == "null"].sort_values("scale_num")
sweet_df = df[df.effect == "SWEET"].sort_values("scale_num")
strong_df = df[df.effect == "strong"].sort_values("scale_num")

ax.plot(null_df.scale_num, null_df.detect_rate * 100, marker="o",
        color="gray", linestyle=":", label="Null False Positive Rate", lw=2)
ax.plot(sweet_df.scale_num, sweet_df.detect_rate * 100, marker="s",
        color="#d62728", label="SWEET-size (Literature) Effect", lw=2)
ax.plot(strong_df.scale_num, strong_df.detect_rate * 100, marker="^",
        color="#2ca02c", label="Strong (2x Literature) Effect", lw=2)

ax.axhline(80, color="k", linestyle="--", lw=1.5, alpha=0.5,
           label="Conventional 80% Power")
ax.axhline(5, color="gray", linestyle="--", lw=1, alpha=0.5)

ax.axvline(1, color="blue", linestyle=":", lw=1.5, alpha=0.7)
ax.axvline(8, color="blue", linestyle=":", lw=1.5, alpha=0.7)

ax.annotate("Current Sample\n(~417 planets)", xy=(1, 95), xytext=(1, 95),
            ha="center", va="top", fontsize=10, color="blue", weight="bold",
            bbox=dict(facecolor="white", alpha=0.9, edgecolor="blue",
                      boxstyle="round,pad=0.3", linewidth=1))
ax.annotate("Largest Scale Tested\n(~3,340 planets)", xy=(8, 95),
            xytext=(8, 95), ha="center", va="top", fontsize=10,
            color="blue", weight="bold",
            bbox=dict(facecolor="white", alpha=0.9, edgecolor="blue",
                      boxstyle="round,pad=0.3", linewidth=1))

ax.set_xscale("log")
ax.set_xticks([0.5, 1, 2, 4, 8])
ax.set_xticklabels(["0.5x", "1x\n(Current)", "2x", "4x", "8x"])
ax.set_xlabel("Sample Size Multiplier")
ax.set_ylabel("Detection Rate (%)")
ax.set_title("Power Grid vs. Sample Size for M-Dwarf Samples (500 replicates/cell)")
ax.set_ylim(0, 100)
ax.legend(loc="center left", framealpha=0.95)
plt.tight_layout()
plt.savefig(MS / "fig4_power.png", dpi=300)
plt.savefig(FIGS / "fig4_power_rep500.png", dpi=300)
plt.close()
print("fig4 regenerated from power_analysis_rep500.json")
