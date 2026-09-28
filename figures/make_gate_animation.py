"""make_gate_animation.py — animated gate convergence for the README.

What it shows (all from committed, frozen inputs — nothing synthetic):
  Left:  null- vs signal-injected beta distributions accumulating
         replicate by replicate (recorded rng order), with the observed
         coefficient fixed throughout.
  Right: running percentile of the observed value inside each distribution,
         converging to the published 45.6% (null) / 81.2% (signal).

Provenance:
  raw_data: results/highrep/task2_gate_rep1000_replicates_raw.csv,
            results/highrep/task2_calibration_kepler_only_trim0.95_rep1000.json
  method: sorted rep order, 40 reps/arm per frame, 25 frames, Pillow GIF fps=5
  style: Okabe-Ito palette, color + linestyle redundancy, constrained layout
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import PillowWriter

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
HIGHREP = ROOT / "results" / "highrep"

NULL = "#999999"
SIGNAL = "#D55E00"
OBS = "black"

gate = pd.read_csv(HIGHREP / "task2_gate_rep1000_replicates_raw.csv",
                   keep_default_na=False)
t2 = json.loads((HIGHREP / "task2_calibration_kepler_only_trim0.95_rep1000.json").read_text())
beta_real = float(t2["beta_real"])
b_null = gate[gate.arm == "null"].sort_values("rep").beta_hat_per_iqr.values.astype(float)
b_sig = gate[gate.arm == "signal"].sort_values("rep").beta_hat_per_iqr.values.astype(float)
assert len(b_null) == 1000 and len(b_sig) == 1000

N_FRAMES = 25
STEP = 1000 // N_FRAMES
HBINS = np.linspace(-0.65, 0.45, 37)
FINAL_NULL = 45.6
FINAL_SIG = 81.2

plt.rcParams.update({
    "font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10,
    "xtick.labelsize": 9, "ytick.labelsize": 9, "legend.fontsize": 8,
    "axes.spines.top": False, "axes.spines.right": False,
})
fig, axes = plt.subplots(1, 2, figsize=(10, 4.4), layout="constrained")

xs = np.arange(STEP, 1001, STEP)
run_null = [100.0 * np.mean(b_null[:k] <= beta_real) for k in xs]
run_sig = [100.0 * np.mean(b_sig[:k] <= beta_real) for k in xs]

# Rebuild as ArtistAnimation on the original figure
from matplotlib.animation import ArtistAnimation

fig.clf()
axes = fig.subplots(1, 2)
artists = []
for f, k in enumerate(xs):
    ax = axes[0]
    _, _, patches_sig = ax.hist(b_sig[:k], bins=HBINS, density=True, alpha=0.45, color=SIGNAL,
                 edgecolor="white", linewidth=0.4)
    _, _, patches_null = ax.hist(b_null[:k], bins=HBINS, density=True, alpha=0.5, color=NULL,
                 edgecolor="white", linewidth=0.4)
    v = ax.axvline(beta_real, color=OBS, linestyle="--", lw=2)
    ax.set(xlabel=r"$\beta_{\mathrm{age}}$ (per IQR of $v_{\mathrm{tan}}$)",
           ylabel="Empirical density", title="(a) Gate distributions accumulate",
           xlim=(HBINS[0], HBINS[-1]), ylim=(0, 4.4))
    t1 = ax.text(0.97, 0.95, f"n = {k:,}/arm", transform=ax.transAxes,
                 ha="right", va="top", fontsize=9,
                 bbox=dict(facecolor="white", alpha=0.9, edgecolor="#BBBBBB",
                           boxstyle="round,pad=0.3"))
    ax = axes[1]
    (l1,) = ax.plot(xs[:f + 1], run_null[:f + 1], marker="o", ms=3, color=NULL,
                    linestyle=":", lw=1.6)
    (l2,) = ax.plot(xs[:f + 1], run_sig[:f + 1], marker="s", ms=3, color=SIGNAL, lw=1.8)
    h1 = ax.axhline(FINAL_NULL, color=NULL, lw=1, alpha=0.5)
    h2 = ax.axhline(FINAL_SIG, color=SIGNAL, lw=1, alpha=0.5)
    ax.set(xlabel="Replicates per arm", ylabel="Running percentile of observed",
           title="(b) Converging to 45.6% / 81.2% (inconclusive)",
           xlim=(0, 1000), ylim=(0, 100))
    st = fig.suptitle(f"Calibrated gate converging — {k:,}/1,000 replicates per arm",
                      fontsize=12, fontweight="bold")
    artists.append([*patches_sig, *patches_null, v, t1, l1, l2, h1, h2, st])

ani = ArtistAnimation(fig, artists, interval=220, repeat_delay=1500)
out = HERE / "gate_convergence.gif"
ani.save(out, writer=PillowWriter(fps=5))
print(f"wrote {out} ({out.stat().st_size / 1e6:.2f} MB), frames={len(artists)}")
