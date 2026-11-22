"""Animate the v2 power-versus-effect curve for the README.

All values come from the committed, authoritative v2 result file
``repro/v2/results/science/power_vs_effect_summary_FINAL.json``. Nothing is
synthetic: each frame reveals one of the ten real Monte Carlo points, sorted by
injected effect size, and shows that the gate verdict is FAIL at every point.

Panels:
  (a) signed arm separation (signal minus null medians, in null standard
      deviations) becoming more negative as the injected effect grows;
  (b) detection rate with Wilson 95 percent intervals, remaining far below the
      80 percent conventional-power reference.

The animation demonstrates that the instrument responds to a larger signal
(the arms separate and detection rises) while the real coefficient's null
placement keeps the verdict inconclusive throughout.

Usage: python figures/make_power_effect_animation.py
"""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
SOURCE = ROOT / "repro" / "v2" / "results" / "science" / "power_vs_effect_summary_FINAL.json"
OUT = HERE / "power_effect_curve.gif"

d = json.loads(SOURCE.read_text())
assert d["status"] == "COMPLETE" and d["points_completed"] == 10, "unexpected curve file"
pts = sorted(d["points"], key=lambda p: p["beta_inj"], reverse=True)

x = np.array([-p["beta_inj"] for p in pts])
sep = np.array([p["arm_separation_nullSD"] for p in pts])
det = np.array([100 * p["detect_rate"] for p in pts])
lo = np.array([100 * p["detect_rate_lo_wilson"] for p in pts])
hi = np.array([100 * p["detect_rate_hi_wilson"] for p in pts])

BLUE, ORANGE, GREY = "#0072B2", "#D55E00", "#777777"
plt.rcParams.update({
    "font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10,
    "xtick.labelsize": 9, "ytick.labelsize": 9, "legend.fontsize": 8,
    "axes.spines.top": False, "axes.spines.right": False,
})
fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.6), layout="constrained")
a, b = axes


def draw(f):
    for ax in axes:
        ax.clear()

    a.axhline(0, color=GREY, lw=0.8, alpha=0.6)
    a.plot(x[:f + 1], sep[:f + 1], "-o", color=BLUE, lw=2, ms=5)
    a.set(xlabel=r"injected effect $|\beta_{\rm inj}|$ per SD",
          ylabel="arm separation (null SD, signed)",
          title="(a) Arms separate as the signal grows",
          xlim=(-0.05, 1.28), ylim=(-2.4, 0.2))
    a.grid(alpha=0.2)

    b.axhline(80, color=GREY, ls=":", lw=1.2)
    b.text(1.24, 82, "80% power", color=GREY, ha="right", fontsize=8)
    b.fill_between(x[:f + 1], lo[:f + 1], hi[:f + 1], color=ORANGE, alpha=0.18)
    b.plot(x[:f + 1], det[:f + 1], "-s", color=ORANGE, lw=2, ms=4.5)
    b.set(xlabel=r"injected effect $|\beta_{\rm inj}|$ per SD",
          ylabel="detection rate [%]",
          title="(b) Detection rises but stays far below 80%",
          xlim=(-0.05, 1.28), ylim=(0, 100))
    b.grid(alpha=0.2)
    b.text(0.03, 0.05,
           "Real coefficient stays null-like;\nno injection changes the verdict.",
           transform=b.transAxes, fontsize=8, color="#333333",
           bbox=dict(facecolor="white", alpha=0.9, edgecolor="#BBBBBB",
                     boxstyle="round,pad=0.3"))

    lab = pts[f]["point"].replace("_", " ")
    fig.suptitle(
        f"Power versus injected effect: {f + 1}/10 points shown, "
        f"verdict = {pts[f]['verdict']} ({lab})",
        fontsize=11, fontweight="bold")
    return []


ani = FuncAnimation(fig, draw, frames=len(pts), interval=520,
                    repeat_delay=2200, blit=False)
ani.save(OUT, writer=PillowWriter(fps=3))
print(f"wrote {OUT} ({OUT.stat().st_size / 1e6:.2f} MB), frames={len(pts)}")
