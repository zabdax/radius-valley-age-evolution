# Figures and animations

Every visual here is rebuilt from committed, frozen data. No value is drawn by
hand. Scripts use the Okabe-Ito colour-blind-safe palette, redundant
colour-plus-marker encoding, and explicit uncertainty intervals.

## Animated (featured in the root README)

| File | Script | What it shows |
|---|---|---|
| `gate_convergence.gif` | `make_gate_animation.py` | The FGK gate converging to the published 45.6% / 81.2% percentiles as 1,000 replicates per arm accumulate. |
| `power_effect_curve.gif` | `make_power_effect_animation.py` | The ten-point v2 power-versus-effect curve building point by point: arms separate and detection rises, but the verdict stays FAIL throughout. |

Rebuild with:

```bash
python figures/make_gate_animation.py
python figures/make_power_effect_animation.py
```

## Static

| File | Purpose |
|---|---|
| `overview.png`, `overview.pdf` | Hero three-panel summary (gate, power grid, frozen samples). Built by `make_overview.py`. |
| `fig1_ess.png` | Effective-sample-size comparison across weighting schemes (earlier analysis). |
| `fig2_gate.png`, `fig3_period.png`, `fig4_power.png` | Earlier publication figures; the current v2 versions live in `repro/v2/figures/`. |

## Current v2 publication figures

`repro/v2/figures/fig1` through `fig7` are the seven figures in the v2
manuscript:

1. radius distributions split by kinematic age
2. TIME Table gyrochronological ages
3. period-mixing diagnostic
4. calibrated gate (FGK and cool hosts)
5. null false-positive / coverage diagnostics
6. power versus sample size
7. power versus injected effect

Their regeneration scripts are in `repro/v2/scripts/` (`regen_fig1_large.py`,
`make_new_figures.py`, `regen_gate_fig.py`, `regen_fig7_sorted.py`).
