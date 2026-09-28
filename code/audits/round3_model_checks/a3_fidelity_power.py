"""Audit A3: generative fidelity of 11_power_analysis synthetic radii vs the
REAL M-sample radius distribution; internal consistency of power_analysis.json;
bincount/GH correctness spot-check."""
import pathlib, sys, json, math
import numpy as np
import pandas as pd
from scipy.stats import ks_2samp

ROOT = pathlib.Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT/"code"))
import importlib.util
spec = importlib.util.spec_from_file_location("pw", ROOT/"code"/"11_power_analysis.py")
pw = importlib.util.module_from_spec(spec); spec.loader.exec_module(pw)

base = pw.load_base()
rng = np.random.default_rng(12345)
real = base.pl_rade.values
VALLEY = 1.85
real_frac = float((real > VALLEY).mean())
print(f"real M: n={len(real)}, mean(pl_rade>1.85)={real_frac:.3f}, "
      f"median err1={np.median(base.pl_radeerr1.abs()):.3f}")

# --- synthetic pools at scale 1 ---
alpha0 = pw.calibrate_alpha(base, 0.0, min(0.495, real_frac/0.90), np.random.default_rng(999))
print(f"calibrated alpha(null)={alpha0:.3f}")
synth_radii, crossing = [], []
for r in range(40):
    d = pw.synth(base, 1.0, 0.0, alpha0, rng)
    synth_radii.append(d.pl_rade.values)
    crossing.append(float((d.pl_rade > VALLEY).mean()))
synth_all = np.concatenate(synth_radii)
ks = ks_2samp(real, synth_all)
print(f"\nKS real vs synthetic (scale1, null): D={ks.statistic:.4f} p={ks.pvalue:.2e} "
      f"(n_real={len(real)}, n_synth={len(synth_all)})")

def kde_modes(x, bw=0.09):
    grid = np.linspace(0.6, 4.2, 1200)
    dens = np.exp(-0.5*((grid[:,None]-x[None,:])/bw)**2).sum(1)
    dens = dens/np.trapezoid(dens, grid)
    # local maxima
    peaks = []
    for i in range(1, len(grid)-1):
        if dens[i] > dens[i-1] and dens[i] >= dens[i+1] and dens[i] > dens.max()*0.15:
            peaks.append((grid[i], dens[i]))
    return grid, dens, peaks

grid, dreal, preal = kde_modes(real)
_, dsyn, psyn = kde_modes(synth_all)
topr = sorted(preal, key=lambda t:-t[1])[:2]
tops = sorted(psyn, key=lambda t:-t[1])[:2]
print(f"KDE modes REAL      : {[f'{m:.2f}' for m,_ in sorted(topr)]}")
print(f"KDE modes SYNTHETIC : {[f'{m:.2f}' for m,_ in sorted(tops)]}")

bw_ = (real > 1.55) & (real < 2.15)
sw_ = (synth_all > 1.55) & (synth_all < 2.15)
print(f"boundary-window [1.55,2.15] mass: real {bw_.mean():.3f} vs synthetic {sw_.mean():.3f}")
print(f"valley-crossing frac: real {real_frac:.3f}; synthetic across reps "
      f"{np.mean(crossing):.3f} +- {np.std(crossing)/np.sqrt(len(crossing)):.3f}")

# analytic class-conditionals vs real modes: overlap across valley
from scipy.stats import norm as Nrm
MU_SE,S_SE,MU_SN,S_SN = 1.30,0.22,2.35,0.45
# reality-matched variant: modes at 1.23 / 2.08, widths from KDE of each side
below, above = real[real<1.6], real[real>2.1]
sd_se_real, sd_sn_real = below.std(), above.std()
print(f"\nreal SE-side scatter sd={sd_se_real:.2f} (x<1.6), SN-side sd={sd_sn_real:.2f} (x>2.1)")
for lab,(mse,ss,mss,ssn) in {"synthetic (1.30/0.22, 2.35/0.45)":(MU_SE,S_SE,MU_SN,S_SN),
                             "reality-match (1.30/0.22*, 2.08/0.45)":(MU_SE,S_SE,2.08,S_SN)}.items():
    # P(SE planet drawn above floor) and P(SN below floor) with typical err 0.10 convolution
    etyp = 0.10
    p_se_up = 1-Nrm.cdf((VALLEY-mse)/math.hypot(ss,etyp))
    p_sn_dn = Nrm.cdf((VALLEY-mss)/math.hypot(ssn,etyp))
    print(f"  {lab}: P(SE crosses up)={p_se_up:.3f}, P(SN crosses down)={p_sn_dn:.3f}, "
          f"sum(misclass proxies)={p_se_up+p_sn_dn:.3f}")

# --- internal consistency of power_analysis.json ---
resj = json.load(open(ROOT/"results"/"power_analysis.json"))
print(f"\npower_analysis.json cells: {len(resj)}")
print(f"{'cell':16s} {'detect':>7s} {'binom95%CI':>13s} {'coverage':>8s} {'med_bhat':>9s}")
for c in resj:
    p = c["detect_rate"]; n = c["reps"]
    se = math.sqrt(max(p*(1-p),1/n)/n); lo=max(0,p-1.96*se); hi=min(1,p+1.96*se)
    print(f"{c['scale']:>5s} {c['effect']:>8s} {p:7.1%}   [{lo:.1%},{hi:.1%}] {c['coverage']:8.1%} {c['median_beta_hat']:+9.3f}")
# monotonicity & attenuation notes
strong = {c["scale"]: c["detect_rate"] for c in resj if c["effect"]=="strong"}
sweet = {c["scale"]: c["detect_rate"] for c in resj if c["effect"]=="SWEET"}
null = {c["scale"]: c["detect_rate"] for c in resj if c["effect"]=="null"}
order = ["0.5x","1x","2x","4x","8x"]
print("\nmonotonicity in scale (detect should be non-decreasing):")
for nm,d in [("null",null),("SWEET",sweet),("strong",strong)]:
    vals=[d[o] for o in order]
    viol=[(order[i],order[i+1]) for i in range(len(vals)-1) if vals[i+1]<vals[i]]
    print(f"  {nm}: {['%.1f%%'%(v*100) for v in vals]} violations={viol}")
print("\nattenuation med_bhat/beta_inj (should ->1 with N):")
for o in order:
    s=[c["median_beta_hat"] for c in resj if c["effect"]=="strong" and c["scale"]==o][0]
    sw=[c["median_beta_hat"] for c in resj if c["effect"]=="SWEET" and c["scale"]==o][0]
    print(f"  {o}: strong {s/-0.28:.2f}, sweet {sw/-0.14:.2f}")
cov=np.array([c["coverage"] for c in resj])
print(f"coverage range [{cov.min():.1%},{cov.max():.1%}], mean {cov.mean():.1%} (nominal 95%)")

# --- bincount host-total spot check vs brute force ---
import numpy.testing as nt
rng2=np.random.default_rng(3)
nh=7; hidx=rng2.integers(0,nh,25); li=rng2.normal(size=25)
bc=np.bincount(hidx,weights=li,minlength=nh)
bf=np.array([li[hidx==h].sum() for h in range(nh)])
assert np.allclose(bc,bf), "bincount mismatch"
print("\nbincount host totals == brute force: OK (incl. single-planet hosts)")
