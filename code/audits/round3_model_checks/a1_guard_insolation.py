"""Audit A1: z-guard impact, insolation formula check, error-column chain,
soft-weight ambiguity census -- replicating build_design exactly."""
import pathlib, sys
import numpy as np
import pandas as pd

D = pathlib.Path(__file__).resolve().parents[3] / "data"

def load(tag, kin_col="vtan"):
    pl = pd.read_csv(D / f"planet_sample_{tag}.csv")
    k = pd.read_csv(D / f"hosts_kinematics_{tag}.csv")
    df = pl.merge(k[["hostname", "vtan", "W"]], on="hostname", how="inner")
    if kin_col == "absW":
        df = df.dropna(subset=["W"]).copy()
        df["absW"] = df.W.abs()
    return df

for tag, kin_col, label in [("M","vtan","M_vtan"), ("M","absW","M_absW"), ("FGK","vtan","FGK_vtan")]:
    df = load(tag, kin_col)
    if tag == "M" and False:
        pass
    n = len(df)
    print(f"\n================ {label}: {n} planets / {df.hostname.nunique()} hosts ================")

    # --- error column chain ---
    err = df.get("pl_radeerr1"); err2 = df.get("pl_radeerr2")
    e1 = np.where(np.isfinite(err.fillna(np.nan)), err.abs(), np.nan)
    used_err1 = np.isfinite(e1)
    e = e1.copy()
    if err2 is not None:
        take2 = ~np.isfinite(e) & np.isfinite(err2.abs())
        e = np.where(np.isfinite(e), e, err2.abs())
    else:
        take2 = np.zeros(n, bool)
    final = np.isfinite(e) & (e > 0.01)
    e_final = np.where(final, e, 0.05 * df.pl_rade)
    print(f"error chain: err1 used {used_err1.sum()} ({100*used_err1.mean():.1f}%), "
          f"fell to err2 {take2.sum()}, fell to 5% radius {(~final).sum()}")
    print(f"  e_final: min {e_final.min():.4f}, med {np.median(e_final):.4f}, "
          f"p95 {np.percentile(e_final,95):.4f}, max {e_final.max():.4f}")
    # asymmetric-error note: |err2| vs err1 when both present
    both = err.notna() & err2.notna()
    ratio = (err2[both].abs() / err[both].abs()).replace([np.inf], np.nan).dropna()
    print(f"  asymmetry |err2|/|err1| median {ratio.median():.2f}; "
          f"|err2|>|err1| in {(ratio>1.05).mean()*100:.0f}% of both-present rows")

    # --- soft weight ---
    valley = 1.85 if tag=="M" else 1.88
    w_sn = 0.5*(1+__import__("scipy.special", fromlist=["erf"]).erf((df.pl_rade.values-valley)/(e_final*np.sqrt(2))))
    amb = (w_sn>0.05)&(w_sn<0.95)
    print(f"soft weight w_sn: ambiguous(0.05<w<0.95) {amb.sum()} ({100*amb.mean():.1f}%); "
          f"w in (0.005,0.995): {((w_sn>0.005)&(w_sn<0.995)).sum()}")

    # --- covariate guard census (exact replication of build_design logic) ---
    mstar = df.st_mass.fillna(0.5).values
    mstar_c = np.clip(mstar, 0.08, 1.6)
    p_yr = df.pl_orbper.values/365.25
    a_au = np.cbrt(mstar_c)*p_yr**(2.0/3.0)
    rstar = df.st_rad.fillna(0.6).values
    tstar = df.st_teff.fillna(4000).values
    lum_buggy = np.clip(rstar,0.1,4)**2 * np.clip(tstar,2500,7500)/5772**4
    lum_correct = np.clip(rstar,0.1,4)**2 * (np.clip(tstar,2500,7500)/5772)**4
    s_in_b = lum_buggy/a_au**2
    s_in_c = lum_correct/a_au**2
    lp = np.log10(df.pl_orbper.values); ls_b = np.log10(s_in_b); ls_c = np.log10(s_in_c)

    def z(v):
        v=np.asarray(v,float); return (v-v.mean())/(v.std()+1e-12)
    vz_lp = z(lp); vz_lsb = z(ls_b); vz_lsc = z(ls_c)
    g_lp = (lp<-9)|(lp>9); g_lsb = (ls_b<-9)|(ls_b>9); g_lsc = (ls_c<-9)|(ls_c>9)
    print(f"logP guard: zeroed {g_lp.sum()}/{n}")
    print(f"BUGGY logS (=R^2*T/5772^4): raw range [{ls_b.min():.2f},{ls_b.max():.2f}], "
          f"zeroed by guard {g_lsb.sum()}/{n} ({100*g_lsb.mean():.1f}%)")
    print(f"CORRECT logS ((T/5772)^4): raw range [{ls_c.min():.2f},{ls_c.max():.2f}], "
          f"zeroed by guard {g_lsc.sum()}/{n}")
    surv = ~g_lsb
    if surv.sum()>0:
        # among surviving (non-zeroed) rows, how does buggy-z relate to correct-z?
        rb = np.corrcoef(vz_lsb[surv], vz_lsc[surv])[0,1]
        print(f"  survivors' buggy-z vs correct-z Pearson r={rb:.3f}; "
              f"survivor correct-S median {10**np.median(ls_c[surv]):.1f} S_Earth "
              f"(buggy-units cut at 10^-9)")
        print(f"  survivor raw buggy-logS >= -9 means true S >= "
              f"{10**np.median(ls_c[surv]):.1f} S_Earth -> guard acts as hard hot/short-P selector")
    # coverage of stellar params (fallback skew)
    print(f"coverage: st_mass {(mstar==df.st_mass).mean()*100:.0f}% real (rest 0.5 fallback), "
          f"st_rad {(rstar==df.st_rad).mean()*100:.0f}%, st_teff {(tstar==df.st_teff).mean()*100:.0f}%")
    if tag=="FGK":
        mm = df.st_mass.dropna()
        print(f"  FGK st_mass median {mm.median():.2f} -> 0.5 fallback understates mass by ~{1-mm.median()/0.5:.1f}x on missing rows")
    # mission family census incl. KIC leakage
    fam = df.hostname.str.extract(r"^(TOI|Kepler|KIC|K2|HIP|HD|GJ|Wolf|LHS|TRAPPIST|EPIC)",
                                  expand=False).fillna("other")
    fammap = fam.map(lambda s: s if s in ("Kepler","K2","TOI") else "other")
    print(f"mission families: {fammap.value_counts().to_dict()}")
    kic = df.hostname.str.match(r"^KIC").sum()
    koioth = df.hostname.str.match(r"^(KOI|TIC)").sum()
    print(f"  KIC-prefix hosts: {kic} (mapped 'other', arguably Kepler-field); KOI/TIC-prefix: {koioth}")

# --- calibrate_alpha offset inputs ---
print("\n================ power-analysis calibration inputs ================")
pl = pd.read_csv(D/"planet_sample_M.csv")
k = pd.read_csv(D/"hosts_kinematics_M.csv")
base = pl.merge(k[["hostname","vtan","W"]], on="hostname", how="inner")
fam = base.hostname.str.extract(r"^(TOI|Kepler|KIC|K2|HIP|HD|GJ|Wolf|LHS|TRAPPIST|EPIC)", expand=False).fillna("other")
fm = fam.map(lambda s: s if s in ("Kepler","K2","TOI") else "other")
fk = (fm=="Kepler").mean(); f2=(fm=="K2").mean()
VALLEY=1.85
real_frac = float((base.pl_rade>VALLEY).mean())
offset = -1.594*fk + 0.893*f2
print(f"M base: {len(base)} planets, frac_Kepler={fk:.3f}, frac_K2={f2:.3f}, other={1-fk-f2:.3f}")
print(f"real_frac mean(pl_rade>1.85)={real_frac:.3f}; calib_frac=min(0.495,{real_frac:.3f}/0.90)={min(0.495,real_frac/0.90):.3f}")
print(f"calibrate_alpha omits: GAMMA_KEP*f+GAMMA_K2*f = {-1.594*fk:+.3f}{0.893*f2:+.3f} = {offset:+.3f} logits")
import math
pt = min(0.495, real_frac/0.90)
eta_t = math.log(pt/(1-pt))
pr = 1/(1+math.exp(-(eta_t+offset)))
print(f"target pre-radius p={pt:.3f} (eta={eta_t:.3f}); realized mean p≈{pr:.3f} "
      f"(odds x{math.exp(offset):.2f})")
# info change for logistic slope detection at p vs p
info = lambda p: p*(1-p)
print(f"Bernoulli info p*(1-p): target {info(pt):.3f} vs realized {info(pr):.3f} "
      f"({100*(info(pr)/info(pt)-1):+.1f}%)")
