"""Audit A2: refit real configurations under corrected GH quadrature and
corrected insolation; quantify shift in beta_age / se / sigma_u."""
import pathlib, sys, time, json
import numpy as np
import pandas as pd
from scipy.special import erf, expit, logsumexp
from scipy.optimize import minimize

ROOT = pathlib.Path("C:/Users/MIT/Downloads/Astro_research_ASGSR/mdwarf_radius_valley")
D = ROOT/"data"; R = ROOT/"results"
sys.path.insert(0, str(ROOT/"code"))
import importlib.util
spec = importlib.util.spec_from_file_location("h9", ROOT/"code"/"09_hierarchical_model.py")
m9 = importlib.util.module_from_spec(spec); spec.loader.exec_module(m9)

GH_X, GH_W = np.polynomial.hermite_e.hermegauss(24)

def gh_grid(su, fixed):
    if fixed:
        return su*GH_X, GH_W/np.sqrt(2*np.pi)
    return m9.gauss_hermite_grid(su)

def build_design_fixed(df, age_col, valley, fe=False, fix_insol=False):
    """Replica of m9.build_design with optional insolation correction."""
    d = df.copy()
    err = df.get("pl_radeerr1"); err2 = df.get("pl_radeerr2")
    e = np.where(np.isfinite(err.fillna(np.nan)), err.abs(), np.nan)
    if err2 is not None:
        e = np.where(np.isfinite(e), e, err2.abs().fillna(np.nan))
    e = np.where(np.isfinite(e) & (e > 0.01), e, 0.05*d.pl_rade)
    w_sn = 0.5*(1+erf((d.pl_rade.values-valley)/(e*np.sqrt(2))))
    mstar = np.clip(d.st_mass.fillna(0.5).values, 0.08, 1.6)
    p_yr = d.pl_orbper.values/365.25
    a_au = np.cbrt(mstar)*p_yr**(2.0/3.0)
    rstar = np.clip(d.st_rad.fillna(0.6).values, 0.1, 4)
    tstar = np.clip(d.st_teff.fillna(4000).values, 2500, 7500)
    lum = rstar**2 * ((tstar/5772)**4 if fix_insol else tstar/5772**4)
    s_in = lum/a_au**2
    lp = np.log10(d.pl_orbper.values); ls = np.log10(s_in)
    fam = d.hostname.str.extract(r"^(TOI|Kepler|KIC|K2|HIP|HD|GJ|Wolf|LHS|TRAPPIST|EPIC)",
                                 expand=False).fillna("other")
    fam = fam.map(lambda s: s if s in ("Kepler","K2","TOI") else "other")
    def z(v):
        v=np.asarray(v,float); return (v-v.mean())/(v.std()+1e-12)
    age_z = z(d[age_col].values)
    iqr = d[age_col].quantile(.75)-d[age_col].quantile(.25)
    X=[np.ones(len(d))]; names=["alpha"]
    X.append(age_z); names.append("beta_age")
    for k,v in (("logP",lp),("logS",ls)):
        vz=z(v); vz[(v<-9)|(v>9)]=0.0
        X.append(vz); names.append(f"gamma_{k}")
    X.append((fam=="Kepler").astype(float).values); names.append("gamma_missionKepler")
    X.append((fam=="K2").astype(float).values); names.append("gamma_missionK2")
    if fe:
        X.append(z(d.st_met.fillna(0.0).values)); names.append("gamma_FeH")
    Xm=np.column_stack(X)
    hosts=d.hostname.values
    uh=pd.unique(hosts)
    hidx=pd.Categorical(hosts,categories=uh).codes
    return {"X":Xm,"names":names,"w":w_sn,"hidx":hidx,"n_hosts":len(uh),
            "age_iqr":iqr,"age_raw":d[age_col].values}

def neg_log_post(theta, mdl, fixed_gh):
    beta=theta[:-1]; log_su=theta[-1]
    if not (-6.0<log_su<2.0): return 1e12
    su=np.exp(log_su)
    eta0=np.clip(mdl["X"]@beta,-30,30)
    w=np.clip(mdl["w"],1e-12,1-1e-12)
    lw,lw1=np.log(w),np.log1p(-w)
    hidx=mdl["hidx"]
    u,gw=gh_grid(su,fixed_gh)
    totals=np.empty((len(u),mdl["n_hosts"]))
    for j,uj in enumerate(u):
        pi=expit(np.clip(eta0+uj,-30,30))
        li=np.logaddexp(lw+np.log(pi),lw1+np.log1p(-pi))
        totals[j]=np.bincount(hidx,weights=li,minlength=mdl["n_hosts"])
    nll=-np.sum(logsumexp(totals+np.log(gw)[:,None],axis=0))
    nll+=0.5*np.sum((beta[1:]/2.5)**2)+0.5*(beta[0]/5.0)**2+0.5*((log_su+1.5)/1.0)**2
    return nll if np.isfinite(nll) else 1e12

def fit(mdl, fixed_gh):
    p=mdl["X"].shape[1]
    th0=np.r_[np.zeros(p),[-1.5]]
    res=minimize(neg_log_post,th0,args=(mdl,fixed_gh),method="L-BFGS-B",
                 bounds=[(-25,25)]*p+[(-6,2)],options={"maxiter":2000})
    f0=neg_log_post(res.x,mdl,fixed_gh)
    n=len(th0); eps=1e-4; H=np.zeros((n,n))
    for i in range(n):
        ei=np.zeros(n); ei[i]=eps
        H[i,i]=(neg_log_post(res.x+ei,mdl,fixed_gh)-2*f0
                +neg_log_post(res.x-ei,mdl,fixed_gh))/eps**2
        for j in range(i+1,n):
            ej=np.zeros(n); ej[j]=eps
            Hij=(neg_log_post(res.x+ei+ej,mdl,fixed_gh)-neg_log_post(res.x+ei-ej,mdl,fixed_gh)
                 -neg_log_post(res.x-ei+ej,mdl,fixed_gh)+neg_log_post(res.x-ei-ej,mdl,fixed_gh))/(4*eps*eps)
            H[i,j]=H[j,i]=Hij
    cov=np.linalg.inv(H+np.eye(n)*1e-6)
    return res.x,np.sqrt(np.diag(cov))

def load(tag,kin_col="vtan"):
    pl=pd.read_csv(D/f"planet_sample_{tag}.csv")
    k=pd.read_csv(D/f"hosts_kinematics_{tag}.csv")
    df=pl.merge(k[["hostname","vtan","W"]],on="hostname",how="inner")
    if kin_col=="absW":
        df=df.dropna(subset=["W"]).copy(); df["absW"]=df.W.abs()
    return df

configs=[
    ("M","vtan","vtan",1.85,False,None,"M_vtan"),
    ("M","absW","absW",1.85,False,None,"M_absW"),
    ("FGK","vtan","vtan",1.88,False,None,"FGK_vtan"),
]
print(f"{'fit':10s} {'variant':16s} {'beta/IQR':>18s} {'z':>6s} {'sig_u(rep)':>10s} {'eff mix SD':>10s}")
for tag,kin_col,acol,valley,fe,dmax,label in configs:
    df=load(tag,kin_col)
    for variant,fix_i,fix_g in [("as-published",False,False),("GH-fixed",False,True),
                                ("insol-fixed",True,False),("both-fixed",True,True)]:
        mdl=build_design_fixed(df,acol,valley,fe=fe,fix_insol=fix_i)
        t0=time.time()
        th,se=fit(mdl,fix_g)
        iqr=mdl["age_iqr"]; sd=mdl["age_raw"].std()
        b=th[1]*(iqr/sd); s=se[1]*(iqr/sd)
        print(f"{label:10s} {variant:16s} {b:+.3f}+-{s:.3f} {b/s:+6.2f} "
              f"{np.exp(th[-1]):10.3f} {np.exp(th[-1])*np.sqrt(2):10.3f}   [{time.time()-t0:.0f}s]")
