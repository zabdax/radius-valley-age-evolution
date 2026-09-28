"""Audit A4: extended selftest for 09 (bias/coverage with more seeds,
realistic host-effect SD) + Laplace Hessian eps stability."""
import pathlib, sys
import numpy as np
import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT/"code"))
import importlib.util
spec = importlib.util.spec_from_file_location("h9", ROOT/"code"/"09_hierarchical_model.py")
m9 = importlib.util.module_from_spec(spec); spec.loader.exec_module(m9)
from scipy.special import expit
from scipy.optimize import minimize

def one_fit(seed, b_true, host_sd):
    rng = np.random.default_rng(100+seed)
    n_host = 600
    hosts = np.repeat([f"H{i}" for i in range(n_host)], rng.integers(1,3,n_host))
    n=len(hosts)
    age=rng.uniform(10,80,n); logP=rng.uniform(-0.5,1.6,n); logS=rng.uniform(-1,2,n)
    kep=(rng.uniform(size=n)<0.7).astype(float)
    eta=(-0.3 + b_true*(age-age.mean())/age.std() + 0.3*logP - 0.2*logS + 0.15*kep
         + rng.normal(0,host_sd,n))
    pi=expit(eta); sn=rng.uniform(size=n)<pi
    r_true=np.where(sn, rng.normal(2.4,0.5,n), rng.normal(1.3,0.25,n))
    r_obs=r_true+rng.normal(0,0.06,n)
    df=pd.DataFrame({"hostname":hosts,"pl_rade":r_obs,"pl_radeerr1":0.06,
                     "pl_orbper":10**logP,"st_mass":1.0,"st_rad":1.0,"st_teff":5700,
                     "sy_dist":300,"st_met":0.0,"vtan":age})
    mdl=m9.build_design(df,"vtan",valley=1.88,fe=False)
    theta,se,cov,_=m9.fit(mdl)
    sd=age.std(); iqr=np.percentile(age,75)-np.percentile(age,25)
    return theta[1]*(iqr/sd), se[1]*(iqr/sd), b_true*(iqr/sd)

def block(name, seeds, b_true, host_sd):
    ests, hits = [], 0
    for s in seeds:
        e_,se_,tr = one_fit(s,b_true,host_sd)
        ests.append(e_)
        hits += (e_-1.96*se_ <= tr) and (tr <= e_+1.96*se_)
    ests=np.array(ests)
    print(f"{name}: n={len(ests)} mean_est={ests.mean():+.3f} sd={ests.std():.3f} "
          f"truth={b_true*(np.percentile(np.arange(10,81),75)-np.percentile(np.arange(10,81),25))/np.arange(10,81).std():+.3f} "
          f"mean_bias={ests.mean()-b_true*1.731:+.3f} coverage={hits}/{len(ests)}")

print("== extended SELFTEST (as-published quadrature/insolation) ==")
block("A: b=-0.8, host_sd=0.35 (paper regime)", range(24), -0.8, 0.35)
block("B: b=-0.8, host_sd=1.5  (real-data regime)", range(24,48), -0.8, 1.5)
block("C: b= 0.0, host_sd=1.5  (null / FPR)",      range(48,72),  0.0, 1.5)

# --- Laplace Hessian eps stability on the real M x vtan fit ---
print("\n== Laplace Hessian eps stability (M x vtan) ==")
D = ROOT/"data"
pl=pd.read_csv(D/"planet_sample_M.csv"); k=pd.read_csv(D/"hosts_kinematics_M.csv")
df=pl.merge(k[["hostname","vtan","W"]],on="hostname",how="inner")
mdl=m9.build_design(df,"vtan",1.85,fe=False)

def hess_cov(mdl, th0v, eps):
    from scipy.optimize import minimize
    res=minimize(m9.neg_log_post,np.r_[np.zeros(mdl["X"].shape[1]),[-1.5]],
                 args=(mdl,),method="L-BFGS-B",
                 bounds=[(-25,25)]*(mdl["X"].shape[1])+[(-6,2)],options={"maxiter":2000})
    x=res.x; f0=m9.neg_log_post(x,mdl); n=len(x); H=np.zeros((n,n))
    for i in range(n):
        ei=np.zeros(n); ei[i]=eps
        H[i,i]=(m9.neg_log_post(x+ei,mdl)-2*f0+m9.neg_log_post(x-ei,mdl))/eps**2
        for j in range(i+1,n):
            ej=np.zeros(n); ej[j]=eps
            Hij=(m9.neg_log_post(x+ei+ej,mdl)-m9.neg_log_post(x+ei-ej,mdl)
                 -m9.neg_log_post(x-ei+ej,mdl)+m9.neg_log_post(x-ei-ej,mdl))/(4*eps*eps)
            H[i,j]=H[j,i]=Hij
    cov=np.linalg.inv(H+np.eye(n)*1e-6)
    return np.sqrt(np.diag(cov)), np.linalg.eigvalsh(H).min()

for eps in (1e-3,1e-4,1e-5):
    se,min_eig = hess_cov(mdl,None,eps)
    print(f"eps={eps:.0e}: se(beta_age)={se[1]:.4f} se(alpha)={se[0]:.4f} "
          f"se(log_su)={se[-1]:.4f} min_eig(H)={min_eig:.2f}")
