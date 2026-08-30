"""audit_v3.py — Evidence-first audit of the high-replicate Monte Carlo.
Raw CSVs = ground truth; JSONs = cross-checks. Old local results loaded
from the project's results/ directory for exact old-vs-new comparison."""
import json
import hashlib
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

BASE = Path(r"C:\Users\MIT\Downloads\Astro_research_ASGSR\mdwarf_radius_valley")
Z1 = BASE / "results" / "colab_audit" / "zip1"
Z2 = BASE / "results" / "colab_audit" / "zip2" / "mdwarf_radius_valley_highrep_fast_v3"
RES = BASE / "results"

def cp_ci(x, n, a=0.05):
    lo = stats.beta.ppf(a / 2, x, n - x + 1) if x > 0 else 0.0
    hi = stats.beta.ppf(1 - a / 2, x + 1, n - x) if x < n else 1.0
    return lo, hi

def wald_pct_ci(p_hat, n):
    return cp_ci(round(p_hat * n), n)

def boot_ci(v, fn=np.median, B=2000, seed=7):
    rng = np.random.default_rng(seed)
    v = np.asarray(v, float)
    out = np.empty(B)
    for i in range(B):
        out[i] = fn(v[rng.integers(0, len(v), len(v))])
    return float(np.quantile(out, .025)), float(np.quantile(out, .975))

def md(x): return float(np.median(x))
def q25(x): return float(np.quantile(x, .025))
def q97(x): return float(np.quantile(x, .975))

print("=" * 78)
print("SECTION 0 - FILE INTEGRITY: zip1 vs zip2 hashes, manifests")
print("=" * 78)
files = ["task2_calibration_kepler_only_trim0.95_rep1000.json",
         "task2_gate_rep1000_replicates_raw.csv",
         "power_analysis_rep500.json",
         "power_grid_rep500_replicates_raw.csv", "manifest.json"]
for f in files:
    h1 = hashlib.sha256((Z1 / f).read_bytes()).hexdigest()[:16]
    h2 = hashlib.sha256((Z2 / f).read_bytes()).hexdigest()[:16]
    print(f"{f}: zip1={h1} zip2={h2} {'MATCH' if h1 == h2 else 'MISMATCH!'}")
man = json.loads((Z1 / "manifest.json").read_text())
acc = json.loads((Z2 / "accelerator_manifest.json").read_text())
print("manifest:", json.dumps(man))
print("accelerator:", json.dumps(acc))

print()
print("=" * 78)
print("SECTION 1 - DATA INTEGRITY (raw CSVs = ground truth)")
print("=" * 78)
gate = pd.read_csv(Z1 / "task2_gate_rep1000_replicates_raw.csv",
                   keep_default_na=False)
grid = pd.read_csv(Z1 / "power_grid_rep500_replicates_raw.csv",
                   keep_default_na=False)
t2 = json.loads((Z1 / "task2_calibration_kepler_only_trim0.95_rep1000.json").read_text())
pj = json.loads((Z1 / "power_analysis_rep500.json").read_text())

print(f"gate CSV: {len(gate)} rows, columns={list(gate.columns)}")
for arm in ("null", "signal"):
    g = gate[gate.arm == arm]
    reps = g.rep.values.astype(int)
    b = g.beta_hat_per_iqr.values.astype(float)
    print(f"[gate {arm}] N={len(g)} | contiguous 0..{len(g)-1}: "
          f"{bool(np.array_equal(np.sort(reps), np.arange(len(g))))} | "
          f"dup reps: {len(reps) - len(np.unique(reps))} | "
          f"NaN/Inf beta: {int(np.sum(~np.isfinite(b)))} | "
          f"alpha unique: {g.alpha.nunique()} ({g.alpha.iloc[0]:+.6f}) | "
          f"beta_inj unique: {g.beta_inj.nunique()} ({g.beta_inj.iloc[0]}) | "
          f"n_planets mean/sd/min/max = {g.n_planets.mean():.1f}/"
          f"{g.n_planets.std():.1f}/{g.n_planets.min()}/{g.n_planets.max()}")

SCALES = ["0.5x", "1x", "2x", "4x", "8x"]
EFFECTS = ["null", "SWEET", "strong"]
SCALE_NUM = {"0.5x": 0.5, "1x": 1, "2x": 2, "4x": 4, "8x": 8}
print(f"grid CSV: {len(grid)} rows, columns={list(grid.columns)}")
cell_rows = {}
for s in SCALES:
    for e in EFFECTS:
        c = grid[grid.cell == f"{s}|{e}"]
        reps = c.rep.values.astype(int)
        b = c.beta_hat_per_sd.values.astype(float)
        se = c.se_per_sd.values.astype(float)
        d, v = c.detect.values.astype(int), c.cover.values.astype(int)
        ok_contig = bool(np.array_equal(np.sort(reps), np.arange(len(c))))
        dup = len(reps) - len(np.unique(reps))
        nan_bad = int(np.sum(~np.isfinite(b)) + np.sum(~np.isfinite(se)))
        flag_ok = bool(np.all(np.isin(d, [0, 1])) and np.all(np.isin(v, [0, 1])))
        cell_rows[(s, e)] = c
        print(f"[grid {s:>4}|{e:>6}] N={len(c):3d} (target 500) | contiguous="
              f"{ok_contig} | dups={dup} | NaN/Inf={nan_bad} | "
              f"flags_in_0/1={flag_ok} | n_planets mean={c.n_planets.mean():6.1f} "
              f"(scale*417={SCALE_NUM[s] * 417:.0f}) sd={c.n_planets.std():5.1f}")

print()
print("=" * 78)
print("SECTION 2 - GATE: independent recomputation from RAW (n=1000/arm)")
print("=" * 78)
b_real, se_real = t2["beta_real"], t2["se_real"]
print(f"real beta (JSON): {b_real:+.6f} +- {se_real:.6f} per IQR")
print(f"[cross-check] local original real fit: -0.0295958711 +- 0.1061026890 "
      f"-> delta = {abs(b_real - (-0.029595871147796216)):.2e} / "
      f"{abs(se_real - 0.10610268898228108):.2e} (fit-tolerance level)")
gate_stats = {}
for arm in ("null", "signal"):
    b = gate[gate.arm == arm].sort_values("rep").beta_hat_per_iqr.values.astype(float)
    pct = float(np.mean(b <= b_real))
    lo, hi = wald_pct_ci(pct, len(b))
    mlo, mhi = boot_ci(b, md)
    s = {"n": len(b), "median": md(b), "q025": q25(b), "q975": q97(b),
         "mean": float(np.mean(b)), "sd": float(np.std(b, ddof=1)),
         "mad": float(np.median(np.abs(b - md(b)))), "min": float(b.min()),
         "max": float(b.max()), "pct_rank": pct, "pct_ci": (lo, hi),
         "median_boot_ci": (mlo, mhi)}
    gate_stats[arm] = s
    print(f"[{arm}] n={s['n']} median={s['median']:+.4f} "
          f"(boot95 [{mlo:+.4f},{mhi:+.4f}]) q025={s['q025']:+.4f} "
          f"q975={s['q975']:+.4f} mean={s['mean']:+.4f} sd={s['sd']:.4f} "
          f"mad={s['mad']:.4f} min={s['min']:+.3f} max={s['max']:+.3f}")
    if arm == "null":
        print(f"      P(null >= real) = {1 - pct:.4f} | pct_rank={pct:.3f} "
              f"(binom95 [{lo:.3f},{hi:.3f}])")
    else:
        print(f"      pct_rank of real in signal = {pct:.3f} "
              f"(binom95 [{lo:.3f},{hi:.3f}])")

print("JSON cross-check (raw-recomputed vs summary JSON):")
for arm in ("null", "signal"):
    j = t2[arm]
    r = gate_stats[arm]
    print(f"  {arm}: median d={abs(r['median'] - j['median']):.2e} "
          f"q025 d={abs(r['q025'] - j['q025']):.2e} "
          f"q975 d={abs(r['q975'] - j['q975']):.2e} "
          f"pct d={abs(r['pct_rank'] - j['pct_rank_of_real']):.2e}")

p_null = 1.0 - gate_stats["null"]["pct_rank"]
inside_sig = gate_stats["signal"]["q025"] <= b_real <= gate_stats["signal"]["q975"]
dist_null = p_null <= 0.05
verdict = "PASS" if (dist_null and inside_sig) else (
    "FAIL (inconclusive)" if not dist_null else "FAIL (null-like)")
print(f"GATE RULE: P(null>=real)={p_null:.4f} -> distinguishable={dist_null}; "
      f"inside signal 95% envelope={inside_sig} -> VERDICT: {verdict}")

print()
print("--- STREAM FINGERPRINT: new first-120 subset vs original 120-run ---")
OLD = dict(null_med=-0.013343129808420381, null_q025=-0.24811345476351288,
           null_q975=0.2213773787135324, null_pct=0.43333333333333335,
           sig_med=-0.17838442092332896, sig_q025=-0.355084009121274,
           sig_q975=0.06806684954795611, sig_pct=0.9083333333333333)
KEYMAP = {"null": ("null_med", "null_q025", "null_q975", "null_pct"),
          "signal": ("sig_med", "sig_q025", "sig_q975", "sig_pct")}
for arm in ("null", "signal"):
    km = KEYMAP[arm]
    b = gate[gate.arm == arm].sort_values("rep")
    sub = b[b.rep < 120].beta_hat_per_iqr.values.astype(float)
    print(f"[{arm} first-120] median={md(sub):+.6f} (old {OLD[km[0]]:+.6f}, "
          f"d={abs(md(sub) - OLD[km[0]]):.2e}) "
          f"q025={q25(sub):+.6f} (old {OLD[km[1]]:+.6f}, "
          f"d={abs(q25(sub) - OLD[km[1]]):.2e}) "
          f"q975={q97(sub):+.6f} (old {OLD[km[2]]:+.6f}, "
          f"d={abs(q97(sub) - OLD[km[2]]):.2e}) "
          f"pct={np.mean(sub <= b_real):.4f} (old {OLD[km[3]]:.4f})")
print("--- fit-level fingerprint vs local smoke run (same seeds/draws) ---")
SMOKE_GATE = {0: -0.12590420884412232, 1: -0.17669393238519407}
b0 = gate[gate.arm == "null"].sort_values("rep")
for r, v in SMOKE_GATE.items():
    newv = float(b0[b0.rep == r].beta_hat_per_iqr.iloc[0])
    print(f"gate null rep{r}: smoke={v:+.9f} new={newv:+.9f} delta={abs(newv - v):.2e}")

print()
print("=" * 78)
print("SECTION 3 - GRID: per-cell recomputation from RAW + JSON cross-check")
print("=" * 78)
OLD_GRID = {
    ("0.5x", "null"): (120, 0.06666666666666667, 0.9333333333333333),
    ("0.5x", "SWEET"): (150, 0.06, 0.9466666666666667),
    ("0.5x", "strong"): (150, 0.16, 0.9133333333333333),
    ("1x", "null"): (120, 0.11666666666666667, 0.8833333333333333),
    ("1x", "SWEET"): (150, 0.09333333333333334, 0.8866666666666667),
    ("1x", "strong"): (150, 0.3, 0.8933333333333333),
    ("2x", "null"): (120, 0.041666666666666664, 0.9583333333333334),
    ("2x", "SWEET"): (150, 0.14666666666666667, 0.9133333333333333),
    ("2x", "strong"): (150, 0.43333333333333335, 0.9),
    ("4x", "null"): (120, 0.025, 0.975),
    ("4x", "SWEET"): (150, 0.14666666666666667, 0.8733333333333333),
    ("4x", "strong"): (150, 0.46, 0.8266666666666667),
    ("8x", "null"): (120, 0.075, 0.925),
    ("8x", "SWEET"): (150, 0.18666666666666668, 0.94),
    ("8x", "strong"): (150, 0.64, 0.88),
}
grid_stats = {}
for s in SCALES:
    for e in EFFECTS:
        c = cell_rows[(s, e)]
        n = len(c)
        d = c.detect.values.astype(int)
        v = c.cover.values.astype(int)
        b = c.beta_hat_per_sd.values.astype(float)
        se = c.se_per_sd.values.astype(float)
        dr, cr = float(d.mean()), float(v.mean())
        dlo, dhi = cp_ci(int(d.sum()), n)
        clo, chi = cp_ci(int(v.sum()), n)
        jrow = [r for r in pj if r["scale"] == s and r["effect"] == e][0]
        K, old_dr, old_cov = OLD_GRID[(s, e)]
        sub = c[c.rep < K]
        sub_dr = float(sub.detect.mean())
        grid_stats[(s, e)] = dict(
            n=n, detect=dr, dci=(dlo, dhi), cover=cr, cci=(clo, chi),
            med_b=md(b), med_se=md(se), json_dr=jrow["detect_rate"],
            json_cov=jrow["coverage"], old_n=K, old_dr=old_dr,
            old_cov=old_cov, sub_dr=sub_dr, sub_n=int(len(sub)))
        print(f"[{s:>4}|{e:>6}] n={n} detect={dr:.3f} CP[{dlo:.3f},{dhi:.3f}] "
              f"cover={cr:.3f} CP[{clo:.3f},{chi:.3f}] med_b={md(b):+.3f} "
              f"med_se={md(se):.3f} | JSON d(dr)={abs(dr - jrow['detect_rate']):.1e} "
              f"d(cov)={abs(cr - jrow['coverage']):.1e} | first-{K}: "
              f"detect {sub_dr:.4f} vs old {old_dr:.4f} (n={len(sub)})")

print()
print("=" * 78)
print("SECTION 4 - NULL FPR AUDIT")
print("=" * 78)
for s in SCALES:
    st = grid_stats[(s, "null")]
    x = round(st["detect"] * st["n"])
    lo, hi = st["dci"]
    pbin = stats.binom.sf(x - 1, st["n"], .05)
    print(f"{s:>4}: FPR={st['detect']:.4f} ({x}/{st['n']}) CP95[{lo:.4f},{hi:.4f}] "
          f"vs nominal 0.05: {'consistent' if lo <= 0.05 <= hi else 'OFF'} "
          f"(binom P(X>={x})={pbin:.4f})")
n8x = round(grid_stats[("8x", "null")]["detect"] * grid_stats[("8x", "null")]["n"])
N8 = grid_stats[("8x", "null")]["n"]
rest_x = sum(round(grid_stats[(s, "null")]["detect"] * grid_stats[(s, "null")]["n"])
             for s in SCALES if s != "8x")
rest_n = sum(grid_stats[(s, "null")]["n"] for s in SCALES if s != "8x")
for s in SCALES:
    if s == "8x":
        continue
    xi = round(grid_stats[(s, "null")]["detect"] * grid_stats[(s, "null")]["n"])
    ni = grid_stats[(s, "null")]["n"]
    p = stats.fisher_exact([[n8x, N8 - n8x], [xi, ni - xi]])[1]
    print(f"Fisher 8x({n8x}/{N8}) vs {s}({xi}/{ni}): p={p:.4f}")
p_pool = stats.fisher_exact([[n8x, N8 - n8x], [rest_x, rest_n - rest_x]])[1]
print(f"Fisher 8x vs pooled others ({rest_x}/{rest_n}): p={p_pool:.4f}")
tot_x, tot_n = n8x + rest_x, N8 + rest_n
print(f"pooled all null cells: {tot_x}/{tot_n} = {tot_x / tot_n:.4f}; "
      f"binom vs 5%: P(X>={tot_x}|{tot_n},.05)={stats.binom.sf(tot_x - 1, tot_n, .05):.4f}")

print()
print("=" * 78)
print("SECTION 5 - COVERAGE AUDIT")
print("=" * 78)
tot_cov = tot_cn = 0
for s in SCALES:
    for e in EFFECTS:
        st = grid_stats[(s, e)]
        x = round(st["cover"] * st["n"])
        tot_cov += x
        tot_cn += st["n"]
        tag = ""
        if st["cci"][1] < 0.95:
            tag = "ANTI-CONSERVATIVE (upper CI < 0.95)"
        elif st["cci"][0] > 0.95:
            tag = "conservative (lower CI > 0.95)"
        print(f"[{s:>4}|{e:>6}] coverage={st['cover']:.3f} ({x}/{st['n']}) "
              f"CP95[{st['cci'][0]:.3f},{st['cci'][1]:.3f}] {tag}")
print(f"pooled coverage all 15 cells: {tot_cov}/{tot_cn} = {tot_cov / tot_cn:.4f} "
      f"(binom P(X<={tot_cov}|{tot_cn},.95)={stats.binom.cdf(tot_cov, tot_cn, .95):.4f})")

print()
print("=" * 78)
print("SECTION 6 - POWER TABLE (raw-derived)")
print("=" * 78)
for s in SCALES:
    npl = cell_rows[(s, "null")].n_planets.mean()
    parts = [f"{s:>4} N~{npl:6.0f}"]
    for e in EFFECTS:
        st = grid_stats[(s, e)]
        parts.append(f"{e}: {st['detect']:.3f} [{st['dci'][0]:.3f},{st['dci'][1]:.3f}] n={st['n']}")
    print(" | ".join(parts))
print("SWEET curve:  ", [round(grid_stats[(s, 'SWEET')]['detect'], 3) for s in SCALES])
print("strong curve: ", [round(grid_stats[(s, 'strong')]['detect'], 3) for s in SCALES])
print("null FPR:     ", [round(grid_stats[(s, 'null')]['detect'], 3) for s in SCALES])
ses = [grid_stats[(s, 'SWEET')]['med_se'] for s in SCALES]
print("median se (SWEET):", [f"{x:.3f}" for x in ses],
      "| ratios vs 8x:", [f"{x / ses[-1]:.2f}" for x in ses],
      "(sqrt-scale expect:", [f"{np.sqrt(8 / sc):.2f}" for sc in (0.5, 1, 2, 4, 8)], ")")

print()
print("=" * 78)
print("SECTION 7 - EXTRAPOLATION CHECK (Gaussian N^-1/2, model-dependent)")
print("=" * 78)
p8 = grid_stats[("8x", "SWEET")]["detect"]
plo, phi = grid_stats[("8x", "SWEET")]["dci"]
for p, tag in [(plo, "CP-lo"), (p8, "point"), (phi, "CP-hi")]:
    x = stats.norm.ppf(p) + 1.96
    f = ((stats.norm.ppf(0.8) + 1.96) / x) ** 2
    print(f"p8={p:.3f} ({tag}): |b|/s={x:.3f} factor={f:.2f} "
          f"N80={3336 * f:,.0f} (~{3336 * f / 417:.0f}x current)")
xo = stats.norm.ppf(0.18666666666666668) + 1.96
fo = ((stats.norm.ppf(0.8) + 1.96) / xo) ** 2
print(f"OLD p8=0.1867 (150 reps): x={xo:.3f} factor={fo:.2f} N80={3336 * fo:,.0f}")

print()
print("=" * 78)
print("SECTION 8 - NUMERICAL / OPTIMIZER DIAGNOSTICS (raw-level)")
print("=" * 78)
for s in SCALES:
    for e in EFFECTS:
        c = cell_rows[(s, e)]
        b = c.beta_hat_per_sd.values.astype(float)
        se = c.se_per_sd.values.astype(float)
        dupb = len(b) - len(np.unique(np.round(b, 12)))
        print(f"[{s:>4}|{e:>6}] |b|>0.5:{int(np.sum(np.abs(b) > 0.5)):3d} "
              f"|b|>1:{int(np.sum(np.abs(b) > 1)):2d} se<0.02:{int(np.sum(se < 0.02)):2d} "
              f"se>1:{int(np.sum(se > 1)):2d} dup_b:{dupb:3d} "
              f"b[{b.min():+.3f},{b.max():+.3f}] se[{se.min():.3f},{se.max():.3f}]")
for arm in ("null", "signal"):
    bb = gate[gate.arm == arm].beta_hat_per_iqr.values.astype(float)
    dupg = len(bb) - len(np.unique(np.round(bb, 12)))
    print(f"gate {arm}: |b|>0.5:{int(np.sum(np.abs(bb) > 0.5))} "
          f"|b|>1:{int(np.sum(np.abs(bb) > 1))} dup_b:{dupg} "
          f"b[{bb.min():+.3f},{bb.max():+.3f}]")

print()
print("=" * 78)
print("SECTION 9 - OLD vs NEW")
print("=" * 78)
for nm, o, n in [("null median", OLD["null_med"], gate_stats["null"]["median"]),
                 ("null q025", OLD["null_q025"], gate_stats["null"]["q025"]),
                 ("null q975", OLD["null_q975"], gate_stats["null"]["q975"]),
                 ("real pct in null", OLD["null_pct"], gate_stats["null"]["pct_rank"]),
                 ("signal median", OLD["sig_med"], gate_stats["signal"]["median"]),
                 ("signal q025", OLD["sig_q025"], gate_stats["signal"]["q025"]),
                 ("signal q975", OLD["sig_q975"], gate_stats["signal"]["q975"]),
                 ("real pct in signal", OLD["sig_pct"], gate_stats["signal"]["pct_rank"])]:
    print(f"{nm:22s} old={o:+.4f} new={n:+.4f} diff={n - o:+.4f}")
for s in SCALES:
    for e in EFFECTS:
        st = grid_stats[(s, e)]
        print(f"{s} {e:>6}: detect {st['old_dr']:.3f}({st['old_n']}) -> "
              f"{st['detect']:.3f}({st['n']}) | cover {st['old_cov']:.3f} -> {st['cover']:.3f}")

summary = {"gate": gate_stats, "verdict": verdict,
           "grid": {f"{s}|{e}": grid_stats[(s, e)] for s in SCALES for e in EFFECTS},
           "manifest": man, "accelerator": acc}
(RES / "colab_audit" / "audit_summary.json").write_text(
    json.dumps(summary, indent=1, default=str))
print("\nsaved audit_summary.json")