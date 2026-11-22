#!/usr/bin/env python3
"""verify_archived_runs_v2.py -- computational verification of two ARCHIVED
Monte-Carlo runs without rerunning them.

Item 2: confirm_independent_production.py output (results/appendix/confirm_*)
Item 3: appendix_B_beta_dist.py output (results/appendix/beta_dist_*)

Strategy:
  * Re-derive every summary statistic from the archived per-replicate RAW CSVs
    and compare it to the archived summary JSON (exact, tight tolerance).
  * Re-generate the beta_injection draws from the exact RNG recipe/seed in
    appendix_B_beta_dist.py and compare them to the archived raw beta_inj column.
  * Confirm run scale (reps/arm), mode == 'independent', and arm separations.
All reads are from the frozen project; all writes go to this audit directory.
"""
import json
import pathlib

import numpy as np
import pandas as pd

ROOT = pathlib.Path(r"REPO")
FROZEN = ROOT / "radius-valley-age-evolution"
APP = FROZEN / "results" / "appendix"
HIGHREP = FROZEN / "results" / "highrep"
OUT = ROOT / "manuscripts" / "q1-paper" / "v2-2026-11" / "results" / "audit"

REAL = {"M": 0.0458, "FGK": -0.0296}


def close(a, b, tol=1e-9):
    try:
        return bool(np.isclose(a, b, rtol=tol, atol=tol))
    except (TypeError, ValueError):
        return bool(a == b)


def arm_stats(b):
    b = np.asarray(b, float)
    b = b[np.isfinite(b)]
    return {
        "n": int(len(b)),
        "median": float(np.median(b)),
        "sd": float(np.std(b, ddof=1)),
        "q025": float(np.quantile(b, 0.025)),
        "q975": float(np.quantile(b, 0.975)),
        "pct_rank_of_real": float(np.mean(b <= np.nan)),  # replaced by caller
    }


def check_confirm():
    summ = json.loads((APP / "confirm_summary.json").read_text())
    # NOTE: the archived CSV writes the literal string "null" in the arm column;
    # pandas' default NA parsing turns that into NaN. Disable default NA parsing
    # and coerce numeric columns explicitly so the "null" arm is preserved.
    raw = pd.read_csv(APP / "confirm_independent_all_raw.csv", keep_default_na=False)
    for col in ("rep", "beta_hat", "se", "n", "beta_inj"):
        raw[col] = pd.to_numeric(raw[col], errors="coerce")
    checks = {}
    for sample in ("M", "FGK"):
        sub = raw[raw["sample"] == sample]
        arms = {}
        for arm in ("null", "signal"):
            b = sub[sub["arm"] == arm]["beta_hat"].values
            b = b[np.isfinite(b)]
            arms[arm] = {
                "n": int(len(b)),
                "median": float(np.median(b)),
                "sd": float(np.std(b, ddof=1)),
                "q025": float(np.quantile(b, 0.025)),
                "q975": float(np.quantile(b, 0.975)),
                "pct_rank_of_real": float(np.mean(b <= REAL[sample])),
            }
        sep = (arms["signal"]["median"] - arms["null"]["median"]) / arms["null"]["sd"]
        stored = summ[sample]
        field_checks = {}
        for arm in ("null", "signal"):
            for k in ("n", "median", "sd", "q025", "q975", "pct_rank_of_real"):
                field_checks[f"{arm}.{k}"] = {
                    "archived": stored["arms"][arm][k],
                    "recomputed": arms[arm][k],
                    "match": close(stored["arms"][arm][k], arms[arm][k]),
                }
        field_checks["arm_separation_nullSD"] = {
            "archived": stored["arm_separation_nullSD"],
            "recomputed": float(sep),
            "match": close(stored["arm_separation_nullSD"], sep),
        }
        checks[sample] = field_checks

    meta = summ["meta"]
    modes = sorted(raw["mode"].unique().tolist())
    per_arm_n = {s: {a: int(((raw["sample"] == s) & (raw["arm"] == a)).sum())
                     for a in ("null", "signal")} for s in ("M", "FGK")}
    return {
        "archived_summary_path": str(APP / "confirm_summary.json"),
        "archived_all_raw_path": str(APP / "confirm_independent_all_raw.csv"),
        "archived_per_arm_raw_paths": [str(APP / f"confirm_independent_{s}_{a}_raw.csv")
                                       for s in ("M", "FGK") for a in ("null", "signal")],
        "script_path": str(FROZEN / "code" / "confirm_independent_production.py"),
        "model_path": str(FROZEN / "code" / "appendix_A_arm_delta.py"),
        "synth_path": str(FROZEN / "code" / "15_stratified_gate.py"),
        "meta": meta,
        "raw_total_rows": int(len(raw)),
        "raw_modes": modes,
        "raw_per_arm_n": per_arm_n,
        "recomputed_vs_archived": checks,
        "all_fields_match": all(
            c["match"] for s in checks.values() for c in s.values()),
        "scale_matches_claim": (per_arm_n["M"] == {"null": 1000, "signal": 1000}
                                and per_arm_n["FGK"] == {"null": 400, "signal": 400}),
        "mode_is_independent": modes == ["independent"],
        "separation_M_approx_minus0.26": close(checks["M"]["arm_separation_nullSD"]["recomputed"], -0.256, 0.005),
        "separation_FGK_approx_minus0.83": close(checks["FGK"]["arm_separation_nullSD"]["recomputed"], -0.835, 0.01),
    }


# ---- Exact copy of the injection recipe in appendix_B_beta_dist.py -----------
R_SD = 0.10


def draw_beta_inj(rng, n):
    r1 = rng.normal(0.51, R_SD, n)
    r2 = rng.normal(0.64, R_SD, n)
    bad = (r1 <= 0.05) | (r2 <= 0.05) | (~np.isfinite(r1)) | (~np.isfinite(r2))
    while bad.any():
        r1[bad] = rng.normal(0.51, R_SD, bad.sum())
        r2[bad] = rng.normal(0.64, R_SD, bad.sum())
        bad = ((r1 <= 0.05) | (r2 <= 0.05) | (~np.isfinite(r1))
               | (~np.isfinite(r2)))
    f1 = 1 / (1 + r1)
    f2 = 1 / (1 + r2)
    logit = lambda f: np.log(f / (1 - f))
    return -(logit(f1) - logit(f2)) / 1.6


def check_betadist():
    summ = json.loads((APP / "beta_dist_summary.json").read_text())
    raw = pd.read_csv(APP / "beta_dist_M_all_raw.csv")
    out = {
        "archived_summary_path": str(APP / "beta_dist_summary.json"),
        "archived_all_raw_path": str(APP / "beta_dist_M_all_raw.csv"),
        "archived_legacy_raw_path": str(APP / "beta_dist_M_legacy_raw.csv"),
        "archived_indep_raw_path": str(APP / "beta_dist_M_indep_raw.csv"),
        "script_path": str(FROZEN / "code" / "appendix_B_beta_dist.py"),
        "recomputed_vs_archived": {},
    }
    seeds = {"legacy": 9001, "independent": 9001 + 100}  # SEED, SEED+100 in code
    reps = {"legacy": 300, "independent": 150}
    for mode in ("legacy", "independent"):
        d = raw[raw["mode"] == mode]
        bh = d["beta_hat"].values.astype(float)
        bi = d["beta_inj"].values.astype(float)
        ok = np.isfinite(bh)
        bh, bi = bh[ok], bi[ok]
        det = d["detect_correct_sign"].values[ok].astype(bool)
        resid = bh - bi
        rec = {
            "n": int(ok.sum()),
            "beta_inj": {
                "median": float(np.median(bi)),
                "sd": float(np.std(bi, ddof=1)),
                "q025": float(np.quantile(bi, 0.025)),
                "q975": float(np.quantile(bi, 0.975)),
                "frac_near_zero_lt0.05": float(np.mean(np.abs(bi) < 0.05)),
                "frac_positive": float(np.mean(bi > 0)),
            },
            "beta_hat": {
                "median": float(np.median(bh)),
                "sd": float(np.std(bh, ddof=1)),
                "q025": float(np.quantile(bh, 0.025)),
                "q975": float(np.quantile(bh, 0.975)),
            },
            "detect_correct_sign_rate": float(np.mean(det)),
            "resid_hat_minus_inj": {
                "median": float(np.median(resid)),
                "sd": float(np.std(resid, ddof=1)),
            },
            "corr_inj_hat": (float(np.corrcoef(bi, bh)[0, 1])
                             if len(bh) > 2 else float("nan")),
        }
        stored = summ[mode]
        fc = {}
        for k in ("n",):
            fc[k] = {"archived": stored[k], "recomputed": rec[k], "match": close(stored[k], rec[k])}
        for grp in ("beta_inj", "beta_hat", "resid_hat_minus_inj"):
            for k in rec[grp]:
                fc[f"{grp}.{k}"] = {"archived": stored[grp][k],
                                    "recomputed": rec[grp][k],
                                    "match": close(stored[grp][k], rec[grp][k])}
        for k in ("detect_correct_sign_rate", "corr_inj_hat"):
            fc[k] = {"archived": stored[k], "recomputed": rec[k],
                     "match": close(stored[k], rec[k])}
        out["recomputed_vs_archived"][mode] = fc

        # Regenerate the injected draws with the exact recipe + seed and compare
        regen = draw_beta_inj(np.random.default_rng(seeds[mode]), reps[mode])
        archived_bi = d["beta_inj"].values.astype(float)
        out.setdefault("injection_regeneration", {})[mode] = {
            "seed": seeds[mode], "reps": reps[mode],
            "max_abs_diff_vs_archived": float(np.max(np.abs(regen - archived_bi))),
            "exact_match": bool(np.allclose(regen, archived_bi, rtol=0, atol=1e-12)),
            "regen_median": float(np.median(regen)),
            "regen_sd": float(np.std(regen, ddof=1)),
            "regen_q025": float(np.quantile(regen, 0.025)),
            "regen_q975": float(np.quantile(regen, 0.975)),
        }

    out["all_fields_match"] = all(
        c["match"] for m in out["recomputed_vs_archived"].values() for c in m.values())
    out["all_injections_regenerated"] = all(
        v["exact_match"] for v in out["injection_regeneration"].values())
    return out


def highrep_manifest():
    m = json.loads((HIGHREP / "manifest.json").read_text())
    return {
        "path": str(HIGHREP / "manifest.json"),
        "gate_reps_per_arm": m.get("gate_reps_per_arm"),
        "grid_reps_per_cell": m.get("grid_reps_per_cell"),
        "elapsed_hours": m.get("elapsed_hours"),
        "model_spec_unchanged": m.get("model_spec_unchanged"),
        "rng_stream_exact": m.get("rng_stream_exact"),
        "generated_by": m.get("generated_by"),
        "note": ("This manifest governs the high-rep GATE/POWER grid run "
                 "(1000 reps/arm gate, 500/cell grid, ~6.03 h). It is the "
                 "production-scale compute context, not the confirm run itself; "
                 "the confirm run is archived separately under results/appendix/."),
    }


def main():
    result = {
        "item2_confirm": check_confirm(),
        "item3_beta_dist": check_betadist(),
        "highrep_manifest": highrep_manifest(),
        "verification_python": "C:/Users/MIT/AppData/Local/Programs/Python/Python312/python.exe",
    }
    p = OUT / "verify_archived_runs_v2.json"
    with open(p, "w") as fh:
        json.dump(result, fh, indent=2)
    # Compact human-readable summary
    c2 = result["item2_confirm"]
    c3 = result["item3_beta_dist"]
    print("=== ITEM 2: confirm_independent_production ===")
    print("  raw rows:", c2["raw_total_rows"], "modes:", c2["raw_modes"],
          "per-arm n:", c2["raw_per_arm_n"])
    print("  scale matches M1000/FGK400:", c2["scale_matches_claim"],
          "| mode independent:", c2["mode_is_independent"])
    for s in ("M", "FGK"):
        sep = c2["recomputed_vs_archived"][s]["arm_separation_nullSD"]
        print(f"  {s} arm separation: archived={sep['archived']:.9f} "
              f"recomputed={sep['recomputed']:.9f} match={sep['match']}")
    print("  ALL confirm fields recomputed from raw match archived:", c2["all_fields_match"])
    print("\n=== ITEM 3: appendix_B_beta_dist ===")
    for mode in ("legacy", "independent"):
        fc = c3["recomputed_vs_archived"][mode]
        med = fc["beta_inj.median"]["recomputed"]
        sd = fc["beta_inj.sd"]["recomputed"]
        q025 = fc["beta_inj.q025"]["recomputed"]
        q975 = fc["beta_inj.q975"]["recomputed"]
        det = fc["detect_correct_sign_rate"]["recomputed"]
        reg = c3["injection_regeneration"][mode]
        print(f"  {mode}: beta_inj median={med:.6f} sd={sd:.6f} "
              f"[{q025:.4f},{q975:.4f}] detect={det:.6f}")
        print(f"       regenerated injection: max|diff|={reg['max_abs_diff_vs_archived']:.3g} "
              f"exact={reg['exact_match']}")
    print("  ALL beta_dist fields match archived:", c3["all_fields_match"],
          "| all injections regenerated:", c3["all_injections_regenerated"])
    print(f"\nWROTE {p}")


if __name__ == "__main__":
    main()
