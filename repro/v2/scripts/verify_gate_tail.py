"""Independently re-derive the gate decision from archived raw replicates.

The injected atmospheric-loss effect is negative, so the one-sided
distinguishability statistic is the LOWER-tail placement of the real coefficient
in the null ensemble, P(null <= real). This script recomputes that placement
from the committed replicate CSVs, evaluates both the correct lower-tail and
the previously used upper-tail convention, and verifies that the PASS/FAIL
verdict is unchanged either way.

Run from anywhere: python repro/v2/scripts/verify_gate_tail.py
"""
import json
from pathlib import Path
import numpy as np
import pandas as pd

from _paths import REPO, SCI

HIGH = REPO / "results" / "highrep"
STRAT = REPO / "results" / "stratified_gate" / "M"
REAL_M = 0.04576344675022682

def placement(null, real):
    """Empirical lower-tail placement P(null <= real)."""
    return float(np.mean(np.asarray(null, float) <= real))

def check(name, null, signal, real, q_sig):
    p_low = placement(null, real)
    p_up = 1.0 - p_low
    inside = q_sig[0] <= real <= q_sig[1]
    verdict_low = "PASS" if (p_low <= 0.05 and inside) else "FAIL"
    verdict_up = "PASS" if (p_up <= 0.05 and inside) else "FAIL"
    row = {
        "sample": name,
        "real_beta": round(float(real), 4),
        "P_null_le_real": round(p_low, 4),
        "P_null_ge_real": round(p_up, 4),
        "real_inside_signal_95pct": bool(inside),
        "verdict_lower_tail_correct": verdict_low,
        "verdict_upper_tail_legacy": verdict_up,
    }
    print(json.dumps(row))
    return row

def main():
    g = pd.read_csv(HIGH / "task2_gate_rep1000_replicates_raw.csv", keep_default_na=False)
    b0 = g[g.arm == "null"].beta_hat_per_iqr.values
    b1 = g[g.arm == "signal"].beta_hat_per_iqr.values
    t2 = json.loads((HIGH / "task2_calibration_kepler_only_trim0.95_rep1000.json").read_text())
    rows = [check("FGK control", b0, b1, float(t2["beta_real"]),
                  (float(np.quantile(b1, .025)), float(np.quantile(b1, .975))))]

    n = pd.read_csv(STRAT / "strat_gate_M_vmag_null_r1000_raw.csv").beta_pooled.dropna().values
    s = pd.read_csv(STRAT / "strat_gate_M_vmag_signal_r1000_raw.csv").beta_pooled.dropna().values
    rows.append(check("Cool hosts", n, s, REAL_M,
                      (float(np.quantile(s, .025)), float(np.quantile(s, .975)))))

    # Cross-check the published placements in the archived gate summaries.
    assert abs(rows[0]["P_null_le_real"] - 0.456) < 0.002, "FGK placement mismatch"
    assert abs(rows[1]["P_null_le_real"] - 0.553) < 0.002, "cool placement mismatch"
    assert all(r["verdict_lower_tail_correct"] == "FAIL" for r in rows)
    assert all(r["verdict_upper_tail_legacy"] == "FAIL" for r in rows)

    out = {
        "status": "VERIFIED",
        "note": "Correct lower-tail test and legacy upper-tail test both return "
                "FAIL for both samples; published verdicts are unchanged.",
        "rows": rows,
    }
    (SCI / "gate_tail_verification.json").write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", SCI / "gate_tail_verification.json")

if __name__ == "__main__":
    main()
