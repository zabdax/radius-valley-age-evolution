"""Rank correlation between the kinematic age proxy and independent TIME Table
gyrochronological ages across the finite-age cool-host crossmatch.

Frozen inputs only. The 65 finite-age matches are exactly those plotted in
Figure 2. The pipeline gates are NOT rerun on these ages: they lie outside the
Angus+19 validity domain and are used only for proxy validation.
"""
import json, os, re
import numpy as np
import pandas as pd
from scipy import stats
from _paths import BUNDLE, REPO, SCI, AUDIT, FIGS, PROV, DATA  # portable paths

SCI = str(SCI)
cm = json.load(open(os.path.join(SCI, 'time_table_crossmatch.json'), encoding='utf-8'))
kin = pd.read_csv(os.path.join(DATA, 'hosts_kinematics_M.csv'))


def norm(s):
    s = re.sub(r'\s+', ' ', str(s).strip().lower())
    s = re.sub(r'^(gliese|gl|gj)\s+', 'gj ', s)
    return re.sub(r'\s+[a-z]$', '', s)


kmap = {norm(r.hostname): float(r.vtan) for _, r in kin.iterrows() if pd.notna(r.vtan)}
rows = []
for m in cm['matched_hosts']:
    for host, st, age in zip(m['sample_hosts'], m['age_statuses'], m['age_Gyr']):
        if st != 'age_estimate' or age is None:
            continue
        key = norm(host)
        if key not in kmap and m.get('time_table_names'):
            key = norm(m['time_table_names'][0])
        rows.append(dict(host=host, age_Gyr=float(age), vtan=kmap.get(key, np.nan),
                         additional=not m['already_one_of_12_gyro']))
df = pd.DataFrame(rows).dropna(subset=['vtan'])
assert len(df) == 65, len(df)

rho, p = stats.spearmanr(df.vtan, df.age_Gyr)
add = df[df.additional]
rho_add, p_add = stats.spearmanr(add.vtan, add.age_Gyr)

out = {
    "status": "COMPLETE",
    "description": "Spearman rank correlation of Gaia tangential velocity with "
                   "independent TIME Table gyrochronological age, cool-host crossmatch",
    "n_matched": len(df),
    "n_additional": int(len(add)),
    "n_overlapping_pipeline_gyro": int((~df.additional).sum()),
    "rho_all": float(rho), "p_all": float(p),
    "rho_additional_only": float(rho_add), "p_additional_only": float(p_add),
    "source_crossmatch": "results/science/time_table_crossmatch.json",
    "source_kinematics": "radius-valley-age-evolution/data/hosts_kinematics_M.csv",
    "note": "Validation only. These ages are not used in the primary gate because "
            "they fall outside the Angus+19 validity domain (3800-6200 K, P_rot <= 45 d).",
}
with open(os.path.join(SCI, 'timetable_spearman.json'), 'w', encoding='utf-8') as f:
    json.dump(out, f, indent=2)
print(json.dumps(out, indent=2))
