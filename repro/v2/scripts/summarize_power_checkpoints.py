import json, pathlib, csv
SCI = pathlib.Path(r'REPO\manuscripts\q1-paper\v2-2026-11\results\science')
CK = SCI / 'checkpoints'
eff = json.load(open(SCI / 'effect_size_derivation.json'))
pts = eff['bracket_points']
names = ['null'] + [n for n in pts if n != 'null']

def load(tag):
    p = CK / f'curve_{tag}.json'
    if not p.exists():
        return None
    d = json.load(open(p))
    return [r for r in d['rows'] if r.get('beta_pooled') is not None]

null = load('null')
null_vals = [r['beta_pooled'] for r in null]
import numpy as np
nmed, nsd = float(np.median(null_vals)), float(np.std(null_vals))
REAL = +0.0458  # cool-sample real beta/IQR from frozen gate JSON
real_pct_null = float(np.mean(np.array(null_vals) <= REAL))

rows, summary = [], []
for i, name in enumerate(names[1:], start=0):
    sig = load(f'sig_{i}')
    if not sig or len(sig) < 50:
        continue
    sv = [r['beta_pooled'] for r in sig]
    smed = float(np.median(sv))
    sep = (smed - nmed) / nsd
    pct_sig = float(np.mean(np.array(sv) <= REAL))
    verdict = 'PASS' if (real_pct_null <= 0.05 and 0.025 <= pct_sig <= 0.975) else 'FAIL'
    rec = dict(point=name, beta_inj=pts[name], n_null=len(null_vals), n_signal=len(sv),
               null_median=nmed, null_sd=nsd, signal_median=smed,
               arm_separation_nullSD=sep, real_pct_null=real_pct_null,
               real_pct_signal=pct_sig, verdict=verdict)
    summary.append(rec)
    for arm, vals in (('null', null_vals), ('signal', sv)):
        for j, v in enumerate(vals):
            rows.append(dict(point=name, beta_inj=pts[name], arm=arm, rep=j,
                             beta_IQR=v))

out = dict(status='partial Monte Carlo from resumable checkpoints; incomplete bracket',
           structural_note='PASS requires real distinguishable from null (p<=0.05), which depends only on null arm and real coefficient; cool real sits at 55th pct of null so verdict is FAIL at every effect size regardless of injected magnitude.',
           n_null=len(null_vals), null_median=nmed, null_sd=nsd,
           real_beta=REAL, real_pct_null=real_pct_null, points=summary)
json.dump(out, open(SCI / 'power_vs_effect_summary.json', 'w'), indent=2)
with open(SCI / 'power_vs_effect_raw.csv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=['point','beta_inj','arm','rep','beta_IQR'])
    w.writeheader(); w.writerows(rows)
print('points:', [(s['point'], round(s['beta_inj'],4), s['n_signal'], round(s['arm_separation_nullSD'],3), s['verdict']) for s in summary])
