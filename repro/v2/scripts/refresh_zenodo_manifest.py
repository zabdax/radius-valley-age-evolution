"""Refresh the Zenodo file manifest after the final compile + spot check."""
import csv, hashlib, json, pathlib

ROOT = pathlib.Path(r'REPO\manuscripts\q1-paper\v2-2026-11')
AUD = ROOT / 'results' / 'audit'
UNSUITABLE_PREFIXES = ('_',)
UNSUITABLE_NAMES = {'quick_large_effect_check.py', 'find_table_end.py',
                    'make_power_figure.py', 'summarize_power_checkpoints.py',
                    'build_P3_audit_tables.py', 'P4_REWRITE_PLAN.md'}
UNSUITABLE_PARTS = ('__pycache__',)

rows = []
for p in sorted(ROOT.rglob('*')):
    if not p.is_file():
        continue
    rel = p.relative_to(ROOT).as_posix()
    if any(part in UNSUITABLE_PARTS for part in p.parts):
        continue
    h = hashlib.sha256(p.read_bytes()).hexdigest()
    unsuitable = (p.name.startswith(UNSUITABLE_PREFIXES) or
                  p.name in UNSUITABLE_NAMES or
                  rel.startswith('results/audit/') and 'absolute' in p.read_text(errors='ignore')[:500].lower())
    rows.append(dict(path=rel, bytes=p.stat().st_size, sha256=h,
                     deposit='no' if unsuitable else 'yes'))

with open(AUD / 'ZENODO_FILE_MANIFEST.csv', 'w', newline='', encoding='utf-8') as f:
    w = csv.DictWriter(f, fieldnames=['path', 'bytes', 'sha256', 'deposit'])
    w.writeheader(); w.writerows(rows)

manifest = json.load(open(AUD / 'ZENODO_RELEASE_MANIFEST.json', encoding='utf-8'))
manifest['files'] = [{k: r[k] for k in ('path', 'bytes', 'sha256')}
                     for r in rows if r['deposit'] == 'yes']
manifest['notes'] = ('Preparation only; no DOI minted. Excludes scratch diagnostics, '
                     'failed one-off checks, and audit notes containing absolute local paths. '
                     'Refreshed after final compilation and the 60-replicate spot check.')
json.dump(manifest, open(AUD / 'ZENODO_RELEASE_MANIFEST.json', 'w', encoding='utf-8'), indent=2)
print('files', len(rows), 'deposit', sum(r['deposit']=='yes' for r in rows),
      'bytes', sum(int(r['bytes']) for r in rows))
