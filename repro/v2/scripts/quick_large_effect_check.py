"""60-replicate large-effect spot check for the calibrated cool-host gate.

FIX (2026-11-21): the original version failed on Windows with
"freeze_support() ... Safe importing of main module".  Cause: all work was at
module level, so the 'spawn' multiprocessing start method re-executed the
script in every child process, recursively spawning workers.  Fix: keep the
dynamic gate import at module level (so 'gate' is importable when unpickling)
and put all executable work behind an `if __name__ == "__main__"` guard.
"""
import importlib.util, json, sys, pathlib
import numpy as np

ROOT = REPO
F = ROOT / 'radius-valley-age-evolution'
spec = importlib.util.spec_from_file_location('gate', F / 'code' / '15_stratified_gate.py')
g = importlib.util.module_from_spec(spec)
sys.modules['gate'] = g          # required for pickling under Windows spawn
spec.loader.exec_module(g)

OUT = ROOT / 'manuscripts' / 'q1-paper' / 'v2-2026-11' / 'results' / 'science'


def main():
    g.CFG = g.SAMPLES['M']
    g.TARGET_FRAC = g.CFG['target_frac']
    g.VALLEY = g.CFG['valley']
    g.GAMMA_LOGP = g.CFG['gamma_logp']
    g.OUT = OUT
    g.CKPT = OUT / 'checkpoints'
    g.CKPT.mkdir(parents=True, exist_ok=True)

    base = g.load_control(g.CFG)
    comp = g.m9.build_design(base, 'vtan', valley=g.VALLEY, fe=False,
                             weights=base.wgt.values)['init_comp']
    cut = float(base['sy_vmag'].median())
    rows = []
    for b, tag in [(0.0, 'quick_null'), (-0.8824, 'quick_gaidos'), (-0.1419, 'quick_primary')]:
        df = g.run_arm(base, comp, 'sy_vmag', cut, b, 60, 777000, tag, workers=7)
        v = df['beta_pooled'].dropna().values
        med = float(np.median(v))
        print(f'{tag}: beta_inj={b}, n={len(v)}, median={med:+.4f}, '
              f'mean={np.mean(v):+.4f}, sd={np.std(v):.4f}', flush=True)
        rows.append({'point': tag, 'beta_inj': b, 'n': len(v),
                     'median': med, 'mean': float(np.mean(v)),
                     'sd': float(np.std(v))})
    nul = next((r for r in rows if r['point'] == 'quick_null'), None)
    if nul is not None:
        nmed = float(nul['median'])
        nsd = float(nul['sd'])
        for r in rows:
            r['arm_separation_nullSD'] = (r['median'] - nmed) / nsd
        # reuse the null SD from the row (populated above)
    json.dump({'status': 'completed', 'reps_per_arm': 60,
               'note': '60-replicate spot check; use only as sanity check, not the production curve',
               'rows': rows},
              open(OUT / 'quick_large_effect_check.json', 'w'), indent=2)
    print('wrote', OUT / 'quick_large_effect_check.json')


if __name__ == '__main__':
    main()
