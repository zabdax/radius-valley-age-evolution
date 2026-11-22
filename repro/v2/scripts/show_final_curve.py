import json
p = r'REPO\manuscripts\q1-paper\v2-2026-11\results\science\power_vs_effect_summary_FINAL.json'
d = json.load(open(p))
print('status', d.get('status'), 'completed', d.get('points_completed'))
print('runtime_s', d.get('runtime_s'), 'reps', d.get('reps_per_arm'))
for x in d['points']:
    print(f"{x['point']:26s} b={x['beta_inj']:+.4f} sep={x['arm_separation_nullSD']:+.3f} "
          f"pct_null={x['real_pct_in_null']:.3f} det={x['detect_rate']:.3f} "
          f"verdict={x['verdict']} n={x['reps_signal']}")
