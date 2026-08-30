# Adversarial audit scripts

The manuscript (Methods, "Adversarial Audits and Analytical Corrections")
reports that the pipeline underwent independent adversarial audits whose
corrections were applied before the final results were reported. The scripts
in this folder are those audits, preserved as run, organized in three rounds:

| Folder | Focus | Scripts |
|---|---|---|
| `round1_data_pipeline/` | Gaia crossmatch integrity, Kepler sample selection, U/V/W kinematics, end-to-end consistency | `audit_crossmatch.py`, `audit_kepler.py`, `audit_uvw.py`, `audit_tests.py`, `audit_final.py` |
| `round2_age_proxies/` | Radius-valley definition, gyrochronology pipeline, `gyrointerp` validity domain, TAP query counts, McQuillan (2014) rotation crossmatch | `audit1_valley.py`, `audit2_gyro.py`, `audit3_gyrointerp.py`, `audit4_tap_counts.py`, `audit5_mcquillan.py` |
| `round3_model_checks/` | Insolation guard, refit variants, gate fidelity/power, Hessian self-test | `a1_guard_insolation.py`, `a2_refit_variants.py`, `a3_fidelity_power.py`, `a4_selftest_hessian.py` |

These are adversarial probes, not part of the production pipeline: each was
written to break a specific assumption, and several succeeded (the mixture
identifiability collapse, the Gauss-Hermite normalization error, the radius
error outlier, component degeneracy, the objective scaling bug, the anchor
matching bug, and the premature "precise null" label described in the
manuscript were all found this way). They are archived unmodified; paths and
minor details reflect the state of the project at the time each audit ran.