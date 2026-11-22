# Orchestrator cross-check: parallax/distance errors from frozen quality table

Computed from `radius-valley-age-evolution/data/gaia_quality_DR3_v2.csv` (2,420 hosts; columns ruwe, parallax, parallax_error):

- Hosts passing RUWE ≤ 1.4 AND parallax SNR > 5: **2,238**.
- RUWE > 1.4: **182** hosts. SNR ≤ 5: **20** hosts.
- Fractional inverse-parallax distance error σ_d/d = σ_π/π: median **0.0088** (all), **0.0082** (passing cuts), 95th percentile **0.0462** (all), **0.0354** (passing).
- For the 20 hosts failing SNR>5: median σ_d/d = **0.350**, 95th percentile 2.42 (unreliable distances).

Conclusion for the disclosure: inverse parallax without zero-point correction is adequate at the median of the adopted sample (<1% distance error) but not for the small low-SNR tail, which the production RUWE/SNR cut removes. Cross-check for the P3 Gaia agent.
