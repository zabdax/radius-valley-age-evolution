# TIME Table (Gaidos et al. 2023) — calibration limits, verbatim

**Primary note:** `research/notes/mnras-000-116-2022.md` (note id `mnras-000-116-2022`,
raw PDF `research/raw/mnras-000-116-2022.pdf`, DOI `arXiv:2301.12109`).
**Full machine-readable table:** Zenodo DOI `10.5281/zenodo.7578269`, file `time_table.mrt`
(CDS format, 249 rows; cached copy at `time_table_mrt_source.txt`, md5 `581d26ae64ae79eb4c9b309bc7d6466a`).

All quotes below are transcribed verbatim from the vault note; line numbers are the
1-based line numbers of that file. Ligatures and dash glyphs are preserved as extracted.
The vault note's printed Table 1 contains only the first 50 of 249 entries — see
note line 988: *"Only the first 50 entries are shown; the complete machine-readable
table is provided on Zenodo (DOI: 10.5281/zenodo.7578269)."*

---

## 1. Sample size and the 227 / 8 / 14 split

**Note lines 57–62 (ABSTRACT):**
> "We curated a catalog of 249 late K-
> and M-type (𝑇eﬀ=3200-4200K) exoplanet host stars with established rotation periods, and
> applied empirical, temperature-dependent rotation-age relations based on relevant published
> gyrochrones, including one derived from observations of the 4 Gyr-old open cluster M67. We
> estimated ages for 227 of these stars, and upper limits for 8 others, excluding 14 which are
> too rapidly rotating or are otherwise outside the valid parameter range of our gyrochronology."

**Note lines 1794–1797 (SUMMARY AND DISCUSSION):**
> "We cataloged 249
> stars with rotation periods, 227 of which we are able to estimate ages
> with a median error of 20% and mode of 14%, and to an additional 8
> we assign upper limits (Table 6)."

**Note lines 21–26 (vault note frontmatter `summary:`):**
> "Gaidos et al. (2023, arXiv:2301.12109) curate the TIME Table: 249 late-K/M
> exoplanet hosts (Teff 3200-4200 K) with rotation periods, applying empirical Teff-dependent
> gyrochrones to estimate ages for 227 stars (median error 20%, mode 14%), upper limits
> for 8, and excluding 14 rapid rotators."

**Verified against the machine-readable table:** 249 rows; 227 with finite age estimates,
8 with upper-limit flags (`ul=Y`), 14 with age-assignment failure flags (`fail=Y`). Exact match.

**Apparent internal inconsistency in the paper (flag, do not silently resolve):**
note lines 1275–1277 read
> "Fourteen stars can only be assigned upper limits
> for ages, and of these 8 are younger than 1 Gyr at the 95-percentile"

This swaps the abstract's 8 upper limits / 14 exclusions. The machine-readable table
supports the abstract's numbers (8 `ul=Y`, 14 `fail=Y`), so treat the line-1275 wording
as a prose inconsistency, not as a different population.

---

## 2. ~20% median age error, ~14% mode

**Note lines 1795–1798:**
> "227 of which we are able to estimate ages
> with a median error of 20% and mode of 14%, and to an additional 8
> we assign upper limits (Table 6). Our fractional error is signiﬁcantly
> higher than the 5-10% estimated by Otani et al. (2022), probably
> because we include additional potential sources of error."

Also stated in the note frontmatter, line 23 (quoted in §1 above).

---

## 3. The <3 Gyr artefact peak

**Note lines 1803–1813:**
> "The age distributions of both Kepler and non-Kepler host stars
> peak at around 3 Gyr with a steady decline to near zero at 10
> Gyr, the age of the Galactic disk (Fig. 7). The resemblance of the
> actual and “mock" populations (the latter with a uniform 0-10 Gyr
> age distribution) shown in Figure 8 indicates that the structure of
> the distribution, particularly at young ages, is partly due to the
> discontinuous and limited coverage of the current gyrochronology
> and dispersion of the distribution due to error. The peaks at <3
> Gyr correspond approximately with the location of the calibration
> ages, and are likely artefacts due to ages that cannot be assigned
> in certain regions of 𝑇eﬀ-𝑃rot space and have only upper limits
> assigned."

**Note lines 1814–1818 (competing selection biases):**
> "Other eﬀects impacting the derived distribution include
> the opposing biases against detection of planets around younger,
> more rapidly rotating, and more active stars (Miyakawa et al. 2022),
> and against detection of rotational variability among older, slowly
> rotating, less active stars (Morris 2020)."

---

## 4. 75 K Teff error → 0–17% age error

**Note lines 614–619 (§5.3 Error in eﬀective temperature):**
> "Eﬀective temperature 𝑇eﬀis usually used as the independent vari-
> able of a gyrochronology since it is set by radiated energy per unit
> area and is thus related to the eddy velocity and magnetodynamo
> strength in the convective envelope of a star. Precise and accurate
> determination of 𝑇eﬀis a classic problem in stellar astrophysics,
> and a particularly acute issue for the coolest stars. Recent calibra-
> tions of the temperature scale of M dwarfs based on interferometric
> measurement of stellar radii, (Gaia) parallaxes, and the Stefan-
> Boltzmann law (Boyajian et al. 2012; Mann et al. 2015) permit
> precision as good as 50 K, but most host stars have 𝑇eﬀwith pre-
> cision no better than 75 K."

**Note lines 619–622:**
> "The M67 gyrochrone of Dungee et al.
> (2022) was ﬁt to the 𝑔-𝑖colors of stars, and converted to 𝑇eﬀusing
> the temperature scale of Mann et al. (2015). This combined with
> the slope of the M67 gyrochrone will produce formal uncertainties"

**Note lines 647–651 (Figure 3 caption and continuation):**
> "Figure 3. Intrinsic dispersion in the 4 Gyr-old M67 gyrochrone (dashed
> black lines) established by Dungee et al. (2022) and dispersion produced by
> errors in 𝑇eﬀof 75K (red lines).
> of 0-17% in age, with the largest value at the cool (3200K) end,
> decreasing to zero near 3900K, and rising to 12% at 4200K."

---

## 5. Untested `n=0.62` spin-down exponent / cool-star extrapolation

**Note lines 470–478 (§4 calibration):**
> "If the only available calibration point was that of M67 (4 Gyr)
> then we calculate the age of the star using the Skumanich-like
> spin-down law that Dungee et al. (2022) found by comparing the
> rotational sequence of M67 with that of 2.7 Gyr-old Ruprecht 147
> (Curtis et al. 2020)
> 𝑡= 𝑡0 (𝑃/𝑃0)1/𝑛,
> (2)
> with 𝑛= 0.62. However, if the age derived in this manner was <2.7
> Gyr (i.e., that of Ruprecht 147) we took this to be an upper limit."

**Note lines 1839–1848 (§7, explicit statement of the untested assumption):**
> "Our gyrochronology assumes that the narrow rotational se-
> quence observed among the late K dwarfs and the warmer M dwarfs
> in the 2.7 Gyr-old Ruprecht 147 cluster (Curtis et al. 2020) extends
> to 3200K by 2.7 Gyr, and that the 𝑛= 0.62 power-law spin-down
> derived by (Dungee et al. 2022) also applies to cooler M dwarfs
> at later times. This assumption could fail if the rotational sequence
> among M67 M dwarfs was formed by stalling, rather than a tran-
> sition from saturated to un-saturated braking laws. Core-envelope
> re-coupling is expected to become weaker and take longer towards
> the fully convective boundary, which could mean that a rotational"

**Note lines 1858–1861 (how it could be resolved):**
> "sequence appears much later, or not at all. This depends in part
> on the 𝑇eﬀor mass dependence of the core-envelope coupling time
> vs. 𝜏𝑐. This issue can only be resolved by deeper monitoring of
> Ruprecht 147 or a cluster of similar age."

**Note lines 1779–1783 (cool edge is an unconstrained extrapolation):**
> "which reaches 120 days at 3250K (Fig. 3) and at cooler temperatures
> is essentially an unconstrained extrapolation. Thus the gyrochronol-
> ogy suggests an age older than M67, as expected, but extension of
> M dwarf gyrochrones into the fully convective area is needed before
> assigning any robust age."

---

## 6. Binary caveat

**Note lines 1206–1208:**
> "Binaries: Twenty-six of the 249 planet host stars are known
> to have stellar companions. This 11% fraction is much lower than
> 26.8±1.4% among ﬁeld M dwarfs (Winters et al. 2019)."

**Note lines 1240–1250 (why binaries distort rotation):**
> "The rotational history of stars in multiple systems can diﬀer
> greatly from that of single stars. Stars with stellar companions are
> more likely to be rapidly rotating relative to single stars of the same
> age/mass (Kraus et al. 2012; Simonian et al. 2019). Very close (≪1
> au) binaries transfer orbital AM to rotational AM via tidal torques
> (Fleming et al. 2019). At moderate separations (∼100 au) a stellar
> companion will truncate a disk and shorten its viscous lifetime
> (Cieza et al. 2009; Kraus et al. 2012; Rosotti & Clarke 2018).
> This removes a sink of stellar angular momentum, allowing the
> star to spin-up unimpeeded during pre-main sequence contraction
> (Messina 2019)."

**Note lines 1251–1256 (the explicit warning):**
> "We computed ages based on rotation without regard to mul-
> tiplicity, but warn that in the case of binaries, such ages could be
> seriously in error. In this sample, at least, the distribution of rotation
> periods does not appear remarkably diﬀerent from single stars (Fig.
> 1) nor does the age distribution of known binaries appear remark-
> able (Fig. 9), although it is greatly limited by the small sample size."

**Note lines 1260–1267 (problematic old ages, 3 of 4 in binaries):**
> "Anomalously old stars: Four host stars (GJ 27.1, GJ 667 C,
> HD 238090, and HIP 70849) are assigned problematic ages that are
> > 2𝜎older than 10 Gyr, the nominal age of the Galactic disk. [...] Three of the stars are in binaries,
> in which stellar companions could directly or indirectly aﬀect the
> rotation evolution."

---

## Implication for our paper ( restrained reading )

- The 20% median / 14% mode error and the 0–17% Teff-driven error are the
  **authors' own** propagated uncertainties; they are not an independent validation of
  the ages against an absolute clock. Cluster calibration errors themselves are
  explicitly declared out of scope at note lines 602–606:
  *"The ages of the clusters used to calibrate the gyrochronology themselves have signiﬁcant standard errors and systematics; full consideration of these is beyond the scope of this work"*.
- The <3 Gyr peaks are flagged by the authors as artefacts, so they must **not** be
  cited as evidence for a genuine young-age excess.
- Any crossmatched host inherits the binary and cool-star-extrapolation caveats above;
  a TIME Table age is a literature gyro age under a different calibration
  (kiauhoku / M67) from our Angus+19 `gyrointerp` ages, not a drop-in replacement.
