# Between-study transfer from published correlate estimates

Run 2026-11-22, following rank 3 of `docs/power-and-framing.md`. Code:
`spikes/summary_pool.py`. Input: `spikes/data/summary_estimates.csv` (published
summary numbers only, with a source for every row). Outputs:
`/mnt/project-files/summary-pool/` (`pooled.json`, `loo_transfer.csv`).

## 1. Answer

- **Neutralising titre transfers across four vaccine platforms with no
  detectable penalty.** COVID-19 hazard ratio per 10-fold ID50 in vaccinees:
  mRNA 0.42, Ad26 0.38, protein 0.39, ChAdOx1 0.28. Between-trial SD of the
  log-HR (τ) is estimated at 0.00, with a 95% upper bound of 0.49 per
  10-fold. The pooled slope is −0.94 per 10-fold (HR 0.39), so the bound means
  a new platform's slope could still differ from the pooled one by up to about
  half its size. Leave-one-trial-out attenuations are all within ±0.38 and none
  stands out (largest |z| = 0.68; null p = 0.88).
- **Binding IgG transfers less well, though the data cannot prove it.** For
  spike IgG, τ = 0.26 (I² = 43%, 95% upper bound 1.34). Taking mRNA (COVE) as
  the target, attenuation is +0.49 ± 0.25: the IgG slope in COVE is about half
  the pooled slope of the other three. The heterogeneity test does not reject
  (p = 0.15). This fits the known fact that binding-antibody levels differ by
  roughly 80-fold between platforms (GM 2,652 BAU/ml in COVE vs 34 in
  ENSEMBLE), so the same BAU value means different things across platforms.
- **Influenza HAI, 37 datasets.** Coudeville et al. 2010 fitted exactly this
  random-effects model to 15 studies (37 strain × study datasets, 5,899
  people, 1,304 cases). The between-dataset SD of the protection-curve slope is
  0.38 [0.10, 0.76], on a mean slope of 1.30 [1.00, 1.69] per natural-log
  titre. So a borrowed HAI slope carries a between-study SD of about 29%
  (8–58%) of its own size. This is the best-determined τ we have, because 37
  datasets is ten times our COVID k.
- **What this means for the project.** For mature antibody correlates
  measured on a harmonised assay, the between-study transfer penalty in the
  slope is modest: roughly 0–50% of the slope's size for COVID nAb, and about
  30% for HAI. That is smaller than the within-study PNG village loss
  (attenuation 0.87 per SD, essentially all of the source effect). On a per-SD
  scale, the COVID ID50 bound is about 0.28 per SD, taking the SD of log10
  ID50 as 0.58 (derived from AZD1222's per-SD and per-10-fold HRs). Whether a
  correlate transfers looks like a property of the marker (nAb vs binding IgG)
  more than of the platform.

## 2. Inputs

All four COVID rows come from one table built by the trial statisticians for
cross-trial comparison: Benkeser et al. 2023, npj Vaccines 8:36, Table 3 (US
sites, harmonised assays, WHO international units). That table was read in
full text (PMC10005913). Each row was checked against its trial's paper where
the paper gives the number in text (AZD1222 Table 2A). COVE's own per-10-fold
values are in a figure only (Gilbert 2022, Fig. 3A), so they rest on Table 3.

| Trial | Platform | Marker day | ID50 HR per 10-fold [95% CI] | Spike IgG HR per 10-fold | Cases with marker |
| --- | --- | --- | --- | --- | --- |
| COVE | mRNA-1273 | D57 | 0.42 [0.27, 0.65] | 0.66 [0.50, 0.88] | 36 |
| ENSEMBLE-US | Ad26.COV2.S | D29 | 0.38 [0.13, 1.12] | 0.62 [0.28, 1.37] | 29 |
| PREVENT-19 | NVX-CoV2373 | D35 | 0.39 [0.19, 0.82] | 0.36 [0.20, 0.64] | 12 |
| AZD1222 | ChAdOx1 | D57 | 0.28 [0.10, 0.77] | 0.32 [0.14, 0.76] | 22 (ID50), 33 (IgG) |

SE = (log hi − log lo) / 3.92. ENSEMBLE all regions (ID50 0.49 [0.29, 0.81],
92 cases, adjusted for region) is kept as a sensitivity row, not pooled, since
it overlaps ENSEMBLE-US.

## 3. Methods, mapped to PR 10

- **Primary metric.** Leave-one-trial-out transfer: source = random-effects
  pool of the other trials, target = the held-out trial. Attenuation =
  sign(β_S)·(β_S − β_T), positive when the correlate is weaker in the target.
  This is the same metric as PR 10, on log-HR per 10-fold rather than per SD,
  because the harmonised IU scale is a common unit across trials and per-SD
  scaling would mix in each trial's titre spread.
- **Null check.** The summary-level stand-in for the random-split check. The
  only thing that can be shuffled at summary level is sampling error, so the
  null is a world where every trial's estimate is drawn from N(μ_fixed, SE_i²)
  (2,000 draws). τ, Q and the largest leave-one-out |z| are recomputed in each
  draw. The observed values are compared with that null.
- **Pooling.** REML τ with a Q-profile CI, Knapp–Hartung CI for the mean, a
  prediction interval, and a grid posterior with half-normal priors on τ
  (SD 0.25, 0.5, 1).

| Family | k | Pooled log-HR per 10-fold [KH 95%] | τ (REML) | τ 95% (Q-profile) | I² | Null p (τ / max LOO z) |
| --- | --- | --- | --- | --- | --- | --- |
| ID50 | 4 | −0.94 [−1.48, −0.39] | 0.00 | [0, 0.49] | 0% | 1.00 / 0.88 |
| Spike IgG | 4 | −0.68 [−1.31, −0.06] | 0.26 | [0, 1.34] | 43% | 0.15 / 0.18 |

Bayesian τ (posterior median [95%]) for ID50 is 0.13 [0.01, 0.45] under a
half-normal(0.25) prior, 0.18 [0.01, 0.72] under (0.5), and 0.22 [0, 1.06]
under (1). With four trials, the upper end of τ is set by the prior, not the
data.

## 4. Caveats

- **Four trials is few.** The ID50 bound of 0.49 is better than typical: by
  simulation, four trials with SE about 0.3 give a median upper bound of 0.94
  when the true τ is 0. These four agree unusually closely (Q = 0.52 on 3 df).
  A typical upper bound below 0.3 needs about 20 trials of this precision.
  Only summary data can supply that many.
- **Platform is confounded with marker timing and follow-up.** Marker days
  are D29, D35 and D57, and follow-up runs 73 to 126 days. The confounding
  makes heterogeneity more likely, not less, so the small ID50 τ is if
  anything conservative.
- **These are correlates of risk in vaccinees, not VE curves.** Exposure
  differences between trials can enter the hazard ratio. All four are US
  sites in the ancestral or early-variant period, which limits that.
- **ENSEMBLE regional HRs are only in a figure** (Extended Data Fig. 7), so
  the cross-region (population) gap could not be added without digitising.
  That is the next cheap addition, and it would be E2 at summary level.
- **Coudeville gives τ from its own Bayesian fit.** The per-dataset counts are
  not in the paper (the original studies date from 1945 to 1997). Their slope
  parameter is on the protection-curve scale, not a hazard ratio, so its
  relative τ (29%) compares with ours and its absolute τ does not.
- **Publication and selection.** Only trials whose correlates papers came out
  are here, and all four come from one harmonised programme. That makes them
  comparable, and also not a random sample of platforms.

## 5. Next steps, cheapest first

1. Digitise ENSEMBLE regional HRs (Extended Data Fig. 7) and the COVE Fig. 3A
   points to cross-check Table 3.
2. Add other harmonised multi-trial sets with HRs in text: COVID booster and
   variant correlates (COVAIL, the Omicron-era CoVPN papers), RSV maternal and
   older-adult programmes, pneumococcal OPA/IgG trials.
3. Pull per-dataset counts for the HAI studies from the 15 source papers, to
   refit Coudeville's model with study-type covariates (challenge vs field,
   vaccinated vs not). That gives τ per gap type, which is the project's real
   question.
