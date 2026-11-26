# Power and framing: are we measuring the transfer penalty the right way?

Spike, 2026-11-22. Code: `spikes/power_rethink.py` (runs in about 5 minutes on
4 cores, no EC2). Outputs: `/mnt/project-files/power-rethink/`
(`sim_precision.csv`, `real_splits.csv`, `meta.json`). Data used: HK H3N2
(Zenodo 5004830), PNG P. vivax (Zenodo 4995749, PR 7 loader), HVTN 505 (CRAN
`vaccine` 1.4.0). AMP was not re-run because its data sits behind the
Dataverse guestbook.

## 1. Short answer

- **Switching metric will not rescue power on its own.** For a single marker,
  the AUC penalty and a likelihood test of the marker coefficient have the same
  power (simulation, section 3). The thing that sets power is **how strong the
  correlate is at home**, not the metric and not the raw event count. A
  50% weakening of a correlate with home AUC 0.6 needs about 600 events per
  side to detect; at AUC 0.7 about 150; at AUC 0.8 about 55. Our HIV, AMP and
  HK sets sit at AUC 0.52 to 0.57, where no feasible sample size works.
- **We should still change the headline estimand**, for reasons other than
  power: the coefficient-scale attenuation stays defined when the source
  signal is weak (calibration slope does not, see AMP's slopes of ±400), it
  separates "the marker means less over there" from "the attack rate is
  different", and it is the one quantity that published papers report, so it
  can be pooled with summary-level data without any DUA.
- **We already have one real, permutation-confirmed penalty.** PNG P. vivax,
  Ilaita to the other villages: the per-SD log-odds coefficient goes from
  −0.84 to about 0, attenuation 0.87 [0.28, 1.46], permutation p = 0.004 with
  the marker selection repeated inside every permutation. This is the first
  result that survives an honest null.
- **The ~100 events per side bar should be replaced** by a rule weighted by
  signal strength (section 4). Dataset hunting should target strong
  correlates, not big event counts on weak ones. Even COVID neutralising titre
  is only moderate per SD, so pooling many studies, most cheaply from published
  summary data, is the realistic route to power.

Update, same day: rank 3 has been run. See `docs/summary-level-pooling.md`.
Across four COVID platforms, the neutralising-titre slope shows no detectable
between-trial heterogeneity (τ 95% upper bound 0.49 log-HR per 10-fold). The
binding-IgG slope may vary more (τ 0.26, not significant).

## 2. Ranked recommendations

| Rank | Do | Cost | What it buys |
| --- | --- | --- | --- |
| 1 | Make **coefficient attenuation** (sign(β_S)·(β_S − β_T), log-OR or log-HR per pooled SD, from one interaction fit on S ∪ T) the primary estimand. Report an upper 95% bound when the test is null. Keep AUC, calibration and Brier as secondary, and suppress calibration slope when β_S's CI covers 0. | 1–2 days: a metric in `penalty/metrics.py` and a report row | A defined, poolable estimand; null results become "attenuation is at most X" instead of "CI covers zero" |
| 2 | Use a **random-split permutation null** for every penalty, with marker selection and CV inside each permutation. | 1 day; code exists in the spike | Calibrated small-n tests; removes the optimism bias of the AUC penalty (+0.03 to +0.08 under the null at 25–130 events here) |
| 3 | **Summary-level meta-analysis** of published per-stratum or per-trial correlate coefficients: USG COVID trials (COVE, ENSEMBLE by region, PREVENT-19, AZD1222 US), influenza HAI cohorts (Coudeville 2010: 15 studies, 1,304 cases), and others. Fit a hierarchical model on log-HR per unit and estimate the between-study SD τ by gap type. | 1–2 weeks, mostly extraction | E1 across platforms with no DUA; many studies, so τ becomes estimable; per-trial power is still modest |
| 4 | **Pick datasets by signal strength** using the rule in section 4. Add RTS,S MALARIA-055 (listed on Vivli) to the Vivli request Imo is already making; keep ENSEMBLE via YODA for its three regions. | No new work, one extra line on the Vivli form | Two large IPD sets with several split axes; check their published β first |
| 5 | **Hierarchical pooling of all splits** (Bayesian random effects on attenuation, with split type as a covariate). Treat τ as the headline: "a geographic gap moves the correlate by about τ". | 2–3 days once 1 and 2 are in | Uses every small split; honest about what is not identifiable (with 3 datasets, τ's CI is 0 to 3.3) |
| 6 | **Challenge studies** as a later E1 source: exposure is fixed, so the attack-rate nuisance disappears and every participant is informative. Influenza HAI/NAI, typhoid Vi, cholera vibriocidal, malaria CHMI. | Uncertain: IPD mostly not public | Clean mechanism signal; small per study, so only useful pooled |
| 7 | Pre-register one primary estimand and one direction per dataset before looking. | An hour per dataset | Stops the 4 metrics × 3 models × 2 directions garden of forking paths |

## 3. What the quick analysis showed

### 3a. Simulation: one standardised marker, 10% attack rate, β_T = β_S/2

| Home AUC | β per SD | Events/side | SD of AUC penalty | SD of attenuation | Power, AUC penalty | Power, attenuation |
| --- | --- | --- | --- | --- | --- | --- |
| 0.6 | 0.37 | 100 | 0.043 | 0.16 | 0.23 | 0.25 |
| 0.6 | 0.37 | 400 | 0.021 | 0.08 | 0.67 | 0.67 |
| 0.7 | 0.77 | 50 | 0.054 | 0.23 | 0.43 | 0.40 |
| 0.7 | 0.77 | 200 | 0.027 | 0.12 | 0.92 | 0.93 |
| 0.8 | 1.29 | 50 | 0.045 | 0.29 | 0.76 | 0.76 |
| 0.8 | 1.29 | 100 | 0.031 | 0.19 | 0.96 | 0.96 |

Full grid in `sim_precision.csv`. Three readings:

1. AUC penalty and attenuation have the same power for a single marker. The
   gain from rank 1 is interpretability and poolability, not power.
2. Calibration slope's null 95% range at AUC 0.6 is [0.35, 2.9] at 100 events
   and [−7, 9] at 25. It is not a usable headline below home AUC 0.7.
3. This simulation scores the source AUC in-sample. The pipeline uses a 5-fold
   CV baseline and ridge, which adds noise, so the pipeline's power is a bit
   lower than this table.

### 3b. Real splits: attenuation with a permutation null (2,000 random splits)

Positive attenuation means the marker's effect is weaker in the target.
Markers are chosen in the source alone where a dataset has several, and the
choice is repeated in each permutation.

| Split | Events S/T | Marker | β_S per SD | Attenuation [95% CI] | p (perm) | AUC penalty | Null mean of AUC penalty |
| --- | --- | --- | --- | --- | --- | --- | --- |
| PNG Pv, Ilaita → other | 75/56 | PvDBPII_P (of 6) | −0.84 | **+0.87 [+0.28, +1.46]** | **0.004** | +0.24 | +0.026 |
| PNG Pv, age < 2.5 → ≥ 2.5 y | 90/41 | PvEBPII (of 6) | −0.42 | −0.18 [−0.80, +0.45] | 0.56 | −0.04 | +0.025 |
| HK H3N2, Perth → Victoria | 135/57 | log2 HAI | −0.53 | +0.23 [−0.23, +0.69] | 0.32 | +0.03 | 0.000 |
| HK H3N2, age < 40 → ≥ 40 | 36/156 | log2 HAI | −0.25 | −0.31 [−0.74, +0.12] | 0.27 | −0.01 | 0.000 |
| HVTN 505, age < 30 → ≥ 30 | 14/11 | IgG V2 (of 3) | −0.30 (SE 0.36) | −0.16 [−1.06, +0.73] | 0.75 | −0.03 | +0.076 |
| HVTN 505, BMI < 25 → ≥ 25 | 12/13 | IgG V2 (of 3) | −0.18 (SE 0.38) | −0.39 [−1.36, +0.58] | 0.46 | −0.07 | +0.072 |

- **PNG village** is a genuine loss of the correlate, not winner's curse: the
  permutation repeats the pick-the-best-of-six step and still gives p = 0.004.
  It is natural immunity, not a vaccine, and the villages differ in
  transmission, so read it as "the same antibody predicts less where exposure
  differs", which is exactly the kind of gap the project is about.
- **HK Perth → Victoria** is not evidence of no penalty. The upper bound,
  0.69 per SD, is larger than the source effect itself (0.53). The data cannot
  rule out total loss.
- **HVTN 505**: the source coefficient is indistinguishable from zero, so any
  attenuation is uninterpretable. The AUC penalty's null mean of +0.07 to +0.08
  means a random split already "shows" a penalty of that size, which is the
  optimism of choosing a marker and its sign in 12–14 cases.
- **Pooling** the three independent primary splits (HK strain, PNG village,
  505 age; DerSimonian–Laird): mean attenuation 0.36 (SE 0.27), τ = 0.35 with a
  Q-profile CI of 0 to 3.3. With three datasets, τ is not identified. This is
  the argument for rank 3: published summary data can supply the dozens of
  studies a τ estimate needs.

## 4. Replacing the ~100 events per side rule

Power to detect a relative attenuation r of a correlate with home coefficient
β (per SD), at two-sided 5% and 80% power, with equal events E per side and a
rare-ish outcome:

    E per side ≈ 22 / (r·β)²        (from the simulation; constant 21–24)

| Home AUC (approx.) | β per SD | E for r = 0.5 | E for r = 0.3 |
| --- | --- | --- | --- |
| 0.60 | 0.37 | ~650 | ~1,800 |
| 0.70 | 0.77 | ~150 | ~410 |
| 0.80 | 1.29 | ~55 | ~150 |

So 100 events per side is plenty for a strong correlate and hopeless for a weak
one. Before requesting a dataset, look up the published per-SD (or per-log10)
hazard ratio and plug it in.

A sobering check: the USG COVID trials report ID50 hazard ratios of about 0.4
per 10-fold (COVE 0.42, ENSEMBLE-US 0.38, PREVENT-19 0.39; from a search
snippet, to confirm against the papers). With a titre SD of roughly 0.5–0.8
log10, that is β ≈ 0.45–0.75 per SD, so detecting a 50% attenuation needs about
150–450 events per side. No single COVID trial has that many vaccinee cases
with markers, so even the strongest correlates we know of only become
measurable when pooled across trials (rank 3). CYD14/CYD15 PRNT and RTS,S
anti-CS should be checked the same way before the Vivli proposal is finalised.

When a source has β near zero (AMP, HVTN 505), there is no correlate to lose,
and the run should say "not estimable" rather than report a penalty.

## 5. Other things in the current framing to fix

- **Correlate of risk vs correlate of protection.** In vaccinee-only data, a
  marker can "transfer" or "fail" because of exposure, not biology. Where a
  placebo arm exists (AMP, HVTN 505, CYD, RTS,S, COVID trials) the penalty
  should be on a placebo-anchored scale, such as the VE(s) curve or the
  marker × arm interaction, so the attack rate cancels by design. The HK and
  PNG cal-intercept shifts are attack-rate differences and should not be
  reported as penalties.
- **Asymmetric baseline.** The current penalty compares a CV estimate on a
  part of S with a full-S fit evaluated on T. The interaction fit uses both
  sides symmetrically, and the permutation null takes care of the CV noise.
- **Two directions are not two results.** A → B and B → A share the same
  data. Report one joint estimate (the interaction) plus the two directional
  rows as description.
- **Multiplicity.** Four metrics, three models and two directions give 24
  intervals per split. One pre-registered primary per dataset keeps a single
  significant row honest.
- **Ridge at the grid floor.** When CV picks the strongest penalty, the fit is
  intercept-only. Calibration slopes in that state (AMP: +422 [−3,024,
  +1,498]) should be suppressed, not printed.

## 6. Datasets and sources

Verified or already in hand are marked **[v]**. Everything else is a lead.

| Source | Type | Why it helps | Route |
| --- | --- | --- | --- |
| USG COVID correlates papers: COVE, ENSEMBLE (US/LatAm/South Africa), PREVENT-19, AZD1222 US | Summary: HR per 10-fold ID50 and bAb, per trial and region, in WHO IU | Cross-platform E1 at summary level, harmonised assays; strong correlates. ENSEMBLE-region and AZD1222 papers checked **[v]** | Published papers and source-data files |
| Coudeville et al. 2010 HAI meta-analysis | Summary: 15 studies, 5,899 adults, 1,304 cases, per-dataset curves **[v]** | τ for influenza HAI across strains and studies; already a random-effects model | Paper (PMC2851702) |
| RTS,S phase 3, MALARIA-055 (NCT00866619) | IPD, 11 sites, two age cohorts | Large E2 with site and age axes; anti-CS titre studied as a correlate. Listed on Vivli **[v]** (previous thread) | Add to Imo's Vivli request |
| ENSEMBLE COV3001 | IPD, three regions | E2 by region; ~92 breakthrough cases, but a strong correlate | YODA (doc/data-access.md) |
| Cross-vaccine COVID titre vs VE (Khoury 2021, Earle 2021, Cromer 2022) | Summary, 7–8 vaccines | Leave-one-vaccine-out curve error is a transfer penalty | Papers. This edges toward efficacy prediction, which the spec puts out of scope; use only as a transfer-error benchmark, never as a VE claim |
| Challenge studies (influenza, typhoid, cholera, CHMI) | IPD, small | Fixed exposure; known to translate poorly across studies | ImmPort, or authors; low priority |

## 7. What this spike did not do

- No AMP rerun (data behind the guestbook; the pipeline outputs already show
  no source signal there).
- No Cox version of attenuation; the same interaction approach applies.
- The HK permutation shuffles rows, not households, and the same people appear
  on both sides. Its SEs are household-clustered; its permutation p is
  approximate.
- None of section 6 has been extracted yet.
