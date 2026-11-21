# Literature review: how far does a correlate of protection travel?

*Prepared 2026-11-21 as the first step of the borrowed-light spec (`SPEC.md`).
Scope follows the spec: methods for moving a correlate inference across studies,
platforms, populations and sample sizes. Nothing here characterises a pathogen or
predicts a vaccine's efficacy.*

Citation status: entries marked **[v]** were checked against an index or
publisher page while writing this. Unmarked entries are standard references
cited from memory; check page numbers before they go into a manuscript.

---

## 1. Summary for the build

1. **The field has the vocabulary but not the number.** Vaccine statistics
   already distinguishes a *specific* correlate (valid inside one trial) from a
   *general* or *bridging* correlate (valid across vaccines, populations or
   settings), and already has formal tools that take a "how imperfect is the
   transport?" quantity as an *input*: the Gilbert–Huang bridging sensitivity
   function and the two bias functions in the 2025 provisional-approval causal
   roadmap. Nobody routinely *estimates* that input from data. That is the gap
   this project fills, and it is the cleanest way to pitch it.
2. **The estimator already exists in another field.** Fit-on-source,
   evaluate-on-target with discrimination plus calibration is external
   validation of a clinical prediction model (Debray, Steyerberg, Van Calster,
   Riley). The spec should borrow its definitions wholesale, especially
   calibration-in-the-large and the calibration hierarchy, rather than
   re-derive them. Reviewers from that community will expect it.
3. **"Calibration breaks before discrimination" is well established, not new.**
   A shift in baseline risk that leaves the ranking alone cannot move AUC much
   and must move the intercept; Van Calster et al. call calibration "the
   Achilles heel of predictive analytics". The publishable part is the
   *vaccine* consequence: threshold-based admissibility decisions inherit the
   miscalibration. Influenza already supplies a real example (the adult HAI
   1:40 threshold gives ~22% protection in children; ~1:110 is needed for 50%).
4. **Be careful which scale you transport.** The spec's model predicts
   absolute risk in vaccinees, which is a *correlate of risk*. Absolute risk
   carries the attack rate of the study, so a large calibration-intercept
   penalty between two trials may only be saying "different epidemic", not
   "the correlate failed". Correlates-of-protection work transports vaccine
   efficacy as a function of the marker, VE(s), which is relative to placebo.
   Recommend adding a placebo-incidence offset (or a VE-curve comparison) so the
   decomposition can separate exposure from mechanism.
5. **Case-cohort weighting is mature and must be used as-is.** IPW and
   calibrated-weight estimators for Cox, logistic, AUC and VE curves exist, are
   what the CoVPN analyses used, and have public code. The bootstrap must
   respect the two-phase design.
6. **Small n is worse than the spike suggests.** Sample-size work on external
   validation says a calibration slope needs on the order of a hundred or more
   events to be estimated usefully. At n = 30 the honest output is "cannot
   tell", and the Bayesian second pass in the spec's open questions is
   justified.
7. **For E2, start today on public HVTN 505 data, then apply for CoVPN
   ENSEMBLE.** HVTN 505 ships in the CRAN `vaccine` package with sampling
   weights and stratifying covariates, so it exercises the full pipeline with
   no paperwork. ENSEMBLE (US / Latin America / South Africa, WHO-calibrated
   assays, published case-cohort SAP) is the best first *scientific* E2 result.
   Details in section 7.

No EC2 is needed at this stage: no data sits in `us-east-1` and the
synthetic suite runs on a laptop.

---

## 2. Definitions and frameworks for correlates

**Nomenclature.** Plotkin (2010) and Plotkin & Gilbert (2012, *Clin Infect Dis*
54:1615) separate a *correlate of risk* (a marker associated with the endpoint
in vaccinees) from a *correlate of protection* (a marker that reliably predicts
vaccine efficacy), and mechanistic from non-mechanistic correlates.

**Validation levels.** Qin, Gilbert, Corey, McElrath & Self (2007, *J Infect
Dis* 196:1304) set out a hierarchy: correlate of risk; surrogate of protection
at level 1 (predicts VE within the same setting); level 2 (predicts VE across
settings, i.e. "bridging"). Gilbert, Qin & Self (2008, *Stat Med*) gave the
statistical frameworks for each level. **[v]** Level 2 is exactly what E1
measures. The literature describes it but has few quantitative criteria for
*how much* level-2 validity a marker has.

**Four frameworks for assessing a correlate inside a trial** (reviewed in
Gilbert et al., *Vaccine* 2024 **[v]**):

| Framework | Key references | What it delivers |
| --- | --- | --- |
| Correlate of risk (Prentice-style regression in vaccinees) | Prentice 1989 *Stat Med*; CoVPN SAPs | Hazard ratio per marker unit; what the spec's ridge models are closest to |
| Principal stratification / principal surrogate | Frangakis & Rubin 2002; Gilbert & Hudgens 2008 *Biometrics*; Follmann 2006 | VE as a function of the marker the person *would* have had if vaccinated |
| Controlled effects / controlled VE curve | Gilbert, Fong, Kenny & Carone 2022 *Biostatistics* | VE(s): efficacy if everyone had marker level s |
| Mediation and stochastic interventions | Hejazi et al. 2021 *Biometrics* (two-phase); COVE mediation estimate 68.5% via Day-29 ID50 **[v]** | Proportion of VE explained; VE under shifted marker distributions |

**Meta-analytic (trial-level) surrogacy.** Daniels & Hughes (1997, *Stat Med*);
Buyse, Molenberghs et al. (2000, *Biostatistics* 1:49) define trial-level R²
between treatment effects on the surrogate and on the outcome across many
trials. This is the established answer to "does it transfer across studies?",
but it needs many trials and says nothing about individual-level calibration.
It is the natural comparator if the project ever has three or more studies
(the spec's leave-one-study-out case).

**Surrogate indices.** Price, Gilbert & van der Laan (2018, *Biometrics*) on
the optimal surrogate from a randomised trial; Athey, Chetty, Imbens & Kang
(2019) on the surrogate index. A multi-marker ridge model is a surrogate index
in this sense, which is a useful framing for the write-up.

## 3. Bridging and transporting a correlate

**Formal transport with a sensitivity function.** Gilbert & Huang (2016,
*Epidemiologic Methods*, PMC5279575) predict overall VE in a new setting by
re-weighting baseline covariates and marker distributions, with a
user-specified *bridging assumption function* for how much the conditional VE
differs between settings. **[v]**

**Provisional approval roadmap.** Gilbert, Peng, Han, … Follmann, Juraska &
Díaz (2025, *Biostatistics* 26(1) kxaf018; arXiv 2407.06350), illustrated with
Group B Streptococcus. Two bias functions: one for the surrogate index being an
imperfect surrogate in the phase 3 study, one for *imperfect transport* of the
surrogate index. The output is a conservative lower bound on the treatment
effect. **[v]** *This paper is the closest neighbour to borrowed-light.* Its
transport bias function is user-supplied; borrowed-light's penalty
distributions, estimated across real study pairs, are a candidate empirical
prior for it. Worth contacting the authors.

**Cross-trial immunogenicity comparison.** "Comparing HIV vaccine
immunogenicity across trials with different populations and study designs"
(arXiv 2410.05594) applies transportability estimators to HVTN trials. **[v]**
A direct methodological template for E1 marker alignment.

**General causal transportability.** Pearl & Bareinboim (2014, *Stat Sci*
29:579); Dahabreh et al. (2020, *Stat Med*) on extending trial inferences to a
target population. These give the identification conditions (which variables
must be measured in both settings) that the marker-mapping file implicitly
asserts.

## 4. What has actually happened when correlates were moved

These are the motivating cases for the write-up. Each is a natural experiment
in transfer.

**COVID-19, the best-instrumented case.** The USG/CoVPN correlates programme
ran harmonised case-cohort analyses with a shared SAP and WHO-calibrated assays
across platforms:

- mRNA-1273 (COVE): Gilbert et al. 2022, *Science* 375:43. **[v]**
- Ad26.COV2.S (ENSEMBLE): Fong et al. 2022, *Nat Microbiol* 7:1996. Compares
  VE-by-ID50 curves across ENSEMBLE-US, COVE and COV002-UK (AZD1222) and
  reports consistency of ID50 as a correlate across trials and platforms. **[v]**
- AZD1222 (US/LatAm): Benkeser et al. 2023, *npj Vaccines*. **[v]**
- NVX-CoV2373 (PREVENT-19): CoVPN/Novavax correlates team, *Nat Commun*, January 2023. **[v journal; authors not checked]**
- COV002 (Oxford): Feng et al. 2021, *Nat Med* 27:2032.
- Programme overview: "A Covid-19 milestone attained — a correlate of
  protection for vaccines", *NEJM* 2022;387:2203. **[v citation; authors not checked]**

Pooled summary-level work: Khoury et al. 2021 (*Nat Med*, 50% protection at
~20% of mean convalescent titre) **[v]**; Earle et al. 2021 (*Vaccine*, seven
vaccines) **[v]**; Cromer et al. 2022 (*Lancet Microbe*, variants and boosting).
Khoury et al. 2023 (*Emerg Infect Dis* 29(2), with Fong and Gilbert) show that
reported protective thresholds disagree across studies until titres are
normalised to a common reference, after which they converge. **[v]** This is
direct evidence that the spec's "assay" dial is real and dominant.

Caveat: "consistent across platforms" in these papers is a visual or
confidence-band statement about VE curves, not a quantified transfer penalty.
That is the opening.

**HIV, the classic non-transfer.** RV144 (Thailand, ALVAC/AIDSVAX B/E)
identified IgG V1V2 as an inverse correlate of risk (Haynes et al. 2012,
*NEJM*). HVTN 702 (South Africa, clade C regimen with MF59) showed no efficacy.
Its correlates analysis (Moodie et al. 2022, *J Infect Dis* 226:246) found a
different immune profile (higher CD4+, lower V1V2 response rate) but still
found V1V2 and polyfunctional CD4+ associations with acquisition. **[v]** So the
*correlate of risk* partly travelled while the *efficacy* did not. Good evidence
for point 4 of the summary.

**HIV bnAbs.** The AMP trials (HVTN 704/HPTN 085 in the Americas and Europe,
HVTN 703/HPTN 081 in sub-Saharan Africa) tested the same product in two
populations; predicted serum neutralisation titre (PT80) tracked prevention
efficacy (AMP correlates analyses, Gilbert et al. 2022, *Nat Med*; trial
results Corey et al. 2021, *NEJM*). **[v trial design and PT80 finding]** One product, two populations:
an E1 design with the platform held fixed.

**Dengue.** CYD14 (Asia, ages 2–14) and CYD15 (Latin America, ages 9–16) were
harmonised twin trials; higher month-13 neutralising titres associated with
higher VE in each (Moodie et al. 2018, *J Infect Dis*, PMC5854020). **[v]**
Serostatus at baseline changes everything, which makes it a natural E2 split.

**Influenza.** The HAI 1:40 threshold (Hobson 1972) is the textbook
borrowed correlate. Coudeville et al. (2010, *BMC Med Res Methodol*) estimated
a 50% protective titre around 1:17 in adults with a Bayesian meta-analysis.
**[v]** Black et al. (2011, *Pediatr Infect Dis J* 30:1081) found 1:40 gave
~22% protection in children and ~1:110 was needed for 50%. **[v]** This is a
calibration failure on transfer across age, with discrimination presumably
intact: the spec's AC3 in the wild.

**Tuberculosis.** MVA85A infant case-control correlates (Fletcher et al. 2016,
*Nat Commun* 7:11290): activated HLA-DR+ CD4+ T cells raise risk; no immune
variable reached AUROC 0.8. **[v]** M72/AS01E: phase 2b 49.7% efficacy; the
phase 3 trial (~20,000 participants) is explicitly intended to validate
correlates from the 2b samples, through the Gates MRI TB Immune Correlates
Program. **[v]** Phase 3 readout is years away, so M72 is a later E1 target,
not a first one.

**Regulatory immunobridging.** ICMRA and the Access Consortium accept
immunobridging (non-inferiority on neutralising GMT ratio and seroresponse,
WHO reference standards) to authorise new COVID-19 vaccines. **[v]** Example:
Taiwan's MVC-COV1901 authorisation used a GMT-ratio lower bound > 0.67 versus
AZD1222. **[v]** These margins are fixed by convention, with no estimate of
the transfer penalty behind them. A decision-flip metric (spec, step 6) speaks
directly to whether such a margin is safe.

## 5. Prediction-model transfer and dataset shift

This is where the estimator's definitions should come from.

- **Reproducibility vs transportability.** Debray et al. (2015, *J Clin
  Epidemiol* 68:279) score how related the validation case-mix is to the
  development case-mix, so a validation result can be read as either. **[v]**
  Recommend reporting a case-mix relatedness statistic alongside every
  penalty, so readers can see how far each transfer stretched.
- **Calibration hierarchy.** Van Calster et al. (2016, *J Clin Epidemiol*
  74:167): mean, weak, moderate and strong calibration. "Calibration: the
  Achilles heel of predictive analytics" (Van Calster et al. 2019, *BMC Med*).
- **Calibration-in-the-large is estimated with the slope fixed at 1** (logistic
  model with the linear predictor as an offset), and the slope is estimated
  separately (Steyerberg, *Clinical Prediction Models*, 2nd ed. 2019). The
  spike fits both jointly, which gives a different intercept whenever the
  slope is not 1. Report both forms, and label them.
- **Recalibration** (intercept-only, then intercept + slope, then refit) is the
  standard ladder for target-specific updating. The step at which a borrowed
  model becomes usable is itself a decision-level summary worth reporting.
- **External-validation sample size.** Riley et al. (2021, *Stat Med*
  40:4230) give minimum sizes for precise C, calibration slope and O/E; Snell
  et al. argue rules like "100 events" are too crude and favour simulation
  **[v Snell]**. This sets realistic expectations for AC4 and argues for adding a
  "target too small to estimate" flag rather than a wide interval alone.
- **Transporting performance without target outcomes.** Steingrimsson,
  Gatsonis, Li & Dahabreh (2022, *Am J Epidemiol* 192:296) identify
  performance measures in a target with covariates but no outcomes; Morrison
  et al. (2024, *Can J Stat*) give doubly robust versions. **[v]** Useful as a
  *predicted* penalty to compare against the *observed* one. If the two agree,
  the penalty can be forecast before an outbreak trial has outcomes, which is
  the project's motivating use.
- **Shift taxonomy.** Covariate shift, prior (label) shift and concept shift
  (Quiñonero-Candela et al. 2009; Moreno-Torres et al. 2012). Prior-shift
  correction by re-estimating the base rate: Saerens, Latinne & Decaestecker
  (2002, *Neural Comput* 14:21); Lipton, Wang & Smola (2018, ICML). The spec's
  four dials map onto this: population ≈ prior/intercept shift, assay and
  noise ≈ covariate shift with measurement error, mechanism ≈ concept shift.
  Using the standard names will help reviewers. Clinical framing: Finlayson et
  al. 2021, *NEJM* 385:283, "The clinician and dataset shift in AI".

## 6. Case-cohort and two-phase sampling

- **Designs and estimators.** Prentice (1986, *Biometrika* 73:1) case-cohort;
  Barlow (1994) robust variance; Borgan et al. (2000) stratified case-cohort;
  Breslow, Lumley et al. (2009, *Am J Epidemiol*) calibrated weights using the
  whole cohort; Saegusa & Wellner (2013, *Ann Stat*) weighted likelihood with
  estimated and calibrated weights **[v]**.
- **Weighted AUC.** IPW AUC estimators for two-phase case-control designs,
  with weights estimated within sampling strata and improved by phase-one
  auxiliaries (Huang et al. 2016, *Biostatistics*, PMC4915610; follow-up in
  *Stat Med* 2019, PMID 30209824). **[v]** The naive empirical AUC "can be
  badly biased" under two-phase sampling.
- **What CoVPN did.** Every COVID correlates paper above sampled a random
  immunogenicity subcohort stratified by baseline covariates and added all
  breakthrough cases, then used IPW with estimated weights. **[v]** Analyses
  were restricted to baseline-negative per-protocol vaccinees.
- **Implications for the build:**
  1. `sampling_weight` should be the inverse of the phase-two sampling
     probability *within strata*, and the schema needs the sampling-stratum
     variable too, because the bootstrap and estimated weights both need it.
  2. Bootstrap by resampling the phase-one cohort within strata and keeping
     each resampled person's phase-two status, or the intervals will be too
     narrow.
  3. Weighted calibration slope and intercept are just weighted logistic
     regressions on the linear predictor. Weighted Brier is a weighted mean.
     Weighted C for time-to-event needs the IPW form of Harrell's or Uno's C.
  4. AC5 can be tested on real public data (HVTN 505, below) as well as
     synthetic data.

## 7. Candidate participant-level datasets for E2

Ranked by how quickly they can produce a first real result. See `docs/data-access.md` for the 2026-11-21 public-data search and step-by-step access routes.

| Dataset | Access | Why it fits E2 | Limitations |
| --- | --- | --- | --- |
| **HVTN 505** (DNA/rAd5 HIV vaccine) | **Public now**, `hvtn505` in the CRAN `vaccine` package (from SCHARP Atlas). About 2,300 rows, IgG env/V2/V3 markers, age, BMI, behavioural risk score, case-control indicator, IPW weights. **[v]** | Zero paperwork. Has real two-phase weights, so AC5 runs on real data. Strata by age, BMI or risk score. | Null efficacy, so it supports correlate-of-risk transfer only. Few cases, so expect wide intervals. Best used to validate plumbing. |
| **ENSEMBLE** (Ad26.COV2.S) | CoVPN proposal process ("submit vaccine study proposals"); data held by CoVPN/SCHARP. **[v]** No public route found. | Three regions (US, Latin America, South Africa) with different variants. WHO-calibrated assays. Case-cohort with published SAP and code. The ideal E2: same trial, real population and variant shift. | Needs a proposal and DUA. Timeline unknown. |
| **COVE, PREVENT-19, AZD1222 US** | Same CoVPN route | Same harmonised pipeline; COVE + ENSEMBLE + AZD1222 is a ready-made E1 with three platforms, which also enables leave-one-study-out. | Same paperwork. |
| **CYD14 / CYD15** (dengue) | Sponsor data-sharing portal (likely Vivli or similar; not verified) | Twin harmonised trials; serostatus and age are natural strata; the pair is also an E1 with platform fixed. | Sponsor approval; serostatus complicates the interpretation. |
| **AMP HVTN 703/704** | HVTN/HPTN request process | Same product in two populations; time-to-event. | bnAb, not a vaccine; few infections per arm. |
| **RV144 / HVTN 702** | CAVD DataSpace (sign-in required) and SCHARP Atlas **[v]** | Classic cross-population pair. | Case-control with small samples; gated. |
| **M72/AS01E** | VISC / Fred Hutch terms, as the spec notes | The parent-programme target. | Phase 3 has not read out. |

**Recommendation.** Use HVTN 505 for build step 4 now, so E2 exists on real
data and AC5 is tested against real weights. In parallel, draft a CoVPN
proposal for ENSEMBLE region-stratified E2, which can grow into a three-trial E1.
Getting that proposal started early matters more than any modelling choice,
since the spec names data access as the binding constraint.

## 8. Notes on the spike, from this reading

- The spike's `bootstrap_interval` re-simulates fresh source and target data
  for each replicate. That is a sampling-distribution interval under a known
  data-generating process, not a participant-level bootstrap of one fixed
  dataset. Real data needs the latter, and AC4 should be re-checked with it,
  because the bootstrap is known to be unreliable at n = 30.
- The calibration intercept is estimated jointly with the slope (section 5).
- The AC3 result is expected from theory; frame it as a known property whose
  regulatory consequence is under-appreciated, and cite the influenza example.
- Consider adding a fifth dial for **attack rate** (exposure) separately from
  baseline susceptibility, since this is the difference between "different
  people" and "different epidemic" and section 4 suggests it dominates
  real cross-trial gaps.

## 9. Open questions this review raises

1. Is the primary estimand transfer of absolute risk (correlate of risk) or of
   VE(s) (correlate of protection)? The second needs placebo arms in both
   studies but is what regulators use.
2. Should the project report a predicted penalty (Steingrimsson-style, no
   target outcomes) next to the observed one? That would make the work usable
   before an outbreak trial has endpoints.
3. Would the Gilbert et al. roadmap authors take borrowed-light penalties as
   empirical input for their transport bias function? A short email could
   settle the framing of the whole write-up.
