# Vivli research proposal: dengue CYD14 and CYD15

**Status:** draft for Imo to review and submit. Nothing has been sent to Vivli
or Sanofi. Fields in `[[double brackets]]` need filling in before submission.

Checked 2026-11-22. **[v]** marks facts confirmed on a live page or paper this
session; anything unmarked is from memory and should be confirmed.

## 0. Read this first: route and feasibility

### CYD14 and CYD15 are not listed on Vivli

A keyword search of the public Vivli study catalogue (search.vivli.org) on
2026-11-22 found no listing for `NCT01373281` (CYD14), `NCT01374516` (CYD15),
"CYD" or "CYD-TDV" **[v]**. "dengue" returns only two studies, neither of them
CYD14 or CYD15: an investigator-led Philippines Dengvaxia serostatus study
(NCT03465254) and a Takeda TDV phase 2 trial **[v]**. Sanofi is a Vivli member,
so the route is still Vivli. Vivli handles unlisted studies through its
**enquiry form for studies not listed on the platform**, which Sanofi reviews
case by case **[v]**. So the order is:

1. Send the enquiry (section 9 holds the text) and ask Sanofi to confirm that
   the CYD14 and CYD15 packages would include the month 0, 7, 13 and 25 PRNT50
   values for the immunogenicity subset and all cases.
2. If Sanofi agrees to list the trials, submit the full data request (sections
   1 to 8) through the Vivli request form.

### Sanofi's sharing criteria, checked against CYD

From Sanofi's Vivli member page **[v]**:

| Criterion | CYD14 / CYD15 |
| --- | --- |
| Interventional phase II–IV study | Phase 3. Meets it. |
| Product approved on or after 1 Jan 2010 by US and/or EU, or development terminated on or after that date | Dengvaxia: EU approval 2018, FDA approval 2019 (from memory). Meets it. |
| Trial completed and primary results accepted for publication | Both published (Capeding 2014 *Lancet*; Villar 2015 *NEJM*), and long-term follow-up is complete (from memory). Meets it. |
| Exceptions (single site, rare disease, legal limits) | Neither is single-site or a rare disease. |
| **Pre-2014 consent: third-party use must address a scientific question in the same disease area** | **Both trials enrolled in 2011, so this applies.** The proposal is framed as a dengue question: how far the dengue neutralising-antibody correlate carries across the populations in these two trials. Section 2 is written to satisfy it. |
| PK data not routinely provided | Not needed. |

Vivli also requires **a qualified statistician on the team** (it is the first
administrative check), a **non-negotiable DUA** signed by the lead investigator
and any team member from another institution, and analysis inside Vivli's
**secure research environment** (free standard tier for 365 days, then
$12/day; non-standard set-ups such as Docker cost extra) **[v]**.

### Case counts: CYD clears ~100 per side on the trial split

From Moodie et al. 2018, *J Infect Dis* 217:742 (PMC5854020) **[v]**:

| | CYD14 (Asia, 2–14 y) | CYD15 (Latin America, 9–16 y) |
| --- | --- | --- |
| DENV-Any cases, month 13–25, with month 13 PRNT50 | 244 (115 vaccine, 129 placebo) | 415 (183 vaccine, 232 placebo) |
| Share of all DENV-Any cases that had month 13 titres | 99.6% | 99.8% |
| Non-case immunogenicity subset with month 13 titres | 1,879 (1,275 vaccine, 604 placebo) | 1,884 (1,275 vaccine, 609 placebo) |

- **Trial split (CYD14 ↔ CYD15): 115 vs 183 vaccinee cases.** Both sides clear
  ~100. This is the primary analysis.
- **Age split within CYD14:** not reported by stratum in the paper. With 115
  vaccinee cases in total, a 2–8 vs 9–14 split will fall below 100 per side.
  Secondary.
- **Pooled age split** (CYD14 2–8 vs CYD14+CYD15 9–16): the 9–16 side will be
  well over 100; the 2–8 side is likely 50–80. Secondary.
- **Baseline serostatus split: limited.** Month 0 titres exist only for the
  random immunogenicity subset, so 86.6% of cases have no baseline sample
  **[v]**. A serostatus split needs inferred serostatus for cases, such as the
  anti-NS1 IgG at month 13 that Sanofi later measured for the Dengvaxia
  serostatus re-analysis (Sridhar et al. 2018 *NEJM*; from memory). We ask for
  it; if it is not in the package, this split is subset-only and exploratory.

Design: the immunogenicity subset was the participants enrolled in the first 2
to 4 months, randomly assigned **[v]**, and titres were measured in all cases.
That is a case-cohort design. Sampling weights are not in the paper, so we ask
for the subset flag and the phase-one counts needed to rebuild them.

### If Sanofi declines: RTS,S phase 3 is listed

GSK's RTS,S/AS01 phase 3 trial (MALARIA-055, NCT00866619, 15,459 enrolled) and
its extension **are** listed on Vivli **[v]**, so they can be requested
directly. Anti-CS antibody titres tracked protection (White et al. 2015
*Lancet Infect Dis*), and the trial has two age cohorts and 11 sites. How many
clinical malaria cases have anti-CS titres is not yet checked; titres were
measured in an immunogenicity subset. Sanofi's high-dose influenza efficacy
trial (FIM12, NCT01427309) is not listed **[v]**.

## 1. Title and lay summary

**Title.** How far does the dengue neutralising-antibody correlate of
protection carry across populations? Transfer penalties for correlate models
between the CYD14 and CYD15 trials.

**Lay summary (about 200 words).** When a new vaccine cannot be tested for
efficacy, regulators sometimes accept immune-response markers instead,
reasoning that a response which protected people in an earlier trial will
protect people now. That reasoning assumes that the relationship between the
marker and protection carries over from the trial where it was measured to a
new population. How much it weakens when it is carried over has rarely been
measured. The two phase 3 trials of the CYD-TDV dengue vaccine, CYD14 in five
Asian countries and CYD15 in five Latin American countries, used the same
vaccine, schedule, assay and endpoint in different populations, ages and
dengue epidemiology. That makes them close to an ideal test. We will fit a
model of dengue risk given month-13 neutralising antibody titres in one trial,
apply it unchanged to the other, and measure how much the model's accuracy
and calibration degrade. We will repeat this across age groups and, where the
data allow, baseline serostatus. The study makes no claim about whether any
vaccine works and no new claim about dengue biology. Its output is a measured
size of the error introduced by borrowing a correlate, with uncertainty, and
openly published code that others can apply to their own trials.

## 2. Background and significance

- Neutralising antibody titres measured by PRNT50 at month 13 were associated
  with lower dengue risk and higher vaccine efficacy in both CYD14 and CYD15
  (Moodie et al. 2018) **[v]**. But the size of the association differed
  between trials. The tertile cut points for a "medium" response were 58–266 in
  CYD14 and 135–631 in CYD15 **[v]**, and baseline serostatus changed the
  picture.
- Dengue is where a borrowed correlate went wrong in public. Efficacy and risk
  depended on prior exposure in a way an antibody threshold alone did not
  capture, which led to the 2017–2018 serostatus restrictions on Dengvaxia.
  Measuring how much a titre-based risk model degrades between the two trials,
  and between age and serostatus groups, is a dengue question with direct
  relevance to second-generation dengue vaccines (TAK-003, Butantan-DV), whose
  correlate analyses will again be borrowed across populations.
- Correlate transfer is usually discussed in terms of discrimination (does the
  marker still rank people?). Our simulation work shows that under a pure
  population shift, discrimination is almost unchanged while calibration moves
  one-for-one with the shift. Since accepting a correlate is a threshold
  decision, calibration is what matters. The CYD trials, with different attack
  rates in the same vaccine, are a direct test of this.
- The analysis code is open, tested against simulations with known answers,
  and already run on public HIV (HVTN 505, AMP) and malaria datasets:
  `[[repository URL: https://github.com/ramics/borrowed-light, if it will be public by submission]]`.

## 3. Specific aims

**Aim 1 (primary): trial-to-trial transfer.** In vaccine recipients, fit a
model of DENV-Any risk from month-13 PRNT50 titres (four serotypes) in CYD14,
apply it unchanged to CYD15, and estimate the transfer penalty in
discrimination, calibration slope, calibration intercept and Brier score.
Repeat in the reverse direction.

**Aim 2: within-trial population transfer.** Repeat Aim 1 across (a) age
groups (CYD14 2–8 y vs 9–14 y; CYD14 2–8 y vs pooled 9–16 y), and (b) baseline
serostatus, using month-0 PRNT50 in the immunogenicity subset and inferred
serostatus for cases if available.

**Aim 3: decision-level consequence.** For each transfer in Aims 1–2, estimate
how often a risk threshold set in the source population gives a different
classification (above or below a target risk) than one set in the target
population.

**Aim 4 (exploratory): which part of the gap is which.** Use placebo-arm
attack rates per stratum to separate a shift in baseline risk from a change in
the titre-risk relationship.

## 4. Analysis plan

**Population.** Per-protocol vaccine recipients in each trial, phase two =
immunogenicity subset plus all DENV-Any cases with month-13 titres. Placebo
recipients are used only for per-stratum attack rates (Aim 4) and as a
sensitivity analysis.

**Endpoint.** Primary: symptomatic virologically confirmed dengue of any
serotype (DENV-Any) from month 13 to month 25, the trials' primary efficacy
endpoint. Analysed as binary and as time to event from month 13. Secondary:
serotype-specific endpoints, hospitalised dengue (likely too few events; will
report counts only if so).

**Markers.** Month-13 log10 PRNT50 for DENV-1 to DENV-4 (values below LLOQ of
10 set to 5 as in the original analysis), and their geometric mean. Month 0
titres as covariates in the serostatus analyses only.

**Models.**
- Ridge-penalised logistic regression (binary endpoint) and ridge-penalised
  Cox regression (time to event).
- The best single marker as a comparator. If the multivariate model does not
  beat it, we report that.
- Penalty strength chosen by cross-validation **within the source population
  only**. The target population is never used to fit or tune anything; the
  code enforces this.

**Sampling weights.** Inverse-probability weights for the case-cohort design,
rebuilt from the immunogenicity-subset flag and the phase-one counts by trial,
arm and age stratum. All fitting and all metrics are weighted.

**Transfer operation.** For a source stratum S and target T:
fit on S; evaluate on held-out S (cross-validated where S is small) for the
home baseline; evaluate on T; penalty = baseline − transfer, per metric.

**Metrics.** Weighted AUC (Harrell's C for time to event), calibration slope,
calibration intercept (calibration-in-the-large), Brier score, and the
decision-flip rate at fixed and target-optimal thresholds.

**Uncertainty.** Participant-level stratified bootstrap (1,000 replicates,
stratified by case status and sampling stratum), 95% percentile intervals for
every penalty. No point estimate is reported without an interval.

**Sample size.** Aim 1 has 115 (CYD14) and 183 (CYD15) vaccinee cases with
month-13 titres. From simulation and from our runs on public data, about 100
cases per side is where calibration-slope intervals become narrow enough to
be informative; at about 50 per side (AMP HIV trials) the calibration
intercept interval was about ±0.6 logit. Aims 2a–2b are expected to fall
below 100 per side and are labelled secondary, with intervals reported as
they come out.

**Pre-specification.** Primary analysis = Aim 1, CYD14 → CYD15, ridge
logistic, DENV-Any, calibration intercept and slope as the headline metrics.
Everything else is secondary or exploratory and labelled as such.

**Software.** Python 3 (numpy, scipy, pandas, scikit-learn) plus our package,
installed from source inside the secure research environment.
`[[Confirm with Vivli that the standard SRE allows Python and installing a
local package; R is the fallback.]]`

## 5. Data elements requested

From the CYD14 (NCT01373281) and CYD15 (NCT01374516) SDTM/ADaM packages, for
all randomised participants:

| Domain | Variables |
| --- | --- |
| Design | Trial, country, site (or site group), randomisation arm, randomisation date, per-protocol and full-analysis flags, immunogenicity-subset flag |
| Demographics | Age at enrolment (years, or the trial's age strata), sex |
| Immunogenicity | PRNT50 titres for DENV-1 to DENV-4 at months 0, 7, 13 and 25, with visit dates and LLOQ flags |
| Serostatus | Baseline serostatus as defined by Sanofi; any month-13 anti-NS1 IgG or other inferred serostatus measured on cases and the subcohort |
| Endpoints | Symptomatic VCD (DENV-Any) with onset date and serotype; hospitalised and severe dengue flags; end of active surveillance date or date of loss to follow-up |
| Exposure | Dose dates (to compute time from dose 3) |

Not requested: names, exact dates of birth, free-text fields, adverse event
narratives, PK data, sequence data, or anything allowing re-identification.
Dates may be offset or replaced by study days.

## 6. Data security and outputs

- Analysis runs only inside Vivli's secure research environment. No
  participant-level data leaves it.
- No attempt at re-identification and no linkage to other datasets.
- Exported outputs are aggregate only: penalty estimates, intervals, and
  summary tables with no cell below `[[Vivli/Sanofi minimum cell size]]`.
- No participant data is ever committed to the public code repository; the
  repository has a test and a `.gitignore` backstop for this.

## 7. Ethics

Secondary analysis of de-identified trial data. We will obtain a
non-human-subjects or exempt determination from
`[[institution]]`'s IRB before data access. `[[Attach determination letter.]]`

## 8. Team, institution, funding, timeline, dissemination

- **Lead investigator:** `[[name, degrees, position]]`
- **Requesting institution (non-commercial):** `[[institution name and
  address]]`. Vivli and Sanofi require a non-commercial academic, research,
  healthcare or government institution to sign the DUA.
- **Statistician (required by Vivli):** `[[name, qualifications]]`
- **Other team members:** `[[names, institutions; each from another
  institution signs a DUA]]`
- **Funding:** `[[source, or "unfunded"]]`
- **Conflicts of interest:** `[[ICMJE disclosure; state no commercial
  interest in dengue vaccines]]`
- **Timeline:** analysis 3 months from data access; preprint within 6 months;
  journal submission within 9 months.
- **Dissemination:** preprint and peer-reviewed paper, crediting Sanofi and
  Vivli as the data source in the form they require. Code public throughout.
- **Commitment to publish:** results will be published whatever they show,
  including a finding that the correlate transfers well.

## 9. Enquiry text for the unlisted studies

> We would like to request individual participant data from Sanofi's
> CYD14 (NCT01373281) and CYD15 (NCT01374516) phase 3 dengue vaccine trials,
> which we could not find on the Vivli platform. Our non-commercial study
> asks how well the month-13 neutralising-antibody correlate of risk for
> dengue carries over between the two trials, and between age and
> serostatus groups within them (methods follow Moodie et al., J Infect Dis
> 2018). We need, for all randomised participants: arm, age, sex, country,
> per-protocol and immunogenicity-subset flags, PRNT50 titres for the four
> serotypes at months 0, 7, 13 and 25, symptomatic virologically confirmed
> dengue with onset date and serotype, and follow-up time. Could you confirm
> (1) whether Sanofi would consider listing these trials, and (2) whether the
> PRNT50 values for the immunogenicity subset and all cases, and any
> inferred baseline serostatus for cases (such as month-13 anti-NS1), would
> be included in the data package? The requesting institution is
> `[[institution]]`.

## Sources

- Vivli study search: https://search.vivli.org/ (keyword searches run 2026-11-22)
- Sanofi on Vivli: https://vivli.org/ourmember/sanofi/
- Vivli request process and costs: https://vivli.org/resources/requestdata/
- Vivli enquiry for unlisted studies: https://vivli.org/members/enquiries-about-studies-not-listed-on-the-vivli-platform/
- Moodie et al. 2018, J Infect Dis 217:742, https://pmc.ncbi.nlm.nih.gov/articles/PMC5854020
- White et al. 2015, Lancet Infect Dis 15:1450 (RTS,S anti-CS antibody and protection)
