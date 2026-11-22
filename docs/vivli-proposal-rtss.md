# Vivli research proposal: RTS,S/AS01 phase 3 (MALARIA-055)

**Status:** draft for Imo to review and submit. Nothing has been sent to Vivli
or GSK. Fields in `[[double brackets]]` need filling in before submission.
This is a **direct request for a listed study**, made in parallel with the CYD
enquiry in [vivli-proposal-cyd.md](vivli-proposal-cyd.md). The two share
methods, so sections 4, 6, 7 and 8 match that document.

Checked 2026-11-22. **[v]** marks facts confirmed on a live page or paper this
session; anything unmarked is from memory or is an estimate, and says so.

## 0. Read this first: route and feasibility

### The trial is listed

A keyword search of the Vivli catalogue for `NCT00866619` on 2026-11-22
returned the phase 3 trial, "Efficacy of GSK Biologicals' Candidate Malaria
Vaccine (257049) Against Malaria Disease Caused by P. Falciparum Infection in
Infants and Children in Africa" (15,459 enrolled), and its long-term extension
(3,084 enrolled) **[v]**. Both can go straight onto a Vivli data request. No
enquiry is needed. Request both: the extension carries follow-up and serology
beyond the main trial.

### How many cases have anti-CS titres

Anti-CS antibodies were measured by ELISA "in the first 200 participants in
each age category at each study site" (RTS,S Clinical Trials Partnership 2015,
*Lancet* 386:31) **[v]**. Measured at enrolment, 1 and 18 months after dose 3,
1, 12 and 24 months after the booster, and at the last visit **[v]**.

From the per-site, per-arm table in White et al. 2015 (*Lancet Infect Dis*
15:1450) **[v]**, the serology subset is about:

| | Control (C3C) | RTS,S, booster (R3R) | RTS,S, no booster (R3C) | RTS,S total |
| --- | --- | --- | --- | --- |
| Children 5–17 months | 810 | 801 | 811 | **1,612** |
| Infants 6–12 weeks | 876 | 873 | 890 | **1,763** |

These totals are our sum of that table, not figures the paper states.

**Neither paper reports how many children in the serology subset had a
clinical malaria episode.** Episode counts are given only for the full trial
populations (for example 8,352 episodes in per-protocol control children over
follow-up) **[v]**. Our estimate: clinical malaria incidence in the control
arm was high enough that, across 11 sites and about 3 to 4 years of follow-up,
well over 100 of the 1,612 vaccinated children and of the 1,763 vaccinated
infants in the subset will have had at least one episode, even within the
first year after dose 3. This is an estimate from the trial-wide incidence,
not a count. It is the first thing to check once the data are open, and the
analysis plan below falls back gracefully if a split is short.

Unlike CYD, where all cases had titres, RTS,S has titres only in this random
subset. So the case count is set by the subset size times the attack rate,
and per-site splits will be small at low-transmission sites.

### GSK requirements

GSK shares through Vivli (moved from ClinicalStudyDataRequest.com in 2023)
**[v]**. Vivli's own rules apply: a qualified statistician on the team, a
non-negotiable DUA, and analysis inside Vivli's secure research environment
**[v]**. `[[Check GSK's Vivli member page for any extra conditions, such as a
disease-area restriction like Sanofi's.]]`

## 1. Title and lay summary

**Title.** How far does the anti-circumsporozoite antibody correlate of
protection for RTS,S/AS01 carry across age groups and transmission settings?
Transfer penalties for correlate models within the phase 3 trial.

**Lay summary.** RTS,S/AS01 is the first malaria vaccine recommended for wide
use, and new malaria vaccines (such as R21/Matrix-M) are judged partly by how
their antibody responses compare with those of RTS,S. That comparison assumes
the relationship between anti-CS antibody levels and protection found in one
group of children holds in others. The RTS,S phase 3 trial gave the same
vaccine to young infants and older children at 11 sites with very different
malaria transmission. We will fit a model of malaria risk given antibody
levels in one group, apply it unchanged to another, and measure how much its
accuracy and calibration degrade. The study makes no claim about whether any
vaccine works and no new claim about malaria biology. Its output is a measured
size of the error introduced by borrowing a correlate, with uncertainty, and
open code others can reuse.

## 2. Background and significance

- Anti-CS titres tracked the magnitude and duration of RTS,S efficacy (White
  et al. 2015) **[v]**. Efficacy against clinical malaria was lower in infants
  than in children, and antibody titres were also lower in infants **[v]**.
  Whether the titre-risk relationship itself differs by age, or only the
  titres do, is exactly the question a transfer penalty answers.
- Second-generation malaria vaccines are compared with RTS,S on anti-CS
  (NANP) antibody levels. If the correlate does not carry across age or
  transmission intensity within one trial, that comparison needs a stated
  margin of error.
- Under a pure shift in baseline risk, such as different transmission
  intensity, our simulations show discrimination barely changes while
  calibration moves one-for-one. MALARIA-055, with sites ranging from low to
  very high transmission, is a direct test.
- Code: `[[repository URL: https://github.com/ramics/borrowed-light, if public
  by submission]]`.

## 3. Specific aims

**Aim 1 (primary): age transfer.** In RTS,S recipients in the serology subset,
fit a model of time to first clinical malaria episode from anti-CS titre at 1
month after dose 3 in children aged 5–17 months, apply it unchanged to infants
aged 6–12 weeks, and estimate the transfer penalty. Repeat in reverse.

**Aim 2: transmission-setting transfer.** Repeat Aim 1 between groups of sites
by transmission intensity (for example the higher- vs lower-incidence half of
sites, grouped by control-arm incidence), within each age category where
numbers allow.

**Aim 3: decision-level consequence.** For each transfer, estimate how often a
risk threshold set in the source gives a different classification than one
set in the target.

**Aim 4 (exploratory): time.** Use the 18-month and post-booster titres to ask
whether a model fitted on peak titres transfers to waned titres.

## 4. Analysis plan

**Population.** Per-protocol RTS,S recipients (R3R and R3C pooled for the
primary schedule) in the serology subset with a titre 1 month after dose 3.
Control-arm participants are used only for per-stratum attack rates.

**Endpoint.** Primary: first clinical malaria episode (temperature ≥37.5 °C
or reported fever, and *P. falciparum* asexual parasitaemia >5,000/µL; the
trial's primary case definition) **[v]** from 1 month after dose 3 to 12
months later. Analysed as time to event and as binary within the window.
Sensitivity: whole follow-up to month 20 (before the booster); all episodes.

**Markers.** log anti-CS (NANP) IgG at 1 month after dose 3. Secondary: anti-HBs
(the vaccine's carrier antigen), if measured in the subset, as a negative-control marker, and enrolment
anti-CS where measured.

**Models, weights, metrics, uncertainty.** As in the CYD proposal: ridge
logistic and ridge Cox with penalty chosen by cross-validation inside the
source only; best single marker as comparator; weighted AUC or Harrell's C,
calibration slope, calibration intercept, Brier score and decision-flip rate;
1,000-replicate participant-level bootstrap stratified by site and event
status, with 95% intervals everywhere. The serology subset was a fixed
first-enrolled sample per site and age category, so we treat it as a random
sample within site × age and weight by the inverse of its sampling fraction
in each stratum.

**Sample size.** About 1,612 vaccinated children and 1,763 vaccinated infants
in the serology subset **[v, our sum]**. Case counts are not published; we
expect well over 100 per side for Aim 1 (estimate). Aim 2 splits may fall
below 100 and are secondary.

**Pre-specification.** Primary = Aim 1, children → infants, ridge Cox,
calibration intercept and slope as headline metrics.

**Software.** Python 3 with our package, installed inside the secure research
environment. `[[Confirm SRE support.]]`

## 5. Data elements requested

From MALARIA-055 (NCT00866619) and its extension, for all randomised
participants:

| Domain | Variables |
| --- | --- |
| Design | Site, age category, randomisation arm, randomisation date, per-protocol and mITT flags, serology-subset flag |
| Demographics | Age at first dose, sex |
| Immunogenicity | Anti-CS (NANP) and anti-HBs titres at every visit measured, with visit dates and assay cut-off flags |
| Endpoints | Every clinical malaria episode with date and parasite density; severe malaria and hospitalisation flags; dates of active surveillance start and end, withdrawal or death |
| Exposure | Dose dates including booster |
| Context | Site-level transmission measure used in the trial, if any; bed-net use if recorded |

Not requested: names, exact dates of birth, free text, adverse-event
narratives, genetic or sequence data. Dates may be offset or given as study
days.

## 6. Data security and outputs

Same as the CYD proposal: analysis only inside Vivli's secure research
environment; no re-identification or linkage; aggregate outputs only, no cell
below `[[minimum cell size]]`; no participant data in the code repository.

## 7. Ethics

Secondary analysis of de-identified trial data. Non-human-subjects or exempt
determination from `[[institution]]`'s IRB before access.

## 8. Team, institution, funding, timeline, dissemination

As in the CYD proposal: `[[lead investigator]]`, `[[non-commercial requesting
institution]]`, `[[statistician]]`, `[[funding]]`, `[[conflicts]]`. Analysis
3 months from access, preprint within 6, journal submission within 9. Results
published whatever they show, crediting GSK and Vivli as the data source.

## Sources

- Vivli study search: https://search.vivli.org/ (keyword search for NCT00866619, 2026-11-22)
- Vivli request process: https://vivli.org/resources/requestdata/
- RTS,S Clinical Trials Partnership 2015, Lancet 386:31, https://pmc.ncbi.nlm.nih.gov/articles/PMC5626001
- White et al. 2015, Lancet Infect Dis 15:1450, https://pmc.ncbi.nlm.nih.gov/articles/PMC4655306
- GSK data sharing move to Vivli: https://www.gsk-studyregister.com/About_GSK_Patient_Level_Data_Sharing_Final_13July2023.pdf
