# Data access: what is public, and how to get the rest

Checked 2026-11-21. Items marked **[v]** were confirmed on the live source page
or API this session; everything else is from papers or search results and
should be confirmed before relying on it. No participant data is committed to
this repository, and none should be.

## 1. Bottom line

- **No public dataset reaches ~100 events per side with markers.** The largest
  public option is the AMP trials (HVTN 704 and HVTN 703), with 60 and 47
  infections in the VRC01 arms. That is roughly double HVTN 505 per side, enough
  for a first E2 that is not purely plumbing, but intervals will still be wide.
- **The best gated route changed.** ENSEMBLE (VAC31518COV3001, 44,325
  enrolled) is now listed on the YODA Project as an available trial with
  "Collected Datasets Available" **[v]**. That is a standard sponsor
  data-sharing route, separate from a CoVPN proposal. The open question is
  whether the Day 29 antibody markers from the correlates case-cohort are in the
  J&J datasets. Ask YODA that before writing a full proposal.
- **One structural check first.** YODA and the HVTN Dataverse licence both
  exclude commercial use, and YODA requires the PI to sit at an academic,
  research, healthcare or government institution that co-signs the DUA. Settle
  which institution is the requester before anything is submitted.

## 2. What a usable dataset needs (from SPEC.md)

1. Participant-level marker values (`marker_*`), on one assay per marker.
2. Binary endpoint, ideally with time to event.
3. Case-cohort or case-control sampling weights, or enough design information
   to rebuild them (subcohort flag, strata, phase-one counts).
4. A split variable for E2 (region, age, sex, risk group, trial).
5. Placebo-arm event counts per stratum for the attack-rate dial.
6. For E1, a second study with the same or a mappable assay.

## 3. Public or registration-only datasets

| Dataset | Where and how | Size that matters | Fits | Catch |
| --- | --- | --- | --- | --- |
| **AMP: HVTN 704/HPTN 085 + HVTN 703/HPTN 081** (VRC01 bnAb infusions) | Harvard Dataverse, HVTN collection **[v]**. Primary results DOI `10.7910/DVN/GDESJI`, PT80 biomarker `10.7910/DVN/50VH2C`, week-104 extension `10.7910/DVN/WEJPUF`. Free download after a guestbook form (name, email, institution, HVTN affiliation) **[v]**. Licence CC BY-NC-ND 4.0 **[v]**. | Infections: 704 placebo 38, VRC01 60; 703 placebo 29, VRC01 47 (Corey et al. NEJM 2021). | **Best public E2.** Split by trial: men and transgender people in the Americas and Europe vs women in sub-Saharan Africa. Same product, two populations. Files include survival data, predicted daily VRC01 concentration, PT80 at infection for cases, and case-control sets. Two doses give a second axis. | Passive antibody, not a vaccine. The marker is a PK-model prediction largely set by dose and body weight. Protection only against sensitive viruses, so the endpoint choice matters. |
| **HVTN 702** (ALVAC-HIV + bivalent gp120, South Africa) | Harvard Dataverse `10.7910/DVN/JI2FUM` **[v]**: `cor.tab` (28 MB), `endpoints.tab`, data dictionary, CoR SAP, protocol. Same guestbook and licence **[v]**. | 60 vaccinee cases and 60 matched non-cases, women only (Moodie et al. JID 2022). | **E1 with a known answer.** It tested the RV144 correlates (IgG V1V2, IgA Env, polyfunctional CD4) on purpose, so a marker mapping already exists. RV144 to HVTN 702 is the textbook correlate that did not transfer. | Small. Case-control, not case-cohort. Needs RV144 (next row). |
| **RV144** (ALVAC-HIV + AIDSVAX B/E, Thailand) | CAVD DataSpace, Haynes et al. NEJM 2012 publication data. Free registration and sign-in. R access via `DataSpaceR`. | 41 vaccinee cases, 205 controls. | Source side of the RV144 to 702 E1. | Small. Registration terms not yet read. |
| **HVTN 505, more markers** | Harvard Dataverse `10.7910/DVN/J7NBBN` (2018 JID: BAMA IgG/IgA + T cells) and `IGHFOI` (2017 JID T cells) **[v]**. | Same 25 vaccinee cases we already have. | Adds markers to the existing E2 run. | Does not fix the case count. |
| **Hong Kong influenza B, children** (TIV vs placebo) | Zenodo record 4972891, CC0 **[v]**. Downloaded and checked: 796 children (479 TIV, 317 placebo), pre- and post-vaccination HAI, time to RT-PCR-confirmed influenza **[v]**. | 31 infections **[v]**. | Clean, tiny, CC0. Good as a non-HIV smoke test and for unit-test fixtures shaped like real data. | Far too few events for a penalty estimate. |
| **Burkina Faso malaria, natural immunity** | Dryad `10.5061/dryad.n1m33qq` / Zenodo 5002068, CC0 **[v]**. | 51 participants. | None in practice. | Not a vaccine; tiny. |
| **ENSEMBLE correlates paper** (Fong et al. Nat Microbiol 2022) | Springer source-data spreadsheets **[v]**. | Aggregate only: curves by region (Latin America, North America, South Africa). | Not usable as input. Useful later as a benchmark for an ENSEMBLE region E2. | Not participant-level. |

Searched with no usable hit: ImmPort (has HAI immunogenicity, e.g. FluPRINT,
but no infection endpoints found), CRAN beyond `vaccine`, Dryad and Zenodo for
COVID antibody-and-breakthrough cohorts (these release aggregates only, citing
ethics restrictions). The CoVPN analysis code is public
(`github.com/CoVPN/correlates_reporting2`), but its data is not.

## 4. Gated routes, in order of payoff per effort

### A. ENSEMBLE through YODA (Johnson & Johnson) — try first

Listing: yoda.yale.edu, NCT04505722 / VAC31518COV3001 **[v]**. ENSEMBLE 2
(VAC31518COV3009, two doses, 31,831 enrolled) is also listed **[v]**; its
published anonymisation report shows an `IS` (immunogenicity) domain with 7,040
rows, which looks like serostatus rather than the correlates panel.

Steps (from the YODA procedures, January 2025 **[v]**):

1. **Pre-submission inquiry.** Ask whether the COV3001 package includes the
   Day 29 spike IgG, RBD IgG and pseudovirus ID50 values for the
   immunogenicity subcohort and cases, the subcohort flag, and the
   case-cohort sampling strata. This decides whether the route is worth it.
2. **Register** on the YODA site and add COV3001 (and COV3009 if useful) to a
   request.
3. **Submit the request**: PI name, degrees and affiliation; key personnel;
   funding and conflict-of-interest forms (modified ICMJE); the research
   proposal (section 5); data specification (trials, IPD); certification of no
   litigation or commercial use.
4. **YODA review** for scientific merit and completeness. Revisions are common.
5. **Data Partner due diligence.** J&J checks it can share the data and that
   the variables exist.
6. **External review** if YODA asks for it (two independent peer reviewers).
7. **DUA** signed by the PI *and* the institution.
8. **Access** in YODA's secure platform. Data cannot be downloaded; the
   analysis runs there. We need to confirm which languages it offers and
   whether we can install this package. Proposal and decision are posted
   publicly. Results must go through peer review or a scientific meeting, and
   any change of aims needs an amendment.

Timeline: YODA does not publish one. Plan for a few months, and confirm in the
inquiry.

### B. HVTN / CoVPN auxiliary-study proposal — for harmonised COVE, ENSEMBLE, AZD1222, PREVENT-19 correlates data

This is the route for the USG/CoVPN correlates datasets themselves, which hold
the harmonised markers, weights and strata, and is the only route to a
three-platform E1. The sponsors own the clinical data for these trials, so
expect CoVPN and sponsor sign-off. Confirm that with SDMC before writing.

Steps (from hvtn.org Auxiliary Studies **[v]**):

1. **Office hours first.** The HVTN Statistics and Data Management Center
   holds drop-in hours (Fridays 9–10 AM PT and Tuesdays 10–11 AM PT, alternate
   weeks) to identify relevant data and walk through access. The contact named
   there is Holly Janes. A sign-up form lets you name the protocol in advance.
   Ask whether CoVPN COVID trials go through this process or a separate CoVPN
   one, and whether the sponsors must also approve.
2. **Fill in the proposal template** downloaded from that page.
3. **Email it** to `vtn.research@hvtn.org`.
4. **Auxiliary Studies Committee review.** It meets on the 1st Tuesday and 3rd
   Thursday. Submit at least 6 business days before a meeting. It checks
   feasibility and classifies the study (exploratory, ancillary or sub-study)
   against the parent consent language.
5. **Approvals and data agreement.** HVTN says 2 weeks to 6 months, longer for
   requesters outside the US, and does not expedite. State any grant deadline
   in the proposal.

The same route covers HVTN 705 (Imbokodo) and HVTN 706 (Mosaico) correlates if
another HIV E1 is wanted.

### C. Vivli — for AZD1222, CYD14/CYD15 dengue, RTS,S malaria

AstraZeneca, Sanofi and GSK share through Vivli (GSK moved there from
ClinicalStudyDataRequest.com in 2023 **[v]**). Search did not confirm listings
for D8110C00001, CYD14 (NCT01373281), CYD15 (NCT01374516) or RTS,S phase 3
(NCT00866619); use Vivli's Quick Study Look-up, or file an enquiry for unlisted
trials.

Steps (Vivli request guide and FAQ):

1. Create an account, find the studies, and fill in the data request form.
2. Vivli checks the form, then the data contributors and any independent
   review panel review it. Each can approve, reject or ask for revisions.
3. Sign Vivli's DUA, which is non-negotiable; the approved research plan is
   attached to it.
4. Data is anonymised and loaded into Vivli's secure research environment.

Timeline: Vivli says a few months on average, faster for one sponsor and a
tight proposal. As with YODA, ask before writing whether the correlates assay
data (PRNT50 for CYD, anti-CSP for RTS,S, the immunogenicity subset for
AZD1222) is in the sponsor package.

## 5. What every proposal must contain

Write one master proposal and adapt it per route.

- **Title and lay summary.** Measuring how much a correlate-of-protection
  model degrades when moved across populations or trials. Stress that it makes
  no efficacy claim and does no pathogen characterisation.
- **Background and significance.** Borrowed correlates are used in
  accelerated approval, and nobody has measured the loss of transfer
  (literature-review.md sections 3–5).
- **Specific aims.** E2: fit on one region, apply to the others. E1, where
  allowed: fit on one trial, apply to another.
- **Exact variable list.** Markers by visit and assay; endpoint and time;
  subcohort flag, sampling strata and weights; split variables; per-protocol
  and baseline-negative flags; placebo-arm events per stratum.
- **Methods.** Ridge logistic and ridge Cox, tuned inside the source only;
  inverse-probability weights; AUC or C, calibration slope and intercept,
  Brier score, and the decision-flip rate; bootstrap intervals. Link this repo,
  since the code is public and tested.
- **Sample-size note.** Events per stratum and why they are enough.
- **Data security.** Analysis stays inside the provider's platform; no
  re-identification; no redistribution; outputs are aggregate only.
- **Ethics.** Most providers expect an IRB exemption or a non-human-subjects
  determination for secondary de-identified analysis. Get it from the
  requester's institution.
- **Team, affiliation, funding, conflicts.** PI and institution signatory.
- **Timeline and dissemination.** Preprint plus peer-reviewed paper; credit
  lines as the provider requires.

## 6. Recommended next steps

1. **Now, free:** download AMP (`GDESJI`, `50VH2C`) and HVTN 702 (`JI2FUM`)
   from the HVTN Dataverse. The guestbook asks for the downloader's own name
   and institution, so Imo has to do this step. Then add an AMP loader to the
   E2 runner, split by trial. Keep the files outside git.
2. **This week:** send YODA the pre-submission inquiry about the COV3001
   correlates variables (section 4A, step 1).
3. **This week:** book an HVTN SDMC office-hours slot about the CoVPN COVID
   correlates datasets (section 4B, step 1).
4. **In parallel:** settle the requesting institution and its IRB route, since
   YODA, Vivli and HVTN all need an institutional signatory.
5. **Optional:** register on CAVD DataSpace for RV144, to run an RV144 to HVTN
   702 E1 whose answer (the correlate did not carry over) is known in advance.

## Sources

- HVTN Auxiliary Studies and public data index: https://www.hvtn.org/scientific-programs/scientific-programs-overview/auxiliary-studies.html
- HVTN Harvard Dataverse: https://dataverse.harvard.edu/dataverse/HVTN
- Atlas HVTN Public Data (pre-April-2025 releases): https://atlas.scharp.org/cpas/project/HVTN%20Public%20Data/begin.view
- YODA ENSEMBLE listing: https://yoda.yale.edu/clinical-trial/nct04505722-a-randomized-double-blind-placebo-controlled-phase-3-study-to-assess-the-efficacy-and-safety-of-ad26-cov2-s-for-the-prevention-of-sars-cov-2-mediated-covid-19-in-adults-aged-18-years-a/
- YODA procedures (Jan 2025): https://yoda.yale.edu/wp-content/uploads/2022/11/YODA-Project-Data-Release-Procedures-January-2025.pdf
- Vivli request process: https://vivli.org/about/data-request-review-process/ and FAQ on timelines
- GSK data sharing move to Vivli: https://www.gsk-studyregister.com/About_GSK_Patient_Level_Data_Sharing_Final_13July2023.pdf
- Corey et al. NEJM 2021 (AMP): https://pmc.ncbi.nlm.nih.gov/articles/PMC8189692
- Moodie et al. JID 2022 (HVTN 702 correlates): https://pmc.ncbi.nlm.nih.gov/articles/PMC9890908
- Fong et al. Nat Microbiol 2022 (ENSEMBLE correlates): https://www.nature.com/articles/s41564-022-01262-1
- Zenodo 4972891 (Hong Kong influenza B): https://zenodo.org/records/4972891
- CAVD DataSpace: https://dataspace.cavd.org
