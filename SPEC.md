# Transfer-penalty estimation for borrowed correlates of protection

**Status:** spec for a first build. A spike validating the core estimator is in
`spike_transfer_penalty.py` and all four acceptance checks pass.

## Why this exists

A countermeasure against a fast-moving or engineered pathogen cannot be
efficacy-tested: there is no time, no stable case pool, and no ethical placebo
arm during a lethal outbreak. The alternative is to argue from immune markers —
"this response looks like the one that protected people in a trial of something
similar." That argument is routine in principle and unquantified in practice.
Nobody can say how much weaker it gets as you stretch it across a different
pathogen, a different vaccine platform, a different population, or a cohort of
thirty instead of thirty thousand.

This project puts numbers on that. The unit of analysis is a *study*, not a
blood sample, so no wet lab is required: the measurements have already been
made by other people.

## Scope

**In scope — the first two experiments.**

- **E1, cross-study transfer.** Fit a correlate model on study A, apply it to
  study B which differs by platform and/or population. Report what was lost.
- **E2, within-study population transfer.** Fit on one stratum of a single
  study, apply to a held-out stratum (age band, prior-exposure status, region,
  site). Same machinery, smaller gap, no cross-study data agreement needed.

**Explicitly not in scope.** New markers. New assays. Any claim about a specific
pathogen's biology. Any prediction of efficacy — the output is a statement about
how much an inference degrades, never a statement that a vaccine works. No
pathogen characterisation of any kind.

## The one operation

E1 and E2 are the same thing with different split axes. Given a source stratum
`S` and target stratum `T`:

```
fit        m = model(S_train)
baseline   b = evaluate(m, S_heldout)     # what the model is worth at home
transfer   t = evaluate(m, T)             # what it is worth across the gap
penalty    = b - t                        # per metric, with uncertainty
```

Everything else in the repo is plumbing around that.

## Metrics, and the one that actually matters

Report all four. They fail in different ways and a regulator needs to see which.

| Metric | What it tells you |
| --- | --- |
| AUC (or Harrell's C for time-to-event) | Does the model still rank people correctly? |
| Calibration slope | Is the model over- or under-confident in the target? |
| Calibration intercept | Is the whole risk scale shifted? |
| Brier score | Overall accuracy, combining both |

**Calibration is the headline, not discrimination.** The spike demonstrates why.
Under a pure population shift, AUC is essentially untouched while the
calibration intercept moves one-for-one with the shift:

| Baseline risk shift | AUC penalty | Calibration intercept |
| --- | --- | --- |
| 0.0 | −0.006 | +0.00 |
| 0.5 | −0.003 | +0.49 |
| 1.0 | −0.002 | +1.00 |
| 1.5 | −0.001 | +1.50 |

A model that ranks perfectly and sits on the wrong scale will put the decision
threshold in the wrong place. Since admissibility is a threshold decision, a
report that leads with AUC hides exactly the failure that matters. This is a
publishable point on its own and should survive into the write-up.

Add a fifth, decision-level metric once the first four work: the fraction of
candidates for which the source-derived threshold gives a different
accept/reject answer than a target-derived threshold would.

## Models

Match the parent programme so results are directly comparable to Repop-TB
model 16:

- Binary endpoint → ridge-penalised logistic regression.
- Time-to-event → ridge-penalised Cox.
- Always report the best single marker as the comparator. If the multivariate
  model cannot beat one well-chosen marker, say so.

Penalty strength chosen by cross-validation **within the source only**. Touching
the target to tune anything invalidates the whole measurement; make this
structurally impossible rather than a rule people remember.

## Data interface

One tidy schema, so synthetic and real data are interchangeable and the real
data never needs the code to change.

| Column | Type | Notes |
| --- | --- | --- |
| `participant_id` | str | |
| `study_id` | str | split axis for E1 |
| `stratum` | str | split axis for E2 |
| `platform` | str | vaccine platform, if known |
| `marker_*` | float | wide numeric columns, one per assay readout |
| `endpoint` | int | 0/1 |
| `time_to_event` | float | optional, for Cox |
| `sampling_weight` | float | **required**, default 1.0 |

`sampling_weight` is not optional plumbing. Correlates datasets are frequently
released as case-cohort or case-control subsamples; ignoring the weights biases
every estimate and the bias is not small. Weighted fitting and weighted metrics
must be in from the first commit, with a test that a weighted case-cohort
subsample reproduces the full-cohort estimate within Monte Carlo error.

Marker alignment across studies is the other quiet killer. Require an explicit
mapping file declaring which marker in study A corresponds to which in study B
and on what scale; refuse to run E1 without one.

## Uncertainty

- Participant-level stratified bootstrap, default 1000 replicates.
- Where three or more studies exist, also report leave-one-study-out.
- Report intervals everywhere. A point estimate of a transfer penalty is
  useless and slightly dishonest.

## Simulation harness

The estimator is validated against data where the true gap is set by a dial, so
the right answer is known. Four independent dials, each zero when source and
target are identical:

- `mechanism` — rotates the true coefficient vector. "The markers mean something
  different over there."
- `assay` — scale/location shift on a subset of markers. "Measured differently."
- `population` — baseline risk shift. "Different people."
- `noise` — extra measurement error in the target.

Dials must be independently settable and composable. Most real gaps are
combinations, and the interesting question is whether the decomposition can tell
them apart.

## Acceptance criteria

Numbers in brackets are what the spike achieved; treat them as the floor.

1. **Null check.** With all dials at zero, the 95% interval for the penalty
   covers zero. [mean −0.0002, CI −0.037 to +0.038]
2. **Monotonicity.** Estimated penalty correlates with the mechanism dial at
   ρ ≥ 0.9 across the range. [0.935]
3. **Calibration sensitivity.** Under population shift alone, AUC penalty moves
   less than 0.02 while the calibration intercept tracks the shift. [confirmed,
   intercept recovers the shift to two decimal places]
4. **Small-n humility.** With a target cohort of 30, the interval is at least
   three times the width at n=1500. [0.279 vs 0.070, a factor of 4]
5. **Sampling weights.** A weighted case-cohort subsample reproduces the
   full-cohort penalty within Monte Carlo error. *(not yet built)*
6. **Runtime.** Full synthetic suite runs in under five minutes on a laptop.
   [13 seconds at spike scale]
7. **No real participant data in the repository, ever.** Synthetic fixtures
   only. Real data is analysed where it lives, under its own data-use terms.

## Repo layout

```
penalty/
  schema.py        tidy-frame validation, weight handling, marker mapping
  models.py        ridge logistic, ridge Cox, single-marker comparator
  metrics.py       AUC/C, calibration slope + intercept, Brier, decision flip rate
  transfer.py      the fit/baseline/transfer/penalty operation
  resample.py      stratified bootstrap, leave-one-study-out
sim/
  generate.py      four-dial simulator
experiments/
  e1_cross_study.py
  e2_within_study.py
  configs/*.yaml
reports/           generated tables and figures, one directory per run
tests/
```

Every experiment writes a run manifest: config, code version, seed, package
versions. If a result cannot be regenerated from its manifest alone, it does not
count.

## Build order

1. Schema, sampling weights end to end (weighted fitting *and* weighted
   metrics), core metrics, tests including AC5. Nothing interesting, everything
   depends on it. Do not defer the weights — retrofitting them means redoing
   every number computed before they landed.
2. Port the spike into `transfer.py` and `sim/generate.py`, with AC1–AC4 as
   automated tests rather than printed output.
3. Add the ridge Cox path for time-to-event endpoints.
4. E2 on a real single-study dataset — the first real-data milestone, and the
   one that needs the least paperwork.
5. E1 across two studies, once data agreements land.
6. Decision-flip metric and the report template.

## Risks

- **Data access is the binding constraint, not compute.** Participant-level data
  needs data-use agreements; summary tables are not enough. Line up two or three
  before relying on any one. The M72/AS01E correlates terms need confirming with
  VISC and Fred Hutch, and that dataset becomes public on publication.
- **Marker alignment between studies** may turn out to be impossible for some
  pairs. That is itself a finding worth reporting.
- **Too few studies** for leave-one-study-out. Fall back to pairwise reporting
  and say so.
- **The penalties may not generalise.** What transfers from this work is the
  estimator, not the numbers it produces on any one pair of studies. The write-up
  must say this plainly; overclaiming here is the fastest way to lose a reviewer.

## Open questions

- Which single-study dataset is easiest to get for E2? That choice sets the
  timeline for the first real result.
- Should the decision-flip metric use a fixed threshold or a target-optimal one?
  Probably both, reported side by side.
- Is a Bayesian version worth it for the small-n case, where the bootstrap is
  known to be unreliable? Likely yes, as a second pass.
