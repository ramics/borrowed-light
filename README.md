# borrowed-light

**How much does a correlate of protection lose when you borrow it across a gap?**

A countermeasure against a fast-moving or engineered pathogen cannot be efficacy-tested: there is no time, no stable case pool, and no ethical placebo arm during a lethal outbreak. The fallback is to argue from immune markers — *this response looks like the one that protected people in a trial of something similar*. That argument is routine in principle and unquantified in practice. Nobody can say how much weaker it gets across a different platform, a different population, or a cohort of thirty instead of thirty thousand.

borrowed-light puts an interval on it. The unit of analysis is a *study* rather than a blood sample, so the measurements have already been made by other people and no wet lab is required.

## The one operation

Given a source stratum `S` and a target stratum `T`:

```
fit        m = model(S_train)
baseline   b = evaluate(m, S_heldout)     # what the model is worth at home
transfer   t = evaluate(m, T)             # what it is worth across the gap
penalty    = b - t                        # per metric, with uncertainty
```

Penalty strength is cross-validated **within the source only**; the target is never touched to tune anything, and that is enforced by the code rather than left to discipline. Everything else in the repo is plumbing around those four lines.

## Calibration is the headline, and that is the point

Four metrics ship with every run: AUC (or Harrell's C), calibration slope, calibration intercept, Brier. They fail in different ways and a regulator needs to see which. Under a pure population shift, discrimination barely moves while calibration tracks the shift one-for-one:

| Baseline risk shift | AUC penalty | Calibration intercept |
| --- | --- | --- |
| 0.0 | −0.006 | +0.00 |
| 0.5 | −0.003 | +0.49 |
| 1.0 | −0.002 | +1.00 |
| 1.5 | −0.001 | +1.50 |

A model that ranks people perfectly and sits on the wrong scale puts the decision threshold in the wrong place. Since admissibility *is* a threshold decision, a report that leads with AUC hides exactly the failure that matters. A fifth, decision-level metric counts it directly: the fraction of candidates whose accept/reject answer changes across the gap, at a fixed threshold and at a source-optimal one, reported side by side.

## What is validated

The estimator is checked against simulated data where the true gap is set by four independent dials — `mechanism` (the markers mean something different over there), `assay` (measured differently), `population` (different people), `noise` — each zero when source and target are identical. The acceptance criteria are automated tests rather than printed output:

| Check | Result |
| --- | --- |
| Null: all dials at zero, the 95% interval covers zero | mean −0.0002, CI −0.037 to +0.038 |
| Monotonicity against the mechanism dial | ρ = 0.935 |
| Calibration sensitivity under population shift alone | intercept recovers the shift to two decimal places |
| Small-n humility: n=30 against n=1500 | intervals 0.279 against 0.070, a factor of 4 |
| Sampling weights: a case-cohort subsample reproduces the full cohort | within Monte Carlo error |
| Runtime, full synthetic suite | 13 seconds |

Sampling weights are load-bearing rather than plumbing. Correlates datasets are routinely released as case-cohort or case-control subsamples, and ignoring the weights biases every estimate by a margin that is not small, so weighted fitting and weighted metrics have been in since the first commit.

## Real data

Experiment configs are committed for the gaps worth measuring first: the AMP bnAb trials (HVTN 703 and 704, split by trial — same product, two populations), HVTN 505 by age band, and the Hong Kong influenza B cohort by strain. Data access is the binding constraint rather than compute, and [`docs/data-access.md`](docs/data-access.md) is the live register of what is public, what needs an agreement, and what each route costs in events per side.

No participant data is committed here, ever. Real data is analysed where it lives, under its own terms; `.gitignore` blocks the common formats as a backstop.

## Layout

| Module | Role |
| --- | --- |
| `penalty/schema.py` | tidy-frame validation, sampling weights, two-phase rows, marker mapping |
| `penalty/metrics.py` | weighted AUC, Brier, calibration-in-the-large, calibration slope, decision flips |
| `penalty/models.py` | weighted ridge logistic and ridge Cox, best single marker as comparator |
| `penalty/transfer.py` | fit / baseline / transfer / penalty |
| `penalty/resample.py` | stratified bootstrap, leave-one-study-out |
| `sim/generate.py` | four-dial simulator and case-cohort sampler |
| `experiments/` | E2 runner, configs, run manifests |

Every experiment writes a manifest with its config, code version, seed and package versions. A result that cannot be regenerated from its manifest alone does not count.

```
pip install -e ".[test]"
pytest        # 73 tests
```

[`SPEC.md`](SPEC.md) is the design and the source of truth for scope. [`docs/literature-review.md`](docs/literature-review.md) places the work: the field already distinguishes a correlate valid inside one trial from one valid across settings, and already has formal tools that take "how imperfect is the transport?" as an *input*. Estimating that input from data is the gap this fills.

## Scope

The output is a statement about how much an inference degrades. It is never a claim that a vaccine works, and the project characterises no pathogen.
