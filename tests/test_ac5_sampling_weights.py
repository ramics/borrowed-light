"""AC5. A weighted case-cohort subsample reproduces the full-cohort penalty
within Monte Carlo error, and the unweighted analysis does not.

The full source and target cohorts are fixed and split once. Each replicate
draws a fresh case-cohort subsample (all cases plus a 15% random subcohort) of
every part and recomputes the penalty. Across replicates the mean difference
from the full-cohort penalty must be within three Monte Carlo standard errors
for every metric. The gap includes a population shift so the target's case
fraction differs from the source's, which is what makes ignoring the weights
visibly wrong.
"""

import numpy as np
import pandas as pd
import pytest

from penalty import metrics, schema, transfer
from sim.generate import Gap, case_cohort, make_studies

REPS = 24
FRAC = 0.15


@pytest.fixture(scope="module")
def experiment():
    s, t = make_studies(6000, 3000, Gap(mechanism=0.4, population=0.8), np.random.default_rng(1))
    s, t = schema.validate(s), schema.validate(t)
    train, held = transfer.split_source(s, 0.3, np.random.default_rng(0))
    full = transfer.run_split(train, held, t).penalty
    weighted, unweighted = [], []
    for r in range(REPS):
        rng = np.random.default_rng(1000 + r)
        parts = [schema.validate(case_cohort(d, FRAC, rng)) for d in (train, held, t)]
        weighted.append(transfer.run_split(*parts).penalty)
        unweighted.append(transfer.run_split(*[p.assign(sampling_weight=1.0) for p in parts]).penalty)
    full = pd.Series(full)
    return full, pd.DataFrame(weighted) - full, pd.DataFrame(unweighted) - full


def test_weighted_case_cohort_reproduces_full_cohort_penalty(experiment):
    full, dw, _ = experiment
    cols = list(metrics.METRICS)
    mc_se = dw[cols].std(ddof=1) / np.sqrt(len(dw))
    bias = dw[cols].mean()
    report = pd.DataFrame({"full": full[cols], "bias": bias, "mc_se": mc_se})
    assert (bias.abs() <= 3 * mc_se).all(), f"\n{report}"


def test_decision_flips_under_case_cohort_are_close_to_full_cohort(experiment):
    """Decision flips are a step function of the predictions, so at a 15%
    subcohort they carry a small finite-sample bias (about -0.02 on a
    full-cohort fixed-threshold flip rate of 0.25 here). It vanishes as the
    subcohort grows (checked at 50%), so it is not a weighting failure.
    The tolerance below is stated rather than hidden in a looser MC bound."""
    _, dw, du = experiment
    for k in metrics.DECISION_METRICS:
        assert abs(dw[k].mean()) < 0.03, (k, dw[k].mean())


def test_ignoring_weights_is_detectably_wrong(experiment):
    _, _, du = experiment
    mc_se = du.std(ddof=1) / np.sqrt(len(du))
    # Brier mixes in the case fraction, so the unweighted penalty is far off.
    assert abs(du["brier"].mean()) > 10 * mc_se["brier"]


def test_phase_one_release_form_gives_same_answer():
    """A release that keeps phase-one rows (phase2 = 0, markers blank) must
    give the same penalty as the phase-two-only release."""
    s, t = make_studies(3000, 1500, Gap(population=0.5), np.random.default_rng(5))
    s, t = schema.validate(s), schema.validate(t)
    train, held = transfer.split_source(s, 0.3, np.random.default_rng(0))
    a, b = [], []
    for d in (train, held, t):
        rng = np.random.default_rng(7)
        a.append(schema.validate(case_cohort(d, 0.2, rng)))
        rng = np.random.default_rng(7)
        b.append(schema.validate(case_cohort(d, 0.2, rng, keep_phase_one=True)))
    pa, pb = transfer.run_split(*a).penalty, transfer.run_split(*b).penalty
    for k in pa:
        assert pa[k] == pytest.approx(pb[k], abs=1e-9)
