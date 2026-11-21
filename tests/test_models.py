import numpy as np
import pytest

from penalty import models, schema, transfer
from sim.generate import BETA, Gap, case_cohort, make_studies


@pytest.fixture(scope="module")
def cohort():
    s, _ = make_studies(20_000, 10, Gap(), np.random.default_rng(10))
    return schema.validate(s)


def test_ridge_recovers_true_coefficients_on_full_cohort(cohort):
    m = models.fit_ridge_logistic(cohort)
    beta = m.coef / m.scale
    assert np.max(np.abs(beta - BETA)) < 0.1


def test_weighted_fit_on_case_cohort_matches_full_cohort(cohort):
    full = models.fit_ridge_logistic(cohort)
    cc = schema.validate(case_cohort(cohort, 0.1, np.random.default_rng(11)))
    w = models.fit_ridge_logistic(cc)
    unw = models.fit_ridge_logistic(cc.assign(sampling_weight=1.0))
    # Risk on the full cohort: weighted model reproduces it, unweighted one
    # is shifted by roughly log(1/f) on the logit scale.
    lp_full, lp_w, lp_u = (m.linear_predictor(cohort) for m in (full, w, unw))
    assert abs(np.mean(lp_w - lp_full)) < 0.1
    assert np.mean(lp_u - lp_full) > 1.5


def test_single_marker_picks_the_strongest(cohort):
    m = models.fit_best_single_marker(cohort)
    assert m.markers == ("marker_0",)


def test_tuning_cannot_see_target():
    s, t = make_studies(800, 400, Gap(), np.random.default_rng(12))
    r = transfer.run(s, t)
    assert not (r.model.fitted_ids & set(t["participant_id"]))
    # Feeding training rows back as a "target" is refused.
    s = schema.validate(s)
    train, held = transfer.split_source(s, 0.3, np.random.default_rng(0))
    m = models.fit_ridge_logistic(train)
    with pytest.raises(RuntimeError, match="used to fit"):
        transfer._evaluate(m, train)
