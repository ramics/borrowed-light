import numpy as np
import pandas as pd
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


@pytest.mark.parametrize("C", [0.001, 0.1, 10.0])
def test_ridge_solver_matches_sklearn_objective(C):
    from sklearn.linear_model import LogisticRegression
    rng = np.random.default_rng(5)
    x = rng.normal(size=(1500, 6))
    y = rng.binomial(1, 0.3, 1500)
    w = rng.uniform(1, 5, 1500)
    ref = LogisticRegression(C=C, max_iter=10_000, tol=1e-10).fit(x, y, sample_weight=w / w.mean())
    b0, b = models._fit_ridge(x, y, w, C)
    np.testing.assert_allclose(b, ref.coef_[0], atol=1e-6)
    assert b0 == pytest.approx(ref.intercept_[0], abs=1e-6)


def test_bootstrap_copies_of_one_person_stay_on_one_side_of_every_split():
    from penalty.resample import _resample
    s, _ = make_studies(400, 10, Gap(), np.random.default_rng(13))
    boot = _resample(schema.validate(s), np.random.default_rng(14))
    base = lambda d: set(d["participant_id"].str.split("#").str[0])  # noqa: E731
    train, held = transfer.split_source(boot, 0.3, np.random.default_rng(15))
    assert not base(train) & base(held)
    fold = transfer.assign_folds(boot, 5, np.random.default_rng(16))
    people = boot["participant_id"].str.split("#").str[0]
    assert (pd.Series(fold).groupby(people.to_numpy()).nunique() == 1).all()


def test_phase_one_bootstrap_lets_the_case_count_vary():
    from penalty.resample import _resample
    s, _ = make_studies(3000, 10, Gap(), np.random.default_rng(17))
    cc = schema.validate(case_cohort(schema.validate(s), 0.1, np.random.default_rng(18), keep_phase_one=True))
    counts = {int(_resample(cc, np.random.default_rng(i))["endpoint"].sum()) for i in range(10)}
    assert len(counts) > 5
