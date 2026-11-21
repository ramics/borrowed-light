import itertools

import numpy as np
import pytest
from sklearn.metrics import brier_score_loss, roc_auc_score
from scipy.special import expit, logit

from penalty import metrics


@pytest.fixture
def data():
    rng = np.random.default_rng(0)
    n = 400
    p = expit(rng.normal(-1, 1.2, n))
    y = rng.binomial(1, p)
    w = rng.uniform(0.5, 6, n)
    return y, p, w


def test_auc_matches_brute_force_pairs_with_ties():
    rng = np.random.default_rng(1)
    y = rng.integers(0, 2, 60)
    s = rng.integers(0, 8, 60).astype(float)  # many ties
    w = rng.uniform(0.2, 4, 60)
    num = den = 0.0
    for i, j in itertools.product(np.where(y == 1)[0], np.where(y == 0)[0]):
        num += w[i] * w[j] * ((s[i] > s[j]) + 0.5 * (s[i] == s[j]))
        den += w[i] * w[j]
    assert metrics.weighted_auc(y, s, w) == pytest.approx(num / den, abs=1e-12)


def test_auc_and_brier_agree_with_sklearn(data):
    y, p, w = data
    assert metrics.weighted_auc(y, p, w) == pytest.approx(roc_auc_score(y, p, sample_weight=w))
    assert metrics.weighted_auc(y, p) == pytest.approx(roc_auc_score(y, p))
    assert metrics.weighted_brier(y, p, w) == pytest.approx(brier_score_loss(y, p, sample_weight=w))


def test_integer_weights_equal_row_replication(data):
    y, p, _ = data
    k = np.random.default_rng(2).integers(1, 4, len(y))
    rep = metrics.evaluate(np.repeat(y, k), np.repeat(p, k))
    wtd = metrics.evaluate(y, p, k.astype(float))
    for m in metrics.METRICS:
        assert wtd[m] == pytest.approx(rep[m], rel=1e-6, abs=1e-8)


def test_calibration_recovers_known_miscalibration():
    rng = np.random.default_rng(3)
    lp = rng.normal(-0.5, 1.5, 200_000)
    y = rng.binomial(1, expit(0.7 + 0.6 * lp))
    slope, joint = metrics.calibration_slope(y, expit(lp))
    assert slope == pytest.approx(0.6, abs=0.02)
    assert joint == pytest.approx(0.7, abs=0.02)
    # Shift-only miscalibration: calibration-in-the-large finds the shift.
    y2 = rng.binomial(1, expit(lp + 0.9))
    assert metrics.calibration_in_the_large(y2, expit(lp)) == pytest.approx(0.9, abs=0.02)


def test_perfect_calibration(data):
    rng = np.random.default_rng(4)
    p = expit(rng.normal(-1, 1, 100_000))
    y = rng.binomial(1, p)
    r = metrics.evaluate(y, p)
    assert r["cal_slope"] == pytest.approx(1, abs=0.03)
    assert r["cal_intercept"] == pytest.approx(0, abs=0.03)


def test_single_class_gives_nan_not_crash():
    r = metrics.evaluate(np.zeros(10), np.full(10, 0.2))
    assert np.isnan(r["auc"]) and np.isnan(r["cal_slope"]) and np.isnan(r["cal_intercept"])


def test_penalty_signs():
    base = {"auc": 0.8, "brier": 0.10, "cal_intercept": 0.0, "cal_slope": 1.0, "cal_intercept_joint": 0.0}
    tgt = {"auc": 0.7, "brier": 0.15, "cal_intercept": 0.5, "cal_slope": 0.8, "cal_intercept_joint": 0.4}
    pen = metrics.penalty(base, tgt)
    assert pen["auc"] == pytest.approx(0.1) and pen["brier"] == pytest.approx(0.05)
    assert pen["cal_intercept"] == pytest.approx(0.5) and pen["cal_slope"] == pytest.approx(-0.2)
