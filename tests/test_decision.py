"""Decision-flip metric (build step 6)."""

import numpy as np
import pandas as pd
import pytest
from scipy.special import expit

from penalty import metrics, transfer
from sim.generate import Gap, make_studies, placebo_incidence


def test_youden_cut_matches_brute_force():
    rng = np.random.default_rng(0)
    y = rng.integers(0, 2, 300)
    p = np.round(rng.random(300), 2)
    w = rng.uniform(0.5, 4, 300)

    def j(c):
        pos = p >= c
        return (w[pos & (y == 1)].sum() / w[y == 1].sum()) - (w[pos & (y == 0)].sum() / w[y == 0].sum())

    best = max(np.unique(p), key=j)
    assert j(metrics.youden_cut(y, p, w)) == pytest.approx(j(best))


def test_calibrated_model_flips_nothing_at_a_fixed_threshold():
    rng = np.random.default_rng(1)
    p = expit(rng.normal(-1, 1.2, 200_000))
    y = rng.binomial(1, p)
    out = metrics.decision_flip(y, p, np.ones_like(p), risk_threshold=0.3, source_cut=0.3)
    assert out["flip_fixed"] < 0.01


def runs(gap, reps=10, seed0=0, placebo=False):
    rows = []
    for s in range(reps):
        rng = np.random.default_rng(seed0 + s)
        src, tgt = make_studies(4000, 3000, gap, rng)
        pl = None
        if placebo:
            n = 20_000
            pl = transfer.PlaceboArms(round(placebo_incidence(n, Gap(), rng) * n), n,
                                      round(placebo_incidence(n, gap, rng) * n), n)
        r = transfer.run(src, tgt, seed=s, placebo=pl)
        rows.append(r.all_penalties())
    return pd.DataFrame(rows).mean()


def test_null_gap_flips_few_decisions():
    pen = runs(Gap(), seed0=100)
    assert abs(pen["flip_fixed"]) < 0.03 and abs(pen["flip_optimal"]) < 0.03


def test_population_shift_flips_decisions_that_auc_cannot_see():
    """The decision-level version of AC3: a risk-scale shift leaves AUC alone
    and moves a large share of fixed-threshold decisions."""
    pen = runs(Gap(population=1.0), seed0=200)
    assert abs(pen["auc"]) < 0.02
    assert pen["flip_fixed"] > 0.15


def test_exposure_offset_removes_attack_rate_flips():
    pen = runs(Gap(attack_rate=1.0), seed0=300, placebo=True)
    assert pen["flip_fixed"] > 0.15
    assert abs(pen["flip_fixed_exposure_adjusted"]) < 0.03


def test_cv_baseline_reports_decision_metrics():
    src, tgt = make_studies(600, 600, Gap(population=1.0), np.random.default_rng(5))
    r = transfer.run(src, tgt, baseline="cv")
    assert {"flip_fixed", "flip_optimal"} <= set(r.penalty)
