"""Ridge Cox path (build step 3): weighted fit, weighted metrics, and the
same acceptance properties as the binary path where they apply."""

import numpy as np
import pandas as pd
import pytest

from penalty import schema, survival, transfer
from sim.generate import BETA, Gap, case_cohort, make_studies, placebo_incidence


def studies(gap, seed, n_source=4000, n_target=2000):
    return make_studies(n_source, n_target, gap, np.random.default_rng(seed), survival=True)


def test_simulator_carries_time_to_event():
    s, _ = studies(Gap(), 0)
    s = schema.validate(s)
    assert (s["time_to_event"] <= 1.0).all() and 0.15 < s["endpoint"].mean() < 0.4


def test_unpenalised_cox_recovers_true_log_hazard_ratios():
    s, _ = make_studies(20_000, 10, Gap(), np.random.default_rng(1), survival=True)
    x, t, d, w = survival._design(schema.validate(s), [f"marker_{j}" for j in range(8)])
    beta = survival.fit_cox(x, t, d, w)
    assert np.max(np.abs(beta - BETA)) < 0.06


def test_weighted_c_index_equals_row_replication():
    rng = np.random.default_rng(2)
    n = 150
    t = rng.exponential(1, n).round(1)   # ties in time
    d = rng.integers(0, 2, n)
    s = rng.normal(size=n).round(1)      # ties in score
    k = rng.integers(1, 4, n)
    a = survival.weighted_c_index(t, d, s, k.astype(float))
    b = survival.weighted_c_index(np.repeat(t, k), np.repeat(d, k), np.repeat(s, k), np.ones(k.sum()))
    assert a == pytest.approx(b, abs=1e-12)


def test_weighted_cox_on_case_cohort_matches_full_cohort():
    s, _ = make_studies(12_000, 10, Gap(), np.random.default_rng(3), survival=True)
    s = schema.validate(s)
    m = [f"marker_{j}" for j in range(8)]
    full = survival.fit_cox(*survival._design(s, m))
    cc = schema.validate(case_cohort(s, 0.1, np.random.default_rng(4)))
    wtd = survival.fit_cox(*survival._design(cc, m))
    assert np.max(np.abs(wtd - full)) < 0.12


def test_cox_null_penalty_is_small():
    rows = [transfer.run(*studies(Gap(), 10 + i), model="ridge_cox", seed=i).penalty for i in range(8)]
    pen = pd.DataFrame(rows).mean()
    assert abs(pen["c_index"]) < 0.01 and abs(pen["cal_intercept"]) < 0.06


def test_cox_population_shift_moves_calibration_not_discrimination():
    """AC3 on the time-to-event path. A log-hazard shift leaves C almost
    unchanged; log(O/E) recovers the shift one for one."""
    res = {}
    for d in (0.0, 1.0):
        rows = [transfer.run(*studies(Gap(population=d), 30 + i), model="ridge_cox", seed=i).penalty
                for i in range(6)]
        res[d] = pd.DataFrame(rows).mean()
    assert abs(res[1.0]["c_index"] - res[0.0]["c_index"]) < 0.02
    assert res[1.0]["cal_intercept"] == pytest.approx(1.0, abs=0.1)


def test_cox_mechanism_gap_costs_discrimination():
    pen = [transfer.run(*studies(Gap(mechanism=d), 50), model="ridge_cox").penalty["c_index"]
           for d in (0.0, 0.5, 1.0)]
    assert pen[0] < pen[1] < pen[2]


def test_cox_weighted_case_cohort_penalty_matches_full():
    """AC5 on the Cox path, a lighter version of the binary-path test."""
    s, t = studies(Gap(population=0.6), 60, 6000, 3000)
    s, t = schema.validate(s), schema.validate(t)
    train, held = transfer.split_source(s, 0.3, np.random.default_rng(0))
    kw = dict(model="ridge_cox", model_kw={"c_grid": (1.0,)}, horizon=0.5)
    full = pd.Series(transfer.run_split(train, held, t, **kw).penalty)
    diffs = []
    for r in range(8):
        rng = np.random.default_rng(700 + r)
        parts = [schema.validate(case_cohort(d, 0.15, rng)) for d in (train, held, t)]
        diffs.append(pd.Series(transfer.run_split(*parts, **kw).penalty) - full)
    dw = pd.DataFrame(diffs)
    se = dw.std(ddof=1) / np.sqrt(len(dw))
    assert (dw.mean().abs() <= 3 * se + 0.01).all(), pd.DataFrame({"bias": dw.mean(), "se": se})


def test_cox_placebo_offset_absorbs_attack_rate():
    gap = Gap(attack_rate=0.8)
    rng = np.random.default_rng(80)
    s, t = make_studies(4000, 3000, gap, rng, survival=True)
    n = 50_000
    pl = transfer.PlaceboArms(round(placebo_incidence(n, Gap(), rng, True) * n), n,
                              round(placebo_incidence(n, gap, rng, True) * n), n)
    r = transfer.run(s, t, model="ridge_cox", placebo=pl)
    assert r.penalty["cal_intercept"] == pytest.approx(0.8, abs=0.12)
    assert r.penalty_adjusted["cal_intercept"] == pytest.approx(0.0, abs=0.12)


def test_cox_needs_time_to_event():
    s, t = make_studies(500, 200, Gap(), np.random.default_rng(9))
    with pytest.raises(schema.SchemaError, match="time_to_event"):
        transfer.run(s, t, model="ridge_cox")
