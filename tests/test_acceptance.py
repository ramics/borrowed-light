"""AC1-AC4 from SPEC.md as automated tests.

Two kinds of interval appear here, and the difference matters:

* *simulation* intervals re-draw fresh source and target data in every
  replicate. They are the sampling distribution under a known truth, and are
  what the spike reported.
* *bootstrap* intervals resample one fixed dataset at participant level. They
  are what real data gets, so AC1 and AC4 are checked with both.

Numbers in SPEC.md are the floor; the assertions below use the spec's
thresholds.
"""

import numpy as np
import pandas as pd
import pytest

from penalty import resample, transfer
from sim.generate import Gap, make_studies


def simulate_runs(gap, reps, seed0, n_source=4000, n_target=1500):
    rows = []
    for s in range(reps):
        src, tgt = make_studies(n_source, n_target, gap, np.random.default_rng(seed0 + s))
        rows.append(_run(src, tgt, s))
    return pd.DataFrame(rows)


def _run(src, tgt, seed):
    r = transfer.run(src, tgt, seed=seed)
    return {**r.penalty, **{"target_" + k: v for k, v in r.target.items()},
            **{"baseline_" + k: v for k, v in r.baseline.items()}}


# AC1 ---------------------------------------------------------------------

def test_ac1_null_simulation_interval_covers_zero():
    runs = simulate_runs(Gap(), reps=150, seed0=0)
    for k in ("auc", "brier", "cal_intercept", "cal_slope"):
        lo, hi = runs[k].quantile([0.025, 0.975])
        assert lo <= 0 <= hi, f"{k}: [{lo:.4f}, {hi:.4f}]"
    assert abs(runs["auc"].mean()) < 0.005


def test_ac1_null_bootstrap_intervals_cover_zero():
    """Across independent null datasets, the participant bootstrap 95%
    interval for the AUC penalty should cover zero about 95% of the time."""
    covered = []
    for s in range(10):
        src, tgt = make_studies(3000, 1000, Gap(), np.random.default_rng(500 + s))
        iv = resample.bootstrap(src, tgt, reps=100, seed=s)["auc"]
        covered.append(iv.covers(0.0))
    assert sum(covered) >= 8, covered


# AC2 ---------------------------------------------------------------------

def test_ac2_penalty_tracks_mechanism_dial():
    dials = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
    means = [simulate_runs(Gap(mechanism=d), reps=30, seed0=100)["auc"].mean() for d in dials]
    rho = np.corrcoef(dials, means)[0, 1]
    assert rho >= 0.9, dict(zip(dials, np.round(means, 4)))
    assert means[-1] > means[0] + 0.05


# AC3 ---------------------------------------------------------------------

def test_ac3_population_shift_moves_calibration_not_discrimination():
    shifts = [0.0, 0.5, 1.0, 1.5]
    runs = {d: simulate_runs(Gap(population=d), reps=30, seed0=300) for d in shifts}
    auc = {d: r["auc"].mean() for d, r in runs.items()}
    citl = {d: r["cal_intercept"].mean() for d, r in runs.items()}
    assert abs(auc[1.5] - auc[0.0]) < 0.02, auc
    for d in shifts:
        # Calibration-in-the-large recovers the shift; the spec's "tracks the
        # shift" made concrete as within 0.1 on the logit scale.
        assert citl[d] == pytest.approx(d, abs=0.1), citl


# AC4 ---------------------------------------------------------------------

@pytest.mark.parametrize("kind", ["simulation", "bootstrap"])
def test_ac4_small_target_widens_interval(kind):
    widths = {}
    for n in (1500, 30):
        if kind == "simulation":
            runs = simulate_runs(Gap(mechanism=0.4), reps=80, seed0=700, n_target=n)
            lo, hi = runs["auc"].dropna().quantile([0.025, 0.975])
            widths[n] = hi - lo
        else:
            src, tgt = make_studies(4000, n, Gap(mechanism=0.4), np.random.default_rng(900))
            widths[n] = resample.bootstrap(src, tgt, reps=150, seed=1)["auc"].width
    assert widths[30] >= 3 * widths[1500], widths


# Attack rate: separating "different epidemic" from "different people" -----

def _placebo(gap, n, rng):
    from sim.generate import placebo_incidence
    ps = placebo_incidence(n, Gap(), rng)
    pt = placebo_incidence(n, gap, rng)
    return transfer.PlaceboArms(round(ps * n), n, round(pt * n), n)


@pytest.mark.parametrize("gap, raw, adjusted", [
    (Gap(attack_rate=1.0), 1.0, 0.0),            # exposure only: offset absorbs it
    (Gap(population=1.0), 1.0, 1.0),             # vaccinee-only shift: offset cannot
    (Gap(attack_rate=1.0, population=0.5), 1.5, 0.5),
])
def test_placebo_offset_separates_exposure_from_population(gap, raw, adjusted):
    rows = []
    for s in range(20):
        rng = np.random.default_rng(1200 + s)
        src, tgt = make_studies(4000, 3000, gap, rng)
        r = transfer.run(src, tgt, seed=s, placebo=_placebo(gap, 20_000, rng))
        rows.append((r.penalty["cal_intercept"], r.penalty_adjusted["cal_intercept"],
                     r.penalty["auc"] - r.penalty_adjusted["auc"]))
    a = np.array(rows).mean(axis=0)
    assert a[0] == pytest.approx(raw, abs=0.1)
    assert a[1] == pytest.approx(adjusted, abs=0.1)
    assert abs(a[2]) < 1e-12  # a constant offset never changes the ranking


def test_bootstrap_reports_exposure_adjusted_intervals():
    rng = np.random.default_rng(3)
    gap = Gap(attack_rate=1.0)
    src, tgt = make_studies(3000, 1500, gap, rng)
    ivs = resample.bootstrap(src, tgt, reps=60, placebo=_placebo(gap, 5000, rng))
    assert not ivs["cal_intercept"].covers(0.0)
    assert ivs["cal_intercept_exposure_adjusted"].covers(0.0)


def test_cv_baseline_null_and_population_shift():
    """The cross-validated baseline used for small sources obeys AC1 and AC3."""
    null = pd.DataFrame([transfer.run(*make_studies(600, 600, Gap(), np.random.default_rng(1500 + s)),
                                      baseline="cv", seed=s).penalty for s in range(12)])
    assert abs(null["auc"].mean()) < 0.02 and abs(null["cal_intercept"].mean()) < 0.15
    shift = pd.DataFrame([transfer.run(*make_studies(600, 600, Gap(population=1.0),
                                                     np.random.default_rng(1600 + s)),
                                       baseline="cv", seed=s).penalty for s in range(12)])
    assert shift["cal_intercept"].mean() == pytest.approx(1.0, abs=0.15)


def test_cv_baseline_auc_is_not_biased_below_half_under_no_signal():
    """With markers unrelated to the outcome, the cross-validated baseline
    AUC must sit near 0.5, not below it (the pooled-CV artefact)."""
    aucs = []
    for s in range(15):
        rng = np.random.default_rng(1700 + s)
        src, tgt = make_studies(300, 300, Gap(), rng)
        src["endpoint"] = rng.binomial(1, 0.15, len(src))  # break the marker link
        aucs.append(transfer.run(src, tgt, baseline="cv", seed=s).baseline["auc"])
    assert abs(np.mean(aucs) - 0.5) < 0.03, np.round(aucs, 3)
