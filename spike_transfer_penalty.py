"""
Spike: can we estimate the penalty a correlate model pays when it is applied
across a gap (different study, platform or population)?

The whole point is that we simulate data where the true gap is set by a dial,
so we know the right answer and can check the estimator against it. No real
participant data is involved, here or ever in this repo.

Run:  python3 spike_transfer_penalty.py
"""

import numpy as np
from dataclasses import dataclass
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, brier_score_loss

RNG = np.random.default_rng(20261121)
P = 8  # number of immune markers


# ----------------------------------------------------------------------------
# Simulation: a source study and a target study separated by a controlled gap
# ----------------------------------------------------------------------------

@dataclass
class Gap:
    """Each dial is 0 when source and target are identical."""
    mechanism: float = 0.0   # radians; rotates the true coefficient vector
    assay: float = 0.0       # multiplicative scale shift on marker 0..2
    population: float = 0.0  # shift in baseline log-odds of the endpoint
    noise: float = 0.0       # extra measurement noise on markers in target


def _rotation(theta, p=P):
    """Rotate in the (0,1) plane. Stands in for 'the markers mean something
    slightly different in the target population'."""
    r = np.eye(p)
    r[0, 0] = r[1, 1] = np.cos(theta)
    r[0, 1] = -np.sin(theta)
    r[1, 0] = np.sin(theta)
    return r


def simulate(n, beta, intercept, gap: Gap, rng):
    cov = 0.3 ** np.abs(np.subtract.outer(np.arange(P), np.arange(P)))
    x = rng.multivariate_normal(np.zeros(P), cov, size=n)
    if gap.noise:
        x = x + rng.normal(0, gap.noise, size=x.shape)
    if gap.assay:
        x[:, :3] *= (1.0 + gap.assay)
    eta = x @ beta + intercept + gap.population
    y = rng.binomial(1, 1 / (1 + np.exp(-eta)))
    return x, y


def make_studies(n_source, n_target, gap: Gap, rng):
    beta_s = np.array([1.1, -0.8, 0.6, 0.0, 0.45, 0.0, -0.3, 0.2])
    beta_t = _rotation(gap.mechanism) @ beta_s
    intercept = -1.2
    xs, ys = simulate(n_source, beta_s, intercept, Gap(), rng)
    xt, yt = simulate(n_target, beta_t, intercept, gap, rng)
    return (xs, ys), (xt, yt)


# ----------------------------------------------------------------------------
# The estimator under test
# ----------------------------------------------------------------------------

def fit_correlate_model(x, y):
    """Ridge-penalised logistic regression, matching the stated approach of
    Repop-TB model 16 so results are comparable to the parent programme."""
    return LogisticRegression(C=1.0, max_iter=2000).fit(x, y)


def calibration(y, p):
    """Slope and intercept of y regressed on the logit of the prediction.
    Perfect calibration is slope 1, intercept 0. This is the part that matters
    for admissibility: a model can rank people correctly and still put the
    decision threshold in the wrong place."""
    p = np.clip(p, 1e-6, 1 - 1e-6)
    lp = np.log(p / (1 - p)).reshape(-1, 1)
    if len(np.unique(y)) < 2:
        return np.nan, np.nan
    m = LogisticRegression(C=np.inf, max_iter=2000).fit(lp, y)
    return float(m.coef_[0][0]), float(m.intercept_[0])


def evaluate(model, x, y):
    p = model.predict_proba(x)[:, 1]
    slope, icpt = calibration(y, p)
    return {
        "auc": roc_auc_score(y, p),
        "brier": brier_score_loss(y, p),
        "cal_slope": slope,
        "cal_intercept": icpt,
    }


def transfer_penalty(gap: Gap, n_source=4000, n_target=1500, rng=None):
    """Fit on source, evaluate on held-out source and on target.
    The penalty is the difference: what you lose by crossing the gap."""
    rng = rng or RNG
    (xs, ys), (xt, yt) = make_studies(n_source, n_target, gap, rng)
    cut = int(0.7 * n_source)
    model = fit_correlate_model(xs[:cut], ys[:cut])
    held_out = evaluate(model, xs[cut:], ys[cut:])
    target = evaluate(model, xt, yt)
    return {
        "auc_penalty": held_out["auc"] - target["auc"],
        "brier_penalty": target["brier"] - held_out["brier"],
        "cal_slope_target": target["cal_slope"],
        "cal_intercept_target": target["cal_intercept"],
        "auc_source": held_out["auc"],
        "auc_target": target["auc"],
    }


def bootstrap_interval(gap: Gap, reps=200, alpha=0.05, **kw):
    vals = [transfer_penalty(gap, rng=np.random.default_rng(s), **kw)["auc_penalty"]
            for s in range(reps)]
    lo, hi = np.quantile(vals, [alpha / 2, 1 - alpha / 2])
    return float(np.mean(vals)), float(lo), float(hi)


# ----------------------------------------------------------------------------
# Acceptance checks
# ----------------------------------------------------------------------------

def check_null():
    """AC1. With no gap the estimated penalty must be indistinguishable from
    zero. If this fails, every later number is a bug, not a finding."""
    mean, lo, hi = bootstrap_interval(Gap(), reps=200)
    covers = lo <= 0 <= hi
    print(f"AC1 null gap        mean penalty {mean:+.4f}  95% [{lo:+.4f}, {hi:+.4f}]  "
          f"covers zero: {covers}")
    return covers


def check_monotone():
    """AC2. As the mechanism gap widens the penalty must grow."""
    dials = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
    means = []
    print("\nAC2 mechanism dial sweep")
    print(f"  {'theta':>6}  {'AUC source':>10}  {'AUC target':>10}  "
          f"{'penalty':>9}  {'cal slope':>9}")
    for d in dials:
        reps = [transfer_penalty(Gap(mechanism=d), rng=np.random.default_rng(100 + s))
                for s in range(40)]
        pen = float(np.mean([r["auc_penalty"] for r in reps]))
        src = float(np.mean([r["auc_source"] for r in reps]))
        tgt = float(np.mean([r["auc_target"] for r in reps]))
        slope = float(np.mean([r["cal_slope_target"] for r in reps]))
        means.append(pen)
        print(f"  {d:6.1f}  {src:10.3f}  {tgt:10.3f}  {pen:+9.4f}  {slope:9.3f}")
    rho = float(np.corrcoef(dials, means)[0, 1])
    print(f"  correlation of dial with estimated penalty: {rho:.3f}")
    return rho >= 0.9


def check_calibration_breaks_first():
    """AC3. The claim worth testing: a population shift can leave ranking
    almost untouched while moving the threshold. If true, discrimination alone
    is the wrong thing to report to a regulator."""
    print("\nAC3 population shift: discrimination vs calibration")
    print(f"  {'shift':>6}  {'penalty(AUC)':>13}  {'cal intercept':>14}")
    rows = []
    for d in [0.0, 0.5, 1.0, 1.5]:
        reps = [transfer_penalty(Gap(population=d), rng=np.random.default_rng(300 + s))
                for s in range(40)]
        pen = float(np.mean([r["auc_penalty"] for r in reps]))
        icpt = float(np.mean([r["cal_intercept_target"] for r in reps]))
        rows.append((pen, icpt))
        print(f"  {d:6.1f}  {pen:+13.4f}  {icpt:+14.3f}")
    auc_stable = abs(rows[-1][0] - rows[0][0]) < 0.02
    cal_moved = abs(rows[-1][1] - rows[0][1]) > 0.5
    print(f"  AUC barely moves: {auc_stable}   calibration moves: {cal_moved}")
    return auc_stable and cal_moved


def check_small_n_humility():
    """AC4. With a first-in-human sized target cohort the interval must widen
    a lot. A method that stays confident at n=30 is lying."""
    print("\nAC4 target cohort size")
    print(f"  {'n':>5}  {'95% interval width':>19}")
    widths = {}
    for n in [1500, 200, 30]:
        _, lo, hi = bootstrap_interval(Gap(mechanism=0.4), reps=120, n_target=n)
        widths[n] = hi - lo
        print(f"  {n:5d}  {hi - lo:19.4f}")
    ok = widths[30] > 3 * widths[1500]
    print(f"  n=30 interval at least 3x the n=1500 interval: {ok}")
    return ok


if __name__ == "__main__":
    print("Transfer-penalty estimator — spike on synthetic data\n")
    results = {
        "AC1 null penalty covers zero": check_null(),
        "AC2 penalty tracks mechanism gap": check_monotone(),
        "AC3 calibration breaks before discrimination": check_calibration_breaks_first(),
        "AC4 small cohorts widen the interval": check_small_n_humility(),
    }
    print("\n" + "=" * 60)
    for k, v in results.items():
        print(f"  {'PASS' if v else 'FAIL'}  {k}")
    print("=" * 60)
