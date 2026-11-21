"""Weighted performance metrics.

Every metric takes ``w``, the inverse-probability sampling weights. With all
weights equal they reduce to the usual unweighted forms. Under two-phase
sampling the unweighted forms estimate the wrong population quantity and the
bias is not small, so there is deliberately no unweighted code path.

Calibration is reported three ways (Steyerberg 2019; Van Calster et al. 2016):

* ``cal_intercept``: calibration-in-the-large, the intercept of a logistic
  model of y with the prediction's logit as an offset (slope fixed at 1).
  0 is perfect; positive means the model under-predicts risk.
* ``cal_slope``: slope of y on the prediction's logit. 1 is perfect; below 1
  means predictions are too extreme (over-confident).
* ``cal_intercept_joint``: the intercept from the same fit as the slope. This
  is what the spike reported. It equals ``cal_intercept`` only when the slope
  is 1, so the two are labelled separately.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import brentq
from scipy.special import expit, logit

EPS = 1e-9
METRICS = ("auc", "brier", "cal_intercept", "cal_slope", "cal_intercept_joint")


def _as_arrays(y, p, w):
    y = np.asarray(y, float)
    p = np.asarray(p, float)
    w = np.ones_like(y) if w is None else np.asarray(w, float)
    if not (y.shape == p.shape == w.shape):
        raise ValueError("y, p and w must have the same shape")
    return y, p, w


def weighted_auc(y, score, w=None) -> float:
    """Weighted probability that a random case outranks a random non-case.

    Pairs are weighted by ``w_case * w_noncase``; ties count one half. This is
    the IPW AUC for two-phase designs. O(n log n).
    """
    y, s, w = _as_arrays(y, score, w)
    cases, ctrls = y == 1, y == 0
    if not cases.any() or not ctrls.any():
        return np.nan
    order = np.argsort(s, kind="mergesort")
    s, y, w = s[order], y[order], w[order]
    # Group tied scores; for each group, control weight strictly below it.
    uniq, start = np.unique(s, return_index=True)
    w_ctrl = np.where(y == 0, w, 0.0)
    w_case = np.where(y == 1, w, 0.0)
    grp_ctrl = np.add.reduceat(w_ctrl, start)
    grp_case = np.add.reduceat(w_case, start)
    below = np.concatenate([[0.0], np.cumsum(grp_ctrl)[:-1]])
    num = np.sum(grp_case * (below + 0.5 * grp_ctrl))
    return float(num / (w_case.sum() * w_ctrl.sum()))


def weighted_brier(y, p, w=None) -> float:
    y, p, w = _as_arrays(y, p, w)
    return float(np.sum(w * (p - y) ** 2) / np.sum(w))


def _lp(p):
    return logit(np.clip(p, EPS, 1 - EPS))


def calibration_in_the_large(y, p, w=None) -> float:
    """Intercept a in logit P(y=1) = a + logit(p), slope fixed at 1."""
    y, p, w = _as_arrays(y, p, w)
    if len(np.unique(y)) < 2:
        return np.nan
    lp = _lp(p)

    def score(a):
        return np.sum(w * (y - expit(a + lp)))

    lo, hi = -20.0, 20.0
    if score(lo) * score(hi) > 0:
        return np.nan
    return float(brentq(score, lo, hi, xtol=1e-10))


def _weighted_logistic(x, y, w, offset=None, iters=50, tol=1e-10):
    """Unpenalised weighted logistic regression by Newton-Raphson.

    ``x`` excludes the intercept column; it is added here. Returns the
    coefficient vector [intercept, *slopes]. Small and dependency-free so the
    calibration metrics do not inherit any library's default penalty.
    """
    X = np.column_stack([np.ones(len(y)), x])
    off = 0.0 if offset is None else offset
    b = np.zeros(X.shape[1])
    for _ in range(iters):
        mu = expit(X @ b + off)
        g = X.T @ (w * (y - mu))
        H = (X * (w * mu * (1 - mu))[:, None]).T @ X
        step = np.linalg.solve(H + 1e-12 * np.eye(len(b)), g)
        b = b + step
        if np.max(np.abs(step)) < tol:
            break
    return b


def calibration_slope(y, p, w=None) -> tuple[float, float]:
    """Return (slope, joint intercept) of y regressed on logit(p)."""
    y, p, w = _as_arrays(y, p, w)
    if len(np.unique(y)) < 2:
        return np.nan, np.nan
    b = _weighted_logistic(_lp(p)[:, None], y, w)
    if not np.all(np.isfinite(b)):
        return np.nan, np.nan
    return float(b[1]), float(b[0])


def evaluate(y, p, w=None) -> dict[str, float]:
    """All binary-endpoint metrics for predictions ``p`` against outcomes ``y``."""
    slope, joint = calibration_slope(y, p, w)
    return {
        "auc": weighted_auc(y, p, w),
        "brier": weighted_brier(y, p, w),
        "cal_intercept": calibration_in_the_large(y, p, w),
        "cal_slope": slope,
        "cal_intercept_joint": joint,
    }


# Sign conventions. For AUC and Brier the penalty is oriented so that positive
# means worse in the target. For the calibration metrics it is the signed shift
# target - baseline, because the direction is the finding (an intercept shift of
# +0.8 says the borrowed model under-predicts target risk) and a signed
# difference keeps the bootstrap well behaved near zero, which |.| does not.

def penalty(baseline: dict[str, float], target: dict[str, float]) -> dict[str, float]:
    return {
        "auc": baseline["auc"] - target["auc"],
        "brier": target["brier"] - baseline["brier"],
        "cal_intercept": target["cal_intercept"] - baseline["cal_intercept"],
        "cal_slope": target["cal_slope"] - baseline["cal_slope"],
        "cal_intercept_joint": target["cal_intercept_joint"] - baseline["cal_intercept_joint"],
    }
