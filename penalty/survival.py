"""Time-to-event path: weighted ridge Cox and its metrics.

Weighted throughout, with the same rules as the binary path: ``sampling_weight``
enters the partial likelihood and every metric, and the penalty is tuned by
cross-validation on the source training frame only.

Ties are handled with Breslow's approximation. The baseline cumulative hazard
is the weighted Breslow estimator. Metrics, all weighted:

* ``c_index``    Harrell's C with IPW pair weights ``w_i * w_j``. Pairs are
  comparable when the earlier time is an event; tied risk scores count 1/2.
* ``cal_slope``  slope of a weighted Cox fit of the outcome on the borrowed
  linear predictor. 1 is perfect.
* ``cal_intercept`` log(observed / expected) events, with expected the sum
  of each person's predicted cumulative hazard at their own follow-up time
  (Crowson, Atkinson & Therneau 2016). It equals the intercept of a Poisson
  model with the log cumulative hazard as offset, so it is the time-to-event
  analogue of calibration-in-the-large: 0 is perfect, positive means the
  model under-predicts, and a pure log-hazard shift is recovered one for one.
* ``brier``      Brier score at ``horizon`` with inverse-probability-of-
  censoring weights (Graf et al. 1999), multiplied by the sampling weights.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from sklearn.model_selection import StratifiedKFold

from . import schema
from .models import DEFAULT_C_GRID, FittedModel, _weighted_standardiser

SURV_METRICS = ("c_index", "brier", "cal_intercept", "cal_slope")


# ----------------------------------------------------------------------------
# Weighted Cox partial likelihood (Breslow ties)
# ----------------------------------------------------------------------------

def _neg_loglik(beta, z, t, d, w, lam):
    eta = z @ beta
    order = np.argsort(-t, kind="mergesort")  # descending time: risk sets are prefixes
    eta_o, t_o, d_o, w_o, z_o = eta[order], t[order], d[order], w[order], z[order]
    m = eta_o.max()
    r = w_o * np.exp(eta_o - m)
    S0 = np.cumsum(r)
    S1 = np.cumsum(r[:, None] * z_o, axis=0)
    # For Breslow ties every subject at time t shares the risk set of all with
    # time >= t, i.e. the cumulative sum up to the *last* index with that time.
    _, first = np.unique(-t_o, return_index=True)
    last = np.append(first[1:] - 1, len(t_o) - 1)
    grp = np.repeat(np.arange(len(first)), np.diff(np.append(first, len(t_o))))
    S0t = S0[last][grp]
    S1t = S1[last][grp]
    ev = d_o == 1
    ll = np.sum(w_o[ev] * (eta_o[ev] - m - np.log(S0t[ev])))
    grad = np.sum(w_o[ev, None] * (z_o[ev] - S1t[ev] / S0t[ev, None]), axis=0)
    return -ll + 0.5 * lam * beta @ beta, -grad + lam * beta


def fit_cox(z, t, d, w, lam=0.0, beta0=None):
    """Weighted (ridge) Cox. Weights are normalised to mean 1 so ``lam`` has
    the same meaning at any sampling fraction, matching the logistic path."""
    w = w / w.mean()
    beta0 = np.zeros(z.shape[1]) if beta0 is None else beta0
    res = minimize(_neg_loglik, beta0, args=(z, t, d, w, lam), jac=True, method="L-BFGS-B",
                   options={"maxiter": 2000, "gtol": 1e-8})
    return res.x


def breslow_cumhaz(eta, t, d, w):
    """Weighted Breslow baseline cumulative hazard; returns (times, H0)."""
    times = np.unique(t[d == 1])
    r = w * np.exp(eta)
    H = np.empty(len(times))
    acc = 0.0
    for i, s in enumerate(times):
        acc += np.sum(w[(t == s) & (d == 1)]) / np.sum(r[t >= s])
        H[i] = acc
    return times, H


def _step(times, values, at, before=0.0):
    idx = np.searchsorted(times, at, side="right") - 1
    return np.where(idx >= 0, values[np.clip(idx, 0, None)], before)


def weighted_km(t, d, w):
    """Weighted Kaplan-Meier survival; returns (times, S)."""
    times = np.unique(t[d == 1])
    S = np.empty(len(times))
    s = 1.0
    for i, u in enumerate(times):
        at_risk = np.sum(w[t >= u])
        s *= 1 - np.sum(w[(t == u) & (d == 1)]) / at_risk
        S[i] = s
    return times, S


# ----------------------------------------------------------------------------
# Model
# ----------------------------------------------------------------------------

def _design(df, markers):
    d = schema.phase_two(df)
    if "time_to_event" not in d.columns:
        raise schema.SchemaError("ridge Cox needs a time_to_event column")
    return (d[list(markers)].to_numpy(float), d["time_to_event"].to_numpy(float),
            d["endpoint"].to_numpy(int), d["sampling_weight"].to_numpy(float))


def _cv_partial_loglik(beta, z_tr, t_tr, d_tr, w_tr, z_all, t_all, d_all, w_all):
    """Verweij & van Houwelingen (1993) CV partial likelihood: the held-out
    fold's contribution is l(all) - l(train) evaluated at the training fit."""
    la, _ = _neg_loglik(beta, z_all, t_all, d_all, w_all / w_all.mean(), 0.0)
    lt, _ = _neg_loglik(beta, z_tr, t_tr, d_tr, w_tr / w_all.mean(), 0.0)
    return la - lt  # negative held-out log-likelihood


def fit_ridge_cox(source: pd.DataFrame, markers=None, c_grid=DEFAULT_C_GRID, folds: int = 5,
                  seed: int = 0) -> FittedModel:
    """Weighted ridge Cox. ``c_grid`` is on the same 1/penalty scale as the
    logistic path (``lam = 1 / C``), and is chosen by weighted CV partial
    likelihood inside the source only."""
    markers = tuple(markers or schema.marker_columns(source))
    x, t, d, w = _design(source, markers)
    centre, scale = _weighted_standardiser(x, w)
    z = (x - centre) / scale
    losses = np.zeros(len(c_grid))
    if len(c_grid) > 1:
        k = int(min(folds, d.sum()))
        if k < 2:
            raise ValueError("need at least two events to tune the penalty")
        for tr, _ in StratifiedKFold(n_splits=k, shuffle=True, random_state=seed).split(z, d):
            for i, C in enumerate(c_grid):
                b = fit_cox(z[tr], t[tr], d[tr], w[tr], lam=1.0 / C)
                losses[i] += _cv_partial_loglik(b, z[tr], t[tr], d[tr], w[tr], z, t, d, w)
    best = int(np.argmin(losses))
    beta = fit_cox(z, t, d, w, lam=1.0 / c_grid[best])
    times, H0 = breslow_cumhaz(z @ beta, t, d, w)
    return FittedModel(
        kind="ridge_cox", markers=markers, intercept=0.0, coef=beta, centre=centre, scale=scale,
        tuning={"C": float(c_grid[best]), "cv_neg_partial_loglik": losses.tolist(),
                "c_grid": [float(c) for c in c_grid], "baseline_times": times, "baseline_cumhaz": H0},
        fitted_ids=frozenset(schema.phase_two(source)["participant_id"]),
    )


def predicted_cumhaz(model: FittedModel, df: pd.DataFrame, at, offset: float = 0.0):
    """Predicted cumulative hazard for the phase-two rows of ``df`` at time(s) ``at``."""
    H0 = _step(model.tuning["baseline_times"], model.tuning["baseline_cumhaz"], at)
    return H0 * np.exp(model.linear_predictor(df) + offset)


def predicted_risk(model: FittedModel, df: pd.DataFrame, horizon: float, offset: float = 0.0):
    return 1 - np.exp(-predicted_cumhaz(model, df, horizon, offset))


# ----------------------------------------------------------------------------
# Metrics
# ----------------------------------------------------------------------------

def weighted_c_index(t, d, score, w) -> float:
    """Harrell's C with pair weights w_i * w_j. O(n^2) in memory-light chunks;
    fine at correlates-study sizes (a few thousand phase-two rows)."""
    t, d, s, w = map(np.asarray, (t, d, score, w))
    num = den = 0.0
    for i in np.where(d == 1)[0]:
        later = t > t[i]
        if not later.any():
            continue
        ww = w[i] * w[later]
        num += np.sum(ww * ((s[i] > s[later]) + 0.5 * (s[i] == s[later])))
        den += np.sum(ww)
    return float(num / den) if den else np.nan


def evaluate_survival(t, d, lp, risk, w, horizon, cumhaz) -> dict:
    """``cumhaz`` is each person's predicted cumulative hazard at their own
    time; ``risk`` is predicted risk by ``horizon``."""
    t, d, lp, risk, w, cumhaz = map(lambda a: np.asarray(a, float), (t, d, lp, risk, w, cumhaz))
    d = d.astype(int)
    out = {"c_index": weighted_c_index(t, d, lp, w)}
    # Calibration slope: weighted Cox of outcome on the borrowed linear predictor.
    if d.sum() >= 2:
        out["cal_slope"] = float(fit_cox(lp[:, None], t, d, w)[0])
    else:
        out["cal_slope"] = np.nan
    observed, expected = np.sum(w * d), np.sum(w * cumhaz)
    out["cal_intercept"] = float(np.log(observed / expected)) if observed > 0 and expected > 0 else np.nan
    # IPCW Brier at the horizon: censoring survival G from reverse KM.
    g_t, g_s = weighted_km(t, 1 - d, w)
    G_t = _step(g_t, g_s, t - 1e-12, before=1.0)
    G_h = _step(g_t, g_s, horizon, before=1.0)
    event_by = (t <= horizon) & (d == 1)
    alive = t > horizon
    ipcw = np.where(event_by, 1 / np.clip(G_t, 1e-8, None), np.where(alive, 1 / max(G_h, 1e-8), 0.0))
    resid = np.where(event_by, (1 - risk) ** 2, risk ** 2)
    out["brier"] = float(np.sum(w * ipcw * resid) / np.sum(w))
    return out


def penalty_survival(b: dict, t: dict) -> dict:
    return {
        "c_index": b["c_index"] - t["c_index"],
        "brier": t["brier"] - b["brier"],
        "cal_intercept": t["cal_intercept"] - b["cal_intercept"],
        "cal_slope": t["cal_slope"] - b["cal_slope"],
    }
