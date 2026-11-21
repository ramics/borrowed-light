"""Correlate models: weighted ridge logistic and the best-single-marker comparator.

All fitting is weighted by ``sampling_weight``. Penalty strength is chosen by
cross-validation inside the frame passed to ``fit``. A fitted model only ever
sees the source frame: there is no argument through which target data could
reach the tuning, and the fitted model records the participants it was fitted
on so :mod:`penalty.transfer` can refuse any overlap with the evaluation sets.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy.special import expit
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold

from . import schema
from .metrics import _weighted_logistic

DEFAULT_C_GRID = tuple(np.logspace(-3, 2, 11))


def _weighted_log_loss(y, p, w):
    p = np.clip(p, 1e-12, 1 - 1e-12)
    return float(-np.sum(w * (y * np.log(p) + (1 - y) * np.log(1 - p))) / np.sum(w))


def _weighted_standardiser(x, w):
    mu = np.average(x, axis=0, weights=w)
    sd = np.sqrt(np.average((x - mu) ** 2, axis=0, weights=w))
    sd[sd == 0] = 1.0
    return mu, sd


@dataclass(frozen=True)
class FittedModel:
    kind: str
    markers: tuple[str, ...]
    intercept: float
    coef: np.ndarray
    centre: np.ndarray
    scale: np.ndarray
    tuning: dict = field(default_factory=dict)
    fitted_ids: frozenset = field(default_factory=frozenset, repr=False)

    def linear_predictor(self, df: pd.DataFrame) -> np.ndarray:
        x = schema.phase_two(df)[list(self.markers)].to_numpy(float)
        return self.intercept + ((x - self.centre) / self.scale) @ self.coef

    def predict(self, df: pd.DataFrame) -> np.ndarray:
        """Predicted risk for the phase-two rows of ``df``, in row order."""
        return expit(self.linear_predictor(df))


def _fit_ridge(x, y, w, C):
    # sklearn's objective is C * sum(w_i * loss_i) + 0.5 ||beta||^2, so the
    # effective penalty depends on the total weight. Normalising the weights to
    # mean 1 makes C mean the same thing whatever the sampling fractions were.
    wn = w / w.mean()
    m = LogisticRegression(C=C, max_iter=5000, tol=1e-8).fit(x, y, sample_weight=wn)
    return float(m.intercept_[0]), m.coef_[0].copy()


def fit_ridge_logistic(source: pd.DataFrame, markers=None, c_grid=DEFAULT_C_GRID,
                       folds: int = 5, seed: int = 0) -> FittedModel:
    """Weighted ridge logistic regression, C chosen by weighted CV log loss."""
    markers = tuple(markers or schema.marker_columns(source))
    x, y, w = schema.design(source, list(markers))
    centre, scale = _weighted_standardiser(x, w)
    z = (x - centre) / scale

    k = int(min(folds, np.bincount(y, minlength=2).min()))
    if k < 2:
        raise ValueError("need at least two cases and two non-cases to tune the penalty")
    cv = StratifiedKFold(n_splits=k, shuffle=True, random_state=seed)
    losses = np.zeros(len(c_grid))
    for tr, te in cv.split(z, y):
        mu, sd = _weighted_standardiser(x[tr], w[tr])
        ztr, zte = (x[tr] - mu) / sd, (x[te] - mu) / sd
        for i, C in enumerate(c_grid):
            b0, b = _fit_ridge(ztr, y[tr], w[tr], C)
            losses[i] += _weighted_log_loss(y[te], expit(b0 + zte @ b), w[te]) / k
    best = int(np.argmin(losses))
    b0, b = _fit_ridge(z, y, w, c_grid[best])
    return FittedModel(
        kind="ridge_logistic", markers=markers, intercept=b0, coef=b,
        centre=centre, scale=scale,
        tuning={"C": float(c_grid[best]), "cv_log_loss": losses.tolist(),
                "c_grid": [float(c) for c in c_grid], "folds": k},
        fitted_ids=frozenset(schema.phase_two(source)["participant_id"]),
    )


def fit_best_single_marker(source: pd.DataFrame, markers=None, folds: int = 5,
                           seed: int = 0) -> FittedModel:
    """The comparator: the one marker with the lowest weighted CV log loss,
    fitted by unpenalised weighted logistic regression."""
    markers = tuple(markers or schema.marker_columns(source))
    x, y, w = schema.design(source, list(markers))
    k = int(min(folds, np.bincount(y, minlength=2).min()))
    if k < 2:
        raise ValueError("need at least two cases and two non-cases to choose a marker")
    cv = StratifiedKFold(n_splits=k, shuffle=True, random_state=seed)
    losses = np.zeros(len(markers))
    for tr, te in cv.split(x, y):
        for j in range(len(markers)):
            b = _weighted_logistic(x[tr, j:j + 1], y[tr], w[tr])
            losses[j] += _weighted_log_loss(y[te], expit(b[0] + b[1] * x[te, j]), w[te]) / k
    j = int(np.argmin(losses))
    b = _weighted_logistic(x[:, j:j + 1], y, w)
    return FittedModel(
        kind="single_marker", markers=(markers[j],), intercept=float(b[0]),
        coef=np.array([b[1]]), centre=np.zeros(1), scale=np.ones(1),
        tuning={"chosen": markers[j], "cv_log_loss": dict(zip(markers, losses.tolist()))},
        fitted_ids=frozenset(schema.phase_two(source)["participant_id"]),
    )


MODELS = {"ridge_logistic": fit_ridge_logistic, "single_marker": fit_best_single_marker}
