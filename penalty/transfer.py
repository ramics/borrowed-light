"""The one operation: fit on source, evaluate at home and across the gap.

```
fit        m = model(S_train)
baseline   b = evaluate(m, S_heldout)
transfer   t = evaluate(m, T)
penalty    = b - t           (per metric; sign conventions in metrics.py)
```

The model function receives ``S_train`` and nothing else, so target data
cannot reach the fit or the penalty tuning. As a second guard, evaluation sets
are checked for any participant the model was fitted on.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.special import expit, logit

from . import metrics, schema, survival
from .models import MODELS, FittedModel

MODELS = {**MODELS, "ridge_cox": survival.fit_ridge_cox}


@dataclass(frozen=True)
class PlaceboArms:
    """Placebo-arm endpoint counts in source and target.

    Placebo incidence measures exposure (the attack rate) without any vaccine
    or marker effect. Its logit difference is added to the borrowed model's
    linear predictor in the target, so the adjusted penalty reports what is
    left once "different epidemic" is taken out. Exact when exposure shifts
    the log-odds additively, as in the simulator; under a multiplicative
    hazard (cloglog) model it is an approximation that is good for rare
    endpoints.
    """

    source_events: int
    source_n: int
    target_events: int
    target_n: int

    def offset(self, link: str = "logit") -> float:
        """Shift in exposure on the model's scale: log-odds for the logistic
        path, log cumulative hazard (cloglog) for Cox."""
        ps = (self.source_events + 0.5) / (self.source_n + 1)
        pt = (self.target_events + 0.5) / (self.target_n + 1)
        if link == "cloglog":
            return float(np.log(-np.log1p(-pt)) - np.log(-np.log1p(-ps)))
        return float(logit(pt) - logit(ps))

    def resample(self, rng) -> "PlaceboArms":
        return PlaceboArms(
            int(rng.binomial(self.source_n, self.source_events / self.source_n)), self.source_n,
            int(rng.binomial(self.target_n, self.target_events / self.target_n)), self.target_n,
        )


@dataclass(frozen=True)
class TransferResult:
    model: FittedModel
    baseline: dict
    target: dict
    penalty: dict
    n: dict
    target_adjusted: dict | None = None
    penalty_adjusted: dict | None = None

    def all_penalties(self) -> dict:
        out = dict(self.penalty)
        if self.penalty_adjusted is not None:
            out.update({f"{k}_exposure_adjusted": v for k, v in self.penalty_adjusted.items()})
        return out


def _units(df: pd.DataFrame) -> tuple[np.ndarray, list[str]]:
    """Resampling units. A bootstrap draw can contain one person several
    times under suffixed ids (``id#k``); all copies must land on the same side
    of any split, or the baseline is measured partly on training rows."""
    base = df["participant_id"].str.split("#", n=1).str[0].to_numpy()
    keys = ["endpoint"] + (["sampling_stratum"] if "sampling_stratum" in df.columns else [])
    return base, keys


def split_source(source: pd.DataFrame, holdout_frac: float, rng) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Random participant-level split, stratified by endpoint and sampling
    stratum. Rows are split at phase one, so each half is itself a valid
    two-phase sample and the original weights stay correct in both."""
    base, keys = _units(source)
    held = np.zeros(len(source), bool)
    for _, idx in source.groupby(keys, sort=False).indices.items():
        people = rng.permutation(np.unique(base[idx]))
        chosen = set(people[: int(round(holdout_frac * len(people)))])
        held[idx] = [b in chosen for b in base[idx]]
    return source.loc[~held], source.loc[held]


def assign_folds(source: pd.DataFrame, k: int, rng) -> np.ndarray:
    """Stratified fold labels, dealt round-robin per person so every fold gets
    its share of cases, with the position carried between strata."""
    base, keys = _units(source)
    fold = np.zeros(len(source), int)
    pos = 0
    for _, idx in source.groupby(keys, sort=False).indices.items():
        people = rng.permutation(np.unique(base[idx]))
        lab = {b: (i + pos) % k for i, b in enumerate(people)}
        fold[idx] = [lab[b] for b in base[idx]]
        pos += len(people)
    return fold


def _predict(model: FittedModel, df: pd.DataFrame, offset: float = 0.0,
             horizon: float | None = None) -> pd.DataFrame:
    """Per-person predictions plus everything the metrics need."""
    overlap = model.fitted_ids & set(schema.phase_two(df)["participant_id"])
    if overlap:
        raise RuntimeError(f"{len(overlap)} evaluation participants were used to fit the model")
    d = schema.phase_two(df)
    out = pd.DataFrame({"y": d["endpoint"].to_numpy(), "w": d["sampling_weight"].to_numpy(),
                        "lp": model.linear_predictor(df) + offset})
    if model.kind == "ridge_cox":
        t = d["time_to_event"].to_numpy()
        out["t"] = t
        out["risk"] = survival.predicted_risk(model, df, horizon, offset=offset)
        out["cumhaz"] = survival.predicted_cumhaz(model, df, t, offset=offset)
    else:
        out["risk"] = expit(out["lp"])
    return out


def _score(pred: pd.DataFrame, cox: bool, horizon: float | None) -> dict:
    if cox:
        return survival.evaluate_survival(pred.t, pred.y, pred.lp, pred.risk, pred.w, horizon, pred.cumhaz)
    return metrics.evaluate(pred.y.to_numpy(), pred.risk.to_numpy(), pred.w.to_numpy())


def _evaluate(model: FittedModel, df: pd.DataFrame, offset: float = 0.0,
              horizon: float | None = None) -> dict:
    return _score(_predict(model, df, offset, horizon), model.kind == "ridge_cox", horizon)


def run(source: pd.DataFrame, target: pd.DataFrame, model: str = "ridge_logistic",
        holdout_frac: float = 0.3, seed: int = 0, marker_map: schema.MarkerMap | None = None,
        markers=None, placebo: PlaceboArms | None = None, model_kw: dict | None = None,
        horizon: float | None = None, baseline: str = "holdout", folds: int = 5) -> TransferResult:
    """Run the operation.

    ``baseline="holdout"`` (the spec's form) fits on a random ``1 - holdout_frac``
    of the source and measures the baseline on the rest. ``baseline="cv"`` is
    for small sources, where a single holdout of a handful of cases makes the
    baseline mostly noise: the baseline is measured on pooled out-of-fold
    predictions from ``folds``-fold cross-validation (each fold model tuned
    by its own inner CV), and the model carried to the target is fitted on
    the whole source. Either way nothing about the target touches any fit.
    """
    source = schema.validate(source)
    target = schema.validate(target)
    schema.require_alignment(source, target, marker_map)
    if marker_map is not None:
        source = schema.validate(marker_map.apply(source))
        target = schema.validate(marker_map.apply(target))
    if baseline == "cv":
        return run_cv(source, target, model=model, seed=seed, markers=markers, placebo=placebo,
                      model_kw=model_kw, horizon=horizon, folds=folds)
    if baseline != "holdout":
        raise ValueError(f"baseline must be 'holdout' or 'cv', not {baseline!r}")
    rng = np.random.default_rng(seed)
    train, held = split_source(source, holdout_frac, rng)
    return run_split(train, held, target, model=model, seed=seed, markers=markers,
                     placebo=placebo, model_kw=model_kw, horizon=horizon)


def run_split(train: pd.DataFrame, held: pd.DataFrame, target: pd.DataFrame,
              model: str = "ridge_logistic", seed: int = 0, markers=None,
              placebo: PlaceboArms | None = None, model_kw: dict | None = None,
              horizon: float | None = None) -> TransferResult:
    """The operation on an already-split source. Frames must be validated and
    on a common marker scale; :func:`run` does both.

    ``horizon`` is the time at which Cox risk, calibration and Brier are read.
    It defaults to the median follow-up time in the source training frame, so
    it too is fixed before the target is looked at."""
    m = MODELS[model](train, markers=markers, seed=seed, **(model_kw or {}))
    cox = m.kind == "ridge_cox"
    if cox and horizon is None:
        horizon = float(np.median(schema.phase_two(train)["time_to_event"]))
    pen = survival.penalty_survival if cox else metrics.penalty
    b = _evaluate(m, held, horizon=horizon)
    t = _evaluate(m, target, horizon=horizon)
    n = {name: int(len(schema.phase_two(d))) for name, d in
         (("train", train), ("heldout", held), ("target", target))}
    ta = pa = None
    if placebo is not None:
        ta = _evaluate(m, target, offset=placebo.offset("cloglog" if cox else "logit"), horizon=horizon)
        pa = pen(b, ta)
    return TransferResult(model=m, baseline=b, target=t, penalty=pen(b, t), n=n,
                          target_adjusted=ta, penalty_adjusted=pa)


def run_cv(source: pd.DataFrame, target: pd.DataFrame, model: str = "ridge_logistic", seed: int = 0,
           markers=None, placebo: PlaceboArms | None = None, model_kw: dict | None = None,
           horizon: float | None = None, folds: int = 5) -> TransferResult:
    """Cross-validated baseline; see :func:`run`. Frames must be validated."""
    fit = lambda df, s: MODELS[model](df, markers=markers, seed=s, **(model_kw or {}))  # noqa: E731
    m = fit(source, seed)
    cox = m.kind == "ridge_cox"
    if cox and horizon is None:
        horizon = float(np.median(schema.phase_two(source)["time_to_event"]))
    pen = survival.penalty_survival if cox else metrics.penalty
    ph2 = schema.phase_two(source)
    k = int(min(folds, ph2["endpoint"].sum(), (1 - ph2["endpoint"]).sum()))
    if k < 2:
        raise ValueError("need at least two cases and two non-cases for a cross-validated baseline")
    fold = assign_folds(source, k, np.random.default_rng(seed))
    oof = []
    for f in range(k):
        mf = fit(source.loc[fold != f], seed + 1 + f)
        oof.append(_predict(mf, source.loc[fold == f], horizon=horizon))
    b = _score(pd.concat(oof, ignore_index=True), cox, horizon)
    # Ranking metrics are not pooled across folds. Each fold model has its own
    # intercept, and a fold with fewer cases in training gets a lower one, so
    # pooled out-of-fold scores rank people by fold and bias AUC below 0.5
    # when the marker signal is weak (Parker et al. 2007; Airola et al. 2011).
    # Within-fold AUC/C, averaged with case x non-case weight, avoids that.
    rank = "c_index" if cox else "auc"
    per_fold, pair_w = [], []
    for pr in oof:
        if cox:
            per_fold.append(survival.weighted_c_index(pr.t, pr.y, pr.lp, pr.w))
        else:
            per_fold.append(metrics.weighted_auc(pr.y, pr.lp, pr.w))
        pair_w.append(pr.w[pr.y == 1].sum() * pr.w[pr.y == 0].sum())
    ok = ~np.isnan(per_fold)
    b[rank] = float(np.average(np.asarray(per_fold)[ok], weights=np.asarray(pair_w)[ok]))
    t = _evaluate(m, target, horizon=horizon)
    n = {"train": int(len(ph2)), "heldout": int(len(ph2)), "target": int(len(schema.phase_two(target))),
         "folds": k}
    ta = pa = None
    if placebo is not None:
        ta = _evaluate(m, target, offset=placebo.offset("cloglog" if cox else "logit"), horizon=horizon)
        pa = pen(b, ta)
    return TransferResult(model=m, baseline=b, target=t, penalty=pen(b, t), n=n,
                          target_adjusted=ta, penalty_adjusted=pa)
