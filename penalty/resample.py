"""Uncertainty for a transfer penalty.

``bootstrap`` is the participant-level bootstrap of one fixed dataset, which
is what real data needs. Each replicate resamples source and target
separately, at phase one when the phase-one cohort is in the frame
(literature review, section 6: resampling phase-two rows only gives
intervals that are too narrow; see ``_resample`` for the fallbacks). Each
resampled participant keeps their phase-two status and weight. The whole
operation is then rerun: split, fit, evaluate. Replicates with too few cases
to fit are dropped and counted, never silently imputed.

By default the penalty strength is fixed at the value chosen on the original
source fit, so replicates cost one fit each. ``retune=True`` reruns the
source-only cross-validation in every replicate, which also propagates tuning
uncertainty and costs about 50 times as much.

``leave_one_study_out`` applies when three or more studies exist.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from . import metrics, schema, transfer


@dataclass(frozen=True)
class Interval:
    estimate: float
    lo: float
    hi: float
    reps: np.ndarray

    @property
    def width(self) -> float:
        return self.hi - self.lo

    def covers(self, value: float) -> bool:
        return self.lo <= value <= self.hi


def _resample(df: pd.DataFrame, rng) -> pd.DataFrame:
    """One bootstrap draw of a frame, choosing the unit by what the frame holds.

    * Whole phase-one cohort present (some ``phase2 == 0``): resample the
      cohort i.i.d. and keep each person's phase-two status and weight. The
      number of cases then varies as it would in a repeat of the study, which
      is most of the uncertainty in a calibration intercept or Brier score.
    * Phase-two rows only, with ``sampling_stratum``: resample within sampling
      strata, the usual approximation when phase one was not released.
    * Otherwise (a simple cohort): resample i.i.d.
    """
    if "phase2" in df.columns and (df["phase2"] == 0).any():
        pick = rng.integers(0, len(df), len(df))
    elif "sampling_stratum" in df.columns:
        parts = [rng.choice(idx, size=len(idx), replace=True)
                 for idx in df.groupby("sampling_stratum", sort=False).indices.values()]
        pick = np.concatenate(parts)
    else:
        pick = rng.integers(0, len(df), len(df))
    out = df.iloc[pick].copy()
    # Copies of one person get distinct ids (id#k) for the overlap guard;
    # transfer.split_source and assign_folds keep the copies together.
    out["participant_id"] = out["participant_id"] + "#" + np.arange(len(out)).astype(str)
    return out


def bootstrap(source: pd.DataFrame, target: pd.DataFrame, reps: int = 1000, alpha: float = 0.05,
              seed: int = 0, retune: bool = False, **kw) -> dict[str, Interval]:
    """Percentile intervals for every penalty metric, plus the exposure-adjusted
    ones when ``placebo`` is given (placebo counts are resampled binomially)."""
    point = transfer.run(source, target, seed=seed, **kw)
    if not retune and point.model.kind == "ridge_logistic":
        kw = {**kw, "model_kw": {**(kw.get("model_kw") or {}), "c_grid": (point.model.tuning["C"],)}}
    source, target = schema.validate(source), schema.validate(target)
    placebo = kw.pop("placebo", None)
    rng = np.random.default_rng(seed)
    rows = []
    for r in range(reps):
        pl = placebo.resample(rng) if placebo is not None else None
        try:
            res = transfer.run(_resample(source, rng), _resample(target, rng), seed=seed + 1 + r,
                               placebo=pl, **kw)
        except ValueError:
            # Too few cases in this draw to fit or to tune; it counts as a
            # failed replicate, reported through Interval.reps' length.
            continue
        rows.append(res.all_penalties())
    draws = pd.DataFrame(rows, columns=list(point.all_penalties()))
    out = {}
    for k, est in point.all_penalties().items():
        v = draws[k].dropna().to_numpy()
        lo, hi = np.quantile(v, [alpha / 2, 1 - alpha / 2]) if len(v) else (np.nan, np.nan)
        out[k] = Interval(float(est), float(lo), float(hi), v)
    return out


def leave_one_study_out(df: pd.DataFrame, marker_map: schema.MarkerMap, **kw) -> pd.DataFrame:
    """Fit on all studies but one, evaluate on the one left out."""
    studies = sorted(df["study_id"].unique())
    if len(studies) < 3:
        raise ValueError("leave-one-study-out needs three or more studies; report pairwise instead")
    rows = []
    for s in studies:
        res = transfer.run(df[df.study_id != s], df[df.study_id == s], marker_map=marker_map, **kw)
        rows.append({"held_out_study": s, **res.all_penalties()})
    return pd.DataFrame(rows)


__all__ = ["Interval", "bootstrap", "leave_one_study_out", "metrics"]
