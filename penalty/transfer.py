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

from . import metrics, schema
from .models import MODELS, FittedModel


@dataclass(frozen=True)
class TransferResult:
    model: FittedModel
    baseline: dict
    target: dict
    penalty: dict
    n: dict


def split_source(source: pd.DataFrame, holdout_frac: float, rng) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Random participant-level split, stratified by endpoint and sampling
    stratum. Rows are split at phase one, so each half is itself a valid
    two-phase sample and the original weights stay correct in both."""
    keys = ["endpoint"] + (["sampling_stratum"] if "sampling_stratum" in source.columns else [])
    held = np.zeros(len(source), bool)
    for _, idx in source.groupby(keys, sort=False).indices.items():
        idx = rng.permutation(idx)
        held[idx[: int(round(holdout_frac * len(idx)))]] = True
    return source.loc[~held], source.loc[held]


def _evaluate(model: FittedModel, df: pd.DataFrame) -> dict:
    overlap = model.fitted_ids & set(schema.phase_two(df)["participant_id"])
    if overlap:
        raise RuntimeError(f"{len(overlap)} evaluation participants were used to fit the model")
    d = schema.phase_two(df)
    return metrics.evaluate(d["endpoint"].to_numpy(), model.predict(df), d["sampling_weight"].to_numpy())


def run(source: pd.DataFrame, target: pd.DataFrame, model: str = "ridge_logistic",
        holdout_frac: float = 0.3, seed: int = 0, marker_map: schema.MarkerMap | None = None,
        markers=None) -> TransferResult:
    source = schema.validate(source)
    target = schema.validate(target)
    schema.require_alignment(source, target, marker_map)
    if marker_map is not None:
        source = schema.validate(marker_map.apply(source))
        target = schema.validate(marker_map.apply(target))
    rng = np.random.default_rng(seed)
    train, held = split_source(source, holdout_frac, rng)
    return run_split(train, held, target, model=model, seed=seed, markers=markers)


def run_split(train: pd.DataFrame, held: pd.DataFrame, target: pd.DataFrame,
              model: str = "ridge_logistic", seed: int = 0, markers=None) -> TransferResult:
    """The operation on an already-split source. Frames must be validated and
    on a common marker scale; :func:`run` does both."""
    m = MODELS[model](train, markers=markers, seed=seed)
    b = _evaluate(m, held)
    t = _evaluate(m, target)
    n = {name: int(len(schema.phase_two(d))) for name, d in
         (("train", train), ("heldout", held), ("target", target))}
    return TransferResult(model=m, baseline=b, target=t, penalty=metrics.penalty(b, t), n=n)
