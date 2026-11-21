"""Loaders that turn a dataset, where it lives, into the tidy schema.

Nothing here contains data. Each loader takes a path outside the repository
(AC7) and returns ``(frame, placebo)`` where ``frame`` is the vaccine arm in
the tidy schema and ``placebo`` is the placebo arm's per-stratum endpoint
counts, used for the exposure (attack-rate) offset.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from penalty import schema


def _read(path: Path) -> pd.DataFrame:
    if path.suffix.lower() in (".rda", ".rdata"):
        try:
            import pyreadr
        except ImportError as e:  # pragma: no cover
            raise ImportError("reading .rda needs `pip install pyreadr`, or export to CSV first") from e
        return next(iter(pyreadr.read_r(str(path)).values()))
    return pd.read_csv(path)


def stratify(values: pd.Series, rule: dict) -> pd.Series:
    """``rule`` = {"cut": x, "labels": [below, at_or_above]}."""
    lo, hi = rule.get("labels", ["below", "at_or_above"])
    return pd.Series(np.where(values.to_numpy(float) >= float(rule["cut"]), hi, lo), index=values.index)


def hvtn505(path: str | Path, split: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    """HVTN 505 (DNA/rAd5 HIV vaccine), as shipped in the CRAN ``vaccine``
    package (``data/hvtn505.rda``) or exported from SCHARP Atlas.

    Endpoint: HIV-1 diagnosis on/after week 28 before unblinding. The
    case-control phase two (``casecontrol == 1``) carries the IgG markers and
    the released IPW weights ``wt``. Phase-one-only vaccinees are kept with
    ``phase2 = 0``. The released file has no sampling-stratum variable, so the
    sampling stratum is taken as (endpoint, weight value): the weights are
    constant within the design's strata, so this recovers them among the
    sampled. Unsampled controls cannot be placed in a stratum and go in their
    own group; the bootstrap therefore resamples phase two within strata and
    phase one as a block, which is an approximation the report states.

    ``split`` = {"column": "age" | "BMI" | "bhvrisk", "cut": ..., "labels": [...]}.
    """
    raw = _read(Path(path))
    raw = raw.assign(stratum=stratify(raw[split["column"]], split))
    vax, pla = raw[raw["trt"] == 1], raw[raw["trt"] == 0]
    ph2 = vax["casecontrol"].astype(int)
    wt = vax["wt"].where(ph2 == 1, 1.0).astype(float)
    y = vax["HIVwk28preunbl"].astype(int)
    df = pd.DataFrame({
        "participant_id": vax["pub_id"].astype(str),
        "study_id": "HVTN505",
        "stratum": vax["stratum"],
        "platform": "DNA/rAd5",
        "marker_IgG_env": vax["IgG_env"].astype(float),
        "marker_IgG_V2": vax["IgG_V2"].astype(float),
        "marker_IgG_V3": vax["IgG_V3"].astype(float),
        "endpoint": y,
        "time_to_event": vax["HIVwk28preunblfu"].astype(float),
        "sampling_weight": wt,
        "phase2": ph2,
        "sampling_stratum": np.where(ph2 == 1, y.astype(str) + ":" + wt.round(6).astype(str), "phase1_only"),
    }).reset_index(drop=True)
    placebo = (pla.groupby("stratum")["HIVwk28preunbl"].agg(events="sum", n="count")
               .astype(int).reset_index())
    return schema.validate(df), placebo


LOADERS = {"hvtn505": hvtn505}
