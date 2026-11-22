"""Tidy-frame validation, sampling-weight handling and marker mapping.

Every dataset, synthetic or real, enters the pipeline as one pandas frame with
the columns below. Real data never needs the code to change; it needs a loader
that produces this frame where the data lives.

| column            | dtype | notes                                              |
| ----------------- | ----- | -------------------------------------------------- |
| participant_id    | str   | unique                                             |
| study_id          | str   | split axis for E1                                  |
| stratum           | str   | split axis for E2                                  |
| platform          | str   | vaccine platform, may be "unknown"                 |
| marker_*          | float | one column per assay readout; NaN only off phase 2 |
| endpoint          | int   | 0/1                                                |
| time_to_event     | float | optional, for Cox                                  |
| sampling_weight   | float | required, > 0, inverse phase-two sampling prob.   |
| sampling_stratum  | str   | optional, two-phase sampling stratum               |
| phase2            | int   | optional 0/1; rows with 0 carry no marker data     |
| cluster_id        | str   | optional; rows sharing it are resampled together   |

Two-phase (case-cohort / case-control) data can be supplied either as the
phase-two rows only, with ``sampling_weight`` = 1 / P(selected | stratum), or
as the whole phase-one cohort with ``phase2`` marking who was measured. In the
second form the phase-one-only rows are kept so the bootstrap can resample the
cohort and keep each person's phase-two status (literature review, section 6),
but they are dropped before any fit or metric.

``cluster_id`` is for data where rows are not independent: one person
observed in several seasons, or members of one household. Rows sharing a
cluster always land on the same side of a holdout split or CV fold, and the
bootstrap resamples whole clusters, jointly across source and target when
the two share clusters (see :mod:`penalty.resample`).
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

MARKER_PREFIX = "marker_"
REQUIRED = ("participant_id", "study_id", "stratum", "platform", "endpoint", "sampling_weight")
OPTIONAL = ("time_to_event", "sampling_stratum", "phase2", "cluster_id")
STRING_COLS = ("participant_id", "study_id", "stratum", "platform", "sampling_stratum", "cluster_id")


class SchemaError(ValueError):
    """The frame does not satisfy the tidy schema."""


def marker_columns(df: pd.DataFrame) -> list[str]:
    return [c for c in df.columns if c.startswith(MARKER_PREFIX)]


def with_default_weights(df: pd.DataFrame) -> pd.DataFrame:
    """Add ``sampling_weight = 1.0`` when absent, loudly.

    The column is required by :func:`validate`. Defaulting is allowed only
    through this explicit call so a loader cannot drop the weights by accident.
    """
    if "sampling_weight" in df.columns:
        return df
    warnings.warn(
        "sampling_weight missing; assuming a simple random sample (all weights 1.0). "
        "If this is a case-cohort or case-control release, this is wrong.",
        stacklevel=2,
    )
    return df.assign(sampling_weight=1.0)


def validate(df: pd.DataFrame) -> pd.DataFrame:
    """Check the frame against the schema and return a normalised copy.

    Raises :class:`SchemaError` describing every problem found, not just the
    first one.
    """
    problems: list[str] = []
    missing = [c for c in REQUIRED if c not in df.columns]
    if missing:
        problems.append(f"missing required columns: {missing}")
    markers = marker_columns(df)
    if not markers:
        problems.append(f"no marker columns (expected names starting {MARKER_PREFIX!r})")
    if problems:
        raise SchemaError("; ".join(problems))

    out = df.copy()
    for c in STRING_COLS:
        if c in out.columns:
            out[c] = out[c].astype(str)

    if out["participant_id"].duplicated().any():
        dupes = out.loc[out["participant_id"].duplicated(), "participant_id"].head(5).tolist()
        problems.append(f"participant_id not unique, e.g. {dupes}")

    y = pd.to_numeric(out["endpoint"], errors="coerce")
    if y.isna().any() or not y.isin([0, 1]).all():
        problems.append("endpoint must be 0/1 with no missing values")
    else:
        out["endpoint"] = y.astype(int)

    w = pd.to_numeric(out["sampling_weight"], errors="coerce")
    if w.isna().any() or not np.isfinite(w).all() or (w <= 0).any():
        problems.append("sampling_weight must be finite and > 0 for every row")
    else:
        out["sampling_weight"] = w.astype(float)

    if "phase2" in out.columns:
        p2 = pd.to_numeric(out["phase2"], errors="coerce")
        if p2.isna().any() or not p2.isin([0, 1]).all():
            problems.append("phase2 must be 0/1")
        else:
            out["phase2"] = p2.astype(int)
    else:
        out["phase2"] = 1

    for c in markers:
        out[c] = pd.to_numeric(out[c], errors="coerce").astype(float)
    if out["phase2"].dtype.kind in "iu":
        measured = out.loc[out["phase2"] == 1, markers]
        if measured.isna().any().any():
            bad = measured.columns[measured.isna().any()].tolist()
            problems.append(f"phase-two rows have missing marker values in {bad}")

    if "time_to_event" in out.columns:
        t = pd.to_numeric(out["time_to_event"], errors="coerce")
        if t.isna().any() or (t < 0).any():
            problems.append("time_to_event must be non-negative with no missing values")
        else:
            out["time_to_event"] = t.astype(float)

    if problems:
        raise SchemaError("; ".join(problems))
    return out


def phase_two(df: pd.DataFrame) -> pd.DataFrame:
    """Rows with marker data: the ones every fit and metric uses."""
    if "phase2" not in df.columns:
        return df
    return df.loc[df["phase2"] == 1]


def design(df: pd.DataFrame, markers: list[str] | None = None):
    """Return (X, y, w) arrays from the phase-two rows of a validated frame."""
    d = phase_two(df)
    markers = markers or marker_columns(d)
    return (
        d[markers].to_numpy(float),
        d["endpoint"].to_numpy(int),
        d["sampling_weight"].to_numpy(float),
    )


# ----------------------------------------------------------------------------
# Marker mapping for E1
# ----------------------------------------------------------------------------

@dataclass(frozen=True)
class MarkerMap:
    """Declares which study column is which canonical marker, and on what scale.

    Read from a CSV with columns
    ``canonical, study_id, column, scale, multiplier, offset``.
    ``scale`` is one of ``linear``, ``log10`` or ``ln`` and says what scale the
    study column is on. Every column is converted to log10 and then put on the
    canonical scale by ``canonical = multiplier * value + offset``
    (multiplier defaults to 1, offset to 0). A study/marker pair absent from the
    file is a refusal, not a silent drop.
    """

    rows: pd.DataFrame = field(repr=False)

    @classmethod
    def from_csv(cls, path: str | Path) -> "MarkerMap":
        return cls.from_frame(pd.read_csv(path))

    @classmethod
    def from_frame(cls, rows: pd.DataFrame) -> "MarkerMap":
        need = {"canonical", "study_id", "column", "scale"}
        if not need <= set(rows.columns):
            raise SchemaError(f"marker map needs columns {sorted(need)}")
        rows = rows.copy()
        rows["study_id"] = rows["study_id"].astype(str)
        rows["multiplier"] = rows.get("multiplier", 1.0)
        rows["offset"] = rows.get("offset", 0.0)
        rows["multiplier"] = rows["multiplier"].fillna(1.0).astype(float)
        rows["offset"] = rows["offset"].fillna(0.0).astype(float)
        bad = set(rows["scale"]) - {"linear", "log10", "ln"}
        if bad:
            raise SchemaError(f"unknown scale(s) in marker map: {sorted(bad)}")
        if rows.duplicated(["canonical", "study_id"]).any():
            raise SchemaError("marker map has a canonical marker mapped twice in one study")
        return cls(rows)

    @property
    def canonical(self) -> list[str]:
        return sorted(self.rows["canonical"].unique())

    def apply(self, df: pd.DataFrame) -> pd.DataFrame:
        """Return ``df`` with its markers replaced by canonical ``marker_*`` columns."""
        out = df.drop(columns=marker_columns(df))
        canon = {c: np.full(len(df), np.nan) for c in self.canonical}
        for study, idx in df.groupby("study_id").groups.items():
            sub = self.rows[self.rows["study_id"] == str(study)]
            absent = set(self.canonical) - set(sub["canonical"])
            if absent:
                raise SchemaError(f"marker map has no entry for study {study!r}: {sorted(absent)}")
            pos = df.index.get_indexer(idx)
            for r in sub.itertuples():
                if r.column not in df.columns:
                    raise SchemaError(f"study {study!r} has no column {r.column!r}")
                v = df.loc[idx, r.column].to_numpy(float)
                if r.scale == "linear":
                    with np.errstate(divide="ignore", invalid="ignore"):
                        v = np.log10(v)
                elif r.scale == "ln":
                    v = v / np.log(10)
                canon[r.canonical][pos] = r.multiplier * v + r.offset
        for c, v in canon.items():
            name = c if c.startswith(MARKER_PREFIX) else MARKER_PREFIX + c
            out[name] = v
        return out


def require_alignment(source: pd.DataFrame, target: pd.DataFrame, marker_map: MarkerMap | None):
    """Refuse to compare across studies without an explicit marker map."""
    s, t = set(source["study_id"]), set(target["study_id"])
    if s != t and marker_map is None:
        raise SchemaError(
            f"source studies {sorted(s)} and target studies {sorted(t)} differ; "
            "E1 needs an explicit MarkerMap declaring how markers correspond"
        )
