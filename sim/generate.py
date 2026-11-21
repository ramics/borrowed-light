"""Dial-driven simulator of a source and a target study, in the tidy schema.

The true gap between source and target is set by independent dials, each zero
when the two are identical, so the right answer is known and the estimator can
be checked against it. Synthetic data only; no real participant data is
involved here or anywhere in the repo.

Dials (composable, any combination):

* ``mechanism``  radians; rotates the true coefficient vector in the plane of
  markers 0 and 1. "The markers mean something different over there."
* ``assay``      multiplicative scale shift on markers 0-2 in the target.
  "Measured differently."
* ``population`` shift in the target's baseline log-odds. "Different people."
* ``noise``      SD of extra measurement error on every target marker.

:func:`case_cohort` turns a full cohort into a two-phase release with the
sampling weights a real correlates dataset would carry.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np
import pandas as pd

P = 8
BETA = np.array([1.1, -0.8, 0.6, 0.0, 0.45, 0.0, -0.3, 0.2])
INTERCEPT = -1.2


@dataclass(frozen=True)
class Gap:
    mechanism: float = 0.0
    assay: float = 0.0
    population: float = 0.0
    noise: float = 0.0

    def as_dict(self):
        return asdict(self)


def rotation(theta: float, p: int = P) -> np.ndarray:
    r = np.eye(p)
    r[0, 0] = r[1, 1] = np.cos(theta)
    r[0, 1] = -np.sin(theta)
    r[1, 0] = np.sin(theta)
    return r


def _cohort(n, beta, intercept, gap: Gap, rng, study_id, stratum, prefix):
    cov = 0.3 ** np.abs(np.subtract.outer(np.arange(P), np.arange(P)))
    x = rng.multivariate_normal(np.zeros(P), cov, size=n)
    eta = x @ beta + intercept + gap.population
    y = rng.binomial(1, 1 / (1 + np.exp(-eta)))
    # Measurement artefacts act on what is recorded, not on the biology that
    # generated the outcome.
    if gap.noise:
        x = x + rng.normal(0, gap.noise, size=x.shape)
    if gap.assay:
        x[:, :3] *= 1.0 + gap.assay
    df = pd.DataFrame(x, columns=[f"marker_{j}" for j in range(P)])
    df.insert(0, "participant_id", [f"{prefix}{i:06d}" for i in range(n)])
    df.insert(1, "study_id", study_id)
    df.insert(2, "stratum", stratum)
    df.insert(3, "platform", "synthetic")
    df["endpoint"] = y
    df["sampling_weight"] = 1.0
    return df


def make_studies(n_source: int, n_target: int, gap: Gap = Gap(), rng=None,
                 same_study: bool = True):
    """Return (source, target) tidy frames separated by ``gap``.

    With ``same_study`` both frames share a study_id and differ by stratum
    (the E2 shape). Otherwise they get distinct study_ids (the E1 shape).
    """
    rng = rng if rng is not None else np.random.default_rng()
    beta_t = rotation(gap.mechanism) @ BETA
    src = _cohort(n_source, BETA, INTERCEPT, Gap(), rng, "SIM-A", "source", "S")
    tgt = _cohort(n_target, beta_t, INTERCEPT, gap, rng,
                  "SIM-A" if same_study else "SIM-B", "target", "T")
    # The noise/assay artefacts above must not apply to the source.
    return src, tgt


def case_cohort(df: pd.DataFrame, subcohort_frac: float, rng=None,
                strata: str | None = None, keep_phase_one: bool = False) -> pd.DataFrame:
    """Two-phase subsample: a random subcohort plus every case.

    ``subcohort_frac`` may be a float or, with ``strata``, a dict from stratum
    value to fraction. Weights are 1 / P(selected): 1 for cases, 1/f for
    non-cases. With ``keep_phase_one`` the unselected rows stay in with
    ``phase2 = 0`` and their markers blanked, the form a full-cohort release
    takes.
    """
    rng = rng if rng is not None else np.random.default_rng()
    if strata is None:
        f = np.full(len(df), float(subcohort_frac))
    else:
        f = df[strata].map(subcohort_frac).to_numpy(float)
        if np.isnan(f).any():
            raise ValueError("subcohort_frac has no entry for some stratum")
    in_sub = rng.random(len(df)) < f
    case = df["endpoint"].to_numpy() == 1
    selected = in_sub | case
    prob = np.where(case, 1.0, f)
    out = df.copy()
    out["sampling_weight"] = 1.0 / prob
    out["sampling_stratum"] = (df[strata].astype(str) if strata else "all")
    out["sampling_stratum"] = out["sampling_stratum"] + np.where(case, ":case", ":noncase")
    if keep_phase_one:
        out["phase2"] = selected.astype(int)
        markers = [c for c in out.columns if c.startswith("marker_")]
        out.loc[~selected, markers] = np.nan
        return out
    out = out.loc[selected].copy()
    out["phase2"] = 1
    return out
