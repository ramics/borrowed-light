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


AMP_TRIALS = {"HVTN 704/HPTN 085": "HVTN704", "HVTN 703/HPTN 081": "HVTN703"}
AMP_DOSE = {"T1": 10.0, "T2": 30.0}
AMP_INFUSION_EVERY = 56  # days; infusions at weeks 0, 8, ..., 72
AMP_INFUSIONS = 10


def amp(path: str | Path, split: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    """AMP trials (HVTN 704/HPTN 085 and HVTN 703/HPTN 081; VRC01 infusions),
    primary-results release on the HVTN Harvard Dataverse (doi:10.7910/DVN/GDESJI).

    ``path`` is the release's ``datasets/`` directory, holding
    ``survival.csv``, ``subject_master.csv`` and ``retention.csv``.

    Endpoint: primary HIV-1 infection (``hiv1event``) with follow-up
    ``hiv1survday``, as in the primary analysis. "Vaccine" arm = pooled VRC01
    arms (T1 10 mg/kg, T2 30 mg/kg); placebo = C3. Full cohort, so every
    weight is 1 and every row is phase two.

    Markers. This release has no assay readout for the whole cohort (serum
    concentrations exist only for a 47-person pilot), so the markers are
    dose-and-adherence proxies, *not* assays:

    - ``marker_dose_mgkg``: assigned dose per infusion.
    - ``marker_adherence``: infusions received before the endpoint or
      censoring day, over infusions scheduled by then.
    - ``marker_delivered_mgkg``: dose x adherence.

    Model-predicted serum concentrations, the marker the AMP correlates work
    uses, are in the separate releases doi:10.7910/DVN/50VH2C and
    doi:10.7910/DVN/WEJPUF.

    ``split`` = {"column": "trial"} (HVTN704 vs HVTN703) or
    {"column": "age" | "sex" | "country", ...}; numeric columns take ``cut``.
    """
    root = Path(path)
    surv = pd.read_csv(root / "survival.csv")
    subj = pd.read_csv(root / "subject_master.csv")
    ret = pd.read_csv(root / "retention.csv")
    key = ["protocol", "uid"]
    d = surv.merge(subj[key + ["age", "sex", "country"]], on=key, how="left", validate="one_to_one")
    d["trial"] = d["protocol"].map(AMP_TRIALS)
    if d["trial"].isna().any():
        raise ValueError(f"unknown protocol(s): {sorted(set(d.loc[d.trial.isna(), 'protocol']))}")
    d = d[d["hiv1survday"] > 0]

    inf = ret[(ret["recv"] == "Y") & ret["visitlabel"].str.contains("Infusion", na=False)
              & ~ret["visitlabel"].str.contains("post", case=False, na=False)]
    inf = inf.merge(d[key + ["hiv1survday"]], on=key)
    received = (inf[inf["visitday"] < inf["hiv1survday"]].groupby(key).size()
                .rename("received").reset_index())
    d = d.merge(received, on=key, how="left").fillna({"received": 0})
    scheduled = np.minimum(AMP_INFUSIONS, np.floor((d["hiv1survday"] - 1) / AMP_INFUSION_EVERY) + 1)
    d["adherence"] = np.clip(d["received"] / scheduled, 0.0, 1.0)

    col = split["column"]
    if "cut" in split:
        d["stratum"] = stratify(d[col], split)
    else:
        d["stratum"] = d[col].astype(str)
    vax, pla = d[d["rx_code"].isin(AMP_DOSE)], d[d["rx_code"] == "C3"]
    dose = vax["rx_code"].map(AMP_DOSE).astype(float)
    df = pd.DataFrame({
        "participant_id": vax["trial"] + "-" + vax["uid"].astype(str),
        # One study_id: the two trials ran as one programme with the same
        # product, schedule, labs and SAP, so a trial split is treated as an
        # E2 population split rather than E1 (which would need a MarkerMap).
        "study_id": "AMP",
        "stratum": vax["stratum"],
        "platform": "VRC01 bnAb infusion",
        "marker_dose_mgkg": dose,
        "marker_adherence": vax["adherence"].astype(float),
        "marker_delivered_mgkg": dose * vax["adherence"].astype(float),
        "endpoint": vax["hiv1event"].astype(int),
        "time_to_event": vax["hiv1survday"].astype(float),
        "sampling_weight": 1.0,
        "phase2": 1,
    }).reset_index(drop=True)
    placebo = (pla.groupby("stratum")["hiv1event"].agg(events="sum", n="count")
               .astype(int).reset_index())
    return schema.validate(df), placebo


LOADERS = {"hvtn505": hvtn505, "amp": amp}
