"""E2, within-study population transfer.

    python -m experiments.e2_within_study experiments/configs/e2_hvtn505_age.yaml \
        --data /path/outside/repo/<file or directory> [--out reports/<run>] [--reps N]

Fits each configured model on the source stratum, measures it at home and on
the target stratum, and bootstraps every penalty. Writes to the output
directory: ``manifest.json``, ``penalties.csv`` (one row per model x metric,
point estimate and interval), ``bootstrap_draws.csv`` and ``report.md`` (see :mod:`experiments.report`).
Only aggregate numbers are written; no participant-level rows leave the
data's location.
"""

from __future__ import annotations

import argparse
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from penalty import resample, schema, transfer

from . import manifest, report
from .loaders import LOADERS

REPO = Path(__file__).resolve().parents[1]


def _refuse_in_repo(path: Path):
    try:
        path.resolve().relative_to(REPO)
    except ValueError:
        return
    raise SystemExit(f"refusing to read data from inside the repository ({path}); AC7")


def _placebo(counts: pd.DataFrame, source: str, target: str) -> transfer.PlaceboArms:
    c = counts.set_index("stratum")
    return transfer.PlaceboArms(int(c.loc[source, "events"]), int(c.loc[source, "n"]),
                                int(c.loc[target, "events"]), int(c.loc[target, "n"]))


def describe(df: pd.DataFrame) -> dict:
    p2 = schema.phase_two(df)
    return {"n_phase1": int(len(df)), "n_phase2": int(len(p2)),
            "cases": int(df["endpoint"].sum()), "cases_phase2": int(p2["endpoint"].sum()),
            "weighted_n": float(p2["sampling_weight"].sum())}


def run(config: dict, data: Path, reps: int | None = None) -> dict:
    _refuse_in_repo(data)
    frame, placebo_counts = LOADERS[config["dataset"]](data, config["split"], **config.get("loader_options", {}))
    src = frame[frame["stratum"] == config["source"]]
    tgt = frame[frame["stratum"] == config["target"]]
    placebo = _placebo(placebo_counts, config["source"], config["target"]) if config.get("placebo_offset") else None
    reps = reps or config.get("bootstrap_reps", 1000)
    rows, draws, timing, fits = [], [], {}, {}
    for model in config["models"]:
        t0 = time.time()
        kw = dict(model=model, baseline=config.get("baseline", "holdout"), folds=config.get("folds", 5),
                  placebo=placebo, risk_threshold=config.get("risk_threshold"))
        point = transfer.run(src, tgt, seed=config["seed"], **kw)
        ivs = resample.bootstrap(src, tgt, reps=reps, alpha=config.get("alpha", 0.05),
                                 seed=config["seed"], **kw)
        timing[model] = round(time.time() - t0, 1)
        tuning = point.model.tuning
        fits[model] = {
            "markers": list(point.model.markers),
            "C": tuning.get("C"),
            "C_at_grid_floor": ("C" in tuning and tuning["C"] == min(tuning["c_grid"])),
            "chosen_marker": tuning.get("chosen"),
            "coef_per_sd": [float(c) for c in point.model.coef],
            "n": point.n,
        }
        for metric, iv in ivs.items():
            base_metric = metric.replace("_exposure_adjusted", "")
            rows.append({
                "model": model, "metric": metric,
                "baseline": point.baseline.get(base_metric),
                "target": (point.target_adjusted if metric.endswith("_adjusted") else point.target).get(base_metric),
                "penalty": iv.estimate, "lo": iv.lo, "hi": iv.hi,
                "covers_zero": iv.covers(0.0), "valid_reps": int(len(iv.reps)),
                "markers": ",".join(point.model.markers),
            })
            draws.append(pd.DataFrame({"model": model, "metric": metric, "draw": iv.reps}))
    return {
        "penalties": pd.DataFrame(rows),
        "draws": pd.concat(draws, ignore_index=True),
        "source": describe(src), "target": describe(tgt),
        "placebo": None if placebo is None else placebo.__dict__ | {"offset_logit": placebo.offset("logit")},
        "reps": reps, "seconds": timing, "fits": fits,
    }


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("config", type=Path)
    ap.add_argument("--data", type=Path, default=os.environ.get("BORROWED_LIGHT_DATA"))
    ap.add_argument("--out", type=Path)
    ap.add_argument("--reps", type=int)
    a = ap.parse_args(argv)
    if a.data is None:
        raise SystemExit("--data (or BORROWED_LIGHT_DATA) is required")
    config = yaml.safe_load(a.config.read_text())
    if a.reps:
        config["bootstrap_reps"] = a.reps
    out = a.out or REPO / "reports" / f"{a.config.stem}-{time.strftime('%Y%m%dT%H%M%S')}"
    m = manifest.build(config, a.data, sys.argv if argv is None else argv)
    res = run(config, a.data, config.get("bootstrap_reps"))
    manifest.write(m | {"seconds": res["seconds"], "fits": res["fits"],
                        "source": res["source"], "target": res["target"], "placebo": res["placebo"]}, out)
    res["penalties"].to_csv(out / "penalties.csv", index=False)
    res["draws"].to_csv(out / "bootstrap_draws.csv", index=False)
    (out / "report.md").write_text(report.render(config, res))
    print((out / "report.md").read_text())
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
