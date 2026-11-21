"""The report template every experiment writes as ``report.md``.

Order is deliberate. Calibration and decision flips come first because
admissibility is a threshold decision and a model can rank perfectly while
sitting on the wrong scale (SPEC.md). Discrimination follows, then overall
accuracy, then the single-marker comparator. Sample-size warnings are placed
before any number they qualify. The closing caveats are fixed text: what
transfers from this work is the estimator, not any one pair's numbers.
"""

from __future__ import annotations

import math

import pandas as pd

# Riley et al. (2021, Stat Med 40:4230): precise calibration-slope estimates
# need on the order of 100+ events in the validation set; with very few, the
# honest answer is "cannot tell".
EVENTS_FOR_SLOPE = 100
EVENTS_FOR_ANYTHING = 10

SECTIONS = [
    ("Calibration (headline)", ["cal_intercept", "cal_slope"],
     "Calibration-in-the-large (`cal_intercept`, 0 is perfect, positive means the borrowed model "
     "under-predicts target risk) and calibration slope (`cal_slope`, 1 is perfect). Rows show "
     "the signed shift target − baseline."),
    ("Decisions", ["flip_fixed", "flip_optimal"],
     "Share of the target population whose call changes. `flip_fixed`: a risk threshold fixed in "
     "the source, borrowed model versus the same model recalibrated to the target. "
     "`flip_optimal`: the source-optimal cut versus the target-optimal cut. Penalty is "
     "target − baseline; positive means more decisions flip across the gap."),
    ("Discrimination", ["auc", "c_index"],
     "Does the model still rank people? Penalty is baseline − target; positive means worse in "
     "the target."),
    ("Overall accuracy", ["brier"],
     "Brier score. Penalty is target − baseline; positive means worse in the target."),
]


def _fmt(x, nd=3, sign=False):
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return "–"
    return f"{x:+.{nd}f}" if sign else f"{x:.{nd}f}"


def _warnings(res: dict) -> list[str]:
    out = []
    for side in ("source", "target"):
        ev = res[side]["cases_phase2"]
        if ev < EVENTS_FOR_ANYTHING:
            out.append(f"**{side.capitalize()} has {ev} events with markers.** Below "
                       f"{EVENTS_FOR_ANYTHING}, no penalty here is interpretable; read every "
                       "interval as \"cannot tell\".")
        elif ev < EVENTS_FOR_SLOPE:
            out.append(f"{side.capitalize()} has {ev} events with markers, under the ~"
                       f"{EVENTS_FOR_SLOPE} needed to estimate a calibration slope usefully "
                       "(Riley et al. 2021). Treat slope rows as not estimable.")
    for model, f in res["fits"].items():
        if f.get("C_at_grid_floor"):
            out.append(f"`{model}`: cross-validation in the source chose the strongest penalty on "
                       "the grid, so the model is close to intercept-only. The source shows no "
                       "marker signal it can separate from noise; its discrimination and slope "
                       "rows describe noise.")
    failed = res["penalties"].query("valid_reps < @res['reps']")
    if len(failed):
        worst = int(res["reps"] - failed["valid_reps"].min())
        out.append(f"Up to {worst} of {res['reps']} bootstrap replicates had too few cases to fit "
                   "and were dropped.")
    return out


# A constant exposure offset cannot change ranking or slope, so adjusted rows
# are shown only where they can differ.
ADJUSTABLE = {"cal_intercept", "flip_fixed", "flip_optimal", "brier"}


def _table(pen: pd.DataFrame, metrics: list[str], adjusted: bool) -> list[str]:
    adj = [m + "_exposure_adjusted" for m in metrics if m in ADJUSTABLE] if adjusted else []
    rows = pen[pen["metric"].isin(metrics + adj)]
    if rows.empty:
        return []
    lines = ["| model | metric | baseline | target | penalty | 95% interval | covers 0 |",
             "| --- | --- | --- | --- | --- | --- | --- |"]
    for r in rows.itertuples():
        label = r.metric.replace("_exposure_adjusted", " (exposure-adjusted)")
        lines.append(f"| {r.model} | {label} | {_fmt(r.baseline)} | {_fmt(r.target)} | "
                     f"{_fmt(r.penalty, sign=True)} | [{_fmt(r.lo, sign=True)}, {_fmt(r.hi, sign=True)}] | "
                     f"{'yes' if r.covers_zero else '**no**'} |")
    return lines


def _comparator(pen: pd.DataFrame) -> list[str]:
    rank = pen[pen["metric"].isin(["auc", "c_index"])].set_index("model")["baseline"]
    if "single_marker" not in rank or "ridge_logistic" not in rank:
        return []
    multi, single = rank["ridge_logistic"], rank["single_marker"]
    verdict = ("The multivariate model does **not** beat the best single marker at home."
               if multi <= single else
               "The multivariate model beats the best single marker at home.")
    return [f"Home AUC: ridge logistic {_fmt(multi)}, best single marker {_fmt(single)}. {verdict} "
            "(Point estimates; see the discrimination table for intervals.)"]


def render(config: dict, res: dict) -> str:
    pen = res["penalties"]
    adjusted = res.get("placebo") is not None
    split = config.get("split", {})
    L = [f"# {config.get('experiment', 'transfer')}: {config['dataset']}, "
         f"{config['source']} → {config['target']}", "",
         "This report measures how much a correlate-of-risk model loses when it is borrowed from "
         "the source population to the target. It makes no claim that any vaccine works, and no "
         "claim about any pathogen.", "",
         "## Set-up", "",
         f"- Split: `{split.get('column')}` at {split.get('cut')}; baseline measured by "
         f"`{config.get('baseline', 'holdout')}`.",
         f"- Uncertainty: {res['reps']} participant-level bootstrap replicates, "
         f"{int(round(100 * (1 - config.get('alpha', 0.05))))}% percentile intervals.",
         "", "| | source | target |", "| --- | --- | --- |"]
    for k, label in (("n_phase1", "participants (phase one)"), ("n_phase2", "with markers (phase two)"),
                     ("cases", "events"), ("cases_phase2", "events with markers")):
        L.append(f"| {label} | {res['source'][k]} | {res['target'][k]} |")
    if adjusted:
        p = res["placebo"]
        L += ["", f"Placebo arms: {p['source_events']}/{p['source_n']} events (source), "
              f"{p['target_events']}/{p['target_n']} (target). Exposure-adjusted rows shift the "
              f"borrowed model by the placebo-incidence difference ({p['offset_logit']:+.2f} on the "
              "logit scale), so what remains is the gap a different epidemic does not explain."]
    warns = _warnings(res)
    if warns:
        L += ["", "## Read this first", ""] + [f"- {w}" for w in warns]
    for title, metrics, blurb in SECTIONS:
        tab = _table(pen, metrics, adjusted)
        if tab:
            L += ["", f"## {title}", "", blurb, ""] + tab
    comp = _comparator(pen)
    if comp:
        L += ["", "## Single-marker comparator", ""] + comp
    L += ["", "## Models", "", "| model | markers | penalty strength |", "| --- | --- | --- |"]
    for model, f in res["fits"].items():
        strength = f"C = {f['C']:.3g}" if f.get("C") is not None else f"chose {f.get('chosen_marker')}"
        L.append(f"| {model} | {', '.join(f['markers'])} | {strength} |")
    L += ["", "## Caveats", "",
          "- What transfers from this work is the estimator, not these numbers. A penalty measured "
          "on one pair of populations says little about another pair.",
          "- Intervals are participant-level bootstrap percentiles, which are known to be "
          "unreliable at small n.",
          "- Penalty strength and every threshold were fixed in the source before the target was "
          "evaluated.", "",
          "## Reproduce", "",
          "Everything needed is in `manifest.json` next to this file: config, commit, seed, "
          "package versions and the data's SHA-256.", ""]
    return "\n".join(L)
