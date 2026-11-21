"""E2 plumbing on a synthetic stand-in shaped like the HVTN 505 release.
No real participant data is used (AC7)."""

import numpy as np
import pandas as pd
import pytest

from experiments import e2_within_study, loaders


def fake_505(path, n=1200, seed=0):
    rng = np.random.default_rng(seed)
    age = rng.integers(18, 51, n)
    trt = rng.integers(0, 2, n)
    env = rng.normal(0.8, 0.3, n)
    y = rng.binomial(1, 1 / (1 + np.exp(-(-2.0 - 1.5 * (env - 0.8) * trt))))
    cc = (y == 1) | (rng.random(n) < 0.2)
    df = pd.DataFrame({
        "pub_id": [f"x-{i}" for i in range(n)], "trt": trt, "HIVwk28preunbl": y,
        "HIVwk28preunblfu": np.where(y == 1, rng.uniform(1, 578, n), 578.0),
        "age": age, "BMI": rng.normal(26, 5, n), "bhvrisk": rng.random(n),
        "casecontrol": cc.astype(int), "wt": np.where(cc, np.where(y == 1, 1.0, 5.0), np.nan),
        "IgG_env": np.where(cc, env, np.nan), "IgG_V2": np.where(cc, rng.normal(1, 0.4, n), np.nan),
        "IgG_V3": np.where(cc, rng.normal(1.5, 0.5, n), np.nan),
    })
    df.to_csv(path, index=False)
    return path


SPLIT = {"column": "age", "cut": 30, "labels": ["young", "old"]}


def test_loader_produces_valid_two_phase_frame(tmp_path):
    frame, placebo = loaders.hvtn505(fake_505(tmp_path / "f.csv"), SPLIT)
    assert set(frame["stratum"]) == {"young", "old"}
    assert set(frame["phase2"]) == {0, 1}
    assert (frame.loc[frame.phase2 == 1, "sampling_weight"] >= 1).all()
    assert set(placebo["stratum"]) == {"young", "old"}


def test_e2_runs_end_to_end_and_writes_manifest(tmp_path):
    data = fake_505(tmp_path / "f.csv")
    cfg = {"dataset": "hvtn505", "split": SPLIT, "source": "young", "target": "old",
           "models": ["ridge_logistic", "single_marker", "ridge_cox"], "baseline": "cv",
           "folds": 5, "seed": 1, "placebo_offset": True, "bootstrap_reps": 5}
    path = tmp_path / "cfg.yaml"
    import yaml
    path.write_text(yaml.safe_dump(cfg))
    out = tmp_path / "out"
    e2_within_study.main([str(path), "--data", str(data), "--out", str(out)])
    pen = pd.read_csv(out / "penalties.csv")
    assert {"auc", "c_index", "cal_intercept_exposure_adjusted"} <= set(pen["metric"])
    rep = (out / "report.md").read_text()
    headings = [l for l in rep.splitlines() if l.startswith("## ")]
    assert headings.index("## Calibration (headline)") < headings.index("## Discrimination")
    assert "## Decisions" in headings and "## Read this first" in headings
    assert "flip_fixed" in set(pen["metric"])
    m = (out / "manifest.json").read_text()
    assert "sha256" in m and str(data) not in m  # fingerprint, not location


def test_refuses_data_inside_repo():
    with pytest.raises(SystemExit, match="AC7"):
        e2_within_study._refuse_in_repo(e2_within_study.REPO / "data" / "x.csv")


def test_report_handles_a_categorical_split():
    from experiments import report
    res = {"penalties": pd.DataFrame(columns=["model", "metric", "baseline", "target", "penalty",
                                              "lo", "hi", "covers_zero", "valid_reps"]),
           "source": dict(n_phase1=10, n_phase2=10, cases=5, cases_phase2=5),
           "target": dict(n_phase1=10, n_phase2=10, cases=5, cases_phase2=5),
           "placebo": None, "reps": 10, "fits": {}}
    md = report.render({"dataset": "amp", "split": {"column": "trial"}, "source": "a", "target": "b"}, res)
    assert "`trial`;" in md and " at " not in md.split("Split:")[1].split("\n")[0]
