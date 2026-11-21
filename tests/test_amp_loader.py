"""AMP loader on a synthetic stand-in shaped like the GDESJI release.
No real participant data is used (AC7)."""

import numpy as np
import pandas as pd

from experiments import e2_within_study, loaders, manifest

PROTOCOLS = list(loaders.AMP_TRIALS)


def fake_amp(root, n=600, seed=0):
    rng = np.random.default_rng(seed)
    root.mkdir(parents=True, exist_ok=True)
    proto = rng.choice(PROTOCOLS, n)
    rx = rng.choice(["T1", "T2", "C3"], n)
    y = rng.binomial(1, 0.08, n)
    t = np.where(y == 1, rng.integers(30, 560, n), 590)
    pd.DataFrame({"protocol": proto, "uid": np.arange(n), "rx_code": rx, "hiv1event": y,
                  "hiv1survday": t}).to_csv(root / "survival.csv", index=False)
    pd.DataFrame({"protocol": proto, "uid": np.arange(n), "age": rng.integers(18, 50, n),
                  "sex": rng.choice(["Male", "Female"], n), "country": "X"}).to_csv(
        root / "subject_master.csv", index=False)
    rows = []
    for i in range(n):
        for k in range(10):
            rows.append((proto[i], i, f"Week {8 * k}/Infusion {k + 1}", "Y",
                         "Y" if rng.random() < 0.9 else "N", 56 * k))
            rows.append((proto[i], i, f"Week {8 * k + 4}/28 Days post Infusion {k + 1}", "Y", "N", 56 * k + 28))
    pd.DataFrame(rows, columns=["protocol", "uid", "visitlabel", "attn", "recv", "visitday"]).to_csv(
        root / "retention.csv", index=False)
    return root, pd.Series(y, index=np.arange(n)), pd.Series(rx), pd.Series(proto)


def test_amp_loader_arms_markers_and_placebo(tmp_path):
    root, y, rx, proto = fake_amp(tmp_path / "datasets")
    frame, placebo = loaders.amp(root, {"column": "trial"})
    assert set(frame["stratum"]) == {"HVTN704", "HVTN703"}
    assert len(frame) == int((rx != "C3").sum())
    assert set(frame["marker_dose_mgkg"]) == {10.0, 30.0}
    assert frame["marker_adherence"].between(0, 1).all()
    assert (frame["sampling_weight"] == 1).all()
    pla = placebo.set_index("stratum")
    for trial, label in loaders.AMP_TRIALS.items():
        m = (rx == "C3") & (proto == trial)
        assert pla.loc[label, "n"] == m.sum() and pla.loc[label, "events"] == y[m].sum()


def test_amp_e2_end_to_end_and_directory_digest(tmp_path):
    root, *_ = fake_amp(tmp_path / "datasets", n=900)
    cfg = {"dataset": "amp", "split": {"column": "trial"}, "source": "HVTN704", "target": "HVTN703",
           "models": ["ridge_logistic", "ridge_cox"], "baseline": "cv", "folds": 3, "seed": 1,
           "placebo_offset": True, "bootstrap_reps": 3}
    res = e2_within_study.run(cfg, root, reps=3)
    assert {"auc", "c_index"} <= set(res["penalties"]["metric"])
    d1 = manifest.file_digest(root)
    (root / "survival.csv").write_text((root / "survival.csv").read_text() + "\n")
    assert manifest.file_digest(root) != d1
