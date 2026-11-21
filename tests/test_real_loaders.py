"""Loaders for released real datasets (AMP, PNG P. vivax), tested on synthetic
stand-ins shaped like each release. No real participant data is used (AC7)."""

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


def fake_pk(root, frame_root, seed=1):
    """Synthetic PT80-release files for a fake_amp cohort: all cases plus 30
    controls per trial, with uids colliding across trials as in the real one."""
    rng = np.random.default_rng(seed)
    surv = pd.read_csv(frame_root / "survival.csv")
    vax = surv[surv.rx_code != "C3"]
    cases = vax[vax.hiv1event == 1]
    ctl = vax[vax.hiv1event == 0].groupby("protocol").head(30)
    rows = []
    for r in pd.concat([cases, ctl]).itertuples():
        for t in range(1, 561):
            rows.append((t, float(np.exp(rng.normal(3, 0.5))), r.uid))
    root.mkdir(parents=True)
    pd.DataFrame(rows, columns=["time", "Cc", "uid"]).to_csv(root / "dat_est_daily_caseCtl.tab", sep="\t", index=False)
    study = cases["protocol"].str.replace(" ", "").str.replace("HVTN", "HVTN", regex=False)
    pd.DataFrame({"uid": cases.uid, "TAFD_MEDIAN": cases.hiv1survday - 10,
                  "study": study, "Cc": 1.0}).to_csv(root / "dat_case_EstInfTime_cc_nab.tab", sep="\t", index=False)
    in703 = pd.concat([cases, ctl])
    in703 = in703[in703.protocol == "HVTN 703/HPTN 081"]
    pd.DataFrame({"ID": range(len(in703)), "uid": in703.uid}).to_csv(root / "id_uid_703.tab", sep="\t", index=False)
    return len(cases), len(ctl)


def test_amp_pk_two_phase_weights_reproduce_phase_one(tmp_path):
    root, *_ = fake_amp(tmp_path / "core", n=900)
    n_cases, n_ctl = fake_pk(tmp_path / "pt80", root)
    frame, _ = loaders.amp(tmp_path, {"column": "trial"}, markers="pk")
    p2 = frame[frame.phase2 == 1]
    assert 0 < len(p2) <= n_cases + n_ctl
    # IPW weights restore the phase-one counts in every sampling stratum.
    assert np.isclose(p2["sampling_weight"].sum(), len(frame))
    assert {"marker_log10_conc_mean", "marker_log10_conc_p10"} <= set(frame.columns)


def fake_png(path, n=200, seed=2):
    """Synthetic file shaped like the PNG P. vivax release (DBP_data_fid)."""
    rng = np.random.default_rng(seed)
    rows = []
    villages = ["Ilaita 1", "Ilaita 4", "Sunuhu 1", "Kamakor"]
    for i in range(n):
        igg = {f"{a}_dil_igg_v1": 10 ** rng.normal(-3.5, 0.5) for a in loaders.PNG_PV_ANTIGENS}
        v, age = villages[i % 4], rng.uniform(1, 4)
        for wk in range(8):
            rows.append({"fid": f"h{i}", "fuwk": wk, "pv_500": int(rng.random() < 0.15),
                         "village": v, "age_yr": age, **igg})
    pd.DataFrame(rows).to_csv(path, index=False)
    return path


def test_png_pv_loader_and_e2(tmp_path):
    path = fake_png(tmp_path / "png.csv")
    frame, placebo = loaders.png_pv(path, {"column": "village_group"})
    assert set(frame["stratum"]) == {"Ilaita", "other"} and len(frame) == 200
    assert len([c for c in frame.columns if c.startswith("marker_")]) == 6
    assert placebo.empty and frame["time_to_event"].between(1, 7).all()
    raw = pd.read_csv(path)
    fu = raw[raw.fuwk >= 1]
    assert frame["endpoint"].sum() == fu[fu.pv_500 > 0]["fid"].nunique()
    cfg = {"dataset": "png_pv", "split": {"column": "village_group"}, "source": "Ilaita",
           "target": "other", "models": ["ridge_logistic", "ridge_cox"], "baseline": "cv",
           "folds": 3, "seed": 1, "placebo_offset": False, "bootstrap_reps": 3}
    res = e2_within_study.run(cfg, path, reps=3)
    assert {"auc", "c_index"} <= set(res["penalties"]["metric"])
