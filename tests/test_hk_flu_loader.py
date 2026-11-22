"""HK serosurvey loader and clustered resampling, on synthetic stand-ins
shaped like the Zenodo release. No real participant data is used (AC7)."""

import numpy as np
import pandas as pd
import pytest

from experiments import loaders
from penalty import resample, transfer


def fake_hk(root, households=120, seed=0):
    rng = np.random.default_rng(seed)
    titres = np.array([5, 10, 20, 40, 80, 160, 320, 640])
    for k in (1, 2, 3, 4):
        rows = []
        for h in range(households):
            for m in range(rng.integers(1, 4)):
                if rng.random() < 0.25:
                    continue  # not in this epidemic's pair of sera
                pre = rng.choice(titres, p=[.6, .05, .1, .08, .07, .05, .03, .02])
                inf = rng.random() < 0.25 / (1 + np.log2(pre / 5))
                post = min(pre * (8 if inf else rng.choice([1, 2])), 1280)
                rows.append({"subjectID": f"K{100000 + h}-{m}", "age.y": int(rng.integers(3, 80)),
                             "sex": 1, "blood.date.ph1": "1/1/2010", "h3.ph1": pre,
                             "blood.date.ph2": "1/1/2011", "h3.ph2": post,
                             "vac": rng.choice([0, 1, np.nan], p=[.75, .2, .05])})
        pd.DataFrame(rows).to_csv(root / f"titer.epi{k}.csv", index=False)
    return root


def test_loader_frame(tmp_path):
    frame, placebo = loaders.hk_flu(fake_hk(tmp_path), {"column": "strain"})
    raw = pd.concat([pd.read_csv(tmp_path / f"titer.epi{k}.csv") for k in (1, 2, 3, 4)])
    assert len(frame) == int((raw["vac"].fillna(0) == 0).sum())  # vaccinees dropped, NA kept
    assert set(frame["stratum"]) == {"Perth", "Victoria"}
    assert frame["participant_id"].is_unique
    assert (frame["marker_log2_hai_pre"] >= 0).all()
    assert frame["cluster_id"].str.fullmatch(r"K\d+").all()  # household
    assert placebo.empty
    by_person = loaders.hk_flu(tmp_path, {"column": "strain"}, cluster="person")[0]
    assert by_person["cluster_id"].str.contains("-").all()


def test_folds_and_holdout_keep_clusters_together(tmp_path):
    frame, _ = loaders.hk_flu(fake_hk(tmp_path), {"column": "strain"})
    src = frame[frame.stratum == "Perth"]
    fold = transfer.assign_folds(src, 5, np.random.default_rng(1))
    assert (pd.Series(fold).groupby(src["cluster_id"].to_numpy()).nunique() == 1).all()
    train, held = transfer.split_source(src, 0.3, np.random.default_rng(2))
    assert not set(train["cluster_id"]) & set(held["cluster_id"])
    assert held["endpoint"].sum() > 0


def test_joint_cluster_bootstrap_draws_one_set_of_clusters(tmp_path):
    frame, _ = loaders.hk_flu(fake_hk(tmp_path), {"column": "strain"})
    src, tgt = frame[frame.stratum == "Perth"], frame[frame.stratum == "Victoria"]
    s, t = resample._resample_clusters([src, tgt], np.random.default_rng(3))
    assert s["participant_id"].is_unique and t["participant_id"].is_unique
    # a household drawn k times contributes its rows k times on each side it appears in
    base = lambda d: d["cluster_id"].str.split("#").str[0]  # noqa: E731
    shared = set(src.cluster_id) & set(tgt.cluster_id)
    copies_s = s.assign(b=base(s)).groupby("b")["cluster_id"].nunique()
    copies_t = t.assign(b=base(t)).groupby("b")["cluster_id"].nunique()
    both = [c for c in shared if c in copies_s.index or c in copies_t.index]
    assert both and all(copies_s.get(c, 0) == copies_t.get(c, 0) for c in both)


def test_clustered_bootstrap_runs(tmp_path):
    frame, _ = loaders.hk_flu(fake_hk(tmp_path), {"column": "strain"})
    src, tgt = frame[frame.stratum == "Perth"], frame[frame.stratum == "Victoria"]
    ivs = resample.bootstrap(src, tgt, reps=5, seed=0, model="ridge_logistic", baseline="cv")
    assert len(ivs["auc"].reps) == 5


def test_bad_cluster_option(tmp_path):
    with pytest.raises(ValueError, match="cluster"):
        loaders.hk_flu(fake_hk(tmp_path), {"column": "strain"}, cluster="village")
