import numpy as np
import pandas as pd
import pytest

from penalty import schema, transfer
from sim.generate import Gap, case_cohort, make_studies


def frame(n=20):
    rng = np.random.default_rng(0)
    return pd.DataFrame({
        "participant_id": [f"p{i}" for i in range(n)],
        "study_id": "A", "stratum": "s", "platform": "x",
        "marker_a": rng.normal(size=n), "endpoint": rng.integers(0, 2, n),
        "sampling_weight": 1.0,
    })


def test_valid_frame_passes_and_gets_phase2():
    out = schema.validate(frame())
    assert (out["phase2"] == 1).all()


def test_weight_column_is_required():
    with pytest.raises(schema.SchemaError, match="sampling_weight"):
        schema.validate(frame().drop(columns="sampling_weight"))


def test_default_weights_only_by_explicit_call_and_warns():
    with pytest.warns(UserWarning, match="simple random sample"):
        out = schema.with_default_weights(frame().drop(columns="sampling_weight"))
    assert (schema.validate(out)["sampling_weight"] == 1.0).all()


@pytest.mark.parametrize("bad", [0.0, -1.0, np.nan, np.inf])
def test_bad_weights_rejected(bad):
    df = frame()
    df.loc[3, "sampling_weight"] = bad
    with pytest.raises(schema.SchemaError, match="sampling_weight"):
        schema.validate(df)


def test_all_problems_reported_together():
    df = frame()
    df.loc[0, "endpoint"] = 2
    df.loc[1, "participant_id"] = "p0"
    with pytest.raises(schema.SchemaError) as e:
        schema.validate(df)
    assert "endpoint" in str(e.value) and "unique" in str(e.value)


def test_phase_one_rows_may_lack_markers_but_phase_two_may_not():
    df = frame()
    df["phase2"] = 1
    df.loc[0, ["phase2", "marker_a"]] = [0, np.nan]
    schema.validate(df)
    df.loc[1, "marker_a"] = np.nan
    with pytest.raises(schema.SchemaError, match="missing marker"):
        schema.validate(df)


def test_e1_refused_without_marker_map():
    s, t = make_studies(300, 300, Gap(), np.random.default_rng(1), same_study=False)
    with pytest.raises(schema.SchemaError, match="MarkerMap"):
        transfer.run(s, t)


def test_marker_map_puts_studies_on_one_scale():
    a = pd.DataFrame({"participant_id": ["a1", "a2"], "study_id": "A", "titre": [100.0, 1000.0]})
    b = pd.DataFrame({"participant_id": ["b1", "b2"], "study_id": "B", "ln_igg": [np.log(100.0), np.log(10.0)]})
    m = schema.MarkerMap.from_frame(pd.DataFrame({
        "canonical": ["igg", "igg"], "study_id": ["A", "B"],
        "column": ["titre", "ln_igg"], "scale": ["linear", "ln"],
        "multiplier": [1.0, 1.0], "offset": [0.0, 0.5],
    }))
    out = m.apply(pd.concat([a, b], ignore_index=True))
    np.testing.assert_allclose(out["marker_igg"], [2.0, 3.0, 2.5, 1.5])


def test_marker_map_refuses_a_study_it_does_not_cover():
    m = schema.MarkerMap.from_frame(pd.DataFrame({
        "canonical": ["igg"], "study_id": ["A"], "column": ["x"], "scale": ["log10"]}))
    df = pd.DataFrame({"participant_id": ["b1"], "study_id": "B", "x": [1.0]})
    with pytest.raises(schema.SchemaError, match="no entry for study 'B'"):
        m.apply(df)


def test_case_cohort_weights_are_inverse_selection_probability():
    s, _ = make_studies(5000, 10, Gap(), np.random.default_rng(2))
    cc = schema.validate(case_cohort(s, 0.2, np.random.default_rng(3)))
    assert (cc.loc[cc.endpoint == 1, "sampling_weight"] == 1).all()
    assert np.allclose(cc.loc[cc.endpoint == 0, "sampling_weight"], 5.0)
    # Weighted count recovers the cohort size within sampling error.
    assert abs(cc["sampling_weight"].sum() - len(s)) < 0.05 * len(s)
    assert (cc.endpoint == 1).sum() == s.endpoint.sum()
