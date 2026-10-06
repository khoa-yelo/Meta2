"""checkgm.healthscore: the arithmetic of the shipped health score, and the two refusals that matter more than it —
an assembly assessment, and a model file carrying the species-presence features whose top coefficient is a batch
artifact. The assessments are written here by hand in the shape cli.write_report writes, so nothing external is needed."""
import json
import os
import warnings

import numpy as np
import pandas as pd
import pytest

from checkgm import healthscore as hs

GRID_COLS = ["analysis_id", "feature_id", "band", "value", "percentile", "call", "call_fdr", "p_two_sided", "fdr_q", "layer"]
SUMMARY_COLS = ["analysis_id", "layer", "n_detected", "n_assessed", "n_low", "n_high", "n_low_raw", "n_high_raw",
                "n_missing", "n_not_assessable", "frac_outside_raw", "layer_status", "quality_band", "bundle_id"]


def write_assessment(d, rows, summaries):
    """rows: (sample, layer, feature_id, percentile, call); summaries: (sample, layer, frac_outside_raw, n_missing)."""
    os.makedirs(d, exist_ok=True)
    sc = pd.DataFrame([dict(analysis_id=s, feature_id=str(f), band="medium", value=0.0, percentile=p, call=c,
                            call_fdr=c, p_two_sided=np.nan, fdr_q=np.nan, layer=l) for s, l, f, p, c in rows],
                      columns=GRID_COLS)
    sc.to_parquet(os.path.join(d, "scores.parquet"), index=False)
    sm = pd.DataFrame([dict(analysis_id=s, layer=l, n_detected=1, n_assessed=1, n_low=0, n_high=0, n_low_raw=0,
                            n_high_raw=0, n_missing=nm, n_not_assessable=0, frac_outside_raw=fo, layer_status="scored",
                            quality_band="medium", bundle_id=b) for s, l, fo, nm, b in summaries], columns=SUMMARY_COLS)
    sm.to_csv(os.path.join(d, "summaries.tsv"), sep="\t", index=False)
    pd.DataFrame(columns=["analysis_id", "reason_code", "detail"]).to_csv(os.path.join(d, "rejections.tsv"), sep="\t", index=False)
    return d


def toy_model(p, features=None, **over):
    m = {"model": "toy", "feature_set": "pct", "intercept": -0.5,
         "features": features if features is not None else [
             {"name": "species:820", "coef": 0.4, "mean": 50.0, "scale": 20.0},
             {"name": "genus:816", "coef": -0.2, "mean": 40.0, "scale": 10.0},
             {"name": "frac_outside:taxonomy_species", "coef": 2.0, "mean": 0.05, "scale": 0.02},
         ]}
    m.update(over)
    json.dump(m, open(p, "w"))
    return p


READS = "gut-reads-adult-global-v0.2-lenient"


def reads_assessment(tmp_path, name="A"):
    return write_assessment(str(tmp_path / name),
                            [("S1", "taxonomy_species", "820", 90.0, "within"),
                             ("S1", "taxonomy_genus", "816", 30.0, "low"),
                             ("S1", "taxonomy_species", "999", 10.0, "within"),       # not a model feature
                             ("S2", "taxonomy_species", "820", np.nan, "expected_but_missing")],
                            [("S1", "taxonomy_species", 0.09, 2, READS), ("S1", "taxonomy_genus", 0.0, 0, READS),
                             ("S2", "taxonomy_species", 0.05, 7, READS)])


# --------------------------------------------------------------------------------------------- the shipped model

def test_shipped_model_is_percentile_only():
    """The percentile-only variant is the one that may ship: the alternative's top coefficient is a batch artifact."""
    m = hs.load_model()
    assert m["feature_set"] == "pct"
    assert not [f for f in m["features"] if str(f["name"]).startswith("present:")]
    assert {k for k, _, _ in m["parsed"]} == {"pct"}
    assert {lay for _, lay, _ in m["parsed"]} <= {"taxonomy_family", "taxonomy_genus", "taxonomy_species"}
    assert all(f["scale"] > 0 for f in m["features"])
    assert np.isfinite(m["intercept"])


def test_load_model_rejects_presence_features(tmp_path):
    p = toy_model(str(tmp_path / "m.json"), features=[{"name": "present:s__[Collinsella]_massiliensis", "coef": -1.4, "mean": 0.3, "scale": 0.4}])
    with pytest.raises(ValueError, match="presence"):
        hs.load_model(p)


@pytest.mark.parametrize("bad, msg", [
    ({"intercept": None}, "non-numeric 'intercept'"),
    ({"features": []}, "no 'features' list"),
    ({"features": [{"name": "species:820", "coef": 0.1, "mean": 1.0, "scale": 0.0}]}, "cannot standardize"),
    ({"features": [{"name": "species:820", "coef": 0.1, "mean": 1.0}]}, "missing scale"),
    ({"features": [{"name": "species:820", "coef": 0.1, "mean": 1.0, "scale": 1.0},
                    {"name": "species:820", "coef": 0.2, "mean": 1.0, "scale": 1.0}]}, "twice"),
    ({"features": [{"name": "nonsense:820", "coef": 0.1, "mean": 1.0, "scale": 1.0}]}, "unknown kind"),
    ({"features": [{"name": "species820", "coef": 0.1, "mean": 1.0, "scale": 1.0}]}, "is not"),
])
def test_load_model_refuses_malformed(tmp_path, bad, msg):
    p = toy_model(str(tmp_path / "m.json"), **bad)
    with pytest.raises(ValueError, match=msg):
        hs.load_model(p)


def test_load_model_missing_and_unparseable(tmp_path):
    with pytest.raises(FileNotFoundError, match="model not found"):
        hs.load_model(str(tmp_path / "nope.json"))
    p = str(tmp_path / "broken.json"); open(p, "w").write("{not json")
    with pytest.raises(ValueError, match="not readable JSON"):
        hs.load_model(p)


# --------------------------------------------------------------------------------------------- the arithmetic

def test_score_is_the_decision_function(tmp_path):
    d = reads_assessment(tmp_path)
    m = hs.load_model(toy_model(str(tmp_path / "m.json")))
    r = hs.score_assessment(d, m).set_index("analysis_id")
    want_s1 = -0.5 + 0.4 * (90 - 50) / 20 + -0.2 * (30 - 40) / 10 + 2.0 * (0.09 - 0.05) / 0.02
    assert r.loc["S1", "score"] == pytest.approx(want_s1)
    assert r.loc["S1", "n_features_present"] == 3 and r.loc["S1", "n_features_absent"] == 0
    # S2 carries one feature, expected-but-missing, so its percentile entered the fit as 0 and does so here
    want_s2 = -0.5 + 0.4 * (hs.ABSENT - 50) / 20 + -0.2 * (hs.ABSENT - 40) / 10 + 2.0 * (0.05 - 0.05) / 0.02
    assert r.loc["S2", "score"] == pytest.approx(want_s2)
    assert r.loc["S2", "n_features_present"] == 1   # only the summary feature; the taxon's percentile is NaN
    assert r.loc["S2", "n_features_absent"] == 2


def test_absent_feature_takes_zero_not_the_mean(tmp_path):
    """0, not the standardization centre: every feature matrix in the training script is built with .fillna(0.0)."""
    d = write_assessment(str(tmp_path / "A"), [("S1", "taxonomy_species", "820", 90.0, "within")],
                         [("S1", "taxonomy_species", 0.05, 0, READS)])
    m = hs.load_model(toy_model(str(tmp_path / "m.json")))
    e = hs.explain(d, "S1", m).set_index("feature")
    assert e.loc["genus:816", "value"] == 0.0 and not e.loc["genus:816", "present"]
    assert e.loc["genus:816", "contribution"] == pytest.approx(-0.2 * (0.0 - 40) / 10)


def test_explain_ranks_by_absolute_contribution_and_sums_to_the_score(tmp_path):
    d = reads_assessment(tmp_path)
    m = hs.load_model(toy_model(str(tmp_path / "m.json")))
    e = hs.explain(d, "S1", m, top=None)
    c = e["contribution"].to_numpy()
    assert (np.abs(c)[:-1] >= np.abs(c)[1:]).all()
    assert e.attrs["intercept"] + c.sum() == pytest.approx(hs.score_assessment(d, m).set_index("analysis_id").loc["S1", "score"])
    assert len(hs.explain(d, "S1", m, top=2)) == 2
    assert e.loc[0, "feature"] == "frac_outside:taxonomy_species"   # the largest single move in this toy model


def test_explain_unknown_sample_names_what_is_there(tmp_path):
    d = reads_assessment(tmp_path)
    m = hs.load_model(toy_model(str(tmp_path / "m.json")))
    with pytest.raises(ValueError, match="was not placed"):
        hs.explain(d, "S9", m)


def test_shipped_model_scores_a_real_shaped_assessment(tmp_path):
    """The whole shipped model against an assessment carrying one of its features, checked against the plain formula."""
    m = hs.load_model()
    f = m["features"][0]; layer = m["parsed"][0][1]; fid = m["parsed"][0][2]
    d = write_assessment(str(tmp_path / "A"), [("S1", layer, fid, 77.0, "within")], [("S1", layer, 0.05, 0, READS)])
    r = hs.score_assessment(d, m)
    rest = sum(float(g["coef"]) * (hs.ABSENT - float(g["mean"])) / float(g["scale"]) for g in m["features"][1:])
    assert r.loc[0, "score"] == pytest.approx(m["intercept"] + f["coef"] * (77.0 - f["mean"]) / f["scale"] + rest)
    assert r.loc[0, "n_features_present"] == 1
    assert len(r.loc[0, "top_features"].split("; ")) == hs.TOP_IN_SUMMARY


# --------------------------------------------------------------------------------------------- the refusals

def test_assembly_assessment_is_refused(tmp_path):
    d = write_assessment(str(tmp_path / "asm"),
                         [("S1", "ko_eggnog", "K00001", 50.0, "within"), ("S1", "taxonomy_genus", "816", 60.0, "within")],
                         [("S1", "ko_eggnog", 0.02, 0, "gut-assembly-adult-global-v0.8-lenient"),
                          ("S1", "taxonomy_genus", 0.04, 1, "gut-assembly-adult-global-v0.8-lenient")])
    for fn in (lambda: hs.score_assessment(d), lambda: hs.explain(d, "S1")):
        with pytest.raises(ValueError, match="read pipeline only"):
            fn()


def test_assembly_is_refused_on_layers_alone(tmp_path):
    """A bundle id that does not follow the gut-assembly-* naming is still caught by the layers it scored."""
    d = write_assessment(str(tmp_path / "asm2"), [("S1", "pfam", "PF00001", 50.0, "within")],
                         [("S1", "pfam", 0.02, 0, "local-rebuild")])
    with pytest.raises(ValueError, match="read pipeline only"):
        hs.score_assessment(d)


def test_unidentifiable_assessment_is_refused(tmp_path):
    """Taxonomy-only and no recognizable bundle id: nothing says this is read-based, so no number is produced."""
    d = write_assessment(str(tmp_path / "amb"), [("S1", "taxonomy_genus", "816", 50.0, "within")],
                         [("S1", "taxonomy_genus", 0.02, 0, "local-rebuild")])
    with pytest.raises(ValueError, match="cannot identify as read-based"):
        hs.score_assessment(d)


def test_no_shared_features_is_refused(tmp_path):
    d = write_assessment(str(tmp_path / "A"), [("S1", "taxonomy_species", "4242", 50.0, "within")],
                         [("S1", "taxonomy_species", 0.05, 0, READS)])
    m = hs.load_model(toy_model(str(tmp_path / "m.json"), features=[{"name": "species:820", "coef": 1.0, "mean": 50.0, "scale": 20.0}]))
    with pytest.raises(ValueError, match="do not share a feature space"):
        hs.score_assessment(d, m)


def test_missing_assessment_files(tmp_path):
    with pytest.raises(FileNotFoundError, match="assessed directory not found"):
        hs.score_assessment(str(tmp_path / "nowhere"))
    d = str(tmp_path / "half"); os.makedirs(d)
    with pytest.raises(FileNotFoundError, match="no scores.parquet"):
        hs.score_assessment(d)
    write_assessment(d, [("S1", "taxonomy_species", "820", 90.0, "within")], [("S1", "taxonomy_species", 0.05, 0, READS)])
    os.remove(os.path.join(d, "summaries.tsv"))
    with pytest.raises(FileNotFoundError, match="no summaries.tsv"):
        hs.score_assessment(d)


def test_no_library_warning_escapes(tmp_path):
    """The one-line stderr contract is load-bearing, so a numpy or pandas warning from here is a failure."""
    d = reads_assessment(tmp_path)
    m = hs.load_model(toy_model(str(tmp_path / "m.json")))
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        hs.score_assessment(d, m)
        hs.explain(d, "S2", m)
        hs.load_model()


def test_raw_json_model_is_accepted(tmp_path):
    """A model handed over as a plain json.load dict is validated on the way in rather than scored blindly."""
    p = toy_model(str(tmp_path / "m.json"))
    d = reads_assessment(tmp_path)
    assert hs.score_assessment(d, json.load(open(p))).loc[0, "score"] == pytest.approx(
        hs.score_assessment(d, hs.load_model(p)).loc[0, "score"])
    with pytest.raises(ValueError, match="presence"):
        hs.score_assessment(d, {"intercept": 0.0, "features": [{"name": "present:s__x", "coef": 1.0, "mean": 0.0, "scale": 1.0}]})
