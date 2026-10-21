"""checkgm.healthscore: the arithmetic of the shipped health scores, the choice of model by the assessment's pipeline,
and the refusals that matter more than either — an assessment that identifies as neither pipeline, a model fitted on
the other pipeline, and a model file carrying the species-presence features whose top coefficient is a batch artifact. The assessments are written here by hand in the shape cli.write_report writes, so nothing external is needed."""
import json
import os
import subprocess
import sys
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
ASSEMBLY = "gut-assembly-adult-global-v0.8-lenient"


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
    assert m["pipeline"] == "reads" and m["path"] == hs.default_model_path("reads")


def test_shipped_assembly_model_reads_assembly_layers_only():
    """The assembly model's features are what an assembly assessment carries: percentiles of family, genus, eggNOG KO
    and KEGG module, and per-layer summary counts of the assembly layers. Nothing of the read pipeline."""
    m = hs.load_model(pipeline="assembly")
    assert m["pipeline"] == "assembly" and m["path"] == hs.default_model_path("assembly")
    assert m["trained_on"]["studies"] == 11
    pct = {lay for k, lay, _ in m["parsed"] if k == "pct"}; summ = {lay for k, lay, _ in m["parsed"] if k == "summary"}
    assert pct <= {"taxonomy_family", "taxonomy_genus", "ko_eggnog", "module"} and pct & hs.ASSEMBLY_LAYERS
    assert summ <= {"taxonomy_family", "taxonomy_genus"} | hs.ASSEMBLY_LAYERS
    assert not [f for f in m["features"] if str(f["name"]).startswith("present:")]
    assert all(f["scale"] > 0 and f["coef"] != 0 for f in m["features"])
    assert np.isfinite(m["intercept"]) and m["thin_features"] > 0


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

def assembly_assessment(tmp_path, name="asm", bundle=ASSEMBLY):
    """One assembly-shaped sample carrying the shipped assembly model's first percentile feature (at 77) and a genus
    fraction outside of 0.2; every other model feature is absent or, for a summary count of a scored layer, 0."""
    m = hs.load_model(pipeline="assembly")
    _, layer, fid = next(t for t in m["parsed"] if t[0] == "pct")
    rows = [("S1", layer, fid, 77.0, "within"), ("S1", "ko_eggnog", "K99999", 50.0, "within")]   # the second is no feature
    summ = [("S1", lay, fo, 0, bundle) for lay, fo in (("taxonomy_genus", 0.2), ("ko_eggnog", 0.02), ("pfam", 0.03))]
    return write_assessment(str(tmp_path / name), rows, summ), m, layer, fid


def test_assembly_assessment_selects_the_assembly_model(tmp_path):
    d, m, _, _ = assembly_assessment(tmp_path)
    assert hs.model_for(d)["pipeline"] == "assembly"
    assert hs.model_for(d)["path"] == hs.default_model_path("assembly")


def test_assembly_assessment_scores_with_the_assembly_model(tmp_path):
    """The default call on an assembly assessment applies the shipped assembly model, checked against the plain formula:
    the one percentile given, genus frac_outside 0.2, every other feature 0 (absent, or a summary count of 0)."""
    d, m, layer, fid = assembly_assessment(tmp_path)
    r = hs.score_assessment(d)
    want = m["intercept"]
    for f, (kind, lay, key) in zip(m["features"], m["parsed"]):
        v = 77.0 if (kind, lay, key) == ("pct", layer, fid) else 0.2 if (kind, lay, key) == ("summary", "taxonomy_genus", "frac_outside_raw") else 0.0
        want += float(f["coef"]) * (v - float(f["mean"])) / float(f["scale"])
    assert r.loc[0, "score"] == pytest.approx(want)
    assert r.loc[0, "score"] == pytest.approx(hs.score_assessment(d, hs.load_model(pipeline="assembly")).loc[0, "score"])
    # present: the percentile, plus each summary feature of the three layers the assessment summarised
    n_summary = sum(1 for k, lay, _ in m["parsed"] if k == "summary" and lay in {"taxonomy_genus", "ko_eggnog", "pfam"})
    assert r.loc[0, "n_features_present"] == 1 + n_summary
    assert r.loc[0, "n_features_absent"] == len(m["features"]) - 1 - n_summary
    e = hs.explain(d, "S1", top=None)
    assert e.attrs["score"] == pytest.approx(want) and len(e) == len(m["features"])


def test_read_model_on_an_assembly_assessment_is_refused(tmp_path):
    d, _, _, _ = assembly_assessment(tmp_path)
    for fn in (lambda: hs.score_assessment(d, hs.load_model()), lambda: hs.explain(d, "S1", hs.load_model()),
               lambda: hs.model_for(d, hs.default_model_path("reads")),
               lambda: hs.score_assessment(d, json.load(open(hs.default_model_path("reads"))))):
        with pytest.raises(ValueError, match="read-pipeline health-score model.*assembly-pipeline baseline"):
            fn()


def test_assembly_model_on_a_read_assessment_is_refused(tmp_path):
    d = reads_assessment(tmp_path)
    with pytest.raises(ValueError, match="assembly-pipeline health-score model.*read-pipeline baseline"):
        hs.score_assessment(d, hs.load_model(pipeline="assembly"))
    assert hs.model_for(d)["pipeline"] == "reads"


def test_assembly_is_identified_on_layers_alone(tmp_path):
    """A bundle id that does not follow the gut-assembly-* naming is still identified by the layers it scored."""
    d = write_assessment(str(tmp_path / "asm2"), [("S1", "pfam", "PF00001", 50.0, "within")],
                         [("S1", "pfam", 0.02, 0, "local-rebuild")])
    assert hs.model_for(d)["pipeline"] == "assembly"
    with pytest.raises(ValueError, match="read-pipeline health-score model"):
        hs.score_assessment(d, hs.load_model())


def test_unidentifiable_assessment_is_refused(tmp_path):
    """Taxonomy-only and no recognizable bundle id: nothing says which pipeline this is, so no number is produced."""
    d = write_assessment(str(tmp_path / "amb"), [("S1", "taxonomy_genus", "816", 50.0, "within")],
                         [("S1", "taxonomy_genus", 0.02, 0, "local-rebuild")])
    for fn in (lambda: hs.score_assessment(d), lambda: hs.model_for(d)):
        with pytest.raises(ValueError, match="cannot identify as read-based or assembly-based"):
            fn()


def test_mixed_assessment_is_refused(tmp_path):
    d = write_assessment(str(tmp_path / "mix"), [("S1", "taxonomy_genus", "816", 50.0, "within")],
                         [("S1", "taxonomy_genus", 0.02, 0, READS), ("S2", "taxonomy_genus", 0.02, 0, ASSEMBLY)])
    with pytest.raises(ValueError, match="both a read-based and an assembly-based baseline"):
        hs.score_assessment(d)


def test_cli_score_on_an_assembly_assessment(tmp_path):
    """`checkgm score` on an assembly assessment uses the assembly model, and `--model` with the read model is a
    one-line refusal (exit 1), not a number."""
    d, m, _, _ = assembly_assessment(tmp_path)
    run = lambda *x: subprocess.run([sys.executable, "-m", "checkgm.cli", "score", "--assessed", d, *x], capture_output=True, text=True)
    r = run("--out", str(tmp_path / "hs"))
    assert r.returncode == 0, r.stderr
    t = pd.read_csv(tmp_path / "hs" / "health_score.tsv", sep="\t")
    assert t.loc[0, "score"] == pytest.approx(hs.score_assessment(d).loc[0, "score"])
    assert f"fewer than {m['thin_features']} of the model's {len(m['features'])}" in r.stderr and "gut-assembly baseline" in r.stderr
    r = run("--out", str(tmp_path / "hs2"), "--model", hs.default_model_path("reads"))
    assert r.returncode == 1 and r.stderr.startswith("checkgm: ") and len(r.stderr.strip().splitlines()) == 1
    assert "read-pipeline health-score model" in r.stderr


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
        hs.load_model(pipeline="assembly")


def test_raw_json_model_is_accepted(tmp_path):
    """A model handed over as a plain json.load dict is validated on the way in rather than scored blindly."""
    p = toy_model(str(tmp_path / "m.json"))
    d = reads_assessment(tmp_path)
    assert hs.score_assessment(d, json.load(open(p))).loc[0, "score"] == pytest.approx(
        hs.score_assessment(d, hs.load_model(p)).loc[0, "score"])
    with pytest.raises(ValueError, match="presence"):
        hs.score_assessment(d, {"intercept": 0.0, "features": [{"name": "present:s__x", "coef": 1.0, "mean": 0.0, "scale": 1.0}]})
