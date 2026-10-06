"""Tests for checkgm.report_figures: a synthetic bundle and a synthetic assess output directory, no external data.

The figures themselves cannot be asserted on pixel by pixel, so what is tested is everything around the drawing: which
features the rule picks, that the files are written in the asked-for format, that a sample which cannot be placed is
skipped rather than raised over, that user errors are single plain messages, and that nothing reaches stderr.
"""
import json
import os

import numpy as np
import pandas as pd
import pytest

from checkgm import report_figures as RF

GRID = [1, 2.5, 5, 10, 25, 50, 75, 90, 95, 97.5, 99]
Z = np.array([-2.33, -1.96, -1.64, -1.28, -0.67, 0.0, 0.67, 1.28, 1.64, 1.96, 2.33])   # normal quantiles of GRID
PCOLS = [f"p{str(p).replace('.', '_')}" for p in GRID]
FAMILIES = [str(1000 + i) for i in range(24)]
KOS = [f"K{i:05d}" for i in range(1, 31)]
MODULES = [f"M{i:05d}" for i in range(1, 13)]
BUNDLE_ID = "gut-assembly-test-v0.8-lenient"


def _feature_rows(ids, mu, sd, banded, prev=0.9):
    rows = []
    for i, f in enumerate(ids):
        d = {"feature_id": f, "n_detected": 90, "n_studies": 3, "prevalence_weighted": prev, "percentiles_available": True,
             "layer": "x", "reliability_tier": "NA"}
        for b in ("low", "medium", "high"):
            d[f"prevalence_band_{b}"] = prev - 0.01 * i   # distinct, so "most prevalent" has one answer
            d[f"n_band_{b}"] = 30
        for c, z in zip(PCOLS, Z):
            d[c] = mu[i] + sd * z
            if banded:
                for b in ("low", "medium", "high"):
                    d[f"{c}_band_{b}"] = mu[i] + sd * z          # same shape per band here; the point is that the columns exist
        rows.append(d)
    return pd.DataFrame(rows)


@pytest.fixture(scope="module")
def bundle_dir(tmp_path_factory):
    """An assembly-shaped bundle: a CLR taxonomy layer (pooled percentiles) and a copies-per-genome layer scored within bands."""
    root = tmp_path_factory.mktemp("bundles")
    b = root / BUNDLE_ID
    (b / "features").mkdir(parents=True)
    (b / "normalization" / "backbone").mkdir(parents=True)
    rng = np.random.default_rng(0)
    _feature_rows(FAMILIES, rng.normal(2.0, 0.5, len(FAMILIES)), 0.8, banded=False).to_parquet(b / "features" / "taxonomy_family.parquet", index=False)
    _feature_rows(KOS, rng.uniform(0.8, 2.0, len(KOS)), 0.25, banded=True).to_parquet(b / "features" / "ko_eggnog.parquet", index=False)
    _feature_rows(MODULES, rng.uniform(0.4, 0.9, len(MODULES)), 0.08, banded=False).to_parquet(b / "features" / "module.parquet", index=False)
    pd.DataFrame({"taxid": [int(f) for f in FAMILIES], "name": [f"Testaceae{i}" for i in range(len(FAMILIES))]}
                 ).to_parquet(b / "normalization" / "backbone" / "ncbi_taxonomy.parquet", index=False)
    json.dump({"bundle_id": BUNDLE_ID, "n_samples": 1234, "n_studies": 17, "percentile_grid": GRID,
               "layers": {"taxonomy_family": {"n_features": len(FAMILIES)}, "ko_eggnog": {"n_features": len(KOS)},
                          "module": {"n_features": len(MODULES)}},
               "profile_type": "assembly", "source_pipeline": {"name": "test", "version": "0"},
               "reference_key": {"body_site": "Gut", "stratum": "adult-global"},
               "scoring": {"detection_floor_cpg": 0.05, "percentiles_within_band": True,
                           "percentiles_within_band_layers": ["ko_eggnog"]}},
              open(b / "manifest.json", "w"))
    return str(root), str(b)


def _scores_for(sample, bundle, n_low=4, n_high=3, n_missing=3, n_na=2):
    """Score rows for one sample: the first n_low features called low, then high, then expected-but-missing, then
    not-assessable, the rest within. Values are placed where the call says they are."""
    rows = []
    for layer, ids in (("taxonomy_family", FAMILIES), ("ko_eggnog", KOS), ("module", MODULES)):
        ref = pd.read_parquet(os.path.join(bundle, "features", f"{layer}.parquet")).set_index("feature_id")
        calls = (["low"] * n_low + ["high"] * n_high + ["expected_but_missing"] * n_missing + ["not_assessable"] * n_na)
        calls += ["within"] * (len(ids) - len(calls))
        for f, call in zip(ids, calls):
            med, lo, hi = ref.loc[f, "p50"], ref.loc[f, "p2_5"], ref.loc[f, "p97_5"]
            v, p = {"low": (lo - 1.0, 0.5), "high": (hi + 1.0, 99.5), "within": (med, 50.0),
                    "expected_but_missing": (np.nan, np.nan), "not_assessable": (np.nan, np.nan)}[call]
            rows.append(dict(analysis_id=sample, feature_id=f, band="medium", value=v, percentile=p, call=call,
                             call_fdr=call, p_two_sided=np.nan, fdr_q=np.nan, layer=layer))
    return pd.DataFrame(rows)


def _summaries_for(samples, scores):
    rows = []
    for (sid, layer), g in scores.groupby(["analysis_id", "layer"]):
        c = g["call"]
        rows.append(dict(analysis_id=sid, layer=layer, n_detected=int((c != "expected_but_missing").sum()),
                         n_assessed=int(c.isin(["low", "high", "within"]).sum()), n_low_raw=int((c == "low").sum()),
                         n_high_raw=int((c == "high").sum()), n_missing=int((c == "expected_but_missing").sum()),
                         n_not_assessable=int((c == "not_assessable").sum()), quality_band="medium",
                         bundle_id=BUNDLE_ID, excess_outside_q=0.01))
    return pd.DataFrame(rows)


@pytest.fixture
def assessed(tmp_path, bundle_dir, monkeypatch):
    """An assess output directory for two samples, with CHECKGM_BUNDLES pointing at the bundle it names."""
    root, bundle = bundle_dir
    monkeypatch.setenv("CHECKGM_BUNDLES", root)
    out = tmp_path / "assessed"
    out.mkdir()
    sc = pd.concat([_scores_for("S1", bundle), _scores_for("S2", bundle, n_low=0, n_high=5)], ignore_index=True)
    sc.to_parquet(out / "scores.parquet", index=False)
    _summaries_for(["S1", "S2"], sc).to_csv(out / "summaries.tsv", sep="\t", index=False)
    pd.DataFrame([{"analysis_id": "S3", "reason_code": "NOT_NORMALIZED", "detail": "mapped fraction below the floor"}]
                 ).to_csv(out / "rejections.tsv", sep="\t", index=False)
    return str(out)


# ------------------------------------------------------------------------------------------------- optional dependency

def test_available_matches_the_import():
    import importlib.util
    assert RF.available() is (importlib.util.find_spec("matplotlib") is not None)


def test_available_is_false_without_matplotlib(monkeypatch):
    """`assess` must be able to skip figures instead of crashing on an install without matplotlib."""
    monkeypatch.setattr(RF.importlib.util, "find_spec", lambda name: None)
    assert RF.available() is False


# ------------------------------------------------------------------------------------------------- the selection rule

def _panel(assessed, layer, sample="S1", cap=10):
    sc = pd.read_parquet(os.path.join(assessed, "scores.parquet")).astype({"analysis_id": str, "feature_id": str})
    sc = sc[(sc["analysis_id"] == sample) & (sc["layer"] == layer) & (sc["call"] != "not_assessable")].set_index("feature_id")
    bundle = RF._bundle_of(assessed, pd.read_csv(os.path.join(assessed, "summaries.tsv"), sep="\t"))
    ref = bundle.features(layer).rename(index=str)
    return RF._panel_rows(sc, ref, bundle, layer, "medium", cap, RF._taxon_names(bundle, sc.index))


def test_quotas_sum_to_the_cap():
    for cap in (4, 5, 10, 20):
        q = RF._quotas(cap)
        assert sum(n for _, n in q) == cap and all(n >= 0 for _, n in q)
    assert dict(RF._quotas(10)) == {"low": 3, "high": 3, "missing": 2, "within": 2}
    assert dict(RF._quotas(20)) == {"low": 6, "high": 6, "missing": 4, "within": 4}


def test_selection_is_symmetric_and_keeps_context(assessed):
    rows = _panel(assessed, "taxonomy_family")
    calls = pd.Series([r["call"] for r in rows]).value_counts().to_dict()
    assert len(rows) == 10
    assert calls == {"low": 3, "high": 3, "expected_but_missing": 2, "within": 2}
    assert [r["call"] for r in rows][:2] == ["expected_but_missing"] * 2        # absent features read first
    xs = [r["x"] for r in rows if np.isfinite(r["x"])]
    assert xs == sorted(xs)                                                     # then lowest to highest, like a report


def test_unused_quota_goes_to_the_other_directions(assessed):
    """S2 has no low calls: the figure shows more high ones, not two rows of context."""
    rows = _panel(assessed, "taxonomy_family", sample="S2")
    calls = pd.Series([r["call"] for r in rows]).value_counts().to_dict()
    assert calls.get("low", 0) == 0 and calls["high"] == 5 and len(rows) == 10


def test_not_assessable_features_are_never_shown(assessed):
    rows = _panel(assessed, "taxonomy_family")
    assert all(r["call"] != "not_assessable" for r in rows)


def test_selection_is_deterministic(assessed):
    assert [r["fid"] for r in _panel(assessed, "taxonomy_family")] == [r["fid"] for r in _panel(assessed, "taxonomy_family")]


def test_rows_carry_names_and_the_right_units(assessed):
    fam = _panel(assessed, "taxonomy_family")
    assert all(r["label"].startswith("Testaceae") for r in fam)                 # taxids resolved through the bundle backbone
    within = next(r for r in fam if r["call"] == "within")
    assert within["x"] == pytest.approx(0.0, abs=1e-9)                          # a value at the median sits on the median tick
    assert within["lo"] == pytest.approx(-1.96 * 0.8 / np.log(2), rel=1e-6)     # CLR: a difference, shown in log2 units
    ko = _panel(assessed, "ko_eggnog")
    kw = next(r for r in ko if r["call"] == "within")
    assert kw["x"] == pytest.approx(0.0, abs=1e-9) and kw["lo"] < 0             # copies per genome: a log2 fold difference


def test_a_bounded_layer_is_plotted_as_a_difference(assessed):
    """Module completeness is a fraction: a fold difference would be meaningless, so the row is the plain difference."""
    rows = _panel(assessed, "module")
    assert RF.MODE_OF_VALUE[RF.LAYER_SPEC["module"][0]] == "linear"
    low = next(r for r in rows if r["call"] == "low")
    assert low["x"] == pytest.approx(low["lo"] - 1.0)       # _scores_for puts a low value one unit below p2_5
    assert "completeness" in RF.X_LABEL["linear"]


def test_band_percentiles_are_preferred_and_fall_back(bundle_dir):
    _, bundle = bundle_dir
    ref = pd.read_parquet(os.path.join(bundle, "features", "ko_eggnog.parquet")).set_index("feature_id")
    row = ref.loc[KOS[0]].copy()
    row["p50_band_medium"] = 99.0
    assert RF._quantiles(row, "medium", banded=True)["med"] == 99.0
    row["p2_5_band_medium"] = np.nan                                            # an incomplete band falls back to the pooled columns
    assert RF._quantiles(row, "medium", banded=True)["med"] == pytest.approx(ref.loc[KOS[0], "p50"])
    assert RF._quantiles(pd.Series({"p50": 1.0}), None, banded=False) is None   # no range edges: not drawable


# ------------------------------------------------------------------------------------------------- writing the figures

@pytest.mark.skipif(not RF.available(), reason="matplotlib is not installed")
def test_writes_one_figure_per_sample(assessed, tmp_path, capsys):
    paths = RF.write_figures(assessed, str(tmp_path / "figs"))
    assert [os.path.basename(p) for p in paths] == ["S1_range_report.png", "S2_range_report.png"]
    for p in paths:
        assert open(p, "rb").read(8) == b"\x89PNG\r\n\x1a\n" and os.path.getsize(p) > 5000
    captured = capsys.readouterr()
    assert captured.out == "" and captured.err == ""      # the one-line stderr contract: no library warning may leak out


@pytest.mark.skipif(not RF.available(), reason="matplotlib is not installed")
def test_samples_argument_selects_and_orders(assessed, tmp_path):
    paths = RF.write_figures(assessed, str(tmp_path / "figs"), samples=["S2"])
    assert len(paths) == 1 and os.path.basename(paths[0]) == "S2_range_report.png"


@pytest.mark.skipif(not RF.available(), reason="matplotlib is not installed")
@pytest.mark.parametrize("fmt,head", [("pdf", b"%PDF"), ("svg", b"<?xml")])
def test_other_formats(assessed, tmp_path, fmt, head):
    paths = RF.write_figures(assessed, str(tmp_path / fmt), samples=["S1"], fmt=fmt)
    assert paths[0].endswith("." + fmt) and open(paths[0], "rb").read(5).startswith(head)


@pytest.mark.skipif(not RF.available(), reason="matplotlib is not installed")
def test_a_sample_that_cannot_be_placed_is_skipped_not_raised(assessed, tmp_path):
    """A rejected sample has no score rows at all, and a sample whose only rows are not assessable has nothing to draw;
    neither is an error, so figures are written for the rest and the skipped sample is simply absent."""
    sc = pd.read_parquet(os.path.join(assessed, "scores.parquet"))
    na = sc[sc["analysis_id"] == "S1"].copy()
    na["analysis_id"] = "S4"; na["call"] = "not_assessable"; na["percentile"] = np.nan
    pd.concat([sc, na], ignore_index=True).to_parquet(os.path.join(assessed, "scores.parquet"), index=False)
    paths = RF.write_figures(assessed, str(tmp_path / "figs"))
    assert [os.path.basename(p) for p in paths] == ["S1_range_report.png", "S2_range_report.png"]
    assert RF.write_figures(assessed, str(tmp_path / "figs2"), samples=["S4"]) == []


@pytest.mark.skipif(not RF.available(), reason="matplotlib is not installed")
def test_single_layer_bundle_gets_the_whole_row_budget(assessed, tmp_path, bundle_dir):
    """With one scorable layer the figure shows up to twenty rows rather than ten and half an empty page."""
    sc = pd.read_parquet(os.path.join(assessed, "scores.parquet"))
    sc[sc["layer"] == "ko_eggnog"].to_parquet(os.path.join(assessed, "scores.parquet"), index=False)
    paths = RF.write_figures(assessed, str(tmp_path / "figs"), samples=["S1"])
    assert len(paths) == 1
    assert RF._quotas(RF.ROW_BUDGET) == [("low", 6), ("high", 6), ("missing", 4), ("within", 4)]


# ------------------------------------------------------------------------------------------------- user errors

def test_bad_format_is_a_plain_error(assessed, tmp_path):
    with pytest.raises(ValueError) as e:
        RF.write_figures(assessed, str(tmp_path / "figs"), fmt="jpeg")
    assert "jpeg" in str(e.value) and "png" in str(e.value) and "\n" not in str(e.value)


def test_unknown_sample_is_a_plain_error(assessed, tmp_path):
    with pytest.raises(ValueError) as e:
        RF.write_figures(assessed, str(tmp_path / "figs"), samples=["S1", "nope"])
    assert "nope" in str(e.value) and "\n" not in str(e.value)


def test_not_an_assess_directory(tmp_path):
    with pytest.raises(FileNotFoundError) as e:
        RF.write_figures(str(tmp_path), str(tmp_path / "figs"))
    assert "scores.parquet" in str(e.value)
    with pytest.raises(FileNotFoundError):
        RF.write_figures(str(tmp_path / "nowhere"), str(tmp_path / "figs"))


def test_unresolvable_bundle_says_where_it_looked(assessed, tmp_path, monkeypatch):
    monkeypatch.setenv("CHECKGM_BUNDLES", str(tmp_path / "empty"))
    with pytest.raises(FileNotFoundError) as e:
        RF.write_figures(assessed, str(tmp_path / "figs"))
    msg = str(e.value)
    assert BUNDLE_ID in msg and "CHECKGM_BUNDLES" in msg and "\n" not in msg


def test_assess_output_without_a_bundle_id(assessed, tmp_path):
    S = pd.read_csv(os.path.join(assessed, "summaries.tsv"), sep="\t").drop(columns=["bundle_id"])
    S.to_csv(os.path.join(assessed, "summaries.tsv"), sep="\t", index=False)
    with pytest.raises(ValueError) as e:
        RF.write_figures(assessed, str(tmp_path / "figs"))
    assert "bundle_id" in str(e.value)


def test_mixed_bundles_are_refused(assessed, tmp_path):
    S = pd.read_csv(os.path.join(assessed, "summaries.tsv"), sep="\t")
    S.loc[S.index[0], "bundle_id"] = "gut-assembly-other-v0.8-lenient"
    S.to_csv(os.path.join(assessed, "summaries.tsv"), sep="\t", index=False)
    with pytest.raises(ValueError) as e:
        RF.write_figures(assessed, str(tmp_path / "figs"))
    assert "2 bundles" in str(e.value)


# ------------------------------------------------------------------------------------------------- palette

def test_palette_is_the_paper_s():
    assert RF._colour("low", 20) == RF.DEV["low_light"] and RF._colour("low", 0.4) == RF.DEV["low_deep"]
    assert RF._colour("high", 98) == RF.DEV["high_light"] and RF._colour("high", 99.6) == RF.DEV["high_deep"]
    assert RF._colour("expected_but_missing", None) == RF.DEV["low_light"]
    assert RF._colour("within", 50) == RF.DEV["neutral"]


def test_value_and_percentile_text():
    assert RF._fold_text("log2_diff", -2.0) == "4.0x lower" and RF._fold_text("log2_ratio", 1.0) == "2.0x higher"
    assert RF._fold_text("linear", -0.25) == "-0.25 vs the median"
    assert RF._fold_text("log2_ratio", 0.0) == "at the reference median"      # a value on the median is not a 1.0x finding
    assert RF._fold_text("linear", 0.0) == "at the reference median"
    assert RF._ordinal(0.5) == "below the 1st percentile" and RF._ordinal(12) == "12th percentile"
    assert RF._ordinal(None) == "" and RF._tick_label("log2_diff", 0) == "typical" and RF._tick_label("log2_diff", -2) == "1/4"


def test_no_feature_is_drawn_twice_when_a_direction_is_empty(assessed, bundle_dir):
    """A panel must never draw one feature on two rows. S2 has no low calls, so the low pool is empty, and an empty
    pool is where this goes wrong: DataFrame.assign of a full-length Series onto an empty frame adopts the Series'
    index instead of aligning down to nothing, which silently turned the empty pool into every feature of the layer
    and drew the three worst high features a second time."""
    import pandas as pd
    from checkgm import report_figures as RF
    from checkgm.score import Bundle

    sc = pd.read_parquet(os.path.join(assessed, "scores.parquet"))
    bundle = Bundle(bundle_dir[1])
    for sid in ("S1", "S2"):
        for layer in ("taxonomy_family", "ko_eggnog"):
            g = sc[(sc.analysis_id == sid) & (sc.layer == layer) & (sc.call != "not_assessable")].set_index("feature_id")
            ref = bundle.features(layer)
            if ref.index.dtype != g.index.dtype:
                ref = ref.rename(index=str)
            sel = RF._select(g, ref, RF._prevalence(ref, "medium"), 10)
            assert len(sel) == len(set(sel)), (sid, layer, [f for f in set(sel) if sel.count(f) > 1])
            assert set(sel) <= set(g.index), (sid, layer)   # and never a feature the sample has no score row for

    rows = _panel(assessed, "taxonomy_family", sample="S2")
    assert len({r["label"] for r in rows}) == len(rows)
