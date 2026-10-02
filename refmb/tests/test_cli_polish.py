"""Command-line and library behaviour added in 0.8.1: --calibration is checked before any work, superseded bundle series
warn, rejected profiles emit no numpy warning, load_config needs a path, and the report drops all-NaN landscape columns."""
import json
import shutil
import subprocess
import sys
import warnings

import numpy as np
import pandas as pd
import pytest

import test_refmb as R
from refmb.cli import bundle_series, write_report
from refmb.normalize_reads import ReadsBundleSpec, normalize_profile
from refmb.score import Bundle


@pytest.fixture(scope="module")
def bundle(tmp_path_factory):
    return R.make_bundle(str(tmp_path_factory.mktemp("bundles")))


def run(*args):
    return subprocess.run([sys.executable, "-m", "refmb.cli", *args], capture_output=True, text=True)


def test_missing_calibration_fails_before_anything_is_written(bundle, tmp_path):
    prof = tmp_path / "S1.txt"; R.write_profile(str(prof)); out = tmp_path / "rep"
    r = run("run-metaphlan", "--bundle", bundle, "--input", str(prof), "--calibration", str(tmp_path / "nonexistent"), "--out", str(out))
    assert r.returncode == 1 and r.stdout == "" and r.stderr.startswith("refmb: --calibration directory not found")
    assert not out.exists()
    r = run("score", "--bundle", bundle, "--normalized", str(tmp_path), "--calibration", str(tmp_path / "nonexistent"), "--out", str(out))
    assert r.returncode == 1 and not out.exists()


def test_superseded_series_warns_and_still_scores(bundle, tmp_path):
    old = tmp_path / "gut-reads-test-v0.1-lenient"; shutil.copytree(bundle, old)
    m = json.load(open(old / "manifest.json")); m["bundle_id"] = old.name; json.dump(m, open(old / "manifest.json", "w"))
    prof = tmp_path / "S1.txt"; R.write_profile(str(prof))
    r = run("run-metaphlan", "--bundle", str(old), "--input", str(prof), "--out", str(tmp_path / "rep"))
    assert r.returncode == 0, r.stderr
    assert r.stderr.strip().splitlines() == [f"refmb: bundle series v0.1 is superseded by v0.2 (bundle '{old.name}'); results from it are not the validated ones"]
    assert (tmp_path / "rep" / "scores.parquet").exists()
    r = run("run-metaphlan", "--bundle", bundle, "--input", str(prof), "--out", str(tmp_path / "rep2"))   # current series: silent
    assert r.returncode == 0 and r.stderr == ""
    assert bundle_series("gut-assembly-adult-global-v0.8-lenient-loso-PRJEB1") == "0.8" and bundle_series("no-version") is None


def test_rejected_profile_emits_no_runtime_warning(bundle):
    spec = ReadsBundleSpec(bundle)
    t = pd.DataFrame({"species_name": ["s__Unknownus_novus"], "taxid": ["424242"], "rel_abundance": [1.0]})
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        n = normalize_profile("odd", t, 1_000_000, spec)
    assert n["qc"]["status"] == "LOW_MAPPED_FRACTION"
    assert not [x for x in w if issubclass(x.category, RuntimeWarning)]


def test_load_config_requires_a_path(tmp_path):
    from refmb.calibrate import load_config
    with pytest.raises(ValueError, match="docs/configs/s11.yaml"):
        load_config()
    p = tmp_path / "c.yaml"; p.write_text("model: {cv_folds: 5}\n")
    assert load_config(str(p)) == {"model": {"cv_folds": 5}}


def test_report_drops_all_nan_landscape_columns(bundle, tmp_path):
    cols = ["analysis_id", "feature_id", "band", "value", "percentile", "call", "call_fdr", "p_two_sided", "fdr_q", "layer"]
    scores = pd.DataFrame([["S1", "815", "high", 0.1, 55.0, "within", "within", 0.9, 0.9, "taxonomy_family"]], columns=cols)
    S = pd.DataFrame([dict(analysis_id="S1", layer="taxonomy_family", n_assessed=6, n_low_raw=0, n_high_raw=0, frac_outside_raw=0.0,
                           excess_outside_p=1.0, excess_outside_q=1.0, excess_outside_ratio=0.0, quality_band="high",
                           genome_equivalents_cov=np.nan, core_distance_pct=48.2, neighbor_studies='{"P1": 5}')])
    write_report({"scores": scores, "summaries": S, "rejections": pd.DataFrame(columns=["analysis_id", "reason_code", "detail"])}, Bundle(bundle), str(tmp_path), "none")
    rep = (tmp_path / "report.md").read_text()
    land = rep.split("## Landscape")[1].split("##")[0]
    assert "core_distance_pct" in land and "genome_equivalents_cov" not in land
