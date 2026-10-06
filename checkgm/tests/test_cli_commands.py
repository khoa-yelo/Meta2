"""The five commands of the public interface: setup, profile, assess, score, end-to-end.

These test the wiring rather than the modules behind it (which have their own files): that each command reaches the
right stage, that `assess` works out the kind of input from what the path holds, that the renamed `score` redirects
instead of reporting an argparse error, and that the two commands which can start a very large or very long job refuse
up front and in one line when their prerequisites are absent.
"""
import json
import os
import subprocess
import sys

import pandas as pd
import pytest

from test_checkgm import make_bundle, write_profile


def run(args, env=None):
    e = dict(os.environ)
    e.pop("CHECKGM_IMAGES", None); e.pop("CHECKGM_DBS", None); e.pop("CHECKGM_IMAGE", None)
    e.update(env or {})
    return subprocess.run([sys.executable, "-m", "checkgm.cli"] + [str(a) for a in args],
                          capture_output=True, text=True, env=e)


@pytest.fixture(scope="module")
def reads(tmp_path_factory):
    root = tmp_path_factory.mktemp("reads")
    b = make_bundle(str(root / "bundles"))
    write_profile(str(root / "S1.txt"))
    return str(root), b


def test_assess_detects_a_metaphlan_profile_and_writes_the_report(reads, tmp_path):
    root, b = reads
    out = tmp_path / "assessed"
    r = run(["assess", "--bundle", os.path.basename(b), "--input", f"{root}/S1.txt", "--out", out],
            {"CHECKGM_BUNDLES": f"{root}/bundles"})
    assert r.returncode == 0, r.stderr
    for f in ("scores.parquet", "summaries.tsv", "rejections.tsv", "report.md"):
        assert (out / f).is_file(), f
    assert (out / "figures" / "S1_range_report.png").is_file()      # figures are on by default
    assert "figures: 1 written" in r.stdout


def test_assess_no_figures_writes_the_tables_only(reads, tmp_path):
    root, b = reads
    out = tmp_path / "assessed"
    r = run(["assess", "--bundle", os.path.basename(b), "--input", f"{root}/S1.txt", "--out", out, "--no-figures"],
            {"CHECKGM_BUNDLES": f"{root}/bundles"})
    assert r.returncode == 0, r.stderr
    assert (out / "scores.parquet").is_file() and not (out / "figures").exists()


def test_assess_rescores_a_directory_it_already_wrote(reads, tmp_path):
    """An assess output is a valid --input: it holds normalized/, so it is scored again without re-normalizing. This is
    how a user re-scores against another baseline or calibration without repeating the expensive half."""
    root, b = reads
    first, again = tmp_path / "a", tmp_path / "b"
    env = {"CHECKGM_BUNDLES": f"{root}/bundles"}
    assert run(["assess", "--bundle", os.path.basename(b), "--input", f"{root}/S1.txt", "--out", first, "--no-figures"], env).returncode == 0
    r = run(["assess", "--bundle", os.path.basename(b), "--input", first, "--out", again, "--no-figures"], env)
    assert r.returncode == 0, r.stderr
    one = pd.read_parquet(first / "scores.parquet").sort_values(["analysis_id", "layer", "feature_id"]).reset_index(drop=True)
    two = pd.read_parquet(again / "scores.parquet").sort_values(["analysis_id", "layer", "feature_id"]).reset_index(drop=True)
    pd.testing.assert_frame_equal(one, two)


def test_assess_refuses_an_input_it_cannot_identify(reads, tmp_path):
    root, b = reads
    d = tmp_path / "junk"; d.mkdir(); (d / "readme.txt").write_text("nothing to score here\n")
    r = run(["assess", "--bundle", os.path.basename(b), "--input", d, "--out", tmp_path / "o"],
            {"CHECKGM_BUNDLES": f"{root}/bundles"})
    assert r.returncode == 1 and r.stdout == ""
    assert len(r.stderr.strip().splitlines()) == 1 and r.stderr.startswith("checkgm: ")


def test_assess_refuses_a_mixture_of_input_kinds(reads, tmp_path):
    """Two kinds in one --out would mix two normalizers' calibrations without saying so."""
    root, b = reads
    q = tmp_path / "q"; q.mkdir()
    for f in ("contigs.tsv", "cds.tsv", "cds_taxonomy.tsv"):
        (q / f).write_text("x\n")
    r = run(["assess", "--bundle", os.path.basename(b), "--input", f"{root}/S1.txt", q, "--out", tmp_path / "o"],
            {"CHECKGM_BUNDLES": f"{root}/bundles"})
    assert r.returncode == 1 and "not all of one kind" in r.stderr
    assert len(r.stderr.strip().splitlines()) == 1


def test_score_redirects_the_flags_it_used_to_take(reads, tmp_path):
    """'score' meant "compare normalized tables with a bundle" before the redesign. Its old flags must produce a
    one-line redirect naming both replacements, not argparse's "unrecognized arguments"."""
    root, b = reads
    r = run(["score", "--bundle", os.path.basename(b), "--normalized", tmp_path, "--out", tmp_path / "o"],
            {"CHECKGM_BUNDLES": f"{root}/bundles"})
    assert r.returncode == 1 and r.stdout == ""
    assert len(r.stderr.strip().splitlines()) == 1
    assert "assess" in r.stderr and "score-normalized" in r.stderr


def test_score_needs_an_assessment(reads, tmp_path):
    r = run(["score", "--out", tmp_path / "o"], {"CHECKGM_BUNDLES": f"{reads[0]}/bundles"})
    assert r.returncode == 1 and "--assessed" in r.stderr


def test_score_runs_on_an_assessment_and_says_when_the_model_barely_overlaps(reads, tmp_path):
    """The test baseline shares almost no features with the shipped model, so the score is mostly the intercept. That
    is legitimate arithmetic and must still produce a number, but it must not pass silently.

    The threshold matters as much as the warning. An absent feature is filled with the value the fit used for one, and
    sparse presence is normal here: over the 9,123 scored cohort samples the median carries 34 of the 104 features and
    97% carry fewer than half, so a "most features absent" rule would fire on nearly every real sample and tell the
    user something false. It is set at the 1st percentile of that distribution instead, and
    test_score_does_not_warn_on_an_ordinary_sample below pins the quiet side."""
    root, b = reads
    out = tmp_path / "assessed"
    env = {"CHECKGM_BUNDLES": f"{root}/bundles"}
    assert run(["assess", "--bundle", os.path.basename(b), "--input", f"{root}/S1.txt", "--out", out, "--no-figures"], env).returncode == 0
    r = run(["score", "--assessed", out, "--out", tmp_path / "hs"], env)
    assert r.returncode == 0, r.stderr
    t = pd.read_csv(tmp_path / "hs" / "health_score.tsv", sep="\t")
    assert len(t) == 1 and t["score"].notna().all()
    assert "fewer than 12 of the model's 104 features" in r.stderr


def test_setup_dry_run_transfers_nothing_and_states_the_size(tmp_path):
    r = run(["setup", "bundles", "--dest", tmp_path / "dl", "--dry-run"])
    assert r.returncode == 0, r.stderr
    assert "would download" in r.stdout and "MB" in r.stdout
    assert not any((tmp_path / "dl").glob("*.tar.gz")) if (tmp_path / "dl").exists() else True


def test_setup_db_states_the_size_before_a_very_large_transfer(tmp_path):
    """The assembly databases are 90 GB over the wire. --dry-run must print that, and the commands it would run."""
    r = run(["setup", "db", "--pipeline", "assembly", "--dest", tmp_path / "dbs", "--dry-run"])
    assert r.returncode == 0, r.stderr
    assert "to download" in r.stdout and "on disk" in r.stdout
    assert "would run:" in r.stdout and "fetch_dbs.sh" in r.stdout
    r2 = run(["setup", "db", "--pipeline", "read", "--dest", tmp_path / "d2", "--dry-run"])
    assert r2.returncode == 0 and "fetch_dbs_B.sh" in r2.stdout


def test_profile_refuses_in_one_line_without_an_image_or_databases(reads, tmp_path):
    """checkgm installs by pip, but the pipelines need a container image and (for assembly) 170 GB of databases. A user
    without them gets one line naming what is missing and how to get it, not a crash from inside a shell script."""
    root, _ = reads
    for pipeline in ("assembly", "read"):
        r = run(["profile", "--pipeline", pipeline, "--sample", "S1", "--out", tmp_path / pipeline, "--r1", f"{root}/S1.txt"])
        assert r.returncode == 1 and r.stdout == ""
        assert len(r.stderr.strip().splitlines()) == 1, r.stderr
        assert r.stderr.startswith("checkgm: ") and "checkgm setup db" in r.stderr
        assert not (tmp_path / pipeline).exists() or not any((tmp_path / pipeline).iterdir())


def test_every_public_command_has_help_and_appears_in_the_top_level_help():
    top = run(["--help"])
    assert top.returncode == 0
    for cmd in ("setup", "profile", "assess", "score", "end-to-end"):
        assert any(l.strip() == cmd or l.strip().startswith(cmd + " ") for l in top.stdout.splitlines()), cmd
        assert run([cmd, "--help"]).returncode == 0, cmd
    for target in ("bundles", "db"):
        assert run(["setup", target, "--help"]).returncode == 0, target


@pytest.mark.skipif(not os.path.isdir("/home/classes/bios/270/khoa/meta2/project/refs/gut-reads-adult-global-v0.2-lenient"),
                    reason="needs the real read-based baseline, which is not redistributed with the source")
def test_score_does_not_warn_on_an_ordinary_sample(tmp_path):
    """The quiet side of the threshold above, on the real baseline and the README's own example: this profile carries
    40 of the 104 features, above the cohort's 1st percentile of 12, so it must be scored without a note. An earlier
    version of the rule used half the features and warned here, which would have been wrong on 97% of real samples."""
    refs = "/home/classes/bios/270/khoa/meta2/project/refs"
    out = tmp_path / "assessed"
    env = {"CHECKGM_BUNDLES": refs}
    a = run(["assess", "--bundle", "gut-reads-adult-global-v0.2-lenient", "--input", "examples/synthetic_healthy_adult.txt",
             "--reads-tsv", "examples/reads.tsv", "--out", out, "--no-figures"], env)
    assert a.returncode == 0, a.stderr
    r = run(["score", "--assessed", out, "--out", tmp_path / "hs"], env)
    assert r.returncode == 0, r.stderr
    assert r.stderr.strip() == "", r.stderr
    t = pd.read_csv(tmp_path / "hs" / "health_score.tsv", sep="\t")
    assert len(t) == 1 and 12 <= int(t["n_features_present"].iloc[0]) <= 104
