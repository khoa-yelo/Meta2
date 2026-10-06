"""User errors must come back as one plain line on stderr and exit status 1, never as a traceback; profiles without species
rows are rejected with EMPTY_PROFILE (exit 0); a header line in --reads-tsv is skipped. Uses the synthetic bundle of
test_checkgm.py plus a tiny synthetic assembly-type bundle (manifest only, enough for the profile_type check)."""
import json
import os
import subprocess
import sys

import pandas as pd
import pytest

from checkgm.cli import main, read_reads_tsv
from checkgm.normalize_reads import ReadsBundleSpec, normalize_profile, read_metaphlan
from test_checkgm import make_bundle, write_profile


@pytest.fixture(scope="module")
def bundles(tmp_path_factory):
    root = str(tmp_path_factory.mktemp("bundles"))
    reads = make_bundle(root)
    asm = os.path.join(root, "gut-assembly-test-v0.8-lenient"); os.makedirs(asm)
    json.dump({"bundle_id": os.path.basename(asm), "profile_type": "assembly", "n_samples": 10, "n_studies": 2, "percentile_grid": [1, 50, 99], "layers": {}},
              open(os.path.join(asm, "manifest.json"), "w"))
    return {"root": root, "reads": reads, "assembly": asm}


def run(args, env_root=None, cwd=None):
    env = dict(os.environ)
    env.pop("REFMB_BUNDLES", None)
    if env_root:
        env["REFMB_BUNDLES"] = env_root
    return subprocess.run([sys.executable, "-m", "checkgm.cli", *args], capture_output=True, text=True, env=env, cwd=cwd)


def assert_one_line_error(res, *fragments):
    assert res.returncode == 1, (res.returncode, res.stderr)
    lines = [l for l in res.stderr.splitlines() if l.strip()]
    assert len(lines) == 1 and lines[0].startswith("checkgm: "), res.stderr
    assert "Traceback" not in res.stderr and res.stdout == ""
    for f in fragments:
        assert f in lines[0], (f, lines[0])


def test_unknown_bundle_without_env(bundles, tmp_path):
    prof = tmp_path / "S.txt"; write_profile(str(prof))
    res = run(["run-metaphlan", "--bundle", "gut-reads-test-v0.2-lenient", "--input", str(prof), "--out", str(tmp_path / "o")])
    assert_one_line_error(res, "not found", "REFMB_BUNDLES")


def test_bundle_typo(bundles, tmp_path):
    prof = tmp_path / "S.txt"; write_profile(str(prof))
    res = run(["run-metaphlan", "--bundle", "gut-reads-test-v0.2-lenieint", "--input", str(prof), "--out", str(tmp_path / "o")], bundles["root"])
    assert_one_line_error(res, "gut-reads-test-v0.2-lenieint", "not found")


def test_missing_input_file(bundles, tmp_path):
    res = run(["run-metaphlan", "--bundle", bundles["reads"], "--input", str(tmp_path / "nope.txt"), "--out", str(tmp_path / "o")])
    assert_one_line_error(res, "--input file not found", "nope.txt")


def test_assembly_bundle_given_to_run_metaphlan(bundles, tmp_path):
    prof = tmp_path / "S.txt"; write_profile(str(prof))
    res = run(["run-metaphlan", "--bundle", bundles["assembly"], "--input", str(prof), "--out", str(tmp_path / "o")])
    assert_one_line_error(res, "assembly-based baseline", "run-metaphlan", "read-based")


def test_reads_bundle_given_to_normalize_and_run_mgnify(bundles, tmp_path):
    q = tmp_path / "q"; q.mkdir()
    res = run(["normalize", "--bundle", bundles["reads"], "--query", str(q), "--out", str(tmp_path / "o")])
    assert_one_line_error(res, "read-based baseline", "normalize", "assembly-based")
    res = run(["run-mgnify", "--bundle", bundles["reads"], "--analysis-dir", str(q), "--out", str(tmp_path / "o")])
    assert_one_line_error(res, "read-based baseline", "run-mgnify")


def test_score_checks_bundle_type_against_normalized_data(bundles, tmp_path):
    prof = tmp_path / "S.txt"; write_profile(str(prof)); out = tmp_path / "norm"
    ok = run(["import-metaphlan", "--bundle", bundles["reads"], "--input", str(prof), "--out", str(out)])
    assert ok.returncode == 0, ok.stderr
    res = run(["score", "--bundle", bundles["assembly"], "--normalized", str(out), "--out", str(tmp_path / "rep")])
    assert_one_line_error(res, "assembly-based baseline", "'score' needs a read-based")


def test_non_metaphlan_input(bundles, tmp_path):
    bad = tmp_path / "genes.tsv"; bad.write_text("gene\tvalue\nabc\t1\n")
    res = run(["run-metaphlan", "--bundle", bundles["reads"], "--input", str(bad), "--out", str(tmp_path / "o")])
    assert_one_line_error(res, "not a MetaPhlAn profile", "clade_name")
    empty = tmp_path / "empty.txt"; empty.write_text("")
    res = run(["run-metaphlan", "--bundle", bundles["reads"], "--input", str(empty), "--out", str(tmp_path / "o2")])
    assert_one_line_error(res, "not a MetaPhlAn profile")


def test_header_only_profile_is_empty_profile(bundles, tmp_path):
    prof = tmp_path / "hdr.txt"
    prof.write_text("#mpa_v30_CHOCOPhlAn_201901\n#12345678 reads processed\n#clade_name\tNCBI_tax_id\trelative_abundance\tadditional_species\n")
    parsed = read_metaphlan(str(prof))
    assert parsed["hdr"]["reads"] == 12345678 and parsed["hdr"]["table"].empty
    n = normalize_profile("hdr", parsed["hdr"]["table"], parsed["hdr"]["reads"], ReadsBundleSpec(bundles["reads"]))
    assert n["qc"]["status"] == "EMPTY_PROFILE" and n["qc"]["n_species_rows"] == 0
    out = tmp_path / "o"
    res = run(["run-metaphlan", "--bundle", bundles["reads"], "--input", str(prof), "--out", str(out)])
    assert res.returncode == 0 and "Traceback" not in res.stderr, res.stderr
    rej = pd.read_csv(out / "rejections.tsv", sep="\t")
    assert list(rej["analysis_id"]) == ["hdr"] and rej["detail"].iloc[0].startswith("EMPTY_PROFILE")


def test_missing_analysis_dir(bundles, tmp_path):
    res = run(["run-mgnify", "--bundle", bundles["assembly"], "--analysis-dir", str(tmp_path / "MGYA0"), "--out", str(tmp_path / "o")])
    assert_one_line_error(res, "--analysis-dir directory not found")


def test_score_on_directory_without_qc(bundles, tmp_path):
    d = tmp_path / "d"; d.mkdir()
    res = run(["score", "--bundle", bundles["reads"], "--normalized", str(d), "--out", str(tmp_path / "o")])
    assert_one_line_error(res, "no qc.tsv")
    res = run(["score", "--bundle", bundles["reads"], "--normalized", str(tmp_path / "absent"), "--out", str(tmp_path / "o")])
    assert_one_line_error(res, "--normalized directory not found")


def test_reads_tsv_header_and_bad_rows(tmp_path, capsys):
    f = tmp_path / "reads.tsv"
    f.write_text("sample\treads\nS1\t12,000,000\nS2\tunknown\nS3\t5000000\n")
    got = read_reads_tsv(str(f))
    assert got == {"S1": 12000000, "S3": 5000000}
    err = capsys.readouterr().err
    assert "1 row(s)" in err and "'S2'" in err and "sample" not in err.split("ignored:")[1]
    g = tmp_path / "plain.tsv"; g.write_text("S1\t100\n")
    assert read_reads_tsv(str(g)) == {"S1": 100} and capsys.readouterr().err == ""
    with pytest.raises(FileNotFoundError):
        read_reads_tsv(str(tmp_path / "none.tsv"))


def test_reads_tsv_header_through_cli(bundles, tmp_path):
    prof = tmp_path / "S9.txt"; write_profile(str(prof), reads=50_000_000)
    rt = tmp_path / "reads.tsv"; rt.write_text("sample\treads\nS9\t5000000\n")   # header line, then a depth that puts S9 in the low band
    out = tmp_path / "o"
    res = run(["run-metaphlan", "--bundle", bundles["reads"], "--input", str(prof), "--reads-tsv", str(rt), "--out", str(out)])
    assert res.returncode == 0 and res.stderr == "", res.stderr
    S = pd.read_csv(out / "summaries.tsv", sep="\t")
    assert S["quality_band"].eq("low").all()


def test_import_then_score_equals_run(bundles, tmp_path):
    """qc.tsv round trip must be exact: the two documented workflows give identical summaries.tsv and scores.parquet."""
    prof = tmp_path / "S5.txt"; write_profile(str(prof))
    one = tmp_path / "one"; two = tmp_path / "two"
    assert run(["run-metaphlan", "--bundle", bundles["reads"], "--input", str(prof), "--out", str(one)]).returncode == 0
    assert run(["import-metaphlan", "--bundle", bundles["reads"], "--input", str(prof), "--out", str(two)]).returncode == 0
    assert run(["score", "--bundle", bundles["reads"], "--normalized", str(two), "--out", str(two)]).returncode == 0
    assert (one / "summaries.tsv").read_bytes() == (two / "summaries.tsv").read_bytes()
    assert pd.read_parquet(one / "scores.parquet").equals(pd.read_parquet(two / "scores.parquet"))
    assert (one / "report.md").read_text() == (two / "report.md").read_text()


def test_main_returns_exit_status(bundles, tmp_path, monkeypatch):
    monkeypatch.delenv("REFMB_BUNDLES", raising=False)
    assert main(["run-metaphlan", "--bundle", "no-such", "--input", "x.txt", "--out", str(tmp_path / "o")]) == 1


def test_help_has_one_sentence_per_argument():
    for cmd in ["score", "run-metaphlan", "import-metaphlan", "run-mgnify", "import-mgnify", "normalize"]:
        res = run([cmd, "--help"])
        assert res.returncode == 0
        text = " ".join(res.stdout.split())
        assert "--bundle BUNDLE bundle directory" in text and "--out OUT output directory" in text, text
    top = run(["--help"]).stdout
    for cmd in ["bundles", "normalize", "import-mgnify", "import-metaphlan", "score", "run-mgnify", "run-metaphlan"]:
        # Python 3.10's argparse prints a long subcommand name on its own line, its help on the next; 3.12 puts both on one line
        assert any(l.strip() == cmd or l.strip().startswith(cmd + " ") for l in top.splitlines()), cmd
