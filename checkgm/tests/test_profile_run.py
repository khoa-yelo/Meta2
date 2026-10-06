"""checkgm.profile_run: the orchestration of the pipeline containers, and above all its refusals.

Nothing here runs a real pipeline. A stub standing in for apptainer records the argv it was called with and writes the
files the real run scripts write, which is enough to pin the invocation (bind mounts, container-side paths, the disabled
scoring step of the read-based script) and the contract the CLI depends on: every user-correctable problem arrives as a
single-line ValueError or FileNotFoundError, a failing pipeline included.
"""
import os
import stat
import subprocess
import sys

import pytest

from checkgm import profile_run as P

PROFILE_HEAD = "#mpa_v30_CHOCOPhlAn_201901\n#45000000 reads processed\n#clade_name\tNCBI_tax_id\trelative_abundance\n"


def fake_dbs(root, pipeline):
    """A database directory holding every file the pipeline's run script opens, so preflight is satisfied."""
    d = os.path.join(str(root), f"dbs_{pipeline}")
    for rel in P.REQUIRED_DBS[pipeline]:
        p = os.path.join(d, rel)
        os.makedirs(os.path.dirname(p), exist_ok=True)
        open(p, "w").write("x")
    return d


def fake_image(root, pipeline):
    p = os.path.join(str(root), P.IMAGE_NAME[pipeline])
    open(p, "w").write("not really a sif")
    return p


def fake_runtime(root, body="", rc=0, stderr=""):
    """An executable stub for apptainer: it logs its argv to argv.txt and then does whatever `body` says."""
    p = os.path.join(str(root), "fake_apptainer")
    # its own PATH, because the fixture empties the one the module (and so the child) runs with
    open(p, "w").write("#!/bin/sh\nPATH=/usr/bin:/bin\nfor a in \"$@\"; do printf '%s\\n' \"$a\" >> \"$ARGV_LOG\"; done\n"
                       + (f"printf '%s\\n' {stderr!r} >&2\n" if stderr else "") + body + f"\nexit {rc}\n")
    os.chmod(p, os.stat(p).st_mode | stat.S_IEXEC)
    return p


def reads(root, *names):
    out = []
    for n in names:
        p = os.path.join(str(root), n)
        os.makedirs(os.path.dirname(p), exist_ok=True)
        open(p, "w").write("@r\nACGT\n+\nIIII\n")
        out.append(p)
    return out


@pytest.fixture(autouse=True)
def clean_env(monkeypatch, tmp_path):
    """No runtime, no image search path and no database search path unless a test asks for one."""
    for v in (P.IMAGES_ENV_VAR, P.DBS_ENV_VAR, P.RUNTIME_ENV_VAR):
        monkeypatch.delenv(v, raising=False)
    monkeypatch.setenv("PATH", str(tmp_path / "empty_path"))
    monkeypatch.setenv("ARGV_LOG", str(tmp_path / "argv.txt"))


def argv_of(tmp_path):
    return open(tmp_path / "argv.txt").read().splitlines()


# --------------------------------------------------------------------------------------------- preflight

def test_preflight_names_all_three_missing_pieces(tmp_path):
    r = P.preflight("assembly", None, None)
    assert len(r) == 3
    assert any(P.IMAGE_NAME["assembly"] in x for x in r)
    assert any("checkgm setup db" in x for x in r)
    assert any("apptainer" in x for x in r)
    assert all("\n" not in x for x in r)


def test_preflight_clean_when_everything_is_there(tmp_path, monkeypatch):
    monkeypatch.setenv(P.RUNTIME_ENV_VAR, fake_runtime(tmp_path))
    assert P.preflight("read", fake_dbs(tmp_path, "read"), fake_image(tmp_path, "read")) == []


def test_preflight_finds_image_and_dbs_through_the_environment(tmp_path, monkeypatch):
    monkeypatch.setenv(P.RUNTIME_ENV_VAR, fake_runtime(tmp_path))
    fake_image(tmp_path, "read")
    monkeypatch.setenv(P.IMAGES_ENV_VAR, str(tmp_path))
    dbs = fake_dbs(tmp_path, "read")
    monkeypatch.setenv(P.DBS_ENV_VAR, dbs)
    assert P.preflight("read", None, None) == []
    assert P.resolve_dbs("read", None)[0] == dbs


def test_half_downloaded_dbs_is_a_different_reason_than_a_missing_one(tmp_path):
    d = fake_dbs(tmp_path, "assembly")
    os.remove(os.path.join(d, "eggnog.db"))
    os.remove(os.path.join(d, "db_kofam.hmm.h3f"))
    why = P.resolve_dbs("assembly", d)[1]
    assert "incomplete" in why and "eggnog.db" in why and "2 of" in why and "checkgm setup db" in why
    assert P.resolve_dbs("assembly", d)[0] is None


def test_a_directory_given_as_the_image_is_said_to_be_one(tmp_path):
    assert "is a directory" in P.resolve_image("read", str(tmp_path))[1]


def test_stale_runtime_variable_reads_as_no_runtime(tmp_path, monkeypatch):
    monkeypatch.setenv(P.RUNTIME_ENV_VAR, str(tmp_path / "gone" / "apptainer"))
    assert P.runtime() is None
    assert any("not executable" in x for x in P.preflight("read", None, None))


def test_unknown_pipeline_is_one_reason(tmp_path):
    assert len(P.preflight("metaphlan", None, None)) == 1


# --------------------------------------------------------------------------------------------- refusals from run()

def one_line(exc):
    msg = str(exc.value)
    assert "\n" not in msg and msg
    return msg


def test_missing_pieces_arrive_as_one_line(tmp_path):
    r1, = reads(tmp_path, "S_1.fq.gz")
    with pytest.raises(ValueError) as e:
        P.run("assembly", "S", str(tmp_path / "out"), r1=r1)
    msg = one_line(e)
    assert "checkgm setup db" in msg and P.IMAGE_NAME["assembly"] in msg
    assert not os.path.exists(tmp_path / "out")   # a refusal creates nothing


def test_missing_read_file_is_reported_by_its_flag(tmp_path):
    with pytest.raises(FileNotFoundError) as e:
        P.run("read", "S", str(tmp_path / "out"), r1=str(tmp_path / "nope.fq.gz"))
    assert "--r1 file not found" in one_line(e)


def test_no_input_at_all_names_the_alternatives(tmp_path):
    with pytest.raises(ValueError) as e:
        P.run("read", "S", str(tmp_path / "out"))
    assert "--profile-in" in one_line(e)
    with pytest.raises(ValueError) as e:
        P.run("assembly", "S", str(tmp_path / "out"))
    assert "--contigs" in one_line(e)


def test_input_belonging_to_the_other_pipeline_says_which(tmp_path):
    f, = reads(tmp_path, "x.fa")
    with pytest.raises(ValueError) as e:
        P.run("read", "S", str(tmp_path / "out"), contigs=f)
    assert "--pipeline assembly" in one_line(e)
    with pytest.raises(ValueError) as e:
        P.run("assembly", "S", str(tmp_path / "out"), profile_in=f)
    assert "--pipeline read" in one_line(e)


def test_mate_without_first_mate_and_reads_with_contigs(tmp_path):
    r1, r2 = reads(tmp_path, "S_1.fq.gz", "S_2.fq.gz")
    with pytest.raises(ValueError) as e:
        P.run("read", "S", str(tmp_path / "out"), r2=r2)
    assert "--r2 was given without --r1" in one_line(e)
    with pytest.raises(ValueError) as e:
        P.run("assembly", "S", str(tmp_path / "out"), r1=r1, contigs=r2)
    assert "alternatives" in one_line(e)


@pytest.mark.parametrize("sample", ["", "a/b", "with space", "-leading"])
def test_sample_ids_that_cannot_be_a_file_name_are_refused(tmp_path, sample):
    with pytest.raises(ValueError) as e:
        P.run("read", sample, str(tmp_path / "out"), r1=reads(tmp_path, "S_1.fq.gz")[0])
    assert "file name" in one_line(e)


def test_threads_must_be_a_positive_integer(tmp_path):
    r1, = reads(tmp_path, "S_1.fq.gz")
    with pytest.raises(ValueError) as e:
        P.run("read", "S", str(tmp_path / "out"), r1=r1, threads=0)
    assert "--threads" in one_line(e)


def test_unknown_pipeline_from_run(tmp_path):
    with pytest.raises(ValueError) as e:
        P.run("humann", "S", str(tmp_path / "out"))
    assert "assembly" in one_line(e)


# --------------------------------------------------------------------------------------------- the invocation

def test_bind_plan_reuses_the_mounts_it_already_has(tmp_path):
    out = str(tmp_path / "out")
    dbs = fake_dbs(tmp_path, "read")
    a, b = reads(tmp_path, "libA/S_1.fq.gz", "libB/S_2.fq.gz")
    inside, = reads(tmp_path, "out/S_contigs.fa")
    binds, paths = P.bind_plan(out, dbs, [a, b, inside, a])
    assert binds[0] == f"{os.path.realpath(out)}:/out" and binds[1] == f"{os.path.realpath(dbs)}:/dbs:ro"
    assert paths[a] == "/in0/S_1.fq.gz" and paths[b] == "/in1/S_2.fq.gz"
    assert paths[inside] == "/out/S_contigs.fa"          # already under --out: addressed, not mounted again
    assert len(binds) == 4 and sum(b.endswith(":ro") for b in binds) == 3


def test_two_reads_in_one_directory_share_one_mount(tmp_path):
    r1, r2 = reads(tmp_path, "lib/S_1.fq.gz", "lib/S_2.fq.gz")
    binds, paths = P.bind_plan(str(tmp_path / "out"), fake_dbs(tmp_path, "read"), [r1, r2])
    assert len(binds) == 3 and paths[r1] == "/in0/S_1.fq.gz" and paths[r2] == "/in0/S_2.fq.gz"


def test_assembly_command_passes_the_pinned_arguments(tmp_path):
    r1, r2 = reads(tmp_path, "lib/S_1.fq.gz", "lib/S_2.fq.gz")
    argv = P.build_command("assembly", "S1", str(tmp_path / "out"), fake_dbs(tmp_path, "assembly"),
                           fake_image(tmp_path, "assembly"), {"r1": r1, "r2": r2}, 16, 120, exe="apptainer")
    assert argv[:3] == ["apptainer", "run", "--cleanenv"]
    assert argv[argv.index("--sample") + 1] == "S1"
    assert argv[argv.index("--out") + 1] == "/out" and argv[argv.index("--dbs") + 1] == "/dbs"
    assert argv[argv.index("--r1") + 1] == "/in0/S_1.fq.gz" and argv[argv.index("--r2") + 1] == "/in0/S_2.fq.gz"
    assert argv[argv.index("--mem-gb") + 1] == "120" and argv[argv.index("--threads") + 1] == "16"
    assert "--bundle" not in argv


def test_read_command_disables_the_scoring_step_of_the_run_script(tmp_path):
    r1, = reads(tmp_path, "lib/S_1.fq.gz")
    argv = P.build_command("read", "S1", str(tmp_path / "out"), fake_dbs(tmp_path, "read"),
                           fake_image(tmp_path, "read"), {"r1": r1}, 8, None)
    # run_pipelineB.sh insists on a bundle and then scores; placement is `checkgm assess`'s job, so its checkgm step is
    # turned into a no-op through the script's own CHECKGM_BIN hook and the bundle argument is a visible placeholder.
    assert argv[argv.index("--env") + 1] == "CHECKGM_BIN=/bin/true"
    assert argv[argv.index("--bundle") + 1] == "not-used-by-checkgm-profile"
    assert "--mem-gb" not in argv


def test_dry_run_prints_the_command_and_changes_nothing(tmp_path, capsys, monkeypatch):
    monkeypatch.setenv(P.RUNTIME_ENV_VAR, fake_runtime(tmp_path))
    r1, = reads(tmp_path, "lib/S_1.fq.gz")
    out = str(tmp_path / "out")
    r = P.run("assembly", "S1", out, r1=r1, dbs=fake_dbs(tmp_path, "assembly"), image=fake_image(tmp_path, "assembly"),
              threads=16, mem_gb=120, dry_run=True)
    printed = capsys.readouterr().out.strip()
    assert printed == " ".join(r["command"]) and "--bind" in printed
    assert r["status"] == "planned" and r["blockers"] == [] and r["assess_input"] is None
    assert not os.path.exists(out) and not os.path.exists(tmp_path / "argv.txt")


def test_dry_run_still_prints_the_command_when_nothing_is_installed(tmp_path, capsys):
    r1, = reads(tmp_path, "lib/S_1.fq.gz")
    r = P.run("read", "S1", str(tmp_path / "out"), r1=r1, dry_run=True)
    cap = capsys.readouterr()
    assert f"/path/to/{P.IMAGE_NAME['read']}" in cap.out and "/path/to/databases:/dbs:ro" in cap.out
    assert len(r["blockers"]) == 3
    assert all(l.startswith("checkgm: ") for l in cap.err.splitlines() if l.strip())


def test_thin_assembly_request_is_a_note_not_a_refusal(tmp_path, capsys, monkeypatch):
    monkeypatch.setenv(P.RUNTIME_ENV_VAR, fake_runtime(tmp_path))
    r1, = reads(tmp_path, "lib/S_1.fq.gz")
    r = P.run("assembly", "S1", str(tmp_path / "out"), r1=r1, dbs=fake_dbs(tmp_path, "assembly"),
              image=fake_image(tmp_path, "assembly"), threads=2, mem_gb=16, dry_run=True)
    assert len(r["notes"]) == 2 and "16 CPUs" in r["notes"][0] and "120 GB" in r["notes"][1]
    assert all(l.startswith("checkgm: ") for l in capsys.readouterr().err.splitlines() if l.strip())


# --------------------------------------------------------------------------------------------- a run that happens

def test_read_pipeline_run_reports_what_the_script_wrote(tmp_path, monkeypatch):
    out = str(tmp_path / "out")
    body = f"mkdir -p {out!r}\nprintf '%s' {PROFILE_HEAD!r} > {out!r}/S1.txt\nprintf 'S1\\t100\\n' > {out!r}/S1.reads.tsv\n" \
           f"echo step > {out!r}/pipelineB.log\n"
    monkeypatch.setenv(P.RUNTIME_ENV_VAR, fake_runtime(tmp_path, body=body))
    r1, = reads(tmp_path, "lib/S_1.fq.gz")
    r = P.run("read", "S1", out, r1=r1, dbs=fake_dbs(tmp_path, "read"), image=fake_image(tmp_path, "read"), threads=4)
    assert r["status"] == "ok" and r["assess_input"] == os.path.join(out, "S1.txt")
    assert r["outputs"]["reads_tsv"].endswith("S1.reads.tsv") and os.path.isfile(r["outputs"]["reads_tsv"])
    assert r["log"].endswith("pipelineB.log") and r["elapsed_s"] >= 0 and r["notes"] == []
    assert argv_of(tmp_path)[:3] == ["run", "--cleanenv", "--bind"]


def test_assembly_pipeline_run_hands_back_the_query_directory(tmp_path, monkeypatch):
    out = str(tmp_path / "out")
    monkeypatch.setenv(P.RUNTIME_ENV_VAR, fake_runtime(tmp_path, body=f"mkdir -p {out!r}/query\n"))
    c, = reads(tmp_path, "asm/contigs.fasta")
    r = P.run("assembly", "S1", out, contigs=c, dbs=fake_dbs(tmp_path, "assembly"),
              image=fake_image(tmp_path, "assembly"), threads=16)
    assert r["status"] == "ok" and r["assess_input"] == os.path.join(out, "query")
    assert "/in0/contigs.fasta" in argv_of(tmp_path)


def test_a_failing_pipeline_is_one_line_with_the_logs_and_the_last_message(tmp_path, monkeypatch):
    monkeypatch.setenv(P.RUNTIME_ENV_VAR, fake_runtime(tmp_path, rc=3, stderr="FATAL: could not open /dbs/eggnog.db"))
    r1, = reads(tmp_path, "lib/S_1.fq.gz")
    out = str(tmp_path / "out")
    with pytest.raises(ValueError) as e:
        P.run("read", "S1", out, r1=r1, dbs=fake_dbs(tmp_path, "read"), image=fake_image(tmp_path, "read"))
    msg = one_line(e)
    assert "exit 3" in msg and P.STDERR_LOG in msg and "could not open /dbs/eggnog.db" in msg
    assert "FATAL" in open(os.path.join(out, P.STDERR_LOG)).read()   # the whole of it is kept, off our stderr


def test_success_without_the_expected_output_is_still_an_error(tmp_path, monkeypatch):
    monkeypatch.setenv(P.RUNTIME_ENV_VAR, fake_runtime(tmp_path))
    r1, = reads(tmp_path, "lib/S_1.fq.gz")
    with pytest.raises(ValueError) as e:
        P.run("read", "S1", str(tmp_path / "out"), r1=r1, dbs=fake_dbs(tmp_path, "read"), image=fake_image(tmp_path, "read"))
    assert "did not write" in one_line(e)


def test_the_pipelines_stderr_never_reaches_ours(tmp_path, monkeypatch, capfd):
    """The one-line contract is load-bearing: a tool's warnings go to the log file, not to our stderr."""
    out = str(tmp_path / "out")
    monkeypatch.setenv(P.RUNTIME_ENV_VAR, fake_runtime(tmp_path, body=f"mkdir -p {out!r}/query\n",
                                                       stderr="RuntimeWarning: invalid value encountered"))
    c, = reads(tmp_path, "asm/contigs.fasta")
    P.run("assembly", "S1", out, contigs=c, dbs=fake_dbs(tmp_path, "assembly"), image=fake_image(tmp_path, "assembly"), threads=16)
    cap = capfd.readouterr()
    assert "RuntimeWarning" not in cap.err and "RuntimeWarning" in open(os.path.join(out, P.STDERR_LOG)).read()


# --------------------------------------------------------------------------------------------- a profile already in hand

def test_profile_in_is_staged_without_image_dbs_or_runtime(tmp_path):
    src = tmp_path / "elsewhere" / "mine.txt"
    os.makedirs(src.parent)
    src.write_text(PROFILE_HEAD + "k__Bacteria|...|s__Escherichia_coli\t562\t1.0\n")
    out = str(tmp_path / "out")
    r = P.run("read", "S1", out, profile_in=str(src))
    assert r["status"] == "staged" and r["assess_input"] == os.path.join(out, "S1.txt")
    assert open(r["outputs"]["profile"]).read() == src.read_text() and r["notes"] == []


def test_profile_without_a_read_count_is_noted(tmp_path, capsys):
    src = tmp_path / "mine.txt"
    src.write_text("#clade_name\tNCBI_tax_id\trelative_abundance\nk__Bacteria|s__E\t562\t1.0\n")
    r = P.run("read", "S1", str(tmp_path / "out"), profile_in=str(src))
    assert len(r["notes"]) == 1 and "--reads-tsv" in r["notes"][0]
    assert capsys.readouterr().err.startswith("checkgm: ")


def test_a_fastq_passed_as_a_profile_is_refused(tmp_path):
    r1, = reads(tmp_path, "S_1.fq.gz")
    with pytest.raises(ValueError) as e:
        P.run("read", "S1", str(tmp_path / "out"), profile_in=r1)
    assert "not a MetaPhlAn profile" in one_line(e)


def test_staging_the_file_already_in_place_is_not_a_copy_onto_itself(tmp_path):
    out = tmp_path / "out"
    out.mkdir()
    (out / "S1.txt").write_text(PROFILE_HEAD)
    r = P.run("read", "S1", str(out), profile_in=str(out / "S1.txt"))
    assert r["command"] == [] and open(r["assess_input"]).read() == PROFILE_HEAD


def test_dry_run_staging_touches_nothing(tmp_path, capsys):
    src = tmp_path / "mine.txt"
    src.write_text(PROFILE_HEAD)
    r = P.run("read", "S1", str(tmp_path / "out"), profile_in=str(src), dry_run=True)
    assert r["status"] == "planned" and "cp" in capsys.readouterr().out
    assert not os.path.exists(tmp_path / "out")


# --------------------------------------------------------------------------------------------- the module itself

def test_module_imports_without_the_optional_world(tmp_path):
    """Importing it must not need numpy, pandas or a container runtime: `checkgm profile --help` imports it."""
    res = subprocess.run([sys.executable, "-c", "import checkgm.profile_run as p; print(p.PIPELINES)"],
                         capture_output=True, text=True)
    assert res.returncode == 0 and res.stderr == "" and "assembly" in res.stdout
