"""checkgm.setup_data: the download/verify/unpack contract, the message for a Zenodo record that is not public yet, and the
delegation to the pinned container/ fetch scripts.

The released archives are 55-65 MB and are not in the repository, so every bundle test registers a synthetic bundle in
setup_data.BUNDLES (one directory with a manifest, tarred the way the released archives are) and checks the mechanism
against it. The recorded figures of the real bundles and databases are asserted separately, against the files they were
read from. HTTP is exercised against a local server that honours Range, including the record lookup, so the "it starts
working once the deposit is published" claim is tested rather than asserted.
"""
import hashlib
import http.server
import json
import os
import stat
import tarfile
import threading
import urllib.parse

import pytest

from checkgm import setup_data


# --------------------------------------------------------------------------------------------- helpers

def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for c in iter(lambda: f.read(1 << 16), b""):
            h.update(c)
    return h.hexdigest()


def make_archive(tmp_path, name="gut-reads-test-v9.9-lenient", into=None, payload="x"):
    """A bundle archive shaped like the released ones: one top-level directory named by the bundle id."""
    root = tmp_path / "build" / name
    (root / "features").mkdir(parents=True, exist_ok=True)
    (root / "manifest.json").write_text(json.dumps({"bundle_id": name, "profile_type": "reads", "percentile_grid": [1, 50, 99],
                                                    "layers": {}, "payload": payload}))
    (root / "features" / "taxonomy_species.parquet").write_text(payload * 100)
    out = (into or (tmp_path / "src"))
    out.mkdir(parents=True, exist_ok=True)
    arch = out / (name + ".tar.gz")
    with tarfile.open(arch, "w:gz") as t:
        t.add(str(root), arcname=name)
    return str(arch)


def register(monkeypatch, archive, name="gut-reads-test-v9.9-lenient", **over):
    """Make `archive` the recorded release of bundle `name`, and the only default bundle."""
    entry = {"bytes": os.path.getsize(archive), "sha256": sha256(archive), "unpacked_bytes": 4096,
             "profile_type": "reads", "tier": "lenient", "deposited": True}
    entry.update(over)
    monkeypatch.setitem(setup_data.BUNDLES, name, entry)
    monkeypatch.setattr(setup_data, "DEFAULT_BUNDLES", (name,))
    return name


class _RangeHandler(http.server.BaseHTTPRequestHandler):
    """Static files with Range support, which http.server's own handler does not have (the resume path needs a 206)."""
    root = ""
    seen_ranges: list = []

    def do_GET(self):
        path = os.path.join(self.root, os.path.basename(urllib.parse.urlparse(self.path).path))
        if not os.path.isfile(path):
            self.send_error(404, "not here"); return
        data = open(path, "rb").read()
        rng = self.headers.get("Range")
        if rng and rng.startswith("bytes="):
            type(self).seen_ranges.append(rng)
            start = int(rng.split("=", 1)[1].split("-", 1)[0])
            body, code = data[start:], 206
        else:
            body, code = data, 200
        self.send_response(code)
        if code == 206:
            self.send_header("Content-Range", f"bytes {len(data) - len(body)}-{len(data) - 1}/{len(data)}")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def log_message(self, *a):
        pass


@pytest.fixture
def server(tmp_path):
    """Serves tmp_path/'src' at the yielded base URL."""
    root = tmp_path / "src"
    root.mkdir(parents=True, exist_ok=True)
    handler = type("H", (_RangeHandler,), {"root": str(root), "seen_ranges": []})
    httpd = http.server.ThreadingHTTPServer(("127.0.0.1", 0), handler)
    threading.Thread(target=httpd.serve_forever, daemon=True).start()
    yield type("S", (), {"base": f"http://127.0.0.1:{httpd.server_address[1]}/", "root": root, "handler": handler})()
    httpd.shutdown()


def stub_container(tmp_path, body, name="fetch_dbs_read.sh"):
    """A container/ directory whose fetch script is `body`. fetch_dbs_assembly.sh must exist for container_dir() to accept it."""
    c = tmp_path / "container"
    c.mkdir(exist_ok=True)
    for n in ("fetch_dbs_assembly.sh", "fetch_dbs_read.sh", "unpack_dbs.sh"):
        (c / n).write_text("#!/usr/bin/env bash\nexit 0\n")
    (c / name).write_text(body)
    (c / name).chmod((c / name).stat().st_mode & ~stat.S_IXUSR)   # the real fetch_dbs_read.sh carries no executable bit
    return str(c)


# --------------------------------------------------------------------------------------------- recorded figures

def test_recorded_bundle_figures_match_the_release():
    """The four current bundles, with the bytes and sha256 of release/RELEASE.json and the README download table."""
    assert set(setup_data.BUNDLES) == {"gut-assembly-adult-global-v0.8-lenient", "gut-assembly-adult-global-v0.8-strict",
                                       "gut-reads-adult-global-v0.2-lenient", "gut-reads-adult-global-v0.2-strict"}
    b = setup_data.BUNDLES["gut-reads-adult-global-v0.2-lenient"]
    assert (b["bytes"], b["sha256"]) == (64799615, "ab6b82bdc62b72dcfa860ea2dbf0e7d6e750c8058723721c3b2e7f25415b5c0e")
    a = setup_data.BUNDLES["gut-assembly-adult-global-v0.8-lenient"]
    assert (a["bytes"], a["sha256"]) == (54519283, "0d31c959dad8ff7607afc0f2baaf04bc2d7e419cdf1dc5617d4ad3eb229d90f8")
    assert all(len(v["sha256"]) == 64 for v in setup_data.BUNDLES.values())
    # Only the lenient tier is deposited; the strict tier is rebuilt, so it cannot be a default.
    assert setup_data.DEFAULT_BUNDLES == ("gut-assembly-adult-global-v0.8-lenient", "gut-reads-adult-global-v0.2-lenient")


def test_plan_without_a_pipeline_is_the_bundles_alone():
    rows = setup_data.plan()
    assert [r["item"] for r in rows] == list(setup_data.DEFAULT_BUNDLES)
    assert all(r["kind"] == "bundle" and r["n_files"] == 1 for r in rows)
    assert sum(r["download_bytes"] for r in rows) == 54519283 + 64799615
    assert setup_data.DOI in rows[0]["source"] and "not public" in rows[0]["note"]


def test_plan_states_the_measured_database_sizes():
    """The point of plan(): the assembly databases are two orders of magnitude bigger than the read one, and the caller
    is told so before anything starts."""
    asm = [r for r in setup_data.plan("assembly") if r["kind"] == "databases"]
    read = [r for r in setup_data.plan("read") if r["kind"] == "databases"]
    assert len(asm) == len(read) == 1
    assert asm[0]["download_bytes"] == 90437315335 and asm[0]["disk_bytes"] == 170047067806 and asm[0]["n_files"] == 24
    assert read[0]["download_bytes"] == 384430080 and read[0]["disk_bytes"] == 2936494131
    assert asm[0]["item"] == "fetch_dbs_assembly.sh" and read[0]["item"] == "fetch_dbs_read.sh"
    assert setup_data.plan("reads")[-1]["pipeline"] == "read"   # the bundles' word for the same pipeline
    assert setup_data.human_bytes(asm[0]["download_bytes"]) == "90.4 GB"


def test_plan_rejects_an_unknown_pipeline():
    with pytest.raises(ValueError) as e:
        setup_data.plan("tierC")
    assert "tierC" in str(e.value) and "assembly" in str(e.value) and "\n" not in str(e.value)


def test_select_takes_ids_and_profile_types():
    assert setup_data.select(None) == list(setup_data.DEFAULT_BUNDLES)
    assert setup_data.select("gut-reads-adult-global-v0.2-strict") == ["gut-reads-adult-global-v0.2-strict"]
    assert setup_data.select("reads") == setup_data.select("read") == ["gut-reads-adult-global-v0.2-lenient"]
    assert setup_data.select(["assembly", "assembly"]) == ["gut-assembly-adult-global-v0.8-lenient"]
    with pytest.raises(ValueError) as e:
        setup_data.select("gut-reads-adult-global-v0.2-lenieint")
    assert "unknown bundle" in str(e.value) and "gut-reads-adult-global-v0.2-lenient" in str(e.value)


# --------------------------------------------------------------------------------------------- bundles()

def test_dry_run_transfers_nothing_and_needs_no_network(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(setup_data, "ZENODO_RECORD_API", "http://127.0.0.1:1/{rid}")   # any call would fail
    dest = tmp_path / "bundles"
    r = setup_data.bundles(str(dest), dry_run=True)
    assert r["dry_run"] and not dest.exists()
    assert r["total_bytes"] == 54519283 + 64799615
    assert all(row["action"] == "would download" and row["verified"] is False for row in r["bundles"])
    assert capsys.readouterr() == ("", "")


def test_local_directory_is_verified_unpacked_and_idempotent(tmp_path, monkeypatch, capsys):
    arch = make_archive(tmp_path)
    name = register(monkeypatch, arch)
    dest = tmp_path / "bundles"

    r = setup_data.bundles(str(dest), url=str(tmp_path / "src"))
    row, = r["bundles"]
    assert (row["name"], row["verified"], row["action"]) == (name, True, "downloaded")
    assert row["sha256"] == sha256(arch) and row["bytes"] == os.path.getsize(arch)
    assert os.path.isfile(os.path.join(row["bundle_dir"], "manifest.json"))
    assert r["downloaded_bytes"] == os.path.getsize(arch)
    out, err = capsys.readouterr()
    assert err == "" and "CHECKGM_BUNDLES" in out   # progress on stdout only; stderr is the error channel

    again = setup_data.bundles(str(dest), url=str(tmp_path / "src"))
    assert again["bundles"][0]["action"] == "skipped" and again["bundles"][0]["verified"] is True
    assert again["downloaded_bytes"] == 0
    assert not os.path.exists(os.path.join(str(dest), name + ".tar.gz.part"))


def test_a_single_archive_path_is_accepted_for_one_bundle(tmp_path, monkeypatch):
    arch = make_archive(tmp_path)
    register(monkeypatch, arch)
    r = setup_data.bundles(str(tmp_path / "b"), url=arch)
    assert r["bundles"][0]["verified"] is True
    with pytest.raises(ValueError) as e:
        setup_data.bundles(str(tmp_path / "b2"), url=arch, only=list(setup_data.BUNDLES)[:2])
    assert "one archive" in str(e.value)


def test_unpacked_bundle_without_its_archive_is_left_alone(tmp_path, monkeypatch):
    arch = make_archive(tmp_path)
    name = register(monkeypatch, arch)
    dest = tmp_path / "bundles"
    setup_data.bundles(str(dest), url=str(tmp_path / "src"))
    os.remove(dest / (name + ".tar.gz"))
    r = setup_data.bundles(str(dest), url=str(tmp_path / "src"))
    assert r["bundles"][0]["verified"] is False and "nothing to verify" in r["bundles"][0]["action"]
    assert r["downloaded_bytes"] == 0


def test_a_corrupt_archive_is_refused_and_not_kept(tmp_path, monkeypatch):
    good = make_archive(tmp_path, into=tmp_path / "good")
    name = register(monkeypatch, good)
    make_archive(tmp_path, payload="y")   # same name under tmp_path/src, different content
    dest = tmp_path / "bundles"
    with pytest.raises(ValueError) as e:
        setup_data.bundles(str(dest), url=str(tmp_path / "src"))
    msg = str(e.value)
    assert "failed verification" in msg and sha256(good) in msg and "\n" not in msg
    assert not os.path.exists(dest / (name + ".tar.gz")) and not os.path.exists(dest / name)


def test_an_archive_that_is_not_a_bundle_is_refused(tmp_path, monkeypatch):
    stray = tmp_path / "src"
    stray.mkdir(parents=True)
    arch = stray / "gut-reads-test-v9.9-lenient.tar.gz"
    (tmp_path / "elsewhere.txt").write_text("hi")
    with tarfile.open(arch, "w:gz") as t:
        t.add(str(tmp_path / "elsewhere.txt"), arcname="../elsewhere.txt")
    register(monkeypatch, str(arch))
    with pytest.raises(ValueError) as e:
        setup_data.bundles(str(tmp_path / "bundles"), url=str(stray))
    assert "not a checkGM bundle archive" in str(e.value)


def test_http_download_and_resume_of_a_part_file(tmp_path, monkeypatch, server):
    arch = make_archive(tmp_path)
    name = register(monkeypatch, arch)
    dest = tmp_path / "bundles"
    dest.mkdir()
    (dest / (name + ".tar.gz.part")).write_bytes(open(arch, "rb").read()[:64])

    r = setup_data.bundles(str(dest), url=server.base)
    assert r["bundles"][0]["verified"] is True and r["source"] == server.base
    assert server.handler.seen_ranges == ["bytes=64-"]
    assert sha256(str(dest / (name + ".tar.gz"))) == sha256(arch)


def test_http_failure_is_one_line(tmp_path, monkeypatch, server):
    arch = make_archive(tmp_path, into=tmp_path / "elsewhere")
    register(monkeypatch, arch)
    with pytest.raises(ValueError) as e:
        setup_data.bundles(str(tmp_path / "bundles"), url=server.base)   # the archive is not under the served root
    assert "HTTP 404" in str(e.value) and "\n" not in str(e.value)


# --------------------------------------------------------------------------------------------- the Zenodo record

def test_the_reserved_doi_fails_with_an_actionable_line(tmp_path, monkeypatch, server):
    """Today's state: the record does not exist, so the lookup 404s."""
    register(monkeypatch, make_archive(tmp_path))
    monkeypatch.setattr(setup_data, "ZENODO_RECORD_API", server.base + "absent_{rid}.json")
    with pytest.raises(ValueError) as e:
        setup_data.bundles(str(tmp_path / "bundles"))
    msg = str(e.value)
    assert "\n" not in msg
    assert setup_data.DOI in msg and "not public yet" in msg and "--url" in msg
    assert "gut-reads-test-v9.9-lenient.tar.gz" in msg   # names the file a request should ask for


def test_the_same_call_works_once_the_record_is_published(tmp_path, monkeypatch, server):
    """The published record is reached through the real records API shape; nothing else about the call changes."""
    arch = make_archive(tmp_path)
    name = register(monkeypatch, arch)
    rid = "23181418"
    (server.root / f"record_{rid}.json").write_text(json.dumps(
        {"files": [{"key": name + ".tar.gz", "size": os.path.getsize(arch),
                    "links": {"self": server.base + "meta", "content": server.base + name + ".tar.gz"}}]}))
    monkeypatch.setattr(setup_data, "ZENODO_RECORD_API", server.base + "record_{rid}.json")

    r = setup_data.bundles(str(tmp_path / "bundles"))
    assert r["source"] == f"zenodo doi {setup_data.DOI}"
    assert r["bundles"][0]["verified"] is True


def test_a_record_whose_file_is_a_different_build_is_refused_before_downloading(tmp_path, monkeypatch, server):
    arch = make_archive(tmp_path)
    name = register(monkeypatch, arch)
    (server.root / "record_23181418.json").write_text(json.dumps(
        {"files": [{"key": name + ".tar.gz", "size": os.path.getsize(arch) + 1,
                    "links": {"content": server.base + name + ".tar.gz"}}]}))
    monkeypatch.setattr(setup_data, "ZENODO_RECORD_API", server.base + "record_{rid}.json")
    with pytest.raises(ValueError) as e:
        setup_data.bundles(str(tmp_path / "bundles"))
    assert "different build" in str(e.value) and not (tmp_path / "bundles" / (name + ".tar.gz")).exists()


def test_a_record_missing_the_strict_tier_says_it_is_rebuilt(tmp_path, monkeypatch, server):
    (server.root / "record_23181418.json").write_text(json.dumps(
        {"files": [{"key": "gut-reads-adult-global-v0.2-lenient.tar.gz", "links": {"content": server.base + "x"}}]}))
    monkeypatch.setattr(setup_data, "ZENODO_RECORD_API", server.base + "record_{rid}.json")
    with pytest.raises(ValueError) as e:
        setup_data.bundles(str(tmp_path / "bundles"), only="gut-reads-adult-global-v0.2-strict")
    assert "not deposited" in str(e.value) and "docs/rebuild.md" in str(e.value)


def test_record_id_comes_from_the_doi():
    assert setup_data._record_id(setup_data.DOI) == "23181418"
    assert setup_data._record_id("https://doi.org/10.5281/zenodo.23181418") == "23181418"
    assert setup_data._record_id("23181418") == "23181418"
    with pytest.raises(ValueError):
        setup_data._record_id("10.1093/nar/gkac1080")


# --------------------------------------------------------------------------------------------- databases()

def test_databases_dry_run_names_the_pinned_scripts(tmp_path, monkeypatch):
    monkeypatch.setenv(setup_data.CONTAINER_ENV, stub_container(tmp_path, "exit 0"))
    r = setup_data.databases("assembly", str(tmp_path / "dbs"), dry_run=True)
    assert [os.path.basename(s["script"]) for s in r["steps"]] == ["fetch_dbs_assembly.sh", "unpack_dbs.sh"]
    assert r["steps"][0]["argv"][:1] == ["bash"] and r["steps"][0]["argv"][-1] == str(tmp_path / "dbs")
    assert r["download_bytes"] == 90437315335 and "missing_tools" in r
    assert not (tmp_path / "dbs").exists() and r["steps"][0]["returncode"] is None
    assert [os.path.basename(s["script"]) for s in setup_data.databases("read", str(tmp_path / "d2"), dry_run=True)["steps"]] \
        == ["fetch_dbs_read.sh"]


def test_databases_delegates_with_dest_as_the_argument(tmp_path, monkeypatch, capsys):
    monkeypatch.setenv(setup_data.CONTAINER_ENV,
                       stub_container(tmp_path, '#!/usr/bin/env bash\necho "[ok] marker"\necho "$1" > "$1/ran.txt"\n'))
    monkeypatch.setattr(setup_data, "_missing_tools", lambda p: [])
    dest = tmp_path / "dbs"
    r = setup_data.databases("read", str(dest), dry_run=False)
    assert r["ok"] is True and r["steps"][0]["returncode"] == 0
    assert (dest / "ran.txt").read_text().strip() == str(dest)
    out, err = capsys.readouterr()
    assert "[ok] marker" in out and err == ""   # the script's output is relayed as progress, on stdout


def test_databases_reports_files_the_script_failed_to_fetch(tmp_path, monkeypatch):
    monkeypatch.setenv(setup_data.CONTAINER_ENV,
                       stub_container(tmp_path, '#!/usr/bin/env bash\necho "[FAIL] eggnog.db"\nexit 0\n', "fetch_dbs_assembly.sh"))
    monkeypatch.setattr(setup_data, "_missing_tools", lambda p: [])
    with pytest.raises(ValueError) as e:
        setup_data.databases("assembly", str(tmp_path / "dbs"))
    msg = str(e.value)
    assert "eggnog.db" in msg and "resumable" in msg and "\n" not in msg


def test_databases_reports_a_nonzero_exit(tmp_path, monkeypatch):
    monkeypatch.setenv(setup_data.CONTAINER_ENV, stub_container(tmp_path, '#!/usr/bin/env bash\nexit 3\n'))
    monkeypatch.setattr(setup_data, "_missing_tools", lambda p: [])
    with pytest.raises(ValueError) as e:
        setup_data.databases("read", str(tmp_path / "dbs"))
    assert "exited 3" in str(e.value) and "run the same command again" in str(e.value)


def test_databases_refuses_a_run_without_its_tools_but_still_plans_one(tmp_path, monkeypatch):
    monkeypatch.setenv(setup_data.CONTAINER_ENV, stub_container(tmp_path, "exit 0"))
    monkeypatch.setattr(setup_data, "_missing_tools", lambda p: ["bowtie2-build"])
    assert setup_data.databases("read", str(tmp_path / "dbs"), dry_run=True)["missing_tools"] == ["bowtie2-build"]
    with pytest.raises(ValueError) as e:
        setup_data.databases("read", str(tmp_path / "dbs"))
    assert "bowtie2-build" in str(e.value) and "checkgm_read.def" in str(e.value)


def test_missing_container_directory_says_where_to_get_it(tmp_path, monkeypatch):
    monkeypatch.delenv(setup_data.CONTAINER_ENV, raising=False)
    monkeypatch.setattr(setup_data, "__file__", str(tmp_path / "nowhere" / "src" / "checkgm" / "setup_data.py"))
    monkeypatch.chdir(tmp_path)
    with pytest.raises(FileNotFoundError) as e:
        setup_data.databases("read", str(tmp_path / "dbs"), dry_run=True)
    assert setup_data.CONTAINER_ENV in str(e.value) and "container/" in str(e.value)


def test_the_repository_ships_the_scripts_databases_delegates_to(monkeypatch):
    """Skipped when the tests run against an installed wheel, which ships no container/."""
    monkeypatch.delenv(setup_data.CONTAINER_ENV, raising=False)
    try:
        cdir = setup_data.container_dir()
    except FileNotFoundError:
        pytest.skip("container/ is not beside the package (wheel install)")
    for d in setup_data.DATABASES.values():
        for s in (d["script"], d["then"]):
            assert s is None or os.path.isfile(os.path.join(cdir, s))
