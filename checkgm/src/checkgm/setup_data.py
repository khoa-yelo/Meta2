"""checkgm.setup_data — fetching the two things a run needs: a reference baseline bundle, and, only for a pipeline run from
raw reads, that pipeline's reference databases.

The two halves are deliberately unlike each other. A bundle is small (55-65 MB), is distributed as one archive with a
recorded sha256, and is fetched here: `bundles()` downloads, verifies and unpacks it, so that the directory it was given
works as ``$CHECKGM_BUNDLES`` (paths.py looks for the directory holding manifest.json, which is what the archive unpacks
to). The pipeline databases are four orders of magnitude larger and their provenance is a pinned list of EBI and cibio
URLs that the container recipes were validated against; reimplementing that list here would let it drift from the
pipeline, so `databases()` delegates to container/fetch_dbs_assembly.sh and container/fetch_dbs_read.sh and only reports what they did.

Nothing large starts without the caller having been told the size: `plan()` answers "what would this download, and how
big is it" from recorded figures alone, without touching the network, and every function takes ``dry_run``.

The Zenodo record that will host the bundles is not published yet. Its DOI is reserved (``DOI`` below, the one printed in
the paper) and resolves to nothing today, so a call with no explicit location fails with one actionable line rather than
an HTTP traceback, and the same call starts working, unchanged, the day the deposit is published.

User errors are raised as ValueError/FileNotFoundError/OSError with a single-sentence message; cli.py turns those into the
one "checkgm: ..." line on stderr and exit 1.
"""
from __future__ import annotations

import contextlib
import hashlib
import json
import os
import re
import shutil
import subprocess
import tarfile
import urllib.error
import urllib.parse
import urllib.request
import warnings

from checkgm import REPO

# Bundle archives: bytes and sha256 from release/RELEASE.json, re-verified with sha256sum against the archives themselves
# and against the download table in README.md. `deposited` marks the two the Zenodo record carries (release/zenodo.json
# "files"): the strict tier is a sensitivity analysis, is not deposited, and is rebuilt from the archived configuration
# (docs/rebuild.md), so asking for it needs an explicit location.
BUNDLES = {
    "gut-assembly-adult-global-v0.8-lenient": {"bytes": 54519283, "unpacked_bytes": 85777710, "profile_type": "assembly",
                                               "tier": "lenient", "deposited": True,
                                               "sha256": "0d31c959dad8ff7607afc0f2baaf04bc2d7e419cdf1dc5617d4ad3eb229d90f8"},
    "gut-assembly-adult-global-v0.8-strict": {"bytes": 54318945, "unpacked_bytes": 85564101, "profile_type": "assembly",
                                              "tier": "strict", "deposited": False,
                                              "sha256": "ddffa36865e3c971b93a826bf47694f72809f935f637a20955bdb9bab363bd5c"},
    "gut-reads-adult-global-v0.2-lenient": {"bytes": 64799615, "unpacked_bytes": 95112954, "profile_type": "reads",
                                            "tier": "lenient", "deposited": True,
                                            "sha256": "ab6b82bdc62b72dcfa860ea2dbf0e7d6e750c8058723721c3b2e7f25415b5c0e"},
    "gut-reads-adult-global-v0.2-strict": {"bytes": 64771253, "unpacked_bytes": 95102595, "profile_type": "reads",
                                           "tier": "strict", "deposited": False,
                                           "sha256": "4742e992125786a8286b030510dd6c45e727fd29620a7e6d4ae5588eaefb8734"},
}
ARCHIVE_SUFFIX = ".tar.gz"
# The default is the two lenient baselines: the tier the paper uses, the only tier deposited, and one bundle per pipeline.
DEFAULT_BUNDLES = tuple(n for n, b in BUNDLES.items() if b["deposited"])

DOI = "10.5281/zenodo.23181418"   # reserved by the draft deposit; printed in paper/main.tex; resolves only once published
ZENODO_RECORD_API = "https://zenodo.org/api/records/{rid}"
HTTP_TIMEOUT = 60                 # applies to the record lookup and to each read of a download, not to a whole transfer
CHUNK = 1 << 20

CONTAINER_ENV = "CHECKGM_CONTAINER"   # where container/ was unpacked, for an install from the wheel (which ships no scripts)

# Database figures measured on the staged copy the baselines were built with (project/staging/dbs and dbs_mpa):
# `download_bytes` sums the files the fetch script's URL list names, `disk_bytes` is what the directory holds after
# unpack_dbs.sh, which keeps both the archive and the gunzipped copy (gunzip -k). Both are larger than the round numbers
# that README.md recorded before 2026-10-06 ("about 70 GB compressed, 110 GB unpacked"), which predate the full eggNOG
# and InterProScan provisioning and have since been corrected there;
# a caller shown too small a figure is exactly the failure plan() exists to prevent, so the measured values are used.
DATABASES = {
    "assembly": {"script": "fetch_dbs_assembly.sh", "then": "unpack_dbs.sh", "n_files": 24,
                 "download_bytes": 90437315335, "disk_bytes": 170047067806,
                 "tools": ("bash", "curl", "xargs", "md5sum", "gunzip", "tar", "python3"),   # unpack_dbs.sh needs the last four
                 "what": "MGnify pipeline v5.0 reference databases (UniRef90 DIAMOND, eggNOG 2.0.0, KOfam, Pfam 32.0, "
                         "InterProScan 5.36-75.0, KEGG module graphs)",
                 "note": "downloaded once, resumable, 4 parallel streams; read QC against hg38 is off by default and its "
                         "index is not fetched"},
    "read": {"script": "fetch_dbs_read.sh", "then": None, "n_files": 1,
             "download_bytes": 384430080, "disk_bytes": 2936494131,
             "tools": ("bash", "curl", "bowtie2-build", "md5sum", "tar", "bunzip2"),   # bunzip2 is its own package on slim images
             "what": "MetaPhlAn 3 marker database mpa_v30_CHOCOPhlAn_201901 (the one curatedMetagenomicData 3 used)",
             "note": "the archive's md5 is pinned in the script; the Bowtie2 index is built locally afterwards, which is "
                     "most of the final size"},
}
PIPELINE_ALIASES = {"reads": "read"}   # bundles say profile_type 'reads'; the pipeline scripts and the CLI say 'read'
PROFILE_TYPE_ALIASES = {"read": "reads"}   # the same two words seen from the bundle side, for `only`


# --------------------------------------------------------------------------------------------- small helpers

def human_bytes(n: int) -> str:
    """A size in decimal units, as the README download table and `plan()` output quote them."""
    for unit, scale in (("TB", 1e12), ("GB", 1e9), ("MB", 1e6), ("kB", 1e3)):
        if n >= scale:
            return f"{n / scale:.1f} {unit}"
    return f"{n} B"


@contextlib.contextmanager
def _quiet():
    """Keep library warnings off stderr while urllib and tarfile are in use.

    The single-line "checkgm: ..." contract is tested, and both modules have emitted DeprecationWarnings on some
    interpreter versions; a user whose download worked must not see one."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        yield


def _say(msg: str):
    """One progress line on stdout, where cli.py puts per-sample progress so that stderr stays the error channel."""
    print(msg, flush=True)


def _sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(CHUNK), b""):
            h.update(chunk)
    return h.hexdigest()


def _is_url(s: str) -> bool:
    return bool(re.match(r"^[a-zA-Z][a-zA-Z0-9+.-]*://", s))


def select(only=None) -> list[str]:
    """The bundle ids a request names: the deposited pair by default, or what `only` asks for.

    `only` takes a bundle id, a profile type ('assembly' or 'reads', the manifest's own word), or a list of either, so
    that a CLI flag can say "both" or "just the read-based one" without a second vocabulary."""
    if only is None:
        return list(DEFAULT_BUNDLES)
    names = [only] if isinstance(only, str) else list(only)
    out: list[str] = []
    for n in names:
        if n in BUNDLES:
            picked = [n]
        else:
            want = PROFILE_TYPE_ALIASES.get(n, n)
            # A profile type means the deposited bundle of that type, not every tier of it: the strict tier is a
            # sensitivity analysis that was never deposited, so including it here would turn `--only assembly` into a
            # request that cannot be served from the record.
            picked = [b for b, m in BUNDLES.items() if m["profile_type"] == want and m["deposited"]]
            if not picked:
                raise ValueError(f"unknown bundle '{n}'; it is a bundle id ({', '.join(sorted(BUNDLES))}) or a profile type "
                                 f"(assembly, reads)")
        out.extend(p for p in picked if p not in out)
    return out


def _not_public_message(doi: str, names) -> str:
    """The message for a DOI that resolves to nothing — the state of the record today.

    It has to name the remedy, because the remedy is not "wait": the archives exist and are handed out by the authors,
    and `--url` takes either a directory URL or a local directory holding them."""
    return (f"the Zenodo record for the checkGM baselines (doi {doi}) is not public yet, so there is no download URL: pass "
            f"--url with a directory URL or a local directory holding {', '.join(n + ARCHIVE_SUFFIX for n in names)}, or "
            f"request the archives from the authors ({REPO}, README section 'Get a reference bundle'); this command needs no "
            f"change once the deposit is published")


# --------------------------------------------------------------------------------------------- where the archives are

def _record_id(doi: str) -> str:
    """The Zenodo record id inside a DOI, a doi.org URL, a record URL or a bare id.

    Zenodo's reserved DOI carries the record id as its suffix, which is what the records API is keyed by; nothing else
    about the DOI is needed here."""
    m = re.search(r"zenodo\.(\d+)", doi) or re.fullmatch(r"\s*(\d+)\s*", doi)
    if not m:
        raise ValueError(f"'{doi}' is not a Zenodo DOI or record id (expected something like {DOI})")
    return m.group(1)


def _zenodo_files(doi: str, names, timeout: int = HTTP_TIMEOUT) -> dict:
    """{file name: (download url, size or None)} for a published Zenodo record.

    A record that does not exist, or exists as an unpublished draft, answers 404 (410 once withdrawn); both mean the same
    thing to a user and get the one message that says what to do instead."""
    rid = _record_id(doi)
    url = ZENODO_RECORD_API.format(rid=rid)
    try:
        with _quiet(), urllib.request.urlopen(url, timeout=timeout) as r:
            rec = json.loads(r.read().decode("utf-8", "replace"))
    except urllib.error.HTTPError as e:
        if e.code in (404, 410):
            raise ValueError(_not_public_message(doi, names)) from None
        raise ValueError(f"zenodo.org answered HTTP {e.code} for record {rid} (doi {doi}); pass --url with a location you can "
                         f"reach, or retry later") from None
    except urllib.error.URLError as e:
        raise ValueError(f"cannot reach zenodo.org to resolve doi {doi} ({e.reason}); pass --url with a directory URL or a local "
                         f"directory holding the bundle archives") from None
    except json.JSONDecodeError:
        raise ValueError(f"zenodo.org did not answer with JSON for record {rid} (doi {doi}); retry, or pass --url") from None
    out = {}
    for f in rec.get("files") or []:
        key = f.get("key") or f.get("filename")
        links = f.get("links") or {}
        # Two shapes are in the wild for the same record: the current records API puts the download under links.content and
        # uses links.self for the file's metadata endpoint, the legacy one had the download under links.self.
        link = links.get("content") or links.get("download") or links.get("self")
        if key and link:
            out[key] = (link, f.get("size"))
    if not out:
        raise ValueError(f"the Zenodo record {rid} (doi {doi}) lists no downloadable files; pass --url with a location holding "
                         f"the bundle archives")
    return out


def _sources(names, url=None, doi=None) -> tuple[dict, str]:
    """({bundle id: (location, size or None)}, a label for the location), without transferring anything.

    `url` wins over `doi` because it is the only thing that works today: it is a local directory, a local archive (when
    one bundle was asked for), or a URL prefix the archive name is appended to. Without either, the reserved DOI is
    tried, which is the call that will start working when the deposit goes public."""
    if url:
        if _is_url(url):
            base = url if url.endswith("/") else url + "/"
            return {n: (urllib.parse.urljoin(base, n + ARCHIVE_SUFFIX), None) for n in names}, base
        if os.path.isfile(url):
            if len(names) != 1:
                raise ValueError(f"--url names one archive ({url}) but {len(names)} bundles were requested; give the directory "
                                 f"holding them, or ask for one bundle")
            return {names[0]: (os.path.abspath(url), os.path.getsize(url))}, os.path.abspath(url)
        if not os.path.isdir(url):
            raise FileNotFoundError(f"--url is neither a URL nor an existing path: {url}")
        return ({n: (os.path.join(os.path.abspath(url), n + ARCHIVE_SUFFIX), None) for n in names},
                os.path.abspath(url))
    doi = doi or DOI
    files = _zenodo_files(doi, names)
    out = {}
    for n in names:
        key = n + ARCHIVE_SUFFIX
        if key not in files:
            rebuilt = "" if BUNDLES[n]["deposited"] else (f" (the {BUNDLES[n]['tier']} tier is not deposited; it is rebuilt from "
                                                          f"the archived configuration, docs/rebuild.md)")
            raise ValueError(f"the Zenodo record for doi {doi} does not hold {key}{rebuilt}; it carries "
                             f"{', '.join(sorted(files))}")
        out[n] = files[key]
    return out, f"zenodo doi {doi}"


# --------------------------------------------------------------------------------------------- transfer and unpack

def _download(src: str, dest_file: str, expect: int):
    """Fetch `src` to `dest_file`, resuming a part-file from an earlier interrupted run where the server allows it.

    The transfer goes to `dest_file`.part and is renamed only once the caller has verified the digest, so an interrupted
    run never leaves something that looks like a complete archive. A server that ignores the Range header answers 200
    with the whole body, in which case the part-file is restarted rather than appended to, which would corrupt it."""
    part = dest_file + ".part"
    if not _is_url(src):
        if not os.path.isfile(src):
            raise FileNotFoundError(f"bundle archive not found: {src}")
        shutil.copyfile(src, part)
    else:
        have = os.path.getsize(part) if os.path.exists(part) else 0
        if have >= expect > 0:
            have = 0   # a part-file at or past the recorded size is from a different archive; start again
        req = urllib.request.Request(src, headers={"Range": f"bytes={have}-"} if have else {})
        try:
            with _quiet(), urllib.request.urlopen(req, timeout=HTTP_TIMEOUT) as r:
                mode = "ab" if have and r.status == 206 else "wb"
                if have and r.status != 206:
                    have = 0
                with open(part, mode) as f:
                    shutil.copyfileobj(r, f, CHUNK)
        except urllib.error.HTTPError as e:
            raise ValueError(f"download of {os.path.basename(dest_file)} failed: HTTP {e.code} for {src}") from None
        except urllib.error.URLError as e:
            raise ValueError(f"download of {os.path.basename(dest_file)} failed: {e.reason} ({src})") from None
    os.replace(part, dest_file)


def _unpack(archive: str, dest: str, name: str) -> str:
    """Unpack a verified bundle archive into `dest` and return the bundle directory it produced.

    Every member is required to sit under the bundle id, which is how the released archives are built; the check keeps an
    archive from scattering files across `dest` (or outside it) and doubles as a statement that what arrived is a bundle.
    The directory is built beside its final name and moved into place, so an interrupted unpack leaves nothing that
    paths.py would accept as a bundle."""
    tmp = os.path.join(dest, "." + name + ".unpacking")
    shutil.rmtree(tmp, ignore_errors=True)
    os.makedirs(tmp)
    try:
        with _quiet(), tarfile.open(archive, "r:gz") as tar:
            members = []
            for m in tar.getmembers():
                p = os.path.normpath(m.name)
                if os.path.isabs(p) or p.startswith("..") or not (p == name or p.startswith(name + os.sep)):
                    raise ValueError(f"{os.path.basename(archive)} is not a checkGM bundle archive: it holds '{m.name}', which is "
                                     f"outside the directory '{name}' a bundle unpacks to")
                if not (m.isfile() or m.isdir()):
                    raise ValueError(f"{os.path.basename(archive)} is not a checkGM bundle archive: '{m.name}' is neither a file "
                                     f"nor a directory")
                members.append(m)
            tar.extractall(tmp, members=members)
        out = os.path.join(dest, name)
        shutil.rmtree(out, ignore_errors=True)
        os.replace(os.path.join(tmp, name), out)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    return out


# --------------------------------------------------------------------------------------------- the three entry points

def plan(pipeline=None) -> list:
    """What a setup would download and roughly how large, without touching the network.

    One row per thing fetched: the default bundles, plus the named pipeline's databases when one is given. The CLI shows
    this and asks before anything starts, which for the assembly databases is the difference between a considered
    download and a surprised one (90 GB over the wire, 170 GB on disk)."""
    rows = []
    for n in select(None):
        b = BUNDLES[n]
        rows.append({"item": n, "kind": "bundle", "pipeline": b["profile_type"], "n_files": 1,
                     "download_bytes": b["bytes"], "disk_bytes": b["bytes"] + b["unpacked_bytes"],
                     "source": f"zenodo doi {DOI}",
                     "what": f"{'assembly-based' if b['profile_type'] == 'assembly' else 'read-based'} adult reference baseline "
                             f"({b['tier']} tier)",
                     "note": "the deposit is not public yet, so a location has to be given with --url until it is"})
    if pipeline is not None:
        p = _pipeline(pipeline)
        d = DATABASES[p]
        rows.append({"item": d["script"], "kind": "databases", "pipeline": p, "n_files": d["n_files"],
                     "download_bytes": d["download_bytes"], "disk_bytes": d["disk_bytes"],
                     "source": f"container/{d['script']} (pinned URL list)", "what": d["what"], "note": d["note"]})
    return rows


def bundles(dest, url=None, doi=None, only=None, dry_run=False) -> dict:
    """Download the baseline bundles into `dest`, verify them against their recorded sha256, and unpack them there.

    `dest` is then usable as ``$CHECKGM_BUNDLES``: each archive unpacks to a directory named by the bundle id, which is
    what paths.resolve_bundle looks for. A bundle whose archive is already present and verifies is not fetched again and
    a bundle directory already in place is not rebuilt, so an interrupted setup is resumed by re-running the same
    command; a partial transfer is resumed from its part-file.

    `url` is a directory URL, a local directory or (for a single bundle) one archive; `doi` overrides the Zenodo record;
    with neither, the reserved DOI is resolved, which today fails with the message _not_public_message writes. `only`
    picks bundles (see `select`). With `dry_run`, nothing is written or transferred and the per-bundle rows say what
    would be.

    Returns {dest, source, dry_run, total_bytes, downloaded_bytes, bundles: [{name, bytes, sha256, verified, action,
    archive, bundle_dir}]}: `bytes` and `sha256` are the recorded ones, and `verified` says the archive was hashed in
    this call and matched."""
    names = select(only)
    dest = os.path.abspath(dest)
    out = {"dest": dest, "source": None, "dry_run": bool(dry_run), "total_bytes": sum(BUNDLES[n]["bytes"] for n in names),
           "downloaded_bytes": 0, "bundles": []}

    if dry_run:
        # No network either: a dry run must be answerable with the record unpublished, which is the whole point of it.
        out["source"] = (os.path.abspath(url) if url and not _is_url(url) else url) or f"zenodo doi {doi or DOI}"
        for n in names:
            b = BUNDLES[n]
            archive = os.path.join(dest, n + ARCHIVE_SUFFIX)
            here = os.path.isfile(archive) and os.path.getsize(archive) == b["bytes"]
            out["bundles"].append({"name": n, "bytes": b["bytes"], "sha256": b["sha256"], "verified": False,
                                   "action": "present, would verify" if here else "would download",
                                   "archive": archive, "bundle_dir": os.path.join(dest, n)})
        return out

    sources, label = _sources(names, url, doi)
    out["source"] = label
    os.makedirs(dest, exist_ok=True)
    for n in names:
        b = BUNDLES[n]
        src, src_size = sources[n]
        archive = os.path.join(dest, n + ARCHIVE_SUFFIX)
        bundle_dir = os.path.join(dest, n)
        row = {"name": n, "bytes": b["bytes"], "sha256": b["sha256"], "verified": False, "action": None,
               "archive": archive, "bundle_dir": bundle_dir}

        if src_size is not None and src_size != b["bytes"]:
            raise ValueError(f"{n}{ARCHIVE_SUFFIX} at the source is {human_bytes(src_size)} ({src_size} bytes) but this release "
                             f"records {human_bytes(b['bytes'])} ({b['bytes']} bytes); the location holds a different build of "
                             f"the bundle, so upgrade checkgm or point --url at the archive matching this one")

        unpacked = os.path.isfile(os.path.join(bundle_dir, "manifest.json"))
        if os.path.isfile(archive) and os.path.getsize(archive) == b["bytes"] and _sha256(archive) == b["sha256"]:
            row.update(verified=True, action="skipped")
            _say(f"{n} archive present and verified ({human_bytes(b['bytes'])})")
        elif os.path.isfile(archive) and unpacked:
            # Present, unpacked, and the wrong bytes: re-fetching over a bundle the user may be mid-analysis with is worse
            # than stopping, and the two remedies (a stale archive to delete, a corrupt one to replace) are different.
            raise ValueError(f"{archive} does not match the recorded sha256 for {n} and {bundle_dir} is already unpacked; delete "
                             f"the archive to fetch it again, or delete both to start over")
        elif unpacked:
            # The bundle is in place and its archive was deleted, which is what a user short of space does. There is
            # nothing to verify and nothing to fetch; saying so is better than spending the download to re-prove it.
            row.update(verified=False, action="skipped (unpacked; archive not kept, so nothing to verify)")
            _say(f"{n} already unpacked in {bundle_dir}; its archive is not here, so its checksum was not re-checked")
        else:
            _say(f"{n} downloading {human_bytes(b['bytes'])} from {src}")
            _download(src, archive, b["bytes"])
            got = _sha256(archive)
            if got != b["sha256"]:
                size = os.path.getsize(archive)
                os.remove(archive)   # keeping it would make the next run skip the download and report the same failure
                raise ValueError(f"{n}{ARCHIVE_SUFFIX} failed verification: sha256 {got} over {size} bytes, expected "
                                 f"{b['sha256']} over {b['bytes']} bytes; the transfer was incomplete or the source holds a "
                                 f"different file, so run the command again")
            out["downloaded_bytes"] += os.path.getsize(archive)
            row.update(verified=True, action="downloaded")
            _say(f"{n} verified sha256 {got}")

        if not unpacked:
            row["bundle_dir"] = _unpack(archive, dest, n)
            _say(f"{n} unpacked to {row['bundle_dir']}")
        out["bundles"].append(row)

    _say(f"bundles ready in {dest}; export CHECKGM_BUNDLES={dest}")
    return out


def _pipeline(pipeline) -> str:
    """The canonical pipeline name ('assembly' or 'read'), accepting the bundles' word 'reads'."""
    p = PIPELINE_ALIASES.get(str(pipeline), str(pipeline))
    if p not in DATABASES:
        raise ValueError(f"unknown pipeline '{pipeline}'; it is 'assembly' (assembly-based, MGnify v5) or 'read' "
                         f"(read-based, MetaPhlAn 3)")
    return p


def container_dir() -> str:
    """The container/ directory holding the pinned fetch scripts.

    The wheel ships the package alone, so the scripts are wherever the repository was cloned or unpacked: $CHECKGM_CONTAINER
    if set, else the checkout this module is running from, else ./container."""
    cands = []
    if os.environ.get(CONTAINER_ENV):
        cands.append(os.environ[CONTAINER_ENV])
    here = os.path.dirname(os.path.abspath(__file__))
    cands.append(os.path.join(os.path.dirname(os.path.dirname(here)), "container"))
    cands.append(os.path.join(os.getcwd(), "container"))
    for c in cands:
        if os.path.isfile(os.path.join(c, DATABASES["assembly"]["script"])):
            return os.path.abspath(c)
    raise FileNotFoundError(f"the pipeline fetch scripts were not found (looked in {', '.join(cands)}): they live in container/ of "
                            f"the repository ({REPO}), which the installed package does not ship, so clone it and set "
                            f"{CONTAINER_ENV}=/path/to/checkgm/container")


def _missing_tools(pipeline: str) -> list[str]:
    """Which of the pipeline's external tools are not on PATH."""
    return [t for t in DATABASES[pipeline]["tools"] if shutil.which(t) is None]


def _run_script(script: str, argv: list) -> tuple[int, list]:
    """Run a fetch script, echoing its output to stdout as it goes, and return (exit status, lines it marked failed).

    The output is relayed rather than captured and printed at the end because these scripts run for hours. fetch_dbs_assembly.sh
    is not `set -e` by design (one failed URL must not abort the other 23, the run being resumable), so its exit status
    alone does not say whether everything arrived: the [FAIL] lines it prints do."""
    failed = []
    proc = subprocess.Popen(argv, cwd=os.path.dirname(script), stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    for line in proc.stdout:
        line = line.rstrip("\n")
        if "[FAIL]" in line or "[md5 FAIL]" in line:
            failed.append(line)
        _say(line)
    return proc.wait(), failed


def databases(pipeline, dest, dry_run=False) -> dict:
    """Fetch the reference databases a pipeline needs, by running the pinned scripts in container/.

    The scripts are the authority on what a pipeline's measurement needs — fetch_dbs_assembly.sh then unpack_dbs.sh for the
    assembly-based pipeline, fetch_dbs_read.sh (which also builds the Bowtie2 index) for the read-based one — and both are
    resumable, so this is re-run after an interruption. They are invoked through bash with `dest` as their argument, not
    executed directly, since not all of them carry the executable bit.

    With `dry_run`, nothing runs: the result carries the size of what would be fetched and the exact argv, and reports
    any missing tool rather than refusing, so the plan can be seen on a host that would not be the one to run it.

    Returns {pipeline, dest, dry_run, download_bytes, disk_bytes, steps: [{script, argv, returncode}], ok}."""
    p = _pipeline(pipeline)
    d = DATABASES[p]
    dest = os.path.abspath(dest)
    cdir = container_dir()
    scripts = [os.path.join(cdir, s) for s in (d["script"], d["then"]) if s]
    for s in scripts:
        if not os.path.isfile(s):
            raise FileNotFoundError(f"{os.path.basename(s)} not found in {cdir}; the {p}-based pipeline needs it, so check out the "
                                    f"repository's container/ directory ({REPO}) or set {CONTAINER_ENV} to it")
    out = {"pipeline": p, "dest": dest, "dry_run": bool(dry_run), "download_bytes": d["download_bytes"],
           "disk_bytes": d["disk_bytes"], "steps": [{"script": s, "argv": ["bash", s, dest], "returncode": None} for s in scripts],
           "ok": None}

    missing = _missing_tools(p)
    if dry_run:
        out["missing_tools"] = missing
        return out
    if missing:
        where = ("container/checkgm_assembly.def" if p == "assembly" else "container/checkgm_read.def")
        raise ValueError(f"the {p}-based database setup needs {', '.join(missing)} on PATH and {'it is' if len(missing) == 1 else 'they are'} "
                         f"not there; run it inside the pipeline image ({where}) or install {'it' if len(missing) == 1 else 'them'} first")
    if shutil.which("bash") is None:
        raise ValueError("bash is not on PATH, and the pipeline database scripts are bash scripts")

    os.makedirs(dest, exist_ok=True)
    _say(f"{p}-based databases: fetching {human_bytes(d['download_bytes'])} into {dest} "
         f"({human_bytes(d['disk_bytes'])} on disk when unpacked)")
    for step in out["steps"]:
        _say(f"running {os.path.basename(step['script'])} {dest}")
        code, failed = _run_script(step["script"], step["argv"])
        step["returncode"] = code
        if code != 0:
            out["ok"] = False
            raise ValueError(f"{os.path.basename(step['script'])} exited {code}; the databases in {dest} are incomplete, and the "
                             f"script is resumable, so fix what it reported and run the same command again")
        if failed:
            out["ok"] = False
            raise ValueError(f"{os.path.basename(step['script'])} could not fetch {len(failed)} of {d['n_files']} files "
                             f"({', '.join(l.split()[-1] for l in failed if len(l.split()) > 1)}); it is resumable, so run the same "
                             f"command again to retry only those")
    out["ok"] = True
    _say(f"{p}-based databases ready in {dest}")
    return out
