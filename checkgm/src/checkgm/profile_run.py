"""checkgm.profile_run — run a measurement pipeline on raw reads (`checkgm profile`).

The measurement itself lives in the container scripts, which pin every tool version the reference was measured with
(`container/run_tierA.sh` for the assembly-based pipeline, `container/run_pipelineB.sh` for the read-based one). Nothing
here reimplements a step: this module decides whether the run can happen at all, rewrites host paths into the bind mounts
the image sees, builds one `apptainer run` invocation, and reports where the output `checkgm assess` consumes ended up.

Two shapes of output, one per pipeline, because the two baselines are measured differently:
  assembly  OUT/query/                     a query directory (docs/query_format.md) -> `checkgm assess` against gut-assembly-*
  read      OUT/<sample>.txt + .reads.tsv  a MetaPhlAn 3 profile and its read count -> `checkgm assess` against gut-reads-*
HUMAnN 3 tables are an input to assess rather than a product of either pipeline (neither run script calls HUMAnN;
run_pipelineB.sh only forwards tables it is handed), so this module never claims to have made them.

Both scripts are resumable in place (tierA stamps `.done_<step>`, pipelineB skips MetaPhlAn when the profile and read
count are already there), so re-running `run()` on the same --out continues rather than starts over.

Three environment variables name what a pip install cannot provide: $CHECKGM_IMAGES (colon-separated directories holding
the .sif files, searched like $CHECKGM_BUNDLES), $CHECKGM_DBS (the database directory, or one holding an `assembly` and a
`read` subdirectory) and $CHECKGM_APPTAINER (the runtime, when it is not on PATH). All three are only defaults for the
explicit arguments, and every one of them missing is a sentence from preflight() rather than a failure inside a tool.
"""
from __future__ import annotations

import os
import re
import shlex
import shutil
import subprocess
import sys
import time

PIPELINES = ("assembly", "read")
WORD = {"assembly": "assembly-based", "read": "read-based"}

# Where this module looks when --image / --dbs are not given. The names are this module's own; whatever
# `checkgm setup` ends up writing has to agree with them, so each is one constant.
IMAGES_ENV_VAR = "CHECKGM_IMAGES"   # colon-separated directories, searched like CHECKGM_BUNDLES
DBS_ENV_VAR = "CHECKGM_DBS"         # one directory per pipeline kind; see DBS_SUBDIR
RUNTIME_ENV_VAR = "CHECKGM_APPTAINER"
RUNTIMES = ("apptainer", "singularity")   # singularity is the pre-rename name of the same runtime and takes the same options

# image and recipe names as container/*.def and the README's distribution table spell them
IMAGE_NAME = {"assembly": "checkgm_tierA.sif", "read": "checkgm_pipelineB.sif"}
DEF_NAME = {"assembly": "checkgm_tierA.def", "read": "checkgm_pipelineB.def"}
# a database directory may hold both pipelines' data side by side, so each kind also gets a conventional subdirectory
DBS_SUBDIR = {"assembly": "assembly", "read": "read"}

# Files the run scripts open by name. Checking them is the difference between "checkgm setup db has not finished" and a
# crash from the middle of a shell script: an incomplete download is the normal failure, not an absent directory.
REQUIRED_DBS = {
    "assembly": ["uniref90_v2019_11_diamond-v0.9.25.dmnd", "db_uniref90_v2019_11.txt.gz", "eggnog.db", "eggnog_proteins.dmnd",
                 "db_kofam.hmm", "db_kofam.hmm.h3f", "graphs-20200805.pkl", "all_pathways_class.txt", "all_pathways_names.txt",
                 "interproscan-5.36-75.0/interproscan.sh"],
    "read": ["mpa_v30_CHOCOPhlAn_201901.pkl", "mpa_v30_CHOCOPhlAn_201901.1.bt2"],
}
# The scripts behind 'checkgm setup db', named for anyone without the new command; their sizes are not repeated here
# because checkgm.setup_data records the measured ones and says the README's round figures are too small.
FETCH_SCRIPT = {"assembly": "container/fetch_dbs.sh then container/unpack_dbs.sh", "read": "container/fetch_dbs_B.sh"}
ASSEMBLY_CPUS, ASSEMBLY_MEM_GB = 16, 128   # README, "From raw reads": roughly this much per sample
TIERA_DEFAULT_MEM_GB = 120                 # run_tierA.sh MEMGB default, i.e. what the validated runs used

# $SAMPLE becomes a file name (run_pipelineB.sh writes $OUT/$SAMPLE.txt) and an unquoted word in places inside the
# scripts, so the ids that cannot survive the trip are refused here rather than producing a file nobody can find.
SAMPLE_RE = re.compile(r"[A-Za-z0-9][A-Za-z0-9._-]*")
READS_PROCESSED_RE = re.compile(r"#\s*([\d,]+)\s+reads processed")   # same header line normalize_reads.py reads
STDERR_LOG = "profile.stderr.log"
PIPELINE_LOG = {"assembly": "tierA.log", "read": "pipelineB.log"}


def _warn(msg: str):
    """One diagnostic on stderr in the package's single form, `checkgm: ...`.

    Not imported from checkgm.cli: cli imports this module, and a notice is not worth an import cycle."""
    print(f"checkgm: {msg}", file=sys.stderr, flush=True)


def _one_line(s: str, limit: int = 240) -> str:
    """Collapse captured tool output to one clause, so that a failure still reports as a single stderr line."""
    s = " ".join(s.split())
    return s if len(s) <= limit else s[:limit] + " ..."


# --------------------------------------------------------------------------------------------- what the run needs

def runtime() -> str | None:
    """The container runtime to invoke, or None when there is none.

    $CHECKGM_APPTAINER wins (a module-built apptainer is often not on a login shell's PATH) and is honoured only if it
    is really executable, so a stale setting is reported as a missing runtime rather than failing at exec time."""
    given = os.environ.get(RUNTIME_ENV_VAR)
    if given:
        return given if os.path.isfile(given) and os.access(given, os.X_OK) else None
    return next((p for p in (shutil.which(n) for n in RUNTIMES) if p), None)


def _search_dirs(var: str) -> list[str]:
    return [d for d in os.environ.get(var, "").split(os.pathsep) if d]


def resolve_image(pipeline: str, image: str | None) -> tuple[str | None, str | None]:
    """(image path, reason it is unusable) — exactly one of the two is None."""
    if image:
        if os.path.isfile(image):
            return os.path.abspath(image), None
        kind = "is a directory, not an image file" if os.path.isdir(image) else "not found"
        return None, f"pipeline image {kind}: {image}"
    for d in _search_dirs(IMAGES_ENV_VAR):
        cand = os.path.join(d, IMAGE_NAME[pipeline])
        if os.path.isfile(cand):
            return os.path.abspath(cand), None
    where = ", ".join(_search_dirs(IMAGES_ENV_VAR))
    # The remedy named here is the build, because it is the only one that always works: the images are distributed beside
    # the repository rather than from a public URL (README's download table, Zenodo DOI still TBD), and checkgm.setup_data
    # fetches bundles and databases, not images.
    return None, (f"no {WORD[pipeline]} pipeline image: " +
                  (f"{IMAGE_NAME[pipeline]} is not under ${IMAGES_ENV_VAR} = {where}" if where else
                   f"neither --image nor ${IMAGES_ENV_VAR} was given") +
                  f"; pass --image, or build it with 'cd container && apptainer build {IMAGE_NAME[pipeline]} {DEF_NAME[pipeline]}'")


def resolve_dbs(pipeline: str, dbs: str | None) -> tuple[str | None, str | None]:
    """(database directory, reason it is unusable) — exactly one of the two is None.

    A directory that exists but lacks files the run script opens is the common state halfway through a download, and it
    is named as such: the shell script would otherwise fail several minutes in, inside a tool, on one missing path."""
    cands = [dbs] if dbs else [os.path.join(d, DBS_SUBDIR[pipeline]) for d in _search_dirs(DBS_ENV_VAR)] + _search_dirs(DBS_ENV_VAR)
    found = next((c for c in cands if c and os.path.isdir(c)), None)
    if not found:
        where = (f"not found: {dbs}" if dbs else
                 f"not given and ${DBS_ENV_VAR} " + ("names no directory that exists" if _search_dirs(DBS_ENV_VAR) else "is not set"))
        return None, (f"database directory for the {WORD[pipeline]} pipeline {where}; get it with 'checkgm setup db' "
                      f"(or with {FETCH_SCRIPT[pipeline]}, which is what it runs)")
    missing = [f for f in REQUIRED_DBS[pipeline] if not os.path.exists(os.path.join(found, f))]
    if missing:
        shown = ", ".join(missing[:3]) + (f" and {len(missing) - 3} more" if len(missing) > 3 else "")
        return None, (f"database directory {found} is incomplete for the {WORD[pipeline]} pipeline: {len(missing)} of "
                      f"{len(REQUIRED_DBS[pipeline])} expected files are missing ({shown}); finish 'checkgm setup db', "
                      f"which is resumable and fetches only what is absent")
    return os.path.abspath(found), None


def preflight(pipeline: str, dbs: str | None, image: str | None) -> list:
    """The reasons this run cannot proceed, as human-readable strings; empty when it can.

    Every reason is one clause with its own remedy, so that a caller can print one of them or join them into the single
    `checkgm: ...` line the CLI promises. Checked here: the pipeline name, the container image, the database directory
    (its contents, not merely its existence) and the container runtime. Per-sample inputs are not visible at this
    signature and are checked by run()."""
    if pipeline not in PIPELINES:
        return [f"unknown pipeline '{pipeline}'; it is one of {', '.join(PIPELINES)}"]
    reasons = []
    for _, why in (resolve_image(pipeline, image), resolve_dbs(pipeline, dbs)):
        if why:
            reasons.append(why)
    if not runtime():
        reasons.append(f"no container runtime: the pipelines only run inside the image, and {' / '.join(RUNTIMES)} is not on "
                       f"PATH" + (f" (${RUNTIME_ENV_VAR} = {os.environ[RUNTIME_ENV_VAR]} is not executable)"
                                  if os.environ.get(RUNTIME_ENV_VAR) else f"; install Apptainer or set ${RUNTIME_ENV_VAR} to it"))
    return reasons


# --------------------------------------------------------------------------------------------- inputs and bind mounts

def _positive_int(v, flag: str) -> int:
    try:
        n = int(v)
    except (TypeError, ValueError):
        raise ValueError(f"{flag} must be a positive integer, not {v!r}") from None
    if n < 1:
        raise ValueError(f"{flag} must be a positive integer, not {v!r}")
    return n


def _check_inputs(pipeline: str, r1, r2, contigs, profile_in) -> dict:
    """{flag: path} of the inputs to pass on, after refusing the combinations neither script can act on.

    The wrong-pipeline cases (contigs for the read pipeline, a profile for the assembly one) are told apart from a
    missing input because the correction differs: one is the other subcommand, the other is a forgotten file."""
    wrong = {"assembly": [("profile_in", profile_in, "a MetaPhlAn profile is read-based input, so use --pipeline read")],
             "read": [("contigs", contigs, "contigs are assembly-based input, so use --pipeline assembly")]}
    for name, val, fix in wrong[pipeline]:
        if val:
            raise ValueError(f"--{name.replace('_', '-')} cannot be used with the {WORD[pipeline]} pipeline: {fix}")
    if r2 and not r1:
        raise ValueError("--r2 was given without --r1; the first mate names the library")
    if profile_in:
        given = {"profile_in": profile_in}
    elif contigs:
        if r1:
            raise ValueError("--contigs and --r1 are alternatives: with contigs the assembly step is skipped, with reads it is run")
        given = {"contigs": contigs}
    elif r1:
        given = {"r1": r1} | ({"r2": r2} if r2 else {})
    else:
        want = "--r1 (with --r2 for paired reads), or --contigs to start from an assembly" if pipeline == "assembly" else \
               "--r1 (with --r2 for paired reads), or --profile-in to use a MetaPhlAn profile you already have"
        raise ValueError(f"the {WORD[pipeline]} pipeline needs reads: give {want}")
    for name, p in given.items():
        if not os.path.isfile(p):
            raise FileNotFoundError(f"--{name.replace('_', '-')} file not found: {p}")
    return given


def bind_plan(out: str, dbs: str, inputs: list) -> tuple[list, dict]:
    """(bind specifications, {host path: path inside the image}) for one invocation.

    The run scripts take absolute paths and see only what is bound, so every input has to be reachable inside the
    image. Output is bound writable at /out and the databases read-only at /dbs; an input already under one of those
    is addressed through it rather than mounted twice (apptainer would otherwise shadow the writable mount with a
    read-only one), and every other input directory gets its own read-only /inN."""
    mounts = [(os.path.realpath(out), "/out", ""), (os.path.realpath(dbs), "/dbs", "ro")]
    paths = {}
    for p in inputs:
        pr = os.path.realpath(p)
        for root, cpath, _ in mounts:
            if pr == root or pr.startswith(root + os.sep):
                paths[p] = cpath if pr == root else cpath + "/" + os.path.relpath(pr, root).replace(os.sep, "/")
                break
        else:
            parent = os.path.dirname(pr)
            cpath = next((c for r, c, _ in mounts if r == parent), None)
            if cpath is None:
                cpath = f"/in{sum(1 for _, c, _ in mounts if c.startswith('/in'))}"
                mounts.append((parent, cpath, "ro"))
            paths[p] = cpath + "/" + os.path.basename(pr)
    return [f"{r}:{c}" + (f":{m}" if m else "") for r, c, m in mounts], paths


def build_command(pipeline: str, sample: str, out: str, dbs: str, image: str, inputs: dict,
                  threads: int, mem_gb: int | None, exe: str = "apptainer") -> list:
    """The single `apptainer run` argv that performs this run, as a list ready for subprocess.

    --cleanenv keeps a host conda, PYTHONPATH or TMPDIR out of the image, whose %environment is what points the scripts
    at the pinned interpreters and tools; it does not touch that %environment. The read-based script was written to
    carry one sample all the way to a score and so insists on a --bundle, but the placement is `checkgm assess`'s job
    now, so its checkgm step is disabled through the hook the script itself provides (CHECKGM_BIN) and the bundle
    argument it never reaches is passed as a visible placeholder."""
    binds, cpath = bind_plan(out, dbs, list(inputs.values()))
    argv = [exe, "run", "--cleanenv", "--bind", ",".join(binds)]
    if pipeline == "read":
        argv += ["--env", "CHECKGM_BIN=/bin/true"]
    argv += [image, "--sample", sample, "--out", "/out", "--dbs", "/dbs"]
    for flag in ("r1", "r2", "contigs"):
        if flag in inputs:
            argv += [f"--{flag}", cpath[inputs[flag]]]
    argv += ["--threads", str(threads)]
    if pipeline == "assembly":
        if mem_gb is not None:
            argv += ["--mem-gb", str(mem_gb)]
    else:
        argv += ["--bundle", "not-used-by-checkgm-profile"]
    return argv


# --------------------------------------------------------------------------------------------- running

def _reads_processed(profile: str) -> int | None:
    """The read count in a MetaPhlAn profile's header, or None when the header carries none.

    Read with the same rule as normalize_reads.py, so that `profile` reports the depth band `assess` will actually see."""
    with open(profile, errors="replace") as fh:
        for ln in fh.readlines()[:50]:
            m = READS_PROCESSED_RE.match(ln)
            if m:
                return int(m.group(1).replace(",", ""))
    return None


def _looks_like_profile(path: str) -> bool:
    """True when the file's first lines mention clade_name, as both MetaPhlAn shapes do (single-sample and merged)."""
    with open(path, errors="replace") as fh:
        return any("clade_name" in ln for ln in fh.readlines()[:50])


def _stage_profile(sample: str, out: str, profile_in: str, dry_run: bool) -> dict:
    """Place an already-measured MetaPhlAn profile under --out, which is all `profile` has left to do for it.

    No image, databases or runtime are involved, so none are demanded: a user who profiled elsewhere can still reach
    `assess` through the same command. The file is copied rather than linked, so that the output directory keeps
    standing on its own once the original moves."""
    if not _looks_like_profile(profile_in):
        raise ValueError(f"--profile-in is not a MetaPhlAn profile (no clade_name column in its first lines): {profile_in}; "
                         f"pass the profiler's table, or --r1/--r2 to measure the reads")
    dest = os.path.join(out, f"{sample}.txt")
    same = os.path.exists(dest) and os.path.samefile(profile_in, dest)
    if dry_run:
        print("already in place" if same else shlex.join(["cp", os.path.abspath(profile_in), dest]), flush=True)
    else:
        os.makedirs(out, exist_ok=True)
        if not same:
            shutil.copyfile(profile_in, dest)
    notes = []
    if not dry_run and _reads_processed(dest) is None:
        notes.append(f"{os.path.basename(dest)} has no '#N reads processed' header line, so its depth band is unknown unless "
                     f"'checkgm assess' is given --reads-tsv (read-based baselines do not use bands for percentiles, so it is "
                     f"still placed)")
    for n in notes:
        _warn(n)
    return {"pipeline": "read", "sample": sample, "out": out, "status": "planned" if dry_run else "staged",
            "command": [] if same else ["cp", os.path.abspath(profile_in), dest], "outputs": {"profile": dest},
            "assess_input": dest, "log": None, "stderr_log": None, "threads": None, "mem_gb": None,
            "elapsed_s": 0.0, "blockers": [], "notes": notes}


def _resource_notes(pipeline: str, threads: int, mem_gb: int | None) -> list:
    """Warnings about a request the pipeline is known not to fit in, which no preflight can turn into a refusal."""
    notes = []
    if pipeline == "assembly":
        if threads < ASSEMBLY_CPUS:
            notes.append(f"assembly needs roughly {ASSEMBLY_CPUS} CPUs and {ASSEMBLY_MEM_GB} GB of memory per sample; this run "
                         f"asks for {threads}")
        if mem_gb is not None and mem_gb < TIERA_DEFAULT_MEM_GB:
            notes.append(f"--mem-gb {mem_gb} is below the {TIERA_DEFAULT_MEM_GB} GB the validated assembly runs used; metaSPAdes "
                         f"may run out of memory on a deep library")
    elif mem_gb is not None:
        notes.append("--mem-gb is ignored by the read-based pipeline, whose run script takes no memory argument")
    return notes


def _outputs(pipeline: str, sample: str, out: str) -> dict:
    """The files `assess` consumes, by the names the run scripts give them."""
    if pipeline == "assembly":
        return {"query": os.path.join(out, "query")}
    return {"profile": os.path.join(out, f"{sample}.txt"), "reads_tsv": os.path.join(out, f"{sample}.reads.tsv")}


def run(pipeline, sample, out, r1=None, r2=None, contigs=None, profile_in=None, dbs=None,
        threads=8, mem_gb=None, image=None, dry_run=False) -> dict:
    """Run one sample through a measurement pipeline and return what it produced.

    pipeline is "assembly" (container/run_tierA.sh: assembly + MGnify v5-equivalent annotation -> OUT/query) or "read"
    (container/run_pipelineB.sh: MetaPhlAn 3 -> OUT/<sample>.txt and OUT/<sample>.reads.tsv). --contigs starts the
    assembly pipeline from an assembly; --profile-in skips the read pipeline altogether and only places the profile.

    Anything the user can correct raises ValueError or FileNotFoundError with a one-line message, which is the form the
    CLI turns into `checkgm: ...` and exit 1 — including the pipeline failing, since a stack of tool output is not a
    user error message. The tools' stderr is kept in OUT/profile.stderr.log and their progress is left on stdout.

    dry_run prints the command (or the copy) that would run, touches nothing, and reports what is missing instead of
    refusing: a dry run on a laptop without the image is exactly where the command line is worth seeing.

    Returns {pipeline, sample, out, status ("ok", "staged" or "planned"), command, outputs, assess_input, log,
    stderr_log, threads, mem_gb, elapsed_s, blockers, notes}; assess_input is the one path to hand to `checkgm assess`."""
    if pipeline not in PIPELINES:
        raise ValueError(f"unknown pipeline '{pipeline}'; it is one of {', '.join(PIPELINES)}")
    sample = str(sample or "").strip()
    if not SAMPLE_RE.fullmatch(sample):
        raise ValueError(f"sample id {sample!r} cannot be used: it becomes a file name inside the container, so it must be "
                         f"letters, digits, dot, dash or underscore and start with a letter or digit")
    if not out:
        raise ValueError("no output directory given")
    out = os.path.abspath(out)
    inputs = _check_inputs(pipeline, r1, r2, contigs, profile_in)
    threads = _positive_int(threads, "--threads")
    mem_gb = None if mem_gb is None else _positive_int(mem_gb, "--mem-gb")
    if "profile_in" in inputs:
        return _stage_profile(sample, out, inputs["profile_in"], dry_run)

    blockers = preflight(pipeline, dbs, image)
    if blockers and not dry_run:
        raise ValueError("; ".join(blockers))
    img = resolve_image(pipeline, image)[0]
    dbdir = resolve_dbs(pipeline, dbs)[0]
    exe = runtime()
    if dry_run:
        # Placeholders stand in for whatever could not be resolved, so the command is still printable and still shows
        # which argument the missing piece would have filled. They are absolute and obviously not real paths, since
        # bind_plan resolves what it is given and would otherwise turn a bracketed word into a path under the cwd.
        img = img or (os.path.abspath(image) if image else f"/path/to/{IMAGE_NAME[pipeline]}")
        dbdir = dbdir or (os.path.abspath(dbs) if dbs else "/path/to/databases")
        exe = exe or RUNTIMES[0]
    argv = build_command(pipeline, sample, out, dbdir, img, inputs, threads, mem_gb, exe)
    notes = _resource_notes(pipeline, threads, mem_gb)
    if dry_run:
        print(shlex.join(argv), flush=True)
        for n in blockers + notes:
            _warn(n)
        return {"pipeline": pipeline, "sample": sample, "out": out, "status": "planned", "command": argv,
                "outputs": _outputs(pipeline, sample, out), "assess_input": None,
                "log": os.path.join(out, PIPELINE_LOG[pipeline]), "stderr_log": os.path.join(out, STDERR_LOG),
                "threads": threads, "mem_gb": mem_gb, "elapsed_s": 0.0, "blockers": blockers, "notes": notes}

    for n in notes:
        _warn(n)
    os.makedirs(out, exist_ok=True)
    err_log = os.path.join(out, STDERR_LOG)
    t0 = time.time()
    with open(err_log, "ab") as fh:   # to a file, not a pipe: a pipeline's stderr is long, and none of it belongs on ours
        fh.write(f"\n### {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime(t0))} {shlex.join(argv)}\n".encode())
        fh.flush()
        rc = subprocess.run(argv, stdout=None, stderr=fh).returncode
    elapsed = time.time() - t0
    log = os.path.join(out, PIPELINE_LOG[pipeline])
    if rc != 0:
        raise ValueError(f"the {WORD[pipeline]} pipeline failed for sample {sample} (exit {rc}, {elapsed:.0f}s); "
                         f"its own log is {log} and the tools' stderr is {err_log}{_last_message(err_log, log)}")
    outputs = _outputs(pipeline, sample, out)
    missing = [p for k, p in outputs.items() if not (os.path.isdir(p) if k == "query" else os.path.isfile(p))]
    if missing:
        raise ValueError(f"the {WORD[pipeline]} pipeline reported success for sample {sample} but did not write "
                         f"{', '.join(os.path.basename(p) for p in missing)} under {out}; see {log}")
    return {"pipeline": pipeline, "sample": sample, "out": out, "status": "ok", "command": argv, "outputs": outputs,
            "assess_input": outputs["query"] if pipeline == "assembly" else outputs["profile"], "log": log,
            "stderr_log": err_log, "threads": threads, "mem_gb": mem_gb, "elapsed_s": round(elapsed, 1),
            "blockers": [], "notes": notes}


def _last_message(err_log: str, log: str) -> str:
    """One clause quoting the last thing the pipeline said, from the tools' stderr or else the script's own log.

    A reader given only an exit status has to open two files to learn whether the run died of a missing database or of a
    full disk; one line of the real message usually settles it, and the paths are in the sentence either way."""
    for p in (err_log, log):
        try:
            lines = [l for l in open(p, errors="replace").read().splitlines() if l.strip() and not l.startswith("### ")]
        except OSError:
            continue
        if lines:
            return f"; last message: {_one_line(lines[-1])}"
    return ""
