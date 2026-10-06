"""Locating reference bundles.

A bundle argument on the command line is either a path to a bundle directory (one containing manifest.json) or a bundle id such
as ``gut-reads-adult-global-v0.2-lenient``. Ids are looked up in the directories listed in the environment variable
``CHECKGM_BUNDLES`` (colon-separated, like PATH). Nothing here depends on where the package was developed.
"""
from __future__ import annotations

import os
import re

ENV_VAR = "CHECKGM_BUNDLES"
LEGACY_ENV_VAR = "REFMB_BUNDLES"   # the name used before the tool was renamed; still read, so existing setups keep working
DEFAULT_ASSEMBLY_VERSION = "0.8"   # current assembly-based bundle series
DEFAULT_READS_VERSION = "0.2"      # current read-based bundle series


def bundle_dirs() -> list[str]:
    """Directories searched for bundle ids, from CHECKGM_BUNDLES (or the legacy REFMB_BUNDLES)."""
    v = os.environ.get(ENV_VAR) or os.environ.get(LEGACY_ENV_VAR, "")
    return [d for d in v.split(os.pathsep) if d]


def bundle_dirs_source() -> str:
    """The variable those directories were read from, so that a message names the one the user actually set.

    CHECKGM_BUNDLES unless only the pre-rename REFMB_BUNDLES carries a value; a user debugging the rename is told the
    wrong thing by a message that attributes their directories to a variable they never exported."""
    return LEGACY_ENV_VAR if not os.environ.get(ENV_VAR) and os.environ.get(LEGACY_ENV_VAR) else ENV_VAR


def legacy_env_note() -> str:
    """One parenthesis appended to a message when the directories came from the pre-rename variable, else ''."""
    return (f" ({LEGACY_ENV_VAR} is the pre-rename name and is still read; {ENV_VAR} is the preferred one)"
            if bundle_dirs_source() == LEGACY_ENV_VAR else "")


def is_bundle(path: str) -> bool:
    """True when `path` is a directory holding a manifest.json."""
    return os.path.isfile(os.path.join(path, "manifest.json"))


def _file_given_note(spec: str) -> str:
    """One clause appended when --bundle names an existing file instead of a bundle directory, else ''.

    Two near misses account for almost all of these, and each has its own correction: a path that reaches one level too
    far into an unpacked bundle, and a release archive that was never unpacked at all, the bundles being distributed as
    .tar.gz. Naming the remedy costs a clause and saves the user a guess."""
    if os.path.basename(spec) == "manifest.json":
        return f"; that is the manifest inside a bundle, so drop the file name and pass {os.path.dirname(spec) or '.'}"
    if re.search(r"\.(tar\.gz|tgz|tar|zip)$", spec):
        return "; that is a bundle archive rather than a bundle, so unpack it first and pass the directory it unpacks to"
    return "; a bundle is the directory holding manifest.json, not a file"


def resolve_bundle(spec: str) -> str:
    """Return the bundle directory for a path or a bundle id; raise FileNotFoundError with a plain explanation otherwise.

    Four mistakes are told apart, because each is corrected differently: a directory that exists but holds no
    manifest.json (one level off, the usual first-run slip), a file where a directory was wanted (an un-unpacked archive,
    or a path reaching down to manifest.json), a path that does not exist, and an id that was searched for under the
    bundle directories and not found. Only the last names those directories, since only it searched them."""
    if is_bundle(spec):
        return os.path.abspath(spec)
    bare = os.sep not in spec
    if bare:
        for d in bundle_dirs():
            cand = os.path.join(d, spec)
            if is_bundle(cand):
                return os.path.abspath(cand)
    if os.path.isdir(spec):
        try:
            inner = next((n for n in sorted(os.listdir(spec)) if is_bundle(os.path.join(spec, n))), None)
        except OSError:
            inner = None   # unreadable directory: the message is still the right one, only without the suggestion
        raise FileNotFoundError(f"bundle '{spec}' is a directory but holds no manifest.json: point --bundle at the directory the "
                                f"bundle archive unpacked to" + (f", which here looks like {os.path.join(spec, inner)}" if inner else
                                                                 " (the one directly containing manifest.json)"))
    if os.path.isfile(spec):
        raise FileNotFoundError(f"bundle directory not found: {spec}{_file_given_note(spec)}")
    if not bare:
        raise FileNotFoundError(f"bundle directory not found: {spec}")
    where = ", ".join(bundle_dirs()) or f"(set {ENV_VAR} to the directory holding your bundles)"
    raise FileNotFoundError(f"bundle '{spec}' not found: not a directory with manifest.json, and not under "
                            f"{bundle_dirs_source()} = {where}{legacy_env_note()}")


def list_bundles() -> list[str]:
    """Bundle directories found under CHECKGM_BUNDLES."""
    out = []
    for d in bundle_dirs():
        if os.path.isdir(d):
            out.extend(os.path.join(d, n) for n in sorted(os.listdir(d)) if is_bundle(os.path.join(d, n)))
    return out

