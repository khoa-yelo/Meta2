"""Locating reference bundles.

A bundle argument on the command line is either a path to a bundle directory (one containing manifest.json) or a bundle id such
as ``gut-reads-adult-global-v0.2-lenient``. Ids are looked up in the directories listed in the environment variable
``REFMB_BUNDLES`` (colon-separated, like PATH). Nothing here depends on where the package was developed.
"""
from __future__ import annotations

import os

ENV_VAR = "REFMB_BUNDLES"
DEFAULT_ASSEMBLY_VERSION = "0.8"   # current assembly-based bundle series
DEFAULT_READS_VERSION = "0.2"      # current read-based bundle series


def bundle_dirs() -> list[str]:
    """Directories searched for bundle ids, from REFMB_BUNDLES."""
    v = os.environ.get(ENV_VAR, "")
    return [d for d in v.split(os.pathsep) if d]


def is_bundle(path: str) -> bool:
    """True when `path` is a directory holding a manifest.json."""
    return os.path.isfile(os.path.join(path, "manifest.json"))


def resolve_bundle(spec: str) -> str:
    """Return the bundle directory for a path or a bundle id; raise FileNotFoundError with a plain explanation otherwise."""
    if is_bundle(spec):
        return os.path.abspath(spec)
    if os.sep not in spec:
        for d in bundle_dirs():
            cand = os.path.join(d, spec)
            if is_bundle(cand):
                return os.path.abspath(cand)
    where = ", ".join(bundle_dirs()) or f"(set {ENV_VAR} to the directory holding your bundles)"
    raise FileNotFoundError(f"bundle '{spec}' not found: not a directory with manifest.json, and not under {ENV_VAR} = {where}")


def list_bundles() -> list[str]:
    """Bundle directories found under REFMB_BUNDLES."""
    out = []
    for d in bundle_dirs():
        if os.path.isdir(d):
            out.extend(os.path.join(d, n) for n in sorted(os.listdir(d)) if is_bundle(os.path.join(d, n)))
    return out

