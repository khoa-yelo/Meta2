"""Locating reference bundles.

A bundle argument on the command line is either a path to a bundle directory (one containing manifest.json) or a bundle id such
as ``gut-reads-adult-global-v0.2-lenient``. Ids are looked up in the directories listed in the environment variable
``REFMB_BUNDLES`` (colon-separated, like PATH). Nothing here depends on where the package was developed.
"""
from __future__ import annotations

import os

ENV_VAR = "REFMB_BUNDLES"
DEFAULT_ASSEMBLY_VERSION = "0.8"   # current pipeline A (assembly) bundle series
DEFAULT_READS_VERSION = "0.2"      # current pipeline B (read-based) bundle series


def bundle_dirs() -> list[str]:
    """Directories searched for bundle ids, from REFMB_BUNDLES."""
    v = os.environ.get(ENV_VAR, "")
    return [d for d in v.split(os.pathsep) if d]


def is_bundle(path: str) -> bool:
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


def bundle_A(tier: str = "lenient", loso: str | None = None, root: str | None = None, version: str = DEFAULT_ASSEMBLY_VERSION) -> str:
    """Name (or path, when root is given) of an assembly bundle: gut-assembly-adult-global-v<version>-<tier>[-loso-<study>]."""
    name = f"gut-assembly-adult-global-v{version}-{tier}" + (f"-loso-{loso}" if loso else "")
    return os.path.join(root, name) if root else name


def bundle_B(tier: str = "lenient", loso: str | None = None, root: str | None = None, version: str = DEFAULT_READS_VERSION) -> str:
    """Name (or path, when root is given) of a read-based bundle: gut-reads-adult-global-v<version>-<tier>[-loso-<study>]."""
    name = f"gut-reads-adult-global-v{version}-{tier}" + (f"-loso-{loso}" if loso else "")
    return os.path.join(root, name) if root else name
