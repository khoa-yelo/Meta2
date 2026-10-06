#!/usr/bin/env python3
"""checkgm command line.

  checkgm normalize --bundle B --query DIR [DIR ...] --out OUT          # assembly-based: query directories -> OUT/normalized/<sample>.parquet + qc.tsv
  checkgm import-mgnify --bundle B --analysis-dir DIR [...] --out OUT   # assembly-based: MGnify v5 analysis folders -> the same
  checkgm score --bundle B --normalized OUT [--calibration CAL] --out REPORT   # -> scores.parquet, summaries.tsv, rejections.tsv, report.md
  checkgm run-normalize --bundle B --query DIR [...] --out REPORT               # normalize + score in one go
  checkgm run-mgnify --bundle B --analysis-dir DIR [...] --out REPORT           # import + score in one go
  checkgm import-metaphlan --bundle B --input PROFILE [...] [--reads-tsv sample<TAB>reads] --out OUT   # read-based: MetaPhlAn output -> normalized
  checkgm run-metaphlan --bundle B --input PROFILE [...] --out REPORT           # import + score against a read-based baseline
        [--humann-genefamilies FILE ...] [--humann-pathabundance FILE ...]   # optional HUMAnN 3 tables -> function layers
  checkgm bundles                                                               # list the bundles found under $CHECKGM_BUNDLES

B is a bundle directory or a bundle id looked up under $CHECKGM_BUNDLES (colon-separated directories).
Sample id = sample.json sample_id, else the directory name.

User errors (missing file, wrong kind of bundle, unreadable input) are reported on stderr as one line starting with
"checkgm:" and the command exits 1. Scientific refusals (a sample that fails a normalization floor) are not errors: the
sample appears in rejections.tsv and the command exits 0.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sys
import time

import pandas as pd

from checkgm import REPO, __version__
from checkgm.normalize import Backbone, BundleSpec, normalize, normalize_mgnify_v5, read_query_dir
from checkgm.paths import (DEFAULT_ASSEMBLY_VERSION, DEFAULT_READS_VERSION, ENV_VAR, bundle_dirs, bundle_dirs_source,
                           legacy_env_note, list_bundles, resolve_bundle)
from checkgm.score import Bundle, score

LAYERS = ["taxonomy_family", "taxonomy_genus", "taxonomy_species", "ko_eggnog", "ko_kofam", "pfam", "module", "ko_humann", "pathway_humann"]
QC_FLOAT_FORMAT = "%.17g"   # 17 significant digits round-trip every double exactly, so import + score == run-* byte for byte

# which kind of baseline each subcommand needs (manifest.json profile_type); `score` infers it from the normalized directory
NEEDS_PROFILE_TYPE = {"normalize": "assembly", "run-normalize": "assembly", "import-mgnify": "assembly", "run-mgnify": "assembly",
                      "import-metaphlan": "reads", "run-metaphlan": "reads"}
PROFILE_TYPE_WORD = {"assembly": "assembly-based", "reads": "read-based"}
# how a normalized directory records which measurement it came from (sample_meta.json "pipeline")
PIPELINE_PROFILE_TYPE = {"metaphlan": "reads", "mgnify_v5_assembly": "assembly"}
CURRENT_SERIES = {"assembly": DEFAULT_ASSEMBLY_VERSION, "reads": DEFAULT_READS_VERSION}   # bundle series this release was validated with


# --------------------------------------------------------------------------------------------- checks

def warn(msg: str):
    """Print one diagnostic on stderr in the single form the docs promise, `checkgm: ...`.

    That prefix is the one recognition rule given to users and to scripts, so every notice goes through here; per-sample
    progress stays on stdout, where a reader redirecting a log expects it."""
    print(f"checkgm: {msg}", file=sys.stderr, flush=True)


def _load_json_object(p: str, what: str, remedy: str) -> dict:
    """A JSON object read from `p`, or ValueError naming the path.

    The bare text of a parse error ("Expecting property name ...", "Unterminated string starting at") tells a user whose
    archive was truncated nothing they can act on, and a file holding a list or a number parses without complaint only to
    fail later inside a comprehension, so both are turned here into one sentence that names the file and what to do."""
    try:
        m = json.load(open(p))
    except json.JSONDecodeError as e:
        raise ValueError(f"{what} is not readable JSON: {p} ({e.msg}: line {e.lineno}, column {e.colno}); {remedy}") from None
    if not isinstance(m, dict):
        raise ValueError(f"{what} is not a JSON object: {p}; {remedy}")
    return m


def _manifest(bundle: str) -> dict:
    """The bundle's manifest.json, reported with its path when it cannot be read (see _load_json_object)."""
    return _load_json_object(os.path.join(bundle, "manifest.json"), "bundle manifest",
                             "the bundle archive may be incomplete, check its sha256 against the one published with it")


def bundle_id_of(bundle: str) -> str:
    """The bundle's own id, falling back to its directory name for a manifest that omits the field.

    Every import step stamps this beside the pipeline name in sample_meta.json, so that `score` can say which bundle
    normalized the tables it is handed rather than assume it was the one now named."""
    return str(_manifest(bundle).get("bundle_id") or os.path.basename(os.path.normpath(bundle)))


def _sample_meta(out: str) -> dict:
    """OUT/sample_meta.json as {sample_id: metadata}, or {} when the directory has none.

    The guards are the manifest's: a half-copied directory and one whose sample_meta.json was written by hand are both
    common enough that the file is never trusted to be an object of objects."""
    p = os.path.join(out, "sample_meta.json")
    if not os.path.isfile(p):
        return {}
    return _load_json_object(p, "sample_meta.json", "it maps sample id -> metadata, so re-run the import step that wrote the directory")


def check_bundle_type(bundle: str, cmd: str, required: str | None = None) -> dict:
    """Read the bundle manifest and raise ValueError when its profile_type does not fit the subcommand.

    `required` overrides the table NEEDS_PROFILE_TYPE (used by `score`, which infers the needed type from the normalized
    data). Returns the manifest so callers do not read it twice. A command for which neither a caller nor the table
    names a type is an internal error rather than a command exempt from the check."""
    m = _manifest(bundle)
    need = required or NEEDS_PROFILE_TYPE.get(cmd)
    if not need:
        raise ValueError(f"no bundle profile_type is registered for '{cmd}' (internal error)")
    have = m.get("profile_type", "")
    warn_superseded(m)
    if have != need:
        art = {"assembly": "an", "reads": "a"}
        raise ValueError(f"bundle '{m.get('bundle_id', os.path.basename(bundle))}' is {art.get(have, 'a')} {PROFILE_TYPE_WORD.get(have, have)} baseline "
                         f"(profile_type '{have}'), but '{cmd}' needs {art[need]} {PROFILE_TYPE_WORD[need]} one "
                         f"({'gut-reads-*' if need == 'reads' else 'gut-assembly-*'}); run 'checkgm bundles' to see what is installed")
    return m


def bundle_series(bundle_id: str) -> str | None:
    """The series part of a bundle id (`gut-reads-adult-global-v0.2-lenient` -> '0.2'), or None when the id has no v<major.minor>."""
    m = re.search(r"-v(\d+\.\d+)(?:-|$)", bundle_id or "")
    return m.group(1) if m else None


def warn_superseded(manifest: dict):
    """Print one stderr line when the bundle belongs to an older series than the one this release was validated with.

    Scoring goes ahead (exit status unchanged); the README lists the superseded series."""
    cur = superseded_by(manifest)
    if cur:
        warn(f"bundle series v{bundle_series(manifest.get('bundle_id', ''))} is superseded by v{cur} "
             f"(bundle '{manifest.get('bundle_id')}'); results from it are not the validated ones")


def superseded_by(manifest: dict) -> str | None:
    """The current series when the bundle belongs to an older one, else None (also None when the id carries no version)."""
    have = bundle_series(manifest.get("bundle_id", "")); cur = CURRENT_SERIES.get(manifest.get("profile_type", ""))
    if have and cur and tuple(int(x) for x in have.split(".")) < tuple(int(x) for x in cur.split(".")):
        return cur
    return None


def load_calibration(a):
    """Check and load `--calibration` (or return None) before any work is done, so a wrong path fails before --out is written."""
    if not getattr(a, "calibration", None):
        return None
    check_exists([a.calibration], "--calibration", "dir")
    check_exists([os.path.join(a.calibration, "calibration.json"), os.path.join(a.calibration, "features.parquet")], "--calibration")
    from checkgm.calibrate import Calibration
    return Calibration.load(a.calibration)


def check_exists(paths, what: str, kind: str = "file"):
    """Raise FileNotFoundError naming the first missing path. kind is 'file' or 'dir'."""
    test = os.path.isdir if kind == "dir" else os.path.isfile
    for p in paths:
        if not test(p):
            raise FileNotFoundError(f"{what} {'directory' if kind == 'dir' else 'file'} not found: {p}")


def read_reads_tsv(path: str) -> dict:
    """Read a two-column sample<TAB>reads table into {sample: reads}.

    A header line is accepted and skipped. Rows whose read count is not a number are reported on stderr and ignored
    rather than aborting the run (the sample then gets depth band 'unknown')."""
    check_exists([path], "--reads-tsv")
    rt = pd.read_csv(path, sep="\t", header=None, names=["sample", "reads"], dtype=str, comment=None, skip_blank_lines=True)
    rt["sample"] = rt["sample"].astype(str).str.strip()
    num = pd.to_numeric(rt["reads"].astype(str).str.replace(",", "", regex=False).str.strip(), errors="coerce")
    bad = rt[num.isna()]
    if len(bad):
        first = bad.index[0]
        looks_like_header = first == 0 and bad.iloc[0]["sample"].lower() in ("sample", "sample_id", "sampleid", "id", "#sample", "#sampleid")
        rest = bad.iloc[1:] if looks_like_header else bad
        if len(rest):
            shown = ", ".join(f"line {i + 1} ({s!r})" for i, s in zip(rest.index[:5], rest["sample"].iloc[:5]))
            warn(f"{len(rest)} row(s) of {path} have no numeric read count and are ignored: {shown}{', ...' if len(rest) > 5 else ''}")
    keep = num.notna()
    return {s: int(r) for s, r in zip(rt.loc[keep, "sample"], num[keep])}


# --------------------------------------------------------------------------------------------- normalized tables on disk

def save_normalized(norm: dict, out: str):
    """Write {sample_id: normalized dict} to OUT/normalized/<sample>.parquet and OUT/qc.tsv (floats at full precision)."""
    nd = os.path.join(out, "normalized"); os.makedirs(nd, exist_ok=True)
    qcs = []
    for sid, d in norm.items():
        parts = []
        for layer in LAYERS:
            if layer in d and len(d[layer]):
                t = d[layer].copy(); t.insert(0, "layer", layer); parts.append(t)
        if parts:
            pd.concat(parts, ignore_index=True).to_parquet(os.path.join(nd, f"{sid}.parquet"), index=False)
        q = dict(d["qc"]); q["rejections"] = "; ".join(f"{c}: {m}" for c, m in q["rejections"]); qcs.append(q)
    pd.DataFrame(qcs).to_csv(os.path.join(out, "qc.tsv"), sep="\t", index=False, float_format=QC_FLOAT_FORMAT)


def check_normalized_dir(out: str):
    """Raise FileNotFoundError when `out` is not an import step's output directory, which qc.tsv is the marker of."""
    if not os.path.isfile(os.path.join(out, "qc.tsv")):
        raise FileNotFoundError(f"--normalized directory has no qc.tsv: {out} (run 'checkgm import-metaphlan', 'import-mgnify' or 'normalize' first)")


def load_normalized(out: str) -> dict:
    """Read back what save_normalized wrote: {sample_id: {layer: DataFrame, ..., "qc": dict}}."""
    nd = os.path.join(out, "normalized"); qcp = os.path.join(out, "qc.tsv")
    check_normalized_dir(out)
    qc = pd.read_csv(qcp, sep="\t", keep_default_na=False, float_precision="round_trip")   # exact: import + score == run-* byte for byte
    norm = {}
    for q in qc.to_dict("records"):
        sid = str(q["sample_id"]); p = os.path.join(nd, f"{sid}.parquet")
        d = {}
        if os.path.exists(p):
            t = pd.read_parquet(p)
            d = {layer: g.drop(columns="layer").reset_index(drop=True) for layer, g in t.groupby("layer")}
        q["rejections"] = [tuple(x.split(": ", 1)) for x in str(q["rejections"]).split("; ") if x]
        for k in ("total_length_all_contigs", "genome_equivalents_cov", "mapped_fraction_family", "n50_cds_contigs"):
            q[k] = float(q[k]) if q.get(k) not in ("", None) else float("nan")
        d["qc"] = q; norm[sid] = d
    return norm


def normalized_profile_type(out: str) -> str | None:
    """'reads' or 'assembly' according to OUT/sample_meta.json, or None when the directory does not say.

    A directory without the file, and one whose entries disagree, both answer None; the caller turns each into its own
    message, the two being corrected differently."""
    def pipeline_name(v: dict) -> str:
        """Pipeline name of one sample_meta.json entry: import steps write a string, a query's sample.json an object {name, version, ...}."""
        pl = v.get("pipeline", "")
        return str(pl.get("name", "")) if isinstance(pl, dict) else str(pl)
    kinds = {PIPELINE_PROFILE_TYPE.get(pipeline_name(v), "assembly") for v in _sample_meta(out).values() if isinstance(v, dict)}
    return kinds.pop() if len(kinds) == 1 else None


def warn_other_bundle(out: str, manifest: dict):
    """Warn when the tables were normalized against a bundle other than the one they are about to be scored against.

    Normalizing once and scoring against several bundles of one series is an advertised use, so the difference is
    reported rather than refused. It is worth reporting because each bundle carries its own CLR basis: a feature that
    left the basis between the two is dropped from the comparison, while the values of the features that remain stay
    centred over the basis they were computed on, and no column of the output records either fact. A directory written
    before this stamp existed carries no bundle id and is passed over in silence."""
    have = str(manifest.get("bundle_id") or "")
    other = sorted({str(v["bundle_id"]) for v in _sample_meta(out).values()
                    if isinstance(v, dict) and v.get("bundle_id") and str(v["bundle_id"]) != have})
    if other:
        warn(f"the tables in {out} were normalized against {' and '.join(other)}, not '{have}'; they are scored as they "
             f"stand, which is what comparing bundles of one series asks for, but every layer's CLR stays centred over "
             f"the basis it was computed on and features outside this bundle's basis are left out of the comparison")


def required_profile_type(out: str) -> str:
    """The kind of bundle the tables in a normalized directory must be scored against, from OUT/sample_meta.json.

    Nothing else in the directory records which measurement produced the tables, so a directory that does not say is
    refused rather than scored against whatever bundle was named: were the check merely skipped, MetaPhlAn species CLR
    values would be read as percentiles of an assembly reference, which is the one comparison the bundle check exists to
    prevent. The file is easily lost, since only qc.tsv is needed to read the tables back, and the corrective action
    differs between a directory that was never written by an import step and one that holds two measurements at once."""
    check_normalized_dir(out)
    t = normalized_profile_type(out)
    if t is not None:
        return t
    p = os.path.join(out, "sample_meta.json")
    if not os.path.isfile(p):
        raise FileNotFoundError(f"--normalized directory has no sample_meta.json, so which measurement it came from cannot be told: "
                                f"{out} (point --normalized at the directory the import step wrote; 'normalize', 'import-mgnify' and "
                                f"'import-metaphlan' all write it beside qc.tsv)")
    raise ValueError(f"--normalized directory does not name one measurement: {out} (sample_meta.json names none or several, and a "
                     f"bundle is a baseline for exactly one; import and score each measurement separately)")


# --------------------------------------------------------------------------------------------- import steps

def cmd_normalize(a):
    """`checkgm normalize`: query directories (docs/query_format.md) -> normalized tables under --out."""
    check_exists(a.query, "--query", "dir")
    bid = bundle_id_of(a.bundle)
    spec = BundleSpec(a.bundle); bb = Backbone(spec.backbone_path); norm = {}; meta = {}
    for qd in a.query:
        t = read_query_dir(qd); sid = t["meta"].get("sample_id") or os.path.basename(os.path.normpath(qd))
        norm[sid] = normalize(sid, t["contigs"], t["cds"], t["cds_tax"], t["modules"], spec, bb); meta[sid] = dict(t["meta"], bundle_id=bid)
        print(sid, norm[sid]["qc"]["status"], norm[sid]["qc"]["quality_band"], flush=True)
    os.makedirs(a.out, exist_ok=True)
    save_normalized(norm, a.out); json.dump(meta, open(os.path.join(a.out, "sample_meta.json"), "w"), indent=1)
    return norm


def cmd_import_mgnify(a):
    """`checkgm import-mgnify`: MGnify v5 assembly analysis folders -> normalized tables under --out."""
    check_exists(a.analysis_dir, "--analysis-dir", "dir")
    spec = BundleSpec(a.bundle); bb = Backbone(spec.backbone_path); norm = {}
    for d in a.analysis_dir:
        sid = os.path.basename(os.path.normpath(d)); t0 = time.time()
        norm[sid] = normalize_mgnify_v5(sid, d, spec, bb)
        print(sid, norm[sid]["qc"]["status"], norm[sid]["qc"]["quality_band"], f"{time.time()-t0:.0f}s", flush=True)
    os.makedirs(a.out, exist_ok=True)
    save_normalized(norm, a.out); json.dump({s: {"pipeline": "mgnify_v5_assembly", "bundle_id": bundle_id_of(a.bundle)} for s in norm}, open(os.path.join(a.out, "sample_meta.json"), "w"), indent=1)
    return norm


def cmd_import_metaphlan(a):
    """`checkgm import-metaphlan`: MetaPhlAn 3 profiles (and optional HUMAnN 3 tables) -> normalized tables under --out."""
    from checkgm.normalize_reads import ReadsBundleSpec, read_metaphlan, normalize_profile
    check_exists(a.input, "--input")
    for f in (getattr(a, "humann_genefamilies", None) or []) + (getattr(a, "humann_pathabundance", None) or []):
        check_exists([f], "HUMAnN")
    reads_of = read_reads_tsv(a.reads_tsv) if getattr(a, "reads_tsv", None) else {}
    spec = ReadsBundleSpec(a.bundle); norm = {}
    for path in a.input:
        for sid, d in read_metaphlan(path).items():
            n = normalize_profile(sid, d["table"], reads_of.get(sid, d["reads"]), spec); norm[sid] = n
            print(sid, n["qc"]["status"], n["qc"]["quality_band"], f"species rows {n['qc']['n_species_rows']}, family mapped {n['qc']['mapped_fraction_family']:.3f}", flush=True)
    # optional function layers from HUMAnN 3 tables (bundles that carry ko_humann / pathway_humann)
    from checkgm.normalize_reads import read_humann, normalize_function
    for layer, files in (("ko_humann", getattr(a, "humann_genefamilies", None)), ("pathway_humann", getattr(a, "humann_pathabundance", None))):
        if not files:
            continue
        if layer not in spec.fbasis:
            warn(f"bundle has no {layer} layer; {len(files)} HUMAnN file(s) ignored"); continue
        tabs, src = {}, {}
        for path in files:
            for sid, v in read_humann(path).items():
                tabs[sid] = v; src[sid] = path
        if len(tabs) == 1 and len(norm) == 1 and next(iter(tabs)) not in norm:
            # One sample on each side: pair them whatever the column is called, since HUMAnN's column names are derived
            # from file names and rarely survive a rename. The pairing is announced, because a table belonging to another
            # sample would otherwise be grafted onto this profile's taxonomy without leaving a trace.
            fsid = next(iter(tabs)); sid = next(iter(norm))
            warn(f"{layer}: HUMAnN sample '{fsid}' ({src[fsid]}) paired with profile '{sid}' (one sample on each side, so the "
                 f"name in the table is ignored); check that the two belong to the same sample")
            tabs = {sid: tabs[fsid]}; src = {sid: src[fsid]}
        for sid, v in tabs.items():
            if sid not in norm:
                warn(f"{layer}: sample {sid} ({src[sid]}) has no MetaPhlAn profile in this run; skipped"); continue
            t, q = normalize_function(v, layer, spec); norm[sid][layer] = t; norm[sid]["qc"].update(q)
            share = q[layer + "_share_in_basis"]
            print(sid, layer, f"detected {int((t.proportion_mapped > 0).sum())} of {len(t)} basis features, share in basis {share:.3f}", flush=True)
            if not share > 0:   # also catches NaN; a layer that mapped nothing is a user error, not a result
                warn(f"{layer}: nothing in {src[sid]} maps to the bundle's {layer} basis for sample {sid}, so the layer carries no "
                     f"information (wrong kind of table, or a gene-family table that is neither UniRef90 nor KO?)")
            elif t["clr"].isna().any():
                # Too few features for the multiplicative zero replacement, which then leaves the CLR undefined. The
                # layer is still written, but a NaN column scores as nothing at all and must not pass unremarked.
                warn(f"{layer}: the CLR of {int(t['clr'].isna().sum())} of {len(t)} basis features is undefined for sample {sid}, too "
                     f"few of them being detected in {src[sid]}; this layer's scores are not usable")
    os.makedirs(a.out, exist_ok=True)
    save_normalized(norm, a.out); json.dump({s: {"pipeline": "metaphlan", "bundle_id": bundle_id_of(a.bundle)} for s in norm}, open(os.path.join(a.out, "sample_meta.json"), "w"), indent=1)
    return norm


# --------------------------------------------------------------------------------------------- report

def taxon_names(bundle: Bundle, taxids) -> dict:
    """{taxid string: scientific name} for the given feature ids, from the bundle's NCBI backbone (empty if the bundle has none)."""
    p = os.path.join(bundle.path, "normalization", "backbone", "ncbi_taxonomy.parquet")
    want = {str(t) for t in taxids}
    if not want or not os.path.isfile(p):
        return {}
    tax = pd.read_parquet(p, columns=["taxid", "name"])
    tax = tax[tax["taxid"].astype(str).isin(want)]
    return dict(zip(tax["taxid"].astype(str), tax["name"]))


def write_report(r: dict, bundle: Bundle, out: str, calibration_id: str):
    """Write scores.parquet, summaries.tsv, rejections.tsv and the readable report.md for a score() result."""
    os.makedirs(out, exist_ok=True)
    r["scores"].to_parquet(os.path.join(out, "scores.parquet"), index=False)
    r["summaries"].to_csv(os.path.join(out, "summaries.tsv"), sep="\t", index=False)
    r["rejections"].to_csv(os.path.join(out, "rejections.tsv"), sep="\t", index=False)
    S = r["summaries"]; m = bundle.manifest; ptype = m.get("profile_type", "")
    sp = m.get("source_pipeline", {}); bid = m.get("bundle_id", "")
    measured = {"gut-assembly-adult-global": "MGnify pipeline v5.0", "gut-reads-adult-global": "curatedMetagenomicData 3: MetaPhlAn 3 + HUMAnN 3"}
    measured = next((v for k, v in measured.items() if bid.startswith(k)), f"{sp.get('name', '')} {sp.get('version', '')}".strip())
    kind = {"assembly": "gut assembly analyses", "reads": "gut read-based profiles"}.get(ptype, "gut samples") + (f" ({measured})" if measured else "")
    remeasure = {"assembly": "re-assemblies of the same reads", "reads": "re-profiling of the same reads"}.get(ptype, "re-measurements of the same reads")
    L = [f"# checkgm report — {m['bundle_id']}\n", f"Reference: {m['n_samples']:,} healthy adult {kind}, {m['n_studies']} studies. Calibration: {calibration_id}. "
         f"Scored {S['analysis_id'].nunique() if not S.empty else 0} sample(s); {len(r['rejections'])} not scored (rejections.tsv).\n",
         "Read this first: a healthy adult gut sample measured like the reference has about 5% of features outside the 2.5–97.5 percentile band. "
         "The sample-level excess test asks whether a sample has more than that (q ≤ 0.05 after BH across samples within a layer). "
         f"Individual feature calls (low/high) are descriptive; binary calls are less reproducible across {remeasure} than percentiles.\n"]
    if not S.empty:
        sc = r["scores"]
        is_tax = sc["layer"].str.startswith("taxonomy_")
        names = taxon_names(bundle, sc.loc[is_tax & sc["call"].isin(["low", "high", "expected_but_missing"]), "feature_id"].unique())

        def with_name(g: pd.DataFrame) -> pd.DataFrame:
            """Insert a `name` column (NCBI name for taxonomy layers, empty otherwise) after feature_id."""
            g = g.copy(); g.insert(g.columns.get_loc("feature_id") + 1, "name", [names.get(str(f), "") if l.startswith("taxonomy_") else "" for f, l in zip(g["feature_id"], g["layer"])])
            return g

        piv = S.pivot_table(index="analysis_id", columns="layer", values="frac_outside_raw")
        L += ["## Fraction of assessed features outside the reference band\n", piv.round(3).to_markdown(floatfmt=".3f"), ""]
        sig = S[S["excess_outside_q"] <= 0.05]
        L += ["## Samples with a significant excess (q ≤ 0.05)\n", (sig[["analysis_id", "layer", "n_assessed", "n_low_raw", "n_high_raw", "excess_outside_ratio", "excess_outside_q"]].round(3).to_markdown(index=False) if len(sig) else "none"), ""]
        if "core_distance_pct" in S.columns:
            ld = S.drop_duplicates("analysis_id")[[c for c in ["analysis_id", "quality_band", "genome_equivalents_cov", "core_distance_pct", "neighbor_studies"] if c in S.columns]]
            ld = ld.dropna(axis=1, how="all")   # e.g. genome_equivalents_cov is not defined for read-based profiles
            L += ["## Landscape (family CLR PCA of the reference)\n", "core_distance_pct = percentile of the sample's distance from the reference centroid among reference samples.\n", ld.round(2).to_markdown(index=False), ""]
        top = sc[sc["call"].isin(["low", "high"])].copy()
        if len(top):
            top["dev"] = (top["percentile"] - 50).abs()
            L += ["## Most extreme features per sample (top 10 by percentile distance)\n", "`name` is the NCBI name of a taxon (taxonomy layers only).\n"]
            for sid, g in top.sort_values("dev", ascending=False).groupby("analysis_id"):
                L += [f"### {sid}\n", with_name(g.head(10))[["layer", "feature_id", "name", "band", "value", "percentile", "call"]].round(3).to_markdown(index=False), ""]
        miss = sc[sc["call"] == "expected_but_missing"]
        ebm = miss.groupby(["analysis_id", "layer"]).size().unstack(fill_value=0)
        L += ["## Expected-but-missing features (reference prevalence ≥ 0.70 in the sample's band, not detected)\n", ebm.to_markdown() if len(ebm) else "none", ""]
        if len(miss):
            L += ["Listed by name for taxonomy layers (up to 40 per layer; gene and pathway layers by identifier, first 40):\n"]
            for (sid, layer), g in miss.groupby(["analysis_id", "layer"]):
                ids = [f"{f} ({names[str(f)]})" if str(f) in names else str(f) for f in g["feature_id"]]
                shown = ", ".join(sorted(ids)[:40]) + (f", ... and {len(ids) - 40} more" if len(ids) > 40 else "")
                L += [f"- **{sid}**, {layer} ({len(ids)}): {shown}"]
            L += [""]
    open(os.path.join(out, "report.md"), "w").write("\n".join(L))


def cmd_score(a, norm=None, cal=None):
    """`checkgm score`: compare normalized tables with the bundle and write the report files to --out.

    norm / cal are passed by the run-* commands, which have already imported the samples and loaded the calibration."""
    if norm is None:
        check_exists([a.normalized], "--normalized", "dir")
        warn_other_bundle(a.normalized, check_bundle_type(a.bundle, "score", required=required_profile_type(a.normalized)))
        cal = load_calibration(a)
        norm = load_normalized(a.normalized)
    bundle = Bundle(a.bundle)
    r = score(norm, bundle, calibration=cal)
    write_report(r, bundle, a.out, "none" if cal is None else cal.id)
    print(json.dumps({"scored": int(r["summaries"]["analysis_id"].nunique()) if not r["summaries"].empty else 0, "rejected": len(r["rejections"]), "rows": len(r["scores"]), "out": a.out}))


def cmd_bundles(a) -> int:
    """`checkgm bundles`: list the bundles found under $CHECKGM_BUNDLES (id, profile type, pool size, path, series).

    This is the command that answers "what do I have and which one should I use", so each row says whether its series is
    the one this release was validated with; the warning would otherwise arrive only once a run is under way. One
    unreadable bundle is reported and passed over rather than aborting the listing, since the bundle a user wants is
    rarely the broken one, and exit 1 is kept for the case where every bundle found was unreadable. A search path that is
    unset, or that holds no bundle, is a question answered rather than an error: one line on stderr and exit 0."""
    var = bundle_dirs_source(); found = list_bundles(); shown = 0
    if not bundle_dirs():
        warn(f"{ENV_VAR} is not set; give --bundle a directory path or export {ENV_VAR}=/dir/with/bundles")
    for b in found:
        try:
            m = _manifest(b)
        except (ValueError, OSError) as e:
            warn(f"{os.path.basename(b)} skipped: {e}"); continue
        cur = superseded_by(m)
        # `current` is a statement about the series this release was validated with, so it can only be made for a profile
        # type that release knows; a bundle of some other kind gets an empty field rather than a reassurance nobody checked.
        known = m.get("profile_type", "") in CURRENT_SERIES and bundle_series(m.get("bundle_id", ""))
        series = f"superseded by v{cur}" if cur else ("current" if known else "")
        print(f"{m.get('bundle_id', os.path.basename(b))}\t{m.get('profile_type', '')}\t"
              f"{m.get('n_samples', '?')} samples, {m.get('n_studies', '?')} studies\t{b}\t{series}")
        shown += 1
    if bundle_dirs() and not found:
        warn(f"no bundles under {var} = {os.pathsep.join(bundle_dirs())}{legacy_env_note()}")
    return 1 if found and not shown else 0


# --------------------------------------------------------------------------------------------- argument parser

H = {  # one sentence per argument, shared by the subcommands that take it (docs/cli.md says the same)
    "bundle": f"bundle directory (one containing manifest.json) or a bundle id looked up under ${ENV_VAR}",
    "out": "output directory; created if needed",
    "calibration": "calibration directory produced with the checkgm.calibrate module, for a measurement pipeline other than the reference's (none is shipped)",
    "input": "MetaPhlAn 3 single-sample profile(s) (#clade_name, NCBI_tax_id, relative_abundance) or a merged table (clade_name, sample columns)",
    "reads_tsv": "two columns sample<TAB>reads processed, used for the depth band when the profile header has no '#N reads processed' line; a header line is skipped",
    "humann_genefamilies": "HUMAnN 3 gene-family table(s) (UniRef90 or KO rows); adds the ko_humann layer when the bundle has it, "
                           "and with one profile and one table the two are paired whatever the table's column is called",
    "humann_pathabundance": "HUMAnN 3 pathway-abundance table(s) (MetaCyc rows); adds the pathway_humann layer when the bundle has it, "
                            "paired with the profile as --humann-genefamilies is",
    "analysis_dir": "MGnify pipeline v5 assembly analysis download folder(s), one per sample, named by the MGYA accession",
    "query": f"one or more query directories in the format described in docs/query_format.md ({REPO}): contigs.tsv, cds.tsv, cds_taxonomy.tsv, optional modules.tsv and sample.json",
    "normalized": "directory written by an import step (normalize, import-mgnify or import-metaphlan): normalized/, qc.tsv and "
                  "sample_meta.json, the last naming the measurement the bundle is checked against and the bundle that normalized "
                  "the tables, so all three are needed",
}
# --out means a normalized directory, a report directory or both, according to the subcommand; the opening words stay the same
OUT_H = {"normalized": "output directory for the normalized tables (normalized/, qc.tsv, sample_meta.json); created if needed",
         "report": "output directory for the report (scores.parquet, summaries.tsv, rejections.tsv, report.md), which may be the "
                   "--normalized directory itself, leaving one directory with both halves; created if needed",
         "both": "output directory for both halves, the normalized tables and the report, which this command writes side by side; created if needed"}


def _add(p, *names, out: str | None = None):
    for n in names:
        flag = "--" + n.replace("_", "-")
        if n == "out":
            p.add_argument(flag, required=True, help=OUT_H[out] if out else H[n])
        elif n == "bundle":
            p.add_argument(flag, required=True, help=H[n])
        elif n in ("input", "analysis_dir", "query"):
            p.add_argument(flag, nargs="+", required=True, help=H[n])
        elif n in ("humann_genefamilies", "humann_pathabundance"):
            p.add_argument(flag, nargs="+", help=H[n])
        elif n == "normalized":
            p.add_argument(flag, required=True, help=H[n])
        else:
            p.add_argument(flag, help=H[n])


def build_parser() -> argparse.ArgumentParser:
    """The checkgm argument parser (subcommands with one-sentence help for every argument)."""
    ap = argparse.ArgumentParser(prog="checkgm", description="Score gut metagenomes against a healthy-adult reference baseline measured the same way. "
                                 "Each baseline (bundle) belongs to one measurement: assembly-based (MGnify v5) or read-based (MetaPhlAn 3 / HUMAnN 3).",
                                 epilog="User errors print one line 'checkgm: ...' on stderr and exit 1; samples that fail a normalization floor are listed in rejections.tsv and exit 0.")
    ap.add_argument("--version", action="version", version=f"checkgm {__version__}")
    sub = ap.add_subparsers(dest="cmd", required=True, metavar="COMMAND")
    sub.add_parser("bundles", help=f"list the bundles found under ${ENV_VAR}",
                   description=f"List the bundles found under ${ENV_VAR}: id, profile type (assembly or reads), pool size, path, and whether "
                               f"the bundle's series is the one this release was validated with.")
    p = sub.add_parser("normalize", help="assembly-based: normalize query directories from any assembly pipeline",
                       description="Normalize one query directory per sample (format: docs/query_format.md, see --query below) with the rules "
                                   "stored in an assembly-based bundle. Writes OUT/normalized/<sample>.parquet, OUT/qc.tsv and "
                                   "OUT/sample_meta.json; score them with 'checkgm score', or use 'run-normalize' to do both at once.")
    _add(p, "bundle", "query", "out", out="normalized")
    p = sub.add_parser("import-mgnify", help="assembly-based: normalize MGnify v5 assembly analyses",
                       description="Normalize MGnify pipeline v5 assembly analysis folders (the reference's own measurement, so no calibration) "
                                   "with the rules of an assembly-based bundle. Writes the same files as 'normalize'; score them with "
                                   "'checkgm score', or use 'run-mgnify' to do both at once.")
    _add(p, "bundle", "analysis_dir", "out", out="normalized")
    p = sub.add_parser("import-metaphlan", help="read-based: normalize MetaPhlAn 3 profiles (and HUMAnN 3 tables)",
                       description="Normalize MetaPhlAn 3 profiles, and optionally HUMAnN 3 gene-family and pathway tables, with the rules of a read-based bundle. "
                                   "Species-level rows are placed on the bundle's taxonomy; proportions over the bundle's basis are CLR-transformed. Writes the same "
                                   "files as 'normalize'; score them with 'checkgm score', or use 'run-metaphlan' to do both at once.")
    _add(p, "bundle", "input", "reads_tsv", "humann_genefamilies", "humann_pathabundance", "out", out="normalized")
    p = sub.add_parser("score", help="score normalized tables against a bundle and write the report",
                       description="Compare normalized tables (from normalize, import-mgnify or import-metaphlan) with the bundle. "
                                   "Writes scores.parquet (one row per sample x layer x feature), summaries.tsv (one row per sample x layer with the excess test), "
                                   "rejections.tsv and report.md.")
    _add(p, "bundle", "normalized", "calibration", "out", out="report")
    p = sub.add_parser("run-normalize", help="assembly-based: normalize and score query directories in one step",
                       description="Run 'normalize' and then 'score' into one output directory.")
    _add(p, "bundle", "query", "calibration", "out", out="both")
    p = sub.add_parser("run-mgnify", help="assembly-based: import-mgnify and score in one step",
                       description="Run 'import-mgnify' and then 'score' into one output directory.")
    _add(p, "bundle", "analysis_dir", "calibration", "out", out="both")
    p = sub.add_parser("run-metaphlan", help="read-based: import-metaphlan and score in one step",
                       description="Run 'import-metaphlan' and then 'score' into one output directory.")
    _add(p, "bundle", "input", "reads_tsv", "humann_genefamilies", "humann_pathabundance", "calibration", "out", out="both")
    return ap


def _dispatch(a) -> int:
    if a.cmd == "bundles":
        return cmd_bundles(a)
    a.bundle = resolve_bundle(a.bundle)
    if a.cmd != "score":
        check_bundle_type(a.bundle, a.cmd)
    if a.cmd == "normalize":
        cmd_normalize(a)
    elif a.cmd == "import-mgnify":
        cmd_import_mgnify(a)
    elif a.cmd == "score":
        cmd_score(a)
    elif a.cmd == "import-metaphlan":
        cmd_import_metaphlan(a)
    elif a.cmd == "run-normalize":
        cal = load_calibration(a); a.normalized = a.out; norm = cmd_normalize(a); cmd_score(a, norm, cal)
    elif a.cmd == "run-metaphlan":
        cal = load_calibration(a); a.normalized = a.out; norm = cmd_import_metaphlan(a); cmd_score(a, norm, cal)
    elif a.cmd == "run-mgnify":
        cal = load_calibration(a); a.normalized = a.out; norm = cmd_import_mgnify(a); cmd_score(a, norm, cal)
    else:
        # Unreachable through argparse, but this chain is the single place a new subcommand has to be registered, and
        # silent success is the worst way to learn that it was not.
        raise ValueError(f"unhandled subcommand '{a.cmd}' (internal error)")
    return 0


def main(argv=None) -> int:
    """Entry point of the `checkgm` console script. Returns the exit status (0 ok, 1 user error, 2 bad arguments)."""
    a = build_parser().parse_args(argv)
    try:
        return _dispatch(a)
    except FileNotFoundError as e:   # must stay ahead of OSError, of which it is a subclass with a better message
        msg = e.args[0] if e.args and isinstance(e.args[0], str) and not e.filename else (f"file not found: {e.filename}" if e.filename else str(e))
        print(f"checkgm: {msg}", file=sys.stderr); return 1
    except OSError as e:
        # Whatever the filesystem refuses — an unreadable manifest, a directory where a file was expected, a full disk —
        # is as much a user error as a missing path, and the contract promises one line for all of them. Python's own
        # str() of these already names the path and the reason ("[Errno 13] Permission denied: '/x/manifest.json'").
        print(f"checkgm: {e}", file=sys.stderr); return 1
    except KeyError as e:
        print(f"checkgm: input is missing the expected column or field {e}", file=sys.stderr); return 1
    except ValueError as e:
        print(f"checkgm: {e}", file=sys.stderr); return 1


if __name__ == "__main__":
    sys.exit(main())
