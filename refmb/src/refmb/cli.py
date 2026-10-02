#!/usr/bin/env python3
"""refmb command line.

  refmb normalize --bundle B --query DIR [DIR ...] --out OUT          # assembly-based: query directories -> OUT/normalized/<sample>.parquet + qc.tsv
  refmb import-mgnify --bundle B --analysis-dir DIR [...] --out OUT   # assembly-based: MGnify v5 analysis folders -> the same
  refmb score --bundle B --normalized OUT [--calibration CAL] --out REPORT   # -> scores.parquet, summaries.tsv, rejections.tsv, report.md
  refmb run-mgnify --bundle B --analysis-dir DIR [...] --out REPORT           # import + score in one go
  refmb import-metaphlan --bundle B --input PROFILE [...] [--reads-tsv sample<TAB>reads] --out OUT   # read-based: MetaPhlAn output -> normalized
  refmb run-metaphlan --bundle B --input PROFILE [...] --out REPORT           # import + score against a read-based baseline
        [--humann-genefamilies FILE ...] [--humann-pathabundance FILE ...]   # optional HUMAnN 3 tables -> function layers
  refmb bundles                                                               # list the bundles found under $REFMB_BUNDLES

B is a bundle directory or a bundle id looked up under $REFMB_BUNDLES (colon-separated directories).
Sample id = sample.json sample_id, else the directory name.

User errors (missing file, wrong kind of bundle, unreadable input) are reported on stderr as one line starting with
"refmb:" and the command exits 1. Scientific refusals (a sample that fails a normalization floor) are not errors: the
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

from refmb import __version__
from refmb.normalize import Backbone, BundleSpec, normalize, normalize_mgnify_v5, read_query_dir
from refmb.paths import DEFAULT_ASSEMBLY_VERSION, DEFAULT_READS_VERSION, ENV_VAR, bundle_dirs, list_bundles, resolve_bundle
from refmb.score import Bundle, score

LAYERS = ["taxonomy_family", "taxonomy_genus", "taxonomy_species", "ko_eggnog", "ko_kofam", "pfam", "module", "ko_humann", "pathway_humann"]
QC_FLOAT_FORMAT = "%.17g"   # 17 significant digits round-trip every double exactly, so import + score == run-* byte for byte

# which kind of baseline each subcommand needs (manifest.json profile_type)
NEEDS_PROFILE_TYPE = {"normalize": "assembly", "import-mgnify": "assembly", "run-mgnify": "assembly",
                      "import-metaphlan": "reads", "run-metaphlan": "reads"}
PROFILE_TYPE_WORD = {"assembly": "assembly-based", "reads": "read-based"}
# how a normalized directory records which measurement it came from (sample_meta.json "pipeline")
PIPELINE_PROFILE_TYPE = {"metaphlan": "reads", "mgnify_v5_assembly": "assembly"}
CURRENT_SERIES = {"assembly": DEFAULT_ASSEMBLY_VERSION, "reads": DEFAULT_READS_VERSION}   # bundle series this release was validated with


# --------------------------------------------------------------------------------------------- checks

def _manifest(bundle: str) -> dict:
    return json.load(open(os.path.join(bundle, "manifest.json")))


def check_bundle_type(bundle: str, cmd: str, required: str | None = None) -> dict:
    """Read the bundle manifest and raise ValueError when its profile_type does not fit the subcommand.

    `required` overrides the table NEEDS_PROFILE_TYPE (used by `score`, which infers the needed type from the normalized
    data). Returns the manifest so callers do not read it twice."""
    m = _manifest(bundle)
    need = required or NEEDS_PROFILE_TYPE.get(cmd)
    have = m.get("profile_type", "")
    warn_superseded(m)
    if need and have != need:
        art = {"assembly": "an", "reads": "a"}
        raise ValueError(f"bundle '{m.get('bundle_id', os.path.basename(bundle))}' is {art.get(have, 'a')} {PROFILE_TYPE_WORD.get(have, have)} baseline "
                         f"(profile_type '{have}'), but '{cmd}' needs {art[need]} {PROFILE_TYPE_WORD[need]} one "
                         f"({'gut-reads-*' if need == 'reads' else 'gut-assembly-*'}); run 'refmb bundles' to see what is installed")
    return m


def bundle_series(bundle_id: str) -> str | None:
    """The series part of a bundle id (`gut-reads-adult-global-v0.2-lenient` -> '0.2'), or None when the id has no v<major.minor>."""
    m = re.search(r"-v(\d+\.\d+)(?:-|$)", bundle_id or "")
    return m.group(1) if m else None


def warn_superseded(manifest: dict):
    """Print one stderr line when the bundle belongs to an older series than the one this release was validated with.

    Scoring goes ahead (exit status unchanged); the README lists the superseded series."""
    have = bundle_series(manifest.get("bundle_id", "")); cur = CURRENT_SERIES.get(manifest.get("profile_type", ""))
    if have and cur and tuple(int(x) for x in have.split(".")) < tuple(int(x) for x in cur.split(".")):
        print(f"refmb: bundle series v{have} is superseded by v{cur} (bundle '{manifest.get('bundle_id')}'); results from it are not the validated ones", file=sys.stderr, flush=True)


def load_calibration(a):
    """Check and load `--calibration` (or return None) before any work is done, so a wrong path fails before --out is written."""
    if not getattr(a, "calibration", None):
        return None
    check_exists([a.calibration], "--calibration", "dir")
    check_exists([os.path.join(a.calibration, "calibration.json"), os.path.join(a.calibration, "features.parquet")], "--calibration")
    from refmb.calibrate import Calibration
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
            print(f"refmb: {len(rest)} row(s) of {path} have no numeric read count and are ignored: {shown}"
                  f"{', ...' if len(rest) > 5 else ''}", file=sys.stderr, flush=True)
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


def load_normalized(out: str) -> dict:
    """Read back what save_normalized wrote: {sample_id: {layer: DataFrame, ..., "qc": dict}}."""
    nd = os.path.join(out, "normalized"); qcp = os.path.join(out, "qc.tsv")
    if not os.path.isfile(qcp):
        raise FileNotFoundError(f"--normalized directory has no qc.tsv: {out} (run 'refmb import-metaphlan', 'import-mgnify' or 'normalize' first)")
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
    """'reads' or 'assembly' according to OUT/sample_meta.json, or None when the directory does not say."""
    p = os.path.join(out, "sample_meta.json")
    if not os.path.isfile(p):
        return None
    def pipeline_name(v: dict) -> str:
        """Pipeline name of one sample_meta.json entry: import steps write a string, a query's sample.json an object {name, version, ...}."""
        pl = v.get("pipeline", "")
        return str(pl.get("name", "")) if isinstance(pl, dict) else str(pl)
    kinds = {PIPELINE_PROFILE_TYPE.get(pipeline_name(v), "assembly") for v in json.load(open(p)).values() if isinstance(v, dict)}
    return kinds.pop() if len(kinds) == 1 else None


# --------------------------------------------------------------------------------------------- import steps

def cmd_normalize(a):
    """`refmb normalize`: query directories (docs/query_format.md) -> normalized tables under --out."""
    check_exists(a.query, "--query", "dir")
    spec = BundleSpec(a.bundle); bb = Backbone(spec.backbone_path); norm = {}; meta = {}
    for qd in a.query:
        t = read_query_dir(qd); sid = t["meta"].get("sample_id") or os.path.basename(os.path.normpath(qd))
        norm[sid] = normalize(sid, t["contigs"], t["cds"], t["cds_tax"], t["modules"], spec, bb); meta[sid] = t["meta"]
        print(sid, norm[sid]["qc"]["status"], norm[sid]["qc"]["quality_band"], flush=True)
    os.makedirs(a.out, exist_ok=True)
    save_normalized(norm, a.out); json.dump(meta, open(os.path.join(a.out, "sample_meta.json"), "w"), indent=1)
    return norm


def cmd_import_mgnify(a):
    """`refmb import-mgnify`: MGnify v5 assembly analysis folders -> normalized tables under --out."""
    check_exists(a.analysis_dir, "--analysis-dir", "dir")
    spec = BundleSpec(a.bundle); bb = Backbone(spec.backbone_path); norm = {}
    for d in a.analysis_dir:
        sid = os.path.basename(os.path.normpath(d)); t0 = time.time()
        norm[sid] = normalize_mgnify_v5(sid, d, spec, bb)
        print(sid, norm[sid]["qc"]["status"], norm[sid]["qc"]["quality_band"], f"{time.time()-t0:.0f}s", flush=True)
    os.makedirs(a.out, exist_ok=True)
    save_normalized(norm, a.out); json.dump({s: {"pipeline": "mgnify_v5_assembly"} for s in norm}, open(os.path.join(a.out, "sample_meta.json"), "w"), indent=1)
    return norm


def cmd_import_metaphlan(a):
    """`refmb import-metaphlan`: MetaPhlAn 3 profiles (and optional HUMAnN 3 tables) -> normalized tables under --out."""
    from refmb.normalize_reads import ReadsBundleSpec, read_metaphlan, normalize_profile
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
    from refmb.normalize_reads import read_humann, normalize_function
    for layer, files in (("ko_humann", getattr(a, "humann_genefamilies", None)), ("pathway_humann", getattr(a, "humann_pathabundance", None))):
        if not files:
            continue
        if layer not in spec.fbasis:
            print(f"bundle has no {layer} layer; {len(files)} HUMAnN file(s) ignored", flush=True); continue
        tabs = {}
        for path in files:
            tabs.update(read_humann(path))
        if len(tabs) == 1 and len(norm) == 1 and next(iter(tabs)) not in norm:   # one sample on each side: pair them whatever the names
            tabs = {next(iter(norm)): next(iter(tabs.values()))}
        for sid, v in tabs.items():
            if sid not in norm:
                print(f"{layer}: sample {sid} has no MetaPhlAn profile in this run; skipped", flush=True); continue
            t, q = normalize_function(v, layer, spec); norm[sid][layer] = t; norm[sid]["qc"].update(q)
            print(sid, layer, f"detected {int((t.proportion_mapped > 0).sum())} of {len(t)} basis features, share in basis {q[layer + '_share_in_basis']:.3f}", flush=True)
    os.makedirs(a.out, exist_ok=True)
    save_normalized(norm, a.out); json.dump({s: {"pipeline": "metaphlan"} for s in norm}, open(os.path.join(a.out, "sample_meta.json"), "w"), indent=1)
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
    L = [f"# refmb report — {m['bundle_id']}\n", f"Reference: {m['n_samples']:,} healthy adult {kind}, {m['n_studies']} studies. Calibration: {calibration_id}. "
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
    """`refmb score`: compare normalized tables with the bundle and write the report files to --out.

    norm / cal are passed by the run-* commands, which have already imported the samples and loaded the calibration."""
    if norm is None:
        check_exists([a.normalized], "--normalized", "dir")
        need = normalized_profile_type(a.normalized)
        check_bundle_type(a.bundle, "score", required=need)
        cal = load_calibration(a)
        norm = load_normalized(a.normalized)
    bundle = Bundle(a.bundle)
    r = score(norm, bundle, calibration=cal)
    write_report(r, bundle, a.out, "none" if cal is None else cal.id)
    print(json.dumps({"scored": int(r["summaries"]["analysis_id"].nunique()) if not r["summaries"].empty else 0, "rejected": len(r["rejections"]), "rows": len(r["scores"]), "out": a.out}))


def cmd_bundles(a):
    """`refmb bundles`: list the bundles found under $REFMB_BUNDLES (id, profile type, pool size, path)."""
    found = list_bundles()
    if not bundle_dirs():
        print(f"{ENV_VAR} is not set; give --bundle a directory path or export {ENV_VAR}=/dir/with/bundles", file=sys.stderr)
    for b in found:
        m = _manifest(b)
        print(f"{m['bundle_id']}\t{m.get('profile_type', '')}\t{m['n_samples']} samples, {m['n_studies']} studies\t{b}")
    if bundle_dirs() and not found:
        print(f"no bundles under {ENV_VAR} = {os.pathsep.join(bundle_dirs())}", file=sys.stderr)


# --------------------------------------------------------------------------------------------- argument parser

H = {  # one sentence per argument, shared by the subcommands that take it (docs/cli.md says the same)
    "bundle": f"bundle directory (one containing manifest.json) or a bundle id looked up under ${ENV_VAR}",
    "out": "output directory; created if needed",
    "calibration": "calibration directory produced with the refmb.calibrate module, for a measurement pipeline other than the reference's (none is shipped)",
    "input": "MetaPhlAn 3 single-sample profile(s) (#clade_name, NCBI_tax_id, relative_abundance) or a merged table (clade_name, sample columns)",
    "reads_tsv": "two columns sample<TAB>reads processed, used for the depth band when the profile header has no '#N reads processed' line; a header line is skipped",
    "humann_genefamilies": "HUMAnN 3 gene-family table(s) (UniRef90 or KO rows); adds the ko_humann layer when the bundle has it",
    "humann_pathabundance": "HUMAnN 3 pathway-abundance table(s) (MetaCyc rows); adds the pathway_humann layer when the bundle has it",
    "analysis_dir": "MGnify pipeline v5 assembly analysis download folder(s), one per sample, named by the MGYA accession",
    "query": "query directory(ies) in the format of docs/query_format.md (contigs.tsv, cds.tsv, cds_taxonomy.tsv, optional modules.tsv and sample.json)",
    "normalized": "directory written by an import step (normalize, import-mgnify or import-metaphlan): normalized/, qc.tsv, sample_meta.json",
}


def _add(p, *names):
    for n in names:
        flag = "--" + n.replace("_", "-")
        if n in ("bundle", "out"):
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
    """The refmb argument parser (subcommands with one-sentence help for every argument)."""
    ap = argparse.ArgumentParser(prog="refmb", description="Score gut metagenomes against a healthy-adult reference baseline measured the same way. "
                                 "Each baseline (bundle) belongs to one measurement: assembly-based (MGnify v5) or read-based (MetaPhlAn 3 / HUMAnN 3).",
                                 epilog="User errors print one line 'refmb: ...' on stderr and exit 1; samples that fail a normalization floor are listed in rejections.tsv and exit 0.")
    ap.add_argument("--version", action="version", version=f"refmb {__version__}")
    sub = ap.add_subparsers(dest="cmd", required=True, metavar="COMMAND")
    sub.add_parser("bundles", help=f"list the bundles found under ${ENV_VAR}",
                   description=f"List the bundles found under ${ENV_VAR}: id, profile type (assembly or reads), pool size and path.")
    p = sub.add_parser("normalize", help="assembly-based: normalize query directories from any assembly pipeline",
                       description="Normalize one query directory per sample (docs/query_format.md) with the rules stored in an assembly-based bundle. "
                                   "Writes OUT/normalized/<sample>.parquet, OUT/qc.tsv and OUT/sample_meta.json; score them with 'refmb score'.")
    _add(p, "bundle", "query", "out")
    p = sub.add_parser("import-mgnify", help="assembly-based: normalize MGnify v5 assembly analyses",
                       description="Normalize MGnify pipeline v5 assembly analysis folders (the reference's own measurement, so no calibration) "
                                   "with the rules of an assembly-based bundle. Writes the same files as 'normalize'.")
    _add(p, "bundle", "analysis_dir", "out")
    p = sub.add_parser("import-metaphlan", help="read-based: normalize MetaPhlAn 3 profiles (and HUMAnN 3 tables)",
                       description="Normalize MetaPhlAn 3 profiles, and optionally HUMAnN 3 gene-family and pathway tables, with the rules of a read-based bundle. "
                                   "Species-level rows are placed on the bundle's taxonomy; proportions over the bundle's basis are CLR-transformed. Writes the same files as 'normalize'.")
    _add(p, "bundle", "input", "reads_tsv", "humann_genefamilies", "humann_pathabundance", "out")
    p = sub.add_parser("score", help="score normalized tables against a bundle and write the report",
                       description="Compare normalized tables (from normalize, import-mgnify or import-metaphlan) with the bundle. "
                                   "Writes scores.parquet (one row per sample x layer x feature), summaries.tsv (one row per sample x layer with the excess test), "
                                   "rejections.tsv and report.md.")
    _add(p, "bundle", "normalized", "calibration", "out")
    p = sub.add_parser("run-mgnify", help="assembly-based: import-mgnify and score in one step",
                       description="Run 'import-mgnify' and then 'score' into one output directory.")
    _add(p, "bundle", "analysis_dir", "calibration", "out")
    p = sub.add_parser("run-metaphlan", help="read-based: import-metaphlan and score in one step",
                       description="Run 'import-metaphlan' and then 'score' into one output directory.")
    _add(p, "bundle", "input", "reads_tsv", "humann_genefamilies", "humann_pathabundance", "calibration", "out")
    return ap


def _dispatch(a) -> int:
    if a.cmd == "bundles":
        cmd_bundles(a); return 0
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
    elif a.cmd == "run-metaphlan":
        cal = load_calibration(a); a.normalized = a.out; norm = cmd_import_metaphlan(a); cmd_score(a, norm, cal)
    elif a.cmd == "run-mgnify":
        cal = load_calibration(a); a.normalized = a.out; norm = cmd_import_mgnify(a); cmd_score(a, norm, cal)
    return 0


def main(argv=None) -> int:
    """Entry point of the `refmb` console script. Returns the exit status (0 ok, 1 user error, 2 bad arguments)."""
    a = build_parser().parse_args(argv)
    try:
        return _dispatch(a)
    except FileNotFoundError as e:
        msg = e.args[0] if e.args and isinstance(e.args[0], str) and not e.filename else (f"file not found: {e.filename}" if e.filename else str(e))
        print(f"refmb: {msg}", file=sys.stderr); return 1
    except KeyError as e:
        print(f"refmb: input is missing the expected column or field {e}", file=sys.stderr); return 1
    except ValueError as e:
        print(f"refmb: {e}", file=sys.stderr); return 1


if __name__ == "__main__":
    sys.exit(main())
