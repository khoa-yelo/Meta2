#!/usr/bin/env python3
"""refmb command line.

  refmb normalize --bundle B --query DIR [DIR ...] --out OUT      # query-format dirs -> OUT/normalized/<sample>.parquet + qc.tsv
  refmb import-mgnify --bundle B --analysis-dir DIR [...] --out OUT  # MGnify v5 folders -> same (tier A reference measurement)
  refmb score --bundle B --normalized OUT [--calibration CAL] --out REPORT   # -> scores.parquet, summaries.tsv, rejections.tsv, report.md
  refmb run-mgnify --bundle B --analysis-dir DIR [...] --out REPORT           # import + score in one go
  refmb import-metaphlan --bundle B --input PROFILE [...] [--reads-tsv sample<TAB>reads] --out OUT   # pipeline B: MetaPhlAn output -> normalized
  refmb run-metaphlan --bundle B --input PROFILE [...] --out REPORT           # import + score against a read-based bundle
        [--humann-genefamilies FILE ...] [--humann-pathabundance FILE ...]   # optional HUMAnN 3 tables -> function layers
  refmb bundles                                                               # list the bundles found under $REFMB_BUNDLES

B is a bundle directory or a bundle id looked up under $REFMB_BUNDLES (colon-separated directories).
Sample id = sample.json sample_id, else the directory name.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import pandas as pd

from refmb import __version__
from refmb.normalize import Backbone, BundleSpec, normalize, normalize_mgnify_v5, read_query_dir
from refmb.paths import ENV_VAR, bundle_dirs, list_bundles, resolve_bundle
from refmb.score import Bundle, score

LAYERS = ["taxonomy_family", "taxonomy_genus", "taxonomy_species", "ko_eggnog", "ko_kofam", "pfam", "module", "ko_humann", "pathway_humann"]


def save_normalized(norm: dict, out: str):
    nd = os.path.join(out, "normalized"); os.makedirs(nd, exist_ok=True)
    qcs = []
    for sid, d in norm.items():
        parts = []
        for layer in LAYERS:
            if layer in d and len(d[layer]):
                t = d[layer].copy(); t.insert(0, "layer", layer); parts.append(t)
        pd.concat(parts, ignore_index=True).to_parquet(os.path.join(nd, f"{sid}.parquet"), index=False)
        q = dict(d["qc"]); q["rejections"] = "; ".join(f"{c}: {m}" for c, m in q["rejections"]); qcs.append(q)
    pd.DataFrame(qcs).to_csv(os.path.join(out, "qc.tsv"), sep="\t", index=False)


def load_normalized(out: str) -> dict:
    nd = os.path.join(out, "normalized"); qc = pd.read_csv(os.path.join(out, "qc.tsv"), sep="\t", keep_default_na=False)
    norm = {}
    for q in qc.to_dict("records"):
        sid = str(q["sample_id"]); t = pd.read_parquet(os.path.join(nd, f"{sid}.parquet"))
        d = {layer: g.drop(columns="layer").reset_index(drop=True) for layer, g in t.groupby("layer")}
        q["rejections"] = [tuple(x.split(": ", 1)) for x in str(q["rejections"]).split("; ") if x]
        for k in ("total_length_all_contigs", "genome_equivalents_cov", "mapped_fraction_family", "n50_cds_contigs"):
            q[k] = float(q[k]) if q.get(k) not in ("", None) else float("nan")
        d["qc"] = q; norm[sid] = d
    return norm


def cmd_normalize(a):
    spec = BundleSpec(a.bundle); bb = Backbone(spec.backbone_path); norm = {}; meta = {}
    for qd in a.query:
        t = read_query_dir(qd); sid = t["meta"].get("sample_id") or os.path.basename(os.path.normpath(qd))
        norm[sid] = normalize(sid, t["contigs"], t["cds"], t["cds_tax"], t["modules"], spec, bb); meta[sid] = t["meta"]
        print(sid, norm[sid]["qc"]["status"], norm[sid]["qc"]["quality_band"], flush=True)
    save_normalized(norm, a.out); json.dump(meta, open(os.path.join(a.out, "sample_meta.json"), "w"), indent=1)
    return norm


def cmd_import_mgnify(a):
    spec = BundleSpec(a.bundle); bb = Backbone(spec.backbone_path); norm = {}
    for d in a.analysis_dir:
        sid = os.path.basename(os.path.normpath(d)); t0 = time.time()
        norm[sid] = normalize_mgnify_v5(sid, d, spec, bb)
        print(sid, norm[sid]["qc"]["status"], norm[sid]["qc"]["quality_band"], f"{time.time()-t0:.0f}s", flush=True)
    save_normalized(norm, a.out); json.dump({s: {"pipeline": "mgnify_v5_assembly"} for s in norm}, open(os.path.join(a.out, "sample_meta.json"), "w"), indent=1)
    return norm


def cmd_import_metaphlan(a):
    from refmb.normalize_reads import ReadsBundleSpec, read_metaphlan, normalize_profile
    spec = ReadsBundleSpec(a.bundle); norm = {}
    reads_of = {}
    if getattr(a, "reads_tsv", None):
        rt = pd.read_csv(a.reads_tsv, sep="\t", header=None, names=["sample", "reads"], dtype={"sample": str}); reads_of = dict(zip(rt["sample"], rt["reads"]))
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
    save_normalized(norm, a.out); json.dump({s: {"pipeline": "metaphlan"} for s in norm}, open(os.path.join(a.out, "sample_meta.json"), "w"), indent=1)
    return norm


def write_report(r: dict, bundle: Bundle, out: str, calibration_id: str):
    os.makedirs(out, exist_ok=True)
    r["scores"].to_parquet(os.path.join(out, "scores.parquet"), index=False)
    r["summaries"].to_csv(os.path.join(out, "summaries.tsv"), sep="\t", index=False)
    r["rejections"].to_csv(os.path.join(out, "rejections.tsv"), sep="\t", index=False)
    S = r["summaries"]; m = bundle.manifest
    kind = {"assembly": "gut assembly analyses", "reads": "gut read-based profiles"}.get(m.get("profile_type", ""), "gut samples"); sp = m.get("source_pipeline", {})
    L = [f"# refmb report — {m['bundle_id']}\n", f"Reference: {m['n_samples']:,} healthy adult {kind} ({sp.get('name', '')} {sp.get('version', '')}), {m['n_studies']} studies. Calibration: {calibration_id}. "
         f"Scored {S['analysis_id'].nunique() if not S.empty else 0} sample(s); {len(r['rejections'])} not scored (rejections.tsv).\n",
         "Read this first: a healthy adult gut sample measured like the reference has about 5% of features outside the 2.5–97.5 percentile band. "
         "The sample-level excess test asks whether a sample has more than that (q ≤ 0.05 after BH across samples within a layer). "
         "Individual feature calls (low/high) are descriptive; binary calls are less reproducible across re-assemblies than percentiles.\n"]
    if not S.empty:
        piv = S.pivot_table(index="analysis_id", columns="layer", values="frac_outside_raw")
        L += ["## Fraction of assessed features outside the reference band\n", piv.round(3).to_markdown(), ""]
        sig = S[S["excess_outside_q"] <= 0.05]
        L += ["## Samples with a significant excess (q ≤ 0.05)\n", (sig[["analysis_id", "layer", "n_assessed", "n_low_raw", "n_high_raw", "excess_outside_ratio", "excess_outside_q"]].round(3).to_markdown(index=False) if len(sig) else "none"), ""]
        if "core_distance_pct" in S.columns:
            ld = S.drop_duplicates("analysis_id")[[c for c in ["analysis_id", "quality_band", "genome_equivalents_cov", "core_distance_pct", "neighbor_studies"] if c in S.columns]]
            L += ["## Landscape (family CLR PCA of the reference)\n", "core_distance_pct = percentile of the sample's distance from the reference centroid among reference samples.\n", ld.round(2).to_markdown(index=False), ""]
        top = r["scores"][r["scores"]["call"].isin(["low", "high"])].copy()
        if len(top):
            top["dev"] = (top["percentile"] - 50).abs()
            L += ["## Most extreme features per sample (top 10 by percentile distance)\n"]
            for sid, g in top.sort_values("dev", ascending=False).groupby("analysis_id"):
                L += [f"### {sid}\n", g.head(10)[["layer", "feature_id", "band", "value", "percentile", "call"]].round(3).to_markdown(index=False), ""]
        ebm = r["scores"][r["scores"]["call"] == "expected_but_missing"].groupby(["analysis_id", "layer"]).size().unstack(fill_value=0)
        L += ["## Expected-but-missing features (reference prevalence ≥ 0.70 in the sample's band, not detected)\n", ebm.to_markdown(), ""]
    open(os.path.join(out, "report.md"), "w").write("\n".join(L))


def cmd_score(a, norm=None):
    bundle = Bundle(a.bundle); norm = norm if norm is not None else load_normalized(a.normalized)
    cal = None
    if getattr(a, "calibration", None):
        from refmb.calibrate import Calibration
        cal = Calibration.load(a.calibration)
    r = score(norm, bundle, calibration=cal)
    write_report(r, bundle, a.out, "none" if cal is None else cal.id)
    print(json.dumps({"scored": int(r["summaries"]["analysis_id"].nunique()) if not r["summaries"].empty else 0, "rejected": len(r["rejections"]), "rows": len(r["scores"]), "out": a.out}))


def cmd_bundles(a):
    found = list_bundles()
    if not bundle_dirs():
        print(f"{ENV_VAR} is not set; give --bundle a directory path or export {ENV_VAR}=/dir/with/bundles", file=sys.stderr)
    for b in found:
        m = json.load(open(os.path.join(b, "manifest.json")))
        print(f"{m['bundle_id']}\t{m.get('profile_type', '')}\t{m['n_samples']} samples, {m['n_studies']} studies\t{b}")
    if bundle_dirs() and not found:
        print(f"no bundles under {ENV_VAR} = {os.pathsep.join(bundle_dirs())}", file=sys.stderr)


def main(argv=None):
    ap = argparse.ArgumentParser(prog="refmb", description="Score gut metagenomes against a healthy-adult reference bundle measured the same way.")
    ap.add_argument("--version", action="version", version=f"refmb {__version__}")
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("bundles", help=f"list bundles found under ${ENV_VAR}")
    p = sub.add_parser("normalize"); p.add_argument("--bundle", required=True); p.add_argument("--query", nargs="+", required=True); p.add_argument("--out", required=True)
    p = sub.add_parser("import-mgnify"); p.add_argument("--bundle", required=True); p.add_argument("--analysis-dir", nargs="+", required=True); p.add_argument("--out", required=True)
    p = sub.add_parser("score"); p.add_argument("--bundle", required=True); p.add_argument("--normalized", required=True); p.add_argument("--calibration"); p.add_argument("--out", required=True)
    p = sub.add_parser("import-metaphlan"); p.add_argument("--bundle", required=True); p.add_argument("--input", nargs="+", required=True); p.add_argument("--reads-tsv"); p.add_argument("--humann-genefamilies", nargs="+"); p.add_argument("--humann-pathabundance", nargs="+"); p.add_argument("--out", required=True)
    p = sub.add_parser("run-metaphlan"); p.add_argument("--bundle", required=True); p.add_argument("--input", nargs="+", required=True); p.add_argument("--reads-tsv"); p.add_argument("--humann-genefamilies", nargs="+"); p.add_argument("--humann-pathabundance", nargs="+"); p.add_argument("--out", required=True); p.add_argument("--calibration")
    p = sub.add_parser("run-mgnify"); p.add_argument("--bundle", required=True); p.add_argument("--analysis-dir", nargs="+", required=True); p.add_argument("--out", required=True); p.add_argument("--calibration")
    a = ap.parse_args(argv)
    if a.cmd == "bundles":
        return cmd_bundles(a)
    a.bundle = resolve_bundle(a.bundle)
    if a.cmd == "normalize":
        cmd_normalize(a)
    elif a.cmd == "import-mgnify":
        cmd_import_mgnify(a)
    elif a.cmd == "score":
        cmd_score(a)
    elif a.cmd == "import-metaphlan":
        cmd_import_metaphlan(a)
    elif a.cmd == "run-metaphlan":
        a.normalized = a.out; norm = cmd_import_metaphlan(a); cmd_score(a, norm)
    elif a.cmd == "run-mgnify":
        a.normalized = a.out; norm = cmd_import_mgnify(a); cmd_score(a, norm)


if __name__ == "__main__":
    main()
