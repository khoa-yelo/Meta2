"""checkgm.normalize — pipeline-agnostic normalization of ONE sample into the reference measurement space.

Input is the query format (checkgm/query_format.md): contigs (id, length, coverage), CDS (id, contig, KO list, Pfam list,
optional KOfam KO), CDS taxonomy (cds_id, NCBI taxid), optional KEGG module completeness. Output is the same set of
normalized tables the reference pool was built from, computed with the rules shipped in the bundle
(bundle/normalization/: rules.json, marker_panel.tsv, clr_basis_<rank>.txt) and the NCBI backbone.

Semantics (identical to the baseline build, docs/configs/s5.yaml and the bundle's rules.json):
  * taxonomy: per-CDS taxid -> lineage; per contig majority vote (>= 50% of resolvable CDS calls) at each rank; contig weight
    = length x coverage; explicit `unmapped` bucket over CDS-bearing contigs; proportions among mapped weight; CLR over the
    bundle's basis with multiplicative replacement (delta = 0.5 x smallest non-zero proportion in the sample).
  * gene layers: per feature cds_count and cov_sum (coverage of the contig carrying the CDS); genome equivalents = median
    cov_sum over the validated marker panel (absent marker = 0); copies_per_genome = cov_sum / GE.
  * quality: total assembled length over ALL contigs -> band; floors from rules.json; N50 over CDS-bearing contigs.
"""
from __future__ import annotations

import csv
import gzip
import json
import os
import re
from collections import Counter, defaultdict

import numpy as np
import pandas as pd

RANKS = ["phylum", "class", "order", "family", "genus"]
MAJORITY = 0.5
CONTIG_RE = re.compile(r"length-(\d+)-cov-([\d.]+)")


class Backbone:
    """The NCBI taxonomy shipped in a bundle (normalization/backbone/): lineage taxid per rank, merged-id table and names.

    `resolve(taxid)` follows merged ids and returns the current taxid, or None when the id is unknown to the backbone;
    `lineage[rank][taxid]` gives the ancestor at that rank; `name[taxid]` the scientific name."""

    def __init__(self, path: str):
        tax = pd.read_parquet(os.path.join(path, "ncbi_taxonomy.parquet"))
        merged = pd.read_parquet(os.path.join(path, "ncbi_merged.parquet"))
        self.lineage = {r: dict(zip(tax["taxid"].to_numpy(), tax[f"{r}_taxid"].to_numpy())) for r in RANKS}
        self.merged = dict(zip(merged["old"], merged["new"]))
        self.valid = set(tax["taxid"].tolist())
        self.name = dict(zip(tax["taxid"], tax["name"]))

    def resolve(self, tid: int):
        """Current taxid for `tid` after following merged ids, or None when the backbone does not know it."""
        tid = self.merged.get(tid, tid)
        return tid if tid in self.valid else None


class BundleSpec:
    """Normalization rules shipped inside a bundle."""

    def __init__(self, bundle: str):
        nd = os.path.join(bundle, "normalization")
        self.rules = json.load(open(os.path.join(nd, "rules.json")))
        self.markers = sorted(pd.read_csv(os.path.join(nd, "marker_panel.tsv"), sep="\t")["ko"])
        self.basis = {r: [l.strip() for l in open(os.path.join(nd, f"clr_basis_{r}.txt")) if l.strip()] for r in self.rules["clr_ranks"]}
        bb = os.path.join(nd, "backbone")
        self.backbone_path = bb if os.path.isdir(bb) else self.rules.get("backbone_fallback")

    def band(self, total_length: float) -> str:
        """Assembly-size class ('low', 'medium', 'high') of a total assembled length in bp from rules.json; 'unbanded' outside all classes."""
        for name, (lo, hi) in self.rules["quality_bands"].items():
            if total_length >= lo and (hi is None or total_length < hi):
                return name
        return "unbanded"


def split_list(s) -> list[str]:
    """Split a `;`- or `,`-separated annotation list (KO ids, Pfam accessions) into a list; empty or NaN -> []."""
    if not isinstance(s, str) or not s:
        return []
    return [x.replace("ko:", "") for x in re.split(r"[;,]", s) if x]


def _n50(lengths: np.ndarray) -> float:
    if not lengths.size:
        return float("nan")
    srt = np.sort(lengths)[::-1]; cs = np.cumsum(srt)
    return float(srt[np.searchsorted(cs, cs[-1] / 2)])


def normalize(sample_id: str, contigs: pd.DataFrame, cds: pd.DataFrame, cds_tax: pd.DataFrame, modules: pd.DataFrame | None,
              spec: BundleSpec, bb: Backbone, total_length_override: float | None = None, kofam: pd.DataFrame | None = None) -> dict:
    """cds_tax: cds_id, taxid and optionally contig_id (used directly when present, otherwise looked up through cds).
    kofam: optional separate table cds_id, contig_id, ko_kofam (one best KOfam hit per CDS); if absent, cds['ko_kofam'] is used."""
    rules = spec.rules
    contigs = contigs.drop_duplicates("contig_id").set_index("contig_id")
    clen = contigs["length"].astype(float); ccov = contigs["coverage"].astype(float)
    n_cds_nocov = int(cds["contig_id"].map(ccov).isna().sum())
    ccov = ccov.fillna(0.0)
    c_of = dict(zip(cds["cds_id"], cds["contig_id"]))
    tax_contig = cds_tax["contig_id"].to_numpy() if "contig_id" in cds_tax.columns else cds_tax["cds_id"].map(c_of).to_numpy()
    cds_contigs = sorted(set(cds["contig_id"]) | set(pd.Series(tax_contig).dropna()))
    cds_contigs = [c for c in cds_contigs if c in clen.index]
    L = clen.reindex(cds_contigs); C = ccov.reindex(cds_contigs); W = L * C
    W_all = float(W.sum())

    # ---- taxonomy: per-contig majority vote per rank -------------------------------------------------------------------
    calls: dict[str, dict[str, list]] = defaultdict(lambda: {r: [] for r in RANKS})
    n_hits = 0; n_bad = 0
    for c, tid in zip(tax_contig, cds_tax["taxid"].to_numpy()):
        try:
            t = bb.resolve(int(tid))
        except (TypeError, ValueError):
            t = None
        if t is None:
            n_bad += 1; continue
        n_hits += 1
        if c is None or c not in clen.index:
            continue
        cc = calls[c]
        for r in RANKS:
            v = bb.lineage[r].get(t, -1)
            if v != -1:
                cc[r].append(v)
    tax_out = {}; mapped_fraction = {}
    for r in RANKS:
        w = Counter(); w_unm = 0.0
        for c in cds_contigs:
            wt = float(W[c])
            cl = calls[c][r] if c in calls else []
            if cl:
                top, k = Counter(cl).most_common(1)[0]
                if k / len(cl) >= MAJORITY:
                    w[top] += wt; continue
            w_unm += wt
        mapped_fraction[r] = (W_all - w_unm) / W_all if W_all else float("nan")
        tax_out[r] = (w, w_unm)

    def clr_table(rank: str) -> pd.DataFrame:
        """Proportions among mapped weight over the bundle's basis at `rank`, and their CLR with multiplicative replacement of zeros."""
        w, _ = tax_out[rank]; basis = spec.basis[rank]
        tot = sum(w.values())
        prop = pd.Series({str(k): v / tot for k, v in w.items()}) if tot > 0 else pd.Series(dtype=float)
        b = prop.reindex(basis).fillna(0.0)
        z = int((b == 0).sum())
        delta = rules["clr_delta_fraction_of_min"] * (b[b > 0].min() if (b > 0).any() else 1.0)
        rep = np.where(b > 0, b * (1 - z * delta), delta)
        clr = np.log(rep) - np.log(rep).mean()
        return pd.DataFrame({"feature_id": b.index, "proportion_mapped": b.to_numpy(), "clr": clr, "rank": rank})

    out = {f"taxonomy_{r}": clr_table(r) for r in spec.basis}

    # ---- gene layers --------------------------------------------------------------------------------------------------
    def gene_layer(tab: pd.DataFrame | None, col: str) -> pd.DataFrame:
        """Per feature of annotation column `col`: number of CDS and sum of the carrying contigs' coverage."""
        if tab is None or col not in tab.columns:
            return pd.DataFrame(columns=["feature_id", "cds_count", "cov_sum"])
        cov_of_cds = tab["contig_id"].map(ccov).fillna(0.0).to_numpy()
        cnt = Counter(); cs = Counter()
        for lst, cv in zip(tab[col].to_numpy(), cov_of_cds):
            for k in sorted(set(split_list(lst))):
                cnt[k] += 1; cs[k] += cv
        return pd.DataFrame({"feature_id": sorted(cnt), "cds_count": [cnt[k] for k in sorted(cnt)], "cov_sum": [cs[k] for k in sorted(cnt)]})

    ko = gene_layer(cds, "ko"); pf = gene_layer(cds, "pfam")
    kf = gene_layer(kofam if kofam is not None else cds, "ko_kofam")
    mk = ko.set_index("feature_id")["cov_sum"].reindex(spec.markers).fillna(0.0)
    ge = float(np.median(mk.to_numpy())); markers_detected = int((mk > 0).sum())
    for name, t in (("ko_eggnog", ko), ("ko_kofam", kf), ("pfam", pf)):
        t = t.copy(); t["copies_per_genome"] = t["cov_sum"] / ge if ge > 0 else np.nan
        out[name] = t
    if modules is not None and len(modules):
        m = modules.rename(columns={"module_id": "feature_id"})[["feature_id", "completeness"]].copy()
        m["completeness"] = m["completeness"].astype(float).clip(0, 1)
        out["module"] = m

    # ---- quality and status -------------------------------------------------------------------------------------------
    total_len = float(total_length_override) if total_length_override is not None and not np.isnan(total_length_override) else float(clen.sum())
    n50 = _n50(L.to_numpy())
    n_cds = int(len(cds))
    prim = rules["taxonomy_primary_rank"]
    status = "retained"; rej = []
    fl = rules["floors"]
    if total_len < fl["min_total_length_bp"]: rej.append(("BELOW_MIN_ASSEMBLY_QUALITY", f"total length < {fl['min_total_length_bp']}"))
    if not (n50 >= fl["min_n50_bp"]): rej.append(("BELOW_MIN_ASSEMBLY_QUALITY", f"N50 < {fl['min_n50_bp']}"))
    if n_cds < fl["min_predicted_cds"]: rej.append(("BELOW_MIN_ASSEMBLY_QUALITY", f"annotated CDS < {fl['min_predicted_cds']}"))
    if markers_detected < fl["min_markers_detected"]: rej.append(("NO_MARKER_PANEL_SIGNAL", f"< {fl['min_markers_detected']} kept marker KOs detected"))
    if not (ge >= fl["min_genome_equivalents"]): rej.append(("NO_MARKER_PANEL_SIGNAL", f"genome equivalents < {fl['min_genome_equivalents']}"))
    if not (mapped_fraction[prim] >= fl["min_mapped_fraction"]): rej.append(("LOW_MAPPED_FRACTION", f"{prim} mapped fraction < {fl['min_mapped_fraction']}"))
    if rej:
        status = rej[0][0]
    qc = {"sample_id": sample_id, "status": status, "rejections": rej, "quality_band": spec.band(total_len) if status == "retained" else "unbanded",
          "total_length_all_contigs": total_len, "total_length_source": "override" if total_length_override is not None else "contigs_table",
          "total_length_cds_contigs": float(L.sum()), "n_contigs_with_cds": len(cds_contigs), "n50_cds_contigs": n50, "n_cds": n_cds,
          "n_cds_without_coverage": n_cds_nocov, "n_taxonomy_hits": n_hits, "n_taxonomy_bad_taxid": n_bad,
          "genome_equivalents_cov": ge, "marker_kos_detected_kept": markers_detected, "n_markers_kept": len(spec.markers),
          **{f"mapped_fraction_{r}": mapped_fraction[r] for r in RANKS}}
    out["qc"] = qc
    return out


# ---- adapters -------------------------------------------------------------------------------------------------------------
def contig_of(cds_id: str) -> str:
    """Contig name of a CDS id of the form <contig>_<start>_<end>_<strand> or <contig>_<n>; the id itself when it has neither suffix."""
    m = re.match(r"^(.*?)(?:_\d+_\d+_[+-]|_\d+)$", cds_id)
    return m.group(1) if m else cds_id


def _find(d: str, suffix: str):
    for f in os.listdir(d):
        if f.endswith(suffix):
            return os.path.join(d, f)
    return None


def from_mgnify_v5(analysis_dir: str) -> dict:
    """MGnify v5 assembly analysis folder -> query-format tables (the reference's own measurement)."""
    gff = _find(analysis_dir, "_annotations.gff.bgz"); dmd = _find(analysis_dir, "_diamond.tsv.gz")
    if not (gff and dmd):
        raise FileNotFoundError(f"{analysis_dir}: missing gff/diamond")
    contigs: dict[str, tuple] = {}
    cds_rows = []
    with gzip.open(gff, "rt") as f:
        for line in f:
            if line[0] == "#":
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[2] != "CDS":
                continue
            c = p[0]
            if c not in contigs:
                m = CONTIG_RE.search(c); contigs[c] = (int(m.group(1)), float(m.group(2))) if m else (None, None)
            attrs = dict(kv.split("=", 1) for kv in p[8].split(";") if "=" in kv)
            cds_rows.append((f"{c}_{p[3]}_{p[4]}_{p[6]}", c, attrs.get("kegg", ""), attrs.get("pfam", "")))   # GFF ID attribute is the contig name
    cds = pd.DataFrame(cds_rows, columns=["cds_id", "contig_id", "ko", "pfam"])
    # taxonomy: DIAMOND rows are per CDS; the CDS id may not match the GFF ID attribute, so map through the contig instead
    tax_rows = []
    with gzip.open(dmd, "rt") as f:
        header = f.readline().rstrip("\n").split("\t"); i_q = header.index("contig_name"); i_t = header.index("tax_id")
        for line in f:
            p = line.rstrip("\n").split("\t")
            s = p[i_t].replace("TaxID=", "")
            q = p[i_q]; c = contig_of(q)
            if c not in contigs:
                m = CONTIG_RE.search(c); contigs[c] = (int(m.group(1)), float(m.group(2))) if m else (None, None)
            tax_rows.append((q, c, int(s) if s.isdigit() else -1))
    cds_tax = pd.DataFrame(tax_rows, columns=["cds_id", "contig_id", "taxid"])
    # KOfam per CDS (best hit by full-sequence score), keyed by the DIAMOND-style CDS id -> contig
    hmm = _find(analysis_dir, "_kofam_hmm.tsv.gz")
    if hmm and os.path.getsize(hmm) > 40:
        best: dict[str, tuple[float, str]] = {}
        with gzip.open(hmm, "rt") as f:
            header = f.readline().rstrip("\n").split("\t")
            i_k = header.index("target_name"); i_q = header.index("query_name"); i_s = header.index("full_sequence_score")
            for line in f:
                p = line.rstrip("\n").split("\t")
                if len(p) <= max(i_k, i_q, i_s):
                    continue
                try:
                    sc = float(p[i_s])
                except ValueError:
                    continue
                if p[i_q] not in best or sc > best[p[i_q]][0]:
                    best[p[i_q]] = (sc, p[i_k])
        kf = pd.DataFrame([(q, contig_of(q), k) for q, (sc, k) in best.items()], columns=["cds_id", "contig_id", "ko_kofam"])
        for c in set(kf["contig_id"]) - set(contigs):   # KOfam may hit CDS on contigs absent from the GFF/DIAMOND tables
            m = CONTIG_RE.search(c); contigs[c] = (int(m.group(1)), float(m.group(2))) if m else (None, None)
    else:
        kf = pd.DataFrame(columns=["cds_id", "contig_id", "ko_kofam"])
    modules = None
    mp = _find(analysis_dir, "_kegg_pathways.csv")
    if mp:
        rows = []
        with open(mp, newline="") as f:
            r = csv.reader(f); next(r, None)
            for row in r:
                if len(row) >= 2 and row[0].startswith("M"):
                    rows.append((row[0], float(row[1]) / 100.0))
        modules = pd.DataFrame(rows, columns=["module_id", "completeness"])
    total_len = None
    qcs = _find(analysis_dir, "_qc_summary.out")
    if qcs and os.path.getsize(qcs) > 0:
        for line in open(qcs):
            p = line.rstrip("\n").split("\t")
            if len(p) == 2 and p[0] == "bp_count":
                total_len = float(p[1])
    ctab = pd.DataFrame([(c, L, cv) for c, (L, cv) in contigs.items() if L is not None], columns=["contig_id", "length", "coverage"])
    return {"contigs": ctab, "cds": cds, "cds_tax": cds_tax, "kofam": kf, "modules": modules, "total_length_override": total_len}


def normalize_mgnify_v5(sample_id: str, analysis_dir: str, spec: BundleSpec, bb: Backbone) -> dict:
    """Normalize one MGnify pipeline v5 assembly analysis folder with the bundle's rules (the reference's own measurement).

    Reads the folder with from_mgnify_v5() and calls normalize(); returns {layer: DataFrame, ..., "qc": dict} as normalize() does."""
    t = from_mgnify_v5(analysis_dir)
    return normalize(sample_id, t["contigs"], t["cds"], t["cds_tax"], t["modules"], spec, bb,
                     total_length_override=t["total_length_override"], kofam=t["kofam"])


def read_query_dir(qdir: str) -> dict:
    """Generic query directory (checkgm/query_format.md) -> tables."""
    rd = lambda n, **kw: pd.read_csv(os.path.join(qdir, n), sep="\t", dtype=str, keep_default_na=False, **kw) if os.path.exists(os.path.join(qdir, n)) else None
    contigs = rd("contigs.tsv"); cds = rd("cds.tsv"); tax = rd("cds_taxonomy.tsv"); mod = rd("modules.tsv")
    if contigs is None or cds is None or tax is None:
        raise FileNotFoundError(f"{qdir}: contigs.tsv, cds.tsv and cds_taxonomy.tsv are required")
    contigs["length"] = contigs["length"].astype(float); contigs["coverage"] = pd.to_numeric(contigs["coverage"], errors="coerce")
    tax["taxid"] = pd.to_numeric(tax["taxid"], errors="coerce").fillna(-1).astype(int)
    if mod is not None:
        mod["completeness"] = mod["completeness"].astype(float)
    meta = json.load(open(os.path.join(qdir, "sample.json"))) if os.path.exists(os.path.join(qdir, "sample.json")) else {}
    keep = [c for c in ("cds_id", "contig_id", "taxid") if c in tax.columns]
    return {"contigs": contigs, "cds": cds, "cds_tax": tax[keep], "modules": mod, "meta": meta}
