"""refmb.normalize_reads — read-based query adapter: MetaPhlAn output -> reference-ready normalized tables.

Accepts MetaPhlAn 3 (or 2/4-style) single-sample profiles and merged tables:
  single-sample : header lines starting with '#', including '#<N> reads processed' when present; columns clade_name,
                  [NCBI_tax_id,] relative_abundance[, additional_species]
  merged table  : optional '#mpa_...' line, then 'clade_name[\tNCBI_tax_id]\t<sample columns>'
Rows at species level (contain 's__', no 't__') are used. Taxid resolution: the file's NCBI_tax_id (last element of the lineage) via
the bundle backbone; else the bundle's species_taxid_map (MetaPhlAn 3 names); else NCBI name lookup on the backbone.
Normalization follows the rules the reference pool was built with (bundle/normalization/rules.json): proportions of the resolved profile, CLR over the bundle's basis per rank with
multiplicative replacement (delta = 0.5 x smallest non-zero proportion), depth band from reads processed, mapped-fraction floor.
"""
from __future__ import annotations

import json
import os
import re

import numpy as np
import pandas as pd

RANKS = ["species", "genus", "family"]
FUNCTION_LAYERS = {"ko_humann": ("UNMAPPED", "UNGROUPED"), "pathway_humann": ("UNMAPPED", "UNINTEGRATED")}   # layer -> QC rows that are not features


class ReadsBundleSpec:
    """Normalization rules and lookup tables of a read-based bundle.

    Loads `normalization/rules.json` (depth bands, CLR delta, mapped-fraction floor), the CLR basis lists per rank and
    per function layer, the bundle's own MetaPhlAn species -> taxid map, and the NCBI backbone (taxonomy, merged ids).
    `resolve()` places a species on the backbone the way the reference did; `band()` gives the depth band of a read count."""

    def __init__(self, bundle: str):
        nd = os.path.join(bundle, "normalization")
        if not os.path.isfile(os.path.join(nd, "rules.json")):
            raise FileNotFoundError(f"bundle has no normalization/rules.json: {bundle}")
        self.rules = json.load(open(os.path.join(nd, "rules.json")))
        self.basis = {r: [l.strip() for l in open(os.path.join(nd, f"clr_basis_{r}.txt")) if l.strip()] for r in RANKS}
        self.fbasis = {l: [x.strip() for x in open(os.path.join(nd, f"clr_basis_{l}.txt")) if x.strip()] for l in FUNCTION_LAYERS if os.path.exists(os.path.join(nd, f"clr_basis_{l}.txt"))}
        self.ko_map_path = os.path.join(nd, "map_ko_uniref90.txt.gz"); self._ko_map = None
        m = pd.read_csv(os.path.join(nd, "species_taxid_map.tsv"), sep="\t", dtype={"species_name": str}).set_index("species_name")
        self.name_map = {n: {r: int(m.loc[n, r]) for r in RANKS} for n in m.index}
        tax = pd.read_parquet(os.path.join(nd, "backbone", "ncbi_taxonomy.parquet"))
        merged = pd.read_parquet(os.path.join(nd, "backbone", "ncbi_merged.parquet"))
        self.merged = dict(zip(merged["old"], merged["new"])); self.valid = set(tax["taxid"].tolist())
        self.rank_of = dict(zip(tax["taxid"], tax["rank"]))
        self.lin = {r: dict(zip(tax["taxid"], tax[f"{r}_taxid"])) for r in ["genus", "family"] if f"{r}_taxid" in tax.columns}
        self.name2 = {}
        for r in tax[tax["rank"].isin(["species", "genus", "family"])].itertuples(index=False):
            self.name2.setdefault(r.name, r.taxid)

    def ko_map(self) -> dict:
        """UniRef90 family -> list of KOs (HUMAnN map_ko_uniref90, shipped in the bundle)."""
        if self._ko_map is None:
            import gzip
            m = {}
            with gzip.open(self.ko_map_path, "rt") as fh:
                for line in fh:
                    f = line.rstrip("\n").split("\t")
                    for u in f[1:]:
                        m.setdefault(u, []).append(f[0])
            self._ko_map = m
        return self._ko_map

    def band(self, reads):
        """Depth band ('low', 'medium', 'high') of a read count from the bundle's depth_bands_reads; 'unknown' when the count is missing."""
        if reads is None or (isinstance(reads, float) and np.isnan(reads)):
            return "unknown"
        for k, (lo, hi) in self.rules["depth_bands_reads"].items():
            if reads >= lo and (hi is None or reads < hi):
                return k
        return "unbanded"

    def resolve(self, species_name: str, taxid=None) -> dict:
        """Place one MetaPhlAn species on the bundle's taxonomy: {species, genus, family} taxids (-1 where it cannot be placed),
        by the bundle's own species map first and the NCBI backbone (with merged ids) second."""
        out = {"species": -1, "genus": -1, "family": -1}
        if species_name in self.name_map:      # the bundle's own map first: a query species is placed exactly where the reference placed it
            return dict(self.name_map[species_name])
        if taxid is not None:
            t = self.merged.get(int(taxid), int(taxid))
            if t in self.valid:
                rk = self.rank_of.get(t)
                if rk == "species":
                    out = {"species": t, "genus": int(self.lin["genus"].get(t, -1)), "family": int(self.lin["family"].get(t, -1))}
                elif rk == "genus":
                    out = {"species": -1, "genus": t, "family": int(self.lin["family"].get(t, -1))}
                elif rk == "family":
                    out = {"species": -1, "genus": -1, "family": t}
                if out["family"] != -1:
                    return out
        if species_name in self.name_map:
            return dict(self.name_map[species_name])
        n = species_name[3:].replace("_", " "); tid = self.name2.get(n)
        if tid is not None and self.rank_of.get(tid) == "species":
            return {"species": tid, "genus": int(self.lin["genus"].get(tid, -1)), "family": int(self.lin["family"].get(tid, -1))}
        g = self.name2.get(n.split(" ")[0])
        if g is not None and self.rank_of.get(g) == "genus":
            return {"species": -1, "genus": g, "family": int(self.lin["family"].get(g, -1))}
        return out


def _species_rows(clades):
    return [c for c in clades if "s__" in c and "t__" not in c]


EMPTY_TABLE_COLUMNS = ["species_name", "taxid", "rel_abundance"]
_CLADE_HEADERS = ("clade_name", "sampleid", "id")   # first header cell of MetaPhlAn 3 / 2 single-sample files and merged tables


def read_metaphlan(path: str) -> dict:
    """Parse a MetaPhlAn profile -> {sample_id: {"table": DataFrame(species_name, taxid, rel_abundance), "reads": int | None}}.

    Accepts a single-sample file (comment header with `#clade_name ... relative_abundance`, optional `#N reads processed`)
    or a merged table (`clade_name [NCBI_tax_id] sample1 sample2 ...`). Only species-level rows are kept; a profile without
    any gets an empty table, which normalize_profile() turns into the EMPTY_PROFILE rejection. A file that does not have
    MetaPhlAn's columns raises ValueError with a plain explanation."""
    reads = None; header = None
    with open(path, errors="replace") as f:
        lines = f.read().splitlines()
    data = []
    for ln in lines:
        if not ln.strip():
            continue
        if ln.startswith("#"):
            m = re.match(r"#\s*([\d,]+)\s+reads processed", ln)
            if m:
                reads = int(m.group(1).replace(",", ""))
            if ln.lower().startswith("#clade_name") or ln.lower().startswith("#sampleid"):
                header = ln[1:].split("\t")
            continue
        data.append(ln.split("\t"))
    if header is None or header[0].strip().lower() not in _CLADE_HEADERS:
        if not data:
            raise ValueError(f"{path}: not a MetaPhlAn profile (no '#clade_name' header line and no table rows)")
        # merged table: first data row is the header
        header = data[0]; data = data[1:]
    cols = [h.strip() for h in header]
    if cols[0].lstrip("#").lower() not in _CLADE_HEADERS:
        raise ValueError(f"{path}: not a MetaPhlAn profile; expected a first column 'clade_name' (got {cols[:3]}). "
                         "Give a MetaPhlAn 3 single-sample output or a merged table")
    clade_col = cols[0]
    has_tid = "NCBI_tax_id" in cols
    single = "relative_abundance" in cols
    sample_cols = [c for c in cols if c not in (clade_col, "NCBI_tax_id", "additional_species", "relative_abundance")]
    if not single and not sample_cols:
        raise ValueError(f"{path}: not a MetaPhlAn profile; no 'relative_abundance' column and no sample columns (columns: {cols})")
    width = len(cols)
    data = [row + [""] * (width - len(row)) if len(row) < width else row[:width] for row in data]
    df = pd.DataFrame(data, columns=cols) if data else pd.DataFrame(columns=cols)
    df = df[df[clade_col].isin(_species_rows(df[clade_col]))].copy()
    df["species_name"] = "s__" + df[clade_col].str.split("s__").str[-1]
    df["taxid"] = df["NCBI_tax_id"].str.split("|").str[-1].where(df["NCBI_tax_id"].str.strip().ne(""), None) if has_tid else None
    out = {}
    if single:                                            # single-sample profile
        sid = os.path.splitext(os.path.basename(path))[0]
        t = df[["species_name", "taxid"]].copy(); t["rel_abundance"] = pd.to_numeric(df["relative_abundance"], errors="coerce").fillna(0.0) / 100.0
        out[sid] = {"table": t[t.rel_abundance > 0], "reads": reads}
    else:                                                  # merged table
        for sc in sample_cols:
            t = df[["species_name", "taxid"]].copy(); t["rel_abundance"] = pd.to_numeric(df[sc], errors="coerce").fillna(0.0) / 100.0
            out[sc] = {"table": t[t.rel_abundance > 0], "reads": None}
    return out


def normalize_profile(sample_id: str, table: pd.DataFrame, reads, spec: ReadsBundleSpec) -> dict:
    """Normalize one parsed MetaPhlAn profile into the reference measurement space.

    table: species_name (MetaPhlAn `s__` name), taxid (string or None), rel_abundance (fraction, not percent) as returned by
    read_metaphlan(); reads: reads processed (int) or None. Each species is placed on the bundle's taxonomy; per rank
    (species, genus, family) the proportions over the bundle's CLR basis are centred-log-ratio transformed with multiplicative
    replacement of zeros. Returns {"taxonomy_<rank>": DataFrame(feature_id, proportion_mapped, clr, rank), "qc": dict};
    qc["status"] is "retained", or a rejection code (EMPTY_PROFILE when there are no species rows, LOW_MAPPED_FRACTION when
    less than the bundle's floor of the abundance could be placed at family level)."""
    tot = float(table["rel_abundance"].sum()) if len(table) else 0.0
    if len(table):
        res = pd.DataFrame([spec.resolve(n, (int(t) if t not in (None, "", "nan") and str(t).lstrip("-").isdigit() else None)) for n, t in zip(table["species_name"], table["taxid"])], index=table.index)
    else:
        res = pd.DataFrame({r: pd.Series(dtype=int) for r in RANKS})
    out = {}; mapped = {}
    for rank in RANKS:
        m = res[rank] != -1
        agg = table.loc[m, "rel_abundance"].groupby(res.loc[m, rank].astype(int).astype(str)).sum()
        mapped[rank] = float(agg.sum()) / tot if tot > 0 else 0.0
        prop = agg / agg.sum() if agg.sum() > 0 else agg
        basis = spec.basis[rank]; b = prop.reindex(basis).fillna(0.0)
        z = int((b == 0).sum()); delta = spec.rules["clr_delta_fraction_of_min"] * (b[b > 0].min() if (b > 0).any() else 1.0)
        with np.errstate(invalid="ignore", divide="ignore"):   # a profile with nothing placed gives NaN CLRs and is rejected below, without a warning
            rep = np.where(b > 0, b * (1 - z * delta), delta); clr = np.log(rep) - np.log(rep).mean()
        out[f"taxonomy_{rank}"] = pd.DataFrame({"feature_id": basis, "proportion_mapped": b.to_numpy(), "clr": clr, "rank": rank})
    rej = []
    if tot <= 0 or len(table) == 0:
        rej.append(("EMPTY_PROFILE", "no species-level abundance"))
    if mapped["family"] < spec.rules["min_mapped_fraction_family"]:
        rej.append(("LOW_MAPPED_FRACTION", f"family mapped fraction < {spec.rules['min_mapped_fraction_family']}"))
    band = spec.band(reads)
    out["qc"] = {"sample_id": sample_id, "status": rej[0][0] if rej else "retained", "rejections": rej, "quality_band": band if not rej else "unbanded",
                 "reads_processed": reads, "n_species_rows": int(len(table)), "mapped_fraction_species": mapped["species"], "mapped_fraction_genus": mapped["genus"],
                 "mapped_fraction_family": mapped["family"], "genome_equivalents_cov": np.nan, "total_length_all_contigs": np.nan, "n50_cds_contigs": np.nan}
    return out


def read_humann(path: str) -> dict:
    """HUMAnN 3 gene-family or pathway table (single sample or joined) -> {sample_id: Series feature -> abundance}. Stratified rows
    (feature|species) are dropped. Units do not matter (RPK, CPM or relative abundance): the layer is renormalized among assigned features."""
    df = pd.read_csv(path, sep="\t", dtype=str); fcol = df.columns[0]
    df = df[~df[fcol].str.contains("|", regex=False)]
    out = {}
    for c in df.columns[1:]:
        sid = re.sub(r"(_Abundance(-RPKs|-CPM|-RELAB)?|_genefamilies|_pathabundance)+$", "", c)
        v = pd.to_numeric(df[c], errors="coerce").fillna(0.0); out[sid] = pd.Series(v.to_numpy(), index=df[fcol].to_numpy())
    return out


def normalize_function(values: pd.Series, layer: str, spec: ReadsBundleSpec) -> tuple:
    """-> (DataFrame feature_id, proportion_mapped, clr over the bundle basis; qc dict). Gene families given as UniRef90 are regrouped
    to KO exactly as the reference was (sum per KO, a family with several KOs counts for each, the rest is UNGROUPED)."""
    v = values[values > 0]; special = FUNCTION_LAYERS[layer]
    if layer == "pathway_humann":
        v = v.groupby(v.index.str.split(":").str[0].str.strip()).sum()
    elif v.index.str.startswith("UniRef").any():
        km = spec.ko_map(); acc = {}
        for f, x in v.items():
            if f == "UNMAPPED":
                acc["UNMAPPED"] = acc.get("UNMAPPED", 0.0) + x; continue
            for k in km.get(f, ["UNGROUPED"]):
                acc[k] = acc.get(k, 0.0) + x
        v = pd.Series(acc, dtype=float)
    else:
        v = v.groupby(v.index.str.split(":").str[0].str.strip()).sum()
    tot = float(v.sum()); qc = {f"{layer}_share_{c.lower()}": (float(v.get(c, 0.0)) / tot if tot > 0 else np.nan) for c in special}
    v = v.drop(labels=[c for c in special if c in v.index]); prop = v / v.sum() if v.sum() > 0 else v
    basis = spec.fbasis[layer]; b = prop.reindex(basis).fillna(0.0); qc[f"{layer}_share_in_basis"] = float(b.sum())
    z = int((b == 0).sum()); delta = spec.rules["clr_delta_fraction_of_min"] * (b[b > 0].min() if (b > 0).any() else 1.0)
    rep = np.where(b > 0, b * (1 - z * delta), delta); clr = np.log(rep) - np.log(rep).mean()
    return pd.DataFrame({"feature_id": basis, "proportion_mapped": b.to_numpy(), "clr": clr, "rank": layer}), qc
