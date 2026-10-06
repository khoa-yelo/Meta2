#!/usr/bin/env python3
"""Helpers for run_tierA.sh: contig filter, Prodigal+FragGeneScan merge (pipeline-v5 CGC semantics), DIAMOND taxonomy join,
KOfam best hits and per-contig KO union, and the final query-directory writer (checkgm/query_format.md)."""
import argparse, gzip, os, re
from collections import defaultdict
import pandas as pd

SPADES_RE = re.compile(r"NODE_(\d+)_length_(\d+)_cov_([\d.]+)")


def read_fasta(path):
    """Yield (header without '>', sequence) pairs from a FASTA file (plain or gzip)."""
    name = None; seq = []
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name = line[1:].strip(); seq = []
            else:
                seq.append(line.strip())
    if name is not None:
        yield name, "".join(seq)


def filter_contigs(a):
    """Copy the contigs of at least --min-length bp to --out (header trimmed to its first word, 80-column lines) and print the counts."""
    n = k = 0
    with open(a.out, "w") as o:
        for name, seq in read_fasta(a.inp):
            n += 1
            if len(seq) >= a.min_length:
                k += 1; o.write(f">{name.split()[0]}\n"); [o.write(seq[i:i + 80] + "\n") for i in range(0, len(seq), 80)]
    print(f"contigs {n} -> {k} (>= {a.min_length} bp)")


def rename_megahit(a):
    """MEGAHIT headers '>k141_12 flag=1 multi=4.0000 len=2164' -> 'NODE_<i>_length_<len>_cov_<multi>' (MGnify's naming; multi = k-mer multiplicity)."""
    n = 0
    with open(a.out, "w") as o:
        for name, seq in read_fasta(a.inp):
            n += 1; m = re.search(r"multi=([\d.]+)", name); L = len(seq)
            o.write(f">NODE_{n}_length_{L}_cov_{m.group(1) if m else '0'}\n"); [o.write(seq[i:i + 80] + "\n") for i in range(0, L, 80)]
    print(f"renamed {n} MEGAHIT contigs")


def merge_cds(a):
    """Prodigal ids <contig>_<n> (header '>contig_n # start # end # strand ...'); FGS ids <contig>_<start>_<end>_<strand>.
    Prodigal has priority; an FGS CDS is kept only if it does not overlap any Prodigal CDS on that contig (any strand)."""
    prod = defaultdict(list); out = []
    for name, seq in read_fasta(a.prodigal_faa):
        parts = name.split(" # "); cid = parts[0]; contig = cid.rsplit("_", 1)[0]
        s, e = int(parts[1]), int(parts[2]); prod[contig].append((s, e)); out.append((cid, seq))
    n_fgs = 0
    if a.fgs_faa and os.path.exists(a.fgs_faa):
        for name, seq in read_fasta(a.fgs_faa):
            cid = name.split()[0]; m = re.match(r"^(.*)_(\d+)_(\d+)_([+-])$", cid)
            if not m:
                continue
            contig, s, e = m.group(1), int(m.group(2)), int(m.group(3))
            if any(not (e < ps or s > pe) for ps, pe in prod.get(contig, [])):
                continue
            out.append((cid, seq)); n_fgs += 1
    with open(a.out, "w") as o:
        for cid, seq in out:
            o.write(f">{cid}\n{seq.rstrip('*')}\n")
    print(f"CDS: prodigal {len(out) - n_fgs}, fgs added {n_fgs}")


def join_taxonomy(a):
    """diamond_raw: sseqid qseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore (best hit per query).
    db_uniref90_result.txt: uniref90_ID <tab> protein_name <tab> num_in_cluster <tab> taxonomy <tab> tax_id <tab> rep_id"""
    d = pd.read_csv(a.diamond, sep="\t", header=None, names=["uniref90_ID", "contig_name", "percentage_of_identical_matches", "lenght", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"], dtype={"uniref90_ID": str})
    d = d.sort_values(["contig_name", "bitscore"], ascending=[True, False]).drop_duplicates("contig_name")
    need = set(d.uniref90_ID)
    rows = {}
    with gzip.open(a.uniref_tax, "rt") as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if p[0] in need:
                rows[p[0]] = p[1:6] if len(p) >= 6 else p[1:] + [""] * (5 - len(p[1:]))
    T = pd.DataFrame.from_dict(rows, orient="index", columns=["protein_name", "num_in_cluster", "taxonomy", "tax_id", "rep_id"]).rename_axis("uniref90_ID").reset_index()
    d = d.merge(T, on="uniref90_ID", how="left"); d["tax_id"] = d["tax_id"].fillna("").str.replace("TaxID=", "", regex=False)
    d["num_in_cluster"] = d["num_in_cluster"].fillna("").str.replace("n=", "", regex=False); d["taxonomy"] = d["taxonomy"].fillna("").str.replace("Tax=", "", regex=False)
    d.to_csv(a.out, sep="\t", index=False); print(f"diamond hits {len(d)}, with taxonomy {d.tax_id.ne('').sum()}")


def kofam_union(a):
    """hmmsearch --domtblout: target = the SEQUENCE (CDS, col 0), query = the PROFILE (KO, col 3); full-sequence score col 7.
    Best KO per CDS by full-sequence score; per-contig union of KOs for the module-completeness tool."""
    best = {}; per = defaultdict(set)
    with open(a.domtbl) as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.split()
            cds, ko, score = p[0], p[3], float(p[7])
            if cds not in best or score > best[cds][0]:
                best[cds] = (score, ko)
            per[contig_of(cds)].add(ko)      # MGnify's parsing_hmmscan unions EVERY KO hit passing --cut_ga per contig, not only the best per CDS
    pd.DataFrame([(c, k, s) for c, (s, k) in best.items()], columns=["cds_id", "ko_kofam", "full_sequence_score"]).to_csv(a.best_hits, sep="\t", index=False)
    with open(a.out, "w") as o:
        for contig in sorted(per):
            o.write("\t".join([contig] + sorted(per[contig])) + "\n")
    print(f"KOfam best hits {len(best)}, contigs with KOs {len(per)}")


def contig_of(cds_id):
    """Contig name of a CDS id (<contig>_<start>_<end>_<strand> from Prodigal/FragGeneScan, or <contig>_<n>)."""
    m = re.match(r"^(.*?)(?:_\d+_\d+_[+-]|_\d+)$", cds_id)
    return m.group(1) if m else cds_id


def build(a):
    """Assemble the query directory (contigs.tsv, cds.tsv, cds_taxonomy.tsv, modules.tsv, sample.json) from the pipeline's annotation files."""
    os.makedirs(a.out, exist_ok=True)
    contigs = []
    for name, seq in read_fasta(a.contigs):
        m = SPADES_RE.search(name); cov = float(m.group(3)) if m else ""
        contigs.append((name.split()[0], len(seq), cov))
    pd.DataFrame(contigs, columns=["contig_id", "length", "coverage"]).to_csv(f"{a.out}/contigs.tsv", sep="\t", index=False)
    cds_ids = [name.split()[0] for name, _ in read_fasta(a.cds)]
    ko = defaultdict(str)
    if a.emapper and os.path.exists(a.emapper):
        with open(a.emapper) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                p = line.rstrip("\n").split("\t")
                if len(p) > 8 and p[8]:
                    ko[p[0]] = ",".join(x.replace("ko:", "") for x in p[8].split(",") if x)
    pf = defaultdict(set)
    if a.ips and os.path.exists(a.ips):
        with open(a.ips) as f:
            for line in f:
                p = line.rstrip("\n").split("\t")
                if len(p) > 4 and re.match(r"PF\d+", p[4]):
                    pf[p[0]].add(p[4])
    kf = {}
    if a.kofam_best and os.path.exists(a.kofam_best):
        t = pd.read_csv(a.kofam_best, sep="\t"); kf = dict(zip(t.cds_id, t.ko_kofam))
    cds = pd.DataFrame({"cds_id": cds_ids, "contig_id": [contig_of(c) for c in cds_ids], "ko": [ko.get(c, "") for c in cds_ids],
                        "pfam": [",".join(sorted(pf.get(c, []))) for c in cds_ids], "ko_kofam": [kf.get(c, "") for c in cds_ids]})
    cds.to_csv(f"{a.out}/cds.tsv", sep="\t", index=False)
    d = pd.read_csv(a.diamond, sep="\t", dtype=str, keep_default_na=False)
    tax = d[d.tax_id.str.isdigit()][["contig_name", "tax_id"]].rename(columns={"contig_name": "cds_id", "tax_id": "taxid"})
    tax["contig_id"] = tax.cds_id.map(contig_of); tax[["cds_id", "contig_id", "taxid"]].to_csv(f"{a.out}/cds_taxonomy.tsv", sep="\t", index=False)
    if a.kegg and os.path.exists(a.kegg):
        k = pd.read_csv(a.kegg, sep="\t", dtype=str, keep_default_na=False)
        col = [c for c in k.columns if c.lower().startswith("module") or c.lower() == "pathway"][0]; comp = [c for c in k.columns if "complet" in c.lower()][0]
        pd.DataFrame({"module_id": k[col], "completeness": pd.to_numeric(k[comp], errors="coerce") / 100.0}).to_csv(f"{a.out}/modules.tsv", sep="\t", index=False)
    import json
    try:
        from checkgm import __version__ as checkgm_version
    except ImportError:
        checkgm_version = "unknown"
    json.dump({"sample_id": a.sample, "body_site": "Gut", "pipeline": {"name": "checkgm-assembly", "version": checkgm_version, "measurement": "MGnify pipeline v5.0 equivalent",
               "tools": {"assembler": os.environ.get("REFMB_ASSEMBLER", "metaSPAdes 3.15.3"), "cds": "Prodigal 2.6.3 + FragGeneScan 1.31", "taxonomy": "DIAMOND 0.9.25 / UniRef90 2019_11", "ko": "eggNOG-mapper 2.0.0 (MGnify v5 eggnog.db)",
                         "pfam": "InterProScan 5.36-75.0 Pfam", "kofam": "HMMER 3.2.1 / KOfam KEGG 90.0", "modules": "MGnify give_pathways"}}}, open(f"{a.out}/sample.json", "w"), indent=1)
    print(f"query dir {a.out}: contigs {len(contigs)}, CDS {len(cds)}, with KO {cds.ko.ne('').sum()}, with Pfam {cds.pfam.ne('').sum()}, taxonomy hits {len(tax)}")


def main():
    """Command line: `pv5_to_query.py <subcommand> ...` with the subcommands used by run_tierA.sh (filter-contigs, rename-megahit, merge-cds, join-taxonomy, kofam-union, build)."""
    ap = argparse.ArgumentParser(); sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("filter-contigs"); p.add_argument("--in", dest="inp", required=True); p.add_argument("--out", required=True); p.add_argument("--min-length", type=int, default=500)
    p = sub.add_parser("rename-megahit"); p.add_argument("--in", dest="inp", required=True); p.add_argument("--out", required=True)
    p = sub.add_parser("merge-cds"); p.add_argument("--prodigal-faa", required=True); p.add_argument("--fgs-faa"); p.add_argument("--out", required=True)
    p = sub.add_parser("join-taxonomy"); p.add_argument("--diamond", required=True); p.add_argument("--uniref-tax", required=True); p.add_argument("--out", required=True)
    p = sub.add_parser("kofam-union"); p.add_argument("--domtbl", required=True); p.add_argument("--out", required=True); p.add_argument("--best-hits", required=True)
    p = sub.add_parser("build"); p.add_argument("--sample", required=True); p.add_argument("--contigs", required=True); p.add_argument("--cds", required=True); p.add_argument("--diamond", required=True)
    p.add_argument("--emapper"); p.add_argument("--ips"); p.add_argument("--kofam-best"); p.add_argument("--kegg"); p.add_argument("--out", required=True)
    a = ap.parse_args()
    {"filter-contigs": filter_contigs, "rename-megahit": rename_megahit, "merge-cds": merge_cds, "join-taxonomy": join_taxonomy, "kofam-union": kofam_union, "build": build}[a.cmd](a)


if __name__ == "__main__":
    main()
