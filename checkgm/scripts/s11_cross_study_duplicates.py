#!/usr/bin/env python3
"""S11.x — audit the read-pipeline sample table for specimens that appear in more than one cMD3 study.

  s11_cross_study_duplicates.py [--out work/s11/pipeB/cross_study_duplicates.tsv] [--report]

Two independent tests are run over every pipeline-B sample (roles reference_pool, heldout_healthy, case_control) and
their results are unioned, because neither test is a superset of the other:

  run_accession  two samples of different studies share at least one ENA/SRA run accession (cMD3 NCBI_accession,
                 ';'-separated).  Catches deposits that reuse the same run records.
  fingerprint    two samples of different studies report an identical (number_reads, number_bases) pair.  Catches
                 deposits of the same sequencing data under *disjoint* accession sets, which the accession test
                 cannot see.

The second test is why this script exists.  LeChatelierE_2013 and NielsenHB_2014 are two cMD3 deposits of the same
MetaHIT Danish cohort: 177 specimens appear in both under the same cMD3 sample_id, subject_id, country and
(number_reads, number_bases), but only MH0001-MH0005 share any run accession, so an accession-only audit records 5
pairs out of 177.  The pairs are not byte-identical copies - the two deposits profile the same specimen from
differently partitioned run sets, so the MetaPhlAn profiles differ slightly - which is why the report also prints the
species-level Bray-Curtis similarity of each matched pair next to the best match that sample has in any other study.

Output columns (symmetric: one row per ordered pair, as the consumers expect)
  sample_key duplicate_of study duplicate_study role duplicate_role match bray_curtis
Consumers: s12_oos_scores.py (a study with a duplicate in the pool is scored against the baseline rebuilt without
that study) and s16_health_score_B.py (studies sharing specimens form one group in the inner cross-validation)."""
import argparse, itertools, os, sys
import numpy as np, pandas as pd

P = "/home/classes/bios/270/khoa/meta2/project"
CMD3 = "/home/classes/bios/270/khoa/meta2/cmd3/cmd3_release_sampleMetadata.csv"
ROLES = ["reference_pool", "heldout_healthy", "case_control"]
# cMD3 renamed two studies between the release this project read and the release metadata table (decisions.yaml,
# pipelineB_extension).  Without the aliases the metadata join silently drops 108 samples and they are never audited.
# HMP_2019 needs no alias: its project sample ids already carry the "ibdmdb_" prefix, so the plain key matches.
ALIAS = {"ThomasAM_2019_a": "ThomasAM_2018a", "ThomasAM_2019_b": "ThomasAM_2018b"}
ap = argparse.ArgumentParser()
ap.add_argument("--out", default=f"{P}/work/s11/pipeB/cross_study_duplicates.tsv")
ap.add_argument("--report", action="store_true", help="print the per-pair Bray-Curtis table")
a = ap.parse_args()

S = pd.read_csv(f"{P}/work/s11/pipeB/samples.tsv", sep="\t", dtype={"sample_key": str}, low_memory=False)
S = S[S.role.isin(ROLES)][["sample_key", "study", "role"]]
C = pd.read_csv(CMD3, dtype=str, keep_default_na=False)
C["cmd3_key"] = C.study_name + "_" + C.sample_id
known = set(C.cmd3_key)   # the plain key first, the renamed key only for the samples the plain key misses
S["cmd3_key"] = [k if k in known else ALIAS.get(st, st) + k[len(st):] for k, st in zip(S.sample_key, S.study)]
M = S.merge(C[["cmd3_key", "sample_id", "subject_id", "country", "NCBI_accession", "number_reads", "number_bases"]], on="cmd3_key", how="left")
unaud = M[M.NCBI_accession.isna()]
if len(unaud):
    print(f"WARNING: {len(unaud)} of {len(S)} pipeline-B samples have no row in {CMD3} and cannot be audited:")
    print(unaud.groupby(["study", "role"]).size().to_string())
M = M[M.NCBI_accession.notna()]

# ---- test 1: shared run accession ------------------------------------------------------------------------------
acc = {}
for k, st, v in zip(M.sample_key, M.study, M.NCBI_accession):
    for r in str(v).split(";"):
        r = r.strip()
        if r:
            acc.setdefault(r, []).append((k, st))
pairs = {}   # frozenset({a, b}) -> set of match kinds
for r, mem in acc.items():
    for (ka, sa), (kb, sb) in itertools.combinations(sorted(set(mem)), 2):
        if sa != sb:
            pairs.setdefault(frozenset((ka, kb)), set()).add("run_accession")

# ---- test 2: identical sequencing fingerprint ------------------------------------------------------------------
F = M[(M.number_reads.str.strip() != "") & (M.number_reads.str.strip().str.upper() != "NA") & (M.number_bases.str.strip() != "")]
for _, g in F.groupby(["number_reads", "number_bases"]):
    if g.study.nunique() < 2:
        continue
    for (ka, sa), (kb, sb) in itertools.combinations(sorted(zip(g.sample_key, g.study)), 2):
        if sa != sb:
            pairs.setdefault(frozenset((ka, kb)), set()).add("fingerprint")

# ---- species-level Bray-Curtis for each matched pair, and the best match elsewhere -----------------------------
sp = pd.read_parquet(f"{P}/work/s11/metaphlan/species_long.parquet")
need = sorted({k for p in pairs for k in p})
W = sp[sp.sample_key.isin(set(need))].pivot_table(index="sample_key", columns="species", values="rel_abundance", fill_value=0.0)
W = W.div(W.sum(axis=1).replace(0, np.nan), axis=0)


def bc_sim(ka, kb):
    if ka not in W.index or kb not in W.index:
        return np.nan
    x, y = W.loc[ka].to_numpy(), W.loc[kb].to_numpy()
    return float(np.minimum(x, y).sum() / ((x.sum() + y.sum()) / 2))


role = dict(zip(M.sample_key, M.role)); study = dict(zip(M.sample_key, M.study))
rows = []
for p, kinds in sorted(pairs.items(), key=lambda kv: sorted(kv[0])):
    ka, kb = sorted(p)
    m = "both" if len(kinds) == 2 else next(iter(kinds))
    s = bc_sim(ka, kb)
    for x, y in ((ka, kb), (kb, ka)):
        rows.append({"sample_key": x, "duplicate_of": y, "study": study[x], "duplicate_study": study[y],
                     "role": role[x], "duplicate_role": role[y], "match": m, "bray_curtis": None if s is None or np.isnan(s) else round(s, 4)})
D = pd.DataFrame(rows).sort_values(["study", "duplicate_study", "sample_key"])
os.makedirs(os.path.dirname(a.out), exist_ok=True)
D.to_csv(a.out, sep="\t", index=False)

print(f"{len(D) // 2} cross-study duplicate pairs over {len(S)} pipeline-B samples -> {a.out}")
t = D.groupby(["study", "duplicate_study", "match"]).agg(pairs=("sample_key", "size"), median_bray_curtis=("bray_curtis", "median")).reset_index()
t["pairs"] = t.pairs  # one row per ordered pair, so this is already per-direction
print(t.to_string(index=False))
print("\nby role of the two members:")
print(D.groupby(["role", "duplicate_role"]).size().to_string())
# how much of each study is a duplicate of another study
own = D.drop_duplicates("sample_key").groupby("study").size().rename("duplicated_samples")
tot = S.groupby("study").size().rename("samples")
print("\nshare of each study that is a duplicate of another study:")
print(pd.concat([tot, own], axis=1).dropna().assign(share=lambda d: (d.duplicated_samples / d.samples).round(3)).to_string())
if a.report:
    print("\nper-pair Bray-Curtis similarity, matched pair vs the best match the sample has in any other study:")
    other = {}
    pool = sorted(set(W.index) | set())
    for k in W.index:
        x = W.loc[k].to_numpy()
        best, who = -1.0, ""
        for j in W.index:
            if j == k or study[j] == study[k]:
                continue
            if frozenset((k, j)) in pairs:
                continue
            y = W.loc[j].to_numpy(); v = float(np.minimum(x, y).sum() / ((x.sum() + y.sum()) / 2))
            if v > best:
                best, who = v, j
        other[k] = (round(best, 4), who)
    R = D.drop_duplicates("sample_key").copy()
    R["best_other"] = R.sample_key.map(lambda k: other[k][0]); R["best_other_sample"] = R.sample_key.map(lambda k: other[k][1])
    print(R[["sample_key", "duplicate_of", "match", "bray_curtis", "best_other", "best_other_sample"]].to_string(index=False))
