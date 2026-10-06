#!/usr/bin/env python3
"""Tidy tables for the v2 figures (drawn in R/ggplot2 by scripts/viz2/figures.R). Every number comes from a results file or
the score tables named here; nothing is typed by hand. Output: figures/v2/data/*.csv."""
import json, os, sys, gzip
import numpy as np, pandas as pd
# the project root is taken from the environment, as the README documents, and falls back to the tree this was written in
P = (os.environ.get("CHECKGM_PROJECT") or os.environ.get("REFMB_PROJECT") or "/home/classes/bios/270/khoa/meta2/project").rstrip("/")
OUT = f"{P}/figures/v2/data"; os.makedirs(OUT, exist_ok=True)
sys.path.insert(0, f"{P}/scripts/viz"); os.environ["REFMB_PROJECT"] = P   # the scripts/viz helpers read this one
import example_case as EC, fig17_range_report as F17
tax = pd.read_parquet(f"{P}/resources/backbone/ncbi_taxonomy.parquet", columns=["taxid", "name"]); NAME = dict(zip(tax.taxid.astype(str), tax.name))
LAYER = {"taxonomy_family": "Family", "taxonomy_genus": "Genus", "taxonomy_species": "Species", "ko_eggnog": "KO (eggNOG)", "ko_kofam": "KO (KOfam)", "pfam": "Pfam",
         "module": "KEGG module", "ko_humann": "KO", "pathway_humann": "Pathway"}

# ---- Fig 1: one sample's report (C. difficile case) and the healthy distribution of one family
r = F17.rows()
keep = ["Lachnospiraceae", "Oscillospiraceae", "Bifidobacteriaceae", "Bacteroidaceae", "Enterobacteriaceae", "Streptococcaceae"]
fam = r[r.group.str.startswith("Bacterial") & r.label.isin(keep)]
miss = r[r.call == "absent"]
ko = r[r.label.str.contains("K20707|K03088")]
f1 = pd.concat([fam, miss, ko]); f1["layer"] = np.where(f1.group.str.startswith("Bacterial"), "Families", "Gene families")
f1.to_csv(f"{OUT}/f1_report.csv", index=False)
ref = EC.reference(); DENS = "Oscillospiraceae"  # the case lies clearly below the range here (Lachnospiraceae sits on its edge)
fid = {v: k for k, v in NAME.items()}[DENS]
B = f"{P}/refs/gut-assembly-adult-global-v0.8-lenient"; pool = set()
inv = pd.read_parquet(f"{P}/work/s0/inventory.parquet", columns=["analysis_id", "study_bioproject", "is_primary_analysis", "body_site_core"])
st = set(pd.read_csv(f"{B}/provenance/studies.tsv", sep="\t").study_bioproject); excl = set(pd.read_csv(f"{B}/provenance/pool_exclusions.tsv", sep="\t").analysis_id)
pool = set(inv[inv.study_bioproject.isin(st - {EC.STUDY}) & inv.is_primary_analysis & (inv.body_site_core == "Gut") & ~inv.analysis_id.isin(excl)].analysis_id)
t = pd.read_parquet(f"{P}/work/s2/taxonomy_family_assembly.parquet", filters=[("feature_id", "==", fid)])
t = t[t.analysis_id.isin(pool) & (t.proportion_mapped > 0)]
med = ref.loc[fid, "p50"]; pd.DataFrame({"x": (t.clr - med) / np.log(2)}).to_csv(f"{OUT}/f1_density.csv", index=False)
json.dump({k: float((ref.loc[fid, f"p{k}"] - med) / np.log(2)) for k in ["1", "2_5", "25", "75", "97_5", "99"]} | {"sample": float(fam[fam.label == DENS].x.iloc[0]), "pct": float(fam[fam.label == DENS].pct.iloc[0]), "n": int(len(t))},
          open(f"{OUT}/f1_density_meta.json", "w"))

# ---- Fig 2: curation counts, calibration, reproduction from raw reads
invr = pd.read_parquet(f"{P}/work/s0/inventory.parquet", columns=["analysis_id", "study_bioproject", "s0_status", "is_primary_analysis", "body_site_core"])
ret = invr[invr.s0_status == "retained"]; gut = ret[ret.is_primary_analysis & (ret.body_site_core == "Gut")]
qc = pd.read_parquet(f"{P}/work/s2/qc_per_analysis.parquet", columns=["analysis_id", "s2_status"]); qok = qc[qc.s2_status == "retained"]
# the QC step is the only funnel row the figure could not label with a study count; it is the gut subset that passed, so
# the studies behind it are the studies of those analyses (every retained analysis is one of the 9,195 gut analyses)
gqc = gut[gut.analysis_id.isin(set(qok.analysis_id))]; assert len(gqc) == len(qok)
ma = json.load(open(f"{B}/manifest.json")); mb = json.load(open(f"{P}/refs/gut-reads-adult-global-v0.2-lenient/manifest.json"))
S = pd.read_csv(f"{P}/work/s11/pipeB/samples.tsv", sep="\t", low_memory=False)
adult = S[S.status.isin(["included", "NOT_CONTROL"])]; ctrl = S[S.status == "included"]
f2 = pd.DataFrame([
    ("Assembly (MGnify v5)", 1, "MGnify v5 analyses", len(invr), invr.study_bioproject.nunique()),
    ("Assembly (MGnify v5)", 2, "human metagenomes", len(ret), ret.study_bioproject.nunique()),
    ("Assembly (MGnify v5)", 3, "gut, one per sample", len(gut), gut.study_bioproject.nunique()),
    ("Assembly (MGnify v5)", 4, "pass assembly QC", len(qok), gqc.study_bioproject.nunique()),
    ("Assembly (MGnify v5)", 5, "healthy-adult baseline", ma["n_samples"], ma["n_studies"]),
    ("Read (cMD3)", 1, "cMD3 profiles", len(S), S.study.nunique()),
    ("Read (cMD3)", 2, "adults, profile resolvable", len(adult), adult.study.nunique()),
    ("Read (cMD3)", 3, "healthy controls", len(ctrl), ctrl.study.nunique()),
    ("Read (cMD3)", 5, "healthy-adult baseline", mb["n_samples"], mb["n_studies"])], columns=["pipeline", "step", "label", "n", "studies"])
f2.to_csv(f"{OUT}/f2_funnel.csv", index=False)
rep = json.load(open(f"{P}/results/s6/report_v08_lenient_base.json")); ps = pd.DataFrame(rep["per_study"])
cal = [pd.DataFrame({"pipeline": "Assembly (MGnify v5)", "kind": "held-out cohort", "layer": ps.layer.map(LAYER), "unit": ps.study_bioproject, "frac": ps.frac_outside_mean})]
lo = pd.read_csv(f"{P}/results/s6/loso_within_pool.tsv", sep="\t"); scol = [c for c in lo.columns if c.startswith("study")][0]
cal.append(pd.DataFrame({"pipeline": "Assembly (MGnify v5)", "kind": "baseline study left out", "layer": lo.layer.map(LAYER), "unit": lo[scol], "frac": lo.frac}))
hb = pd.read_csv(f"{P}/results/s11/pipelineB_heldout.tsv", sep="\t"); lb = pd.read_csv(f"{P}/results/s11/pipelineB_loso_all.tsv", sep="\t")
cal.append(pd.DataFrame({"pipeline": "Read (cMD3)", "kind": "held-out cohort", "layer": hb.layer.map(LAYER), "unit": hb.study, "frac": hb["mean"]}))
cal.append(pd.DataFrame({"pipeline": "Read (cMD3)", "kind": "baseline study left out", "layer": lb.layer.map(LAYER), "unit": lb.study, "frac": lb["median"]}))
pd.concat(cal).dropna(subset=["layer"]).to_csv(f"{OUT}/f2_calibration.csv", index=False)
ta = pd.read_csv(f"{P}/results/s11/tierA_anchor_compare.tsv", sep="\t"); ta = ta[ta.variant == "default"]
tb = pd.read_csv(f"{P}/results/s11/pipelineB_fidelity.tsv", sep="\t"); tb = tb[~tb.deposit_mismatch.astype(bool)]
fid_ = pd.concat([pd.DataFrame({"pipeline": "Assembly (MGnify v5)", "layer": ta.layer.map(LAYER), "sample": ta.anchor, "rho": ta.spearman_pct}),
                  pd.DataFrame({"pipeline": "Read (cMD3)", "layer": tb.layer.map(LAYER), "sample": tb.anchor, "rho": tb.spearman_pct})])
fid_.dropna().to_csv(f"{OUT}/f2_fidelity.csv", index=False)

# ---- Fig 3: AUROC per study and method
A = pd.read_csv(f"{P}/results/s8/loso_auroc.tsv", sep="\t"); A = A[A.study != "POOLED"]
# the assembly reference_relative set spans family, genus, eggNOG KO and module layers (s8_evaluate.LAYERS_RR), so it is
# "both"; s8 also stores the taxonomy-only and genes-only ablations, which give the same decomposition as the read pipeline
mapA = {"reference_relative": "checkGM percentiles, taxa + genes", "reference_relative_taxonomy_only": "checkGM percentiles",
        "reference_relative_genes_only": "checkGM percentiles, genes only", "raw_clr": "Raw abundances",
        "alpha_diversity": "Alpha diversity", "health_index": "GMHI (genus approx.)"}
A = A[A.feature_set.isin(mapA)]; fa = pd.DataFrame({"pipeline": "Assembly pipeline (11 studies)", "method": A.feature_set.map(mapA), "study": A.study, "auroc": A.auroc_loso, "unseen": True})
H = pd.read_csv(f"{P}/results/s15/health_indices_B.tsv", sep="\t"); base = H.drop_duplicates("study").set_index("study")
rows = [("checkGM percentiles", base.percentiles), ("Raw abundances", base.raw_abundances)]
for n, lab in [("GMHI", "GMHI (published)"), ("Shannon", "Alpha diversity"), ("GMWI2", "GMWI2 (published)")]:
    rows.append((lab, H[H["index"] == n].set_index("study").auroc_direct))
HS = pd.read_csv(f"{P}/results/s16/health_score_B.tsv", sep="\t"); rows.append(("checkGM health score", HS[HS.feature_set == "pct+summary+presence"].set_index("study").auroc))
# the same L1 model trained on raw CLR instead of percentiles, so the transformation is the only difference (s16 LAYERS=rawcmp)
RC = f"{P}/results/s16/health_score_B_rawcmp.tsv"
if os.path.exists(RC):
    rows.append(("same score on raw abundances", pd.read_csv(RC, sep="\t").query("feature_set == 'raw+summary+presence'").set_index("study").auroc))
# do gene families and pathways add anything? the same leave-one-study-out random forest on the function percentiles, and on
# taxa + function together (results/s11/pipelineB_loso_auroc.tsv, written by s13)
LF = pd.read_csv(f"{P}/results/s11/pipelineB_loso_auroc.tsv", sep="\t")
LF = LF[(LF.training == "all_studies") & (LF.features == "reference_relative")]
for fs, lab in [("function", "checkGM percentiles, genes only"), ("all", "checkGM percentiles, taxa + genes")]:
    rows.append((lab, LF[LF.feature_set == fs].set_index("study").auroc))
fb = pd.concat([pd.DataFrame({"method": lab, "study": s_.index, "auroc": s_.values}) for lab, s_ in rows])
fb["unseen"] = ~fb.study.map(base.in_gmwi2_training).astype(bool)
f3 = pd.concat([fa, fb.assign(pipeline="Read pipeline (14 studies)"), fb[fb.unseen].assign(pipeline="Read pipeline, 4 studies GMWI2 never saw")])
f3.to_csv(f"{OUT}/f3_auroc.csv", index=False)

# ---- Fig 4: C. difficile cohort (assembly) and colorectal cancer cohorts (read)
sc = pd.read_parquet(f"{P}/work/s8/scores/{EC.STUDY}.parquet")
inv2 = pd.read_parquet(f"{P}/work/s0/inventory.parquet", columns=["analysis_id", "cur_Health_status_group", "is_primary_analysis"]).set_index("analysis_id")
sc = sc[sc.analysis_id.map(inv2.is_primary_analysis).astype(bool)]; sc["group"] = sc.analysis_id.map(inv2.cur_Health_status_group).map({"Diseased": "C. difficile cases", "Healthy": "controls", "Control": "controls"})
want = {"taxonomy_family": ["Lachnospiraceae", "Oscillospiraceae", "Bifidobacteriaceae", "Enterobacteriaceae"], "taxonomy_genus": ["Lactobacillus"]}
rows = []
for layer, names in want.items():
    for n in names:
        f_ = {v: k for k, v in NAME.items()}[n]; ids = sc[(sc.layer == layer)].analysis_id.unique()
        g = sc[(sc.layer == layer) & (sc.feature_id.astype(str) == f_)].set_index("analysis_id")
        for a in ids:
            rows.append({"taxon": n, "sample": a, "group": sc[sc.analysis_id == a].group.iloc[0], "pct": float(g.percentile.get(a, 0.0)) if pd.notna(g.percentile.get(a, 0.0)) else 0.0})
f4a = pd.DataFrame(rows).dropna(subset=["group"]); f4a.to_csv(f"{OUT}/f4_cdi.csv", index=False)
F = pd.read_csv(f"{P}/results/s12/disease_atlas_B_features.tsv", sep="\t", dtype={"feature_id": str}); F = F[F.condition == "CRC"]
C = pd.read_csv(f"{P}/results/s12/disease_atlas_B_consistent.tsv", sep="\t", dtype={"feature_id": str}); C = C[(C.condition == "CRC") & C.layer.str.startswith("taxonomy")]
F = F.merge(C[["layer", "feature_id", "n_studies", "median_shift"]], on=["layer", "feature_id"])
F["name"] = F.feature_id.map(NAME); F["rank"] = F.layer.map(LAYER); F["sig"] = (F.q <= 0.05) & (F["shift"].abs() >= 10)
F[["study", "rank", "name", "shift", "sig", "n_studies", "median_shift"]].to_csv(f"{OUT}/f4_crc.csv", index=False)

# ---- Fig 5: Hadza hunter-gatherers and US participants of the same study vs the industrialized baseline (out-of-sample, assembly)
mA = pd.read_parquet(f"{P}/work/s12/A/meta.parquet").drop_duplicates("sample").set_index("sample")
SA = pd.read_parquet(f"{P}/work/s12/A/oos_summaries.parquet"); SA = SA[SA["sample"].map(mA.role).isin(["reference_pool", "heldout_healthy"])]
SA["group"] = np.where(SA["sample"].map(mA.study) != "PRJEB49206", "Other baseline studies",
                np.where(SA["sample"].map(mA.country) == "Tanzania", "Hadza (Tanzania)", np.where(SA["sample"].map(mA.country) == "United States", "US participants, same study", "Nepal")))
SA = SA[SA.group != "Nepal"]; SA["layer"] = SA.layer.map(LAYER)
SA[SA.layer.isin(["Family", "KO (eggNOG)", "Pfam"])][["sample", "group", "layer", "frac_outside_raw", "n_missing"]].to_csv(f"{OUT}/f5_frac.csv", index=False)
# genera carried by at least half of healthy adults in some size class (common gut genera; contig-level assignments of rare
# environmental genera are left out), ranked by the share of Hadza samples above minus below the healthy range
rg = pd.read_parquet(f"{B}/features/taxonomy_genus.parquet").set_index("feature_id"); pv = rg[[c for c in rg.columns if c.startswith("prevalence_band_")]].max(axis=1)
common = set(pv.index[pv >= 0.5].astype(str))
X = pd.read_parquet(f"{P}/work/s12/A/oos_scores.parquet", columns=["sample", "feature_id", "call"], filters=[("layer", "==", "taxonomy_genus")])
X = X[X["sample"].isin(set(SA["sample"])) & X.feature_id.astype(str).isin(common)]; X["group"] = X["sample"].map(SA.drop_duplicates("sample").set_index("sample").group)
n_by = SA.drop_duplicates("sample").groupby("group").size()
agg = X.assign(low=X.call == "low", high=X.call == "high", missing=X.call == "expected_but_missing").groupby(["group", "feature_id"])[["low", "high", "missing"]].sum()
agg = agg.div(n_by, level="group", axis=0).reset_index(); agg["name"] = agg.feature_id.astype(str).map(NAME); agg["n_group"] = agg.group.map(n_by)
hz = agg[agg.group == "Hadza (Tanzania)"].assign(score=lambda d: d.high - d.low - d.missing).set_index("feature_id")
top = list(hz.sort_values("score").tail(10).index) + list(hz.sort_values("score").head(5).index)
agg[agg.feature_id.isin(top)].assign(direction=lambda d: np.where(d.feature_id.isin(hz.sort_values("score").tail(10).index), "above", "below")).to_csv(f"{OUT}/f5_genera.csv", index=False)
print("written", sorted(os.listdir(OUT)))

# ---- Fig 4a: per-patient calls in the C. difficile study (assembly pipeline, scored against the baseline without the study).
# Families are chosen by a fixed rule: the 8 with the largest case-minus-control difference in the share below the range or
# missing, and the 4 with the largest difference in the share above it, among bacterial families (domain 2; contig-level hits to
# eukaryotic families such as Trichuridae are assignment artifacts).
BACT = set(pd.read_parquet(f"{P}/resources/backbone/ncbi_taxonomy.parquet", columns=["taxid", "domain_taxid"]).query("domain_taxid == 2").taxid.astype(str))
sc = pd.read_parquet(f"{P}/work/s8/scores/{EC.STUDY}.parquet", columns=["analysis_id", "layer", "feature_id", "call"], filters=[("layer", "==", "taxonomy_family")])
inv3 = pd.read_parquet(f"{P}/work/s0/inventory.parquet", columns=["analysis_id", "cur_Health_status_group", "is_primary_analysis"]).set_index("analysis_id")
sc = sc[sc.analysis_id.map(inv3.is_primary_analysis).astype(bool) & sc.feature_id.astype(str).isin(BACT)]
sc["group"] = sc.analysis_id.map(inv3.cur_Health_status_group).map({"Diseased": "cases", "Healthy": "controls"}); sc = sc.dropna(subset=["group"])
ng = sc.drop_duplicates("analysis_id").groupby("group").size()
t = sc.groupby(["feature_id", "group"]).call.value_counts().unstack(fill_value=0).reindex(columns=["low", "expected_but_missing", "high"], fill_value=0)
t = t.div(t.index.get_level_values("group").map(ng), axis=0).reset_index(); t["name"] = t.feature_id.astype(str).map(NAME)
w = t.pivot(index="name", columns="group", values=["low", "expected_but_missing", "high"]).fillna(0)
dlow = (w["low"]["cases"] + w["expected_but_missing"]["cases"]) - (w["low"]["controls"] + w["expected_but_missing"]["controls"]); dhigh = w["high"]["cases"] - w["high"]["controls"]
pick = list(dlow.sort_values().tail(8).index) + list(dhigh.sort_values().tail(4).index)
t[t.name.isin(pick)].assign(n_group=lambda d: d.group.map(ng), side=lambda d: np.where(d.name.isin(dhigh.sort_values().tail(4).index), "gained", "lost")).to_csv(f"{OUT}/f4_cdi_calls.csv", index=False)
print("cdi picks", pick, ng.to_dict())

# ---- Fig 4c: the gene families and pathways that shift in several colorectal cancer cohorts (read pipeline). MetaCyc names
# come from the HUMAnN tables themselves ("PWY-xxxx: name"); KO names from HUMAnN's map_ko_name.
import glob, gzip
PNAME = {}
for f_ in sorted(glob.glob(f"{P}/staging/cmd3_humann/pathway_abundance/*.parquet"))[:4]:
    for v in pd.read_parquet(f_, columns=["feature_id"]).feature_id.unique():
        if ":" in str(v):
            PNAME.setdefault(str(v).split(":")[0].strip(), str(v).split(":", 1)[1].strip())
KNAME = {}
with gzip.open(f"{P}/staging/cmd3_humann/mapping/map_ko_name.txt.gz", "rt") as fh:
    for line in fh:
        q = line.rstrip("\n").split("\t"); KNAME[q[0]] = q[1] if len(q) > 1 else q[0]
Cf = pd.read_csv(f"{P}/results/s12/disease_atlas_B_consistent.tsv", sep="\t", dtype={"feature_id": str})
Cf = Cf[(Cf.condition == "CRC") & Cf.layer.str.endswith("_humann")]
Ff = pd.read_csv(f"{P}/results/s12/disease_atlas_B_features.tsv", sep="\t", dtype={"feature_id": str})
Ff = Ff[(Ff.condition == "CRC") & Ff.layer.str.endswith("_humann")].merge(Cf[["layer", "feature_id", "n_studies", "median_shift"]], on=["layer", "feature_id"])
pick = []
for lay, n_min, k in [("pathway_humann", 3, 7), ("ko_humann", 2, 5)]:   # pathways replicate more widely than single KOs
    g = Cf[(Cf.layer == lay) & (Cf.n_studies >= n_min)].copy()
    pick += list(g.reindex(g.median_shift.abs().sort_values(ascending=False).index).head(k).feature_id)
Ff = Ff[Ff.feature_id.isin(pick)].copy()
Ff["name"] = [PNAME.get(f_, f_) if lay == "pathway_humann" else KNAME.get(f_, f_) for f_, lay in zip(Ff.feature_id, Ff.layer)]
Ff["rank"] = Ff.layer.map({"pathway_humann": "pathways", "ko_humann": "gene families"})
Ff["sig"] = (Ff.q <= 0.05) & (Ff["shift"].abs() >= 10)
Ff[["study", "rank", "name", "shift", "sig", "n_studies", "median_shift"]].to_csv(f"{OUT}/f4_crc_function.csv", index=False)

# ---- Fig 5c: the gene families that place Hadza samples outside the range (assembly pipeline, KOs carried by most healthy adults)
rk = pd.read_parquet(f"{B}/features/ko_eggnog.parquet").set_index("feature_id")
pvk = rk[[c for c in rk.columns if c.startswith("prevalence_band_")]].max(axis=1); commonk = set(pvk.index[pvk >= 0.9].astype(str))
K = pd.read_parquet(f"{P}/work/s12/A/oos_scores.parquet", columns=["sample", "feature_id", "call"], filters=[("layer", "==", "ko_eggnog")])
K = K[K["sample"].isin(set(SA["sample"])) & K.feature_id.astype(str).isin(commonk)]
K["group"] = K["sample"].map(SA.drop_duplicates("sample").set_index("sample").group)
aggk = K.assign(low=K.call == "low", high=K.call == "high", missing=K.call == "expected_but_missing").groupby(["group", "feature_id"])[["low", "high", "missing"]].sum()
aggk = aggk.div(n_by, level="group", axis=0).reset_index(); aggk["name"] = aggk.feature_id.astype(str).map(lambda k_: KNAME.get(k_, k_))
aggk = aggk[~aggk.name.str.contains("uncharacterized|hypothetical", case=False)]   # unnamed orthologs carry no readable label
hk = aggk[aggk.group == "Hadza (Tanzania)"].assign(score=lambda d: d.high - d.low - d.missing).set_index("feature_id")
topk = list(hk.sort_values("score").tail(4).index) + list(hk.sort_values("score").head(3).index)
aggk[aggk.feature_id.isin(topk)].assign(direction=lambda d: np.where(d.feature_id.isin(hk.sort_values("score").tail(4).index), "above", "below")).to_csv(f"{OUT}/f5_ko.csv", index=False)
print("function panels written")
