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
# The cohort behind this panel is the baseline rebuilt WITHOUT the patient's own study (EC.BUNDLE), not the shipped
# baseline: 1,927 adults from 25 studies, of whom n_detected carry the family. The shipped baseline is a different and
# larger cohort (1,941 adults, 26 studies, 1,936 of them carrying it), so the two must not be paired in a caption.
# n_drawn is the number of points the density is estimated from and should track n_detected.
mref = json.load(open(f"{EC.BUNDLE}/manifest.json"))
json.dump({k: float((ref.loc[fid, f"p{k}"] - med) / np.log(2)) for k in ["1", "2_5", "25", "75", "97_5", "99"]} | {"sample": float(fam[fam.label == DENS].x.iloc[0]), "pct": float(fam[fam.label == DENS].pct.iloc[0]), "n": int(len(t)),
           "n_drawn": int(len(t)), "n_detected": int(ref.loc[fid, "n_detected"]), "n_pool": int(mref["n_samples"]), "n_pool_studies": int(mref["n_studies"])},
          open(f"{OUT}/f1_density_meta.json", "w"))
print(f"fig1 density: {DENS} carried by {int(ref.loc[fid, 'n_detected'])} of {mref['n_samples']} adults in the LOSO baseline ({mref['n_studies']} studies); {len(t)} points drawn")

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
# All four series in this panel are the same per-study statistic: the MEAN over the study's samples of the per-sample share
# of features outside the range. It has to be the mean, because that is the only summary the assembly tables carry (s6
# frac_outside_mean; s6_loso_from_s12.py frac=("frac_outside_raw","mean")); the read tables carry both (s13_pipelineB_eval.py
# frac_table) and lb["median"] was taken here until 2026-10-06, which mixed medians into a panel of means. The per-sample
# distribution is right-skewed (mean > median on 93/93 read leave-one-study-out rows), so the mix displaced the read
# "baseline study left out" points left of every other series: family 3.92 % against 5.60 % on the matched statistic.
cal.append(pd.DataFrame({"pipeline": "Read (cMD3)", "kind": "held-out cohort", "layer": hb.layer.map(LAYER), "unit": hb.study, "frac": hb["mean"]}))
cal.append(pd.DataFrame({"pipeline": "Read (cMD3)", "kind": "baseline study left out", "layer": lb.layer.map(LAYER), "unit": lb.study, "frac": lb["mean"]}))
pd.concat(cal).dropna(subset=["layer"]).to_csv(f"{OUT}/f2_calibration.csv", index=False)
ta = pd.read_csv(f"{P}/results/s11/tierA_anchor_compare.tsv", sep="\t"); ta = ta[ta.variant == "default"]
tb = pd.read_csv(f"{P}/results/s11/pipelineB_fidelity.tsv", sep="\t"); tb = tb[~tb.deposit_mismatch.astype(bool)]
fid_ = pd.concat([pd.DataFrame({"pipeline": "Assembly (MGnify v5)", "layer": ta.layer.map(LAYER), "sample": ta.anchor, "rho": ta.spearman_pct}),
                  pd.DataFrame({"pipeline": "Read (cMD3)", "layer": tb.layer.map(LAYER), "sample": tb.anchor, "rho": tb.spearman_pct})])
fid_.dropna().to_csv(f"{OUT}/f2_fidelity.csv", index=False)

# ---- Fig 3: AUROC per study and method
# The assembly pipeline is not plotted here. checkGM ships no assembly deviation score --- the shipped score is defined
# on the read pipeline's species percentiles and refuses an assembly assessment --- so every assembly bar would be a
# model nobody can run. Its representation check (percentiles against raw centred log-ratios, same logistic regression)
# is reported in the text from results/s8/report_loso_l1.md.
H = pd.read_csv(f"{P}/results/s15/health_indices_B.tsv", sep="\t"); base = H.drop_duplicates("study").set_index("study")
# Every bar is a logistic regression or a published score evaluated as published: no random forest anywhere, and no
# checkGM variant that the software does not distribute. The deviation score plotted is therefore the percentile-only
# model that `checkgm score` applies, not the pre-specified variant whose species-presence features carry the
# [Collinsella] batch artifact; that one is reported in the text. Both it and the raw-abundance comparator are taken
# from the one rawcmp run, so the pair is matched by construction rather than assembled from two runs.
RC = pd.read_csv(f"{P}/results/s16/health_score_B_rawcmp.tsv", sep="\t")
rows = [("checkGM deviation score", RC.query("feature_set == 'pct'").set_index("study").auroc),
        ("raw abundances, same model", RC.query("feature_set == 'raw'").set_index("study").auroc)]
for n, lab in [("GMWI2", "GMWI2 (published)"), ("GMHI", "GMHI (published)"), ("Shannon", "Alpha diversity")]:
    rows.append((lab, H[H["index"] == n].set_index("study").auroc_direct))
# The fraction outside, the one deviation summary that needs no labels at all, scored the same way the published indices
# are: directly, higher meaning more case-like. Species is the layer it does best on (family 0.52, genus 0.53), so the
# comparison is not stacked against it. It is in the figure because its distance from the score is the paper's argument
# for per-feature output: the same percentiles, collapsed to one number, lose almost all of their discrimination.
from sklearn.metrics import roc_auc_score
_elig = S[S.role.isin(["reference_pool", "heldout_healthy", "case_control"]) & S.status.isin(["included", "NOT_CONTROL"])]
_lab = _elig.set_index("sample_key")[["condition", "study"]].assign(y=lambda d: d.condition.ne("control").astype(int))
_fo = (pd.read_parquet(f"{P}/work/s12/B/oos_summaries.parquet")
         .query("layer == 'taxonomy_species'").set_index("sample").frac_outside_raw)
_rows = {}
for _st in base.index:
    _ix = _lab.index[_lab.study == _st]
    _y, _x = _lab.loc[_ix, "y"], _fo.reindex(_ix)
    _ok = _x.notna()
    if _ok.sum() >= 20 and _y[_ok].nunique() == 2:
        _rows[_st] = roc_auc_score(_y[_ok], _x[_ok])
rows.append(("fraction outside", pd.Series(_rows)))
fb = pd.concat([pd.DataFrame({"method": lab, "study": s_.index, "auroc": s_.values}) for lab, s_ in rows])
fb["unseen"] = ~fb.study.map(base.in_gmwi2_training).astype(bool)
f3 = pd.concat([fb.assign(pipeline="Read pipeline (14 studies)"),
                fb[fb.unseen].assign(pipeline="Read pipeline, 4 studies GMWI2 never saw")])
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

# ---- Figs 3 and 4 as report views: the same geometry as Fig. 1c, carried into the two case studies.
# In percentile space the reference range is 2.5--97.5 for every feature, so one shaded band serves every row and a
# point outside it is outside the reference range by construction. That is what these two files are for: a percentile
# and a call per sample per feature, taxa and function together, rather than the shares and shifts the panels used
# before, which showed the deviation without ever showing the range it was a deviation from.

# read the cohort's scores again rather than reuse `sc`, which is rebound further up to a reduced frame carrying
# neither percentiles nor groups; a panel that silently depends on a name defined 100 lines away is how this broke once
_cdi = pd.read_parquet(f"{P}/work/s8/scores/{EC.STUDY}.parquet",
                       columns=["analysis_id", "layer", "feature_id", "percentile", "call"])
_cdi = _cdi[_cdi.analysis_id.map(inv2.is_primary_analysis).astype(bool)].copy()
_cdi["group"] = _cdi.analysis_id.map(inv2.cur_Health_status_group).map(
    {"Diseased": "C. difficile cases", "Healthy": "controls", "Control": "controls"})
_cdi = _cdi[_cdi.group.notna()]
assert {"percentile", "call", "group"} <= set(_cdi.columns) and len(_cdi), "C. difficile scores did not load"
# Restricted to features the reference population actually carries, as everything else in the paper is. Without it the
# ranking is won by artifacts: a rare environmental family, a bacteriophage (Herelleviridae) and a whipworm
# (Trichuridae) all separate cases from controls perfectly, as do single-hit effector genes at thin coverage, because a
# feature almost nobody carries is "outside the range" for almost everybody. The thresholds are the ones the Hadza
# panels use: half the reference adults in some size class for taxa, nine tenths for gene families.
# Two restrictions, each for a reason the data forced. Bacteria only (domain_taxid 2), because the contig-level
# assignments include a bacteriophage family and Trichuridae, the whipworms, which a panel headed "bacterial families"
# must not contain --- Trichuridae is not even rare, carried by 1,829 of the 1,941 reference adults, which is homology
# at the contig level rather than helminth DNA. And study-weighted prevalence rather than the most permissive size
# band, since taking the maximum over bands admitted Herelleviridae at 0.34 weighted on the strength of 0.51 in one
# band alone.
_dom = pd.read_parquet(f"{P}/resources/backbone/ncbi_taxonomy.parquet", columns=["taxid", "domain_taxid"])
_bact = set(_dom.query("domain_taxid == 2").taxid.astype(str))
_rf = pd.read_parquet(f"{B}/features/taxonomy_family.parquet").set_index("feature_id")
_commonf = set(_rf.index[_rf.prevalence_weighted >= 0.5].astype(str)) & _bact
_cdi = _cdi[((_cdi.layer == "taxonomy_family") & _cdi.feature_id.astype(str).isin(_commonf))
            | ((_cdi.layer == "ko_eggnog") & _cdi.feature_id.astype(str).isin(commonk))]
# Features are chosen by the rule the share panel already used and the text already describes: the eight whose loss
# (below the range or absent) most separates patients from controls, and the four most gained. Re-deriving a selection
# instead produced a different, unnamed set each time the rule was tweaked, and the prose would then have described a
# panel that was not there.
# Gene families have no published list, so they are ranked the way the families were: most lost, then most gained.
_k = _cdi[(_cdi.layer == "ko_eggnog") & _cdi.feature_id.astype(str).isin(commonk)]
_ng = _k.drop_duplicates("analysis_id").groupby("group").size()
_c = _k.groupby(["feature_id", "group"]).call.value_counts().unstack(fill_value=0) \
       .reindex(columns=["low", "expected_but_missing", "high"], fill_value=0)
_c = _c.div(_c.index.get_level_values("group").map(_ng), axis=0).reset_index()
_w = _c.pivot(index="feature_id", columns="group", values=["low", "expected_but_missing", "high"]).fillna(0)
_named = [f for f in _w.index if KNAME.get(str(f)) and "uncharacterized" not in KNAME[str(f)].lower()]
_w = _w.loc[_named]
_lost = ((_w["low"]["C. difficile cases"] + _w["expected_but_missing"]["C. difficile cases"])
         - (_w["low"]["controls"] + _w["expected_but_missing"]["controls"])).sort_values()
_gain = (_w["high"]["C. difficile cases"] - _w["high"]["controls"]).sort_values()
_pick_ko = list(_lost.tail(3).index) + list(_gain.tail(2).index)

# the published panel already records which side each family was picked for; carrying it through means the row's
# label states the reason the row is there, rather than whichever share happens to be larger
_pub = pd.read_csv(f"{OUT}/f4_cdi_calls.csv").drop_duplicates("name")
_side_of = dict(zip(_pub["name"], _pub["side"]))
_picked = _pub["name"].tolist()
_byname = {v: k for k, v in NAME.items()}
_rows = []
_ko_side = {str(f): ("lost" if i < 3 else "gained") for i, f in enumerate(_pick_ko)}
for layer, fids, blk in (("taxonomy_family", [_byname[n] for n in _picked if n in _byname], "Bacterial families"),
                         ("ko_eggnog", [str(x) for x in _pick_ko], "Gene families")):
    for fid in dict.fromkeys(str(x) for x in fids):
        g = _cdi[(_cdi.layer == layer) & (_cdi.feature_id.astype(str) == fid)]
        nm = NAME.get(fid) or KNAME.get(fid, fid)
        for r in g.itertuples():
            _rows.append({"block": blk, "name": nm, "sample": r.analysis_id, "group": r.group,
                          "pct": r.percentile, "call": r.call,
                          "side": _side_of.get(nm) or _ko_side.get(fid, "lost")})
f3r = pd.DataFrame(_rows)
f3r.to_csv(f"{OUT}/f3_range_cdi.csv", index=False)
print(f"f3 report view: {f3r['name'].nunique()} features, {f3r['sample'].nunique()} samples, "
      f"{f3r.groupby('block')['name'].nunique().to_dict()}")

_sg = SA.drop_duplicates("sample").set_index("sample").group
# The features are the ones the share panels and the paper's text already name, ranked by net directional deviation
# (high - low - missing) rather than by total outside-ness: ranking by the total picks a different, unnamed set and the
# figure would then contradict the sentence describing it.
_pick = {"taxonomy_genus": (list(hz.sort_values("score").tail(4).index) + list(hz.sort_values("score").head(3).index),
                            "Common gut genera"),
         "ko_eggnog": (list(hk.sort_values("score").tail(3).index) + list(hk.sort_values("score").head(2).index),
                       "Gene families")}
_hz = []
for layer, (_ids, blk) in _pick.items():
    _ids = [str(x) for x in _ids]
    Z = pd.read_parquet(f"{P}/work/s12/A/oos_scores.parquet", columns=["sample", "feature_id", "percentile", "call"],
                        filters=[("layer", "==", layer)])
    Z = Z[Z["sample"].isin(set(SA["sample"])) & Z.feature_id.astype(str).isin(set(_ids))].copy()
    Z["group"] = Z["sample"].map(_sg)
    Z["name"] = Z.feature_id.astype(str).map(lambda f_: NAME.get(f_) or KNAME.get(f_, f_))
    Z["block"] = blk
    _hz.append(Z[Z.group.notna()][["block", "name", "sample", "group", "pct" if "pct" in Z else "percentile", "call"]]
               .rename(columns={"percentile": "pct"}))
f4r = pd.concat(_hz, ignore_index=True)
f4r.to_csv(f"{OUT}/f4_range_hadza.csv", index=False)
print(f"f4 report view: {f4r.groupby('block')['name'].nunique().to_dict()}, {f4r['sample'].nunique()} samples")

# The colorectal panel on the same axis as the rest: a cohort's case median and control median are both percentiles of
# one baseline, so they can be drawn against the same reference band instead of being reduced to the difference between
# them. Taxa and function together, restricted as before to features the reference population carries.
_A = pd.read_csv(f"{P}/results/s12/disease_atlas_B_features.tsv", sep="\t", dtype={"feature_id": str})
_A = _A[_A.condition == "CRC"]
_Ac = pd.read_csv(f"{P}/results/s12/disease_atlas_B_consistent.tsv", sep="\t", dtype={"feature_id": str})
_Ac = _Ac[(_Ac.condition == "CRC") & (_Ac.n_studies >= 3)]
_A = _A.merge(_Ac[["layer", "feature_id", "n_studies", "median_shift"]], on=["layer", "feature_id"])
_A["block"] = np.where(_A.layer.str.startswith("taxonomy"), "Taxa", "Function")
_A["name"] = np.where(_A.layer.str.startswith("taxonomy"), _A.feature_id.map(NAME),
                      _A.feature_id.map(lambda k_: PNAME.get(k_) or KNAME.get(k_) or k_))
_A = _A[_A["name"].notna()]
_keep = (_A.groupby(["block", "name"]).median_shift.first().abs().sort_values(ascending=False)
           .groupby(level=0).head(6).index)
_A = _A.set_index(["block", "name"]).loc[_keep.unique()].reset_index()
_A["sig"] = (_A.q <= 0.05) & (_A["shift"].abs() >= 10)
# median_shift and n_studies come from the consistent-shift table and are the pair the running text quotes
# ("Lachnospiraceae in four cohorts at a median of -35 points"); a median taken over all nine cohorts here instead
# would put a different number beside the same name.
_A[["block", "name", "study", "median_pct_case", "median_pct_control", "shift", "sig", "n_studies", "median_shift"]] \
  .to_csv(f"{OUT}/f3_range_crc.csv", index=False)
print(f"f3 crc view: {_A.groupby('block')['name'].nunique().to_dict()}, {_A.study.nunique()} cohorts")
