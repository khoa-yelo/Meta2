#!/usr/bin/env python3
"""A supervised health-versus-disease score built on refmb output (read pipeline), evaluated leave-one-study-out.

Samples: all cMD3 adults of the read pipeline — baseline healthy adults (6,494), held-out healthy adults (317), and the cases
and controls of the disease studies (label 1 = case). Features (per sample):
  pct      percentiles of family, genus and species against the healthy baseline (out-of-sample: a baseline study's samples
           are scored against a baseline rebuilt without it; work/s12/B/oos_scores.parquet); undetected = 0
  summary  fraction outside the healthy range and number of expected-but-missing features, per taxonomy layer
  presence presence (relative abundance > 1e-5) of each MetaPhlAn 3 species detected in >= 0.5 % of samples — this carries
           the rare, mostly oral-origin taxa that are below the baseline's prevalence filter
Model: L1-penalised logistic regression (liblinear, balanced class weights) on standardized features; the penalty C is chosen
inside each training set by grouped 5-fold CV over studies. Outer loop: each of the 14 case/control studies used in the paper
is left out; the model is trained on every other study (healthy-only, case-only and case/control studies alike).
Feature sets: presence | pct | pct+summary | pct+summary+presence. With LAYERS=taxfun the run is restricted to the samples
that also carry HUMAnN tables (8,465 of 9,097; 1,375 of 1,376 evaluation samples) and each set is run twice, with and without
the gene-family and pathway percentiles, so that adding function is the only difference between the two. Writes
results/s16/health_score_B{suffix}.{tsv,md} and the coefficients of the final model (…_coefficients.tsv)."""
import os, sys, time
import numpy as np, pandas as pd
from joblib import Parallel, delayed
from scipy.stats import wilcoxon
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

P = "/home/classes/bios/270/khoa/meta2/project"; OUT = f"{P}/results/s16"; os.makedirs(OUT, exist_ok=True); t0 = time.time()
NJ = int(os.environ.get("SLURM_CPUS_PER_TASK", 8)); CS = [0.003, 0.01, 0.03, 0.1]
S = pd.read_csv(f"{P}/work/s11/pipeB/samples.tsv", sep="\t", low_memory=False)
S = S[S.role.isin(["reference_pool", "heldout_healthy", "case_control"]) & S.status.isin(["included", "NOT_CONTROL"])].copy()
S["y"] = S.condition.ne("control").astype(int)
TRAIN = os.environ.get("TRAIN", "all"); LAY = os.environ.get("LAYERS", "tax")   # tax = taxonomy only (primary) | taxfun = also KO and pathway percentiles | rawcmp = percentiles vs raw CLR, same model
SUF = ("" if TRAIN == "all" else f"_{TRAIN}") + ("" if LAY == "tax" else f"_{LAY}")
if LAY == "taxfun":   # same sample set for every feature set, so that "with function" vs "without" is the only contrast
    _f = pd.read_parquet(f"{P}/work/s12/B/oos_summaries.parquet", columns=["sample", "layer"])
    _ok = set(_f[_f.layer == "ko_humann"]["sample"]) & set(_f[_f.layer == "pathway_humann"]["sample"])
    S = S[S.sample_key.isin(_ok)].copy(); print(f"restricted to {len(S)} samples with HUMAnN tables", flush=True)
if TRAIN == "mixed":   # only studies with >= 10 cases and >= 10 controls, so study identity cannot predict the label
    cnt = S.groupby("study").y.agg(["sum", "size"]); ok = cnt.index[(cnt["sum"] >= 10) & (cnt["size"] - cnt["sum"] >= 10)]
    EV = pd.read_csv(f"{P}/results/s11/pipelineB_case_control.tsv", sep="\t").study
    S = S[S.study.isin(set(ok) | set(EV))].copy(); print("training studies with both classes:", len(ok), sorted(ok), flush=True)   # evaluation studies always kept
EVAL = pd.read_csv(f"{P}/results/s11/pipelineB_case_control.tsv", sep="\t").study.tolist()
X = pd.read_parquet(f"{P}/work/s12/B/oos_scores.parquet", columns=["sample", "layer", "feature_id", "percentile", "call"])
X = X[X.layer.str.startswith("taxonomy") & X["sample"].isin(set(S.sample_key))]
X["feature"] = X.layer.str.replace("taxonomy_", "") + ":" + X.feature_id.astype(str); X["v"] = X.percentile.fillna(0.0)
PCT = X.pivot_table(index="sample", columns="feature", values="v", aggfunc="first").reindex(S.sample_key).fillna(0.0)
M = pd.read_parquet(f"{P}/work/s12/B/oos_summaries.parquet"); M = M[M.layer.str.startswith("taxonomy")]
SUM = pd.concat([M.pivot_table(index="sample", columns="layer", values="frac_outside_raw").add_prefix("frac_outside:"),
                 M.pivot_table(index="sample", columns="layer", values="n_missing").add_prefix("n_missing:")], axis=1).reindex(S.sample_key).fillna(0.0)
sp = pd.read_parquet(f"{P}/work/s11/metaphlan/species_long.parquet"); sp = sp[sp.sample_key.isin(set(S.sample_key)) & (sp.rel_abundance > 1e-5)]
PRES = pd.crosstab(sp.sample_key, sp.species).clip(upper=1).reindex(S.sample_key).fillna(0).add_prefix("present:")   # species kept per fold (>= 0.5 % of TRAINING samples)
SETS = {"presence": [PRES], "pct": [PCT], "pct+summary": [PCT, SUM], "pct+summary+presence": [PCT, SUM, PRES]}
if LAY == "rawcmp":   # the same model class on raw CLR instead of percentiles, so that the transformation is the only difference
    XR = pd.read_parquet(f"{P}/work/s12/B/oos_scores.parquet", columns=["sample", "layer", "feature_id", "value"])
    XR = XR[XR.layer.str.startswith("taxonomy") & XR["sample"].isin(set(S.sample_key))]
    XR["feature"] = XR.layer.str.replace("taxonomy_", "") + ":" + XR.feature_id.astype(str)
    RAW = XR.pivot_table(index="sample", columns="feature", values="value", aggfunc="first").reindex(S.sample_key).fillna(0.0).astype(np.float32)
    del XR
    SETS = {"pct": [PCT], "raw": [RAW], "pct+summary+presence": [PCT, SUM, PRES], "raw+summary+presence": [RAW, SUM, PRES]}
    print(f"raw CLR features: {RAW.shape[1]}", flush=True)
if LAY == "taxfun":
    F = pd.read_parquet(f"{P}/work/s12/B/oos_scores.parquet", columns=["sample", "layer", "feature_id", "percentile"],
                        filters=[("layer", "in", ["ko_humann", "pathway_humann"])])
    F = F[F["sample"].isin(set(S.sample_key))]
    F["feature"] = F.layer.str.replace("_humann", "") + ":" + F.feature_id.astype(str)
    PCTF = F.pivot_table(index="sample", columns="feature", values="percentile", aggfunc="first").reindex(S.sample_key).fillna(0.0).astype(np.float32)
    del F
    MF = pd.read_parquet(f"{P}/work/s12/B/oos_summaries.parquet"); MF = MF[MF.layer.str.endswith("_humann")]
    SUMF = pd.concat([MF.pivot_table(index="sample", columns="layer", values="frac_outside_raw").add_prefix("frac_outside:"),
                      MF.pivot_table(index="sample", columns="layer", values="n_missing").add_prefix("n_missing:")], axis=1).reindex(S.sample_key).fillna(0.0)
    SETS = {"pct": [PCT], "pct+function": [PCT, PCTF],
            "pct+summary+presence": [PCT, SUM, PRES], "pct+summary+presence+function": [PCT, SUM, PRES, PCTF, SUMF]}
    print(f"function features: pct {PCTF.shape[1]}, summary {SUMF.shape[1]}", flush=True)
y = S.y.to_numpy(); st = S.study.to_numpy()
# studies that share samples (same run or identical profile; work/s11/pipeB/cross_study_duplicates.tsv) form one group in the inner CV,
# and a training sample that duplicates a test-study sample is dropped from that fold's training set
DUP = pd.read_csv(f"{P}/work/s11/pipeB/cross_study_duplicates.tsv", sep="\t"); grp = pd.Series(st, index=S.sample_key)
for a_, b_ in DUP[["study", "duplicate_study"]].drop_duplicates().itertuples(index=False):
    grp[grp == b_] = a_ if a_ < b_ else grp[grp == b_]; grp[grp == a_] = min(a_, b_)
grp = grp.to_numpy(); dup_of_study = DUP.groupby("duplicate_study").sample_key.apply(set).to_dict()
print(f"samples {len(S)} (cases {y.sum()}), studies {len(set(st))}; features pct {PCT.shape[1]}, summary {SUM.shape[1]}, presence {PRES.shape[1]}; {time.time()-t0:.0f}s", flush=True)


def model(C):
    return make_pipeline(StandardScaler(), LogisticRegression(penalty="l1", C=C, solver="liblinear", class_weight="balanced", max_iter=2000))


def fold(name, Xm, cols, test_study):
    te = st == test_study; tr = ~te & ~S.sample_key.isin(dup_of_study.get(test_study, set())).to_numpy(); g = grp[tr]
    keep = np.array([not c.startswith("present:") or Xm[tr, k].mean() >= 0.005 for k, c in enumerate(cols)]); Xm = Xm[:, keep]; best = (None, -1)
    for C in CS:   # inner grouped CV on the training studies only
        sc = []
        for a, b in GroupKFold(5).split(Xm[tr], y[tr], g):
            if len(set(y[tr][b])) < 2:
                continue
            m = model(C).fit(Xm[tr][a], y[tr][a]); sc.append(roc_auc_score(y[tr][b], m.decision_function(Xm[tr][b])))
        if np.mean(sc) > best[1]:
            best = (C, np.mean(sc))
    m = model(best[0]).fit(Xm[tr], y[tr])
    return {"feature_set": name, "study": test_study, "C": best[0], "inner_auroc": best[1], "auroc": roc_auc_score(y[te], m.decision_function(Xm[te])),
            "n_nonzero": int((m[-1].coef_ != 0).sum())}


rows = []
for name, parts in SETS.items():
    Xd = pd.concat(parts, axis=1); Xm = Xd.to_numpy(dtype=np.float32)
    rows += Parallel(n_jobs=NJ)(delayed(fold)(name, Xm, list(Xd.columns), s) for s in EVAL)
    print(name, f"{time.time()-t0:.0f}s", flush=True)
R = pd.DataFrame(rows)
H = pd.read_csv(f"{P}/results/s15/health_indices_B.tsv", sep="\t"); seen = H.drop_duplicates("study").set_index("study")["in_gmwi2_training"]
G = H[H["index"] == "GMWI2"].set_index("study")["auroc_direct"]; base = H.drop_duplicates("study").set_index("study")
R["gmwi2"] = R.study.map(G); R["percentiles_rf"] = R.study.map(base["percentiles"]); R["in_gmwi2_training"] = R.study.map(seen)
R.to_csv(f"{OUT}/health_score_B{SUF}.tsv", sep="\t", index=False)
summ = []
for name, g in R.groupby("feature_set", sort=False):
    for lab, h in [("all 14", g), ("4 not seen by GMWI2", g[~g.in_gmwi2_training.astype(bool)])]:
        d = h.auroc - h.gmwi2
        summ.append({"feature_set": name, "studies": lab, "score_auroc": h.auroc.mean(), "gmwi2_auroc": h.gmwi2.mean(), "score_better": int((d > 0).sum()), "n": len(h),
                     "wilcoxon_p_vs_gmwi2": wilcoxon(d).pvalue if len(h) >= 6 else np.nan, "median_nonzero_features": h.n_nonzero.median()})
Sm = pd.DataFrame(summ)
# final model on all samples, best feature set, C by grouped CV, for the coefficient table
best = "pct+summary+presence" if LAY in ("tax", "rawcmp") else "pct+summary+presence+function"   # primary model, fixed before results (configs/decisions.yaml health_score_prespec)
Xb = pd.concat(SETS[best], axis=1); Xb = Xb[[c for c in Xb.columns if not c.startswith("present:") or Xb[c].mean() >= 0.005]]; cv = {C: np.mean([roc_auc_score(y[b], model(C).fit(Xb.values[a], y[a]).decision_function(Xb.values[b])) for a, b in GroupKFold(5).split(Xb, y, grp)]) for C in CS}
mf = model(max(cv, key=cv.get)).fit(Xb.values, y); coef = pd.Series(mf[-1].coef_[0], index=Xb.columns); coef = coef[coef != 0].sort_values()
coef.rename("coefficient").to_csv(f"{OUT}/health_score_B{SUF}_coefficients.tsv", sep="\t")
L = [f"# Supervised health-vs-disease score on refmb output (read pipeline), leave one study out; training studies: {TRAIN}; layers: {LAY}\n",
     f"{len(S):,} samples ({y.sum():,} cases) from {len(set(st))} studies; evaluated on the 14 case/control studies, each scored by a model trained without it. GMWI2 is scored as published (its training set includes 10 of the 14 studies).\n",
     Sm.round(3).to_markdown(index=False), f"\nFinal model ({best}, all samples): {len(coef)} non-zero features (results/s16/health_score_B_coefficients.tsv).\n",
     "## Per study\n", R.round(3).to_markdown(index=False), f"\nWall {time.time()-t0:.0f}s."]
open(f"{OUT}/health_score_B{SUF}.md", "w").write("\n".join(L)); print("\n".join(L[:4]))
