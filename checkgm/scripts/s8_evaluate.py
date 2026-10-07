#!/usr/bin/env python3
"""S8 — case/control evaluation, leave-one-study-out (plan §7 S8 step 3, plus ablations).

Folds: the case/control studies in work/s8_fold_studies.txt. Every study's samples were scored (S7) against a
reference built WITHOUT that study (work/s8/scores/<study>.parquet), so no sample ever sees a reference containing
its own study. Feature sets, same classifier (random forest, 500 trees, seeded), same folds:
  reference_relative  — per-feature percentile (undetected -> 0; expected_but_missing -> 0) for family, genus,
                        eggNOG KO and module layers, plus per-layer deviation counts
  raw_clr             — family + genus CLR, and log10(copies-per-genome + 1e-3) for eggNOG KO, module completeness
  alpha_diversity     — Shannon, richness (family proportions), Shannon on KO copies
  health_index        — genus-level approximation of GMHI (Gupta 2020) from family/genus proportions (stated as
                        an approximation: GMHI is defined on 50 MetaPhlAn species)
Task: case (Diseased) vs control (Healthy/Control/cMD3 control), in-stratum adults only. Train on all other folds
(cases + controls), test on the held-out study; AUROC per study, pooled over studies (macro), and per condition
group where >= 2 studies carry it. Within-study AUROC (5-fold CV inside each study) is reported as the ceiling.
"""
from __future__ import annotations

import json
import os
import re
import time

import numpy as np
import pandas as pd
import yaml
from scipy.stats import entropy
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold

P = "/home/classes/bios/270/khoa/meta2/project"
SEED = 20261020
CFG5 = yaml.safe_load(open(f"{P}/configs/s5.yaml"))
LAYERS_RR = ["taxonomy_family", "taxonomy_genus", "ko_eggnog", "module"]
MIN_PREV = 0.10   # features kept for modelling: detected in >= 10% of evaluation samples
GMHI_GOOD = {"Alistipes", "Bifidobacterium", "Sutterella"}
GMHI_BAD = {"Anaerotruncus", "Atopobium", "Blautia", "Clostridium", "Eggerthella", "Flavonifractor", "Fusobacterium", "Gemella", "Granulicatella",
            "Holdemania", "Klebsiella", "Lactobacillus", "Peptostreptococcus", "Mediterraneibacter", "Solobacterium", "Streptococcus", "Subdoligranulum", "Veillonella"}


SCORES_DIR = f"{P}/work/s8/scores"; TAG = "" if os.environ.get("MODEL", "rf") == "rf" else "_" + os.environ["MODEL"]


def load():
    folds = [l.strip() for l in open(f"{P}/work/s8_fold_studies.txt") if l.strip()]
    inv = pd.read_parquet(f"{P}/work/s0/inventory.parquet")
    st = yaml.safe_load(open(f"{P}/configs/inclusion.yaml"))["stratum"]
    meta = inv[inv["study_bioproject"].isin(folds) & inv["is_primary_analysis"]].copy()
    age = pd.to_numeric(meta["cur_Age_years"], errors="coerce"); age = age.where(age <= 110)
    oos = (age.notna() & (age < 18)) | meta["cur_Age_category"].fillna("").str.contains("|".join(st["age_category_disallow"]))
    meta = meta[~oos]
    meta["y"] = np.where(meta["cur_Health_status_group"].eq("Diseased"), 1,
                         np.where(meta["cur_Health_status_group"].isin(["Healthy", "Control"]) | meta["cmd3_study_condition"].eq("control"), 0, -1))
    meta = meta[meta["y"] >= 0]
    scores = pd.concat([pd.read_parquet(f"{SCORES_DIR}/{s}.parquet") for s in folds if os.path.exists(f"{SCORES_DIR}/{s}.parquet")], ignore_index=True)
    summ = pd.concat([pd.read_parquet(f"{SCORES_DIR}/{s}_summaries.parquet") for s in folds if os.path.exists(f"{SCORES_DIR}/{s}_summaries.parquet")], ignore_index=True)
    meta = meta[meta["analysis_id"].isin(scores["analysis_id"])]
    return folds, meta.set_index("analysis_id"), scores, summ


def feature_sets(meta, scores, summ):
    ids = meta.index
    X = {}
    # reference-relative
    blocks = []
    for layer in LAYERS_RR:
        t = scores[(scores["layer"] == layer) & scores["analysis_id"].isin(ids)]
        pct = t.pivot_table(index="analysis_id", columns="feature_id", values="percentile", aggfunc="first")
        ebm = t[t["call"] == "expected_but_missing"].pivot_table(index="analysis_id", columns="feature_id", values="percentile", aggfunc="size")
        pct = pct.reindex(ids)
        pct = pct.loc[:, pct.notna().mean() >= MIN_PREV]
        pct = pct.fillna(0.0)                     # undetected / not assessable / expected-but-missing -> 0th percentile
        pct.columns = [f"{layer}|{c}|pct" for c in pct.columns]
        blocks.append(pct)
    dev = summ[summ["analysis_id"].isin(ids)].pivot_table(index="analysis_id", columns="layer", values=["n_low_raw", "n_high_raw", "n_missing", "frac_outside_raw"]).reindex(ids).fillna(0)
    dev.columns = [f"dev|{a}|{b}" for a, b in dev.columns]
    X["reference_relative"] = pd.concat(blocks + [dev], axis=1)
    X["reference_relative_taxonomy_only"] = pd.concat([b for b in blocks if b.columns[0].startswith("taxonomy")] + [dev.filter(like="taxonomy")], axis=1)
    X["reference_relative_genes_only"] = pd.concat([b for b in blocks if not b.columns[0].startswith("taxonomy")] + [dev.loc[:, ~dev.columns.str.contains("taxonomy")]], axis=1)
    # raw
    fam = pd.read_parquet(f"{P}/work/s2/taxonomy_family_assembly.parquet").query("analysis_id in @ids")
    gen = pd.read_parquet(f"{P}/work/s2/taxonomy_genus_assembly.parquet").query("analysis_id in @ids")
    ko = pd.read_parquet(f"{P}/work/s2/ko_eggnog_assembly.parquet", columns=["analysis_id", "feature_id", "copies_per_genome"]).query("analysis_id in @ids")
    mod = pd.read_parquet(f"{P}/work/s2/module_assembly.parquet", columns=["analysis_id", "feature_id", "completeness"]).query("analysis_id in @ids")
    def piv(t, v, prefix, fill):
        m = t.pivot_table(index="analysis_id", columns="feature_id", values=v, aggfunc="first").reindex(ids)
        m = m.loc[:, m.notna().mean() >= MIN_PREV].fillna(fill); m.columns = [f"{prefix}|{c}" for c in m.columns]; return m
    fam_clr = piv(fam, "clr", "fam_clr", fam["clr"].min()); gen_clr = piv(gen, "clr", "gen_clr", gen["clr"].min())
    ko_log = piv(ko.assign(v=np.log10(ko["copies_per_genome"] + 1e-3)), "v", "ko_log", np.log10(1e-3))
    mod_c = piv(mod, "completeness", "module", 0.0)
    X["raw_clr"] = pd.concat([fam_clr, gen_clr, ko_log, mod_c], axis=1)
    # alpha diversity
    fp = fam.pivot_table(index="analysis_id", columns="feature_id", values="proportion_mapped", aggfunc="first").reindex(ids).fillna(0)
    kp = ko.pivot_table(index="analysis_id", columns="feature_id", values="copies_per_genome", aggfunc="first").reindex(ids).fillna(0)
    alpha = pd.DataFrame({"shannon_family": entropy(fp.to_numpy().T), "richness_family": (fp > 0.001).sum(axis=1),
                          "shannon_ko": entropy((kp.div(kp.sum(axis=1), axis=0)).to_numpy().T), "richness_ko": (kp > 0.05).sum(axis=1)}, index=ids)
    X["alpha_diversity"] = alpha
    # genus-level GMHI approximation
    tax = pd.read_parquet(f"{P}/resources/backbone/ncbi_taxonomy.parquet", columns=["taxid", "name"]); name = dict(zip(tax["taxid"].astype(str), tax["name"]))
    gp = gen.pivot_table(index="analysis_id", columns="feature_id", values="proportion_mapped", aggfunc="first").reindex(ids).fillna(0)
    gp.columns = [name.get(c, c) for c in gp.columns]
    good = [c for c in gp.columns if c in GMHI_GOOD]; bad = [c for c in gp.columns if c in GMHI_BAD]
    theta = 1e-5
    psi_g = (gp[good] > theta).sum(axis=1) / max(len(good), 1) * np.log((gp[good].clip(lower=theta)).mean(axis=1) / theta + 1e-9)
    psi_b = (gp[bad] > theta).sum(axis=1) / max(len(bad), 1) * np.log((gp[bad].clip(lower=theta)).mean(axis=1) / theta + 1e-9)
    X["health_index"] = pd.DataFrame({"gmhi_genus_approx": np.log10((psi_g + 1e-5) / (psi_b + 1e-5)), "n_good_genera": (gp[good] > theta).sum(axis=1), "n_bad_genera": (gp[bad] > theta).sum(axis=1)}, index=ids)
    return X


# MODEL=l1 (the paper's classifier) or MODEL=rf (the original run, kept so the earlier tables stay reproducible). Every
# representation is scored by the same model class, so a difference between feature sets cannot be a difference of
# classifier; the L1 configuration is identical to s16's, including its C grid and its grouped inner selection.
MODEL = os.environ.get("MODEL", "rf")
CS = [0.003, 0.01, 0.03, 0.1]


def rf():
    return RandomForestClassifier(n_estimators=500, random_state=SEED, n_jobs=8, class_weight="balanced", min_samples_leaf=2)


def l1(C):
    return make_pipeline(StandardScaler(), LogisticRegression(penalty="l1", C=C, solver="liblinear",
                                                              class_weight="balanced", max_iter=2000))


def fit_predict(Mtr, ytr, gtr, Mte):
    """Train on one fold and return P(case) for its test rows."""
    if MODEL == "rf":
        return rf().fit(Mtr, ytr).predict_proba(Mte)[:, 1]
    k = min(5, len(set(gtr)))
    if k < 2:
        best = 0.03
    else:
        cv = {}
        for C in CS:
            sc = []
            for a, b in GroupKFold(k).split(Mtr, ytr, gtr):
                if len(set(ytr[a])) < 2 or len(set(ytr[b])) < 2:
                    continue
                sc.append(roc_auc_score(ytr[b], l1(C).fit(Mtr[a], ytr[a]).decision_function(Mtr[b])))
            if sc:
                cv[C] = np.mean(sc)
        best = max(cv, key=cv.get) if cv else 0.03
    return l1(best).fit(Mtr, ytr).decision_function(Mte)


def main():
    t0 = time.time()
    folds, meta, scores, summ = load()
    X = feature_sets(meta, scores, summ)
    y = meta["y"].to_numpy(); study = meta["study_bioproject"].to_numpy(); cond = meta["cur_Condition_group"].fillna("").to_numpy()
    rows = []; preds = []
    for fs, M in X.items():
        M = M.reindex(meta.index).fillna(0).to_numpy(dtype=float)
        pooled_scores = np.full(len(y), np.nan)
        for s in folds:
            te = study == s; tr = ~te
            if te.sum() == 0 or len(set(y[te])) < 2 or len(set(y[tr])) < 2:
                continue
            p = fit_predict(M[tr], y[tr], study[tr], M[te]); pooled_scores[te] = p
            rng = np.random.default_rng(SEED); yt = y[te]; boots = []
            for _ in range(300):
                ix = rng.integers(0, len(yt), len(yt))
                if len(set(yt[ix])) == 2: boots.append(roc_auc_score(yt[ix], p[ix]))
            rows.append(dict(feature_set=fs, study=s, n_test=int(te.sum()), n_cases=int(y[te].sum()), condition=pd.Series(cond[te][y[te] == 1]).mode().iat[0] if (y[te] == 1).any() else "",
                             auroc_loso=roc_auc_score(y[te], p), ci_low=np.percentile(boots, 2.5), ci_high=np.percentile(boots, 97.5), n_features=M.shape[1]))
            # within-study ceiling
            if te.sum() >= 20 and min((y[te] == 1).sum(), (y[te] == 0).sum()) >= 5:
                skf = StratifiedKFold(5, shuffle=True, random_state=SEED); pw = np.zeros(te.sum())
                Mt, yt = M[te], y[te]
                for a, b in skf.split(Mt, yt):
                    pw[b] = rf().fit(Mt[a], yt[a]).predict_proba(Mt[b])[:, 1]
                rows[-1]["auroc_within_study_cv"] = roc_auc_score(yt, pw)
        ok = ~np.isnan(pooled_scores)
        rows.append(dict(feature_set=fs, study="POOLED", n_test=int(ok.sum()), n_cases=int(y[ok].sum()), condition="all", auroc_loso=roc_auc_score(y[ok], pooled_scores[ok]), n_features=M.shape[1]))
        preds.append(pd.DataFrame({"analysis_id": meta.index[ok], "study": study[ok], "y": y[ok], "score": pooled_scores[ok], "feature_set": fs}))
        print(fs, f"{time.time()-t0:.0f}s", flush=True)
    R = pd.DataFrame(rows)
    os.makedirs(f"{P}/results/s8", exist_ok=True)
    pd.concat(preds, ignore_index=True).to_parquet(f"{P}/results/s8/loso_predictions{TAG}.parquet", index=False)
    R.to_csv(f"{P}/results/s8/loso_auroc{TAG}.tsv", sep="\t", index=False)
    per_study = R[R["study"] != "POOLED"].pivot(index=["study", "condition", "n_test", "n_cases"], columns="feature_set", values="auroc_loso").round(3)
    macro = R[R["study"] != "POOLED"].groupby("feature_set")["auroc_loso"].agg(["mean", "median"]).round(3)
    pooled = R[R["study"] == "POOLED"].set_index("feature_set")["auroc_loso"].round(3)
    within = R[R["study"] != "POOLED"].groupby("feature_set")["auroc_within_study_cv"].median().round(3)
    # per condition with >= 2 studies
    cg = R[R["study"] != "POOLED"].groupby(["condition", "feature_set"]).agg(n_studies=("study", "nunique"), auroc_mean=("auroc_loso", "mean")).reset_index()
    cg = cg[cg["n_studies"] >= 2].pivot(index=["condition", "n_studies"], columns="feature_set", values="auroc_mean").round(3)
    from scipy.stats import wilcoxon
    ps = R[R["study"] != "POOLED"].pivot(index="study", columns="feature_set", values="auroc_loso")
    paired = []
    for other in [c for c in ps.columns if c != "reference_relative"]:
        d = ps["reference_relative"] - ps[other]
        paired.append(dict(comparison=f"reference_relative − {other}", mean_diff=round(d.mean(), 3), wins=int((d > 0).sum()), losses=int((d < 0).sum()),
                           wilcoxon_p=round(wilcoxon(d, alternative="greater").pvalue, 3) if d.abs().sum() > 0 else np.nan))
    paired = pd.DataFrame(paired)
    L = [f"# S8 — leave-one-study-out case/control AUROC ({len(folds)} studies, {len(meta):,} in-stratum samples, {int(y.sum()):,} cases)\n",
         "Same random forest (500 trees, balanced class weights, seed 20261020), same folds, for every feature set. Each study's features were computed against a reference built without that study.\n",
         "## Per held-out study\n", per_study.to_markdown(), "\n\n## Summary\n",
         pd.DataFrame({"macro_mean": macro["mean"], "macro_median": macro["median"], "pooled_auroc": pooled, "within_study_cv_median (ceiling)": within}).to_markdown(),
         f"\n\n## Paired comparison across the {len(folds)} held-out studies (one-sided Wilcoxon, reference_relative greater)\n", paired.to_markdown(index=False),
         "\n\n## Per-study 95% bootstrap CI (300 resamples of the test set)\n",
         R[R["study"] != "POOLED"].assign(ci=lambda d: d.apply(lambda r: f"{r.auroc_loso:.3f} [{r.ci_low:.2f}, {r.ci_high:.2f}]", axis=1)).pivot(index="study", columns="feature_set", values="ci").to_markdown(),
         "\n\n## Per condition group (≥ 2 studies)\n", cg.to_markdown() if len(cg) else "none", "\n"]
    open(f"{P}/results/s8/report_loso{TAG}.md", "w").write("\n".join(L)); print("\n".join(L))
    json.dump({"created": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "folds": folds, "n_samples": int(len(meta)), "n_cases": int(y.sum()), "seed": SEED,
               "feature_counts": {k: int(v.shape[1]) for k, v in X.items()}, "wall_s": round(time.time() - t0)}, open(f"{P}/results/s8/provenance_loso{TAG}.json", "w"), indent=1)


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(); ap.add_argument("--scores-dir", default=SCORES_DIR); ap.add_argument("--tag", default=TAG)   # default comes from MODEL; an explicit --tag still wins
    a = ap.parse_args(); SCORES_DIR = a.scores_dir; TAG = a.tag
    main()
