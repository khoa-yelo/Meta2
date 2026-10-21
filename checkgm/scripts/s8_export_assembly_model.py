#!/usr/bin/env python3
"""Export the assembly-pipeline health score: the s8 reference_relative L1 model, refitted on all 11 case/control studies.

The model is exactly the one s8_evaluate.py (MODEL=l1) evaluates leave-one-study-out: the same samples (load()), the same
reference_relative feature matrix (feature_sets(): per-feature percentiles of family, genus, eggNOG KO and KEGG module,
undetected / expected-but-missing / not-assessable -> 0, features kept when assessable in >= 10 % of the evaluation
samples, plus n_low_raw / n_high_raw / n_missing / frac_outside_raw per layer), the same StandardScaler + L1 logistic
regression (l1()), and the same C selection (GroupKFold over studies on the training samples, C grid CS). Every
sample's percentiles were computed against the bundle built without its own study (work/s8/scores), so the fit sees
out-of-sample percentiles, as a new sample assessed against the shipped bundle does.

Checks before writing: (1) two leave-one-study-out folds are refitted here and must reproduce results/s8/loso_auroc_l1.tsv;
(2) the exported JSON is applied with checkgm.healthscore itself to the same score tables and must reproduce the sklearn
decision function on every sample, which pins the feature-name mapping and the absent -> 0 convention end to end.

Writes results/s8/health_score_A_model.json and copies it to the package as src/checkgm/health_score_model_assembly.json.
Run: MODEL=l1 python scripts/s8_export_assembly_model.py [--package-dir Meta2/checkgm/src/checkgm]
"""
import argparse
import json
import os
import shutil
import sys
import time

os.environ["MODEL"] = "l1"
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s8_evaluate as S8   # noqa: E402

P = S8.P
PKG = "/home/classes/bios/270/khoa/meta2/Meta2/checkgm/src/checkgm"
PREFIX = {"taxonomy_family": "family", "taxonomy_genus": "genus", "ko_eggnog": "ko_eggnog", "module": "module"}
STAT = {"frac_outside_raw": "frac_outside", "n_missing": "n_missing", "n_low_raw": "n_low", "n_high_raw": "n_high"}


def healthscore_name(col: str) -> str:
    """s8 column name -> checkgm.healthscore feature name ('<layer>|<id>|pct' or 'dev|<stat>|<layer>')."""
    a, b, c = col.split("|")
    if a == "dev":
        return f"{STAT[b]}:{c}"
    assert c == "pct", col
    return f"{PREFIX[a]}:{b}"


def select_C(M, y, g):
    """s8_evaluate.fit_predict's inner selection, on the given training rows."""
    k = min(5, len(set(g))); cv = {}
    for C in S8.CS:
        sc = []
        for a, b in GroupKFold(k).split(M, y, g):
            if len(set(y[a])) < 2 or len(set(y[b])) < 2:
                continue
            sc.append(roc_auc_score(y[b], S8.l1(C).fit(M[a], y[a]).decision_function(M[b])))
        if sc:
            cv[C] = float(np.mean(sc))
    return (max(cv, key=cv.get) if cv else 0.03), cv


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--package-dir", default=PKG); a = ap.parse_args()
    t0 = time.time()
    folds, meta, scores, summ = S8.load()
    X = S8.feature_sets(meta, scores, summ)["reference_relative"]
    Xd = X.reindex(meta.index).fillna(0); cols = list(Xd.columns); M = Xd.to_numpy(dtype=float)
    y = meta["y"].to_numpy(); study = meta["study_bioproject"].to_numpy()
    print(f"{len(y)} samples, {int(y.sum())} cases, {len(set(study))} studies, {M.shape[1]} features; {time.time()-t0:.0f}s", flush=True)

    # (1) the pipeline here is s8's: two folds reproduce the published leave-one-study-out AUROC
    ref = pd.read_csv(f"{P}/results/s8/loso_auroc_l1.tsv", sep="\t").query("feature_set == 'reference_relative'").set_index("study").auroc_loso
    for s in ("PRJEB26165", "PRJEB62821"):
        te = study == s; tr = ~te
        C, _ = select_C(M[tr], y[tr], study[tr])
        au = roc_auc_score(y[te], S8.l1(C).fit(M[tr], y[tr]).decision_function(M[te]))
        print(f"reproduce {s}: {au:.6f} vs s8 {ref[s]:.6f}", flush=True)
        assert abs(au - ref[s]) < 1e-9, (s, au, ref[s])

    # the final model: all 11 studies, C by the same grouped selection
    C, cv = select_C(M, y, study)
    mf = S8.l1(C).fit(M, y); sc_, lr = mf[0], mf[-1]
    coef = lr.coef_[0]; nz = np.flatnonzero(coef)
    feats = [{"name": healthscore_name(cols[i]), "coef": float(coef[i]), "mean": float(sc_.mean_[i]), "scale": float(sc_.scale_[i])} for i in nz]
    assert len({f["name"] for f in feats}) == len(feats)
    dec = mf.decision_function(M)
    print(f"C={C} (grouped CV {cv}); {len(nz)} non-zero of {len(cols)}; training AUROC {roc_auc_score(y, dec):.3f}", flush=True)

    model = {
        "model": "L1-penalised logistic regression on checkGM output, standardized features",
        "pipeline": "assembly",
        "feature_set": "reference_relative",
        "baseline": "gut-assembly-adult-global-v0.8-lenient (training percentiles: each study against its leave-one-study-out bundle, work/s8/scores)",
        "trained_on": {"samples": int(len(y)), "cases": int(y.sum()), "studies": int(len(set(study))), "training_studies": folds},
        "C": float(C), "C_grid": list(S8.CS), "C_selection": {"method": "GroupKFold(5) over studies, mean AUROC", "cv_auroc": {str(k): v for k, v in cv.items()}},
        "n_candidate_features": int(len(cols)), "min_prevalence": S8.MIN_PREV,
        "intercept": float(lr.intercept_[0]),
        "orientation": "higher score = more case-like; decision_function of the disease class",
        "validation": "leave-one-study-out AUROC over the 11 studies: results/s8/loso_auroc_l1.tsv (feature_set reference_relative)",
        "source": "scripts/s8_export_assembly_model.py (pipeline of scripts/s8_evaluate.py, MODEL=l1)",
        "created": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "features": feats,
    }

    # (2) checkgm.healthscore applies the export to the same tables and must agree with sklearn on every sample
    sys.path.insert(0, os.path.dirname(a.package_dir))
    from checkgm import healthscore as hs
    tmp = f"{P}/results/s8/health_score_A_model.json"
    json.dump(model, open(tmp, "w"), indent=1)
    m = hs.load_model(tmp); assert m["pipeline"] == "assembly"
    sm = summ[summ.analysis_id.isin(meta.index)]; sc = scores[scores.analysis_id.isin(meta.index)]
    assert hs.assessment_pipeline("s8", sc, sm) == "assembly"
    V, present = hs.matrix(sc, sm, m)
    _, s_hs = hs._contributions(V, m)
    s_hs = pd.Series(s_hs, index=V.index).reindex(meta.index.astype(str))
    assert s_hs.notna().all()
    err = float(np.max(np.abs(s_hs.to_numpy() - dec)))
    print(f"checkgm.healthscore vs sklearn: max |diff| {err:.2e} over {len(dec)} samples", flush=True)
    assert err < 1e-8, err
    npres = present.sum(axis=1)
    model["features_present_in_training"] = {"median": float(np.median(npres)), "p1": float(np.percentile(npres, 1)), "n_features": len(feats)}
    model["thin_features"] = int(np.floor(np.percentile(npres, 1)))
    json.dump(model, open(tmp, "w"), indent=1)
    dst = os.path.join(a.package_dir, "health_score_model_assembly.json"); shutil.copyfile(tmp, dst)
    print(f"features present per training sample: median {np.median(npres):.0f}, p1 {np.percentile(npres, 1):.1f} of {len(feats)}")
    print(f"wrote {tmp} and {dst}; {time.time()-t0:.0f}s")
    top = sorted(feats, key=lambda f: -abs(f["coef"]))[:15]
    print(pd.DataFrame(top)[["name", "coef"]].to_string(index=False))
    kinds = pd.Series([f["name"].split(":")[0] for f in feats]).value_counts(); print(kinds.to_string())


if __name__ == "__main__":
    main()
