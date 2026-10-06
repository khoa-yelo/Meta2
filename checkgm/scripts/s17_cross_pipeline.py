#!/usr/bin/env python3
"""Do the two pipelines assess the same sample the same way?

Every other comparison in this project asks whether a pipeline reproduces its own source (s11 fidelity). This one asks
the question a reader of a two-pipeline framework actually has: take one biological sample, measure it both ways, score
each measurement against its own baseline, and see whether the two reports agree.

The bridge is work/s11/pipeB/samples.tsv, which carries the MGnify analysis_id beside the cMD3 sample_key, so a sample
present in both can be matched. Both sides are scored out of sample (s12 scores pool studies against a baseline rebuilt
without them), so neither report is inflated by having seen the sample. Four of the bridge samples need the stronger
form of that guarantee, and get it: NielsenHB_2014_MH0050 and _MH0054 also sit in the pool as LeChatelierE_2013
specimens (the two studies are two cMD3 deposits of one MetaHIT cohort), and GuptaA_2019_GupDM_HCQ and _HDF also sit in
it as DhakanDB_2019 specimens, so s12 scores all four against a baseline rebuilt without the duplicate's study as well
as without their own (work/s11/pipeB/cross_study_duplicates.tsv, s11_cross_study_duplicates.py).

Two things are measured, and they answer differently:
  ranking  — Spearman rho between the two sets of percentiles, over the features both pipelines assess
  calls    — whether the same feature is placed outside the healthy range by both

The second needs a chance baseline. The two pipelines flag very different numbers of features, so a per-sample
agreement rate is misleading: when one side flags one feature out of twenty-five, most samples score zero overlap by
construction. The pooled counts and the overlap expected under independence are reported instead.

Writes results/s11/cross_pipeline_agreement.tsv (per sample x layer) and prints the summary.
"""
import os
import numpy as np, pandas as pd
from scipy.stats import spearmanr

P = os.environ.get("CHECKGM_PROJECT") or os.environ.get("REFMB_PROJECT") or "/home/classes/bios/270/khoa/meta2/project"
OUT_SET = ["low", "high", "expected_but_missing"]
LAYERS = ["taxonomy_family", "taxonomy_genus"]


def bridge() -> pd.DataFrame:
    """Samples carrying both an MGnify analysis and a cMD3 profile."""
    S = pd.read_csv(f"{P}/work/s11/pipeB/samples.tsv", sep="\t", low_memory=False)
    S["analysis_id"] = S.analysis_id.astype(str)
    anchors = set(pd.read_csv(f"{P}/work/s11/anchor_panel.tsv", sep="\t").analysis_id)
    return S[S.analysis_id.isin(anchors)][["sample_key", "analysis_id", "study"]].reset_index(drop=True)


def main():
    b = bridge()
    print(f"{len(b)} samples measured by both pipelines, from {b.study.nunique()} studies: {b.study.value_counts().to_dict()}\n")
    rows, pooled = [], []
    for layer in LAYERS:
        A = pd.read_parquet(f"{P}/work/s12/A/oos_scores.parquet", columns=["sample", "feature_id", "percentile", "call"], filters=[("layer", "==", layer)])
        B = pd.read_parquet(f"{P}/work/s12/B/oos_scores.parquet", columns=["sample", "feature_id", "percentile", "call"], filters=[("layer", "==", layer)])
        A["feature_id"] = A.feature_id.astype(str); B["feature_id"] = B.feature_id.astype(str)
        n_a = n_b = n_both = n_shared = 0; expected = 0.0
        for r in b.itertuples():
            x = A[A["sample"] == r.analysis_id].set_index("feature_id")
            y = B[B["sample"] == r.sample_key].set_index("feature_id")
            common = x.index.intersection(y.index)
            if len(common) < 10:
                continue
            xx, yy = x.loc[common], y.loc[common]
            ok = xx.percentile.notna() & yy.percentile.notna()   # an expected-but-missing feature has no percentile
            rho = spearmanr(xx.percentile[ok], yy.percentile[ok]).statistic if ok.sum() >= 10 else np.nan
            fa = set(common[xx.call.isin(OUT_SET)]); fb = set(common[yy.call.isin(OUT_SET)])
            n_a += len(fa); n_b += len(fb); n_both += len(fa & fb); n_shared += len(common)
            expected += len(fa) * len(fb) / len(common)
            rows.append(dict(layer=layer, sample=r.sample_key, study=r.study, n_assessed_assembly=len(x), n_assessed_read=len(y),
                             n_common=len(common), n_usable=int(ok.sum()), spearman_pct=rho,
                             call_agreement=float((xx.call == yy.call).mean()),
                             n_flagged_assembly=len(fa), n_flagged_read=len(fb), n_flagged_both=len(fa & fb)))
        pooled.append(dict(layer=layer, shared_assessments=n_shared, flagged_assembly=n_a, flagged_read=n_b,
                           flagged_both=n_both, expected_if_independent=round(expected, 1),
                           jaccard=round(n_both / max(n_a + n_b - n_both, 1), 3)))
    R = pd.DataFrame(rows); R.to_csv(f"{P}/results/s11/cross_pipeline_agreement.tsv", sep="\t", index=False)
    print("ranking (per sample, over the features both pipelines assess)")
    print(R.groupby("layer").agg(samples=("sample", "size"), assessed_assembly=("n_assessed_assembly", "median"),
                                 assessed_read=("n_assessed_read", "median"), common=("n_common", "median"),
                                 median_rho=("spearman_pct", "median"), call_agreement=("call_agreement", "median")).round(3).to_string())
    print("\ncalls (pooled over samples, with the overlap expected if the two flagged independently)")
    print(pd.DataFrame(pooled).to_string(index=False))


if __name__ == "__main__":
    main()
