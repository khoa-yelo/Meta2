#!/usr/bin/env python3
"""Assembly pipeline: alpha diversity scored directly, as the read pipeline scores it (s15: AUROC of -Shannon per study,
lower diversity meaning disease, no model), so that Fig. 2 compares the two pipelines on one definition. The assembly
output has no species layer, so Shannon is taken on family proportions, the same quantity s8's alpha_diversity feature
set uses. Samples and studies are those of the s8 leave-one-study-out evaluation. Writes results/s8/alpha_direct_A.tsv."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pandas as pd
from sklearn.metrics import roc_auc_score
import s8_evaluate as e

folds, meta, scores, summ = e.load()
a = e.feature_sets(meta, scores, summ)["alpha_diversity"].reindex(meta.index)
y, st = meta["y"], meta["study_bioproject"]
r = pd.DataFrame([{"study": s, "auroc_direct": roc_auc_score(y[st == s], -a.shannon_family[st == s])} for s in sorted(st.unique())])
r.to_csv(f"{e.P}/results/s8/alpha_direct_A.tsv", sep="\t", index=False)
print(r.round(3).to_string(index=False), "\nmean", round(r.auroc_direct.mean(), 3))
