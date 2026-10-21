#!/usr/bin/env python3
"""What does ageing look like read against the reference range?

The disease sections ask which features of a patient sit outside a range built on reference adults. This asks the
same of a variable that is not a disease at all and that every clinical reference interval is stratified by: age.

The design is the disease atlas's, so the numbers are comparable to it. Within one study, the shift of a feature is
the median percentile of its reference adults aged 65 and over minus that of its adults aged 18 to 39, both scored
against a baseline rebuilt without that study. Only the six studies carrying at least twenty of each contribute, so
age is contrasted within a study and never across them: the confound that made the population analysis equivocal is
removed by construction here rather than argued about.

A shift is significant at ten points and q <= 0.05 (Mann-Whitney, Benjamini-Hochberg) and consistent when it holds in
the same direction in at least three of the six studies.

Writes results/s19/age_signature.{tsv,md} and the figure data.
"""
import os
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

P = os.environ.get("CHECKGM_PROJECT") or "/home/classes/bios/270/khoa/meta2/project"
OUT, FIG = f"{P}/results/s19", os.environ.get("CHECKGM_FIGDATA") or f"{P}/figures/v2/data"
LAYERS = {"taxonomy_family": "Family", "taxonomy_genus": "Genus", "taxonomy_species": "Species"}
YOUNG, OLD, MIN_ARM, MIN_STUDIES, MIN_SHIFT = (18, 39), (65, 200), 20, 3, 10.0


def cohort() -> pd.DataFrame:
    S = pd.read_csv(f"{P}/work/s11/pipeB/samples.tsv", sep="\t", low_memory=False)
    S = S[S.role.isin(["reference_pool", "heldout_healthy"])].copy()
    S["age"] = pd.to_numeric(S.age, errors="coerce")
    S = S.dropna(subset=["age"])
    S["arm"] = np.where(S.age.between(*YOUNG), "young", np.where(S.age.between(*OLD), "old", None))
    S = S[S.arm.notna()]
    ok = S.groupby(["study", "arm"]).size().unstack(fill_value=0)
    ok = ok.index[(ok.get("young", 0) >= MIN_ARM) & (ok.get("old", 0) >= MIN_ARM)]
    return S[S.study.isin(set(ok))]


def main():
    os.makedirs(OUT, exist_ok=True)
    S = cohort()
    tax = pd.read_parquet(f"{P}/resources/backbone/ncbi_taxonomy.parquet", columns=["taxid", "name"])
    NAME = dict(zip(tax.taxid.astype(str), tax.name))
    arm = S.set_index("sample_key").arm
    rows = []
    for layer, lab in LAYERS.items():
        X = pd.read_parquet(f"{P}/work/s12/B/oos_scores.parquet", columns=["sample", "feature_id", "percentile"],
                            filters=[("layer", "==", layer)])
        X = X[X["sample"].isin(set(S.sample_key)) & X.percentile.notna()]
        X["arm"] = X["sample"].map(arm)
        X["study"] = X["sample"].map(S.set_index("sample_key").study)
        for (study, fid), g in X.groupby(["study", "feature_id"]):
            y, o = g.percentile[g.arm == "young"], g.percentile[g.arm == "old"]
            if len(y) < MIN_ARM or len(o) < MIN_ARM:
                continue
            rows.append({"layer": lab, "feature_id": str(fid), "study": study,
                         "n_young": len(y), "n_old": len(o),
                         "shift": float(o.median() - y.median()),
                         "p": float(mannwhitneyu(o, y, alternative="two-sided").pvalue)})
    F = pd.DataFrame(rows)
    F["q"] = np.nan
    for (l, s), idx in F.groupby(["layer", "study"]).groups.items():
        F.loc[idx, "q"] = multipletests(F.loc[idx, "p"], method="fdr_bh")[1]
    F["sig"] = (F.q <= 0.05) & (F["shift"].abs() >= MIN_SHIFT)
    F["name"] = F.feature_id.map(NAME)
    F.to_csv(f"{OUT}/age_signature_per_study.tsv", sep="\t", index=False)

    agg = F.groupby(["layer", "feature_id"]).agg(
        studies=("study", "nunique"), n_up=("shift", lambda s: int((s >= MIN_SHIFT).sum())),
        n_down=("shift", lambda s: int((s <= -MIN_SHIFT).sum())),
        n_sig=("sig", "sum"), median_shift=("shift", "median")).reset_index()
    agg["name"] = agg.feature_id.map(NAME)
    agg["consistent"] = (agg[["n_up", "n_down"]].max(axis=1) >= MIN_STUDIES) & \
                        (agg[["n_up", "n_down"]].min(axis=1) == 0) & (agg.n_sig >= 2)
    agg = agg[agg["name"].notna()].sort_values("median_shift")
    agg.to_csv(f"{OUT}/age_signature.tsv", sep="\t", index=False)

    C = agg[agg.consistent]
    top = pd.concat([C.nsmallest(8, "median_shift"), C.nlargest(8, "median_shift")]).drop_duplicates("feature_id")
    F[F.feature_id.isin(set(top.feature_id))].merge(top[["feature_id", "median_shift", "consistent"]],
                                                    on="feature_id").to_csv(f"{FIG}/f6_age.csv", index=False)
    L = [f"# Ageing read against the reference range (read pipeline)\n",
         f"{len(S):,} reference adults from {S.study.nunique()} studies carrying at least {MIN_ARM} aged 18-39 and "
         f"{MIN_ARM} aged 65+: {(S.arm=='young').sum():,} young, {(S.arm=='old').sum():,} old. Age is contrasted "
         f"within a study, never across studies. Shift is the old median percentile minus the young, in points.\n",
         f"Consistent (same direction and at least {MIN_SHIFT:.0f} points in {MIN_STUDIES} of "
         f"{S.study.nunique()} studies, significant in at least two): **{int(C.shape[0])} features** "
         f"({int((C.median_shift<0).sum())} down with age, {int((C.median_shift>0).sum())} up).\n",
         "## The largest consistent shifts\n",
         top.sort_values("median_shift")[["layer", "name", "studies", "n_up", "n_down", "n_sig", "median_shift"]]
            .round(1).to_markdown(index=False)]
    open(f"{OUT}/age_signature.md", "w").write("\n".join(L))
    print("\n".join(L))


if __name__ == "__main__":
    main()
