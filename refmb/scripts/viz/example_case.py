"""The worked example shared by fig1_schematic, fig5_tree, fig14_sample_card and make_summary: one C. difficile infection case
(analysis MGYA00694698 of study PRJEB26165, assembly pipeline), scored against a baseline built without its study.
Every count the figures print (assessed, low, high, expected-but-missing) is computed here from the scores table over the
232 "leaf" families (reference prevalence >= 50 % in some band and percentiles available), never typed by hand."""
import os, sys
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import tokens as T
sys.path.insert(0, T.P); from refmb.paths import bundle_A

STUDY, CASE = "PRJEB26165", "MGYA00694698"
MIN_PREV = 0.5
BUNDLE = bundle_A("lenient", loso=STUDY, root=os.path.join(T.P, "work", "s8", "refs"))


def reference():
    return pd.read_parquet(os.path.join(BUNDLE, "features", "taxonomy_family.parquet")).set_index("feature_id")


def leaves(ref=None):
    ref = reference() if ref is None else ref
    pc = [c for c in ref.columns if c.startswith("prevalence_band_")]
    return ref.index[(ref[pc].max(axis=1) >= MIN_PREV) & ref["percentiles_available"]]


def names():
    tax = pd.read_parquet(os.path.join(T.P, "resources", "backbone", "ncbi_taxonomy.parquet"), columns=["taxid", "name"])
    return dict(zip(tax["taxid"].astype(str), tax["name"]))


def case_families(case=CASE):
    """One row per leaf family for the case: percentile, call ('within', 'low', 'high', 'expected_but_missing' or
    'unscored' when the family is neither detected nor expected), name."""
    sc = pd.read_parquet(os.path.join(T.P, "work", "s8", "scores", f"{STUDY}.parquet"),
                         filters=[("analysis_id", "==", case), ("layer", "==", "taxonomy_family")]).set_index("feature_id")
    lv = leaves(); nm = names()
    d = pd.DataFrame({"percentile": sc["percentile"].reindex(lv), "call": sc["call"].reindex(lv).fillna("unscored")})
    d["name"] = [nm.get(i, i) for i in d.index]
    return d


def counts(d=None):
    d = case_families() if d is None else d
    c = d["call"]
    n_ass = int(c.isin(["within", "low", "high"]).sum()); n_low = int((c == "low").sum()); n_high = int((c == "high").sum())
    return {"n_families": int(len(d)), "n_assessed": n_ass, "n_low": n_low, "n_high": n_high, "n_outside": n_low + n_high,
            "n_missing": int((c == "expected_but_missing").sum()), "frac_outside": (n_low + n_high) / n_ass if n_ass else float("nan")}


LAYER_ORDER = ["taxonomy_family", "taxonomy_genus", "ko_eggnog", "ko_kofam", "pfam", "module"]


def layer_counts(case=CASE):
    """Sample-level numbers of the case on every layer, from the same scores table as case_families() (work/s8/scores,
    baseline built without the study), over all assessed features of the layer (not only the leaf families): n_assessed,
    n_outside, frac_outside (share of assessed features outside the healthy range) and n_missing (expected but missing).
    This is the source of the gene-layer counts in results/s14/example_reports.tsv (e.g. 179 expected-but-missing KOs);
    work/s12/A/oos_summaries.parquet holds a second scoring of the same sample with slightly different counts (184)."""
    sc = pd.read_parquet(os.path.join(T.P, "work", "s8", "scores", f"{STUDY}.parquet"), filters=[("analysis_id", "==", case)], columns=["layer", "call"])
    rows = []
    for layer, g in sc.groupby("layer"):
        c = g["call"]; n_ass = int(c.isin(["within", "low", "high"]).sum()); n_out = int(c.isin(["low", "high"]).sum())
        rows.append({"layer": layer, "n_assessed": n_ass, "n_outside": n_out, "frac_outside_raw": n_out / n_ass if n_ass else float("nan"), "n_missing": int((c == "expected_but_missing").sum())})
    return pd.DataFrame(rows).set_index("layer").reindex([l for l in LAYER_ORDER if l in set(sc["layer"])])


def none_or(n):
    return "none" if n == 0 else f"{n:,}"


if __name__ == "__main__":
    print(counts()); print(layer_counts())
