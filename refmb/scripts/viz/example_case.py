"""The worked example shared by fig1_schematic, fig5_tree and make_summary: one C. difficile infection case
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


def none_or(n):
    return "none" if n == 0 else f"{n:,}"


if __name__ == "__main__":
    print(counts())
