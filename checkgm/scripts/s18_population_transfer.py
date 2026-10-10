#!/usr/bin/env python3
"""Does a reference range transfer to the populations it was not dominated by?

The paper tests this on one cohort, the Hadza, and the result is strong but single: 20.6% of families outside the
range against 4.4% elsewhere. One study cannot separate a population from the protocol that measured it, which the
paper says. This asks the same question of every constituent population of the read baseline, where 19 studies span
17 countries and five of them are non-industrialised or transitional.

Every sample here is a reference adult scored against a baseline rebuilt without its own study (s12 leave-one-study-
out), so a country's excess is a statement about a baseline that never saw it.

The statistic is the paper's: the mean per-sample fraction outside within a study, then the median over the studies of
a country. Per-sample medians must not be substituted -- with about 22 assessed families a 5% expectation is roughly
one feature, so more than half of the samples of several countries have exactly zero and the median is 0.0 for
reasons that have nothing to do with fit.

The decisive column is the last. Country and study are confounded wherever a country contributes one study, so the
only internal evidence that a population rather than a protocol is being measured comes from the six countries with
two or more studies: if studies within a country disagree as much as countries do, the population reading fails.

Writes results/s18/population_transfer.{tsv,md}.
"""
import os
import numpy as np, pandas as pd

P = os.environ.get("CHECKGM_PROJECT") or "/home/classes/bios/270/khoa/meta2/project"
OUT = f"{P}/results/s18"
LAYERS = {"taxonomy_family": "family", "taxonomy_genus": "genus", "taxonomy_species": "species"}
MIN_N = 40
IND_COUNTRIES = ["USA", "GBR", "FRA", "DEU", "DNK", "NLD", "ISR", "ESP", "SWE", "IRL"]


def load() -> pd.DataFrame:
    S = pd.read_csv(f"{P}/work/s11/pipeB/samples.tsv", sep="\t", low_memory=False)
    S = S[S.role.isin(["reference_pool", "heldout_healthy"])]
    M = pd.read_parquet(f"{P}/work/s12/B/oos_summaries.parquet")
    rows = []
    for layer, lab in LAYERS.items():
        m = M[M.layer == layer].set_index("sample")
        d = S.assign(frac=S.sample_key.map(m.frac_outside_raw), n_assessed=S.sample_key.map(m.n_assessed))
        rows.append(d.dropna(subset=["frac"]).assign(layer=lab))
    return pd.concat(rows, ignore_index=True)


def main():
    os.makedirs(OUT, exist_ok=True)
    D = load()
    # mean per sample within a study, then median over a country's studies
    per_study = D.groupby(["layer", "country", "study"]).agg(
        n=("frac", "size"), mean_frac=("frac", "mean"), n_assessed=("n_assessed", "median")).reset_index()
    per_country = per_study.groupby(["layer", "country"]).agg(
        studies=("study", "nunique"), n=("n", "sum"), pct=("mean_frac", "median"),
        spread=("mean_frac", lambda s: (s.max() - s.min()) * 100 if len(s) > 1 else np.nan),
        feats=("n_assessed", "median")).reset_index()
    per_country["pct"] = 100 * per_country.pct
    keep = per_country[per_country.n >= MIN_N]

    wide = keep.pivot_table(index="country", columns="layer", values="pct").round(1)
    meta = keep[keep.layer == "family"].set_index("country")[["studies", "n", "feats"]]
    sprd = keep.pivot_table(index="country", columns="layer", values="spread").round(1)
    T = meta.join(wide).join(sprd.add_suffix("_spread")).sort_values("family", ascending=False)
    T.to_csv(f"{OUT}/population_transfer.tsv", sep="\t")

    # Ranges are not comparable across groups of different size -- the range over 17 countries is wider than the range
    # over two studies for arithmetic reasons alone. Standard deviations are, so the confound is judged on those: the
    # spread of study means within a country against the spread of country means.
    multi = keep[(keep.studies >= 2) & keep.spread.notna()]
    dec = []
    for layer in LAYERS.values():
        st = per_study[(per_study.layer == layer) & (per_study.n >= 20)]
        mm = st.groupby("country").filter(lambda g: len(g) > 1)
        within = float(np.sqrt(np.mean([g.mean_frac.var(ddof=1) for _, g in mm.groupby("country")])) * 100)
        betw = float(st.groupby("country").mean_frac.mean().std(ddof=1) * 100)
        ind = st[st.country.isin(IND_COUNTRIES)].mean_frac.median() * 100
        top = keep[(keep.layer == layer)].nlargest(3, "pct")
        dec.append({"layer": layer, "sd_between_countries": round(betw, 2), "sd_within_country": round(within, 2),
                    "industrialised_median": round(ind, 1),
                    "top3_excess_in_within_SDs": round(float((top.pct.mean() - ind) / within), 1)})
    DEC = pd.DataFrame(dec)
    L = ["# Does the reference range transfer across the baseline's own populations? (read pipeline)\n",
         f"Reference adults only, each scored against a baseline rebuilt without its own study. "
         f"{int(meta.n.sum())} samples, {int(per_study[per_study.layer=='family'].study.nunique())} studies, "
         f"{len(T)} countries with at least {MIN_N} samples. 5% outside is expected by construction.\n",
         "Statistic: mean per-sample fraction outside within a study, then the median over a country's studies. "
         "The per-sample median is not usable here: with a median of "
         f"{int(meta.feats.median())} assessed families a 5% expectation is about one feature, so over half the "
         "samples of several countries have exactly zero outside and the median is 0.0 regardless of fit.\n",
         T.round(1).to_markdown(), "\n## Is it the population or the study?\n",
         "Country and study coincide wherever a country contributes one study. Where a country contributes two or "
         "more, the spread between its own studies measures what a protocol alone can do:\n",
         multi.pivot_table(index="country", columns="layer", values="spread").round(1).to_markdown(),
         "\nJudged on standard deviations rather than ranges, which depend on how many groups are compared:\n",
         DEC.to_markdown(index=False),
         "\nStudy-to-study variation within a country is almost as large as variation between countries "
         f"({DEC.sd_within_country.iloc[0]:.2f} against {DEC.sd_between_countries.iloc[0]:.2f} points on families), so "
         "most of the spread in the table above cannot be attributed to population rather than protocol. The exception "
         "is the top of it: Madagascar, Mongolia and India sit about "
         f"{DEC.top3_excess_in_within_SDs.iloc[0]:.0f} within-country standard deviations above the industrialised "
         "median, which protocol alone is not observed to do.\n"]
    # A single layer can be elevated for a reason that has nothing to do with the population: India's families are
    # scored on 16 features, where one feature is six points, and its elevation does not survive at the species layer.
    # Requiring the excess to appear on all three is the cheapest guard against reading such a cell as a finding.
    base = {l: T[T.index.isin(IND_COUNTRIES)][l].median() for l in LAYERS.values()}
    exc = T[list(LAYERS.values())].sub(pd.Series(base), axis=1).round(1)
    exc["replicates"] = np.where(exc.min(axis=1) >= 2.0, "all three layers",
                          np.where(exc.species >= 2.0, "genus and species", "one layer only"))
    exc.insert(0, "assessed_families", T.feats)
    L += ["\n## Does the excess replicate across layers?\n",
          "Excess over the industrialised median (" +
          ", ".join(f"{k} {v:.1f}%" for k, v in base.items()) + "), in points:\n",
          exc.sort_values("family", ascending=False).head(8).to_markdown(),
          "\nOnly Madagascar and Mongolia are elevated on all three layers. India is elevated on families alone, on "
          f"{int(T.loc['IND', 'feats'])} assessed features where one feature is six points, and sits at the "
          "industrialised median on species; it is not counted. Japan and China run the other way, at the "
          "industrialised level on families and two to five points above it on genus and species, which is the more "
          "surprising result of the two: a reference range can transfer at one rank and not at another.\n"]
    # figure data: one row per country per layer, plus the span of that country's own study means, which is the
    # protocol yardstick the panel is read against
    FIG = os.environ.get("CHECKGM_FIGDATA") or f"{P}/figures/v2/data"
    span = per_study[per_study.layer == "family"].groupby("country").mean_frac.agg(["min", "max", "size"]) * 100
    span["size"] = span["size"] / 100
    f = keep[["layer", "country", "pct", "studies", "n"]].copy()
    f["replicates"] = f.country.map(exc["replicates"])
    f = f.join(span[["min", "max"]].rename(columns={"min": "study_lo", "max": "study_hi"}), on="country")
    f["expected"] = 5.0
    f.to_csv(f"{FIG}/f5_countries.csv", index=False)
    print(f"wrote {FIG}/f5_countries.csv: {f.country.nunique()} countries x {f.layer.nunique()} layers")
    open(f"{OUT}/population_transfer.md", "w").write("\n".join(L))
    print("\n".join(L))


if __name__ == "__main__":
    main()
