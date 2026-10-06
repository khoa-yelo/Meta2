"""Figure 1 — the tool (double column). (a) Left-to-right schematic: a user's raw reads enter one of two pipelines
(assembly pipeline, read pipeline); each is normalized with the rules stored in its own healthy baseline and scored
against that baseline; the report is the same for both. Baseline sizes are read from the bundle manifests. (b) The user's
view: nine rows of the family-level report of one real sample (the C. difficile case of fig5_tree), each placed at its
percentile among healthy adults with the healthy range shaded; counts in the title are computed (example_case.py).
Pipelines are told apart by fill and outline (filled light grey = assembly, outlined dashed = read), never by hue:
blue and orange are reserved for the direction of a deviation (tokens.py, colour conventions)."""
import json, os, sys
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import tokens as T, example_case as EC
sys.path.insert(0, T.P); from checkgm.paths import bundle_A

P = T.P


def box(ax, x, y, w, h, title, lines, style, sub=None, title_col=T.INK["primary"], lw=0.8):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.015,rounding_size=0.05", facecolor=style["facecolor"], edgecolor=style["edgecolor"],
                                linewidth=lw, linestyle=style.get("linestyle", "-")))
    ax.text(x + 0.07, y + h - 0.07, title, fontsize=T.PT_BODY, weight="bold", va="top", color=title_col)
    ax.text(x + 0.07, y + h - 0.3, "\n".join(lines), fontsize=T.PT_MIN, va="top", color=T.INK["secondary"], linespacing=1.3)
    if sub:
        ax.text(x + 0.07, y + 0.06, sub, fontsize=T.PT_MIN, va="bottom", color=T.INK["muted"], style="italic", linespacing=1.2)


def arrow(ax, x0, y0, x1, y1, ls="-"):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=6, linewidth=0.7, color=T.INK["primary"], shrinkA=0, shrinkB=0, linestyle=ls))


def _numbers():
    a = json.load(open(f"{bundle_A()}/manifest.json")); b = json.load(open(f"{P}/refs/gut-reads-adult-global-v0.2-lenient/manifest.json"))
    return (a["n_samples"], a["n_studies"]), (b["n_samples"], b["n_studies"])


def _rows(d):
    """Nine rows of the case's family report, chosen by rule so the strip shows every kind of call: the two low families
    with the highest and the two with the lowest percentiles, the within families nearest the 10th, 50th, 90th and 96th
    percentiles, and the first expected-but-missing family (alphabetical), among families carried by at least 90 % of healthy adults in some size class, most prevalent first."""
    d = d.drop_duplicates("name").sort_values("name")
    # rows come from families that at least 90 % of healthy adults carry in some size class, so the strip shows familiar gut families
    ref = EC.reference(); pc = [c for c in ref.columns if c.startswith("prevalence_band_")]
    prev = ref[pc].max(axis=1); common = set(prev.index[prev >= 0.9]); d = d[d.index.isin(common)].copy(); d["prev"] = prev.reindex(d.index).to_numpy()
    lows = d[d["call"] == "low"].sort_values("prev", ascending=False).head(4)
    within = d[d["call"] == "within"].sort_values("prev", ascending=False).head(12)
    pick = [within.iloc[(within["percentile"] - q).abs().argsort().iloc[0]] for q in (10, 50, 75, 96)]
    miss = d[d["call"] == "expected_but_missing"].sort_values("prev", ascending=False).head(1)
    return pd.concat([lows, pd.DataFrame(pick), miss]).drop_duplicates("name").sort_values("percentile", na_position="first")
    lows = d[d["call"] == "low"].sort_values(["percentile", "name"], ascending=[False, True]); low = pd.concat([lows.head(2), lows.tail(2)])
    within = d[d["call"] == "within"]; pick = [within.iloc[(within["percentile"] - q).abs().argsort().iloc[0]] for q in (10, 50, 90, 96)]
    miss = d[d["call"] == "expected_but_missing"].head(1)
    out = pd.concat([low, pd.DataFrame(pick), miss]).drop_duplicates("name")
    return out.sort_values("percentile", na_position="first")


def make(out_dir):
    T.print_profile(); (nA, sA), (nB, sB) = _numbers()
    fam = EC.case_families(); c = EC.counts(fam); smp = _rows(fam)
    fig = plt.figure(figsize=(T.DOUBLE_IN, 3.6))
    # ---------------- a: flow ----------------
    ax = fig.add_axes([0, 0.42, 1, 0.58]); ax.set_xlim(0, 10); ax.set_ylim(0, 2.7); ax.axis("off")
    ax.text(0.08, 2.68, "a  one sample, two pipelines, one kind of report", fontsize=T.PT_TITLE, va="top")
    yA, yB, h = 1.3, 0.04, 1.2
    SA, SB = T.PIPE["assembly"], T.PIPE["read"]
    box(ax, 0.06, yB, 1.14, yA + h - yB, "Your sample", ["gut metagenome,", "raw reads (FASTQ)", "", "one sample is", "enough: the", "comparison", "group ships", "with the tool"], {"facecolor": T.INK["surface"], "edgecolor": T.INK["axis"]})
    box(ax, 1.55, yA, 2.15, h, "Assembly pipeline", ["pinned container; reproduces", "MGnify v5: metaSPAdes/MEGAHIT,", "Prodigal, FragGeneScan; eggNOG,", "KOfam, Pfam. Hours per sample"], SA)
    box(ax, 1.55, yB, 2.15, h, "Read pipeline", ["MetaPhlAn 3.0.14; reproduces", "curatedMetagenomicData 3;", "HUMAnN 3 for genes, pathways.", "Minutes per sample"], SB)
    box(ax, 4.0, yA, 1.6, h, "Normalize", ["with the baseline's rules:", "genome equivalents, CLR,", "assembly-size class"], SA)
    box(ax, 4.0, yB, 1.6, h, "Normalize", ["with the baseline's rules:", "relative abundance, CLR,", "read-depth class"], SB)
    box(ax, 5.85, yA, 2.25, h, "Healthy baseline, assembly", [f"{nA:,} healthy adults, {sA} studies", "family, genus, KO, Pfam, module;", "per feature: prevalence,", "percentiles, bootstrap CIs"], SA, sub="statistics only,\nno sample-level data")
    box(ax, 5.85, yB, 2.25, h, "Healthy baseline, read", [f"{nB:,} healthy adults, {sB} studies", "family, genus, species, KO,", "pathway; per feature: prevalence,", "percentiles, bootstrap CIs"], SB, sub="statistics only,\nno sample-level data")
    # report box: two short lists with hanging indents
    rx, rw = 8.4, 1.55
    ax.add_patch(FancyBboxPatch((rx, yB), rw, yA + h - yB, boxstyle="round,pad=0.015,rounding_size=0.05", facecolor=T.INK["surface"], edgecolor=T.INK["primary"], linewidth=0.8))
    ax.text(rx + 0.07, yA + h - 0.07, "Report", fontsize=T.PT_BODY, weight="bold", va="top")
    ax.text(rx + 0.07, yA + h - 0.3, "per feature", fontsize=T.PT_MIN, va="top", color=T.INK["primary"])
    ax.text(rx + 0.19, yA + h - 0.46, "percentile among\nhealthy adults\nwithin, low or high\nexpected but missing", fontsize=T.PT_MIN, va="top", color=T.INK["secondary"], linespacing=1.3)
    ax.text(rx + 0.07, yA + h - 1.18, "per sample", fontsize=T.PT_MIN, va="top", color=T.INK["primary"])
    ax.text(rx + 0.19, yA + h - 1.34, "share outside the\nhealthy range, with\na test against the\nexpected ~5 %", fontsize=T.PT_MIN, va="top", color=T.INK["secondary"], linespacing=1.3)
    for y, ls in ((yA + h / 2, "-"), (yB + h / 2, SB["linestyle"])):
        for x0, x1 in [(3.7, 4.0), (5.6, 5.85), (8.1, 8.4)]:
            arrow(ax, x0, y, x1, y, ls)
    arrow(ax, 1.2, (yA + yB + h) / 2, 1.55, yA + h / 2, "-"); arrow(ax, 1.2, (yA + yB + h) / 2, 1.55, yB + h / 2, SB["linestyle"])
    # ---------------- b: user's view (full lower width) ----------------
    cx = fig.add_axes([0.17, 0.115, 0.56, 0.215])
    cx.axvspan(2.5, 97.5, color=T.INK["grid"], lw=0, zorder=0); cx.axvline(50, color=T.INK["axis"], lw=0.5, zorder=1)
    y = np.arange(len(smp))
    for yi, (_, r) in zip(y, smp.iterrows()):
        if r["call"] == "expected_but_missing":
            cx.scatter([0], [yi], s=18, zorder=3, **T.MISSING_RING); cx.text(103, yi, "expected but missing", va="center", ha="left", fontsize=T.PT_MIN, color=T.INK["muted"])
            continue
        p = r["percentile"]; col = T.DEV["neutral"] if r["call"] == "within" else (T.DEV["low_deep"] if p < 1 else T.DEV["low_light"]) if r["call"] == "low" else (T.DEV["high_deep"] if p > 99 else T.DEV["high_light"])
        cx.scatter([p], [yi], s=18, color=col, edgecolor=T.INK["surface"], linewidth=0.3, zorder=3)
        cx.text(103, yi, f"{r['call']}  ({p:.2f})" if p < 2.5 or p > 97.5 else f"within  ({p:.0f})", va="center", ha="left", fontsize=T.PT_MIN, color=T.INK["secondary"])
    cx.set_yticks(y); cx.set_yticklabels(smp["name"], style="italic"); cx.tick_params(axis="y", length=0); cx.set_ylim(len(smp) - 0.5, -0.5)
    cx.set_xlim(-1, 101); cx.set_xticks([0, 25, 50, 75, 100]); cx.tick_params(axis="x", length=2, width=0.4)
    cx.set_xlabel("percentile among healthy adults (shaded = healthy range, 2.5th–97.5th)")
    cx.text(-0.02, -1.0, "family", transform=cx.get_yaxis_transform(), ha="right", va="center", fontsize=T.PT_MIN, color=T.INK["muted"])
    cx.text(103, -1.0, "call (percentile)", ha="left", va="center", fontsize=T.PT_MIN, color=T.INK["muted"])
    cx.set_title(f"b  the user's view: nine rows of one C. difficile report ({c['n_assessed']} assessed, {c['n_missing']} expected but missing)", loc="left", pad=8)
    T.save(fig, "fig1_schematic", out_dir); plt.close(fig)


if __name__ == "__main__":
    make(sys.argv[1] if len(sys.argv) > 1 else f"{P}/figures")
