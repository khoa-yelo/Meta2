"""Figure 12 — disease deviations within age, sex, BMI and region (double column). Colorectal cancer, pipeline B.
a: median shift of the reference percentile (cases minus same-study controls, median over studies) of the taxa that
are consistent in the disease atlas, within each level; five-step diverging scale (tokens), as in Figure 10.
b: share of all atlas-consistent features (taxa, KOs, pathways; 203) whose shift within the level has the atlas direction.
c: median absolute shift of those features within the level. Levels are labelled with the number of studies."""
import os, sys
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import tokens as T

P = T.P
BOUNDS = [-100, -20, -5, 5, 20, 100]; STEP_LABELS = ["≤ −20", "−20 to −5", "−5 to 5", "5 to 20", "≥ 20"]
CONDITION = "CRC"
LEVELS = [("age_bin", "40-64"), ("age_bin", "65+"), ("sex", "female"), ("sex", "male"), ("bmi_class", "18.5-25"), ("bmi_class", "25-30"),
          ("region", "Europe"), ("region", "N. America"), ("region", "E. Asia"), ("region", "S. Asia")]
VAR_LABEL = {"age_bin": "age", "sex": "sex", "bmi_class": "BMI", "region": "region"}
RANK = {"taxonomy_family": 0, "taxonomy_genus": 1, "taxonomy_species": 2}


def load():
    L = pd.read_csv(f"{P}/results/s14/disease_strata_{'B'}_condition_levels.tsv", sep="\t", dtype={"feature_id": str})
    L = L[L["condition"] == CONDITION].copy()
    L["same"] = np.sign(L["median_shift"]) == L["consistent"].map({"up": 1, "down": -1})
    L["key"] = list(zip(L["variable"], L["level"]))
    summ = L.groupby("key").agg(features=("same", "size"), share=("same", "mean"), mabs=("median_shift", lambda s: s.abs().median()), studies=("n_studies", "max")).reindex(LEVELS)
    tax = pd.read_parquet(f"{P}/resources/backbone/ncbi_taxonomy.parquet", columns=["taxid", "name"]); names = dict(zip(tax["taxid"].astype(str), tax["name"]))
    X = L[L["layer"].isin(RANK)].copy(); X["name"] = X["feature_id"].map(lambda i: names.get(i, i)); X["rank"] = X["layer"].map(RANK)
    M = X.pivot_table(index=["rank", "name"], columns="key", values="median_shift").reindex(columns=LEVELS)
    order = M.assign(m=M.mean(axis=1)).reset_index().sort_values(["rank", "m"], ascending=[False, False]).set_index(["rank", "name"]).index
    return M.loc[order], summ, int(L["feature_id"].nunique())


def make(out_dir):
    T.print_profile()
    M, summ, nfeat = load()
    cmap = ListedColormap([T.DEV[k] for k in T.DEV_ORDER]); norm = BoundaryNorm(BOUNDS, cmap.N)
    labels = [f"{VAR_LABEL[v]}  {l.replace('-', '–')} ({int(summ.loc[[(v, l)], 'studies'].iloc[0])})" for v, l in LEVELS]
    breaks = [k for k in range(1, len(LEVELS)) if LEVELS[k][0] != LEVELS[k - 1][0]]
    fig = plt.figure(figsize=(T.DOUBLE_IN, 3.3))
    gs = fig.add_gridspec(1, 3, width_ratios=[10, 5.2, 5.2], wspace=0.12, left=0.245, right=0.985, top=0.92, bottom=0.39)
    # a: heatmap of the atlas-consistent taxa
    ax = fig.add_subplot(gs[0])
    ax.pcolormesh(np.arange(M.shape[1] + 1), np.arange(M.shape[0] + 1), np.ma.masked_invalid(M.values), cmap=cmap, norm=norm, edgecolors=T.INK["surface"], linewidth=0.8)
    rank_lab = {0: "family", 1: "genus", 2: "species"}
    ax.set_yticks(np.arange(M.shape[0]) + 0.5); ax.set_yticklabels([n for _, n in M.index], style="italic", color=T.INK["secondary"])
    ax.set_xticks(np.arange(M.shape[1]) + 0.5); ax.set_xticklabels(labels, rotation=90, color=T.INK["secondary"])
    ax.tick_params(length=0); ax.set_xlim(0, M.shape[1]); ax.set_ylim(0, M.shape[0])
    for sp in ax.spines.values():
        sp.set_visible(False)
    for x in breaks:
        ax.axvline(x, color=T.INK["secondary"], lw=0.6, zorder=4)
    ranks = np.array([r for r, _ in M.index])
    for y in np.where(ranks[1:] != ranks[:-1])[0] + 1:
        ax.axhline(y, color=T.INK["secondary"], lw=0.6, zorder=4)
    for r in sorted(set(ranks)):
        idx = np.where(ranks == r)[0]
        fig.text(0.012, ax.get_position().y0 + ax.get_position().height * (idx.mean() + 0.5) / M.shape[0], rank_lab[r], rotation=90, ha="left", va="center", fontsize=T.PT_MIN, color=T.INK["muted"])
    ax.set_title("a  the 18 taxa consistent across studies", loc="left")
    # b, c: all atlas-consistent features, per level (top to bottom in the column order of a)
    y = np.arange(len(LEVELS))[::-1]
    bx = fig.add_subplot(gs[1]); cx = fig.add_subplot(gs[2])
    for a_, vals, lab, xlim, ticks, fmt in [(bx, summ["share"].values * 100, "features keeping the\ncross-study direction (%)", (50, 112), [50, 75, 100], "{:.0f}"),
                                           (cx, summ["mabs"].values, "median |shift|\n(percentile points)", (0, 31), [0, 10, 20], "{:.0f}")]:
        a_.hlines(y, xlim[0], vals, color=T.INK["axis"], lw=0.8, zorder=2)
        a_.scatter(vals, y, s=14, color=T.INK["primary"], edgecolor=T.INK["surface"], linewidth=0.3, zorder=3)
        for v, yy in zip(vals, y):
            a_.text(v + (xlim[1] - xlim[0]) * 0.035, yy, fmt.format(v), va="center", ha="left", fontsize=T.PT_MIN, color=T.INK["primary"])
        for k in breaks:
            a_.axhline(len(LEVELS) - k - 0.5, color=T.INK["grid"], lw=0.6, zorder=1)
        a_.set_xlim(*xlim); a_.set_ylim(-0.6, len(LEVELS) - 0.4); a_.set_xticks(ticks); a_.set_xlabel(lab)
        a_.set_yticks(y); a_.tick_params(axis="y", length=0); a_.tick_params(axis="x", length=2, width=0.4)
    bx.set_yticklabels(labels, color=T.INK["secondary"]); cx.set_yticklabels([])
    bx.axvline(100, color=T.INK["grid"], lw=0.6, zorder=1)
    bx.set_title(f"b  all {nfeat} features", loc="left"); cx.set_title("c  size of shift", loc="left")
    # panel a sits left of b's tick labels: leave room by shrinking a
    pa = ax.get_position(); ax.set_position([pa.x0, pa.y0, pa.width * 0.6, pa.height])
    h = [Patch(facecolor=T.DEV[k], edgecolor="none", label=l) for k, l in zip(T.DEV_ORDER, STEP_LABELS)]
    pa = ax.get_position()   # colour legend directly under panel a, left-aligned with it
    fig.legend(handles=h, loc="lower left", bbox_to_anchor=(pa.x0 - 0.2, 0.005), ncol=5, frameon=False, handlelength=1.2, handletextpad=0.4, columnspacing=1.0, borderaxespad=0,
               title="shift in percentile points, cases minus same-study controls; blank = not assessed", title_fontsize=T.PT_MIN, alignment="left")
    T.save(fig, "fig12_disease_strata", out_dir); plt.close(fig)


if __name__ == "__main__":
    make(sys.argv[1] if len(sys.argv) > 1 else f"{P}/figures")
