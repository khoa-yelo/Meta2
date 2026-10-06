"""Figure 8 — case/control performance (double column, 2 x 2). Top row, assembly pipeline: (a) leave-one-study-out AUROC in
each held-out study with 95 % bootstrap CIs for the four feature sets; (b) mean AUROC over the studies with the within-study
cross-validated AUROC as a ceiling. Bottom row, read pipeline (taxonomy features): (c) percentiles vs raw CLR in each held-out
study; (d) mean AUROC by condition of the held-out study. Feature sets carry the four categorical colours (tokens.CAT: violet,
green, magenta, yellow), which are never used for deviation direction. The pooled ROC over all held-out predictions, drawn in an
earlier version, is not shown: it mixed studies of different size and prevalence and added nothing to (a) and (b)."""
import os, sys
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import tokens as T

P = T.P
SETS = ["reference_relative", "raw_clr", "alpha_diversity", "health_index"]
COL = dict(zip(SETS, T.CAT))


def make(out_dir):
    T.print_profile()
    R = pd.read_csv(f"{P}/results/s8/loso_auroc.tsv", sep="\t")
    R = R[R["feature_set"].isin(SETS)]
    per = R[R["study"] != "POOLED"].copy()
    per["label"] = per["study"] + "  " + per["condition"].str.replace("Inflammatory bowel disease", "IBD").str.replace("Parkinson's disease", "Parkinson's").str.replace("Rheumatoid arthritis", "RA").str.replace("C. difficile infection", "CDI").str.replace("Colorectal neoplasia", "CRC") + " (" + per["n_test"].astype(str) + ")"
    order = per[per["feature_set"] == "reference_relative"].sort_values("auroc_loso", ascending=False)["study"].tolist()
    fig = plt.figure(figsize=(T.DOUBLE_IN, 5.6))
    gs = fig.add_gridspec(2, 2, width_ratios=[2.1, 1.4], height_ratios=[1, 1.0], wspace=0.55, hspace=0.62, left=0.25, right=0.985, top=0.935, bottom=0.1)
    ax = fig.add_subplot(gs[0, 0])
    off = np.linspace(-0.27, 0.27, len(SETS))
    for k, fs in enumerate(SETS):
        d = per[per["feature_set"] == fs].set_index("study").reindex(order)
        yy = np.arange(len(order)) + off[k]
        ax.hlines(yy, d["ci_low"], d["ci_high"], color=COL[fs], lw=0.9, zorder=2)
        ax.scatter(d["auroc_loso"], yy, s=9, color=COL[fs], edgecolor=T.INK["surface"], linewidth=0.3, zorder=3, label=T.CAT_LABELS[fs])
    ax.axvline(0.5, color=T.INK["axis"], lw=0.6, zorder=1)
    ax.set_yticks(np.arange(len(order))); ax.set_yticklabels(per.drop_duplicates("study").set_index("study").reindex(order)["label"]); ax.invert_yaxis()
    ax.set_xlim(0.15, 1.0); ax.set_xlabel("AUROC in the held-out study (95 % bootstrap CI)")
    ax.set_title("a  assembly pipeline:\neach study held out", loc="left")
    h, l = ax.get_legend_handles_labels(); h.append(plt.Line2D([], [], marker="|", ls="", ms=5, markeredgewidth=0.9, color=T.INK["primary"])); l.append("within-study ceiling (b)")
    ax.legend(h, l, loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=3, frameon=False, handletextpad=0.3, columnspacing=1.0, borderaxespad=0)
    # b summary: horizontal dot plot (labels never collide), within-study ceiling as a bar
    ax = fig.add_subplot(gs[0, 1])
    mac = per.groupby("feature_set")["auroc_loso"].mean().reindex(SETS); ceil = per.groupby("feature_set")["auroc_within_study_cv"].median().reindex(SETS)
    yy = np.arange(len(SETS))
    ax.hlines(yy, 0.5, mac, color=T.INK["axis"], lw=0.8, zorder=2)
    ax.scatter(mac, yy, s=22, color=[COL[s_] for s_ in SETS], edgecolor=T.INK["surface"], linewidth=0.3, zorder=3)
    ax.scatter(ceil, yy, marker="|", s=40, color=T.INK["primary"], linewidth=0.9, zorder=3)
    for yi, v in zip(yy, mac):
        ax.text(v, yi - 0.3, f"{v:.2f}", ha="center", va="bottom", fontsize=T.PT_MIN, color=T.INK["primary"])
    ax.axvline(0.5, color=T.INK["axis"], lw=0.6, zorder=1)
    ax.set_yticks(yy); ax.set_yticklabels(["percentiles", "raw abundances", "alpha diversity", "GMHI (genus approx.)"]); ax.invert_yaxis(); ax.tick_params(axis="y", length=0)
    ax.set_xlim(0.4, 0.8); ax.set_xticks([0.4, 0.5, 0.6, 0.7, 0.8]); ax.set_ylim(len(SETS) - 0.5, -0.8); ax.set_xlabel("mean AUROC, 11 studies"); ax.set_title("b  assembly pipeline:\nmean over studies", loc="left")
    # c, d — read pipeline: cross-study leave-one-study-out, reference-relative vs raw CLR (family + genus + species)
    cc = pd.read_csv(f"{P}/results/s11/pipelineB_case_control.tsv", sep="\t").set_index("study")
    b = pd.read_csv(f"{P}/results/s11/pipelineB_loso_ref_relative.tsv", sep="\t").set_index("study")["auroc"].to_frame("reference_relative")
    b["raw_clr"] = pd.read_csv(f"{P}/results/s11/pipelineB_loso_raw_clr.tsv", sep="\t").set_index("study")["auroc"]
    b = b.sort_values("reference_relative", ascending=False); cond = cc["condition"].reindex(b.index).replace({"schizofrenia": "schizophrenia"})
    n = (cc["n_case"] + cc["n_control"]).reindex(b.index)
    ax = fig.add_subplot(gs[1, 0]); yy = np.arange(len(b))
    ax.hlines(yy, b.min(axis=1), b.max(axis=1), color=T.INK["axis"], lw=0.9, zorder=2)
    for fs in ["raw_clr", "reference_relative"]:
        ax.scatter(b[fs], yy, s=11, color=COL[fs], edgecolor=T.INK["surface"], linewidth=0.3, zorder=3, label=T.CAT_LABELS[fs])
    ax.axvline(0.5, color=T.INK["axis"], lw=0.6, zorder=1)
    ax.set_yticks(yy); ax.set_yticklabels([f"{s}  {c} ({int(k)})" for s, c, k in zip(b.index, cond, n)]); ax.invert_yaxis()
    ax.set_xlim(0.15, 1.0); ax.set_xlabel("AUROC in the held-out study"); ax.set_title("c  read pipeline (taxonomy):\neach study held out", loc="left")
    ax.legend(loc="upper center", bbox_to_anchor=(0.45, -0.17), ncol=2, frameon=False, handletextpad=0.3, columnspacing=1.0, borderaxespad=0)
    # d — read pipeline: published health indices used as they are (results/s15/health_indices_B.tsv) against percentiles;
    # GMWI2 was trained on samples of 7 of the 14 studies, so its fair comparison is the 7 studies it never saw
    H = pd.read_csv(f"{P}/results/s15/health_indices_B.tsv", sep="\t"); Hs = H.pivot(index="study", columns="index", values="auroc_direct")
    base = H.drop_duplicates("study").set_index("study"); unseen = ~base["in_gmwi2_training"].astype(bool)
    ax = fig.add_subplot(gs[1, 1])
    rows = [("percentiles", base["percentiles"], COL["reference_relative"], "o"), ("raw abundances", base["raw_abundances"], COL["raw_clr"], "o"),
            ("Shannon diversity", Hs["Shannon"], T.CAT[2], "o"), ("GMHI (50 species)", Hs["GMHI"], T.CAT[3], "o"), ("GMWI2", Hs["GMWI2"], T.INK["primary"], "D")]
    # supervised health score on refmb output (pre-specified primary model; results/s16/health_score_B.tsv), every study scored by a model trained without it
    HS = pd.read_csv(f"{P}/results/s16/health_score_B.tsv", sep="\t"); HS = HS[HS.feature_set == "pct+summary+presence"].set_index("study")["auroc"]
    rows.insert(1, ("refmb health score", HS, COL["reference_relative"], "s"))
    yl = []; y = 0
    for gi, (glab, m) in enumerate([(f"all {len(base)} studies", pd.Series(True, index=base.index)), (f"{int(unseen.sum())} studies not used to train GMWI2", unseen)]):
        ax.text(0.405, y - 0.75, glab, fontsize=T.PT_MIN, color=T.INK["primary"], weight="bold", va="center")
        for lab, ser, col, mk in rows:
            v = ser[m.reindex(ser.index).fillna(False).astype(bool)].mean()
            ax.hlines(y, 0.5, v, color=T.INK["axis"], lw=0.8, zorder=2); ax.scatter(v, y, s=20, color=col, marker=mk, edgecolor=T.INK["surface"], linewidth=0.3, zorder=3)
            ax.text(v + 0.012, y, f"{v:.2f}", ha="left", va="center", fontsize=T.PT_MIN, color=T.INK["primary"]); yl.append((y, lab)); y += 1
        y += 1.2
    ax.axvline(0.5, color=T.INK["axis"], lw=0.6, zorder=1); ax.set_yticks([a for a, _ in yl]); ax.set_yticklabels([b for _, b in yl]); ax.tick_params(axis="y", length=0)
    ax.set_xlim(0.4, 0.9); ax.set_xticks([0.4, 0.5, 0.6, 0.7, 0.8, 0.9]); ax.set_ylim(y - 1.5, -1.4); ax.set_xlabel("mean AUROC over studies"); ax.set_title("d  read pipeline: health\nindices vs percentiles", loc="left")
    T.save(fig, "fig8_performance", out_dir); plt.close(fig)


if __name__ == "__main__":
    make(sys.argv[1] if len(sys.argv) > 1 else f"{P}/figures")
