"""Figure 8 — case/control performance (double column). a: leave-one-study-out AUROC per held-out study with 95%
bootstrap CIs for the four feature sets (forest). b: pooled ROC curves over all held-out predictions. c: summary
(macro-mean AUROC) with the within-study ceiling."""
import os, sys
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, roc_auc_score
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
    fig = plt.figure(figsize=(T.DOUBLE_IN, 5.9))
    gs = fig.add_gridspec(2, 3, width_ratios=[2.1, 1.0, 1.3], height_ratios=[1, 1.0], wspace=0.7, hspace=0.6, left=0.25, right=0.985, top=0.94, bottom=0.1)
    ax = fig.add_subplot(gs[0, 0])
    off = np.linspace(-0.27, 0.27, len(SETS))
    for k, fs in enumerate(SETS):
        d = per[per["feature_set"] == fs].set_index("study").reindex(order)
        yy = np.arange(len(order)) + off[k]
        ax.hlines(yy, d["ci_low"], d["ci_high"], color=COL[fs], lw=0.9, zorder=2)
        ax.scatter(d["auroc_loso"], yy, s=9, color=COL[fs], edgecolor=T.INK["surface"], linewidth=0.3, zorder=3, label=T.CAT_LABELS[fs].replace("reference-relative", "percentiles (reference-relative)"))
    ax.axvline(0.5, color=T.INK["axis"], lw=0.6, zorder=1)
    ax.set_yticks(np.arange(len(order))); ax.set_yticklabels(per.drop_duplicates("study").set_index("study").reindex(order)["label"]); ax.invert_yaxis()
    ax.set_xlim(0.15, 1.0); ax.set_xlabel("AUROC in the held-out study (95 % bootstrap CI)")
    ax.set_title("a  assembly pipeline:\neach study held out", loc="left")
    h, l = ax.get_legend_handles_labels(); h.append(plt.Line2D([], [], marker="|", ls="", ms=5, markeredgewidth=0.9, color=T.INK["primary"])); l.append("within-study ceiling (c)")
    ax.legend(h, l, loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=3, frameon=False, handletextpad=0.3, columnspacing=1.0, borderaxespad=0)
    # b ROC
    ax = fig.add_subplot(gs[0, 1])
    pf = f"{P}/results/s8/loso_predictions.parquet"
    if os.path.exists(pf):
        pr = pd.read_parquet(pf)
        for fs in SETS:
            d = pr[pr["feature_set"] == fs]
            fpr, tpr, _ = roc_curve(d["y"], d["score"]); a = roc_auc_score(d["y"], d["score"])
            ax.plot(fpr, tpr, color=COL[fs], lw=1.0, zorder=3)
            ax.text(1.0, {"reference_relative": -0.52, "raw_clr": -0.64, "alpha_diversity": -0.76, "health_index": -0.88}[fs], f"{T.CAT_LABELS[fs].replace('reference-relative', 'percentiles')} {a:.2f}", transform=ax.transAxes, ha="right", va="center", fontsize=T.PT_MIN, color=T.INK["secondary"])
    else:
        ax.text(0.5, 0.5, "pooled predictions\nnot yet available", ha="center", va="center", fontsize=T.PT_MIN, color=T.INK["muted"])
    ax.plot([0, 1], [0, 1], color=T.INK["axis"], lw=0.6, zorder=1)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_xticks([0, 0.5, 1]); ax.set_yticks([0, 0.5, 1]); ax.set_aspect("equal")
    ax.set_xlabel("false positive rate"); ax.set_ylabel("true positive rate"); ax.set_title("b  assembly pipeline:\npooled ROC", loc="left"); ax.set_anchor("N")
    # c summary: horizontal dot plot (labels never collide), within-study ceiling as a bar
    ax = fig.add_subplot(gs[0, 2])
    mac = per.groupby("feature_set")["auroc_loso"].mean().reindex(SETS); ceil = per.groupby("feature_set")["auroc_within_study_cv"].median().reindex(SETS)
    yy = np.arange(len(SETS))
    ax.hlines(yy, 0.5, mac, color=T.INK["axis"], lw=0.8, zorder=2)
    ax.scatter(mac, yy, s=22, color=[COL[s_] for s_ in SETS], edgecolor=T.INK["surface"], linewidth=0.3, zorder=3)
    ax.scatter(ceil, yy, marker="|", s=40, color=T.INK["primary"], linewidth=0.9, zorder=3)
    for yi, v in zip(yy, mac):
        ax.text(v, yi - 0.3, f"{v:.2f}", ha="center", va="bottom", fontsize=T.PT_MIN, color=T.INK["primary"])
    ax.axvline(0.5, color=T.INK["axis"], lw=0.6, zorder=1)
    ax.set_yticks(yy); ax.set_yticklabels(["percentiles", "raw CLR", "alpha diversity", "GMHI (genus approx.)"]); ax.invert_yaxis(); ax.tick_params(axis="y", length=0)
    ax.set_xlim(0.4, 0.8); ax.set_xticks([0.4, 0.5, 0.6, 0.7, 0.8]); ax.set_ylim(len(SETS) - 0.5, -0.8); ax.set_xlabel("mean AUROC, 11 studies"); ax.set_title("c  assembly pipeline:\nmean over studies", loc="left")
    # d, e — pipeline B: cross-study leave-one-study-out, reference-relative vs raw CLR (family + genus + species)
    cc = pd.read_csv(f"{P}/results/s11/pipelineB_case_control.tsv", sep="\t").set_index("study")
    b = pd.read_csv(f"{P}/results/s11/pipelineB_loso_ref_relative.tsv", sep="\t").set_index("study")["auroc"].to_frame("reference_relative")
    b["raw_clr"] = pd.read_csv(f"{P}/results/s11/pipelineB_loso_raw_clr.tsv", sep="\t").set_index("study")["auroc"]
    b = b.sort_values("reference_relative", ascending=False); cond = cc["condition"].reindex(b.index).replace({"schizofrenia": "schizophrenia"})
    n = (cc["n_case"] + cc["n_control"]).reindex(b.index)
    ax = fig.add_subplot(gs[1, 0]); yy = np.arange(len(b))
    ax.hlines(yy, b.min(axis=1), b.max(axis=1), color=T.INK["axis"], lw=0.9, zorder=2)
    for fs in ["raw_clr", "reference_relative"]:
        ax.scatter(b[fs], yy, s=11, color=COL[fs], edgecolor=T.INK["surface"], linewidth=0.3, zorder=3, label=T.CAT_LABELS[fs].replace("reference-relative", "percentiles (reference-relative)"))
    ax.axvline(0.5, color=T.INK["axis"], lw=0.6, zorder=1)
    ax.set_yticks(yy); ax.set_yticklabels([f"{s}  {c} ({int(k)})" for s, c, k in zip(b.index, cond, n)]); ax.invert_yaxis()
    ax.set_xlim(0.15, 1.0); ax.set_xlabel("AUROC in the held-out study"); ax.set_title("d  read pipeline (taxonomy):\neach study held out", loc="left")
    ax.legend(loc="upper center", bbox_to_anchor=(0.45, -0.17), ncol=2, frameon=False, handletextpad=0.3, columnspacing=1.0, borderaxespad=0)
    ax = fig.add_subplot(gs[1, 1:]); crc = cond.eq("CRC").to_numpy()
    groups = [(f"all ({len(b)})", np.ones(len(b), bool)), (f"colorectal cancer ({crc.sum()})", crc), (f"other conditions ({(~crc).sum()})", ~crc)]
    for gi, (lab, m) in enumerate(groups):
        for k, fs in enumerate(["reference_relative", "raw_clr"]):
            v = b.loc[m, fs].mean(); y = gi + (-0.17 if k == 0 else 0.17)
            ax.hlines(y, 0.5, v, color=T.INK["axis"], lw=0.8, zorder=2); ax.scatter(v, y, s=22, color=COL[fs], edgecolor=T.INK["surface"], linewidth=0.3, zorder=3)
            ax.text(v + 0.012, y, f"{v:.2f}", ha="left", va="center", fontsize=T.PT_MIN, color=T.INK["primary"])
    ax.axvline(0.5, color=T.INK["axis"], lw=0.6, zorder=1); ax.set_yticks(range(len(groups))); ax.set_yticklabels([g for g, _ in groups]); ax.invert_yaxis(); ax.tick_params(axis="y", length=0)
    ax.set_xlim(0.4, 0.9); ax.set_xticks([0.4, 0.5, 0.6, 0.7, 0.8, 0.9]); ax.set_ylim(len(groups) - 0.5, -0.5); ax.set_xlabel("mean AUROC over held-out studies"); ax.set_title("e  read pipeline:\nmean by condition", loc="left")
    T.save(fig, "fig8_performance", out_dir); plt.close(fig)


if __name__ == "__main__":
    make(sys.argv[1] if len(sys.argv) > 1 else f"{P}/figures")
