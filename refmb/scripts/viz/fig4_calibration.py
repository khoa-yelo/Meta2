"""Figure 4 — calibration check (double column, 2 x 2). Rows: assembly pipeline (a, b), read pipeline (c, d).
Left: mean share of a healthy sample's assessed features outside the healthy range for each independent healthy cohort
(never part of the baseline), per layer, against the 5 % expectation and the pre-set 8 % limit; a cohort above the limit is
drawn in the deviation-high orange (an exceedance). Right: each baseline study scored against a baseline rebuilt without it.
Cohorts are told apart by marker shape (one hue); the cohort legend sits in the empty right half of a and c."""
import json, os, sys
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import tokens as T
sys.path.insert(0, T.P); from refmb.paths import VTAG as _VTAG

P = T.P
LAYER_LABEL = {"taxonomy_family": "family", "taxonomy_genus": "genus", "ko_eggnog": "KO (eggNOG)", "ko_kofam": "KO (KOfam)", "pfam": "Pfam", "module": "module",
               "taxonomy_species": "species", "ko_humann": "KO (HUMAnN)", "pathway_humann": "pathway"}
LAYERS = ["taxonomy_family", "taxonomy_genus", "ko_eggnog", "ko_kofam", "pfam", "module"]
LAYERS_B = ["taxonomy_family", "taxonomy_genus", "taxonomy_species", "ko_humann", "pathway_humann"]
XT = dict(xlim=(0, 0.17), xticks=[0, 0.05, 0.10, 0.15], xticklabels=["0", "5 %", "10 %", "15 %"])


def heldout(ax, d, layers, value, title):
    studies = sorted(d["study"].unique()); y = np.arange(len(layers)); mk = dict(zip(studies, ["o", "s", "D", "^", "v"]))
    ax.axvspan(0.03, 0.08, color=T.INK["grid"], lw=0, zorder=0); ax.axvline(0.05, color=T.INK["axis"], lw=0.6, zorder=1)
    for s in studies:
        v = d[d["study"] == s].set_index("layer").reindex(layers)[value]
        ax.scatter(v, y, marker=mk[s], s=14, facecolor=np.where(v > 0.08, T.DEV["high_deep"], T.INK["primary"]), edgecolor=T.INK["surface"], linewidth=0.4, zorder=3)
        ax.scatter([], [], marker=mk[s], s=14, color=T.INK["primary"], label=f"{s} (n = {int(d[d.study == s].n.iloc[0])})")
    ax.legend(loc="upper right", bbox_to_anchor=(1.0, 1.0), ncol=1, frameon=False, handletextpad=0.2, borderaxespad=0, labelspacing=0.35)
    ax.set_yticks(y); ax.set_yticklabels([LAYER_LABEL[l] for l in layers]); ax.set_ylim(len(layers) - 0.4, -0.9); ax.set(**XT)
    ax.set_xlabel("share of assessed features outside the healthy range"); ax.set_title(title, loc="left")
    ax.text(0.049, -0.55, "5 %", fontsize=T.PT_MIN, color=T.INK["secondary"], va="bottom", ha="right"); ax.text(0.083, -0.55, "8 % limit", fontsize=T.PT_MIN, color=T.INK["secondary"], va="bottom", ha="left")


def loso_panel(ax, loso, layers, col, title):
    n = loso["study"].nunique()
    for i, l in enumerate(layers):
        v = np.sort(loso.loc[loso["layer"] == l, col].to_numpy())
        if not len(v):
            continue
        q1, q2, q3 = np.percentile(v, [25, 50, 75])
        ax.plot([q1, q3], [i, i], color=T.INK["axis"], lw=2.5, solid_capstyle="butt", zorder=2); ax.scatter([q2], [i], s=12, color=T.INK["primary"], zorder=3)
        yy = np.full(len(v), i) + np.linspace(-0.22, 0.22, len(v)); inside = v <= 0.168
        ax.scatter(v[inside], yy[inside], s=3, color=T.INK["muted"], zorder=1, linewidth=0)
        ax.scatter(np.full((~inside).sum(), 0.166), yy[~inside], s=7, marker=">", color=T.INK["muted"], zorder=1, linewidth=0)   # off-scale studies (> 16.8 %) drawn as arrows at the edge
    ax.axvspan(0.03, 0.08, color=T.INK["grid"], lw=0, zorder=0); ax.axvline(0.05, color=T.INK["axis"], lw=0.6)
    ax.set_yticks(np.arange(len(layers))); ax.set_yticklabels([]); ax.set_ylim(len(layers) - 0.4, -0.9); ax.set(**XT)
    ax.set_xlabel(f"median share per study ({n} studies, each left out)"); ax.set_title(title, loc="left")
    ax.text(0.166, -0.55, "beyond axis", fontsize=T.PT_MIN, color=T.INK["muted"], va="bottom", ha="right")


def make(out_dir):
    T.print_profile()
    rep = json.load(open(f"{P}/results/s6/report_{_VTAG}_lenient_base.json"))
    ps = pd.DataFrame(rep["per_study"]).rename(columns={"study_bioproject": "study", "n_samples": "n"})
    la = pd.read_csv(f"{P}/results/s6/loso_within_pool.tsv", sep="\t")
    hb = pd.read_csv(f"{P}/results/s11/pipelineB_heldout.tsv", sep="\t"); lb = pd.read_csv(f"{P}/results/s11/pipelineB_loso_all.tsv", sep="\t")
    fig = plt.figure(figsize=(T.DOUBLE_IN, 4.3))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.35, 1], wspace=0.12, hspace=0.75, left=0.1, right=0.975, top=0.905, bottom=0.1)
    heldout(fig.add_subplot(gs[0, 0]), ps, LAYERS, "frac_outside_mean", "a  assembly pipeline:\nindependent healthy cohorts")
    loso_panel(fig.add_subplot(gs[0, 1]), la, LAYERS, "frac", "b  assembly pipeline:\neach baseline study left out")
    heldout(fig.add_subplot(gs[1, 0]), hb, LAYERS_B, "mean", "c  read pipeline:\nindependent healthy cohorts")
    loso_panel(fig.add_subplot(gs[1, 1]), lb, LAYERS_B, "median", "d  read pipeline:\neach baseline study left out")
    T.save(fig, "fig4_calibration", out_dir); plt.close(fig)


if __name__ == "__main__":
    make(sys.argv[1] if len(sys.argv) > 1 else f"{P}/figures")
