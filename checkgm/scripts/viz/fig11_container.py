"""Figure 11 — container fidelity (double column), default-profile anchors. a: assembled length, container vs native.
b: container-vs-native percentile Spearman per layer (one point per anchor, median bar). c: family and KO Spearman by
native assembler group."""
import os, sys
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import tokens as T

P = T.P
LAYERS = [("ko_eggnog", "KO\neggNOG"), ("ko_kofam", "KO\nKOfam"), ("pfam", "Pfam"), ("module", "module"), ("taxonomy_family", "family"), ("taxonomy_genus", "genus")]
ERAS = ["single-end\nMEGAHIT", "paired\n3.15.3", "paired\n3.12–3.14", "paired\nversion\nunknown"]


def load():
    C = pd.read_csv(f"{P}/results/s11/tierA_anchor_compare.tsv", sep="\t"); C = C[C["variant"] == "default"].copy()
    Q = pd.read_csv(f"{P}/results/s11/tierA_anchor_qc.tsv", sep="\t"); Q = Q[Q["variant"] == "default"].copy()
    A = pd.read_csv(f"{P}/work/s11/anchor_panel.tsv", sep="\t").set_index("anchor_id"); E = pd.read_csv(f"{P}/work/s11/anchor_reads.tsv", sep="\t").set_index("anchor_id")
    inv = pd.read_parquet(f"{P}/work/s0/inventory.parquet", columns=["analysis_id", "assembler"]).set_index("analysis_id")["assembler"].fillna("")
    for D in (C, Q):
        D["base"] = D["anchor"].str.replace(r"\+.*", "", regex=True)
        lay = D["base"].map(E["library_layout"]); asm = D["base"].map(A["analysis_id"]).map(inv).str.lower()
        D["era"] = np.where(lay == "SINGLE", ERAS[0], np.where(asm.str.contains("3.15", regex=False), ERAS[1], np.where(asm.eq(""), ERAS[3], ERAS[2])))
    return C, Q


def jitter(n, seed):
    return np.random.default_rng(seed).uniform(-0.22, 0.22, n)


def make(out_dir):
    T.print_profile(); C, Q = load()
    fig = plt.figure(figsize=(T.DOUBLE_IN, 3.0))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.85, 1.3], wspace=0.4, left=0.075, right=0.965, top=0.9, bottom=0.32)
    # a assembled length
    ax = fig.add_subplot(gs[0]); single = Q["era"] == ERAS[0]
    lo, hi = 5, 2000
    ax.plot([lo, hi], [lo, hi], color=T.INK["axis"], lw=0.6, zorder=1)
    for m, mk, lab in [(~single, T.LAYOUT_MARK["paired"], "paired (metaSPAdes)"), (single, T.LAYOUT_MARK["single"], "single-end (MEGAHIT)")]:
        ax.scatter(Q.loc[m, "assembled_Mb_native"], Q.loc[m, "assembled_Mb_container"], s=10, marker=mk, color=T.INK["primary"], edgecolor=T.INK["surface"], linewidth=0.3, zorder=3, label=lab)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(lo, hi); ax.set_ylim(lo, hi); ax.set_aspect("equal"); ax.set_anchor("N")
    ax.set_xticks([10, 100, 1000]); ax.set_yticks([10, 100, 1000]); ax.set_xticklabels(["10", "100", "1,000"]); ax.set_yticklabels(["10", "100", "1,000"]); ax.minorticks_off()
    r = Q["assembled_Mb_container"] / Q["assembled_Mb_native"]
    ax.text(0.05, 0.95, f"n = {len(Q)}\n{int(((r - 1).abs() <= 0.02).sum())} within 2 %", transform=ax.transAxes, ha="left", va="top", fontsize=T.PT_MIN, color=T.INK["secondary"])
    ax.set_xlabel("MGnify assembly (Mb)"); ax.set_ylabel("container assembly (Mb)"); ax.set_title("a  assembled length", loc="left")
    ax.legend(loc="upper center", bbox_to_anchor=(0.45, -0.3), ncol=1, frameon=False, handletextpad=0.2, borderaxespad=0)
    # b per layer
    ax = fig.add_subplot(gs[1])
    for i, (l, _) in enumerate(LAYERS):
        v = C.loc[C["layer"] == l, "spearman_pct"].dropna().values
        ax.scatter(i + jitter(len(v), i), v, s=5, color=T.INK["muted"], linewidth=0, zorder=2)
        ax.hlines(np.median(v), i - 0.32, i + 0.32, color=T.INK["primary"], lw=1.0, zorder=3)
        ax.text(i, 1.035, f"{np.median(v):.2f}", ha="center", va="bottom", fontsize=T.PT_MIN, color=T.INK["primary"])
    ax.set_xticks(range(len(LAYERS))); ax.set_xticklabels([n for _, n in LAYERS]); ax.set_xlim(-0.6, len(LAYERS) - 0.4); ax.set_ylim(0.3, 1.1); ax.set_yticks([0.4, 0.6, 0.8, 1.0])
    ax.grid(axis="y", zorder=0); ax.set_axisbelow(True)
    ax.set_ylabel("Spearman correlation of percentiles,\ncontainer vs MGnify"); ax.set_title("b  per layer (bar = median)", loc="left")
    # c by era
    ax = fig.add_subplot(gs[2]); off = {"ko_eggnog": -0.2, "taxonomy_family": 0.2}; mk = {"ko_eggnog": T.LAYER_MARK["KO"], "taxonomy_family": T.LAYER_MARK["family"]}; fc = {"ko_eggnog": T.INK["surface"], "taxonomy_family": T.INK["primary"]}
    for i, e in enumerate(ERAS):
        for l in off:
            v = C.loc[(C["era"] == e) & (C["layer"] == l), "spearman_pct"].dropna().values
            ax.scatter(i + off[l] + jitter(len(v), i) * 0.5, v, s=6, marker=mk[l], facecolor=fc[l], edgecolor=T.INK["primary"], linewidth=0.4, zorder=2)
            if len(v):
                ax.hlines(np.median(v), i + off[l] - 0.17, i + off[l] + 0.17, color=T.INK["primary"], lw=1.0, zorder=3)
        n = C.loc[C["era"] == e, "anchor"].nunique(); ax.text(i, 1.035, f"n = {n}", ha="center", va="bottom", fontsize=T.PT_MIN, color=T.INK["secondary"])
    ax.set_xticks(range(len(ERAS))); ax.set_xticklabels(ERAS); ax.set_ylim(0.3, 1.1); ax.set_yticks([0.4, 0.6, 0.8, 1.0]); ax.grid(axis="y", zorder=0); ax.set_axisbelow(True)
    ax.set_title("c  by MGnify assembler", loc="left"); ax.set_xlabel("read layout; MGnify metaSPAdes version")
    h = [plt.Line2D([], [], marker=mk[l], ls="", ms=3.2, markerfacecolor=fc[l], markeredgecolor=T.INK["primary"], markeredgewidth=0.5, label=n) for l, n in [("ko_eggnog", "KO (eggNOG)"), ("taxonomy_family", "family")]]
    ax.legend(handles=h, loc="upper center", bbox_to_anchor=(0.5, -0.4), ncol=2, frameon=False, handletextpad=0.2, borderaxespad=0)
    T.save(fig, "fig11_container", out_dir); plt.close(fig)


if __name__ == "__main__":
    make(sys.argv[1] if len(sys.argv) > 1 else f"{P}/figures")
