"""Figure 14 — what the user gets at sample level (double column): the two sample-level numbers of a report, for the two
example samples of the paper, each against the distribution of the same number in healthy adults scored out of sample.
Rows: (a, b) the C. difficile infection case of Fig. 5 (assembly pipeline); (c, d) the colorectal cancer case of Fig. 13
(read pipeline). Columns: share of assessed features outside the healthy range; number of expected-but-missing features.
Healthy distribution = healthy adults scored against a baseline built without their study (pool studies left out in turn,
plus the held-out healthy cohorts; work/s12/<pipeline>/oos_summaries.parquet): grey bar p10–p90, tick = median, light dots =
individual healthy samples (thinned).

Where the two example samples' numbers come from (the footnote of the figure states it). The C. difficile case is read from
scripts/viz/example_case.layer_counts(), the work/s8 scores table that also feeds Fig. 5, Fig. 13 and
results/s14/example_reports.tsv (179 expected-but-missing KOs on the eggNOG layer); work/s12/A/oos_summaries.parquet holds a
second scoring of the same sample with slightly different counts (184) and is used here only for the healthy distribution.
The colorectal cancer case exists only in the work/s12/B scoring, which is also the source of its counts in
example_reports.tsv (88 missing KOs). Value labels sit to the right of whichever is further right, the sample's diamond or
the healthy median tick, on a white box, so they never cross the tick.

Takeaway (read from the rendered data, not assumed): the per-sample share outside is a weak signal on its own, because
healthy samples scored out of sample spread widely (p90 about 20 % on taxonomy); the C. difficile case is above the
healthy median on both taxonomy layers but inside the healthy p10–p90, while its count of expected-but-missing genes
(179 KOs, 376 Pfam families) is beyond the healthy 90th percentile. The colorectal cancer case is inside the healthy
spread on every layer and on both numbers, and the report says so."""
import os, sys
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import tokens as T, example_case as EC

P = T.P
A_CASE, B_CASE = EC.CASE, "FengQ_2015_SID31223"
LAYERS = {"A": [("taxonomy_family", "family"), ("taxonomy_genus", "genus"), ("ko_eggnog", "KO (eggNOG)"), ("ko_kofam", "KO (KOfam)"), ("pfam", "Pfam"), ("module", "KEGG module")],
          "B": [("taxonomy_family", "family"), ("taxonomy_genus", "genus"), ("taxonomy_species", "species"), ("ko_humann", "KO (HUMAnN)"), ("pathway_humann", "pathway")]}
HEALTHY_ROLES = {"reference_pool", "heldout_healthy"}
LABEL_GAP = 0.03   # value label offset, as a share of the axis range, right of max(diamond, median tick)


def load(pl):
    S = pd.read_parquet(f"{P}/work/s12/{pl}/oos_summaries.parquet"); M = pd.read_parquet(f"{P}/work/s12/{pl}/meta.parquet")
    S = S.merge(M[["sample", "is_healthy", "role", "study"]], on="sample", how="left")
    H = S[S["is_healthy"].fillna(False) & S["role"].isin(HEALTHY_ROLES)]
    return S, H


def case_table(pl, S, case_id):
    """The example sample's per-layer numbers: the C. difficile case from example_case.layer_counts() (the scoring behind
    Fig. 5, Fig. 13 and results/s14/example_reports.tsv), the colorectal cancer case from the out-of-sample scoring table."""
    if pl == "A":
        return EC.layer_counts(case_id).reset_index()[["layer", "frac_outside_raw", "n_missing"]]
    return S.loc[S["sample"] == case_id, ["layer", "frac_outside_raw", "n_missing"]]


def panel(ax, H, case, layers, col, xlab, title, xlim, fmt, seed):
    rng = np.random.default_rng(seed); y = np.arange(len(layers)); span = xlim[1] - xlim[0]
    for i, (layer, _) in enumerate(layers):
        v = H.loc[H["layer"] == layer, col].dropna().to_numpy()
        if not len(v):
            continue
        p10, p50, p90 = np.percentile(v, [10, 50, 90])
        thin = rng.choice(v, size=min(len(v), 160), replace=False)
        ax.scatter(np.clip(thin, xlim[0], xlim[1]), i + rng.uniform(-0.22, 0.22, len(thin)), s=2, color=T.INK["axis"], linewidth=0, zorder=1)
        ax.plot([p10, p90], [i, i], color=T.INK["muted"], lw=3, solid_capstyle="butt", zorder=2); ax.plot([p50, p50], [i - 0.2, i + 0.2], color=T.INK["primary"], lw=0.8, zorder=3)
        c = case.loc[case["layer"] == layer, col]
        if len(c):
            cv = float(c.iloc[0]); out = cv > p90; x = min(cv, xlim[1])
            ax.scatter([x], [i], s=26, color=T.DEV["high_deep"] if out else T.DEV["neutral"], edgecolor=T.INK["surface"], linewidth=0.4, zorder=4, marker="D")
            xl = max(x, p50) + LABEL_GAP * span; ha = "left"
            if xl > xlim[1] - 0.06 * span:   # no room on the right: label left of the diamond and tick instead
                xl = min(x, p50) - LABEL_GAP * span; ha = "right"
            ax.text(xl, i - 0.02, fmt(cv), va="center", ha=ha, fontsize=T.PT_MIN, color=T.INK["primary"], zorder=6,
                    bbox=dict(facecolor=T.INK["surface"], edgecolor="none", pad=0.6))
    ax.set_yticks(y); ax.set_yticklabels([l for _, l in layers]); ax.set_ylim(len(layers) - 0.5, -0.6); ax.tick_params(axis="y", length=0); ax.tick_params(axis="x", length=2, width=0.4)
    ax.set_xlim(*xlim); ax.set_xlabel(xlab); ax.set_title(title, loc="left"); ax.grid(axis="x", zorder=0); ax.set_axisbelow(True)


def make(out_dir):
    T.print_profile()
    fig = plt.figure(figsize=(T.DOUBLE_IN, 4.1))
    gs = fig.add_gridspec(2, 2, width_ratios=[1, 1], wspace=0.45, hspace=0.6, left=0.11, right=0.985, top=0.93, bottom=0.2)
    pct = lambda v: f"{100 * v:.0f} %"; cnt = lambda v: f"{int(v)}"
    n_h = {}
    for r, (pl, case_id, label) in enumerate([("A", A_CASE, "C. difficile infection case, assembly pipeline"), ("B", B_CASE, "colorectal cancer case, read pipeline")]):
        S, H = load(pl); case = case_table(pl, S, case_id); n_h[pl] = (H["sample"].nunique(), H["study"].nunique())
        la, lb = ("a", "b") if r == 0 else ("c", "d")
        panel(fig.add_subplot(gs[r, 0]), H, case, LAYERS[pl], "frac_outside_raw", "share of assessed features outside the healthy range", f"{la}  {label}", (0, 0.42), pct, 1 + r)
        panel(fig.add_subplot(gs[r, 1]), H, case, LAYERS[pl], "n_missing", "expected-but-missing features (count)", f"{lb}  same sample", (0, 420), cnt, 11 + r)
    for ax in fig.axes[::2]:
        ax.set_xticks([0, 0.1, 0.2, 0.3, 0.4]); ax.set_xticklabels(["0", "10 %", "20 %", "30 %", "40 %"])
    T.save(fig, "fig14_sample_card", out_dir); plt.close(fig)


if __name__ == "__main__":
    make(sys.argv[1] if len(sys.argv) > 1 else f"{P}/figures")
