"""Figure 16 — function features that shift consistently across colorectal cancer cohorts (read pipeline, double column).
Rows: MetaCyc pathways and KEGG orthologs shifted in the same direction (>= 10 points, q <= 0.05) in at least three of the
nine CRC cohorts, plus the five strongest shifted in two (results/s12/disease_atlas_B_consistent.tsv); cells: shift in median
percentile, cases minus same-study controls (results/s12/disease_atlas_B_features.tsv); right margin: cohorts agreeing."""
import os, sys, glob, gzip
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import tokens as T

P = T.P
BOUNDS = [-100, -20, -5, 5, 20, 100]; STEP_LABELS = ["≤ −20", "−20 to −5", "−5 to 5", "5 to 20", "≥ 20"]


def names():
    nm = {}
    for f in sorted(glob.glob(f"{P}/staging/cmd3_humann/pathway_abundance/*.parquet"))[:8]:
        for x in pd.read_parquet(f, columns=["feature_id"]).feature_id.unique():
            if ":" in x:
                nm[x.split(":")[0].strip()] = x.split(":", 1)[1].strip().replace("&gamma;", "γ").replace("&beta;", "β").replace("&alpha;", "α")
    with gzip.open(f"{P}/staging/cmd3_humann/mapping/map_ko_name.txt.gz", "rt") as fh:
        for line in fh:
            k, v = line.rstrip("\n").split("\t", 1); nm[k] = v.split(" [EC")[0]
    return nm


def rows():
    C = pd.read_csv(f"{P}/results/s12/disease_atlas_B_consistent.tsv", sep="\t", dtype={"feature_id": str})
    C = C[(C.condition == "CRC") & C.layer.isin(["pathway_humann", "ko_humann"])]
    top = pd.concat([C[C.n_studies >= 3], C[C.n_studies == 2].assign(a=lambda d: d.median_shift.abs()).sort_values("a", ascending=False).head(5)])
    top = top.assign(a=top.median_shift.abs()); return top.sort_values(["layer", "n_studies", "a"], ascending=[False, False, False])


def make(out_dir):
    T.print_profile(); R = rows(); nm = names()
    F = pd.read_csv(f"{P}/results/s12/disease_atlas_B_features.tsv", sep="\t", dtype={"feature_id": str}); F = F[F.condition == "CRC"]
    F["sig"] = (F["q"] <= 0.05) & (F["shift"].abs() >= 10); studies = sorted(F.study.unique())
    key = list(zip(R.layer, R.feature_id)); F = F.set_index(["layer", "feature_id"]).sort_index()
    M = np.full((len(key), len(studies)), np.nan); S = np.zeros_like(M, dtype=bool)
    for i, k in enumerate(key):
        if k in F.index:
            g = F.loc[[k]].set_index("study")
            for j, st in enumerate(studies):
                if st in g.index:
                    M[i, j] = g.loc[st, "shift"]; S[i, j] = bool(g.loc[st, "sig"])
    lab = []
    for (layer, fid) in key:
        n = nm.get(fid, fid); n = n if len(n) <= 52 else n[:50].rstrip() + "…"
        lab.append(f"{n}  ({fid})" if layer == "ko_humann" else n)
    cmap = ListedColormap([T.DEV[k] for k in T.DEV_ORDER]); norm = BoundaryNorm(BOUNDS, cmap.N)
    fig = plt.figure(figsize=(T.DOUBLE_IN, 0.6 + 0.115 * len(key) + 1.3))
    ax = fig.add_axes([0.44, 0.24, 0.42, 0.69])
    ax.pcolormesh(np.arange(len(studies) + 1), np.arange(len(key) + 1), np.ma.masked_invalid(M), cmap=cmap, norm=norm, edgecolors=T.INK["surface"], linewidth=0.8)
    yy, xx = np.where(S); deep = np.abs(M[yy, xx]) >= 20
    ax.scatter(xx[deep] + 0.5, yy[deep] + 0.5, s=2.2, color=T.INK["surface"], linewidth=0, zorder=3)
    ax.scatter(xx[~deep] + 0.5, yy[~deep] + 0.5, s=2.2, color=T.INK["primary"], linewidth=0, zorder=3)
    ax.set_yticks(np.arange(len(key)) + 0.5); ax.set_yticklabels(lab, color=T.INK["secondary"])
    ax.set_xticks(np.arange(len(studies)) + 0.5); ax.set_xticklabels(studies, rotation=90, color=T.INK["secondary"])
    ax.tick_params(length=0); ax.set_xlim(0, len(studies)); ax.set_ylim(len(key), 0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    nko = int((R.layer == "ko_humann").sum())
    ax.axhline(len(key) - nko, color=T.INK["secondary"], lw=0.6)
    for i, (n, d) in enumerate(zip(R.n_studies, R.consistent)):
        ax.text(len(studies) + 0.3, i + 0.5, f"   {n} of 9 {'lower' if d == 'down' else 'higher'}", fontsize=T.PT_MIN, va="center", color=T.INK["secondary"])
    fig.text(0.02, 0.985, "MetaCyc pathways above the line, KEGG orthologs below it (KO identifier in brackets)", ha="left", va="top", fontsize=T.PT_BODY, color=T.INK["primary"])
    ax.text(len(studies) + 0.3, -0.25, "   cohorts shifted", fontsize=T.PT_MIN, color=T.INK["muted"], va="bottom")
    h = [Patch(facecolor=T.DEV[k], edgecolor="none", label=l) for k, l in zip(T.DEV_ORDER, STEP_LABELS)]
    h.append(plt.Line2D([], [], marker="o", ls="", ms=1.8, color=T.INK["primary"], label="q ≤ 0.05 and |shift| ≥ 10"))
    fig.legend(handles=h, loc="lower center", bbox_to_anchor=(0.5, 0.0), ncol=6, frameon=False, handlelength=1.2, handletextpad=0.4, columnspacing=1.2, borderaxespad=0,
               title="shift in median percentile, colorectal cancer cases minus same-study controls (percentile points); blank = not assessed", title_fontsize=T.PT_MIN)
    T.save(fig, "fig16_function_shifts", out_dir); plt.close(fig)


if __name__ == "__main__":
    make(sys.argv[1] if len(sys.argv) > 1 else f"{P}/figures")
