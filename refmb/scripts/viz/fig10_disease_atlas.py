"""Figure 10 — disease atlas (double column). Shift of the median reference percentile, cases minus same-study controls,
for genera, per case/control study. a: pipeline A, b: pipeline B. Five-step diverging scale (tokens); a dot marks
q <= 0.05. Genera are those with the most studies significant, ties by absolute median shift; studies grouped by condition."""
import os, sys
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import tokens as T

P = T.P
BOUNDS = [-100, -20, -5, 5, 20, 100]; STEP_LABELS = ["≤ −20", "−20 to −5", "−5 to 5", "5 to 20", "≥ 20"]
SHORT = {"Inflammatory bowel disease": "IBD", "Parkinson's disease": "PD", "Rheumatoid arthritis": "RA", "Colorectal neoplasia": "CRC", "Seasonal allergies": "allergy",
         "Hemodialysis": "dialysis", "Hypertension": "hypert.", "Kidney stones": "k. stones", "C. difficile infection": "CDI", "schizofrenia": "SCZ"}
NTOP = 22


def matrix(pl, names):
    F = pd.read_csv(f"{P}/results/s12/disease_atlas_{pl}_features.tsv", sep="\t", dtype={"feature_id": str})
    F = F[F["layer"] == "taxonomy_genus"].copy(); F["sig"] = (F["q"] <= 0.05) & (F["shift"].abs() >= 10)
    F["cond"] = F["condition"].map(lambda c: SHORT.get(c, c))
    rank = F.groupby("feature_id").agg(ns=("sig", "sum"), ms=("shift", lambda s: s.abs().median()), n=("shift", "size"))
    rank = rank[rank["n"] >= 2].sort_values(["ns", "ms"], ascending=False).head(NTOP)
    cols = F.drop_duplicates("study").sort_values(["cond", "study"])[["study", "cond"]]
    M = F.pivot_table(index="feature_id", columns="study", values="shift").reindex(index=rank.index, columns=cols["study"])
    S = F.pivot_table(index="feature_id", columns="study", values="sig", aggfunc="max").reindex(index=rank.index, columns=cols["study"]).fillna(False).astype(bool)
    order = M.mean(axis=1).sort_values().index
    return M.loc[order], S.loc[order], cols, [names.get(i, i) for i in order]


def make(out_dir):
    T.print_profile()
    tax = pd.read_parquet(f"{P}/resources/backbone/ncbi_taxonomy.parquet", columns=["taxid", "name"]); names = dict(zip(tax["taxid"].astype(str), tax["name"]))
    cmap = ListedColormap([T.DEV[k] for k in T.DEV_ORDER]); norm = BoundaryNorm(BOUNDS, cmap.N)
    fig = plt.figure(figsize=(T.DOUBLE_IN, 4.4))
    gs = fig.add_gridspec(1, 2, width_ratios=[10, 15], wspace=0.5, left=0.19, right=0.995, top=0.95, bottom=0.33)
    for k, pl in enumerate("AB"):
        M, S, cols, rl = matrix(pl, names); ax = fig.add_subplot(gs[k])
        ax.pcolormesh(np.arange(M.shape[1] + 1), np.arange(M.shape[0] + 1), np.ma.masked_invalid(M.values), cmap=cmap, norm=norm, edgecolors=T.INK["surface"], linewidth=0.8)
        yy, xx = np.where(S.values); deep = np.abs(M.values[yy, xx]) >= 20   # white dot on the two deep steps, ink on the light ones
        ax.scatter(xx[deep] + 0.5, yy[deep] + 0.5, s=2.2, color=T.INK["surface"], linewidth=0, zorder=3)
        ax.scatter(xx[~deep] + 0.5, yy[~deep] + 0.5, s=2.2, color=T.INK["primary"], linewidth=0, zorder=3)
        ax.set_yticks(np.arange(len(rl)) + 0.5); ax.set_yticklabels(rl, style="italic", color=T.INK["secondary"])
        ax.set_xticks(np.arange(M.shape[1]) + 0.5); ax.set_xticklabels([f"{c}  {s}" for s, c in zip(cols["study"], cols["cond"])], rotation=90, color=T.INK["secondary"])
        ax.tick_params(length=0); ax.set_xlim(0, M.shape[1]); ax.set_ylim(0, M.shape[0])
        for sp in ax.spines.values():
            sp.set_visible(False)
        cv = cols["cond"].astype(str).to_numpy(); b = [k for k in range(1, len(cv)) if cv[k] != cv[k - 1]]   # one separator per condition group
        for x in b:
            ax.axvline(x, color=T.INK["secondary"], lw=0.6, zorder=4)
        ax.set_title({"A": "a  assembly pipeline (A), genera", "B": "b  read pipeline (B), genera"}[pl], loc="left")
    h = [Patch(facecolor=T.DEV[k], edgecolor="none", label=l) for k, l in zip(T.DEV_ORDER, STEP_LABELS)]
    h.append(plt.Line2D([], [], marker="o", ls="", ms=1.8, color=T.INK["primary"], label="q ≤ 0.05 and |shift| ≥ 10"))
    fig.legend(handles=h, loc="lower center", bbox_to_anchor=(0.5, 0.0), ncol=6, frameon=False, handlelength=1.2, handletextpad=0.4, columnspacing=1.2, borderaxespad=0,
               title="shift in median percentile, cases minus same-study controls (percentile points); blank = genus not assessed in the study", title_fontsize=T.PT_MIN)
    T.save(fig, "fig10_disease_atlas", out_dir); plt.close(fig)


if __name__ == "__main__":
    make(sys.argv[1] if len(sys.argv) > 1 else f"{P}/figures")
