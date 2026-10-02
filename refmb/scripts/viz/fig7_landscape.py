"""Figure 7 — the healthy map (single column). Healthy baseline samples as a hexbin density on the first two principal
components of family CLR abundance, with the baseline's 50 % and 95 % density contours; one independent healthy cohort and
one disease cohort projected as points. All coordinates, baseline and cohorts alike, are the ones the scoring code writes
into a report (results/s7/summaries_gut.parquet, PC1/PC2), so the three groups are on one footing: the bundle's own
reference_coords.parquet is on a different scale from the scored projection (see the viz round-2 report) and is not drawn.
PC1 tracks richness: its Spearman correlation with the number of assessed families is printed in the axis label."""
import json, os, sys
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from scipy.stats import spearmanr, gaussian_kde
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import tokens as T
sys.path.insert(0, T.P); from refmb.paths import bundle_A as _bundle_A

P = T.P
BUNDLE = _bundle_A("lenient")
HEALTHY_Q = "PRJEB62833"   # independent healthy adults (never in the baseline)
DISEASE_Q = "PRJEB26165"   # C. difficile infection cases


def density_levels(x, y, masses=(0.5, 0.95), grid=160):
    """KDE on a grid; the contour level enclosing each mass fraction."""
    k = gaussian_kde(np.vstack([x, y]))
    xs = np.linspace(x.min() - 3, x.max() + 3, grid); ys = np.linspace(y.min() - 3, y.max() + 3, grid)
    X, Y = np.meshgrid(xs, ys); Z = k(np.vstack([X.ravel(), Y.ravel()])).reshape(X.shape)
    z = np.sort(Z.ravel())[::-1]; cum = np.cumsum(z) / z.sum()
    levels = [float(z[np.searchsorted(cum, m)]) for m in masses]
    return X, Y, Z, sorted(levels)


def load():
    S = pd.read_parquet(f"{P}/results/s7/summaries_gut.parquet")
    S = S[S["layer"] == "taxonomy_family"].drop_duplicates("analysis_id")
    base = S[S["in_reference_pool"] == True]
    hq = S[S["study_bioproject"] == HEALTHY_Q]
    dq = S[(S["study_bioproject"] == DISEASE_Q) & (S["cur_Health_status_group"] == "Diseased")]
    return base, hq, dq


def stats():
    """Spearman correlation of PC1 and PC2 with the number of assessed families, over the baseline samples (for captions)."""
    base, hq, dq = load()
    return {"rho_PC1_n_assessed": float(spearmanr(base["PC1"], base["n_assessed"]).statistic), "rho_PC2_n_assessed": float(spearmanr(base["PC2"], base["n_assessed"]).statistic),
            "n_base": int(len(base)), "n_heldout": int(len(hq)), "n_cdi": int(len(dq))}


def make(out_dir):
    T.print_profile()
    ve = json.load(open(f"{BUNDLE}/landscape/core_distance.json"))["variance_explained"]
    base, hq, dq = load(); st = stats(); rho1, rho2 = st["rho_PC1_n_assessed"], st["rho_PC2_n_assessed"]
    fig, ax = plt.subplots(figsize=(T.SINGLE_IN, 3.6))
    cmap = ListedColormap(T.SEQ); norm = BoundaryNorm([1, 2, 4, 8, 16, 32], cmap.N, extend="max")
    hb = ax.hexbin(base["PC1"], base["PC2"], gridsize=26, cmap=cmap, norm=norm, mincnt=1, linewidths=0.15, edgecolors=T.INK["surface"], zorder=1)
    X, Y, Z, lev = density_levels(base["PC1"].to_numpy(), base["PC2"].to_numpy())
    cs = ax.contour(X, Y, Z, levels=lev, colors=[T.INK["primary"]], linewidths=[0.5, 0.8], linestyles=["dashed", "solid"], zorder=2)
    ax.scatter(hq["PC1"], hq["PC2"], s=7, facecolor="none", edgecolor=T.INK["primary"], linewidth=0.5, zorder=3, label=f"independent healthy cohort ({HEALTHY_Q}, n = {len(hq)})")
    ax.scatter(dq["PC1"], dq["PC2"], s=9, marker="D", facecolor=T.INK["primary"], edgecolor=T.INK["surface"], linewidth=0.3, zorder=4, label=f"C. difficile infection cases ({DISEASE_Q}, n = {len(dq)})")
    ax.plot([], [], color=T.INK["primary"], lw=0.8, label="densest 50 % of baseline samples"); ax.plot([], [], color=T.INK["primary"], lw=0.5, ls="--", label="densest 95 % of baseline samples")
    cb = fig.colorbar(hb, ax=ax, fraction=0.045, pad=0.02, ticks=[1, 2, 4, 8, 16, 32]); cb.ax.tick_params(labelsize=T.PT_MIN, length=2, width=0.4); cb.outline.set_linewidth(0.4)
    cb.set_label(f"healthy baseline samples per hex (n = {len(base):,})", fontsize=T.PT_MIN, color=T.INK["secondary"])
    ax.set_xlabel(f"PC1 ({ve[0]*100:.1f} % of variance; Spearman ρ = {rho1:.2f} with assessed families)")
    ax.set_ylabel(f"PC2 ({ve[1]*100:.1f} %; ρ = {rho2:.2f})")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=1, frameon=False, handletextpad=0.3, borderaxespad=0, markerscale=1.2, labelspacing=0.3)
    ax.tick_params(length=2, width=0.4)
    ax.set_title("Position on the map mostly reflects richness; the independent\ncohort lies in the dense half of the healthy cloud, the\nC. difficile cases at its edge", loc="left", fontsize=T.PT_BODY)
    fig.subplots_adjust(left=0.15, right=0.86, top=0.85, bottom=0.33)
    T.save(fig, "fig7_landscape", out_dir); plt.close(fig)
    return st


if __name__ == "__main__":
    print(make(sys.argv[1] if len(sys.argv) > 1 else f"{P}/figures"))
