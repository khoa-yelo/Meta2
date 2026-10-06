"""Figure 7 — the healthy map (single column). The healthy baseline samples as a hexbin density on the first two principal
components of family CLR abundance, with contours enclosing the densest 50 % and 95 % of baseline samples; one independent
healthy cohort and one disease cohort projected as points on the same axes.

Coordinates. The baseline cloud and its contours are the coordinates stored in the bundle
(refs/<bundle>/landscape/reference_coords.parquet, PC1/PC2 of the 1,941 pool samples as written by the builder; the pool
is not re-scored). The two cohorts were scored with the packaged checkgm 0.8.1 scorer (Meta2/checkgm/src/checkgm/score.py, which
fixed the 0.8.0 projection bug: undetected basis families are now filled the way the builder filled them, so an in-pool
sample re-projects onto its stored coordinate to within 1e-6) into figures/_data/landscape_v081_<study>.parquet, one row per
sample with PC1-3, core distance and the family-layer summary. The earlier draft of this figure drew results/s7
(0.8.0 scorer) coordinates for all three groups; those are not used here.

Richness. stats() recomputes the Spearman correlation of PC1 and PC2 with the number of assessed families over the
re-scored cohort samples only (the stored reference coordinates carry no sample id to join on). The axis label states the
correlation only when |rho| >= RHO_LABEL both over the pooled cohort samples and within each cohort (a pooled correlation can be
carried by the offset between two cohorts); otherwise the axes read PC1 / PC2. stats() also reports the share of each cohort
outside the baseline's 95 % density contour (KDE on the stored coordinates), quoted in the caption."""
import json, os, sys
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from scipy.stats import spearmanr, gaussian_kde
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import tokens as T
sys.path.insert(0, T.P); from checkgm.paths import bundle_A as _bundle_A

P = T.P
BUNDLE = _bundle_A("lenient")
HEALTHY_Q = "PRJEB62833"   # independent healthy adults (never in the baseline)
DISEASE_Q = "PRJEB26165"   # C. difficile infection cases
DATA = f"{P}/figures/_data"
RHO_LABEL = 0.7            # |Spearman rho| needed before the axis label mentions richness


def density(x, y, masses=(0.5, 0.95), grid=160):
    """KDE of the baseline coordinates on a grid, the contour level enclosing each share of the baseline samples (so a 95 %
    contour has 5 % of baseline samples outside it by construction), and the KDE itself for point-wise evaluation."""
    k = gaussian_kde(np.vstack([x, y]))
    xs = np.linspace(x.min() - 3, x.max() + 3, grid); ys = np.linspace(y.min() - 3, y.max() + 3, grid)
    X, Y = np.meshgrid(xs, ys); Z = k(np.vstack([X.ravel(), Y.ravel()])).reshape(X.shape)
    at_pts = k(np.vstack([x, y]))   # level enclosing a mass m = the density below which the (1 - m) least dense baseline samples lie
    levels = {m: float(np.quantile(at_pts, 1 - m)) for m in masses}
    return X, Y, Z, levels, k


def load():
    base = pd.read_parquet(f"{BUNDLE}/landscape/reference_coords.parquet")
    hq = pd.read_parquet(f"{DATA}/landscape_v081_{HEALTHY_Q}.parquet")
    dq = pd.read_parquet(f"{DATA}/landscape_v081_{DISEASE_Q}.parquet")
    return base, hq, dq


def stats():
    """Numbers quoted in the caption: Spearman rho of PC1/PC2 with assessed families over the re-scored cohort samples
    (both cohorts pooled, and per cohort), cohort sizes, and the share of each cohort outside the 95 % contour."""
    base, hq, dq = load(); q = pd.concat([hq, dq], ignore_index=True)
    _, _, _, lev, k = density(base["PC1"].to_numpy(), base["PC2"].to_numpy())
    outside = lambda d: float((k(np.vstack([d["PC1"], d["PC2"]])) < lev[0.95]).mean())
    rho = lambda d, pc: float(spearmanr(d[pc], d["n_assessed"]).statistic)
    return {"rho_PC1_n_assessed": rho(q, "PC1"), "rho_PC2_n_assessed": rho(q, "PC2"),
            "rho_PC1_heldout": rho(hq, "PC1"), "rho_PC1_cdi": rho(dq, "PC1"), "rho_PC2_heldout": rho(hq, "PC2"), "rho_PC2_cdi": rho(dq, "PC2"),
            "n_base": int(len(base)), "n_heldout": int(len(hq)), "n_cdi": int(len(dq)),
            "frac_outside95_heldout": outside(hq), "frac_outside95_cdi": outside(dq), "frac_outside95_base": outside(base),
            "median_core_pct_heldout": float(hq["core_distance_pct"].median()), "median_core_pct_cdi": float(dq["core_distance_pct"].median())}


def make(out_dir):
    T.print_profile()
    ve = json.load(open(f"{BUNDLE}/landscape/core_distance.json"))["variance_explained"]
    base, hq, dq = load(); st = stats(); rho1, rho2 = st["rho_PC1_n_assessed"], st["rho_PC2_n_assessed"]
    fig, ax = plt.subplots(figsize=(T.SINGLE_IN, 3.3))
    cmap = ListedColormap(T.SEQ); norm = BoundaryNorm([1, 2, 4, 8, 16, 32], cmap.N, extend="max")
    hb = ax.hexbin(base["PC1"], base["PC2"], gridsize=26, cmap=cmap, norm=norm, mincnt=1, linewidths=0.15, edgecolors=T.INK["surface"], zorder=1)
    X, Y, Z, lev, _ = density(base["PC1"].to_numpy(), base["PC2"].to_numpy())
    ax.contour(X, Y, Z, levels=[lev[0.95], lev[0.5]], colors=[T.INK["primary"]], linewidths=[0.5, 0.8], linestyles=["dashed", "solid"], zorder=2)
    ax.scatter(hq["PC1"], hq["PC2"], s=7, facecolor="none", edgecolor=T.INK["primary"], linewidth=0.5, zorder=3, label=f"independent healthy cohort ({HEALTHY_Q}, n = {len(hq)})")
    ax.scatter(dq["PC1"], dq["PC2"], s=9, marker="D", facecolor=T.INK["primary"], edgecolor=T.INK["surface"], linewidth=0.3, zorder=4, label=f"C. difficile infection cases ({DISEASE_Q}, n = {len(dq)})")
    ax.plot([], [], color=T.INK["primary"], lw=0.8, label="densest 50 % of baseline samples"); ax.plot([], [], color=T.INK["primary"], lw=0.5, ls="--", label="densest 95 % of baseline samples")
    cb = fig.colorbar(hb, ax=ax, fraction=0.045, pad=0.02, ticks=[1, 2, 4, 8, 16, 32]); cb.ax.tick_params(labelsize=T.PT_MIN, length=2, width=0.4); cb.outline.set_linewidth(0.4)
    cb.set_label(f"healthy baseline samples per hex (n = {len(base):,})", fontsize=T.PT_MIN, color=T.INK["secondary"])
    rich1 = min(abs(rho1), abs(st["rho_PC1_heldout"]), abs(st["rho_PC1_cdi"])) >= RHO_LABEL   # pooled and within each cohort
    rich2 = min(abs(rho2), abs(st["rho_PC2_heldout"]), abs(st["rho_PC2_cdi"])) >= RHO_LABEL
    x_rich = f"; Spearman ρ = {rho1:.2f} with assessed families" if rich1 else ""
    y_rich = f"; ρ = {rho2:.2f}" if rich2 else ""
    ax.set_xlabel(f"PC1 ({ve[0]*100:.1f} % of variance{x_rich})")
    ax.set_ylabel(f"PC2 ({ve[1]*100:.1f} %{y_rich})")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=1, frameon=False, handletextpad=0.3, borderaxespad=0, markerscale=1.2, labelspacing=0.3)
    ax.tick_params(length=2, width=0.4)
    fig.subplots_adjust(left=0.15, right=0.86, top=0.97, bottom=0.34)
    T.save(fig, "fig7_landscape", out_dir); plt.close(fig)
    return st


if __name__ == "__main__":
    print(make(sys.argv[1] if len(sys.argv) > 1 else f"{P}/figures"))
