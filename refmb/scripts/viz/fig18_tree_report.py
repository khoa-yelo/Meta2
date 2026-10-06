"""Figure 18 — one sample on a radial taxonomy tree (double column). The C. difficile infection case MGYA00694698 (assembly
pipeline), scored against the baseline rebuilt without its study. Inside: the taxonomy of the 232 families carried by at
least half of healthy adults in some assembly-size class, drawn as a radial dendrogram (phylum, class, order, family).
Outside: a percentile track (inner edge 0, outer edge 100) with the healthy range (2.5th–97.5th) shaded; each family's dot
is the sample's percentile (blue below the range, orange above, grey within; open ring = expected but missing). Outermost
ring: families outside the range in at least 25 % (light) or 50 % (deep) of the study's 56 cases. Labelled: families out of
range both here and in at least 25 % of cases. Data: fig5_tree.load()."""
import os, sys
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge, Arc
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import tokens as T, fig5_tree as F5

P = T.P
R_RANK = {"root": 0.0, "phylum_taxid": 0.32, "class_taxid": 0.52, "order_taxid": 0.72, "leaf": 0.9}
R_IN, R_OUT = 1.0, 1.75


def make(out_dir):
    T.print_profile(); d, ref, ncases, cnt = F5.load()
    d = d.copy(); n = len(d); theta = np.linspace(0, 2 * np.pi, n, endpoint=False) + np.pi / 2 + np.pi / n; d["theta"] = theta
    rp = lambda p: R_IN + (R_OUT - R_IN) * np.clip(p, 0, 100) / 100
    fig = plt.figure(figsize=(T.DOUBLE_IN, 4.6)); ax = fig.add_axes([0.02, 0.01, 0.78, 0.98]); ax.set_aspect("equal"); ax.axis("off")
    ax.set_xlim(-3.0, 2.4); ax.set_ylim(-2.2, 2.2)
    line = dict(color=T.INK["secondary"], lw=0.45, solid_capstyle="round", zorder=2)
    # radial dendrogram: each rank's node spans its leaves' angles
    ranks = ["phylum_taxid", "class_taxid", "order_taxid"]
    def node_angle(sub):
        return float(sub["theta"].mean())
    for i, rk in enumerate(ranks):
        r = R_RANK[rk]; rparent = R_RANK[ranks[i - 1]] if i else 0.0
        for key, sub in d.groupby(ranks[: i + 1], sort=False):
            a = node_angle(sub); t0, t1 = sub["theta"].min(), sub["theta"].max()
            ax.plot([rparent * np.cos(a), r * np.cos(a)], [rparent * np.sin(a), r * np.sin(a)], **line)
            if t1 > t0:
                ax.add_patch(Arc((0, 0), 2 * r, 2 * r, theta1=np.degrees(t0), theta2=np.degrees(t1), color=T.INK["secondary"], lw=0.45, zorder=2))
    for (ph, cl), sub in d.groupby(["phylum_taxid", "class_taxid"], sort=False):
        pass
    # phylum-level arc joining phyla at the root
    a_ph = d.groupby("phylum_taxid", sort=False)["theta"].mean()
    ax.scatter([0], [0], s=6, color=T.INK["secondary"], zorder=3)
    for _, r in d.iterrows():   # leaf stems
        c, s_ = np.cos(r.theta), np.sin(r.theta)
        ax.plot([R_RANK["order_taxid"] * c, R_RANK["leaf"] * c], [R_RANK["order_taxid"] * s_, R_RANK["leaf"] * s_], color=T.INK["axis"], lw=0.35, zorder=1)
    # percentile track
    ax.add_patch(Wedge((0, 0), rp(97.5), 0, 360, width=rp(97.5) - rp(2.5), facecolor=T.INK["grid"], edgecolor="none", zorder=0))
    for p_, lw in ((0, 0.5), (50, 0.4), (100, 0.5)):
        ax.add_patch(plt.Circle((0, 0), rp(p_), fill=False, edgecolor=T.INK["axis"], lw=lw, ls="-" if p_ != 50 else (0, (2, 2)), zorder=1))
    for _, r in d.iterrows():
        c, s_ = np.cos(r.theta), np.sin(r.theta)
        if r.call == "expected_but_missing":
            ax.scatter([rp(0) * c], [rp(0) * s_], s=13, zorder=5, **T.MISSING_RING); continue
        if pd.isna(r.pct):
            continue
        p = r.pct; out = p < 2.5 or p > 97.5
        col = (T.DEV["low_deep"] if p < 1 else T.DEV["low_light"]) if p < 2.5 else ((T.DEV["high_deep"] if p > 99 else T.DEV["high_light"]) if p > 97.5 else T.DEV["neutral"])
        ax.scatter([rp(p) * c], [rp(p) * s_], s=15 if out else 5, color=col, edgecolor=T.INK["surface"] if out else "none", linewidth=0.4, zorder=6 if out else 4)
    # cohort ring
    w = 180 / n
    for _, r in d.iterrows():
        fh, fl = r.cohort_high, r.cohort_low
        if max(fh, fl) >= 0.25:
            col = (T.DEV["high_deep"] if fh >= 0.5 else T.DEV["high_light"]) if fh >= fl else (T.DEV["low_deep"] if fl >= 0.5 else T.DEV["low_light"])
            ax.add_patch(Wedge((0, 0), R_OUT + 0.13, np.degrees(r.theta) - w * 0.9, np.degrees(r.theta) + w * 0.9, width=0.08, facecolor=col, edgecolor="none", zorder=3))
    # phylum labels (curved position, horizontal text) for phyla with >= 8 families, with a thin arc
    for ph, sub in d.groupby("phylum_taxid", sort=False):
        t0, t1 = sub["theta"].min() - np.pi / n, sub["theta"].max() + np.pi / n
        gap = np.pi / n * 0.8
        ax.add_patch(Arc((0, 0), 2 * (R_OUT + 0.24), 2 * (R_OUT + 0.24), theta1=np.degrees(t0 + gap), theta2=np.degrees(t1 - gap), color=T.INK["axis"], lw=1.4, zorder=2))
        if len(sub) >= 8:
            a = (t0 + t1) / 2; rr = R_OUT + 0.4
            ax.text(rr * np.cos(a), rr * np.sin(a), ph, ha="center" if abs(np.cos(a)) < 0.3 else ("left" if np.cos(a) > 0 else "right"), va="center", fontsize=T.PT_MIN, color=T.INK["secondary"])
    # labels: out of range here and in >= 25 % of cases
    lab = d[(d.call.isin(["low", "high"])) & (np.maximum(d.cohort_low, d.cohort_high) >= 0.25)]
    for side in (-1, 1):
        col = lab[(np.cos(lab.theta) >= 0) == (side > 0)].sort_values("theta", key=lambda t: -np.sin(t))
        ys = [2.0 * np.sin(t) for t in col.theta]
        for k in range(1, len(ys)):
            ys[k] = min(ys[k], ys[k - 1] - 0.2)
        for (_, r), yl in zip(col.iterrows(), ys):
            x0, y0 = rp(r.pct) * np.cos(r.theta), rp(r.pct) * np.sin(r.theta); x1 = side * 2.1
            ax.plot([x0, x1, x1 + side * 0.1], [y0, yl, yl], color=T.INK["muted"], lw=0.35, zorder=3)
            ax.text(x1 + side * 0.14, yl, r["name"], ha="left" if side > 0 else "right", va="center", fontsize=T.PT_MIN, style="italic", color=T.INK["primary"])
    ang = np.radians(-28)   # percentile scale along one radius
    for p_, lab in ((0, "0"), (50, "50"), (100, "100th")):
        rr = rp(p_)
        ax.text(rr * np.cos(ang) + 0.04, rr * np.sin(ang) - 0.04, lab, fontsize=T.PT_MIN, color=T.INK["muted"], ha="left", va="top", zorder=7)
    fig.text(0.01, 0.975, "a  the sample on the taxonomy tree", fontsize=T.PT_TITLE, va="top")
    # legend
    lg = fig.add_axes([0.78, 0.2, 0.2, 0.3]); lg.axis("off")
    items = [(T.DEV["low_deep"], "o", "below the 1st percentile"), (T.DEV["low_light"], "o", "1st–2.5th percentile"), (T.DEV["neutral"], "o", "within the healthy range"),
             (T.DEV["high_light"], "o", "above the 97.5th (none here)")]
    for k, (c, m, t) in enumerate(items):
        lg.scatter([0.06], [0.9 - 0.13 * k], s=14, color=c, marker=m, transform=lg.transAxes, clip_on=False); lg.text(0.14, 0.9 - 0.13 * k, t, va="center", fontsize=T.PT_MIN, transform=lg.transAxes)
    lg.scatter([0.06], [0.9 - 0.13 * 4], s=13, transform=lg.transAxes, clip_on=False, **T.MISSING_RING); lg.text(0.14, 0.9 - 0.13 * 4, "expected but missing", va="center", fontsize=T.PT_MIN, transform=lg.transAxes)
    lg.text(0.0, 0.9 - 0.13 * 5.3, "radius = percentile among\nhealthy adults; grey ring =\nhealthy range (2.5th–97.5th);\noutermost ticks = out of\nrange in ≥ 25 % of the cases", va="top", fontsize=T.PT_MIN, color=T.INK["secondary"], transform=lg.transAxes)
    T.save(fig, "fig18_tree_report", out_dir); plt.close(fig)


if __name__ == "__main__":
    make(sys.argv[1] if len(sys.argv) > 1 else f"{P}/figures")
