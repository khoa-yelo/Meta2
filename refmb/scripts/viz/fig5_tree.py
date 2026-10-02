"""Figure 5 — radial family report of one sample (double column). Leaves = the 232 gut families carried by at least half
of healthy samples (reference prevalence >= 50 %), ordered phylum > class > order so relatives are adjacent. The percentile
track runs from the inner edge (0) to the outer edge (100) with the healthy range (2.5–97.5) shaded; each leaf's mark is
this case's percentile against a baseline built without its study. Families outside the range that are also outside in
>= 25 % of the study's cases are labelled; expected-but-missing families are hollow rings. (a) wheel; (b) inset: one
labelled low family's full baseline distribution with this case marked, joined to its leaf by a leader line.
Every count in the title and centre is computed from the scores table (example_case.py), never typed."""
import os, sys
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge, ConnectionPatch
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import tokens as T, example_case as EC

P = T.P
STUDY, CASE = EC.STUDY, EC.CASE
INSET_FAMILY = "Lachnospiraceae"   # a labelled low family, so the inset and the wheel tell one story


def ordinal(p):
    s = f"{p:.1f}" if p < 10 else f"{p:.0f}"
    n = int(round(p)); suf = "th" if 10 <= n % 100 <= 20 else {1: "st", 2: "nd", 3: "rd"}.get(n % 10, "th")
    return s + ("th" if "." in s else suf)


def load():
    inv = pd.read_parquet(f"{P}/work/s0/inventory.parquet", columns=["analysis_id", "study_bioproject", "cur_Health_status_group", "is_primary_analysis"]).set_index("analysis_id")
    sc = pd.read_parquet(f"{P}/work/s8/scores/{STUDY}.parquet"); sc = sc[sc["layer"] == "taxonomy_family"]
    cases = [a for a in sc["analysis_id"].unique() if inv.loc[a, "study_bioproject"] == STUDY and inv.loc[a, "cur_Health_status_group"] == "Diseased" and inv.loc[a, "is_primary_analysis"]]
    sc = sc[sc["analysis_id"].isin(cases)]
    ref = EC.reference(); leaves = EC.leaves(ref)
    tax = pd.read_parquet(f"{P}/resources/backbone/ncbi_taxonomy.parquet", columns=["taxid", "name", "phylum_taxid", "class_taxid", "order_taxid"])
    tax["taxid"] = tax["taxid"].astype(str); tx = tax.set_index("taxid"); names = dict(zip(tax["taxid"], tax["name"]))
    # representative case: the sample whose out-of-range families agree most with the cohort-frequent deviations
    freq0 = sc[sc["feature_id"].isin(leaves)].assign(out=lambda t: t["call"].isin(["low", "high"])).groupby("feature_id")["out"].mean()
    frequent = set(freq0.index[freq0 >= 0.25])
    agree = sc[sc["feature_id"].isin(frequent) & sc["call"].isin(["low", "high"])].groupby("analysis_id").size()
    rep = agree.idxmax() if len(agree) else sc["analysis_id"].iloc[0]
    if rep != CASE:
        print(f"fig5: representative case by rule is {rep}, but the shared example (fig1, summary) is {CASE}; drawing {CASE}", file=sys.stderr)
    one = EC.case_families(CASE)
    freq_out = sc[sc["feature_id"].isin(leaves)].assign(hi=lambda t: t["call"] == "high", lo=lambda t: t["call"] == "low").groupby("feature_id")[["hi", "lo"]].mean()
    d = pd.DataFrame({"pct": one["percentile"], "call": one["call"], "name": one["name"],
                      "cohort_high": freq_out["hi"].reindex(leaves).fillna(0), "cohort_low": freq_out["lo"].reindex(leaves).fillna(0)})
    for r in ["phylum_taxid", "class_taxid", "order_taxid"]:
        d[r] = [names.get(str(int(tx.loc[i, r])), "?") if i in tx.index and tx.loc[i, r] != -1 else "?" for i in d.index]
    d = d.sort_values(["phylum_taxid", "class_taxid", "order_taxid", "name"])
    return d, ref, len(cases), EC.counts(one)


def make(out_dir):
    T.print_profile()
    d, ref, ncases, c = load()
    n = len(d); theta = np.linspace(0, 2 * np.pi, n, endpoint=False) + np.pi / 2
    r_in, r_out = 1.0, 2.0                      # percentile 0 -> r_in, 100 -> r_out
    rp = lambda p: r_in + (r_out - r_in) * np.clip(p, 0, 100) / 100
    fig = plt.figure(figsize=(T.DOUBLE_IN, 4.6))
    ax = fig.add_axes([0.0, 0.0, 0.66, 0.93]); ax.set_aspect("equal"); ax.axis("off")
    ax.set_xlim(-4.4, 3.4); ax.set_ylim(-3.0, 3.0)   # left-heavy so the left column of family labels is not clipped
    ax.add_patch(Wedge((0, 0), rp(97.5), 0, 360, width=rp(97.5) - rp(2.5), facecolor=T.INK["grid"], edgecolor="none", zorder=0))
    for p in (2.5, 50, 97.5):
        ax.add_patch(plt.Circle((0, 0), rp(p), fill=False, edgecolor=T.INK["axis"], lw=0.4, zorder=1))
    ax.add_patch(plt.Circle((0, 0), r_in, fill=False, edgecolor=T.INK["axis"], lw=0.4)); ax.add_patch(plt.Circle((0, 0), r_out, fill=False, edgecolor=T.INK["axis"], lw=0.4))
    # phylum arcs (outer ring); labels anchored on their outward side so long names extend away from the ring
    for ph, g in d.groupby("phylum_taxid", sort=False):
        idx = np.where(d["phylum_taxid"].to_numpy() == ph)[0]; a0, a1 = np.degrees(theta[idx.min()]) - 180 / n, np.degrees(theta[idx.max()]) + 180 / n
        ax.add_patch(Wedge((0, 0), r_out + 0.2, a0, a1, width=0.03, facecolor=T.INK["axis"], edgecolor="none"))
        am = np.radians((a0 + a1) / 2); rr = r_out + 0.3; ca, sa = np.cos(am), np.sin(am)
        if len(idx) >= 8:
            ax.text(rr * ca, rr * sa, ph, ha="left" if ca > 0.3 else "right" if ca < -0.3 else "center", va="bottom" if sa > 0.3 else "top" if sa < -0.3 else "center", fontsize=T.PT_MIN, color=T.INK["secondary"])
    labels = []
    # outer ring: cohort deviation frequency (share of cases high = orange, low = blue), light >= 25 %, deep >= 50 %
    for i, (fid, row) in enumerate(d.iterrows()):
        a0, a1 = np.degrees(theta[i]) - 180 / n * 0.9, np.degrees(theta[i]) + 180 / n * 0.9
        fh, fl = row["cohort_high"], row["cohort_low"]
        if max(fh, fl) >= 0.25:
            col = (T.DEV["high_deep"] if fh >= 0.5 else T.DEV["high_light"]) if fh >= fl else (T.DEV["low_deep"] if fl >= 0.5 else T.DEV["low_light"])
            ax.add_patch(Wedge((0, 0), r_out + 0.12, a0, a1, width=0.09, facecolor=col, edgecolor="none", zorder=2))
    # spokes + per-sample marks
    leaf_xy = {}
    for i, (fid, row) in enumerate(d.iterrows()):
        cs, sn = np.cos(theta[i]), np.sin(theta[i])
        ax.plot([r_in * cs, r_out * cs], [r_in * sn, r_out * sn], color=T.INK["grid"], lw=0.3, zorder=1)
        p = row["pct"]
        if row["call"] == "expected_but_missing":
            ax.scatter([rp(0) * cs], [rp(0) * sn], s=14, zorder=4, **T.MISSING_RING); continue
        if np.isnan(p):
            continue
        if p < 2.5: col, hatch = (T.DEV["low_deep"] if p < 1 else T.DEV["low_light"]), "////"
        elif p > 97.5: col, hatch = (T.DEV["high_deep"] if p > 99 else T.DEV["high_light"]), "\\\\"
        else: col, hatch = T.DEV["neutral"], None
        out = p < 2.5 or p > 97.5
        ax.scatter([rp(p) * cs], [rp(p) * sn], s=16 if out else 6, color=col, edgecolor=T.INK["surface"], linewidth=0.3, zorder=5 if out else 3, hatch=hatch)
        leaf_xy[row["name"]] = (rp(p) * cs, rp(p) * sn, p)
        if out and max(row["cohort_high"], row["cohort_low"]) >= 0.25:
            labels.append((theta[i], rp(p), row["name"]))
    # leader labels: horizontal text in a left and a right column, spaced >= 0.17 apart vertically
    for side in (-1, 1):
        col = sorted([(t, r, nm) for t, r, nm in labels if (np.cos(t) >= 0) == (side > 0)], key=lambda x: -np.sin(x[0]))
        ys = [2.55 * np.sin(t) for t, _, _ in col]
        for k in range(1, len(ys)):
            ys[k] = min(ys[k], ys[k - 1] - 0.17)
        for (t, r, nm), yl in zip(col, ys):
            x0, y0 = r * np.cos(t), r * np.sin(t); x1 = side * 2.45; x2 = side * 2.62
            ax.plot([x0, x1, x2], [y0, yl, yl], color=T.INK["muted"], lw=0.35, zorder=4)
            ax.text(x2 + 0.04 * side, yl, nm, ha="left" if side > 0 else "right", va="center", fontsize=T.PT_MIN, style="italic", color=T.INK["primary"])
    ax.text(0, 0, f"{c['n_families']} families\ncase {CASE}\n{c['n_outside']}/{c['n_assessed']} assessed\nfamilies outside\n({c['frac_outside']:.0%})", ha="center", va="center", fontsize=T.PT_MIN, color=T.INK["secondary"], linespacing=1.25)
    ang = np.radians(-28)   # radius labels along one direction clear of the expected-but-missing rings at the bottom
    for p_, lab, dx in [(0, "0", -0.05), (2.5, "2.5", 0.09), (50, "50", 0), (97.5, "97.5", -0.1), (100, "100", 0.16)]:
        rr = rp(p_) + (0.09 if p_ == 100 else -0.05 if p_ == 97.5 else -0.11 if p_ == 0 else 0.07)
        ax.text(rr * np.cos(ang) + dx, rr * np.sin(ang), lab, fontsize=T.PT_MIN, color=T.INK["muted"], ha="center", va="center")
    fig.text(0.02, 0.925, "a", fontsize=T.PT_TITLE, va="top", ha="left")
    fig.text(0.02, 0.975, f"One C. difficile infection case: {c['n_low']} of {c['n_assessed']} assessed families below the healthy range, {EC.none_or(c['n_high'])} above", fontsize=T.PT_TITLE, va="top")
    # legend (right column)
    lg = fig.add_axes([0.66, 0.5, 0.33, 0.44]); lg.axis("off")
    none_hi = "  (none in this case)" if c["n_high"] == 0 else ""
    items = [(T.DEV["high_deep"], "high: > 99th percentile" + none_hi), (T.DEV["high_light"], "high: 97.5th–99th" + none_hi), (T.DEV["neutral"], "within the healthy range (2.5th–97.5th)"),
             (T.DEV["low_light"], "low: 1st–2.5th"), (T.DEV["low_deep"], "low: < 1st percentile")]
    for k, (col, lab) in enumerate(items):
        lg.scatter([0.05], [0.95 - 0.11 * k], s=16, color=col, edgecolor=T.INK["surface"], linewidth=0.3, transform=lg.transAxes); lg.text(0.12, 0.95 - 0.11 * k, lab, va="center", fontsize=T.PT_MIN, transform=lg.transAxes)
    lg.scatter([0.05], [0.95 - 0.11 * 5], s=14, transform=lg.transAxes, **T.MISSING_RING); lg.text(0.12, 0.95 - 0.11 * 5, "expected but missing\n(carried by > 50 % of healthy adults)", va="center", fontsize=T.PT_MIN, transform=lg.transAxes, linespacing=1.2)
    lg.add_patch(plt.Rectangle((0.03, 0.95 - 0.125 * 6 - 0.012), 0.04, 0.024, facecolor=T.DEV["low_deep"], edgecolor="none", transform=lg.transAxes))
    lg.text(0.12, 0.95 - 0.125 * 6, f"outer ring: share of the {ncases} cases outside\n(light ≥ 25 %, deep ≥ 50 %)", va="center", fontsize=T.PT_MIN, transform=lg.transAxes, linespacing=1.2)
    # inset (b): one labelled low family's baseline distribution, joined to its leaf
    fid = d.index[d["name"] == INSET_FAMILY][0] if (d["name"] == INSET_FAMILY).any() else d["pct"].idxmin()
    ins = fig.add_axes([0.73, 0.1, 0.25, 0.3])
    grid = [1, 2.5, 5, 10, 25, 50, 75, 90, 95, 97.5, 99]; vals = [ref.loc[fid, f"p{str(g).replace('.', '_')}"] for g in grid]
    ins.step(vals, grid, where="post", color=T.INK["primary"], lw=1.0); ins.axhspan(2.5, 97.5, color=T.INK["grid"], lw=0, zorder=0)
    pm = d.loc[fid, "pct"]; vm = np.interp(pm, grid, vals)
    mcol = (T.DEV["low_deep"] if pm < 1 else T.DEV["low_light"]) if pm < 2.5 else (T.DEV["high_deep"] if pm > 99 else T.DEV["high_light"]) if pm > 97.5 else T.DEV["neutral"]
    ins.scatter([vm], [pm], s=20, color=mcol, edgecolor=T.INK["surface"], linewidth=0.3, zorder=3)
    ins.annotate(f"this case: {ordinal(pm)} percentile", xy=(vm, pm), xytext=(0.06, 0.78), textcoords="axes fraction", ha="left", va="bottom", fontsize=T.PT_MIN, color=T.INK["secondary"],
                 arrowprops=dict(arrowstyle="-", color=T.INK["muted"], lw=0.4, shrinkA=0, shrinkB=3))
    ins.set_yticks([2.5, 50, 97.5]); ins.set_yticklabels(["2.5", "50", "97.5"]); ins.set_ylim(0, 100)
    ins.set_xlabel(f"CLR abundance of {d.loc[fid, 'name']}", fontsize=T.PT_MIN, style="italic"); ins.set_ylabel("percentile in healthy adults", fontsize=T.PT_MIN); ins.tick_params(labelsize=T.PT_MIN, length=2, width=0.4)
    ins.set_title("b  how a percentile is read", loc="left", fontsize=T.PT_BODY, pad=3)
    nm = d.loc[fid, "name"]
    if nm in leaf_xy:
        lx, ly, _ = leaf_xy[nm]
        fig.add_artist(ConnectionPatch(xyA=(vm, pm), coordsA=ins.transData, xyB=(lx, ly), coordsB=ax.transData, color=T.INK["muted"], lw=0.4, ls=(0, (2, 2)), zorder=1))
    T.save(fig, "fig5_tree", out_dir); plt.close(fig)


if __name__ == "__main__":
    make(sys.argv[1] if len(sys.argv) > 1 else f"{P}/figures")
