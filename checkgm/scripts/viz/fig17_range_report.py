"""Figure 17 — the healthy range and where one sample lands, read like a laboratory report (double column).
Example: the C. difficile infection case MGYA00694698 (assembly pipeline), scored against the baseline rebuilt without its
study. Each row is one feature: the grey bar is the healthy range (2.5th–97.5th percentile of healthy adults; darker
inner bar 25th–75th; tick = median; thin line 1st–99th), on an axis of fold difference from the typical (median) healthy
adult. Families: difference in centred log-ratio (log2 units); gene families (KEGG orthologs): log2 ratio of copies per
genome, against healthy adults of the same assembly-size class. The dot is the sample (blue below, orange above, grey
within the range; open ring at the left edge = expected but absent). Right column: fold difference and percentile.
Rows: the families below the range in at least a quarter of the study's cases, Enterobacteriaceae (expected to rise in
this infection), two typical families within the range and one expected-but-missing family; KOs out of range here and in
at least a quarter of the cases."""
import os, sys, gzip
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import tokens as T, example_case as EC

P = T.P; LN2 = np.log(2); XL = (-7.5, 5.5)
FAM = ["Lachnospiraceae", "Oscillospiraceae", "Butyricicoccaceae", "Eubacteriaceae", "Erysipelotrichaceae", "Bifidobacteriaceae", "Bacteroidaceae", "Enterobacteriaceae"]
KOS = ["K20707", "K01879", "K00205", "K07084", "K01006", "K03088"]


def ko_names():
    nm = {}
    with gzip.open(f"{P}/staging/cmd3_humann/mapping/map_ko_name.txt.gz", "rt") as fh:
        for line in fh:
            k, v = line.rstrip("\n").split("\t", 1); nm[k] = v.split(" [EC")[0]
    return nm


def rows():
    sc = pd.read_parquet(f"{P}/work/s8/scores/{EC.STUDY}.parquet", filters=[("analysis_id", "==", EC.CASE)])
    band = sc.band.iloc[0]; ref_f = EC.reference(); ref_k = pd.read_parquet(f"{EC.BUNDLE}/features/ko_eggnog.parquet").set_index("feature_id")
    names = EC.names(); inv = {v: k for k, v in names.items()}; out = []
    fam = sc[sc.layer == "taxonomy_family"].set_index("feature_id")
    miss = fam[fam.call == "expected_but_missing"].index.tolist()
    prev = ref_f[[c for c in ref_f.columns if c.startswith("prevalence_band_")]].max(axis=1)
    within = [names.get(i, i) for i in prev.reindex(fam.index[fam.call == "within"]).sort_values(ascending=False).index[:2]]
    miss = sorted(miss, key=lambda i: -prev.get(i, 0))
    for n in FAM + within + [names.get(miss[0], miss[0])]:
        fid = inv.get(n, n); r = ref_f.loc[fid]; q = {k: r[f"p{k}"] for k in ["1", "2_5", "25", "50", "75", "97_5", "99"]}
        row = fam.loc[fid] if fid in fam.index else None
        v = None if row is None or pd.isna(row["value"]) else (row["value"] - q["50"]) / LN2
        out.append({"group": "Bacterial families", "label": n, "lo99": (q["1"] - q["50"]) / LN2, "lo": (q["2_5"] - q["50"]) / LN2, "q1": (q["25"] - q["50"]) / LN2, "q3": (q["75"] - q["50"]) / LN2,
                    "hi": (q["97_5"] - q["50"]) / LN2, "hi99": (q["99"] - q["50"]) / LN2, "x": v, "pct": None if row is None else row["percentile"], "call": "absent" if row is None or row["call"] == "expected_but_missing" else row["call"]})
    ko = sc[sc.layer == "ko_eggnog"].set_index("feature_id"); kn = ko_names()
    for k in KOS:
        r = ref_k.loc[k]; g = lambda p: r.get(f"p{p}_band_{band}", r[f"p{p}"]); med = g("50"); f = lambda p: np.log2(g(p) / med)
        row = ko.loc[k]; lab = kn.get(k, k).split(";")[-1].strip(); lab = {"RNA polymerase sigma-70 factor, ECF subfamily": "ECF sigma factor", "pyruvate, orthophosphate dikinase": "pyruvate phosphate dikinase", "putative amino acid transporter": "amino acid transporter", "glycyl-tRNA synthetase beta chain": "glycyl-tRNA synthetase β"}.get(lab, lab)
        out.append({"group": "Gene families (KEGG orthologs)", "label": f"{lab} ({k})", "lo99": f("1"), "lo": f("2_5"), "q1": f("25"), "q3": f("75"), "hi": f("97_5"), "hi99": f("99"),
                    "x": np.log2(row["value"] / med), "pct": row["percentile"], "call": row["call"]})
    return pd.DataFrame(out)


def colour(call, pct):
    if call == "low":
        return T.DEV["low_deep"] if pct < 1 else T.DEV["low_light"]
    if call == "high":
        return T.DEV["high_deep"] if pct > 99 else T.DEV["high_light"]
    return T.DEV["neutral"]


def fold_txt(x):
    f = 2 ** x
    if f >= 1:
        return f"{f:.1f}× higher" if f < 10 else f"{f:.0f}× higher"
    g = 1 / f
    return f"{g:.1f}× lower" if g < 10 else f"{g:,.0f}× lower"


def ordinal(p):
    if p < 1:
        return "below the 1st"
    if p > 99:
        return "above the 99th"
    n = int(round(p)); suf = "th" if 10 <= n % 100 <= 20 else {1: "st", 2: "nd", 3: "rd"}.get(n % 10, "th")
    return f"{n}{suf}"


def tick_lab(t):
    if t == 0:
        return "typical"
    f = 2 ** t
    return f"{f:g}×" if f >= 1 else f"1/{1 / f:g}"


def panel(ax, d, xl, ticks, title):
    y = 0; ys = []
    for _, r in d.iterrows():
        lo99, hi99 = max(r.lo99, xl[0]), min(r.hi99, xl[1])
        ax.plot([lo99, hi99], [y, y], color=T.INK["axis"], lw=0.6, zorder=1, solid_capstyle="butt")
        ax.add_patch(plt.Rectangle((max(r.lo, xl[0]), y - 0.3), min(r.hi, xl[1]) - max(r.lo, xl[0]), 0.6, facecolor=T.INK["grid"], edgecolor="none", zorder=2))
        ax.add_patch(plt.Rectangle((r.q1, y - 0.3), r.q3 - r.q1, 0.6, facecolor=T.INK["axis"], edgecolor="none", zorder=3))
        ax.plot([0, 0], [y - 0.34, y + 0.34], color=T.INK["secondary"], lw=0.8, zorder=4)
        if r.call == "absent":
            ax.scatter([xl[0] + 0.3], [y], s=22, zorder=5, **T.MISSING_RING); txt = "absent; carried by ≥ 70 % of healthy adults"
        else:
            out_lo, out_hi = r.x < xl[0] + 0.3, r.x > xl[1] - 0.3; xs = float(np.clip(r.x, xl[0] + 0.3, xl[1] - 0.3))
            ax.scatter([xs], [y], s=30, marker="<" if out_lo else (">" if out_hi else "o"), color=colour(r.call, r.pct), edgecolor=T.INK["surface"], linewidth=0.5, zorder=6)
            txt = f"{fold_txt(r.x)} · {ordinal(r.pct)} percentile"
        ax.text(1.02, y, txt, va="center", ha="left", fontsize=T.PT_MIN, color=T.INK["secondary"], transform=ax.get_yaxis_transform(), clip_on=False)
        ys.append(y); y += 1
    ax.set_yticks(ys); ax.set_yticklabels(d.label, color=T.INK["primary"]); ax.tick_params(axis="y", length=0, labelsize=T.PT_MIN)
    ax.set_ylim(y - 0.4, -0.7); ax.set_xlim(*xl); ax.set_xticks(ticks); ax.set_xticklabels([tick_lab(t) for t in ticks])
    ax.grid(axis="x", color=T.INK["grid"], lw=0.4, zorder=0); ax.set_axisbelow(True)
    for sp in ["left", "top", "right"]:
        ax.spines[sp].set_visible(False)
    ax.set_title(title, loc="left", fontsize=T.PT_TITLE, pad=4)


def make(out_dir):
    T.print_profile(); d = rows()
    fam = d[d.group.str.startswith("Bacterial")]; ko = d[~d.group.str.startswith("Bacterial")]
    fig = plt.figure(figsize=(T.DOUBLE_IN, 4.3))
    gs = fig.add_gridspec(2, 1, height_ratios=[len(fam) + 1.5, len(ko) + 1.5], hspace=0.45, left=0.25, right=0.71, top=0.94, bottom=0.09)
    a = fig.add_subplot(gs[0]); panel(a, fam, (-12, 6), [-10, -8, -6, -4, -2, 0, 2, 4, 6], "b  bacterial families: healthy range and this sample")
    for lab in a.get_yticklabels():
        lab.set_style("italic")
    b = fig.add_subplot(gs[1]); panel(b, ko, (-5, 4), [-4, -2, 0, 2, 4], "c  gene families (KEGG orthologs): healthy range and this sample")
    b.set_xlabel("difference from the typical healthy adult (median); grey bar = healthy range (2.5th–97.5th percentile), darker = middle half", fontsize=T.PT_MIN)
    T.save(fig, "fig17_range_report", out_dir); plt.close(fig)


if __name__ == "__main__":
    make(sys.argv[1] if len(sys.argv) > 1 else f"{P}/figures")
