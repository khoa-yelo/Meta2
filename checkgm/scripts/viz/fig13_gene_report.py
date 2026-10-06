"""Figure 13 — what a user gets on the gene layers (double column, two radial reports side by side, about 100 mm tall): the
radial deviation map of Fig. 5 applied to (a) KEGG orthologs (eggNOG annotation) for the same C. difficile case (assembly
pipeline) and to (b) KEGG orthologs (HUMAnN 3) for one colorectal cancer case (read pipeline, FengQ_2015). The sparse
layers of the same two samples, KEGG module completeness (assembly) and MetaCyc pathways (read), are not drawn: their rings
carried almost no information at print size, so each is summarised in one text line under the ring of its sample
("1 of 188 modules outside ...", "5 of 320 pathways outside ..."). Leaves are the features with reference prevalence >= 50 % in
the sample's class, ordered by functional class (KEGG BRITE ko00001 level 2). Mark radius = the sample's percentile; grey
ring = 2.5–97.5; outer ring = histogram per 5° sector of features called high (orange) or low (blue) in >= 25 % of the
cohort's cases; numbered = out of range in this sample and in >= 25 % of the cohort, the MAX_LAB most frequent listed in full.
Per-layer counts of both samples on all four layers are compared with results/s14/example_reports.tsv (the data columns:
features, assessed, low, high, expected_but_missing, cohort_cases); a differing table is written beside the figure for review.
The TSV's 'labelled' column is a presentation count from an earlier version and is not compared."""
import json, os, re, sys
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import tokens as T
sys.path.insert(0, T.P); from checkgm.paths import bundle_A as _bundle_A

P = T.P
A_STUDY, A_CASE = "PRJEB26165", "MGYA00694698"          # the case shown in Fig. 5
B_STUDY = "FengQ_2015"
MIN_PREV = 0.5; FREQ = 0.25; MAX_LAB = 8   # the eight most frequent shared deviations are numbered and listed in full
DATA_COLS = ["features", "assessed", "low", "high", "expected_but_missing", "cohort_cases"]   # compared with results/s14/example_reports.tsv
SHORT = {"Carbohydrate metabolism": "Carbohydrate", "Energy metabolism": "Energy", "Lipid metabolism": "Lipid", "Nucleotide metabolism": "Nucleotide",
         "Amino acid metabolism": "Amino acid", "Metabolism of other amino acids": "Other amino acids", "Glycan biosynthesis and metabolism": "Glycan",
         "Metabolism of cofactors and vitamins": "Cofactors, vitamins", "Metabolism of terpenoids and polyketides": "Terpenoids", "Biosynthesis of other secondary metabolites": "Secondary metab.",
         "Xenobiotics biodegradation and metabolism": "Xenobiotics", "Genetic Information Processing": "Genetic information", "Environmental Information Processing": "Environmental information",
         "Cellular Processes": "Cellular processes", "Organismal Systems": "Organismal", "Human Diseases": "Disease-associated", "Not Included in Pathway or Brite": "Unclassified",
         "Brite Hierarchies": "Unclassified", "Biosynthesis of terpenoids and polyketides": "Terpenoids", "Glycan metabolism": "Glycan", "Xenobiotics biodegradation": "Xenobiotics"}


def brite():
    """KO -> (level-2 class, short name) from KEGG BRITE ko00001; metabolism classes take precedence, then the first other class."""
    j = json.load(open(f"{P}/staging/kegg/ko00001.json")); cls = {}; names = {}
    for l1 in j["children"]:
        top = re.sub(r"^\d+\s+", "", l1["name"])
        for l2 in l1.get("children", []):
            c2 = re.sub(r"^\d+\s+", "", l2["name"]); c = SHORT.get(c2 if top == "Metabolism" else top, c2 if top == "Metabolism" else top)
            for l3 in l2.get("children", []):
                for l4 in l3.get("children", []):
                    m = re.match(r"(K\d{5})\s+(.*)", l4["name"])
                    if not m:
                        continue
                    k, nm = m.groups(); names.setdefault(k, nm)
                    if k not in cls or (top == "Metabolism" and cls[k] not in METAB):
                        cls[k] = c
    return cls, names


METAB = {"Carbohydrate", "Energy", "Lipid", "Nucleotide", "Amino acid", "Other amino acids", "Glycan", "Cofactors, vitamins", "Terpenoids", "Secondary metab.", "Xenobiotics"}
ORDER = ["Carbohydrate", "Energy", "Lipid", "Nucleotide", "Amino acid", "Other amino acids", "Glycan", "Cofactors, vitamins", "Terpenoids", "Secondary metab.", "Xenobiotics",
         "Genetic information", "Environmental information", "Cellular processes", "Organismal", "Disease-associated", "Unclassified"]


def ko_label(fid, nm):
    """KO id plus the full KEGG short name (EC numbers dropped; of alternative names separated by ' / ' the first is kept)."""
    sym, _, desc = nm.partition(";"); desc = re.sub(r"\s*\[EC:[^\]]*\]", "", desc).strip().split(" / ")[0].strip(); sym = sym.split(",")[0].strip()
    return f"{fid}  {desc or sym}"


def table(scores, ref, cases, sample, names, classes, name_short):
    pc = [c for c in ref.columns if c.startswith("prevalence_band_")]; prev = ref[pc].max(axis=1)
    leaves = ref.index[(prev >= MIN_PREV) & ref["percentiles_available"]]
    sc = scores[scores["feature_id"].isin(leaves)]; one = sc[sc["analysis_id"] == sample].set_index("feature_id")
    fr = sc.assign(hi=sc["call"] == "high", lo=sc["call"] == "low").groupby("feature_id")[["hi", "lo"]].mean()
    d = pd.DataFrame({"pct": one["percentile"].reindex(leaves), "call": one["call"].reindex(leaves).fillna("undetected"), "cohort_high": fr["hi"].reindex(leaves).fillna(0), "cohort_low": fr["lo"].reindex(leaves).fillna(0)})
    d["name"] = [name_short(i, names.get(i, i)) for i in d.index]; d["cls"] = [classes.get(i, "Unclassified") for i in d.index]
    d["ord"] = d["cls"].map({c: k for k, c in enumerate(ORDER)}).fillna(len(ORDER)); d = d.sort_values(["ord", "cls", "name"]); return d


def summary(d, ncases):
    """Counts of one layer for one sample: features (leaves), assessed, low, high, expected-but-missing, cohort cases and the
    number of features meeting the labelling criterion (outside in this sample and in >= FREQ of the cohort)."""
    p = d["pct"].to_numpy(); call = d["call"].to_numpy()
    ok = ~np.isnan(p) & np.isin(call, ["within", "low", "high"]); lo_ = ok & (call == "low"); hi_ = ok & (call == "high")
    freq = np.maximum(d["cohort_high"].to_numpy(), d["cohort_low"].to_numpy())
    return {"features": int(len(d)), "assessed": int(ok.sum()), "low": int(lo_.sum()), "high": int(hi_.sum()), "expected_but_missing": int((call == "expected_but_missing").sum()),
            "cohort_cases": int(ncases), "labelled": int(((lo_ | hi_) & (freq >= FREQ)).sum())}


def radial(ax, d, title, ncases, cls_min=12):
    n = len(d); theta = np.linspace(0, 2 * np.pi, n, endpoint=False) + np.pi / 2; r_in, r_out = 1.25, 2.0
    rp = lambda p: r_in + (r_out - r_in) * np.clip(p, 0, 100) / 100
    ax.set_aspect("equal"); ax.axis("off"); ax.set_xlim(-4.8, 4.8); ax.set_ylim(-3.72, 3.72); k_cls = 0   # 9.6 x 7.44 data units; room for the class labels
    ax.add_patch(Wedge((0, 0), rp(97.5), 0, 360, width=rp(97.5) - rp(2.5), facecolor=T.INK["grid"], edgecolor="none", zorder=0))
    for p in (2.5, 50, 97.5):
        ax.add_patch(plt.Circle((0, 0), rp(p), fill=False, edgecolor=T.INK["axis"], lw=0.4, zorder=1))
    for r in (r_in, r_out):
        ax.add_patch(plt.Circle((0, 0), r, fill=False, edgecolor=T.INK["axis"], lw=0.4))
    for c in d["cls"].unique():
        idx = np.where(d["cls"].to_numpy() == c)[0]; a0, a1 = np.degrees(theta[idx.min()]) - 180 / n, np.degrees(theta[idx.max()]) + 180 / n
        ring = 0.2 if n <= 400 else 0.5
        ax.add_patch(Wedge((0, 0), r_out + ring, a0, a1, width=0.03, facecolor=T.INK["axis"], edgecolor="none"))
        am = np.radians((a0 + a1) / 2)
        if len(idx) >= max(cls_min, 0.07 * n):
            rr = r_out + ring + ((0.3, 0.58) if n <= 400 else (0.55, 0.85))[k_cls % 2]; k_cls += 1; c = c.replace(" information", "\ninformation").replace("Cofactors, vitamins", "Cofactors,\nvitamins").replace("Cellular processes", "Cellular\nprocesses")
            ax.text(rr * np.cos(am), rr * np.sin(am), c, ha="center" if abs(np.cos(am)) < 0.35 else ("left" if np.cos(am) > 0 else "right"), va="center", fontsize=T.PT_MIN, color=T.INK["secondary"])
    fh_all, fl_all = d["cohort_high"].to_numpy(), d["cohort_low"].to_numpy()
    if n <= 400:   # sparse layer: one wedge per feature, as in Fig. 5
        for i in range(n):
            fh, fl = fh_all[i], fl_all[i]
            if max(fh, fl) >= FREQ:
                a0, a1 = np.degrees(theta[i]) - 180 / n * 0.9, np.degrees(theta[i]) + 180 / n * 0.9
                col = (T.DEV["high_deep"] if fh >= 0.5 else T.DEV["high_light"]) if fh >= fl else (T.DEV["low_deep"] if fl >= 0.5 else T.DEV["low_light"])
                ax.add_patch(Wedge((0, 0), r_out + 0.12, a0, a1, width=0.09, facecolor=col, edgecolor="none", zorder=2))
    else:          # dense layer: per-feature wedges merge at print size, so bin the features into NBIN sectors and stack counts (low inside, high outside)
        NBIN = 72; sector = (np.arange(n) * NBIN) // n
        hi_c = np.bincount(sector[(fh_all >= FREQ) & (fh_all >= fl_all)], minlength=NBIN); lo_c = np.bincount(sector[(fl_all >= FREQ) & (fl_all > fh_all)], minlength=NBIN)
        scale = 0.42 / max(1, (hi_c + lo_c).max()); deg = 360 / NBIN
        for b in range(NBIN):
            a0 = 90 + b * deg * 1.0 + deg * 0.08; a1 = a0 + deg * 0.84   # sectors follow the leaf order (theta starts at 90°)
            r0 = r_out + 0.04
            if lo_c[b]:
                ax.add_patch(Wedge((0, 0), r0 + lo_c[b] * scale, a0, a1, width=lo_c[b] * scale, facecolor=T.DEV["low_deep"], edgecolor="none", zorder=2)); r0 += lo_c[b] * scale
            if hi_c[b]:
                ax.add_patch(Wedge((0, 0), r0 + hi_c[b] * scale, a0, a1, width=hi_c[b] * scale, facecolor=T.DEV["high_deep"], edgecolor="none", zorder=2))
        ax.add_patch(plt.Circle((0, 0), r_out + 0.04, fill=False, edgecolor=T.INK["axis"], lw=0.3))
    c_, s_ = np.cos(theta), np.sin(theta); p = d["pct"].to_numpy(); call = d["call"].to_numpy()
    miss = call == "expected_but_missing"; ax.scatter(rp(0) * c_[miss], rp(0) * s_[miss], s=9, zorder=4, **T.MISSING_RING)
    ok = ~np.isnan(p) & np.isin(call, ["within", "low", "high"]); lo_ = ok & (call == "low"); hi_ = ok & (call == "high"); within = ok & ~lo_ & ~hi_
    ax.scatter(rp(p[within]) * c_[within], rp(p[within]) * s_[within], s=1.4 if n > 400 else 6, color=T.DEV["neutral"], linewidth=0, zorder=3)
    for m, col in [(lo_ & (p < 1), T.DEV["low_deep"]), (lo_ & (p >= 1), T.DEV["low_light"]), (hi_ & (p <= 99), T.DEV["high_light"]), (hi_ & (p > 99), T.DEV["high_deep"])]:
        ax.scatter(rp(p[m]) * c_[m], rp(p[m]) * s_[m], s=11, color=col, edgecolor=T.INK["surface"], linewidth=0.3, zorder=5)
    out = lo_ | hi_; freq = np.maximum(d["cohort_high"].to_numpy(), d["cohort_low"].to_numpy())
    cand = sorted([(theta[i], rp(p[i]), d["name"].iloc[i], freq[i], d["cohort_high"].iloc[i] >= d["cohort_low"].iloc[i]) for i in np.where(out & (freq >= FREQ))[0]], key=lambda x: -x[3])[:MAX_LAB]
    cand = sorted(cand, key=lambda x: -x[0]); listed = []
    r_lab = r_out + (0.3 if n <= 400 else 0.7); min_gap = np.radians(7)
    ang = np.array([t for t, *_ in cand], float)
    for _ in range(60):   # relax: push neighbours apart along the circle until every gap >= min_gap
        moved = False
        for j in range(len(ang) - 1):
            gap = ang[j] - ang[j + 1]
            if gap < min_gap:
                push = (min_gap - gap) / 2; ang[j] += push; ang[j + 1] -= push; moved = True
        if not moved:
            break
    for j, ((t, r, nm, fq, hi), ta) in enumerate(zip(cand, ang)):
        ax.plot([r * np.cos(t), (r_lab - 0.12) * np.cos(ta)], [r * np.sin(t), (r_lab - 0.12) * np.sin(ta)], color=T.INK["muted"], lw=0.35, zorder=5)
        ax.text(r_lab * np.cos(ta), r_lab * np.sin(ta), str(j + 1), ha="center", va="center", fontsize=T.PT_MIN, color=T.INK["primary"], zorder=6)
        listed.append(f"{j + 1}  {nm} · {fq:.0%} {'high' if hi else 'low'}")
    n_out = int(out.sum()); n_ass = int(ok.sum())
    ax.text(0, 0, f"{n:,} features\n{n_out} of {n_ass:,}\noutside ({n_out / max(1, n_ass):.0%})\n{int(miss.sum())} expected,\nmissing", ha="center", va="center", fontsize=T.PT_MIN, color=T.INK["secondary"])
    ax.set_title(title, loc="left", fontsize=T.PT_TITLE, pad=2)
    return listed, summary(d, ncases)


def make(out_dir):
    T.print_profile(); kcls, knames = brite(); rows = []
    # pipeline A: same case as Fig. 5
    inv = pd.read_parquet(f"{P}/work/s0/inventory.parquet", columns=["analysis_id", "study_bioproject", "cur_Health_status_group", "is_primary_analysis"]).set_index("analysis_id")
    scA = pd.read_parquet(f"{P}/work/s8/scores/{A_STUDY}.parquet"); cases = [a for a in scA["analysis_id"].unique() if inv.loc[a, "study_bioproject"] == A_STUDY and inv.loc[a, "cur_Health_status_group"] == "Diseased" and inv.loc[a, "is_primary_analysis"]]
    scA = scA[scA["analysis_id"].isin(cases)]; bA = _bundle_A("lenient", loso=A_STUDY, root=f"{P}/work/s8/refs")
    mf = pd.read_csv(f"{P}/resources/module_features.tsv", sep="\t").set_index("feature_id"); mcls = {i: SHORT.get(c.split(";")[1].strip(), c.split(";")[1].strip()) if c.startswith("Pathway") else "Signature" for i, c in mf["class"].items()}
    mnames = mf["name"].to_dict(); mshort = lambda fid, s: f"{fid}  {s}"
    # pipeline B: representative CRC case of B_STUDY: fraction of KOs outside near the cohort median, then most agreement with cohort-frequent deviations
    meta = pd.read_parquet(f"{P}/work/s12/B/meta.parquet").set_index("sample"); casesB = list(meta[(meta.study == B_STUDY) & (~meta.is_healthy)].index)
    X = pd.read_parquet(f"{P}/work/s12/B/oos_scores.parquet"); scB = X[X["sample"].isin(casesB)].rename(columns={"sample": "analysis_id"})
    k = scB[scB.layer == "ko_humann"]; fr = k.assign(o=k.call.isin(["low", "high"])).groupby("feature_id").o.mean(); frequent = set(fr.index[fr >= FREQ])
    outc = k.assign(o=k.call.isin(["low", "high"])).groupby("analysis_id").o.mean(); mid = outc[(outc >= outc.quantile(0.4)) & (outc <= outc.quantile(0.6))].index
    agree = k[k.analysis_id.isin(mid) & k.feature_id.isin(frequent) & k.call.isin(["low", "high"])].groupby("analysis_id").size(); repB = agree.idxmax() if len(agree) else outc.sub(outc.median()).abs().idxmin()
    bB = f"{P}/refs/gut-reads-adult-global-v0.2-lenient"; pnames = {}
    for f in sorted(os.listdir(f"{P}/staging/cmd3_humann/pathway_abundance"))[:8]:
        for x in pd.read_parquet(f"{P}/staging/cmd3_humann/pathway_abundance/{f}", columns=["feature_id"]).feature_id.unique():
            if ":" in x:
                pnames[x.split(":")[0].strip()] = x.split(":", 1)[1].strip().replace("&gamma;", "γ").replace("&beta;", "β").replace("&alpha;", "α")
    KW = [("biosynthesis", "Biosynthesis"), ("degradation", "Degradation"), ("fermentation", "Fermentation"), ("salvage", "Salvage"), ("glycolysis", "Energy"), ("tca", "Energy"), ("oxidat", "Energy"), ("respiration", "Energy"), ("photosynth", "Energy")]
    pcls = {i: next((c for kw, c in KW if kw in n.lower()), "Other") for i, n in pnames.items()}
    pshort = lambda fid, s: f"{fid}  {s}"
    refA = lambda layer: pd.read_parquet(f"{bA}/features/{layer}.parquet").set_index("feature_id")
    refB = lambda layer: pd.read_parquet(f"{bB}/features/{layer}.parquet").set_index("feature_id")
    # the two radial panels (dense KO layers) and, under each, the one-line summary of the sample's sparse layer
    panels = [("a  assembly pipeline, KEGG orthologs: C. difficile infection case", scA[scA.layer == "ko_eggnog"], refA("ko_eggnog"), cases, A_CASE, knames, kcls, ko_label, "ko_eggnog", "A",
               ("KEGG modules, same case", scA[scA.layer == "module"], refA("module"), mnames, mcls, mshort, "module")),
              ("b  read pipeline, KEGG orthologs: colorectal cancer case", scB[scB.layer == "ko_humann"], refB("ko_humann"), casesB, repB, knames, kcls, ko_label, "ko_humann", "B",
               ("MetaCyc pathways, same case", scB[scB.layer == "pathway_humann"], refB("pathway_humann"), pnames, pcls, pshort, "pathway_humann"))]
    fig = plt.figure(figsize=(T.DOUBLE_IN, 4.25))   # 108 mm tall
    top, rh, lh = 0.955, 0.64, 0.225   # radial axes top and height (equal aspect, 9.6 x 7.44 data units -> 2.72 in), list height; the sparse-layer summary sits between
    for k_, (title, sc, ref, cs, sample, names, classes, short, layer, pipe, sparse) in enumerate(panels):
        x0 = 0.005 + 0.5 * k_; y_rad = top - rh
        ax = fig.add_axes([x0, y_rad, 0.49, rh])
        d = table(sc, ref, cs, sample, names, classes, short); listed, r = radial(ax, d, title, len(cs), cls_min=12 if len(d) > 400 else 6)
        rows.append({"pipeline": pipe, "study": A_STUDY if pipe == "A" else B_STUDY, "sample": sample, "layer": layer, **r})
        s_title, s_sc, s_ref, s_names, s_cls, s_short, s_layer = sparse
        m = summary(table(s_sc, s_ref, cs, sample, s_names, s_cls, s_short), len(cs))
        rows.append({"pipeline": pipe, "study": A_STUDY if pipe == "A" else B_STUDY, "sample": sample, "layer": s_layer, **m})
        fig.text(x0 + 0.015, y_rad - 0.006, f"{s_title}: {m['low'] + m['high']} of {m['assessed']} outside, {m['expected_but_missing']} expected but missing;\n"
                 f"{'no deviation' if m['labelled'] == 0 else str(m['labelled']) + ' deviations'} shared with ≥ 25 % of the cohort's cases",
                 ha="left", va="top", fontsize=T.PT_MIN, color=T.INK["secondary"], linespacing=1.25)
        lx = fig.add_axes([x0 + 0.015, 0.012, 0.48, lh]); lx.axis("off")
        if listed:   # one column across the full panel width: id and full name, never truncated
            lx.text(0.0, 1.0, "\n".join(listed), ha="left", va="top", fontsize=T.PT_MIN, color=T.INK["primary"], transform=lx.transAxes, linespacing=1.3)
        else:
            lx.text(0.0, 1.0, "no KEGG ortholog is outside the healthy range\nin both this sample and ≥ 25 % of the cohort's cases", ha="left", va="top", fontsize=T.PT_MIN, color=T.INK["muted"], transform=lx.transAxes)
    T.save(fig, "fig13_gene_report", out_dir); plt.close(fig)
    new = pd.DataFrame(rows); rec = f"{P}/results/s14/example_reports.tsv"
    if os.path.exists(rec):
        old = pd.read_csv(rec, sep="\t")
        key = ["pipeline", "sample", "layer"]
        a = old.set_index(key)[DATA_COLS].sort_index(); b = new.set_index(key)[DATA_COLS].sort_index()
        if a.equals(b):
            return
    # results/ is frozen: a differing table is written beside the figure for review, never over the recorded file
    new.to_csv(os.path.join(out_dir, "example_reports_CHECK.tsv"), sep="\t", index=False); print("fig13: example counts differ from results/s14/example_reports.tsv; see example_reports_CHECK.tsv")


if __name__ == "__main__":
    make(sys.argv[1] if len(sys.argv) > 1 else f"{P}/figures")
