"""Figure 15 — how each healthy baseline is curated and built (double column). Two columns, one per resource: the
attrition from the public archive to the healthy-adult baseline (counts computed from the inventory tables), the held-out
sets reserved for evaluation, and a shared bottom row with the baseline-build steps. Pipelines by fill/outline (tokens.PIPE)."""
import os, sys, json
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import tokens as T

P = T.P


def counts():
    inv = pd.read_parquet(f"{P}/work/s0/inventory.parquet", columns=["analysis_id", "study_bioproject", "s0_status", "is_primary_analysis", "body_site_core"])
    ret = inv[inv.s0_status == "retained"]; gut = ret[ret.is_primary_analysis & (ret.body_site_core == "Gut")]
    qc = pd.read_parquet(f"{P}/work/s2/qc_per_analysis.parquet", columns=["analysis_id", "s2_status"]); qcok = qc[qc.s2_status == "retained"]
    ma = json.load(open(f"{P}/refs/gut-assembly-adult-global-v0.8-lenient/manifest.json")); mb = json.load(open(f"{P}/refs/gut-reads-adult-global-v0.2-lenient/manifest.json"))
    mA = pd.read_parquet(f"{P}/work/s12/A/meta.parquet"); hA = mA[mA.role == "heldout_healthy"].drop_duplicates("sample")
    S = pd.read_csv(f"{P}/work/s11/pipeB/samples.tsv", sep="\t", low_memory=False)
    inc = S[S.status == "included"]; hB = S[S.role == "heldout_healthy"]; cc = S[S.role == "case_control"]
    bal = pd.read_csv(f"{P}/results/s14/disease_strata_B_balance.tsv", sep="\t"); ev = bal.study.nunique(); nev = int((bal.n_case + bal.n_control).sum())   # studies with >= 15 cases and >= 15 controls scored
    first = open(f"{P}/results/s8/report_loso.md").readline()   # '# S8 — ... (11 studies, 2,294 in-stratum samples, 1,407 cases)'
    import re; a_cc = re.search(r"\((\d+) studies, ([\d,]+) in-stratum samples", first).groups()
    return dict(
        A=[f"{len(inv):,} MGnify v5 analyses, {inv.study_bioproject.nunique()} studies",
           f"{len(ret):,} human metagenomes with data\n(excluded: metatranscriptomes, non-human hosts,\nnegative and synthetic controls)",
           f"{len(gut):,} gut analyses, {gut.study_bioproject.nunique()} studies\n(one primary analysis per sample)",
           f"{len(qcok):,} pass assembly quality filters\n(≥ 5 Mb assembled, ≥ 20 of 52 marker genes)"],
        B=[f"{len(S):,} cMD3 MetaPhlAn 3 profiles, {S.study.nunique()} studies",
           f"adult stool samples: {len(inc):,} controls and {(S.status == 'NOT_CONTROL').sum():,} cases\n(excluded: {(S.status == 'OUT_OF_STRATUM').sum()} under 18, {(S.status == 'LOW_MAPPED_FRACTION').sum()} with < 90 % of the\nprofile resolvable to a family)"],
        A_out=[f"baseline:\n{ma['n_samples']:,} healthy adults,\n{ma['n_studies']} studies", f"held out:\n{len(hA)} healthy adults,\n{hA.study.nunique()} studies", f"case/control:\n{a_cc[0]} studies,\n{a_cc[1]} samples"],
        B_out=[f"baseline:\n{mb['n_samples']:,} healthy adults,\n{mb['n_studies']} studies", f"held out:\n{len(hB)} healthy adults,\n{hB.study.nunique()} studies", f"case/control:\n{ev} studies,\n{nev:,} samples"])


def box(ax, x, y, w, h, text, style, fs=T.PT_MIN, bold=False):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=0.08", linewidth=0.6, **style))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs, color=T.INK["primary"], weight="bold" if bold else "normal", linespacing=1.25)


def arrow(ax, x0, y0, x1, y1):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=6, lw=0.6, color=T.INK["secondary"], shrinkA=0, shrinkB=0))


def make(out_dir):
    T.print_profile(); c = counts()
    fig = plt.figure(figsize=(T.DOUBLE_IN, 4.3)); ax = fig.add_axes([0.005, 0.005, 0.99, 0.99]); ax.set_xlim(0, 20); ax.set_ylim(0, 12); ax.axis("off")
    plain = {"facecolor": T.INK["surface"], "edgecolor": T.INK["axis"]}
    for k, (key, title, x0, style) in enumerate([("A", "a  assembly pipeline: MGnify v5 cohort", 0.2, T.PIPE["assembly"]), ("B", "b  read pipeline: curatedMetagenomicData 3 cohort", 10.2, T.PIPE["read"])]):
        ax.text(x0, 11.7, title, fontsize=T.PT_TITLE, va="center", ha="left")
        steps = c[key]; y = 10.75; hh = 0.95 if key == "A" else 1.15
        ys = []
        for i, s in enumerate(steps):
            h = hh if i == 0 else (1.3 if key == "B" else 1.05)
            box(ax, x0 + 0.6, y - h, 8.4, h, s, plain); ys.append((y - h, y)); y = y - h - 0.45
        for (b0, t0), (b1, t1) in zip(ys[:-1], ys[1:]):
            arrow(ax, x0 + 4.8, b0, x0 + 4.8, t1)
        yo = 3.75; bw = 2.7
        for j, s in enumerate(c[key + "_out"]):
            bx = x0 + 0.6 + j * (bw + 0.15); box(ax, bx, yo, bw, 1.45, s, style if j == 0 else plain, bold=False)
            arrow(ax, x0 + 4.8, ys[-1][0], bx + bw / 2, yo + 1.45)
    # shared baseline build row
    ax.text(0.2, 3.05, "c  baseline build and validation (same steps for both pipelines)", fontsize=T.PT_TITLE, va="center", ha="left")
    steps = ["normalize with fixed rules:\nCLR over features in ≥ 5 %\nof samples; genes as\ncopies per genome",
             "per feature: prevalence by\nassembly-size or read-depth\nclass; percentiles 1–99,\nequal weight per study",
             "statistics only, no samples;\nbootstrap CIs over studies;\nversioned bundle",
             "validate: held-out cohorts and\neach study left out:\n~5 % outside the 2.5–97.5 range"]
    w = 4.55
    for i, s in enumerate(steps):
        x = 0.8 + i * (w + 0.25); box(ax, x, 0.35, w, 2.1, s, plain)
        if i:
            arrow(ax, x - 0.25, 1.4, x, 1.4)
    T.save(fig, "fig15_curation", out_dir); plt.close(fig)


if __name__ == "__main__":
    make(sys.argv[1] if len(sys.argv) > 1 else f"{P}/figures")
