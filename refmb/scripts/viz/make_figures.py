#!/usr/bin/env python3
"""Figure driver — render every figure twice (determinism), run lint on each PDF/SVG pair, write figures/LINT.json and
figures/README.md (lint table, file-to-manuscript mapping, colour conventions, captions). Exit non-zero if any figure fails
lint. The project root comes from REFMB_PROJECT (tokens.py); every number in a caption that describes the worked example is
computed (example_case.py, fig7_landscape.stats), the rest are quoted from the result files named in the figure scripts."""
import glob, importlib, json, os, re, shutil, sys, time
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lint, tokens as T

P = T.P; OUT = f"{P}/figures"; RERUN = f"{OUT}/_rerun"
FIGS = ["fig1_schematic", "fig2_composition", "fig4_calibration", "fig5_tree", "fig7_landscape", "fig8_performance", "fig9_strata", "fig10_disease_atlas", "fig11_container", "fig12_disease_strata", "fig13_gene_report", "fig14_sample_card", "fig15_curation", "fig16_function_shifts"]


def newest_draft():
    """The newest manuscript draft by version number (results/paper/refmb_draft_v*.md), project-relative; None if absent."""
    ds = glob.glob(f"{P}/results/paper/refmb_draft_v*.md")
    key = lambda p: tuple(int(x) for x in re.search(r"_v([0-9.]+)\.md$", p).group(1).split("."))
    return os.path.relpath(max(ds, key=key), P) if ds else None


def manuscript_mapping():
    """figure file -> figure number, parsed from the image lines '![Figure N](figures_vX/<file>.png)' of the newest draft;
    a figure the draft does not include is labelled so."""
    draft = newest_draft(); m = {}
    if draft:
        for num, stem in re.findall(r"!\[Figure ([0-9S]+)\]\([^)]*?/(fig\d+_[a-z_]+)\.png\)", open(f"{P}/{draft}").read()):
            m.setdefault(stem, f"Fig. {num}")
    return {n: m.get(n, "not in the manuscript") for n in FIGS}, draft


def manuscript_order(mapping):
    """FIGS sorted by manuscript number: main figures, then supplementary (S1, S2, ...), then figures not in the draft."""
    def key(n):
        m = re.match(r"Fig\. (S?)(\d+)$", mapping[n])
        return (2, 0, n) if not m else (1 if m.group(1) else 0, int(m.group(2)), n)
    return sorted(FIGS, key=key)


MANUSCRIPT, DRAFT = manuscript_mapping(); ORDERED = manuscript_order(MANUSCRIPT)
COLOUR_CONVENTIONS = [
    ("blue / orange (deviation scale)", "direction of a deviation from the healthy range only: blue = below (low), orange = above (high); deep step beyond the 1st/99th percentile or a shift of 20 points or more, light step the milder class; grey = within. fig1 (b), fig5, fig10, fig12, fig13; fig14 (orange diamond = above the healthy 90th percentile); fig4 (orange = cohort above the pre-set 8 % limit)."),
    ("four categorical colours", "the four feature sets of the classification benchmark only (fig8): percentiles, raw CLR, alpha diversity, GMHI approximation."),
    ("sequential blues", "an ordinal magnitude: baseline samples per hex (fig7) and the assembly-size or read-depth class of a study's samples (fig2 a, d)."),
    ("greys and marker shape", "everything else: assembly pipeline = filled light-grey box or filled bar, read pipeline = outlined (dashed) box or outlined bar (fig1, fig2); family layer = filled circle, KO layer = open square (fig9, fig11); paired-end = circle, single-end = triangle (fig11); medians, quantile bars and reference lines in ink greys."),
]


def captions():
    import example_case as EC, fig7_landscape as F7
    sys.path.insert(0, P); from refmb.paths import bundle_A
    c = EC.counts(); mA = json.load(open(f"{bundle_A()}/manifest.json")); mB = json.load(open(f"{P}/refs/gut-reads-adult-global-v0.2-lenient/manifest.json")); s7 = F7.stats()
    nA, sA, nB, sB = mA["n_samples"], mA["n_studies"], mB["n_samples"], mB["n_studies"]
    # fig13 / fig14 example counts from the recorded table (results/s14/example_reports.tsv), never typed
    ex = pd.read_csv(f"{P}/results/s14/example_reports.tsv", sep="\t").set_index("layer"); ex["outside"] = ex["low"] + ex["high"]
    koA, mod, koB, pw = ex.loc["ko_eggnog"], ex.loc["module"], ex.loc["ko_humann"], ex.loc["pathway_humann"]
    # fig7: outside-contour counts and the richness statement (only when |rho| >= RHO_LABEL pooled and within each cohort)
    out_h, out_c = round(s7["frac_outside95_heldout"] * s7["n_heldout"]), round(s7["frac_outside95_cdi"] * s7["n_cdi"])
    rich1 = min(abs(s7["rho_PC1_n_assessed"]), abs(s7["rho_PC1_heldout"]), abs(s7["rho_PC1_cdi"])) >= F7.RHO_LABEL
    rich2 = min(abs(s7["rho_PC2_n_assessed"]), abs(s7["rho_PC2_heldout"]), abs(s7["rho_PC2_cdi"])) >= F7.RHO_LABEL
    rho_txt = (f"Over the re-scored cohort samples PC1 {'tracks' if rich1 else 'does not consistently track'} the number of assessed families (Spearman ρ = {s7['rho_PC1_n_assessed']:.2f} pooled; "
               f"{s7['rho_PC1_heldout']:.2f} and {s7['rho_PC1_cdi']:.2f} within the two cohorts) and PC2 {'tracks it' if rich2 else 'does not consistently'} (ρ = {s7['rho_PC2_n_assessed']:.2f} pooled; "
               f"{s7['rho_PC2_heldout']:.2f} and {s7['rho_PC2_cdi']:.2f}).")
    return {
        "fig1_schematic": f"**One sample, two pipelines, one kind of report.** (a) A user's raw reads go through the assembly pipeline (filled boxes; pinned container reproducing MGnify v5) or the read pipeline (outlined boxes; MetaPhlAn 3.0.14 reproducing curatedMetagenomicData 3, HUMAnN 3 optional), are normalized with the rules stored in the matching healthy baseline ({nA:,} healthy adults from {sA} studies; {nB:,} from {sB}) and are scored against it. The report gives every feature's percentile among healthy adults with a call (within / low / high), lists expected-but-missing features, and gives the sample's share of features outside the healthy range, tested against the about 5 % expected. Baselines ship as statistics only, without sample-level data. (b) The user's view: nine rows of the family-level report of one C. difficile infection sample ({c['n_assessed']} assessed families, {c['n_missing']} expected but missing), each placed at its percentile among healthy adults with the healthy range (2.5th–97.5th) shaded; blue = below the range, grey = within, hollow = expected but missing.",
        "fig2_composition": f"**Both baselines are dominated by a few large European and North American studies, plus two Tanzanian studies (372 adults) in the assembly baseline; curated metadata adds health status, age and sex for many samples.** (a, d) Healthy adults per study, stacked by assembly-size class (assembly pipeline, {nA:,} adults, {sA} studies) or read-depth class (read pipeline, {nB:,} adults, {sB} studies). (b, e) Samples by country (filled bar = assembly pipeline, outlined bar = read pipeline; country codes resolved to names, so Tanzania is one bar). (c) Share of 9,195 gut samples with each metadata field after curation (light grey) and after joining curatedMetagenomicData 3 metadata (dark grey).",
        "fig4_calibration": "**Independent healthy cohorts fall near the expected 5 % of features outside the healthy range on most layers, with named exceptions.** (a, c) Mean share of assessed features outside the healthy range (2.5th–97.5th percentile) for three healthy cohorts per pipeline that never entered the baseline (marker shape = cohort, sample counts in the legend); orange = cohort above the pre-set 8 % limit (assembly: PRJEB62833 on the gene layers; read: YeZ_2018, n = 43, on every layer). (b, d) Each baseline study scored against a baseline built without it (assembly: 26 studies; read: 19); small dots = studies, large dot = median, bar = interquartile range, arrow = study beyond the axis. Grey band 3–8 %, line = 5 %.",
        "fig5_tree": f"**One C. difficile infection case has {c['n_low']} of {c['n_assessed']} assessed families below the healthy range and {EC.none_or(c['n_high'])} above.** (a) Each point is one of the {c['n_families']} bacterial families carried by at least half of healthy adults; its distance from the centre is this case's percentile among healthy adults (baseline built without the study), the grey ring is the healthy range (2.5th–97.5th); blue = low, orange = high, hollow = expected but missing ({c['n_missing']} families). The outer ring marks families outside the range in at least 25 % (light) or 50 % (deep) of the study's 56 cases; labelled families are outside both in this case and in at least 25 % of cases. (b) How a percentile is read from one family's baseline distribution (Lachnospiraceae, joined to its leaf by the dotted line). Families ordered phylum > class > order.",
        "fig7_landscape": f"**The independent healthy cohort lies inside the healthy cloud ({'none' if out_h == 0 else out_h} of its {s7['n_heldout']} samples outside the 95 % contour); {out_c} of the {s7['n_cdi']} C. difficile infection cases ({100 * s7['frac_outside95_cdi']:.0f} %) lie outside it.** Hexbin density of the {s7['n_base']:,} baseline samples on the first two principal components of family CLR abundance, drawn from the coordinates stored in the bundle, with contours enclosing the densest 50 % (solid) and 95 % (dashed) of baseline samples; an independent healthy cohort (PRJEB62833, n = {s7['n_heldout']}, open circles) and the C. difficile infection cases (PRJEB26165, n = {s7['n_cdi']}, filled diamonds) were scored with refmb 0.8.1 and projected onto the same axes (an in-pool sample re-projects onto its stored coordinate to within 1e-6). {rho_txt}",
        "fig8_performance": "**Percentiles transfer between studies modestly better than raw abundances; the gain is clearest for C. difficile infection and colorectal cancer.** (a) Assembly pipeline: AUROC in each held-out study for a random forest trained on the other 10, with 95 % bootstrap CI, for four feature sets; the baseline is rebuilt without the held-out study. (b) Mean AUROC over the 11 studies (percentiles 0.618 vs raw CLR 0.581, 9 of 11 studies won, p = 0.021), with the within-study cross-validated AUROC as a ceiling (black tick). (c) Read pipeline, taxonomy features: percentiles vs raw CLR in 14 held-out studies. (d) Mean AUROC by condition of the held-out study (all: 0.690 vs 0.664, p = 0.052; colorectal cancer: 0.769 vs 0.715, p = 0.004). Vertical line = AUROC 0.5. GMHI (genus approx.) = genus-level approximation of the Gut Microbiome Health Index.",
        "fig9_strata": "**The healthy baseline holds across age, sex and BMI groups; it does not fit African cohorts.** Median share of features outside the healthy range for healthy adults by age, sex, BMI, lifestyle and region; (a) assembly pipeline, (b) read pipeline; filled circle = family layer, open square = KO layer. Every sample is scored against a baseline built without its study. Dashed line = the pre-set 8 % limit; sample and study counts at right. African samples reach 20 % of families outside in the assembly pipeline.",
        "fig10_disease_atlas": "**Colorectal cancer studies share a block of genera that are lower in cases; the two inflammatory bowel disease studies share a block that is higher.** Each cell is the shift in median percentile of a genus, cases minus controls of the same study, in percentile points (blue lower in cases, orange higher; five-step scale). Dots (white on the deep steps) mark shifts of at least 10 points with q ≤ 0.05 (Mann–Whitney, Benjamini–Hochberg within layer). Columns are studies grouped by condition, separated by lines; blank = genus not assessed. (a) Assembly pipeline, 10 studies. (b) Read pipeline, 15 studies. Genera shown are the 22 with the most significant studies.",
        "fig11_container": "**The assembly pipeline's container reproduces MGnify's analysis of the same reads: gene layers agree at Spearman 0.98–1.00, taxonomy at 0.86–0.88.** (a) Assembled length, container vs MGnify, for 61 public samples (51 within 2 %; circle = paired-end reads, triangle = single-end). (b) Spearman correlation of feature percentiles between the two analyses, per layer (bar = median). (c) KO (open square) and family (filled circle) correlation by read layout and the metaSPAdes version MGnify used; agreement is lowest where that version is unknown.",
        "fig12_disease_strata": "**Colorectal cancer shifts keep their direction within age, sex and BMI groups; North America and East Asia show weaker shifts, each from one or two studies.** Read pipeline. (a) Median shift in percentile, cases minus same-study controls, of the 18 taxa consistent across studies, within each group. (b) Share of all 203 cross-study-consistent features (taxa, KOs, pathways) that keep their direction within the group. (c) Median absolute shift of those features. Numbers in brackets = contributing studies.",
        "fig13_gene_report": f"**The gene layers are read the same way as the family layer, and say both when a sample deviates and when it does not.** Radial reports of KEGG orthologs for (a) the C. difficile case of the family report (assembly pipeline, eggNOG annotation) and (b) one colorectal cancer case (read pipeline, FengQ_2015, HUMAnN 3). Conventions as in the family report; features ordered by functional class (arcs); the outer ring is a histogram of features per 5° sector that are low (blue) or high (orange) in at least 25 % of the cohort's cases. Numbered features are outside the healthy range both in this sample and in at least 25 % of cases; the eight most frequent are listed with KO identifier, full name and the share of cases. The line under each ring summarises the same sample's sparse layer (KEGG modules; MetaCyc pathways). In (a) {koA['outside']} of {koA['assessed']:,} assessed orthologs ({100 * koA['outside'] / koA['assessed']:.0f} %) are outside and {koA['expected_but_missing']} expected genes are missing; {mod['outside']} of {mod['assessed']} modules is outside. In (b) {koB['outside']} of {koB['assessed']:,} orthologs ({100 * koB['outside'] / koB['assessed']:.0f} %) and {pw['outside']} of {pw['assessed']} pathways are outside, and no deviation is shared with a quarter of the cohort.",
        "fig15_curation": "**Each baseline is curated from a public resource by fixed rules, and held-out sets are reserved before any evaluation.** (a) MGnify v5 analyses to the assembly baseline; (b) curatedMetagenomicData 3 profiles to the read baseline; counts at each step; held-out healthy cohorts and case/control studies never enter a baseline. (c) Build and validation steps shared by both pipelines.",
        "fig16_function_shifts": "**Function features shifted in the same direction in several colorectal cancer cohorts are candidates beyond the known taxa.** Read pipeline; MetaCyc pathways and KEGG orthologs shifted (at least 10 percentile points, q ≤ 0.05) in the same direction in at least three of nine cohorts, plus the five strongest shifted in two (branched-chain amino acid biosynthesis); cell colour = shift in median percentile, cases minus same-study controls.",
        "fig14_sample_card": f"**The sample-level numbers of a report, read against healthy adults: the C. difficile case stands out mainly by its missing genes, the colorectal cancer case is unremarkable.** (a, b) The C. difficile infection case of the family report (assembly pipeline); (c, d) the colorectal cancer case of the gene report (read pipeline). Left: share of all assessed features outside the healthy range (over all assessed features, so smaller than the common-family share quoted in the text); right: number of expected-but-missing features. Grey bar = 10th–90th percentile and tick = median of healthy adults scored against a baseline built without their study (assembly: 2,174 samples, 29 studies; read: 6,811 samples, 22 studies); diamond = the example sample, orange when above the healthy 90th percentile; small grey dots, 160 of those healthy samples per layer; values beyond the axis are drawn at its edge. The C. difficile case's values come from the same scoring as the family and gene reports ({koA['expected_but_missing']} missing KOs on the eggNOG layer), the colorectal cancer case's from the out-of-sample scoring. Healthy samples spread widely on the share outside (90th percentile about 20 % on taxonomy), so this number alone separates little; the C. difficile case is above the healthy median on taxonomy but inside the spread, while its counts of missing KOs and Pfam families exceed the healthy 90th percentile.",
    }


CAPTIONS = captions()


def readme(results, ok):
    L = [f"# Figures — generated {time.strftime('%Y-%m-%d %H:%M UTC', time.gmtime())}; lint {'PASS' if ok else 'FAIL'}\n",
         "Print profile: exact IEEE widths (88.9 mm single, 181.9 mm double), DejaVu Sans shipped from resources/fonts, 7 pt body / 6 pt floor, fonts embedded (Type 42), real text in SVG, token palette only, deterministic (two renders byte-identical). Validation of the palette: scripts/viz/tokens_validation.txt. Scripts read the project root from `REFMB_PROJECT` (when unset, the directory two levels above scripts/viz/) and fail with a plain error if it is not a project directory. The text-in-page lint check needs the pymupdf package.\n",
         "| figure | width (mm) | min font (pt) | fonts embedded | deterministic | lint |", "|---|---:|---:|---|---|---|"]
    for r in results:
        L.append(f"| {r['file']} | {r['width_mm']} | {r['min_font_pt']} | {'yes' if all(r['fonts'].values()) else 'NO'} | {'yes' if not any('non-deterministic' in p for p in r['problems']) else 'NO'} | {'PASS' if r['pass'] else 'FAIL: ' + '; '.join(r['problems'])} |")
    L += ["\n## Figure files and manuscript numbers\n", f"Numbering is read from the image references of the newest manuscript draft ({DRAFT}); rows are in manuscript order. Script numbers fig3 and fig6 are not produced (cross-pipeline calibration; metabolic network), both listed as future work.\n",
          "| figure file | manuscript figure |", "|---|---|"] + [f"| {n} | {MANUSCRIPT[n]} |" for n in ORDERED]
    L += ["\n## Colour conventions\n", "One meaning per colour across the set (the same table is a comment block in scripts/viz/tokens.py):\n", "| colour | meaning and where it is used |", "|---|---|"] + [f"| {a} | {b} |" for a, b in COLOUR_CONVENTIONS]
    L.append("\n## Captions\n")
    for n in ORDERED:
        L.append(f"- **{n}** ({MANUSCRIPT[n]}) — {CAPTIONS[n]}")
    L.append("\nCaptions open with the bold takeaway sentence; counts for the worked example are computed by scripts/viz/example_case.py. The HTML summary is figures/summary.html (scripts/viz/make_summary.py).\n")
    return L


def main():
    if "--readme-only" in sys.argv:   # rebuild README.md and LINT.json-derived table without re-rendering
        results = json.load(open(f"{OUT}/LINT.json")); ok = all(r["pass"] for r in results)
        open(f"{OUT}/README.md", "w").write("\n".join(readme(results, ok))); print("README rebuilt from LINT.json"); return
    t0 = time.time()
    if os.path.isdir(RERUN):
        shutil.rmtree(RERUN)
    for name in FIGS:
        m = importlib.import_module(name)
        m.make(OUT); m.make(RERUN)
        print(f"rendered {name} ({time.time()-t0:.0f}s)", flush=True)
    results = [lint.lint(f"{OUT}/{n}.pdf", f"{OUT}/{n}.svg", f"{RERUN}/{n}.svg") for n in ORDERED]
    json.dump(results, open(f"{OUT}/LINT.json", "w"), indent=1)
    ok = all(r["pass"] for r in results)
    L = readme(results, ok)
    open(f"{OUT}/README.md", "w").write("\n".join(L)); print("\n".join(L[:4 + len(results)]))
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
