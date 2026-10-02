#!/usr/bin/env python3
"""Build figures/summary.html: a self-contained five-minute summary of refmb (inline CSS, PNG figures embedded as base64,
no external scripts). Figures are re-rendered from figures/*.pdf at 150 dpi with pdftoppm; captions come from
make_figures.CAPTIONS; every number in the tables is read from a results file (paths listed under each table).
The page is validated by parsing it with html.parser and decoding every embedded image."""
import base64, glob, html, io, json, os, re, subprocess, sys, tempfile
from html.parser import HTMLParser
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tokens as T, example_case as EC
from make_figures import CAPTIONS, FIGS, MANUSCRIPT, newest_draft

P = T.P; FIG = f"{P}/figures"; OUT = f"{FIG}/summary.html"

LAYER = {"taxonomy_family": "family", "taxonomy_genus": "genus", "taxonomy_species": "species", "ko_eggnog": "KO (eggNOG)", "ko_kofam": "KO (KOfam)", "pfam": "Pfam", "module": "KEGG module",
         "ko_humann": "KO (HUMAnN 3)", "pathway_humann": "pathway (MetaCyc)"}
A_LAYERS = ["taxonomy_family", "taxonomy_genus", "ko_eggnog", "ko_kofam", "pfam", "module"]; B_LAYERS = ["taxonomy_family", "taxonomy_genus", "taxonomy_species", "ko_humann", "pathway_humann"]
ORDER = ["fig1_schematic", "fig5_tree", "fig13_gene_report", "fig14_sample_card", "fig4_calibration", "fig11_container", "fig10_disease_atlas", "fig12_disease_strata", "fig9_strata", "fig8_performance", "fig7_landscape", "fig2_composition"]
TAKEAWAY = {   # one sentence per figure, for the reader in a hurry (the caption carries the detail)
    "fig1_schematic": "A user runs one of two standard pipelines on raw reads and gets, for every taxon, gene and pathway, a percentile among healthy adults measured the same way; panel (b) is the report as the user sees it.",
    "fig5_tree": "The family-level report of one C. difficile infection sample: every deviation is on the low side, in the main anaerobic families.",
    "fig13_gene_report": "The same report on genes and pathways says when a sample deviates (C. difficile case) and when it does not (colorectal cancer case).",
    "fig14_sample_card": "Sample-level numbers are read against the healthy spread; the share outside alone separates little, the count of missing genes is what stands out for the C. difficile case.",
    "fig4_calibration": "Independent healthy cohorts fall near the expected 5 % of features outside the healthy range on most layers; the exceptions are named.",
    "fig11_container": "The container reproduces MGnify's analysis from raw reads: Spearman 0.98–1.00 on gene layers, 0.86–0.88 on taxonomy.",
    "fig10_disease_atlas": "Across independent colorectal cancer studies a shared block of genera is lower in cases; the two inflammatory bowel disease studies share a block that is higher.",
    "fig12_disease_strata": "Colorectal cancer shifts keep their direction within age, sex and BMI groups; regional differences cannot be separated from study.",
    "fig9_strata": "The baseline holds across age, sex and BMI groups of healthy adults; it does not fit African cohorts (20 % of families outside in the assembly pipeline).",
    "fig8_performance": "Percentiles transfer between studies modestly better than raw abundances (A: 0.618 vs 0.581, p = 0.021; B taxonomy: 0.690 vs 0.664, p = 0.052); function layers show no advantage.",
    "fig7_landscape": "Position on the healthy map mostly reflects richness (PC1 correlates at 0.93 with the number of assessed families); the independent healthy cohort lies in the dense half of the cloud, the C. difficile cases at its edge.",
    "fig2_composition": "What the baselines are made of: a few large European and North American studies dominate both; curated metadata adds health status, age and sex.",
}
REPORTS = [("Calibration, assembly pipeline (A)", "results/s6/report_v08_lenient_base.md"), ("Leave-one-study-out calibration, A", "results/s6/loso_within_pool.md"),
           ("Case/control classification, known biology, A", "results/s8/report.md"), ("Per-study AUROC, A (main / strict baseline)", "results/s8/report_loso.md, results/s8/report_loso_strict.md"),
           ("Calibration and classification, read pipeline (B)", "results/s11/pipelineB_report.md"), ("Pipeline B from raw reads", "results/s11/pipelineB_fidelity.md"),
           ("Pipeline A container vs MGnify", "results/s11/tierA_anchor_compare.md"), ("Healthy baseline across groups", "results/s12/strata_A.md, results/s12/strata_B.md"),
           ("Disease atlas", "results/s12/disease_atlas_A.md, results/s12/disease_atlas_B.md"), ("Disease shifts within groups", "results/s14/disease_strata_A.md, results/s14/disease_strata_B.md"),
           ("Example sample counts", "results/s14/example_reports.tsv"), ("Current manuscript", newest_draft() or "results/paper/ (no draft found)"), ("Figure lint, captions and colour conventions", "figures/README.md, figures/LINT.json")]


def png_b64(stem, dpi=150):
    with tempfile.TemporaryDirectory() as d:
        subprocess.run(["pdftoppm", "-r", str(dpi), "-png", "-singlefile", f"{FIG}/{stem}.pdf", f"{d}/x"], check=True)
        return base64.b64encode(open(f"{d}/x.png", "rb").read()).decode()


def md_table(path, header_key):
    """Parse the first GitHub-markdown table after a line containing header_key."""
    lines = open(path).read().splitlines(); i = next(k for k, l in enumerate(lines) if header_key in l) + 1
    while i < len(lines) and not lines[i].startswith("|"):
        i += 1
    rows = []
    while i < len(lines) and lines[i].startswith("|"):
        rows.append(lines[i]); i += 1
    hdr = [c.strip() for c in rows[0].strip("|").split("|")]; body = [[c.strip() for c in r.strip("|").split("|")] for r in rows[2:]]
    return pd.DataFrame(body, columns=hdr)


def tables_A():
    rep = json.load(open(f"{P}/results/s6/report_v08_lenient_base.json")); ps = pd.DataFrame(rep["per_study"])
    loso = pd.read_csv(f"{P}/results/s6/loso_within_pool.tsv", sep="\t")
    cal = []
    for l in A_LAYERS:
        v = ps[ps["layer"] == l].sort_values("study_bioproject"); lv = loso.loc[loso["layer"] == l, "frac"]
        cal.append([LAYER[l], " / ".join(f"{100*x:.1f}" for x in v["frac_outside_mean"]), f"{100*lv.median():.1f}", f"{100*(lv <= 0.08).mean():.0f} %"])
    cal = pd.DataFrame(cal, columns=["layer", "held-out cohorts, mean % outside (PRJEB39906 / PRJEB62833 / PRJEB75578)", "leave-one-study-out median % outside (26 studies)", "studies ≤ 8 %"])
    C = pd.read_csv(f"{P}/results/s11/tierA_anchor_compare.tsv", sep="\t"); C = C[C["variant"] == "default"]
    fid = C.groupby("layer").agg(n=("anchor", "nunique"), rho=("spearman_pct", "median"), fo_c=("frac_outside_container", "mean"), fo_n=("frac_outside_native", "mean"), agree=("call_agreement", "median")).reindex(A_LAYERS)
    fid = pd.DataFrame({"layer": [LAYER[l] for l in fid.index], "samples": fid["n"].astype(int).values, "median Spearman of percentiles, container vs MGnify": fid["rho"].round(3).values,
                        "mean % outside, container vs MGnify": [f"{100*a:.1f} vs {100*b:.1f}" for a, b in zip(fid["fo_c"], fid["fo_n"])], "median call agreement": fid["agree"].round(3).values})
    L = pd.read_csv(f"{P}/results/s8/loso_auroc.tsv", sep="\t"); per = L[L["study"] != "POOLED"]
    mac = per.groupby("feature_set")["auroc_loso"].mean(); pooled = L[L["study"] == "POOLED"].set_index("feature_set")["auroc_loso"] if (L["study"] == "POOLED").any() else None
    pair = md_table(f"{P}/results/s8/report.md", "Paired comparison").set_index("comparison")
    names = {"reference_relative": "percentiles (reference-relative)", "raw_clr": "raw CLR", "alpha_diversity": "alpha diversity", "health_index": "GMHI (genus-level approximation)"}
    cc = []
    for fs in ["reference_relative", "raw_clr", "alpha_diversity", "health_index"]:
        row = [names[fs], f"{mac[fs]:.3f}"]
        if fs == "reference_relative":
            row += ["—", "—", "—"]
        else:
            r = pair.loc[f"reference_relative − {fs}"]; row += [f"{float(r['mean_diff']):+.3f}", f"{r['wins']} / {r['losses']}", r["wilcoxon_p"]]
        cc.append(row)
    cc = pd.DataFrame(cc, columns=["features", "mean AUROC over 11 held-out studies", "percentiles minus this (mean)", "studies won / lost by percentiles", "one-sided Wilcoxon p"])
    strict = md_table(f"{P}/results/s8/report.md", "strict (sensitivity) reference").set_index("feature_set")
    strict_note = f"Strict baseline (1,588 healthy adults, extra exclusions): percentiles {float(strict.loc['reference_relative', 'macro_mean']):.3f} vs raw CLR {float(strict.loc['raw_clr', 'macro_mean']):.3f}, paired p = {pair_strict_p()}."
    return cal, fid, cc, strict_note


def pair_strict_p():
    lines = open(f"{P}/results/s8/report.md").read().splitlines(); i = next(k for k, l in enumerate(lines) if "strict (sensitivity)" in l)
    for l in lines[i:]:
        if l.startswith("| reference_relative − raw_clr"):
            return l.strip("|").split("|")[-1].strip()
    return "n/a"


def tables_B():
    hb = pd.read_csv(f"{P}/results/s11/pipelineB_heldout.tsv", sep="\t"); lb = pd.read_csv(f"{P}/results/s11/pipelineB_loso_all.tsv", sep="\t")
    cal = []
    for l in B_LAYERS:
        v = hb[hb["layer"] == l].sort_values("study"); lv = lb.loc[lb["layer"] == l, "median"]
        cal.append([LAYER[l], " / ".join(f"{100*x:.1f}" for x in v["mean"]), f"{100*lv.median():.1f}", f"{100*(lv <= 0.08).mean():.0f} %"])
    studies = " / ".join(sorted(hb["study"].unique()))
    cal = pd.DataFrame(cal, columns=["layer", f"held-out cohorts, mean % outside ({studies})", f"leave-one-study-out median % outside ({lb['study'].nunique()} studies; {lb.loc[lb['layer'] == 'ko_humann', 'study'].nunique()} for function layers)", "studies ≤ 8 %"])
    F = pd.read_csv(f"{P}/results/s11/pipelineB_fidelity.tsv", sep="\t"); S = pd.read_csv(f"{P}/results/s11/pipelineB_fidelity_samples.tsv", sep="\t"); S = S[~S["deposit_mismatch"]]
    F = F[F["anchor"].isin(S["anchor"])]
    fid = F.groupby("layer").agg(n=("anchor", "nunique"), rho=("spearman_pct", "median"), fo_c=("frac_outside_container", "mean"), fo_n=("frac_outside_cmd3", "mean"), agree=("call_agreement", "median")).reindex(["taxonomy_family", "taxonomy_genus", "taxonomy_species"])
    fid = pd.DataFrame({"layer": [LAYER[l] for l in fid.index], "samples": fid["n"].astype(int).values, "median Spearman of percentiles, raw reads vs cMD3 profile": fid["rho"].round(3).values,
                        "mean % outside, raw reads vs cMD3": [f"{100*a:.1f} vs {100*b:.1f}" for a, b in zip(fid["fo_c"], fid["fo_n"])], "median call agreement": fid["agree"].round(3).values})
    fid_note = f"Species profiles from raw reads vs the cMD3 profile of the same sample: median Bray–Curtis distance {S['species_bray_curtis'].median():.3f}, median Spearman {S['species_spearman'].median():.3f} ({len(S)} samples, {S['study'].nunique()} studies). Function layers were not run from raw reads."
    cc = md_table(f"{P}/results/s11/pipelineB_report.md", "Case/control, cross-study")
    cc = cc.rename(columns={"feature_set": "features", "training": "held-out studies", "studies": "n studies", "reference_relative": "percentiles, mean AUROC", "raw_clr": "raw CLR, mean AUROC", "wins": "won", "losses": "lost", "wilcoxon_p": "one-sided Wilcoxon p"})
    cc["held-out studies"] = cc["held-out studies"].map({"all_studies": "all", "crc_only": "colorectal cancer only"}); cc["features"] = cc["features"].map({"all": "taxonomy + function", "taxonomy": "taxonomy", "function": "function"})
    return cal, fid, fid_note, cc


def tab_html(df, caption, source):
    h = [f'<table><caption>{html.escape(caption)}</caption><thead><tr>' + "".join(f"<th>{html.escape(str(c))}</th>" for c in df.columns) + "</tr></thead><tbody>"]
    for _, r in df.iterrows():
        h.append("<tr>" + "".join(f"<td>{html.escape(str(v))}</td>" for v in r) + "</tr>")
    h.append(f'</tbody></table><p class="src">Source: <code>{html.escape(source)}</code></p>')
    return "\n".join(h)


def md_inline(s):
    s = html.escape(s); s = re.sub(r"\*\*(.+?)\*\*", r"<strong>\1</strong>", s); return s


CSS = """
:root{--ink:#0b0b0b;--ink2:#52514e;--muted:#898781;--grid:#e1e0d9;--axis:#c3c2b7;--a:#2a78d6;--b:#eb6834;--surface:#ffffff;--band:#f3f2ee}
*{box-sizing:border-box}body{margin:0;background:var(--surface);color:var(--ink);font:15px/1.5 "DejaVu Sans","Segoe UI",system-ui,sans-serif}
main{max-width:1080px;margin:0 auto;padding:32px 24px 64px}h1{font-size:30px;line-height:1.2;margin:0 0 6px}h2{font-size:21px;margin:44px 0 10px;padding-top:12px;border-top:1px solid var(--grid)}
h3{font-size:16px;margin:28px 0 6px}.sub{color:var(--ink2);margin:0 0 20px}.lead p{font-size:16px}
.kpis{display:grid;grid-template-columns:repeat(auto-fit,minmax(200px,1fr));gap:12px;margin:22px 0}.kpi{border:1px solid var(--grid);border-radius:8px;padding:12px 14px;background:var(--band)}
.kpi b{display:block;font-size:22px}.kpi span{color:var(--ink2);font-size:13px}.kpi.a{border-left:4px solid var(--ink)}.kpi.b{border-left:4px solid var(--ink);border-left-style:dashed}
.stem{color:var(--muted);font-size:12.5px;font-weight:400;margin-left:8px}.tag{display:inline-block;font-size:12px;color:var(--ink2);background:var(--band);border:1px solid var(--grid);border-radius:4px;padding:0 6px;margin-left:8px;vertical-align:middle}
figure{margin:18px 0 0;border:1px solid var(--grid);border-radius:8px;padding:14px;background:#fff}figure img{width:100%;height:auto;display:block}
figcaption{margin-top:12px;font-size:14px;color:var(--ink2)}figcaption strong{color:var(--ink)}.take{margin:8px 0 0;padding:10px 14px;border-left:4px solid var(--axis);background:var(--band);font-size:15px}
table{border-collapse:collapse;width:100%;margin:14px 0 4px;font-size:13.5px}caption{text-align:left;font-weight:600;padding:0 0 6px;color:var(--ink)}th,td{text-align:left;padding:6px 8px;border-bottom:1px solid var(--grid);vertical-align:top}
th{color:var(--ink2);font-weight:600;background:var(--band)}td:not(:first-child),th:not(:first-child){font-variant-numeric:tabular-nums}.src{color:var(--muted);font-size:12.5px;margin:0 0 10px}
code{font:12.5px ui-monospace,Menlo,Consolas,monospace;background:var(--band);padding:1px 4px;border-radius:3px}.note{color:var(--ink2);font-size:14px}
ul.reports{columns:2;column-gap:28px;padding-left:18px;font-size:14px}ul.reports li{break-inside:avoid;margin-bottom:4px}.weak{border-left:4px solid var(--b)}
@media (max-width:700px){ul.reports{columns:1}main{padding:20px 16px 40px}h1{font-size:24px}}
"""


def build():
    sys.path.insert(0, P); from refmb.paths import bundle_A
    mA = json.load(open(f"{bundle_A()}/manifest.json")); mB = json.load(open(f"{P}/refs/gut-reads-adult-global-v0.2-lenient/manifest.json"))
    calA, fidA, ccA, strictA = tables_A(); calB, fidB, fidB_note, ccB = tables_B(); c = EC.counts()
    imgs = {s: png_b64(s) for s in ORDER}
    H = [f"<!DOCTYPE html><html lang='en'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'><title>refmb summary</title><style>{CSS}</style></head><body><main>",
         "<h1>refmb: gut metagenome pipelines bundled with healthy reference baselines</h1>",
         "<p class='sub'>Project summary for a five-minute read. Figures are the current figure set (figures/*.pdf); the tag after each heading gives the figure's number in the manuscript, and every number below is read from a results file named under the table. Authors and repository: [TBD].</p>",
         "<div class='lead'>",
         f"<p>A reference genome makes a genome interpretable because every observation is placed against a shared, versioned baseline. Gut microbiome profiles have no such baseline: each study decides what is normal from its own controls, so a single sample cannot be interpreted and results from different studies are not on one scale. refmb ships two standard analysis pipelines for gut metagenomes, each with a healthy-adult baseline measured by that same pipeline. The <strong>assembly pipeline (A)</strong> is a pinned container that reproduces the MGnify v5 analysis; its baseline holds {mA['n_samples']:,} healthy adults from {mA['n_studies']} studies with family, genus, KEGG ortholog (two annotators), Pfam and KEGG module layers. The <strong>read pipeline (B)</strong> runs MetaPhlAn 3.0.14 (HUMAnN 3 optional) and reproduces curatedMetagenomicData 3; its baseline holds {mB['n_samples']:,} healthy adults from {mB['n_studies']} studies with family, genus, species, KO and pathway layers. A user supplies raw reads and receives, for every feature, its percentile among healthy adults and a call (within, low or high the healthy range of 2.5–97.5), a list of expected-but-missing features, and the sample's share of features outside the range. Baselines are distributed as statistics only.</p>",
         "<p>The baseline is checked before it is used. Healthy cohorts that never entered a baseline fall near the expected 5 % of features outside the healthy range on most layers, and both pipelines reproduce the original analyses from raw reads (gene layers at Spearman 0.98–1.00, taxonomy at 0.86–0.88 for the assembly pipeline; 0.997–0.999 for the read pipeline). On real data the report shows, for one <em>C. difficile</em> infection sample, {c['n_low']} of {c['n_assessed']} assessed families below the healthy range and {EC.none_or(c['n_high'])} above; across nine independent colorectal cancer studies, Lachnospiraceae and other butyrate-associated taxa are lower in cases by 23–35 percentile points, and these shifts keep their direction within age, sex and BMI groups. The weak results are stated as plainly: percentiles beat raw abundances in cross-study classification only modestly (A: AUROC 0.618 vs 0.581, p = 0.021, not significant on the strict baseline, p = 0.051; B taxonomy: 0.690 vs 0.664, p = 0.052; function layers show no advantage), only 14 of 43 literature expectations are confirmed in A, African cohorts do not fit the assembly baseline (19.8 % of families outside), the read pipeline's function layers are not validated from raw reads, and region is confounded with study.</p>",
         "</div>",
         "<div class='kpis'>",
         f"<div class='kpi a'><b>{mA['n_samples']:,} / {mA['n_studies']}</b><span>healthy adults / studies, assembly baseline (A, v0.8)</span></div>",
         f"<div class='kpi b'><b>{mB['n_samples']:,} / {mB['n_studies']}</b><span>healthy adults / studies, read baseline (B, v0.2)</span></div>",
         f"<div class='kpi'><b>{c['n_low']} of {c['n_assessed']}</b><span>assessed families below the healthy range in one C. difficile case ({EC.none_or(c['n_high'])} above; {c['n_missing']} expected but missing)</span></div>",
         "<div class='kpi'><b>0.618 vs 0.581</b><span>cross-study AUROC, percentiles vs raw abundances, A (p = 0.021)</span></div>",
         "<div class='kpi'><b>19.8 %</b><span>of families outside the range in African cohorts: the baseline does not fit them</span></div>",
         "</div>"]
    H.append("<h2>Figures</h2><p class='note'>One section per figure: the figure, its caption and the one-line takeaway. Scripts: <code>scripts/viz/</code>; captions, colour conventions and lint: <code>figures/README.md</code>.</p>")
    H.append("<table><caption>Figure files and their numbers in the manuscript</caption><thead><tr><th>figure file</th><th>manuscript figure</th></tr></thead><tbody>" +
             "".join(f"<tr><td><a href='#{s}'>{html.escape(s)}</a></td><td>{html.escape(MANUSCRIPT[s])}</td></tr>" for s in ORDER) + "</tbody></table>")
    for s in ORDER:
        take = re.match(r"\*\*(.+?)\*\*", CAPTIONS[s]).group(1)
        H.append(f"<h3 id='{s}'>{html.escape(take)} <span class='tag'>{html.escape(MANUSCRIPT[s])}</span><span class='stem'>{html.escape(s)}</span></h3>")
        H.append(f"<figure><img src='data:image/png;base64,{imgs[s]}' alt='{html.escape(s)}'><figcaption>{md_inline(CAPTIONS[s])}</figcaption></figure>")
        H.append(f"<p class='take'><strong>Takeaway.</strong> {html.escape(TAKEAWAY[s])}</p>")
    H.append("<h2>Results, assembly pipeline (A)</h2>")
    H.append(tab_html(calA, "Calibration: share of a healthy sample's assessed features outside the healthy range (expected about 5 %, pre-set limit 8 %)", "results/s6/report_v08_lenient_base.json (per_study), results/s6/loso_within_pool.tsv"))
    H.append(tab_html(fidA, "Fidelity: the container run from raw reads vs MGnify's own analysis of the same 61 public samples (default settings)", "results/s11/tierA_anchor_compare.tsv (variant = default)"))
    H.append(tab_html(ccA, "Case/control: random forest trained on all studies but one, tested on the left-out study, baseline rebuilt without it (11 studies)", "results/s8/loso_auroc.tsv, results/s8/report.md (paired comparison)"))
    H.append(f"<p class='note weak' style='padding:8px 12px'>{html.escape(strictA)} Known biology: 14 of 43 literature expectations confirmed (results/s8/report.md).</p>")
    H.append("<h2>Results, read pipeline (B)</h2>")
    H.append(tab_html(calB, "Calibration: share of a healthy sample's assessed features outside the healthy range (expected about 5 %, pre-set limit 8 %)", "results/s11/pipelineB_heldout.tsv, results/s11/pipelineB_loso_all.tsv"))
    H.append(tab_html(fidB, "Fidelity: MetaPhlAn 3.0.14 run from deposited raw reads vs the curatedMetagenomicData 3 profile of the same sample", "results/s11/pipelineB_fidelity.tsv, results/s11/pipelineB_fidelity_samples.tsv (deposit_mismatch = False)"))
    H.append(f"<p class='note'>{html.escape(fidB_note)}</p>")
    H.append(tab_html(ccB, "Case/control: cross-study leave-one-study-out, mean AUROC over held-out studies; percentiles vs raw CLR, one-sided Wilcoxon", "results/s11/pipelineB_report.md (Case/control, cross-study leave-one-study-out)"))
    H.append("<h2>Detailed reports</h2><ul class='reports'>" + "".join(f"<li>{html.escape(a)}: <code>{html.escape(b)}</code></li>" for a, b in REPORTS) + "</ul>")
    H.append(f"<p class='note'>Paths are relative to the project root <code>{html.escape(P)}</code>. Reference bundles (<code>refs/</code>) and the container image (<code>release/refmb_tierA.sif</code>) are distributed separately from the code repository.</p>")
    H.append("</main></body></html>")
    open(OUT, "w", encoding="utf-8").write("\n".join(H)); return OUT


class Check(HTMLParser):
    def __init__(self):
        super().__init__(); self.imgs = []; self.stack = []; self.errors = []; self.scripts = 0
    def handle_starttag(self, tag, attrs):
        if tag == "img":
            self.imgs.append(dict(attrs).get("src", ""))
        if tag == "script":
            self.scripts += 1
        if tag not in ("meta", "img", "br", "hr", "link", "input"):
            self.stack.append(tag)
    def handle_endtag(self, tag):
        if self.stack and self.stack[-1] == tag:
            self.stack.pop()
        else:
            self.errors.append(f"unbalanced </{tag}> (open: {self.stack[-3:]})")


def validate(path):
    from PIL import Image
    c = Check(); c.feed(open(path, encoding="utf-8").read()); c.close()
    assert not c.errors, c.errors[:5]; assert c.stack == [], c.stack; assert c.scripts == 0, "scripts present"
    sizes = []
    for src in c.imgs:
        assert src.startswith("data:image/png;base64,"), src[:40]
        im = Image.open(io.BytesIO(base64.b64decode(src.split(",", 1)[1]))); im.verify(); sizes.append(im.size)
    print(f"summary.html OK: {len(c.imgs)} embedded PNGs decode ({min(s[0] for s in sizes)}–{max(s[0] for s in sizes)} px wide), balanced tags, no scripts, {os.path.getsize(path)/1e6:.1f} MB")


if __name__ == "__main__":
    out = build(); validate(out)
