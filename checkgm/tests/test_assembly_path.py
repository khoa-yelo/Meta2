"""End-to-end test of the assembly-based path: a minimal query directory in the format of docs/query_format.md and a tiny
synthetic assembly bundle (built with checkgm.build_reference, plus a normalization/ directory) go through `checkgm normalize`
and `checkgm score` on the command line. Floors and size classes are the released bundle's; no external data; a few seconds."""
import json
import os
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest

from checkgm.build_reference import build

GRID = [1, 2.5, 5, 10, 25, 50, 75, 90, 95, 97.5, 99]
# species taxid, genus taxid, family taxid, plus dummy order/class/phylum ids shared by all
TAXA = [(820, 816, 815), (853, 216851, 216572), (1680, 1678, 31953), (239935, 239934, 1647988), (562, 561, 543), (1301, 1300, 1300001)]
EXTRA_FAMILIES = [186803, 186806]          # in the basis and the pool, never detected in the query -> expected_but_missing
KOS = [f"K{i:05d}" for i in range(2860, 2890)]   # 30 KOs; the first 25 form the marker panel
PFAMS = [f"PF{i:05d}" for i in range(1, 21)]
MODULES = [f"M{i:05d}" for i in range(1, 11)]
RULES = {"taxonomy_primary_rank": "family", "clr_ranks": ["family", "genus"], "clr_delta_fraction_of_min": 0.5, "contig_majority": 0.5,
         "floors": {"min_total_length_bp": 5_000_000, "min_n50_bp": 1000, "min_predicted_cds": 5000, "min_markers_detected": 20,
                    "min_genome_equivalents": 10, "min_mapped_fraction": 0.5},
         "quality_bands": {"low": [5_000_000, 90_000_000], "medium": [90_000_000, 180_000_000], "high": [180_000_000, None]}}


def write_backbone(bb: str):
    rows = [dict(taxid=1, parent=1, rank="no rank", name="root", phylum_taxid=-1, class_taxid=-1, order_taxid=-1, family_taxid=-1, genus_taxid=-1, species_taxid=-1, domain_taxid=-1)]
    for sp, g, f in TAXA:
        base = dict(phylum_taxid=1224, class_taxid=1236, order_taxid=91347, domain_taxid=2)
        rows.append(dict(taxid=sp, parent=g, rank="species", name=f"species {sp}", family_taxid=f, genus_taxid=g, species_taxid=sp, **base))
        rows.append(dict(taxid=g, parent=f, rank="genus", name=f"genus {g}", family_taxid=f, genus_taxid=g, species_taxid=-1, **base))
        rows.append(dict(taxid=f, parent=91347, rank="family", name=f"family {f}", family_taxid=f, genus_taxid=-1, species_taxid=-1, **base))
    for f in EXTRA_FAMILIES:
        rows.append(dict(taxid=f, parent=91347, rank="family", name=f"family {f}", phylum_taxid=1224, class_taxid=1236, order_taxid=91347, family_taxid=f, genus_taxid=-1, species_taxid=-1, domain_taxid=2))
    os.makedirs(bb, exist_ok=True)
    pd.DataFrame(rows).to_parquet(os.path.join(bb, "ncbi_taxonomy.parquet"), index=False)
    pd.DataFrame({"old": [999999], "new": [562]}).to_parquet(os.path.join(bb, "ncbi_merged.parquet"), index=False)


def make_bundle(root: str, seed: int = 0) -> str:
    """A tiny assembly-based bundle: features built with checkgm.build_reference from a synthetic pool, plus normalization/."""
    rng = np.random.default_rng(seed)
    out = os.path.join(root, "gut-assembly-test-v0.8-lenient")
    n = 36
    samples = pd.DataFrame({"sample_id": [f"P{i:03d}" for i in range(n)], "study": [f"PRJ{i % 3}" for i in range(n)],
                            "band": [("low", "medium", "high")[i % 3] for i in range(n)]})
    fam_basis = sorted(str(f) for _, _, f in TAXA) + sorted(str(f) for f in EXTRA_FAMILIES)
    gen_basis = sorted(str(g) for _, g, _ in TAXA)
    layers = {k: [] for k in ("taxonomy_family", "taxonomy_genus", "ko_eggnog", "ko_kofam", "pfam", "module")}
    for sid in samples["sample_id"]:
        for rank, basis in (("family", fam_basis), ("genus", gen_basis)):
            p = rng.dirichlet(np.ones(len(basis))); p[rng.random(len(basis)) < 0.15] = 0; p = p / p.sum()
            z = int((p == 0).sum()); delta = 0.5 * p[p > 0].min(); rep = np.where(p > 0, p * (1 - z * delta), delta)
            clr = np.log(rep) - np.log(rep).mean()
            layers[f"taxonomy_{rank}"].append(pd.DataFrame({"sample_id": sid, "feature_id": basis, "value": clr, "detect": p}))
        for layer, feats in (("ko_eggnog", KOS), ("ko_kofam", KOS), ("pfam", PFAMS)):
            v = np.exp(rng.normal(0, 0.5, len(feats)))
            layers[layer].append(pd.DataFrame({"sample_id": sid, "feature_id": feats, "value": v, "detect": v}))
        c = rng.uniform(0.3, 1.0, len(MODULES))
        layers["module"].append(pd.DataFrame({"sample_id": sid, "feature_id": MODULES, "value": c, "detect": c}))
    layers = {k: pd.concat(v, ignore_index=True) for k, v in layers.items()}
    spec = {"percentiles": GRID, "min_detected_for_percentiles": 5, "bootstrap_n": 0, "band_layers": ["ko_eggnog", "ko_kofam", "pfam", "module"],
            "gene_layers": ["ko_eggnog", "ko_kofam", "pfam"], "detection_floor": 0.05,
            "landscape": {"layer": "taxonomy_family", "n_components": 3, "min_prevalence": 0.1, "k_neighbors": 5}}
    build(out, samples, layers, spec, {"profile_type": "assembly", "reference_key": {"body_site": "Gut", "stratum": "adult-global"},
                                        "source_pipeline": {"name": "mgnify", "version": "5.0 (test)"}})
    nd = os.path.join(out, "normalization"); os.makedirs(nd, exist_ok=True)
    json.dump(RULES, open(os.path.join(nd, "rules.json"), "w"))
    pd.DataFrame({"ko": KOS[:25], "name": ["marker"] * 25}).to_csv(os.path.join(nd, "marker_panel.tsv"), sep="\t", index=False)
    open(os.path.join(nd, "clr_basis_family.txt"), "w").write("\n".join(fam_basis) + "\n")
    open(os.path.join(nd, "clr_basis_genus.txt"), "w").write("\n".join(gen_basis) + "\n")
    write_backbone(os.path.join(nd, "backbone"))
    return out


def write_query(qdir: str, seed: int = 1):
    """Minimal query directory of docs/query_format.md: 120 contigs (96 Mb in all, size class medium), 5,400 CDS."""
    rng = np.random.default_rng(seed)
    os.makedirs(qdir, exist_ok=True)
    n_contigs = 120
    contigs = pd.DataFrame({"contig_id": [f"c{i}" for i in range(n_contigs)], "length": 800_000, "coverage": rng.uniform(12, 40, n_contigs).round(2)})
    contigs.to_csv(os.path.join(qdir, "contigs.tsv"), sep="\t", index=False)
    taxon_of_contig = [TAXA[i % len(TAXA)][0] if i % 10 else -1 for i in range(n_contigs)]   # every 10th contig has no taxonomy
    cds = []; tax = []
    k = 0
    for i in range(n_contigs):
        for j in range(45):
            cid = f"c{i}_{j}"
            ko = ""
            if j % 3:                                   # two CDS in three carry a KO; the counter walks through all 30 KOs
                ko = KOS[k % len(KOS)]; k += 1
            cds.append((cid, f"c{i}", ko, PFAMS[j % len(PFAMS)] if j % 2 else "", ko if j % 4 else ""))
            if taxon_of_contig[i] != -1 and j % 5:
                tax.append((cid, taxon_of_contig[i]))
    pd.DataFrame(cds, columns=["cds_id", "contig_id", "ko", "pfam", "ko_kofam"]).to_csv(os.path.join(qdir, "cds.tsv"), sep="\t", index=False)
    pd.DataFrame(tax, columns=["cds_id", "taxid"]).to_csv(os.path.join(qdir, "cds_taxonomy.tsv"), sep="\t", index=False)
    pd.DataFrame({"module_id": MODULES, "completeness": np.linspace(0.2, 1.0, len(MODULES))}).to_csv(os.path.join(qdir, "modules.tsv"), sep="\t", index=False)
    json.dump({"sample_id": "Q1", "body_site": "Gut", "pipeline": {"name": "test-assembly", "version": "0"}}, open(os.path.join(qdir, "sample.json"), "w"))


def run(*args):
    return subprocess.run([sys.executable, "-m", "checkgm.cli", *args], capture_output=True, text=True)


def test_normalize_then_score_through_the_cli(tmp_path):
    bundle = make_bundle(str(tmp_path)); q = tmp_path / "Q1"; write_query(str(q))
    r = run("normalize", "--bundle", bundle, "--query", str(q), "--out", str(tmp_path / "norm"))
    assert r.returncode == 0, r.stderr
    assert r.stdout.split()[:3] == ["Q1", "retained", "medium"]
    qc = pd.read_csv(tmp_path / "norm" / "qc.tsv", sep="\t")
    assert qc.loc[0, "quality_band"] == "medium" and qc.loc[0, "total_length_all_contigs"] == pytest.approx(96_000_000)
    assert qc.loc[0, "n_cds"] == 5400 and qc.loc[0, "marker_kos_detected_kept"] == 25 and qc.loc[0, "genome_equivalents_cov"] > 10
    assert 0.5 <= qc.loc[0, "mapped_fraction_family"] < 1.0        # the contigs without taxonomy sit in the unmapped bucket
    t = pd.read_parquet(tmp_path / "norm" / "normalized" / "Q1.parquet")
    for rank in ("family", "genus"):
        g = t[t["layer"] == f"taxonomy_{rank}"]
        assert abs(g["clr"].sum()) < 1e-9 and g["proportion_mapped"].sum() == pytest.approx(1.0)
    assert set(t["layer"]) == {"taxonomy_family", "taxonomy_genus", "ko_eggnog", "ko_kofam", "pfam", "module"}
    assert json.load(open(tmp_path / "norm" / "sample_meta.json"))["Q1"]["pipeline"]["name"] == "test-assembly"

    r = run("score", "--bundle", bundle, "--normalized", str(tmp_path / "norm"), "--out", str(tmp_path / "rep"))
    assert r.returncode == 0, r.stderr
    sc = pd.read_parquet(tmp_path / "rep" / "scores.parquet")
    assert set(sc["layer"]) == {"taxonomy_family", "taxonomy_genus", "ko_eggnog", "ko_kofam", "pfam", "module"}
    assert sc["band"].eq("medium").all()
    assert sc.loc[sc["call"] != "expected_but_missing", "percentile"].between(0, 100).all()
    assert set(sc.loc[sc["call"] == "expected_but_missing", "feature_id"]) >= {str(f) for f in EXTRA_FAMILIES}
    S = pd.read_csv(tmp_path / "rep" / "summaries.tsv", sep="\t")
    assert sorted(S["layer"]) == sorted({"taxonomy_family", "taxonomy_genus", "ko_eggnog", "ko_kofam", "pfam", "module"})
    assert (S["layer_status"] == "scored").all() and S["core_distance_pct"].between(0, 100).all()
    assert pd.read_csv(tmp_path / "rep" / "rejections.tsv", sep="\t").empty
    rep = (tmp_path / "rep" / "report.md").read_text()
    assert "# checkgm report" in rep and "genome_equivalents_cov" in rep       # the assembly report keeps the coverage column


def test_query_below_a_floor_is_rejected_not_an_error(tmp_path):
    bundle = make_bundle(str(tmp_path)); q = tmp_path / "small"; write_query(str(q))
    c = pd.read_csv(q / "contigs.tsv", sep="\t"); c["length"] = 20_000; c.to_csv(q / "contigs.tsv", sep="\t", index=False)   # 2.4 Mb assembled
    r = run("normalize", "--bundle", bundle, "--query", str(q), "--out", str(tmp_path / "norm"))
    assert r.returncode == 0 and "BELOW_MIN_ASSEMBLY_QUALITY unbanded" in r.stdout
    r = run("score", "--bundle", bundle, "--normalized", str(tmp_path / "norm"), "--out", str(tmp_path / "rep"))
    assert r.returncode == 0, r.stderr
    rej = pd.read_csv(tmp_path / "rep" / "rejections.tsv", sep="\t")
    assert list(rej["reason_code"]) == ["NOT_NORMALIZED"] and "BELOW_MIN_ASSEMBLY_QUALITY" in rej.loc[0, "detail"]


def test_run_normalize_matches_normalize_then_score(tmp_path):
    """`run-normalize` is the one-step form of `normalize` followed by `score`, so it must produce the same report
    as the two steps it replaces; it was added without a test, and an untested shortcut is how the two paths drift."""
    bundle = make_bundle(str(tmp_path / "b")); q = tmp_path / "q"; write_query(str(q))

    two = tmp_path / "two"
    assert run("normalize", "--bundle", bundle, "--query", str(q), "--out", str(two / "norm")).returncode == 0
    assert run("score", "--bundle", bundle, "--normalized", str(two / "norm"), "--out", str(two / "rep")).returncode == 0

    one = tmp_path / "one"
    r = run("run-normalize", "--bundle", bundle, "--query", str(q), "--out", str(one))
    assert r.returncode == 0, r.stderr
    for f in ("scores.parquet", "summaries.tsv", "rejections.tsv", "report.md", "qc.tsv"):
        assert (one / f).exists(), f

    a = pd.read_parquet(two / "rep" / "scores.parquet").sort_values(["analysis_id", "layer", "feature_id"]).reset_index(drop=True)
    b = pd.read_parquet(one / "scores.parquet").sort_values(["analysis_id", "layer", "feature_id"]).reset_index(drop=True)
    pd.testing.assert_frame_equal(a, b)


def test_too_few_taxa_is_rejected_not_scored_as_nan(tmp_path):
    """A query whose taxonomy spans only one basis family passes every documented floor, but the multiplicative
    replacement of zeros is then undefined (z * delta exceeds the observed mass) and every CLR is NaN. That used to be
    reported as a successful scoring: NaN percentiles, a landscape position derived from the NaN vector, and an empty
    rejections.tsv. It must be a rejection with a named reason, as the read path already does for an unplaceable
    profile, and no bare numpy warning may reach stderr."""
    bundle = make_bundle(str(tmp_path / "b")); q = tmp_path / "q"; write_query(str(q))
    # collapse the taxonomy onto a single taxon, leaving every other floor satisfied
    t = pd.read_csv(q / "cds_taxonomy.tsv", sep="\t")
    t["taxid"] = TAXA[0][0]
    t.to_csv(q / "cds_taxonomy.tsv", sep="\t", index=False)

    r = run("normalize", "--bundle", bundle, "--query", str(q), "--out", str(tmp_path / "norm"))
    assert r.returncode == 0, r.stderr
    assert "RuntimeWarning" not in r.stderr and "invalid value" not in r.stderr, r.stderr

    qc = pd.read_csv(tmp_path / "norm" / "qc.tsv", sep="\t")
    assert qc["status"].eq("TOO_FEW_TAXA_FOR_CLR").all(), qc[["sample_id", "status"]].to_dict("records")

    r = run("score", "--bundle", bundle, "--normalized", str(tmp_path / "norm"), "--out", str(tmp_path / "rep"))
    assert r.returncode == 0, r.stderr
    rej = pd.read_csv(tmp_path / "rep" / "rejections.tsv", sep="\t")
    assert len(rej) == 1 and "centred log-ratio" in rej.iloc[0]["detail"], rej.to_dict("records")
    # and nothing is passed off as a score
    sc = pd.read_parquet(tmp_path / "rep" / "scores.parquet")
    assert sc.empty or not sc["layer"].str.startswith("taxonomy").any(), sc.head().to_dict("records")
