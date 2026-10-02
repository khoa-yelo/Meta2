"""Minimal end-to-end test: build a tiny synthetic read-based bundle, parse a MetaPhlAn 3 profile, normalize and score it,
through the library and through the command line. No external data; runs in seconds."""
import json
import os
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest

import refmb
from refmb.normalize_reads import ReadsBundleSpec, normalize_profile, read_metaphlan
from refmb.paths import resolve_bundle
from refmb.score import Bundle, score

GRID = [1, 2.5, 5, 10, 25, 50, 75, 90, 95, 97.5, 99]
# species (taxid, genus taxid, family taxid, MetaPhlAn name)
SPECIES = [
    (820, 816, 815, "s__Bacteroides_uniformis"),
    (821, 816, 815, "s__Bacteroides_vulgatus"),
    (853, 216851, 216572, "s__Faecalibacterium_prausnitzii"),
    (1680, 1678, 31953, "s__Bifidobacterium_adolescentis"),
    (239935, 239934, 1647988, "s__Akkermansia_muciniphila"),
    (562, 561, 543, "s__Escherichia_coli"),
]


def make_bundle(root: str, seed: int = 0) -> str:
    rng = np.random.default_rng(seed)
    b = os.path.join(root, "gut-reads-test-v0.0-lenient")
    for sub in ("features", "normalization/backbone"):
        os.makedirs(os.path.join(b, sub), exist_ok=True)
    # NCBI-like backbone: species, genus and family rows
    rows = []
    for sp, g, f, name in SPECIES:
        rows.append(dict(taxid=sp, parent=g, rank="species", name=name[3:].replace("_", " "), family_taxid=f, genus_taxid=g, species_taxid=sp))
        rows.append(dict(taxid=g, parent=f, rank="genus", name=f"genus{g}", family_taxid=f, genus_taxid=g, species_taxid=-1))
        rows.append(dict(taxid=f, parent=1, rank="family", name=f"family{f}", family_taxid=f, genus_taxid=-1, species_taxid=-1))
    tax = pd.DataFrame(rows).drop_duplicates("taxid")
    tax.to_parquet(os.path.join(b, "normalization", "backbone", "ncbi_taxonomy.parquet"), index=False)
    pd.DataFrame({"old": [999999], "new": [562]}).to_parquet(os.path.join(b, "normalization", "backbone", "ncbi_merged.parquet"), index=False)
    pd.DataFrame([dict(species_name=n, species=s, genus=g, family=f) for s, g, f, n in SPECIES]).to_csv(os.path.join(b, "normalization", "species_taxid_map.tsv"), sep="\t", index=False)
    basis = {"species": sorted({str(s) for s, _, _, _ in SPECIES}), "genus": sorted({str(g) for _, g, _, _ in SPECIES}), "family": sorted({str(f) for _, _, f, _ in SPECIES})}
    for r, feats in basis.items():
        open(os.path.join(b, "normalization", f"clr_basis_{r}.txt"), "w").write("\n".join(feats) + "\n")
    json.dump({"depth_bands_reads": {"low": [0, 1e7], "medium": [1e7, 3e7], "high": [3e7, None]}, "clr_delta_fraction_of_min": 0.5,
               "min_mapped_fraction_family": 0.9}, open(os.path.join(b, "normalization", "rules.json"), "w"))
    layers = {}
    for r, feats in basis.items():
        rec = []
        for f in feats:
            mu = rng.normal(0, 1); sd = 0.8
            q = mu + sd * np.array([-2.33, -1.96, -1.64, -1.28, -0.67, 0, 0.67, 1.28, 1.64, 1.96, 2.33])
            d = dict(feature_id=f, n_detected=90, n_studies=3, prevalence_weighted=0.9, percentiles_available=True, layer=f"taxonomy_{r}", reliability_tier="NA")
            for band in ("low", "medium", "high"):
                d[f"prevalence_band_{band}"] = 0.9; d[f"n_band_{band}"] = 30
            for p, v in zip(GRID, q):
                d[f"p{str(p).replace('.', '_')}"] = v
            rec.append(d)
        pd.DataFrame(rec).to_parquet(os.path.join(b, "features", f"taxonomy_{r}.parquet"), index=False)
        layers[f"taxonomy_{r}"] = {"n_features": len(feats)}
    json.dump({"bundle_id": os.path.basename(b), "n_samples": 100, "n_studies": 3, "percentile_grid": GRID, "layers": layers,
               "profile_type": "reads", "source_pipeline": {"name": "MetaPhlAn", "version": "3 (test)"},
               "reference_key": {"body_site": "Gut", "stratum": "adult-global"},
               "scoring": {"detection_floor_cpg": 0.0, "percentiles_within_band": False, "percentiles_within_band_layers": []}},
              open(os.path.join(b, "manifest.json"), "w"))
    return b


def write_profile(path: str, reads: int = 50_000_000):
    lines = ["#mpa_v30_CHOCOPhlAn_201901", "#/usr/bin/metaphlan sample.fq --input_type fastq", f"#{reads} reads processed",
             "#SampleID\tMetaphlan_Analysis", "#clade_name\tNCBI_tax_id\trelative_abundance\tadditional_species"]
    ab = [40.0, 25.0, 20.0, 10.0, 4.0, 1.0]
    for (sp, g, f, name), a in zip(SPECIES, ab):
        lines.append(f"k__Bacteria|p__X|c__X|o__X|f__family{f}|g__genus{g}|{name}\t2|1|1|1|{f}|{g}|{sp}\t{a:.5f}\t")
    lines.append("k__Bacteria|p__X|c__X|o__X|f__family815|g__genus816\t2|1|1|1|815|816\t65.0\t")   # genus row: must be ignored
    open(path, "w").write("\n".join(lines) + "\n")


@pytest.fixture(scope="module")
def bundle(tmp_path_factory):
    return make_bundle(str(tmp_path_factory.mktemp("bundles")))


def test_version():
    assert refmb.__version__ == "0.8.0"


def test_parse_normalize_score(bundle, tmp_path):
    prof = tmp_path / "S1.txt"; write_profile(str(prof))
    parsed = read_metaphlan(str(prof))
    assert list(parsed) == ["S1"] and parsed["S1"]["reads"] == 50_000_000 and len(parsed["S1"]["table"]) == 6
    spec = ReadsBundleSpec(bundle)
    n = normalize_profile("S1", parsed["S1"]["table"], parsed["S1"]["reads"], spec)
    assert n["qc"]["status"] == "retained" and n["qc"]["quality_band"] == "high"
    assert n["qc"]["mapped_fraction_family"] == pytest.approx(1.0)
    fam = n["taxonomy_family"]
    assert abs(fam["clr"].sum()) < 1e-9 and fam["proportion_mapped"].sum() == pytest.approx(1.0)
    r = score({"S1": n}, Bundle(bundle))
    assert r["rejections"].empty
    assert set(r["scores"]["layer"]) == {"taxonomy_family", "taxonomy_genus", "taxonomy_species"}
    s = r["scores"]; assert s["percentile"].between(0, 100).all()
    assert set(s["call"]) <= {"low", "high", "within"}
    S = r["summaries"]; assert len(S) == 3 and (S["layer_status"] == "scored").all()
    assert (S["n_assessed"] == S["n_detected"]).all() and S["excess_outside_q"].notna().all()


def test_percentile_of_median_is_50(bundle, tmp_path):
    """A value equal to the reference median must score at the 50th percentile."""
    from refmb.score import interp_percentile
    grid = np.array(GRID, dtype=float); q = np.array([-2, -1.5, -1, -0.5, -0.2, 0, 0.2, 0.5, 1, 1.5, 2.0])
    assert interp_percentile(np.array([0.0]), np.array([q]), grid)[0] == pytest.approx(50.0)
    assert interp_percentile(np.array([-9.0]), np.array([q]), grid)[0] == pytest.approx(0.5)


def test_unresolvable_profile_is_rejected(bundle):
    """Species the bundle cannot place on its backbone leave the family mapped fraction below the floor: not scored, with a reason."""
    spec = ReadsBundleSpec(bundle)
    t = pd.DataFrame({"species_name": ["s__Unknownus_novus", "s__Bacteroides_uniformis"], "taxid": ["424242", "820"], "rel_abundance": [0.95, 0.05]})
    n = normalize_profile("odd", t, 1_000_000, spec)
    assert n["qc"]["status"] == "LOW_MAPPED_FRACTION" and n["qc"]["quality_band"] == "unbanded"
    assert n["qc"]["mapped_fraction_family"] == pytest.approx(0.05)
    r = score({"odd": n}, Bundle(bundle))
    assert list(r["rejections"]["reason_code"]) == ["NOT_NORMALIZED"] and r["summaries"].empty and r["scores"].empty


def test_cli_run_metaphlan_with_env_lookup(bundle, tmp_path, monkeypatch):
    prof = tmp_path / "S2.txt"; write_profile(str(prof), reads=5_000_000)
    monkeypatch.setenv("REFMB_BUNDLES", os.path.dirname(bundle))
    assert resolve_bundle(os.path.basename(bundle)) == bundle
    out = tmp_path / "report"
    env = dict(os.environ, REFMB_BUNDLES=os.path.dirname(bundle))
    res = subprocess.run([sys.executable, "-m", "refmb.cli", "run-metaphlan", "--bundle", os.path.basename(bundle), "--input", str(prof), "--out", str(out)],
                         capture_output=True, text=True, env=env)
    assert res.returncode == 0, res.stderr
    for f in ("scores.parquet", "summaries.tsv", "rejections.tsv", "report.md", "qc.tsv", "normalized/S2.parquet"):
        assert (out / f).exists(), f
    S = pd.read_csv(out / "summaries.tsv", sep="\t")
    assert S["quality_band"].eq("low").all() and S["bundle_id"].eq(os.path.basename(bundle)).all()
    assert "# refmb report" in (out / "report.md").read_text()
    bad = subprocess.run([sys.executable, "-m", "refmb.cli", "score", "--bundle", "no-such-bundle", "--normalized", str(out), "--out", str(tmp_path / "x")], capture_output=True, text=True, env=env)
    assert bad.returncode != 0 and "not found" in (bad.stderr + bad.stdout)
