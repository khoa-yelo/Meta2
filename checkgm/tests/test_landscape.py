"""The healthy map: samples of the reference pool, re-projected through checkgm.score, must land exactly on the coordinates the
bundle stores (landscape/reference_coords.parquet), and their core-distance percentiles must be spread evenly.

The pool matrix the builder fits the PCA on holds a CLR only for detected features and fills the rest with the pool minimum;
a normalized query carries a multiplicative-replacement CLR for every basis feature. The scorer has to mask undetected
features before filling, or every sample is displaced along PC1 (release 0.8.0 reported healthy samples near the 90th
core-distance percentile for this reason)."""
import json
import os

import numpy as np
import pandas as pd
import pytest

from checkgm.build_reference import build
from checkgm.score import Bundle, score

GRID = [1, 2.5, 5, 10, 25, 50, 75, 90, 95, 97.5, 99]
FAMILIES = [str(f) for f in (815, 216572, 31953, 1647988, 543, 186803, 186806, 541000, 1570339, 128827)]


def clr_table(prop: pd.Series, delta_frac: float = 0.5) -> pd.DataFrame:
    """The CLR convention of checkgm.normalize / normalize_reads: multiplicative replacement of zeros over the whole basis."""
    b = prop.reindex(FAMILIES).fillna(0.0)
    z = int((b == 0).sum()); delta = delta_frac * (b[b > 0].min() if (b > 0).any() else 1.0)
    rep = np.where(b > 0, b * (1 - z * delta), delta)
    clr = np.log(rep) - np.log(rep).mean()
    return pd.DataFrame({"feature_id": b.index, "proportion_mapped": b.to_numpy(), "clr": clr, "rank": "family"})


def make_pool(n: int = 60, seed: int = 3):
    rng = np.random.default_rng(seed)
    samples = pd.DataFrame({"sample_id": [f"S{i:03d}" for i in range(n)], "study": [f"P{i % 3}" for i in range(n)],
                            "band": [("low", "medium", "high")[i % 3] for i in range(n)]})
    norm = {}; rows = []
    for sid, band in zip(samples["sample_id"], samples["band"]):
        p = rng.dirichlet(np.full(len(FAMILIES), 0.6))
        p[rng.random(len(FAMILIES)) < 0.3] = 0.0            # undetected families, as in real assemblies
        if (p > 0).sum() < 3:
            p[:3] = 0.2
        p = p / p.sum()
        t = clr_table(pd.Series(p, index=FAMILIES))
        norm[sid] = {"taxonomy_family": t, "qc": {"status": "retained", "quality_band": band, "mapped_fraction_family": 1.0,
                                                 "genome_equivalents_cov": np.nan, "rejections": []}}
        rows.append(pd.DataFrame({"sample_id": sid, "feature_id": t["feature_id"], "value": t["clr"], "detect": t["proportion_mapped"]}))
    return samples, {"taxonomy_family": pd.concat(rows, ignore_index=True)}, norm


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    samples, layers, norm = make_pool()
    out = os.path.join(str(tmp_path_factory.mktemp("bundles")), "gut-assembly-test-v0.8-lenient")
    spec = {"percentiles": GRID, "min_detected_for_percentiles": 5, "bootstrap_n": 0, "band_layers": [], "gene_layers": [],
            "detection_floor": 0.0, "landscape": {"layer": "taxonomy_family", "n_components": 3, "min_prevalence": 0.1, "k_neighbors": 5}}
    build(out, samples, layers, spec, {"profile_type": "assembly", "reference_key": {"body_site": "Gut", "stratum": "adult-global"},
                                        "source_pipeline": {"name": "test", "version": "0"}})
    return out, norm


def test_pool_samples_reproject_onto_stored_coordinates(built):
    out, norm = built
    ref = pd.read_parquet(os.path.join(out, "landscape", "reference_coords.parquet"))
    meta = json.load(open(os.path.join(out, "landscape", "core_distance.json")))
    assert meta["feature_fill_value"] == pytest.approx(min(t["taxonomy_family"].loc[t["taxonomy_family"]["proportion_mapped"] > 0, "clr"].min() for t in norm.values()))
    S = score(norm, Bundle(out))["summaries"].drop_duplicates("analysis_id").sort_values("analysis_id").reset_index(drop=True)
    assert len(S) == len(ref)                       # build() sorts the pool by sample id; the pivot in score() does the same
    for c in ("PC1", "PC2", "PC3", "core_distance"):
        assert np.abs(S[c].to_numpy() - ref[c].to_numpy()).max() < 1e-6, c
    assert S["landscape_centring"].eq("reference").all()


def test_pool_core_distance_percentiles_are_uniform(built):
    out, norm = built
    S = score(norm, Bundle(out))["summaries"].drop_duplicates("analysis_id")
    pct = S["core_distance_pct"].to_numpy()
    assert 40 <= np.median(pct) <= 60
    assert pct.min() < 10 and pct.max() > 90
    assert S["neighbor_studies"].map(lambda s: sum(json.loads(s).values())).eq(5).all()   # k nearest reference samples


def test_undetected_features_are_masked_not_taken_at_face_value(built):
    """Replacing the undetected features' CLR by any other number must not move a sample: only detected values count."""
    out, norm = built
    sid = next(iter(norm)); t = norm[sid]["taxonomy_family"].copy()
    t.loc[t["proportion_mapped"] == 0, "clr"] = -40.0
    a = score({sid: norm[sid]}, Bundle(out))["summaries"].iloc[0]
    b = score({sid: {**norm[sid], "taxonomy_family": t}}, Bundle(out))["summaries"].iloc[0]
    assert a["PC1"] == pytest.approx(b["PC1"]) and a["core_distance"] == pytest.approx(b["core_distance"])
