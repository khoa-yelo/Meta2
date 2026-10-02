# Rebuilding a baseline

The released bundles were built from public data with the project's build scripts. This page records what a rebuild
involves so the baselines can be reproduced or extended; the per-stage scripts live in the project workspace and are not
part of the installable package (they depend on cluster storage holding the staged analyses). The bundle builder itself
is in the package: `refmb.build_reference.build`.

## Assembly-based baseline (MGnify v5)

Inputs: the curated sample metadata (health status, age, body site, study) and the MGnify v5 assembly analyses of the
candidate samples. The metadata drives labels, stratum membership (adult, stool, healthy), study splits and pool
membership; it does not change how an analysis is normalized or how the container measures a sample.

Stages, in order, with typical wall time when no new analyses need staging:

| stage | reads metadata? | time | output |
|---|---|---|---|
| inventory and splits | yes | 5 min | which analyses are candidates, held-out, case/control; printed role changes |
| staging and per-analysis parsing | only new analysis ids | hours per ~1,000 new analyses | normalized long tables per analysis |
| matrices | no | 15 min | per-layer long tables (`sample_id`, `feature_id`, `value`, `detect`) |
| bundles (lenient, strict) | yes | 10 min each | `refs/gut-assembly-adult-global-v<V>-<tier>` via `refmb.build_reference.build`, study bootstrap 500 |
| held-out calibration check | yes | 10 min | fraction outside the band on held-out healthy cohorts (gate 8 % per study) |
| scoring of all gut analyses | labels only | 1–2 h | out-of-sample scores for every candidate |
| leave-one-study-out fold bundles and evaluation | yes | 1–2 h | case/control AUROC with the held-out study's reference rebuilt without it |
| packaging and round trip | no | 15 min | tarballs, sha256, re-score check |
| strata and disease atlas | yes | 1–2 h | healthy baseline by age, sex, BMI, region; disease deviations |
| figures | no | 5 min | |

A version bump writes new bundles beside the old ones; nothing is overwritten. Every stage writes a stamp keyed by the
sha256 of its inputs, so a second run re-does only what changed, and every stage's report keeps `in = out + rejected`.

Before trusting a rebuilt baseline:

1. Read the inventory report and the split role changes.
2. Held-out cohorts within the 8 % gate on the layers that passed before.
3. Leave-one-study-out AUROC and the paired comparison; strata and disease atlas.
4. Round trip identical in the release manifest.
5. Record the decision with the input file hashes.

## Read-based baseline (curatedMetagenomicData 3)

Inputs: the curatedMetagenomicData 3 MetaPhlAn 3 species profiles and metadata (Zenodo release of cMD3), and its HUMAnN 3
tables for the function layers. Stages: parse the profiles → place species on the NCBI backbone and form the CLR tables
(`refmb.normalize_reads` rules) → build the bundle (`refmb.build_reference.build`) → check held-out healthy studies
(CosteaPI_2017, XieH_2016, YeZ_2018 are never in the pool) → package. A refresh of cMD does not touch the assembly baseline.

## `refmb.build_reference.build`

```python
from refmb.build_reference import build
manifest = build(out="refs/my-bundle-v0.1-lenient",
                 samples=samples,            # DataFrame: sample_id, study, band
                 layers=layers,              # {layer: DataFrame sample_id, feature_id, value, detect}
                 spec={"percentiles": [1, 2.5, 5, 10, 25, 50, 75, 90, 95, 97.5, 99], "min_detected_for_percentiles": 20,
                       "bootstrap_n": 500, "band_layers": [...], "gene_layers": [...], "detection_floor": 0.05,
                       "landscape": {"layer": "taxonomy_family", "n_components": 10, "min_prevalence": 0.05, "k_neighbors": 10}},
                 manifest_extra={...}, exclusions=exclusions, seed=0)
```

Samples are weighted so that each study carries equal total weight; percentiles are weighted quantiles among detected
samples; bands get their own prevalence and (for `band_layers`) their own percentiles; a study-level bootstrap gives
intervals on the 2.5th and 97.5th percentiles. The `normalization/` directory (rules, basis lists, backbone, marker panel
or species map) must be added beside the output so that queries can be normalized the way the pool was.
