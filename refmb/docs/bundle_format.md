# Bundle format

A bundle is a directory named by its id, for example `gut-reads-adult-global-v0.2-lenient`. It holds summary statistics
of the reference pool and the rules needed to normalize a new sample the same way. It contains no sample-level values.

```
<bundle_id>/
  manifest.json
  README.md                      human-readable summary (pool, layers, studies, exclusions)
  CHECKSUMS.sha256
  features/<layer>.parquet       one file per layer
  landscape/pca_loadings.parquet, reference_coords.parquet, core_distance.json
  normalization/                 rules and lookup tables for query normalization
  provenance/studies.tsv, pool_exclusions.tsv, build_log.json, inclusion_criteria.yaml (assembly bundles)
```

## `manifest.json`

| field | meaning |
|---|---|
| `bundle_id` | directory name; `gut-<assembly|reads>-adult-global-v<version>-<lenient|strict>[-loso-<study>]` |
| `created`, `builder_version`, `wall_seconds` | build provenance |
| `profile_type` | `assembly` or `reads` |
| `source_pipeline` | `{name, version, workflow}` of the measurement the pool was produced with (MGnify 5.0 assembly; MetaPhlAn + HUMAnN 3 of curatedMetagenomicData 3) |
| `reference_key` | `{body_site: Gut, stratum: adult-global}`; the scorer refuses samples declaring another body site |
| `inclusion_tier` | `lenient` or `strict` |
| `n_samples`, `n_studies`, `band_counts` | pool size and quality-band composition |
| `study_weighting` | `equal_per_study`: every study contributes the same total weight to prevalence and percentiles |
| `percentile_grid` | `[1, 2.5, 5, 10, 25, 50, 75, 90, 95, 97.5, 99]` |
| `bootstrap` | `{n: 500, unit: study}`: study-level bootstrap behind the confidence intervals on the tail percentiles |
| `assembly_quality_bands` | assembly bundles: total assembled length (bp) bounds of the `low`, `medium`, `high` bands |
| `depth_bands_reads` | read-based bundles: reads-processed bounds of the bands |
| `scoring` | `detection_floor_cpg` (0.05 on assembly gene layers, 0 otherwise), `percentiles_within_band_layers` (layers whose percentiles are taken within the sample's band) |
| `layers` | per layer: `n_samples`, `n_studies`, `n_features`, `n_with_percentiles`, `non_monotone_percentile_vectors`, `prevalence_bands` |
| `landscape` | `layer` (taxonomy_family), `k_neighbors`, `variance_explained` |
| `normalization` | how values were formed (taxonomy weighting, CLR basis, gene value definition) and a pointer to `normalization/` |
| `annotation_versions`, `config_sha256` | database versions and hashes of the configuration the pool was processed with |
| `checksums` | sha256 of every `features/*.parquet` |
| `loso_excluded_studies`, `fallback_bundle`, `provisional` | leave-one-study-out variants list the excluded study; the released bundles have none |

## `features/<layer>.parquet`

Layers: `taxonomy_family`, `taxonomy_genus`, `ko_eggnog`, `ko_kofam`, `pfam`, `module` (assembly);
`taxonomy_family`, `taxonomy_genus`, `taxonomy_species`, `ko_humann`, `pathway_humann` (reads). One row per feature.

| column | meaning |
|---|---|
| `feature_id` | NCBI taxid as a string for taxonomy layers; KO id (`K00001`), Pfam accession, KEGG module id, or MetaCyc pathway id otherwise |
| `n_detected`, `n_studies` | reference samples and studies in which the feature was detected |
| `prevalence_weighted` | study-weighted detection prevalence over the pool |
| `prevalence_band_<b>`, `prevalence_ci_low_<b>`, `prevalence_ci_high_<b>`, `n_band_<b>` | the same within quality band `b` (Wilson 95 % interval); `prevalence_band_<b> >= 0.70` defines "expected" for the expected-but-missing call |
| `p1 … p99` | study-weighted percentiles of the value among samples in which the feature was detected (`p2_5` = 2.5th, `p97_5` = 97.5th) |
| `p<q>_band_<b>`, `n_detected_band_<b>` | assembly gene layers only: the same percentiles within band `b`, used when `percentiles_within_band_layers` lists the layer |
| `p2_5_ci_low/high`, `p97_5_ci_low/high` | study-bootstrap 95 % intervals on the two band edges |
| `percentiles_available` | false when fewer than the minimum number of detections (percentile columns are then null and the feature is `not_assessable`) |
| `layer`, `reliability_tier` | layer name; reserved |

The value column that the percentiles describe is per layer: `clr` for taxonomy and the HUMAnN layers, `copies_per_genome`
for KO and Pfam in assembly bundles, `completeness` (0–1) for KEGG modules.

## `landscape/`

`pca_loadings.parquet` (`feature_id`, `PC1..PC10`) holds the loadings of a PCA on the reference's family-level CLR values
(features with prevalence ≥ the landscape threshold); `reference_coords.parquet` (`study_bioproject`, `PC1..PC10`,
`core_distance`, `knn_distance`, `quality_band`) the coordinates of every reference sample, without sample ids;
`core_distance.json` the centroid, column means, fill value for undetected features and percentile tables of the distances.
A new sample is projected with the reference's centring, and its distance from the centroid is reported as a percentile among
reference samples.

## `normalization/`

Assembly bundles:

| file | content |
|---|---|
| `rules.json` | taxonomy rule (contig majority 0.5, length × coverage weighting, primary rank family), CLR ranks and the multiplicative-replacement delta, genome-equivalent definition, gene value (`copies_per_genome = coverage sum / genome equivalents`), KO source, quality floors (total length, N50, CDS count, marker KOs, genome equivalents, mapped fraction) and the band bounds |
| `marker_panel.tsv` | the 52 validated single-copy marker KOs whose median coverage gives genome equivalents |
| `clr_basis_family.txt`, `clr_basis_genus.txt` | the feature sets the CLR is taken over (584 families, 1,690 genera) |
| `backbone/ncbi_taxonomy.parquet`, `backbone/ncbi_merged.parquet` | the NCBI taxonomy (taxid, parent, rank, name, lineage taxids per rank) and merged-id table the reference used |

Read-based bundles:

| file | content |
|---|---|
| `rules.json` | depth bands, CLR delta, family mapped-fraction floor (0.9), basis sizes (254 species, 119 genera, 49 families), held-out healthy studies |
| `rules_function.json` | the same for the HUMAnN layers (4,518 KO, 423 pathways) and the KO regrouping rule |
| `clr_basis_<species|genus|family|ko_humann|pathway_humann>.txt` | basis feature lists |
| `species_taxid_map.tsv` | `species_name` (MetaPhlAn 3 `s__` name), `species`, `genus`, `family` taxids as the reference placed them |
| `map_ko_uniref90.txt.gz` | HUMAnN's UniRef90 → KO map used to regroup gene families |
| `backbone/` | as above |

## `provenance/`

`studies.tsv` lists the pool studies with sample counts per band; `pool_exclusions.tsv` the candidate samples left out
and why (`CASE_CONTROL_STUDY`, `HELD_OUT`, `LOW_MAPPED_FRACTION`, `OUT_OF_STRATUM`, quality floors). Sample identifiers in
`pool_exclusions.tsv` are public run or analysis accessions.

## Versioning

The bundle id carries the version. Scores are only comparable across samples scored against the same bundle id;
`summaries.tsv` records it. The released tarballs round-trip: unpacking and re-scoring the packaged test queries
reproduces `scores.parquet` exactly (recorded in the release manifest).
