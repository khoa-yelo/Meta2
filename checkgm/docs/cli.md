# Command line

```
checkgm --version
checkgm bundles
checkgm normalize        --bundle B --query DIR [DIR ...] --out OUT
checkgm import-mgnify    --bundle B --analysis-dir DIR [DIR ...] --out OUT
checkgm import-metaphlan --bundle B --input PROFILE [PROFILE ...] [--reads-tsv FILE] [--humann-genefamilies FILE ...] [--humann-pathabundance FILE ...] --out OUT
checkgm score            --bundle B --normalized OUT [--calibration DIR] --out REPORT
checkgm run-mgnify       --bundle B --analysis-dir DIR [DIR ...] [--calibration DIR] --out REPORT
checkgm run-metaphlan    --bundle B --input PROFILE [PROFILE ...] [--reads-tsv FILE] [--humann-genefamilies FILE ...] [--humann-pathabundance FILE ...] [--calibration DIR] --out REPORT
```

`B` is either a bundle directory (one containing `manifest.json`) or a bundle id such as
`gut-reads-adult-global-v0.2-lenient`, which is looked up in the directories listed in the environment variable
`CHECKGM_BUNDLES` (colon-separated). An unknown bundle is reported with the places that were searched. Every bundle is
either assembly-based (`gut-assembly-*`, `profile_type: assembly`) or read-based (`gut-reads-*`, `profile_type: reads`);
`normalize`, `import-mgnify` and `run-mgnify` need an assembly-based bundle, `import-metaphlan` and `run-metaphlan` a
read-based one, and `score` checks the bundle against what the normalized directory was made from. A mismatch is
reported and the command exits 1.

`checkgm --help` lists the subcommands and `checkgm <command> --help` explains every argument in one sentence; the text
below says the same at greater length.

Every command has two halves: an **import** step that turns a measurement into normalized tables (`OUT/normalized/<sample>.parquet`,
`OUT/qc.tsv`, `OUT/sample_meta.json`) and a **score** step that compares those tables with the bundle. The `run-*` commands do
both into one directory.

## `checkgm bundles`

Lists the bundles found under `$CHECKGM_BUNDLES`: id, profile type (`assembly` or `reads`), pool size and path.

## `checkgm normalize` (assembly-based, any pipeline)

Input: one *query directory* per sample in the format of `docs/query_format.md` (`contigs.tsv`, `cds.tsv`, `cds_taxonomy.tsv`,
optional `modules.tsv` and `sample.json`). The sample id is `sample.json: sample_id`, else the directory name. The
normalization rules come from the bundle (`normalization/rules.json`): contig-majority taxonomy on the bundle's NCBI
backbone, length × coverage weighting, CLR over the bundle's family and genus basis, copies per genome for KO and Pfam from
the bundle's marker panel, quality band from total assembled length. Prints `sample status band` per sample.

Only meaningful against an assembly bundle (`gut-assembly-adult-global-*`).

## `checkgm import-mgnify` (assembly-based, MGnify v5 analyses)

Input: MGnify pipeline v5 assembly analysis download folders (one per sample, named by the MGYA accession). The same
normalization as above, read from MGnify's own output files. This is exactly how the reference pool was measured, so no
calibration is involved.

## `checkgm import-metaphlan` (read-based)

Input: MetaPhlAn 3 single-sample profiles (header lines starting with `#`, columns `clade_name`, `NCBI_tax_id`,
`relative_abundance`) or a merged table (`clade_name [NCBI_tax_id] sample1 sample2 ...`). Only species-level rows are used.
Species are placed on the bundle's taxonomy using the bundle's own species map first, then the file's NCBI taxid, then a
name lookup. Per rank, proportions over the bundle's basis are CLR-transformed with multiplicative replacement of zeros
(the rule set in `normalization/rules.json`).

- `--reads-tsv FILE`: two columns, `sample<TAB>reads processed`, used for the depth band when the profile header has no
  `#N reads processed` line (merged tables never do). A header line is accepted and skipped; rows whose read count is not
  a number are reported on stderr and ignored. If the depth is unknown the sample gets band `unknown`. Read-based bundles
  do not use bands for percentiles, so such a sample is still scored.
- `--humann-genefamilies FILE ...` and `--humann-pathabundance FILE ...`: HUMAnN 3 tables (single sample or joined;
  stratified `feature|species` rows are dropped). Gene families given as UniRef90 are regrouped to KO with the map shipped in
  the bundle; the `UNMAPPED`, `UNGROUPED` and `UNINTEGRATED` rows are removed before the CLR. These add the `ko_humann` and
  `pathway_humann` layers when the bundle has them. These layers have not been validated from raw reads (see
  `docs/validation.md`).

A profile whose species cannot be placed on the backbone to at least 90 % of its abundance (family level) fails the floor
`LOW_MAPPED_FRACTION`; a profile without species rows (for example a header-only file) fails `EMPTY_PROFILE`. Both appear in
`rejections.tsv` with `reason_code` `NOT_NORMALIZED` and the floor name in the `detail` column; the command still exits 0. A file that is not a MetaPhlAn table at all (no `clade_name` column) is a
user error: one line on stderr, exit 1.

## `checkgm score`

Reads the normalized tables from `--normalized OUT` and writes to `--out REPORT`:

| file | content |
|---|---|
| `scores.parquet` | one row per sample × layer × feature: `analysis_id`, `layer`, `feature_id`, `band`, `value`, `percentile`, `call`, `call_fdr`, `p_two_sided`, `fdr_q` |
| `summaries.tsv` | one row per sample × layer, 31 columns: `analysis_id`, `layer`; counts (`n_detected`, `n_assessed`, `n_low_raw`, `n_high_raw`, `n_low`, `n_high`, `n_missing`, `n_not_assessable`); `frac_outside_raw`; `expected_false_positives_at_alpha` = 0.05 × `n_assessed`, the number of features a healthy sample is expected to show outside the band by chance; `weighted_deviation_score`; the excess test `excess_outside_p`, `excess_outside_ratio` (observed fraction outside / 0.05) and `excess_outside_q`; `layer_status` = `scored`, or `PIPELINE_NOT_CALIBRATED` when a calibration was given that does not cover the layer (counts are then 0 and fractions empty); `quality_band` (depth band or assembly-size class the sample was scored in); `mapped_fraction_family` = share of the profile's abundance that maps to a family of the bundle's basis (read-based profiles; the normalization floor `LOW_MAPPED_FRACTION` uses it); `genome_equivalents_cov` = median coverage of the single-copy marker KOs in genome equivalents (assembly-based profiles; empty for read-based bundles); `bundle_id`; `calibration_applied` (`none` or the calibration id); `fallback_used` (always `False` in this version: no fallback scoring path exists, the column is reserved); `landscape_centring` = `reference` when the family CLR is centred with the column means stored in the bundle, `query_batch` when an older bundle stores none and the query batch's own means are used; landscape columns `core_distance`, `core_distance_pct`, `neighbor_studies` (JSON: study → count among the k nearest reference samples), `neighbor_bands` (JSON: depth band or size class → count among the same neighbours), and the sample's coordinates `PC1`, `PC2`, `PC3` |
| `rejections.tsv` | `analysis_id`, `reason_code`, `detail` for samples not scored. `reason_code` is `NOT_NORMALIZED` (a normalization floor failed; `detail` names the floor: `EMPTY_PROFILE`, `LOW_MAPPED_FRACTION`, `BELOW_MIN_ASSEMBLY_QUALITY`, `NO_MARKER_PANEL_SIGNAL`) or `NO_MATCHING_REFERENCE` (body site or depth band without a reference) |
| `report.md` | the readable summary: fraction outside the band per layer, samples with a significant excess, landscape position, the ten most extreme features per sample (taxa with their NCBI name next to the taxid), and the expected-but-missing features listed by name |

Definitions:

- `percentile`: position of the sample's value within the reference distribution of that feature (linear interpolation on
  the bundle's percentile grid 1, 2.5, 5, 10, 25, 50, 75, 90, 95, 97.5, 99; values beyond the grid are placed at 0.5 or 99.5).
  For assembly bundles the gene layers use the percentile grid of the sample's quality band.
- `call`: `low` below the reference 2.5th percentile, `high` above the 97.5th, `within` otherwise; `not_assessable` when the
  feature has no percentiles in the bundle or the value is under the detection floor (0.05 copies per genome on gene layers
  of assembly bundles); `expected_but_missing` when a feature present in at least 70 % of reference samples of the same
  quality band was not detected.
- `call_fdr`: `low`/`high` only if the two-sided tail probability survives Benjamini–Hochberg across the sample's features
  in that layer; otherwise `within_after_fdr`.
- `excess_outside_p`: binomial test that the number of features outside the band exceeds 5 % of the assessed features;
  `excess_outside_q` is this p-value adjusted across samples within a layer. This is the primary statement about a sample.
- Landscape: the sample's coordinates on the reference's PCA of family CLR values, its distance from the reference centroid
  and the percentile of that distance among reference samples (`core_distance_pct`), and the studies of its ten nearest
  reference samples.

`--calibration DIR` applies a calibration object produced with the `checkgm.calibrate` module (for a measurement pipeline
other than the reference's; see the module docstring). No calibration is shipped with this release; the two pipelines in
`container/` reproduce the reference measurement and need none.

## Exit status and errors

User errors exit 1 with a single line on stderr of the form `checkgm: <what is wrong>` and no traceback: a bundle that
cannot be found (the searched places are named), a bundle of the wrong kind for the command (assembly-based vs read-based),
a missing `--input` file, `--analysis-dir`, `--query` or `--calibration` directory, a `--normalized` directory without
`qc.tsv`, or an `--input` file that does not have MetaPhlAn's columns. Bad command-line syntax exits 2 (argparse).

Scientific refusals are not errors: a sample that fails a normalization floor (`EMPTY_PROFILE`, `LOW_MAPPED_FRACTION`,
assembly quality floors, named in the `detail` column) or has no matching reference appears in `rejections.tsv` with its
`reason_code` (`NOT_NORMALIZED` or `NO_MATCHING_REFERENCE`) and the command exits 0. A bundle from a superseded series
(read-based v0.1, assembly-based v0.7 and below) is still scored, with one warning line on stderr. Unparseable rows of `--reads-tsv` are reported on stderr and skipped, also with exit 0.

The two-step workflow (`import-metaphlan` then `score`) and the one-step `run-metaphlan` write identical files: `qc.tsv`
stores floats with 17 significant digits and is read back exactly.

## Library use

```python
from checkgm.normalize_reads import ReadsBundleSpec, read_metaphlan, normalize_profile
from checkgm.score import Bundle, score
from checkgm.paths import resolve_bundle

b = resolve_bundle("gut-reads-adult-global-v0.2-lenient")
spec = ReadsBundleSpec(b)
parsed = read_metaphlan("S1.txt")                      # {sample_id: {"table": DataFrame, "reads": int | None}}
norm = {sid: normalize_profile(sid, d["table"], d["reads"], spec) for sid, d in parsed.items()}
result = score(norm, Bundle(b))                        # {"scores", "summaries", "landscape", "rejections"} DataFrames
```

`checkgm.build_reference.build(...)` builds a bundle from long tables (see `docs/rebuild.md`); `checkgm.calibrate` fits a
calibration between a user's pipeline and the reference's using samples measured both ways.
