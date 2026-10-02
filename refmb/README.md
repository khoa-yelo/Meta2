# refmb — a healthy reference baseline for gut metagenomes

refmb reports a new human gut metagenome as percentiles against a reference of healthy adults that were measured in
exactly the same way. The idea is the one genomics takes for granted: a variant is only interpretable against a
reference genome. In microbiome analysis there has been no equivalent, so a taxon, gene or pathway abundance in one
sample has had no scale to be read on. refmb supplies that scale for two standard measurement pipelines and ships each
pipeline together with its own baseline.

| pipeline | measurement | reference baseline (bundle) | reference pool | layers scored |
|---|---|---|---|---|
| **assembly-based** | reproduces MGnify pipeline v5.0 (metaSPAdes/MEGAHIT, Prodigal + FragGeneScan, DIAMOND vs UniRef90, eggNOG-mapper KO, Pfam, KOfam, KEGG module completeness) | `gut-assembly-adult-global-v0.8-{lenient,strict}` | 1,941 healthy adults, 26 studies | family, genus, KO (eggNOG), KO (KOfam), Pfam, KEGG module |
| **read-based** | reproduces curatedMetagenomicData 3 (MetaPhlAn 3.0.14, marker database mpa_v30_CHOCOPhlAn_201901; optional HUMAnN 3) | `gut-reads-adult-global-v0.2-{lenient,strict}` | 6,494 healthy adults, 19 studies | family, genus, species; KO and pathway when HUMAnN 3 tables are given |

A *bundle* is a directory of summary statistics (per-feature prevalence and percentile grids, a PCA landscape of the
reference) with the normalization rules the reference was built with. It contains no sample-level data. *Lenient* and
*strict* differ in how strictly "healthy adult" was defined when the pool was assembled; lenient is the primary baseline,
strict is the sensitivity analysis.

The output for a sample is, per layer, the percentile of every detected feature within the reference distribution, a
call (`within`, `low`, `high`, or `expected_but_missing`), and one sample-level test: whether the fraction of features
outside the central 95 % band exceeds the 5 % a healthy sample is expected to show.

## Install

```bash
pip install git+https://github.com/khoa-yelo/Meta2.git@refmb#subdirectory=refmb
refmb --version          # refmb 0.8.0
```

Requirements: Python 3.10 or later; numpy, pandas, scipy, pyarrow, pyyaml and tabulate are installed automatically.
To work from a clone: `cd refmb && pip install -e .[test] && pytest`.

The package only normalizes and scores. Producing the measurement from raw reads needs the pipeline containers below.

## Get a reference bundle

Bundles and the pipeline A container image are distributed outside GitHub (download location to be announced [TBD]).
File names, sizes and checksums:

| file | bytes | sha256 |
|---|---:|---|
| `gut-assembly-adult-global-v0.8-lenient.tar.gz` | 54,519,283 | `0d31c959dad8ff7607afc0f2baaf04bc2d7e419cdf1dc5617d4ad3eb229d90f8` |
| `gut-assembly-adult-global-v0.8-strict.tar.gz` | 54,318,945 | `ddffa36865e3c971b93a826bf47694f72809f935f637a20955bdb9bab363bd5c` |
| `gut-reads-adult-global-v0.2-lenient.tar.gz` | 64,799,615 | `ab6b82bdc62b72dcfa860ea2dbf0e7d6e750c8058723721c3b2e7f25415b5c0e` |
| `gut-reads-adult-global-v0.2-strict.tar.gz` | 64,771,253 | `4742e992125786a8286b030510dd6c45e727fd29620a7e6d4ae5588eaefb8734` |
| `refmb_tierA.sif` (pipeline A container, Apptainer) | 783,478,784 | `e7b0de48eef988a8d9f420918a90da89e1315637850ada6fe93d037e890075a6` |

```bash
mkdir -p ~/refmb_bundles && cd ~/refmb_bundles
sha256sum -c <<< "ab6b82bdc62b72dcfa860ea2dbf0e7d6e750c8058723721c3b2e7f25415b5c0e  gut-reads-adult-global-v0.2-lenient.tar.gz"
tar xzf gut-reads-adult-global-v0.2-lenient.tar.gz
export REFMB_BUNDLES=~/refmb_bundles     # bundle ids are looked up here; a directory path works without it
refmb bundles
```

Earlier bundle versions (assembly v0.7 and below, reads v0.1) are superseded and should not be used.

## Five-minute example (read-based pipeline)

`examples/CosteaPI_2017_alien2-11-0-0.txt` is a MetaPhlAn 3 profile of a public stool sample from a healthy adult
(CosteaPI_2017, a study held out of the reference pool), and `examples/reads.tsv` gives its read count (the depth band
comes from reads processed; MetaPhlAn single-sample files carry it in their header, merged tables do not).

```bash
refmb run-metaphlan --bundle gut-reads-adult-global-v0.2-lenient \
    --input examples/CosteaPI_2017_alien2-11-0-0.txt --reads-tsv examples/reads.tsv --out example_report
cat example_report/report.md
```

Expected: the sample is scored on the family, genus and species layers (quality band `high`, 70 million reads); the fraction
of features outside the reference band is 0.111 (family), 0.060 (genus) and 0.073 (species), and no layer shows a
significant excess. That is what a healthy sample looks like: a few features in the tails, no more than chance allows. If you have
HUMAnN 3 tables for the same sample, add `--humann-genefamilies S_genefamilies.tsv --humann-pathabundance S_pathabundance.tsv`
to score the KO and pathway layers as well.

Already have a MGnify v5 assembly analysis (the MGYA download folder)? Score it against the assembly bundle:

```bash
refmb run-mgnify --bundle gut-assembly-adult-global-v0.8-lenient --analysis-dir MGYA00XXXXXX --out mgya_report
```

## From raw reads

- **Read-based (pipeline B):** `container/refmb_pipelineB.def` builds an Apptainer image with MetaPhlAn 3.0.14, Bowtie2 and
  refmb. `container/fetch_dbs_B.sh` fetches the marker database (0.4 GB download, about 3 GB with the index);
  `container/run_pipelineB.sh` runs profile → normalize → score for one sample. Both mates are given to MetaPhlAn as
  unpaired reads and there is no read QC step, which is how the reference profiles were produced.
- **Assembly-based (pipeline A):** `container/refmb_tierA.def` (image `refmb_tierA.sif`, see the table above) pins the
  MGnify v5 tool versions. `container/fetch_dbs.sh` and `unpack_dbs.sh` fetch the reference databases (about 70 GB
  compressed, 110 GB unpacked; downloaded once). `container/run_tierA.sh` goes from reads to a *query directory*
  (`docs/query_format.md`), which `refmb normalize` + `refmb score` turn into a report. Assembly needs roughly 16 CPUs and
  128 GB of memory per sample.
- **Scoring-only image:** `container/Dockerfile` installs just the package (`podman build -t refmb -f container/Dockerfile .`
  from this directory).

## Reading a report

`refmb score` and the `run-*` commands write to `--out`:

| file | content |
|---|---|
| `report.md` | human-readable summary: fraction outside the band per layer, samples with a significant excess, landscape position, the ten most extreme features per sample, expected-but-missing counts |
| `scores.parquet` | one row per sample × layer × feature: `value`, `percentile`, `call` (`within`/`low`/`high`/`expected_but_missing`/`not_assessable`), `call_fdr`, `p_two_sided`, `fdr_q`, `band` |
| `summaries.tsv` | one row per sample × layer: `n_assessed`, `n_low_raw`, `n_high_raw`, `n_missing`, `frac_outside_raw`, `excess_outside_p`, `excess_outside_q` (Benjamini–Hochberg across samples within a layer), landscape coordinates, quality band, bundle id |
| `rejections.tsv` | samples not scored, with a reason code (`NOT_NORMALIZED`, `NO_MATCHING_REFERENCE`) |
| `qc.tsv`, `normalized/` | the normalized tables that were scored (written by the import step) |

How to read it: a healthy adult sample measured like the reference has about 5 % of its features outside the 2.5–97.5
percentile band; the held-out healthy cohorts used for validation sit between 3 % and 10 %. The sample-level excess test
(`excess_outside_q ≤ 0.05`) is the primary statement. Percentiles are the robust quantity; individual low/high calls are
descriptive and less reproducible between re-measurements of the same reads than percentiles are.

## What has been validated, and what has not

All numbers below come from the reports in `docs/reports/` (see `docs/validation.md` for the map).

- Re-measuring public samples from raw reads with the pipeline A container and scoring them reproduces scores obtained
  from the original MGnify analyses (61 samples; percentile Spearman 0.98–0.99 on gene layers, 0.86–0.88 on taxonomy).
  Pipeline B reproduces the curatedMetagenomicData 3 profiles from raw reads (34 samples, median species Bray–Curtis 0.011,
  score percentile Spearman ≥ 0.997). The pipeline B function layers (HUMAnN 3) have **not** been validated from raw reads.
- Held-out healthy cohorts mostly stay near the expected 5 % outside the band, but not all do: cohorts from Africa do not
  fit the assembly baseline (19.8 % of families outside), and location is confounded with study throughout.
- Reference-relative features separate cases from controls only modestly better than raw abundances (pipeline A
  leave-one-study-out AUROC 0.618 vs 0.581, p = 0.021; not significant on the strict baseline, p = 0.051; pipeline B taxonomy
  0.690 vs 0.664, p = 0.052; colorectal cancer only 0.769 vs 0.715, p = 0.004). Function layers show no advantage. Of 43
  known-biology expectations tested on pipeline A, 14 were confirmed.
- The baselines are adult, stool, global. Other body sites and infants are refused (`NO_MATCHING_REFERENCE`).

## Documentation

- `docs/cli.md` — every subcommand, its arguments and outputs
- `docs/query_format.md` — the per-sample query directory any assembly pipeline can produce
- `docs/bundle_format.md` — manifest fields, `features/*.parquet` columns, `normalization/`
- `docs/rebuild.md` — how a baseline is rebuilt
- `docs/validation.md` — how validation was done, with the report files in `docs/reports/`
- `CHANGELOG.md`

## Citation

Paper, authors and DOI to be announced [TBD]; see `CITATION.cff`. Until then please cite this repository.

## Licence

MIT (see `LICENSE`). `container/emapper2/` is eggNOG-mapper 2.0.0 (GPL v2) as vendored by MGnify pipeline v5 and is
used only inside the pipeline A container; `container/pv5_ref/` holds MGnify pipeline-v5 workflow files kept for provenance.
