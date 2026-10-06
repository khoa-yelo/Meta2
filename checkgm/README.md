# checkGM — a healthy reference baseline for gut metagenomes

checkgm reports a new human gut metagenome as percentiles against a reference of healthy adults that were measured in
exactly the same way. The idea is the one genomics takes for granted: a variant is only interpretable against a
reference genome. In microbiome analysis there has been no equivalent, so a taxon, gene or pathway abundance in one
sample has had no scale to be read on. checkgm supplies that scale for two standard measurement pipelines and ships each
pipeline together with its own baseline.

| pipeline | measurement | reference baseline (bundle) | reference pool | layers scored |
|---|---|---|---|---|
| **assembly-based** | reproduces MGnify pipeline v5.0 (metaSPAdes/MEGAHIT, Prodigal + FragGeneScan, DIAMOND vs UniRef90, eggNOG-mapper KO, Pfam, KOfam, KEGG module completeness) | `gut-assembly-adult-global-v0.8-{lenient,strict}` | 1,941 healthy adults, 26 studies | family, genus, KO (eggNOG), KO (KOfam), Pfam, KEGG module |
| **read-based** | reproduces curatedMetagenomicData 3 (MetaPhlAn 3.0.14, marker database mpa_v30_CHOCOPhlAn_201901; optional HUMAnN 3) | `gut-reads-adult-global-v0.2-{lenient,strict}` | 6,494 healthy adults, 19 studies | family, genus, species; KO and pathway when HUMAnN 3 tables are given |

A *bundle* is a directory of summary statistics (per-feature prevalence and percentile grids, a PCA landscape of the
reference) with the normalization rules the reference was built with. It contains no sample-level data. Each baseline
comes in two variants: the **standard baseline** (bundle name `lenient`), used throughout and the one to start with, and
the **strict baseline**, which further excludes reported antibiotic use, age above 65, BMI outside 18.5-30 and pregnancy
and serves as a sensitivity analysis.

The output for a sample is, per layer, the percentile of every detected feature within the reference distribution, a
call (`within`, `low`, `high`, or `expected_but_missing`), and one sample-level test: whether the fraction of features
outside the central 95 % band exceeds the 5 % a healthy sample is expected to show.

## Install

```bash
pip install git+https://github.com/khoa-yelo/Meta2.git@checkgm#subdirectory=checkgm
checkgm --version          # checkGM 0.8.2
```

Requirements: Python 3.10 or later; numpy, pandas, scipy, pyarrow, pyyaml and tabulate are installed automatically.
To work from a clone: `cd checkgm && pip install -e .[test] && pytest` (30 tests, about 40 s).

The package only normalizes and scores. Producing the measurement from raw reads needs the pipeline containers below.

## Get a reference bundle

Bundles, the assembly-based pipeline image and the baseline pool accession lists are distributed outside GitHub
(to be deposited on Zenodo together with the container image `checkgm_tierA.sif`, DOI [TBD]). File names, sizes and checksums:

| file | bytes | sha256 |
|---|---:|---|
| `gut-assembly-adult-global-v0.8-lenient.tar.gz` (standard assembly-based baseline) | 54,519,283 | `0d31c959dad8ff7607afc0f2baaf04bc2d7e419cdf1dc5617d4ad3eb229d90f8` |
| `gut-assembly-adult-global-v0.8-strict.tar.gz` (strict) | 54,318,945 | `ddffa36865e3c971b93a826bf47694f72809f935f637a20955bdb9bab363bd5c` |
| `gut-reads-adult-global-v0.2-lenient.tar.gz` (standard read-based baseline) | 64,799,615 | `ab6b82bdc62b72dcfa860ea2dbf0e7d6e750c8058723721c3b2e7f25415b5c0e` |
| `gut-reads-adult-global-v0.2-strict.tar.gz` (strict) | 64,771,253 | `4742e992125786a8286b030510dd6c45e727fd29620a7e6d4ae5588eaefb8734` |
| `checkgm_tierA.sif` (the assembly-based pipeline image, Apptainer) | 783,478,784 | `e7b0de48eef988a8d9f420918a90da89e1315637850ada6fe93d037e890075a6` |
| `pool_gut-assembly-adult-global-v0.8-lenient.tsv` (1,941 public accessions: study, MGnify analysis and assembly, ENA runs, sample, BioSample, size class) | 169,843 | `a1b3ca2bfda3261b74c9c8fad05c1a0fa1c469866c31e956f2dec1c00aaa79b7` |
| `pool_gut-assembly-adult-global-v0.8-strict.tsv` (1,588 accessions) | 136,627 | `75889428751fa5fdb623dc917572325639d00fb0c28ea1ff721bbc81d070a5fe` |
| `heldout_healthy_gut-assembly-adult-global-v0.8.tsv` (233 accessions of the three healthy cohorts kept out for validation) | 29,678 | `45e26202d15b20ddca3c4aac2714ed205b15719bab74168190976df81d777de7` |
| `pool_gut-reads-adult-global-v0.2-lenient.tsv` (6,494 curatedMetagenomicData 3 sample ids with NCBI run accessions where cMD3 records them, depth band) | 468,292 | `095d7747f8d02ae00f1af07c006c102430f6234b5cb9f35cb9b032eb32ba4b90` |
| `pool_gut-reads-adult-global-v0.2-strict.tsv` (6,235 sample ids) | 440,077 | `e20d399da81a7bf805885a421615d9fb53b13a92feba918d84f4f8309d9ee762` |
| `heldout_healthy_gut-reads-adult-global-v0.2.tsv` (317 sample ids of the three held-out healthy cohorts) | 24,269 | `033dc67590bd77bcaa9730294d451ee0c4e5c86006566591cc2b1f078f578473` |

```bash
mkdir -p ~/checkgm_bundles && cd ~/checkgm_bundles
sha256sum -c <<< "ab6b82bdc62b72dcfa860ea2dbf0e7d6e750c8058723721c3b2e7f25415b5c0e  gut-reads-adult-global-v0.2-lenient.tar.gz"
tar xzf gut-reads-adult-global-v0.2-lenient.tar.gz
export CHECKGM_BUNDLES=~/checkgm_bundles     # bundle ids are looked up here; a directory path works without it
checkgm bundles
```

Earlier bundle versions (assembly v0.7 and below, reads v0.1) are superseded and should not be used. The pool lists
hold public identifiers only (no measurements); `docs/bundle_format.md` describes them.

## Five-minute example (read-based pipeline)

The example files live in the repository, so clone it first (the `pip install` above does not copy them):

```bash
git clone -b checkgm https://github.com/khoa-yelo/Meta2.git && cd Meta2/checkgm
```

`examples/synthetic_healthy_adult.txt` is a MetaPhlAn 3 profile in the exact output format, but of no real person: it
was generated from the read-based bundle's own summary statistics (every basis species with prevalence at least 0.3,
at its median abundance), so a healthy sample's profile can be shown without redistributing anyone's data.
`examples/reads.tsv` gives its read count the way a merged table would need it (single-sample MetaPhlAn files carry the
count in their header, as this one does, so the option is redundant here).

```bash
checkgm run-metaphlan --bundle gut-reads-adult-global-v0.2-lenient \
    --input examples/synthetic_healthy_adult.txt --reads-tsv examples/reads.tsv --out example_report
cat example_report/report.md
```

The bundle `gut-reads-adult-global-v0.2-lenient` must already be unpacked under `$CHECKGM_BUNDLES` (previous section);
the baselines and the container image will be deposited at Zenodo, DOI [TBD]. Expected: the sample is scored on the family, genus and species layers (quality band `high`, 45 million reads), with
28, 56 and 102 features assessed and 187 score rows in total; the fraction of features outside the reference band is
0.000 on every layer and no layer shows a significant excess, as it must be for a profile sitting at the reference
medians. One family (Coprobacillaceae) is listed as expected but missing, because the synthetic profile carries no
species of that family although at least 70 % of healthy adults in the same depth band do. The command took 16 s wall time on one CPU core
(`docs/validation.md`, section 6). A real sample has a few features in the tails; the held-out healthy cohorts of
`docs/validation.md` show what that looks like. If you have HUMAnN 3 tables for the same sample, add
`--humann-genefamilies S_genefamilies.tsv --humann-pathabundance S_pathabundance.tsv` to score the KO and pathway
layers as well.

Already have a MGnify v5 assembly analysis (the MGYA download folder)? Score it against the assembly-based bundle:

```bash
checkgm run-mgnify --bundle gut-assembly-adult-global-v0.8-lenient --analysis-dir MGYA00XXXXXX --out mgya_report
```

Mistakes are reported as one line (`checkgm: ...`) with exit status 1: a bundle that cannot be found, a read-based bundle
given to an assembly command or the other way round, a missing input file, or a file that is not a MetaPhlAn table.

## From raw reads

- **Read-based pipeline:** `container/checkgm_pipelineB.def` builds an Apptainer image with MetaPhlAn 3.0.14, Bowtie2 and
  checkGM (`cd container && apptainer build checkgm_pipelineB.sif checkgm_pipelineB.def`). `container/fetch_dbs_B.sh` fetches
  the marker database (0.4 GB download, about 3 GB with the index; the archive's md5 is pinned in the script);
  `container/run_pipelineB.sh` runs profile, normalize and score for one sample. Both mates are given to MetaPhlAn as
  unpaired reads and there is no read QC step, which is how the reference profiles were produced.
- **Assembly-based pipeline:** `container/checkgm_tierA.def` (image `checkgm_tierA.sif`, see the table above; build with
  `cd container && apptainer build checkgm_tierA.sif checkgm_tierA.def`) pins the MGnify v5 tool versions.
  `container/fetch_dbs.sh` and `unpack_dbs.sh` fetch the reference databases (about 70 GB compressed, 110 GB unpacked;
  downloaded once). `container/run_tierA.sh` goes from reads to a *query directory* (`docs/query_format.md`), which
  `checkgm normalize` + `checkgm score` turn into a report. Assembly needs roughly 16 CPUs and 128 GB of memory per sample.
- **Scoring-only image:** `container/Dockerfile` installs just the package (`podman build -t checkGM -f container/Dockerfile .`
  from this directory).

## Reading a report

`checkgm score` and the `run-*` commands write to `--out`:

| file | content |
|---|---|
| `report.md` | human-readable summary: fraction outside the band per layer, samples with a significant excess, landscape position, the ten most extreme features per sample (taxa with their name next to the NCBI taxid), expected-but-missing features by name |
| `scores.parquet` | one row per sample × layer × feature: `value`, `percentile`, `call` (`within`/`low`/`high`/`expected_but_missing`/`not_assessable`), `call_fdr`, `p_two_sided`, `fdr_q`, `band` |
| `summaries.tsv` | one row per sample × layer: `n_assessed`, `n_low_raw`, `n_high_raw`, `n_missing`, `frac_outside_raw`, `excess_outside_p`, `excess_outside_q` (Benjamini–Hochberg across samples within a layer), landscape coordinates, quality band, bundle id |
| `rejections.tsv` | samples not scored: `reason_code` is `NOT_NORMALIZED` (a normalization floor failed) or `NO_MATCHING_REFERENCE` (body site or depth band without a reference); the `detail` column names the floor, for example `EMPTY_PROFILE`, `LOW_MAPPED_FRACTION` or `BELOW_MIN_ASSEMBLY_QUALITY` |
| `qc.tsv`, `normalized/` | the normalized tables that were scored (written by the import step) |

How to read it: a healthy adult sample measured like the reference has about 5 % of its features outside the 2.5–97.5
percentile band; the held-out healthy cohorts used for validation sit between 2 % and 13 % depending on layer and
cohort (most 3-8 %; `docs/validation.md`, section 2). The sample-level excess test (`excess_outside_q ≤ 0.05`) is the
primary statement. Percentiles are the robust quantity; individual low/high calls are descriptive and less reproducible
between re-measurements of the same reads than percentiles are.

## What has been validated, and what has not

All numbers below come from the reports in `docs/reports/` (see `docs/validation.md` for the map).

- Re-measuring public samples from raw reads with the assembly-based pipeline container and scoring them reproduces
  scores obtained from the original MGnify analyses (61 samples; percentile Spearman 0.98–0.99 on gene layers, 0.86–0.88
  on taxonomy). The read-based pipeline reproduces the curatedMetagenomicData 3 profiles from raw reads (34 samples,
  median species Bray–Curtis 0.011, score percentile Spearman ≥ 0.997). The read-based function layers (HUMAnN 3) have
  **not** been validated from raw reads.
- Held-out healthy cohorts mostly stay near the expected 5 % outside the band, but not all do: one assembly cohort
  reaches 13 % on the gene layers, cohorts from Africa do not fit the assembly-based baseline (19.8 % of families
  outside), and location is confounded with study throughout.
- Reference-relative features separate cases from controls only modestly better than raw abundances (assembly-based
  pipeline, leave-one-study-out AUROC 0.618 vs 0.581, p = 0.021 uncorrected; not significant on the strict baseline,
  p = 0.051; read-based taxonomy 0.690 vs 0.664, p = 0.052; colorectal cancer only 0.769 vs 0.715, p = 0.004). Function
  layers show no advantage. Of 43 known-biology expectations tested on the assembly-based pipeline, 14 were confirmed.
- The baselines are adult, stool, global. Other body sites and infants are refused (`NO_MATCHING_REFERENCE`).

## Documentation

- `docs/cli.md` — every subcommand, its arguments, outputs and exit status
- `docs/query_format.md` — the per-sample query directory any assembly pipeline can produce
- `docs/bundle_format.md` — manifest fields, `features/*.parquet` columns, `normalization/`
- `docs/rebuild.md` — how a baseline is rebuilt
- `docs/validation.md` — how validation was done, with the report files in `docs/reports/`
- `docs/configs/` — read-only copies of the baseline-build and evaluation settings (thresholds, seeds)
- `scripts/viz/` — the scripts that draw every figure of the paper (they read the project workspace, see its README)
- `CHANGELOG.md`

## Development notes

- The early example profile `checkgm/examples/CosteaPI_2017_alien2-11-0-0.txt` (a public curatedMetagenomicData 3 profile,
  commits `97cba6b`, `92c85ce`) stays in the history of branch `checkgm` by the owner's decision of 2026-10-02 (`CHANGELOG.md`);
  the example shipped since 0.8.0 is synthetic.
- Tests: `pip install -e .[test] && pytest` (30 tests, about 40 s). CI runs them on Python 3.10 and 3.12
  (`.github/workflows/test.yml` at the repository root).
- `scripts/viz/` holds byte-identical copies of the project workspace's figure scripts; they are not edited here.

## Citation

Paper, authors and DOI to be announced [TBD]; see `CITATION.cff`. Until then please cite this repository.

## Licence

MIT (see `LICENSE`). `container/emapper2/` is eggNOG-mapper 2.0.0 (GPL v2) as vendored by MGnify pipeline v5 and is
used only inside the assembly-based pipeline container; `container/pv5_ref/` holds MGnify pipeline-v5 workflow files kept
for provenance.
