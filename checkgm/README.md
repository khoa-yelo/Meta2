# checkGM — a reference baseline for gut metagenomes

checkgm reports a new human gut metagenome as percentiles against a reference of healthy adults that were measured in
exactly the same way. The idea is the one genomics takes for granted: a variant is only interpretable against a
reference genome. In microbiome analysis there has been no equivalent, so a taxon, gene or pathway abundance in one
sample has had no scale to be read on. checkgm supplies that scale for two standard measurement pipelines and ships each
pipeline together with its own baseline.

| pipeline | measurement | reference baseline (bundle) | reference pool | layers scored |
|---|---|---|---|---|
| **assembly-based** | reproduces MGnify pipeline v5.0 (metaSPAdes/MEGAHIT, Prodigal + FragGeneScan, DIAMOND vs UniRef90, eggNOG-mapper KO, Pfam, KOfam, KEGG module completeness) | `gut-assembly-adult-global-v0.8-{lenient,strict}` | 1,941 healthy adults, 26 studies | family, genus, KO (eggNOG), KO (KOfam), Pfam, KEGG module |
| **read-based** | reproduces curatedMetagenomicData 3 (MetaPhlAn 3.0.14, marker database mpa_v30_CHOCOPhlAn_201901; optional HUMAnN 3) | `gut-reads-adult-global-v0.2-{lenient,strict}` | 6,494 healthy adults, 19 studies (6,323 distinct specimens: see below) | family, genus, species; KO and pathway when HUMAnN 3 tables are given |

A *bundle* is a directory of summary statistics (per-feature prevalence and percentile grids, a PCA landscape of the
reference) with the normalization rules the reference was built with. It contains no sample-level data. Each baseline
comes in two variants, of which the **standard baseline** (bundle name `lenient`) is used throughout and is the one to
start with. The **strict baseline** serves as a sensitivity analysis. Where the standard tier keeps a sample whose
health status was never recorded, the strict tier requires the health label itself, and that requirement bites, health
status being stated for only 64 % of gut samples. Reported antibiotic use, age above 65, BMI outside 18.5-30 and
pregnancy are excluded on top of it (`docs/configs/inclusion.yaml`).

One fact about the read-based pool is worth stating next to its size. 342 of its 6,494 samples are 171 pairs in which
the same specimen was deposited twice in curatedMetagenomicData 3, as LeChatelierE_2013 and as NielsenHB_2014, two
deposits of one MetaHIT Danish cohort. The pool therefore describes **6,323 distinct specimens**, and those 171 donors
carry twice the weight of the rest in every percentile the read-based bundle stores. The pairs are not byte-identical
copies — each deposit profiles the specimen from a differently partitioned set of sequencing runs — so they were not
de-duplicated, but they are detected and they are held out together whenever one of the two studies is left out
(`docs/validation.md` section 2).

The output for a sample is, per layer, the percentile of every detected feature within the reference distribution, a
call (`within`, `low`, `high`, or `expected_but_missing`), and one sample-level test: whether the fraction of features
outside the central 95 % band exceeds the 5 % a healthy sample is expected to show.

## Install

```bash
pip install git+https://github.com/khoa-yelo/Meta2.git@checkgm#subdirectory=checkgm
checkgm --version          # checkgm 0.8.2
```

Two spellings of one name are used deliberately: checkGM for the project and the paper, `checkgm` for the package, the
console script and the container tag.

Requirements: Python 3.10 or later; numpy, pandas, scipy, pyarrow, pyyaml and tabulate are installed automatically.
To work from a clone: `cd checkgm && pip install -e ".[test]" && pytest` (155 tests, about 100 s). The quotes matter
in zsh, which would otherwise read `.[test]` as a glob. `pip install "checkgm[figures]"` adds matplotlib, which
`checkgm assess` uses to draw the per-sample report figure; without it the tables are written and the figure is
skipped with a note.

The package itself assesses and scores. Producing the measurement from raw reads needs the pipeline containers below.

## Commands

Five commands cover the whole workflow, and `checkgm <command> --help` documents every argument:

| command | what it does |
|---|---|
| `setup` | downloads the reference baselines (`setup bundles`), or a pipeline's reference databases (`setup db`) |
| `profile` | runs a measurement pipeline on raw reads: assembly, or MetaPhlAn 3 |
| `assess` | places a sample against its baseline — percentiles, calls, report and figure |
| `score` | the supervised health score, computed from an `assess` output |
| `end-to-end` | `profile`, then `assess`, then `score`, under one output directory |

`assess` works out what kind of input it was given from what the path holds — an assembly query directory, an MGnify v5
analysis directory, a MetaPhlAn profile, or a directory `assess` itself wrote — so one command covers every entry
point. Re-running it on its own output re-scores without repeating the normalization, which is how a sample is compared
against a second baseline or calibration.

Both commands that can start a very large or very long job (`setup db`, `profile`) check their prerequisites first and
take `--dry-run`, so the size of a 90 GB download or the exact pipeline invocation can be seen before committing to it.

The stages these are built from stay available under their own names (`normalize`, `import-mgnify`, `import-metaphlan`,
`score-normalized`, `run-normalize`, `run-mgnify`, `run-metaphlan`) and are unchanged; the paper's methods refer to
them. One rename matters: `score` now means the health score, and what `score` used to do is `assess`, or
`score-normalized` for the stage on its own. Passing the old flags to `score` prints that redirect.

## Get a reference bundle

Bundles, the assembly-based pipeline image and the baseline pool accession lists are too large for GitHub and are
distributed beside it; a Zenodo deposit covering them and the container image `checkgm_assembly.sif` is prepared, and its
DOI is reserved as `10.5281/zenodo.23181418` and will resolve once the record is public. **That deposit is not public yet**, so until it is, the archives have to be requested
from the authors and given to `setup` directly:

```bash
checkgm setup bundles --dest ~/checkgm-bundles --url /path/holding/the/archives
export CHECKGM_BUNDLES=~/checkgm-bundles
```

`setup bundles` verifies every archive against the checksum recorded below, skips what is already present and correct,
and prints the sizes first with `--dry-run`. Once the deposit is public the same command works without `--url`. The
table is also what a request to the authors should yield, so that whatever arrives can be checked by hand:

| file | bytes | sha256 |
|---|---:|---|
| `gut-assembly-adult-global-v0.8-lenient.tar.gz` (standard assembly-based baseline) | 54,519,283 | `0d31c959dad8ff7607afc0f2baaf04bc2d7e419cdf1dc5617d4ad3eb229d90f8` |
| `gut-assembly-adult-global-v0.8-strict.tar.gz` (strict) | 54,318,945 | `ddffa36865e3c971b93a826bf47694f72809f935f637a20955bdb9bab363bd5c` |
| `gut-reads-adult-global-v0.2-lenient.tar.gz` (standard read-based baseline) | 64,799,615 | `ab6b82bdc62b72dcfa860ea2dbf0e7d6e750c8058723721c3b2e7f25415b5c0e` |
| `gut-reads-adult-global-v0.2-strict.tar.gz` (strict) | 64,771,253 | `4742e992125786a8286b030510dd6c45e727fd29620a7e6d4ae5588eaefb8734` |
| `checkgm_assembly.sif` (the assembly-based pipeline image, Apptainer) | 783,478,784 | `e7b0de48eef988a8d9f420918a90da89e1315637850ada6fe93d037e890075a6` |
| `pool_gut-assembly-adult-global-v0.8-lenient.tsv` (1,941 public accessions: study, MGnify analysis and assembly, ENA runs, sample, BioSample, size class) | 169,843 | `a1b3ca2bfda3261b74c9c8fad05c1a0fa1c469866c31e956f2dec1c00aaa79b7` |
| `pool_gut-assembly-adult-global-v0.8-strict.tsv` (1,588 accessions) | 136,627 | `75889428751fa5fdb623dc917572325639d00fb0c28ea1ff721bbc81d070a5fe` |
| `heldout_healthy_gut-assembly-adult-global-v0.8.tsv` (233 accessions of the three healthy cohorts kept out for validation) | 29,678 | `45e26202d15b20ddca3c4aac2714ed205b15719bab74168190976df81d777de7` |
| `pool_gut-reads-adult-global-v0.2-lenient.tsv` (6,494 curatedMetagenomicData 3 sample ids with NCBI run accessions where cMD3 records them, depth band) | 468,292 | `095d7747f8d02ae00f1af07c006c102430f6234b5cb9f35cb9b032eb32ba4b90` |
| `pool_gut-reads-adult-global-v0.2-strict.tsv` (6,235 sample ids) | 440,077 | `e20d399da81a7bf805885a421615d9fb53b13a92feba918d84f4f8309d9ee762` |
| `heldout_healthy_gut-reads-adult-global-v0.2.tsv` (317 sample ids of the three held-out healthy cohorts) | 24,269 | `033dc67590bd77bcaa9730294d451ee0c4e5c86006566591cc2b1f078f578473` |

```bash
mkdir -p ~/checkgm_bundles && cd ~/checkgm_bundles
# obtain gut-reads-adult-global-v0.2-lenient.tar.gz into this directory first: there is no public URL yet (see above)
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
at its median abundance), so a reference-median profile can be shown without redistributing anyone's data.
`examples/reads.tsv` gives its read count the way a merged table would need it (single-sample MetaPhlAn files carry the
count in their header, as this one does, so the option is redundant here).

```bash
checkgm assess --bundle gut-reads-adult-global-v0.2-lenient \
    --input examples/synthetic_healthy_adult.txt --reads-tsv examples/reads.tsv --out example_report
cat example_report/report.md
```

The bundle `gut-reads-adult-global-v0.2-lenient` must already be unpacked under `$CHECKGM_BUNDLES`, which the previous
section explains how to arrange. The sample is then scored on the family, genus and species layers, in quality band
`high` at 45 million reads, with 28, 56 and 102 features assessed and 187 score rows in total. On every layer the
fraction of features outside the reference band is 0.000 and no layer shows a significant excess, as it must be for a
profile sitting at the reference medians. One family (Coprobacillaceae) is listed as expected but missing, because the
synthetic profile carries no species of that family although at least 70 % of healthy adults in the same depth band do,
and the whole command takes 16 s of wall time on one CPU core (`docs/validation.md`, section 6). A real sample has a
few features in the tails; the held-out healthy cohorts of
`docs/validation.md` show what that looks like. If you have HUMAnN 3 tables for the same sample, add
`--humann-genefamilies S_genefamilies.tsv --humann-pathabundance S_pathabundance.tsv` to score the KO and pathway
layers as well.

Already have a MGnify v5 assembly analysis (the MGYA download folder)? Score it against the assembly-based bundle:

```bash
checkgm assess --bundle gut-assembly-adult-global-v0.8-lenient --input MGYA00XXXXXX --out mgya_report
```

Mistakes are reported as one line (`checkgm: ...`) with exit status 1: a bundle that cannot be found, a read-based bundle
given to an assembly command or the other way round, a missing input file, or a file that is not a MetaPhlAn table.

## From raw reads

- **Read-based pipeline:** `container/checkgm_read.def` builds an Apptainer image with MetaPhlAn 3.0.14, Bowtie2 and
  checkGM (`cd container && apptainer build checkgm_read.sif checkgm_read.def`). `container/fetch_dbs_read.sh` fetches
  the marker database (0.4 GB download, about 3 GB with the index; the archive's md5 is pinned in the script);
  `container/run_read.sh` runs profile, normalize and score for one sample. Both mates are given to MetaPhlAn as
  unpaired reads and there is no read QC step, which is how the reference profiles were produced.
- **Assembly-based pipeline:** `container/checkgm_assembly.def` (image `checkgm_assembly.sif`, see the table above; build with
  `cd container && apptainer build checkgm_assembly.sif checkgm_assembly.def`) pins the MGnify v5 tool versions.
  `container/fetch_dbs_assembly.sh` and `unpack_dbs.sh` fetch the reference databases, downloaded once: about 90 GB over
  the wire and 170 GB on disk, which `checkgm setup db --pipeline assembly --dry-run` prints before transferring
  anything. `container/run_assembly.sh` goes from reads to a *query directory* (`docs/query_format.md`), which
  `checkgm assess --input` turns into a report. `checkgm profile --pipeline assembly` runs this for you and checks the
  image and databases are present first; `checkgm end-to-end` continues straight into the assessment and the score.
  Assembly needs roughly 16 CPUs and 128 GB of memory per sample.
- **Scoring-only image:** `container/Dockerfile` installs just the package, and is built from this directory with
  `podman build -t checkgm:0.8.2 -f container/Dockerfile .` (the tag has to be lowercase, which both podman and docker
  insist on). Its entrypoint is the installed `checkgm` console script, so `podman run checkgm:0.8.2 --help` prints
  the command list.

## Reading a report

`checkgm score` and the `run-*` commands write to `--out`:

| file | content |
|---|---|
| `report.md` | human-readable summary: fraction outside the band per layer, samples with a significant excess, landscape position, the ten most extreme features per sample (taxa with their name next to the NCBI taxid), expected-but-missing features by name |
| `scores.parquet` | one row per sample × layer × feature: `value`, `percentile`, `call` (`within`/`low`/`high`/`expected_but_missing`/`not_assessable`), `call_fdr`, `p_two_sided`, `fdr_q`, `band` |
| `summaries.tsv` | one row per sample × layer: `n_assessed`, `n_low_raw`, `n_high_raw`, `n_missing`, `frac_outside_raw`, `excess_outside_p`, `excess_outside_q` (Benjamini–Hochberg across samples within a layer), landscape coordinates, quality band, bundle id |
| `rejections.tsv` | samples not scored: `reason_code` is `NOT_NORMALIZED` (a normalization floor failed) or `NO_MATCHING_REFERENCE` (body site or depth band without a reference); the `detail` column names the floor, for example `EMPTY_PROFILE`, `LOW_MAPPED_FRACTION` or `BELOW_MIN_ASSEMBLY_QUALITY` |
| `qc.tsv`, `normalized/` | the normalized tables that were scored (written by the import step) |
| `sample_meta.json` | which measurement each normalized sample came from and which bundle normalized it, so that `score` can both refuse a bundle of the wrong kind and warn when the tables are scored against a different bundle than the one they were centred on. The import and `run-*` commands write those two fields alone, as `{"S1": {"pipeline": "metaphlan", "bundle_id": "gut-reads-adult-global-v0.2-lenient"}}`; `checkgm normalize` stamps the bundle id onto a query directory's whole `sample.json` entry instead, so anything beyond the two fields is there because that file carried it |

The table covers the files a reader will open, and the columns of `scores.parquet`, `summaries.tsv` and `rejections.tsv`
are given in full in `docs/cli.md`.

How to read it: a healthy adult sample measured like the reference has about 5 % of its features outside the 2.5–97.5
percentile band. The held-out healthy cohorts used for validation sit between 2 % and 13 % on the taxonomy and gene
layers, most of them between 3 % and 8 %, and far lower on the module layer, where the three cohorts run from 0.2 % to
0.9 % (`docs/validation.md`, section 2). The sample-level excess test (`excess_outside_q ≤ 0.05`) is the
primary statement. Percentiles are the robust quantity; individual low/high calls are descriptive and less reproducible
between re-measurements of the same reads than percentiles are.

## What has been validated, and what has not

The numbers below come from the validation runs that `docs/validation.md` describes report by report, each run's output
being copied verbatim into `docs/reports/`, and they agree with the manuscript in `paper/`. The four read-based reports
(`reads_calibration_and_evaluation.md`, `reads_disease_atlas.md`, `reads_healthy_strata.md`,
`reads_disease_by_stratum.md`) were re-copied on 2026-10-06 from the run of that date, which supersedes every earlier
read-based figure quoted anywhere in this repository, including the 0.690 at p = 0.052 of 2026-09-30 and the 0.693 at
p = 0.039 of earlier on 2026-10-06. That run differs from the one before it for one reason: the out-of-sample scoring
was changed so that no study is split across two baselines and so that two curatedMetagenomicData 3 deposits of one
donor cohort are held out together (`docs/validation.md` section 2).

- Re-measuring public samples from raw reads with the assembly-based pipeline container and scoring them reproduces
  scores obtained from the original MGnify analyses (61 samples; percentile Spearman 0.98–0.99 on gene layers, 0.86–0.88
  on taxonomy). The read-based pipeline reproduces the curatedMetagenomicData 3 profiles from raw reads (34 samples,
  median species Bray–Curtis 0.011, score percentile Spearman ≥ 0.997). The read-based function layers (HUMAnN 3) have
  **not** been validated from raw reads.
- Held-out healthy cohorts mostly stay near the expected 5 % outside the band, but not all do, one assembly cohort
  reaching 13 % on the gene layers. Cohorts from Africa do not fit the assembly-based baseline at all: a median of
  19.8 % of families per sample falls outside the band across 372 samples from two African studies
  (`docs/reports/assembly_healthy_strata.md`, the region × taxonomy_family row). Location is confounded with study
  throughout.
- Reference-relative features separate cases from controls only modestly better than raw abundances. With a random
  forest trained leave-one-study-out, the assembly-based pipeline reaches an AUROC of 0.618 against 0.581 for raw
  abundances (9 of 11 studies, p = 0.021 uncorrected), and 0.615 against 0.581 on the strict baseline (p = 0.051). The
  read-based taxonomy reaches 0.688 against 0.664 over 14 studies (10 of them, p = 0.034), and 0.770 against 0.715 in
  colorectal cancer alone (8 of 9, p = 0.004). That advantage is specific to the classifier and nothing more is claimed
  for it, because a matched L1-penalised model does not reproduce it (0.720 against 0.705, two-sided p = 0.30); what the
  percentiles buy is an output that can be read, at about half the features for the same accuracy (a median of 104.5
  non-zero coefficients against 204). Nothing is added by the function layers on top of the taxa, so their place is
  earned by what they make readable rather than by what they predict.
  They are not weak in themselves: on the assembly-based pipeline, whose taxonomy stops at genus, the gene layers alone
  are the stronger ones, reaching 0.631 against 0.585 for the taxa. Of 43 known-biology expectations tested on that
  pipeline, 14 were confirmed.
- The published comparators are matched rather than beaten. Over the same 14 read-based studies GMHI reaches 0.62 and
  GMWI2 0.72, against 0.75 for the manuscript's supervised health score. That score is an analysis reported in `paper/`
  and not a feature of this package: no subcommand fits or applies one. Neither comparison is made on wholly equal
  terms. Only 33 of GMHI's 50 species are matched by name in MetaPhlAn 3, so its figure should be read as a lower bound
  on what that index achieves on its native input. GMWI2's training set, meanwhile, contains samples from 10 of those
  14 studies, whereas every checkGM figure here comes from a model trained without the study it scored. On the four
  studies neither model has seen the two are level, at 0.62 against 0.61, which is too few studies to separate them and
  no separation is claimed.
- The baselines are adult, stool, global. Other body sites and infants are refused (`NO_MATCHING_REFERENCE`).

## Documentation

- `docs/cli.md` — every subcommand, its arguments, outputs and exit status
- `docs/query_format.md` — the per-sample query directory any assembly pipeline can produce
- `docs/bundle_format.md` — manifest fields, `features/*.parquet` columns, `normalization/`
- `docs/rebuild.md` — how a baseline is rebuilt
- `docs/validation.md` — how validation was done, with the report files in `docs/reports/`
- `docs/configs/` — read-only copies of the baseline-build and evaluation settings (thresholds, seeds)
- `paper/` — a snapshot of the manuscript and the five figures it includes (`main.tex`, IEEEtran conference class;
  `tectonic main.tex` builds it; the project workspace holds the draft under revision)
- `scripts/viz2/` — the scripts that draw those five figures: `prep.py` writes one tidy CSV per panel and `figures.R`
  draws from it, so no number in a figure is typed by hand (see its README)
- `scripts/viz/` — the superseded Python set behind the earlier drafts, sixteen figures numbered up to 18, kept for
  provenance and drawn on by no figure of the current manuscript
- `CHANGELOG.md`

## Development notes

- The early example profile `checkgm/examples/CosteaPI_2017_alien2-11-0-0.txt` (a public curatedMetagenomicData 3 profile,
  commits `97cba6b`, `92c85ce`) stays in the history of branch `checkgm` by the owner's decision of 2026-10-02 (`CHANGELOG.md`);
  the example shipped since 0.8.0 is synthetic.
- Tests: `pip install -e ".[test]" && pytest` (155 tests, about 100 s). CI runs them on Python 3.10 and 3.12
  (`.github/workflows/test.yml` at the repository root).
- `scripts/viz2/` and `scripts/viz/` are dated snapshots of the project workspace's figure scripts and are not edited
  here, the workspace copy being the authority, so neither should be expected to diff clean against it. The v2 pair
  was taken on 2026-10-06 and the workspace has moved on within the same day; the older `scripts/viz/` snapshot has
  drifted further, 14 of its 21 scripts now differing, and is best read as a record of how the superseded figures were
  drawn.

## Citation

Khoa Hoang, Stanford University. The paper is not yet published; see `CITATION.cff` for the intended venue, and cite
this repository meanwhile. The baselines and container image carry the reserved DOI `10.5281/zenodo.23181418`.

## Licence

MIT (see `LICENSE`). `container/emapper2/` is eggNOG-mapper 2.0.0 (GPL v2) as vendored by MGnify pipeline v5 and is
used only inside the assembly-based pipeline container; `container/pv5_ref/` holds reference copies of MGnify
pipeline-v5 workflow files from `EBI-Metagenomics/pipeline-v5` (Apache License 2.0), kept for provenance and executed
by nothing here.
