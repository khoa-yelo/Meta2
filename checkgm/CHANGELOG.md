# Changelog

## Unreleased

- 2026-10-06, **the command-line interface is now five commands**: `setup`, `profile`, `assess`, `score` and
  `end-to-end`. `assess` determines the kind of input from what the path holds (an assembly query directory, an MGnify
  v5 analysis, a MetaPhlAn profile, or a directory `assess` itself wrote), so one command serves every entry point, and
  it draws the per-sample report figure as well as the tables. `profile` runs a measurement pipeline and `setup`
  fetches the baselines or a pipeline's databases; both check their prerequisites and refuse in one line before
  starting a job that would otherwise fail after a 90 GB download. One rename affects existing command lines: `score`
  now means the supervised health score, computed from an `assess` output, and what `score` used to do is `assess` or
  the unchanged stage command `score-normalized`. The old flags given to `score` print that redirect. The eight stage
  commands remain registered and unchanged, since the paper's methods name them.
  - The shipped health-score model is the percentile-only fit (104 features, no species-presence term and so none of
    the batch artifact the paper describes). It is exported complete — intercept, and per feature its coefficient and
    the standardizer's centre and scale — because a coefficient table alone cannot score a sample; the previous export
    could not have been applied by anyone. An assembly assessment is refused rather than scored, and a sample carrying
    fewer features than the cohort's 1st percentile (12 of 104) is flagged as uninformative.
  - `matplotlib` is an optional extra (`pip install "checkgm[figures]"`). Without it `assess` writes the tables and
    says the figure was skipped, since the percentiles are the product and the figure is one rendering of them.
  - Fixed, found while wiring the above: an empty call direction in a report panel drew the three worst features of
    the other direction a second time (`DataFrame.assign` of a full-length Series onto an empty frame adopts the
    Series' index instead of aligning to nothing); the panel ranked its extremes by percentile, which saturates at
    0.5/99.5 so every badly deviant feature tied; the print profile was written into the process-global `rcParams`,
    restyling later plots of any program embedding checkgm; an absent `matplotlib` raised `ModuleNotFoundError` past
    the one-line error contract; and an assessment that rejected every sample crashed the figure writer.
  - `README.md` corrected: the assembly databases are about 90 GB to download and 170 GB on disk, not the 70/110 GB
    recorded before eggNOG and InterProScan were fully provisioned (the figures now reconcile to the byte against the
    staged copy, excluding the hg38 index the pipeline does not use by default), and the note that the scoring image
    still carried an entrypoint under the tool's former name was stale — it does not.

- 2026-10-06, **the read-based out-of-sample scoring was wrong in two ways and has been corrected; every read-based
  number in this repository and in `paper/` comes from the re-run.** Both defects were in how the project decides which
  baseline scores a sample (`project/scripts/s12_oos_scores.py`), and both were found by auditing the one condition the
  manuscript rests on: that no sample was ever compared with a baseline that contained it.
  - **A study was split across two baselines.** The exclusion of a sample that also sits in the pool under another
    study's name was applied sample by sample. In GuptaA_2019 all 30 controls are DhakanDB_2019 duplicates and none of
    the 30 colorectal cancer cases is, so the controls were scored against the baseline rebuilt without DhakanDB_2019
    and the cases against the full baseline. The one baseline difference inside that study was therefore perfectly
    confounded with the case/control label — the exact artifact the manuscript's leakage paragraph rules out — and the
    study's reported shifts were differences between two baselines rather than percentiles of one. The exclusion is now
    a property of the study: if any sample of a study duplicates a pool study, the whole study is scored against the
    baseline rebuilt without it. Checked sample by sample, no study is now split (one was).
  - **A donor cohort deposited twice was only half held out.** The duplicate audit matched ENA run-accession strings,
    and LeChatelierE_2013 and NielsenHB_2014 are two curatedMetagenomicData 3 deposits of one MetaHIT Danish cohort
    under largely disjoint accession sets: 177 specimens appear in both, of which only MH0001–MH0005 share any
    accession. The audit recorded 5 of 171 pool-internal pairs, so when one deposit was left out for the
    leave-one-study-out calibration the other deposit's profile of the same specimen stayed in the baseline that scored
    it — for 342 of the 6,494 pool samples. The audit is now a committed script,
    `project/scripts/s11_cross_study_duplicates.py`, which runs an accession test and a sequencing-fingerprint test
    (identical `number_reads` and `number_bases`) and unions them, because neither is a superset of the other; it finds
    201 pairs where the old table held 35. A leave-one-study-out baseline now holds out every study sharing a specimen
    with the one being scored, so the two deposits are held out together. Checked sample by sample, 0 of the 9,123
    scored samples now sees itself or another deposit of itself in its baseline (342 did).
  - **What moved.** The held-out-cohort calibration is unchanged to the last digit, no sample of those cohorts being
    duplicated anywhere. In-pool leave-one-study-out figures moved for the two MetaHIT deposits and for no other study;
    of the per-layer summaries in `reads_calibration_and_evaluation.md` (median over studies of the study median) only
    genus changed, 4.7 % to 4.8 %, and of the median-over-studies of the study *mean* that the manuscript's calibration
    panel plots, likewise only genus, 5.53 % to 5.63 %. The read-based taxonomy case/control comparison is now
    0.688 against 0.664 over 14 studies (10 of them, one-sided p = 0.034, from 0.693 and p = 0.039), and 0.770 against
    0.715 in colorectal cancer alone (8 of 9, p = 0.004). GMHI (0.615), GMWI2 (0.725) and Shannon (0.476) are
    unchanged, being computed from raw abundances. The pre-specified health score is 0.747 as before. The matched
    sparse-model control moved most: percentiles against raw abundances with summaries and presence added is now higher
    in 7 of 14 studies at two-sided p = 0.81, where it was 10 of 14 at p = 0.30, which strengthens rather than weakens
    the manuscript's claim that the random forest's advantage is classifier-specific. In the colorectal cancer atlas
    the consistent-feature count went from 198 to 199. The six families and ten genera the manuscript names are
    unchanged in membership and direction, and in median shift but for one two-cohort genus that moved 0.2 points; the
    function layers moved, to 99 MetaCyc pathways (86 down) and 82 KEGG orthologs (72 down) from 100 (85 down) and 80
    (70 down). Only GuptaA_2019's own feature rows changed at all, by a mean 1.1 and at most 7.3 points of shift.
  - The four read-based validation reports in `docs/reports/` were re-copied from the new run, which also closes the
    open item recorded below: the repository no longer carries two different verdicts on the taxonomy comparison.
  - Not fixed, because it is a composition fact rather than a leakage one and correcting it would mean rebuilding the
    released bundle: the read-based pool's 6,494 samples are 6,323 distinct specimens, so 171 donors carry twice the
    weight of the rest. This is now stated in the README and in `docs/validation.md` section 2.

- 2026-10-06, `checkgm run-normalize` added: the one-step form of `normalize` followed by `score`, which the
  `run-mgnify` and `run-metaphlan` pairs already had and the generic assembly path did not. A test asserts that
  its `scores.parquet` is identical to the two-step result, so the shortcut cannot drift from the path it
  replaces. The README and the tier-A container recipe now point assembly users at it.

- 2026-10-06, documentation review round 6. Four claims that round 5 either introduced or strengthened past what the
  repository's own files support have been withdrawn. The README no longer asserts that `scripts/viz2/` matches the
  project workspace exactly, an equality that could not hold of a tree still being edited and did not hold when
  checked: both scripts differ, the repository's copies dating from 00:15 and the workspace's from 01:02 and 01:09 of
  the same day. Both snapshots are therefore described by their date alone. The `sample_meta.json` row no longer
  promises a sample id, a body site and a pipeline version, since only `checkgm normalize` copies a query's
  `sample.json` wholesale, while the import and `run-*` commands write the pipeline name by itself; the row now says
  which is which. "The function layers add no predictive advantage at all" is replaced by the narrower claim the
  manuscript defends, that nothing is added on top of the taxa, because on the assembly-based pipeline the gene layers
  alone are in fact the stronger ones at 0.631 against 0.585 (`docs/reports/assembly_case_control_evaluation.md`).
  `scripts/viz2/README.md` gains the two R dependencies that no `library()` call names, `ragg` and a cairo-enabled R,
  both reached through `pkg::fn` inside the `save_fig` helper and both required for a fresh install to get through the
  first figure.
- Corrected with them, the README having been less careful about a competitor than the paper is:
  - GMHI now carries the caveat the manuscript gives it, that only 33 of its 50 species are matched by name in
    MetaPhlAn 3, so its 0.62 is a lower bound on what the index achieves on its native input;
  - the 0.75 beside it is marked as a manuscript analysis rather than a shipped feature, nothing in the package
    fitting or applying a health score;
  - the superseded figure set is sixteen figures numbered up to 18, not eighteen figures, in both places that counted
    it;
  - the unfinished rename of `REFMB_PROJECT` to `CHECKGM_PROJECT` is documented as unfinished, the shipped
    `tokens.py` reading the new name while `prep.py` exports the old one, and both are now set in the figure command;
  - the scoring image's leftover `REFMB_BUNDLES` default is noted beside the entrypoint caveat it belongs with.
- ~~**Open item.**~~ *Resolved by the entry at the top of this file on 2026-10-06.* `docs/reports/reads_calibration_and_evaluation.md`
  was still the run of 2026-09-30, and section 4 of `docs/validation.md` still quoted it, so the repository carried the
  superseded 0.690 at p = 0.052 beside the README's 0.693 at p = 0.039. Both now carry the run of 2026-10-06 16:12 UTC,
  0.688 at p = 0.034 and 0.770, and that run supersedes the 0.693 as well.
- 2026-10-06, documentation review round 5. The read-based case/control figures in the README now follow the project
  workspace's `results/s11/pipelineB_report.md` as regenerated on 2026-10-06, which gives taxonomy 0.693 against 0.664
  over 14 studies (10 wins, p = 0.039) and 0.771 against 0.715 in colorectal cancer alone (p = 0.004). The earlier
  run's 0.690 at p = 0.052 is superseded, and because 0.052 and 0.039 fall on opposite sides of 0.05 the repository
  and the manuscript had been giving a reader two different verdicts on one comparison. Stated beside them is the
  caveat the paper makes, that a matched L1-penalised model does not reproduce the margin (0.719 against 0.705,
  p = 0.33), so the random-forest result is not read as a general one. GMHI and GMWI2 are now named in the README with
  the training-overlap caveat that makes the GMWI2 comparison an indirect one.
- Corrected with them, each against the file it has to agree with:
  - the test count (31, about 50 s, where 0.8.2 recorded 30 and 40 s);
  - `pip install -e ".[test]"`, quoted as the CI workflow quotes it, because zsh reads the bare form as a glob;
  - the genus CLR basis in `docs/query_format.md` (1,687 ids, as a bundle's `normalization/clr_basis_genus.txt` and
    its manifest both have it, not 1,690);
  - the version banner's capitalisation, with one sentence settling the convention that checkGM names the project and
    `checkgm` the package, the command and the image tag;
  - the scoring image's build command, whose `checkGM` tag podman and docker both reject;
  - the strict tier's description, which omitted the health-label requirement that distinguishes it;
  - a `sample_meta.json` row in the output table;
  - the module layer's far lower fraction outside the band;
  - a citation for the 19.8 % of families outside in the African cohorts, which sits in
    `assembly_healthy_strata.md` rather than in either report the text had pointed at.
- `docs/rebuild.md` no longer implies that `checkgm.build_reference.build` produced the released assembly baseline. It
  did not: those manifests record `refmb 0.1.0 (project/scripts/s5_build_reference.py)`, a standalone script that
  never imports the package, while the packaged builder wrote the read-based bundles only. Whether the two agree value
  for value is untested, and the page now says so.
- The documentation map gains `paper/` with its build line, and names `scripts/viz2/` as the set that draws the
  manuscript's five figures, `scripts/viz/` being the superseded set that no current `\includegraphics` reaches. The
  claim that `scripts/viz/` holds byte-identical copies of the workspace scripts is withdrawn, here and for 0.8.2,
  having held at neither release: 14 of its 21 scripts now differ, and both directories are described as dated
  snapshots instead. `scripts/viz2/README.md` also drops the suggestion that `CHECKGM_PROJECT` alone makes the scripts
  portable, since the project root is hard-coded in both of them.
- The *Get a reference bundle* section says plainly that no public download exists yet and that the tarballs have to
  be requested from the authors, and the shell block carries that missing step as its first line. Verifying and
  unpacking a file nothing had fetched was the largest usability gap in the page, because every command after it
  depends on the bundle being there.
- 2026-10-02: the repository owner decided to leave the early example profile `checkgm/examples/CosteaPI_2017_alien2-11-0-0.txt`
  (a public curatedMetagenomicData 3 MetaPhlAn 3 profile, commits `97cba6b` and `92c85ce`) in the pushed history of branch
  `checkgm` rather than rewrite it; the example profile shipped since 0.8.0, `examples/synthetic_healthy_adult.txt`, is synthetic.

## 0.8.2 — 2026-10-02, packaging round 4

Documentation and build hygiene after the third external review; percentiles, calls, the excess test and the landscape
columns are unchanged (`scores.parquet` and `summaries.tsv` are identical to 0.8.1 on the packaged test inputs).

- The report header names the measurement plainly: `6,494 healthy adult gut read-based profiles (curatedMetagenomicData 3:
  MetaPhlAn 3 + HUMAnN 3)` and `1,941 healthy adult gut assembly analyses (MGnify pipeline v5.0)` for the released bundle
  series; other bundles keep the manifest's `source_pipeline` name and version. The fraction-outside table prints three
  decimals (`0.000`, as the README says) instead of dropping trailing zeros.
- `docs/cli.md` documents all 31 columns of `summaries.tsv` (added: `expected_false_positives_at_alpha`, `layer_status`,
  `mapped_fraction_family`, `genome_equivalents_cov`, `fallback_used`, `landscape_centring`, `neighbor_bands`).
- `paths.py` drops `bundle_A`/`bundle_B`, which nothing in the package,
  tests or containers called; the figure scripts import them from the project workspace's own copy of `paths.py`.
- `container/Dockerfile` no longer copies `tests/` into the scoring image (they were never run there);
  `container/checkgm_tierA.def` installs `python=3.8` once instead of twice. The assembly image does not install the checkGM
  package (its Python 3.8 predates checkGM's requirement of 3.10), so `pv5_to_query.py` records the checkGM version in
  `sample.json` only where checkGM is importable and writes `unknown` otherwise; 0.8.1's note is corrected below and in
  `docs/query_format.md`.
- README: 30 tests (about 40 s); the five-minute example says that the bundle must already be unpacked under
  `$REFMB_BUNDLES` and that the baselines and the container image will be deposited on Zenodo (DOI [TBD]); a *Development
  notes* section records the owner's decision above about the early example profile.
- `docs/reports/assembly_all_gut_analyses_scoring.md` is the verbatim copy of the survey re-scored with the 0.8.1 scorer
  (2026-10-02); only its landscape table (distance-from-centre percentiles by health label) differs from the 0.8.0 copy.
- `docs/validation.md` and `docs/configs/README.md` define the internal labels that the verbatim report copies and
  configuration files carry (S6, S7, S8, S11, S12/S14, anchors, tier A/B).
- The default branch's README (`main`) now points to `checkgm/` on branch `checkgm`. `scripts/viz/` re-copied byte-identical
  from the project workspace after the fourth figure revision (10 of 17 scripts changed; see the project's figure notes).
  The byte-identical half of that claim is withdrawn under *Unreleased*, not having held at this release either.

## 0.8.1 — 2026-10-02, packaging round 3

Bug fixes after the second external review. Percentiles, calls and the excess test are unchanged (`scores.parquet` is
identical to 0.8.0 on the packaged test inputs); the healthy-map columns change for every sample.

- **Healthy map fixed.** `score.py` projected a query's family CLR with the multiplicative-replacement values of undetected
  families, while the baseline's principal components were fitted on a matrix in which undetected families are missing
  and filled with the pool minimum. Every sample was displaced along PC1, and healthy pool samples re-scored at a median
  core-distance percentile of about 90 instead of 50. The scorer now masks undetected basis features before the fill.
  Re-projecting pool samples through `checkgm.score` now reproduces the bundle's `landscape/reference_coords.parquet`
  (40 samples of the released assembly bundle: largest difference 5.6e-7, from the stored column means being rounded to six
  decimals; synthetic bundles: exact). `core_distance_pct`, `PC1`-`PC3`, `neighbor_studies` and `neighbor_bands` in
  `summaries.tsv` and the report's landscape table are therefore different from 0.8.0 for every sample; 0.8.0 values
  should not be used. `build_reference` now stores column means and the centroid at full precision and records the fill
  convention in `core_distance.json` (`tests/test_landscape.py`).
- `checkgm score` after `checkgm normalize` failed with a traceback when the query's `sample.json` carried a `pipeline`
  object (`{"name": ..., "version": ...}`, as the container writes it). Found by the new assembly-path test.
- `--calibration` is checked (directory, `calibration.json`, `features.parquet`) and loaded before the import step, so a
  wrong path fails before anything is written to `--out`.
- A bundle from a superseded series (read-based v0.1, assembly-based v0.7 and below) prints one stderr line
  `checkgm: bundle series v0.1 is superseded by v0.2 ...` and is still scored (exit 0).
- The report's landscape table drops columns that are entirely NaN (`genome_equivalents_cov` for read-based bundles).
- Rejected read-based profiles no longer emit a numpy `RuntimeWarning` from the CLR.
- `checkgm.calibrate.load_config()` requires a path (ValueError naming `docs/configs/s11.yaml`); the old default pointed
  into the development tree.
- Tests: 20 -> 30 (`tests/test_landscape.py`, `tests/test_assembly_path.py`: the documented query directory through
  `checkgm normalize` and `checkgm score` against a synthetic assembly bundle built with `checkgm.build_reference`;
  `tests/test_cli_polish.py`).
- Continuous integration: the workflow moved from `checkgm/.github/workflows/` to `.github/workflows/test.yml` at the
  repository root, where GitHub reads it (0.8.0 claimed CI runs before any run existed). The first run (commit 7fbc355,
  run 37029618083) passed on Python 3.12 and failed on 3.10 on one test's assumption about argparse's help layout, fixed
  here; the suite then passes locally on CPython 3.10.22 (pandas 2.3.3, numpy 2.2.6, pyarrow 25.0.1) and 3.12.14
  (pandas 3.0.6). The run on this commit is the first green one, if the Actions page agrees.
- Wording: internal stage labels removed from code comments and shell headers (`fetch_dbs.sh`, `unpack_dbs.sh`,
  `normalize.py`, `build_reference.py`, `score.py`, `calibrate.py`, which now defines anchor samples as samples measured by
  both pipelines); `run_tierA.sh` header says metaSPAdes 3.15.3; `pv5_to_query.py` writes pipeline name `checkgm-assembly`,
  the checkgm version when the package is importable (it is not inside the assembly image, whose Python 3.8 cannot import
  checkGM, so the field reads `unknown` there; see 0.8.2) and the assembler actually used into `sample.json`; one-line docstrings on the public names
  that lacked one. `docs/configs/evaluation_settings.md` states the inclusion rule as applied per pipeline (assembly-based
  at least 20 cases and 10 controls, 11 studies; read-based at least 20 and 20, 14 studies) and that `HMP_2019_ibdmdb` is a
  largely adolescent cohort (median age 16.5) included in the read-pipeline evaluation; `docs/validation.md` gives the
  package versions the timing was measured with (pandas 3.0.6, numpy 2.5.3, pyarrow 25.0.0); README and `docs/cli.md`
  say that `reason_code` is `NOT_NORMALIZED` or `NO_MATCHING_REFERENCE` and that floor names appear in `detail`.
- **Open item for the repository owner.** The public cMD3 per-sample profile `checkgm/examples/CosteaPI_2017_alien2-11-0-0.txt`
  that 0.8.0 removed from the tree is still in the pushed history of branch `checkgm` (commits `97cba6b` and `92c85ce`,
  before `7ceaba7`). Removing it requires rewriting and force-pushing the branch, which needs the owner's decision; the
  exact command is prepared in the project workspace (`distribution/HISTORY_CLEANUP.md`, a workspace file that is not
  shipped in this repository). Until then the history holds one
  public individual's MetaPhlAn profile.

## 0.8.0 — 2026-10-02, packaging round 2 (version number unchanged)

Packaging fixes after the first external review; scoring is unchanged (`scores.parquet` byte-identical on the packaged
test inputs for both pipelines).

- Command line: user errors (bundle not found, bundle of the wrong kind for the command, missing input file or directory,
  a file that is not a MetaPhlAn table) exit 1 with one line `checkgm: ...` on stderr instead of a traceback. Every
  subcommand and argument has help text (`checkgm <command> --help`). The bundle's `profile_type` is checked against the
  subcommand (`run-metaphlan`/`import-metaphlan` need a read-based bundle, `normalize`/`import-mgnify`/`run-mgnify` an
  assembly-based one, `score` checks against the normalized directory). A profile without species rows is rejected as
  `EMPTY_PROFILE` in `rejections.tsv`; `--reads-tsv` accepts a header line and reports unparseable rows.
- `report.md` prints taxon names next to NCBI taxids and lists expected-but-missing features by name; the
  reproducibility sentence says "re-profiling" for read-based bundles. `qc.tsv` stores floats with 17 significant digits
  and is read back exactly, so `import-metaphlan` + `score` now equals `run-metaphlan` byte for byte.
- Example: `examples/CosteaPI_2017_alien2-11-0-0.txt`, a public individual's MetaPhlAn profile, is replaced by
  `examples/synthetic_healthy_adult.txt`, generated from the read-based bundle's percentile medians. No per-sample
  public profile is redistributed in the repository.
- Repository: figure scripts (`scripts/viz/`), evaluation and baseline-build configuration files (`docs/configs/`) and
  two more validation reports (`docs/reports/assembly_all_gut_analyses_scoring.md`; `assembly_known_biology.md` was
  already present) added; baseline pool accession lists (public identifiers only) are distributed beside the bundles,
  not on GitHub (README download table). `CITATION.cff` validates against the 1.2.0 schema (`date-released` and `doi`
  omitted until known). GitHub Actions runs the tests on Python 3.10 and 3.12. The sdist now includes `docs/` and
  `examples/`. Wording: "assembly-based" and "read-based" pipelines throughout; "standard baseline (bundle name lenient)"
  and "strict baseline". `fetch_dbs_B.sh` pins the md5 of the marker database archive. Docstrings on the public
  functions; unused names removed. Timing of the README example recorded in `docs/validation.md`.

## 0.8.0 — 2026-10-02

First public packaging of the tool. The version number follows the current assembly-based bundle series (v0.8); the
read-based bundle series is at v0.2.

- Installable package (`pip install .`) with the console script `checkgm` (subcommands `normalize`, `import-mgnify`, `score`,
  `run-mgnify`, `import-metaphlan`, `run-metaphlan`, `bundles`).
- Bundles are given as a directory path or as a bundle id resolved under `$REFMB_BUNDLES`; the development-machine paths
  that `checkgm/paths.py` used to hold are gone. Scoring, normalization and import code are unchanged numerically (verified
  on the packaged test inputs: identical `scores.parquet` against the pre-packaging code for both pipelines).
- `checkgm --version` and `checkgm bundles` added.
- Container recipes: `container/Dockerfile` (scoring only), `container/checkgm_tierA.def` (the assembly-based pipeline;
  the file name is kept so the released image checksum stays valid), `container/checkgm_pipelineB.def` (the read-based
  pipeline, MetaPhlAn 3.0.14 + bowtie2 + checkGM). `run_pipelineB.sh` now calls the
  installed `checkgm` command; `run_tierA.sh` lost its development-host fallback paths (set `SPADES_BIN`, `EMAPPER_PY`,
  `JAVA8_BIN` outside the container).
- Documentation in `docs/`, validation reports copied to `docs/reports/`, a worked example in `examples/`, and a small
  pytest suite that scores a MetaPhlAn profile against a synthetic bundle.
- Not packaged: the project-internal script `bundle_add_normalization.py` (it rewrote development bundles in place and
  depended on project configuration files). Bundles already carry their `normalization/` directory.

### Reference bundles released alongside (distributed separately, see README)

- gut-assembly-adult-global-v0.8-{lenient,strict}: 1,941 healthy adults, 26 studies, MGnify v5 assembly measurement.
- gut-reads-adult-global-v0.2-{lenient,strict}: 6,494 healthy adults, 19 studies, MetaPhlAn 3 + HUMAnN 3 measurement.
  The v0.1 read-based bundles are obsolete and must not be used.
