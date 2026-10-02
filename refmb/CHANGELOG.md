# Changelog

## 0.8.0 — 2026-10-02, packaging round 2 (version number unchanged)

Packaging fixes after the first external review; scoring is unchanged (`scores.parquet` byte-identical on the packaged
test inputs for both pipelines).

- Command line: user errors (bundle not found, bundle of the wrong kind for the command, missing input file or directory,
  a file that is not a MetaPhlAn table) exit 1 with one line `refmb: ...` on stderr instead of a traceback. Every
  subcommand and argument has help text (`refmb <command> --help`). The bundle's `profile_type` is checked against the
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

- Installable package (`pip install .`) with the console script `refmb` (subcommands `normalize`, `import-mgnify`, `score`,
  `run-mgnify`, `import-metaphlan`, `run-metaphlan`, `bundles`).
- Bundles are given as a directory path or as a bundle id resolved under `$REFMB_BUNDLES`; the development-machine paths
  that `refmb/paths.py` used to hold are gone. Scoring, normalization and import code are unchanged numerically (verified
  on the packaged test inputs: identical `scores.parquet` against the pre-packaging code for both pipelines).
- `refmb --version` and `refmb bundles` added.
- Container recipes: `container/Dockerfile` (scoring only), `container/refmb_tierA.def` (the assembly-based pipeline;
  the file name is kept so the released image checksum stays valid), `container/refmb_pipelineB.def` (the read-based
  pipeline, MetaPhlAn 3.0.14 + bowtie2 + refmb). `run_pipelineB.sh` now calls the
  installed `refmb` command; `run_tierA.sh` lost its development-host fallback paths (set `SPADES_BIN`, `EMAPPER_PY`,
  `JAVA8_BIN` outside the container).
- Documentation in `docs/`, validation reports copied to `docs/reports/`, a worked example in `examples/`, and a small
  pytest suite that scores a MetaPhlAn profile against a synthetic bundle.
- Not packaged: the project-internal script `bundle_add_normalization.py` (it rewrote development bundles in place and
  depended on project configuration files). Bundles already carry their `normalization/` directory.

### Reference bundles released alongside (distributed separately, see README)

- gut-assembly-adult-global-v0.8-{lenient,strict}: 1,941 healthy adults, 26 studies, MGnify v5 assembly measurement.
- gut-reads-adult-global-v0.2-{lenient,strict}: 6,494 healthy adults, 19 studies, MetaPhlAn 3 + HUMAnN 3 measurement.
  The v0.1 read-based bundles are obsolete and must not be used.
