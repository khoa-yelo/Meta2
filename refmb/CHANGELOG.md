# Changelog

## 0.8.0 — 2026-10-02

First public packaging of the tool. The version number follows the current assembly-based bundle series (v0.8); the
read-based bundle series is at v0.2.

- Installable package (`pip install .`) with the console script `refmb` (subcommands `normalize`, `import-mgnify`, `score`,
  `run-mgnify`, `import-metaphlan`, `run-metaphlan`, `bundles`).
- Bundles are given as a directory path or as a bundle id resolved under `$REFMB_BUNDLES`; the development-machine paths
  that `refmb/paths.py` used to hold are gone. Scoring, normalization and import code are unchanged numerically (verified
  on the packaged test inputs: identical `scores.parquet` against the pre-packaging code for both pipelines).
- `refmb --version` and `refmb bundles` added.
- Container recipes: `container/Dockerfile` (scoring only), `container/refmb_tierA.def` (pipeline A, assembly-based),
  `container/refmb_pipelineB.def` (pipeline B, MetaPhlAn 3.0.14 + bowtie2 + refmb). `run_pipelineB.sh` now calls the
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
