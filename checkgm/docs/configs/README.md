# Evaluation and baseline-build settings

Read-only copies of the configuration files the released baselines and their evaluation were run with, so that every
threshold and seed is inspectable. The scripts that consume them live in the project workspace (see `docs/rebuild.md`);
nothing in the installable package reads these files. The copies are verbatim and keep the workspace's internal labels:
S5 = the baseline build, S6 = the held-out calibration check, S7 = the excess-test survey over 8,742 MGnify gut analyses,
S8 = the case/control evaluation, S11 = the reproduction from raw reads, S12 and S14 = the healthy strata and the disease
atlas, "anchors" = samples measured by both pipelines, and "tier A" / "tier B" = the assembly-based / read-based pipeline.

| file | used by | what it fixes |
|---|---|---|
| `s5.yaml` | baseline build (assembly-based) | seed 20261020, bundle version and body site, pool exclusion roles and the 3-samples-per-study minimum, equal-per-study weighting, the percentile grid, study bootstrap (500), minimum detections for percentiles (20), the six layers and their value columns, the landscape PCA settings, the detection floor (0.05 copies per genome) and which layers take percentiles within the assembly-size class |
| `inclusion.yaml` | pool selection | what "adult" (stratum) means, and the standard (`lenient`) and `strict` definitions of a healthy adult (health label, antibiotics, age 18-65, BMI 18.5-30, pregnancy) |
| `splits_rules.yaml` | study-level roles (assembly-based) | the pinned held-out healthy studies (PRJEB62833, PRJEB75578, PRJEB39906), studies needing review, and the case/control fold rule (at least 20 cases and 10 controls) |
| `s11.yaml` | `checkgm.calibrate` | thresholds of the calibration gate, fixed before any result was seen; no calibration is shipped, the file documents the method |
| `evaluation_settings.md` | case/control evaluation | the settings of the cross-study classification that are written in the evaluation script itself rather than in a YAML file |

The read-based baseline (`gut-reads-adult-global-v0.2`) carries its own rules inside the bundle
(`normalization/rules.json`: depth bands, CLR delta, mapped-fraction floor, held-out studies, seed) and in its
`manifest.json`; `docs/bundle_format.md` describes them.
