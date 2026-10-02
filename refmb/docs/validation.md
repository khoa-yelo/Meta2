# How the baselines were validated

The reports in `docs/reports/` are verbatim copies of the result files the numbers in the paper and in the README come from.
This page says what each one tested and where its headline numbers are.

## 1. Does the pipeline reproduce the reference measurement from raw reads?

A baseline is only usable if a new sample can be measured the same way the pool was.

- **Assembly-based pipeline** — `reports/assembly_fidelity_from_raw_reads.md`. 61 public samples that exist both as MGnify v5
  analyses and as deposited reads were re-assembled and re-annotated with the container, then scored against the same
  bundle as the original analyses. Per layer the report gives the fraction outside the band for container and native
  measurement, the Spearman correlation of percentiles, the Jaccard overlap of detected features and the call agreement.
  Gene layers (KO, Pfam) reproduce to percentile Spearman 0.98–0.99 and call agreement 0.99; family and genus taxonomy to
  0.86–0.88 and 0.97–0.98 (taxonomy is the layer most sensitive to re-assembly). All layers pass the gate.
- **Read-based pipeline** — `reports/reads_fidelity_from_raw_reads.md`. 34 samples (17 studies) were profiled from deposited
  reads with MetaPhlAn 3.0.14 and compared with the curatedMetagenomicData 3 profile: median species Bray–Curtis 0.011,
  Spearman 0.999; scores against the baseline agree with percentile Spearman ≥ 0.997 and identical calls; 30 of 34 land in the
  same depth band. One further sample was set aside because the deposited run has less than half the reads cMD3 reports.
  The HUMAnN 3 function layers were **not** re-run from raw reads; their fidelity is unverified.

## 2. Do healthy samples that were never in the pool fit the baseline?

Expectation: about 5 % of features outside the 2.5–97.5 percentile band; the gate used during development was 8 % per study.

- **Assembly** — `reports/assembly_heldout_calibration_check.md` (three held-out healthy studies, 233 samples; per layer and
  study the mean fraction outside and whether it passes) and `reports/assembly_loso_within_pool.md` (each pool study scored
  against a reference rebuilt without it: medians 0.05–0.055 on all layers except module; 73–85 % of studies within 8 %).
  Cohorts from Africa do not fit (19.8 % of families outside), and location cannot be separated from study.
- **Reads** — `reports/reads_calibration_and_evaluation.md` (three held-out healthy studies and in-pool leave-one-study-out,
  all five layers).

## 3. Does the healthy baseline move with age, sex, BMI or region?

`reports/assembly_healthy_strata.md` and `reports/reads_healthy_strata.md`: out-of-sample scores of healthy samples by age
bin, sex, BMI class and region, with the median fraction outside, the share of samples over the gate, the share with a
significant excess and the landscape distance. Region effects are confounded with study.

## 4. Do reference-relative features carry disease signal?

`reports/assembly_case_control_evaluation.md` with per-study tables in `reports/assembly_case_control_per_study_lenient.md`
and `..._strict.md`: leave-one-study-out AUROC of a classifier on reference-relative features versus raw CLR abundances,
alpha diversity and a published health index, on 11 held-out case/control studies. Reference-relative features win
modestly: macro-mean AUROC 0.618 vs 0.581 for raw abundances (one-sided Wilcoxon p = 0.021); on the strict baseline the
difference is not significant (p = 0.051). Gene layers alone are at least as good as all layers (0.631). For the read-based
baseline the taxonomy comparison is 0.690 vs 0.664 (p = 0.052), and 0.769 vs 0.715 (p = 0.004) when restricted to
colorectal cancer studies (`reports/reads_calibration_and_evaluation.md`). Function layers show no advantage in either
pipeline.

`reports/assembly_disease_atlas.md` and `reports/reads_disease_atlas.md`: for every case/control study, the share of cases
and controls flagged by the sample-level excess test, the median fraction outside, and per-feature deviations.
`reports/assembly_known_biology.md` tests 43 directions expected from the literature (case versus same-study control
percentile shifts) on the assembly pipeline; 14 are confirmed. `reports/assembly_disease_by_stratum.md`
and `reports/reads_disease_by_stratum.md` repeat the case/control comparison within age, sex, BMI and location strata where
a study has at least 10 cases and 10 controls, and tabulate the balance of cases and controls first.

`reports/example_reports.tsv` lists the individual sample reports used as worked examples.

## 5. Packaging round trip

The release manifest records, for each tarball, that unpacking it and re-scoring 40 test queries reproduces
`scores.parquet` and `summaries.tsv` exactly. The packaged tool was checked the same way: scoring the test inputs with
the installed package and with the pre-packaging code gives identical `scores.parquet`, `summaries.tsv` and `report.md`
(18,733 score rows for the assembly test, 531 for the read-based test).

## What is not validated

- HUMAnN 3 layers of the read-based pipeline from raw reads.
- Any body site other than stool, infants, or populations outside the pool's regions (African cohorts fail the fit).
- Individual feature calls as diagnostic statements: they are descriptive; only the sample-level excess test is a
  significance statement.
