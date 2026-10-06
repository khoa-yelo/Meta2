# How the baselines were validated

The reports in `docs/reports/` are verbatim copies of the result files the numbers in the paper and in the README come from.
This page says what each one tested and where its headline numbers are. Because the copies are verbatim, they still carry
the development workspace's internal labels: S6 = the held-out calibration check (section 2), S7 = the excess-test survey over
8,742 MGnify gut analyses (`assembly_all_gut_analyses_scoring.md`), S8 = the case/control evaluation (section 4), S11 = the
reproduction from raw reads (section 1), S12 and S14 = the reference strata and the disease atlas (sections 3 and 4), "anchors"
= samples measured by both pipelines. Outputs and workspace scripts produced before 2026-10-06 call the two pipelines
"tier A" and "tier B", and some file names still carry those codenames; they mean the assembly-based and the read-based
pipeline. The software itself now says only "assembly" and "read".

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

## 2. Do reference samples that were never in the pool fit the baseline?

Expectation: about 5 % of features outside the 2.5–97.5 percentile band; the gate used during development was 8 % per study.

Every figure in sections 2 to 4 rests on one condition — that no sample was scored against a baseline holding it — and
that condition is only as good as the test for whether two deposits are the same specimen. In the read-based pipeline
they are tested two ways and the results unioned (`project/scripts/s11_cross_study_duplicates.py`, output
`work/s11/pipeB/cross_study_duplicates.tsv`): a shared ENA or SRA run accession, and an identical
(`number_reads`, `number_bases`) pair in the cMD3 metadata. Neither test alone is enough. The audit finds 201 pairs
among the 9,123 scored samples. 30 are DhakanDB_2019 / GuptaA_2019, caught by both tests, and the profiles are
effectively identical (species Bray–Curtis similarity 0.9999 to 1.00). The other 171 are LeChatelierE_2013 /
NielsenHB_2014, two curatedMetagenomicData 3 deposits of one MetaHIT Danish cohort: all 171 agree on cMD3 sample id,
subject id, country, sequencing platform, `number_reads` and `number_bases`, but each deposit profiles the specimen
from a differently partitioned set of runs (MH0006, for instance, from 6 runs in one deposit and 12 in the other, with
no accession in common), so the two profiles differ — Bray–Curtis similarity median 0.90, minimum 0.69 — and only
MH0001 to MH0005 share any run accession at all. An accession-only audit therefore records 5 of those 171 pairs, which
is what this project did until 2026-10-06.

Two consequences follow, and both are now handled in `project/scripts/s12_oos_scores.py`. A leave-one-study-out
baseline holds out every study that shares a specimen with the one being scored, not only that study, so the two
MetaHIT deposits are held out together; and the exclusion applies to a whole study rather than to its duplicated
samples alone, so a study's cases and its controls are never scored against two different baselines. Checked sample by
sample, 0 of the 9,123 scored samples now sees itself or another deposit of itself in its baseline, against 342 under
the per-sample rule, and no study is split across two baselines, against GuptaA_2019 under the per-sample rule.

One fact about the pool composition is not a leakage question and remains: 342 of the 6,494 read-based pool samples are
171 mirrored pairs, so the pool holds **6,323 distinct specimens**, and those 171 donors carry twice the weight of the
rest. Two MetaCardis_2020_a pool samples have no row in the cMD3 metadata table and cannot be audited either way.

- **Assembly** — `reports/assembly_heldout_calibration_check.md` (three held-out reference studies, 233 samples; per layer and
  study the mean fraction outside and whether it passes) and `reports/assembly_loso_within_pool.md` (each pool study scored
  against a reference rebuilt without it: medians 0.05–0.055 on all layers except module; 73–85 % of studies within 8 %).
  Cohorts from Africa do not fit (19.8 % of families outside), and location cannot be separated from study.
- **Reads** — `reports/reads_calibration_and_evaluation.md` (three held-out reference studies and in-pool leave-one-study-out,
  all five layers).

## 3. Does the reference baseline move with age, sex, BMI or region?

`reports/assembly_healthy_strata.md` and `reports/reads_healthy_strata.md`: out-of-sample scores of reference samples by age
bin, sex, BMI class and region, with the median fraction outside, the share of samples over the gate, the share with a
significant excess and the landscape distance. Region effects are confounded with study.

## 4. Do reference-relative features carry disease signal?

`reports/assembly_case_control_evaluation.md` with per-study tables in `reports/assembly_case_control_per_study_lenient.md`
and `..._strict.md`: leave-one-study-out AUROC of a classifier on reference-relative features versus raw CLR abundances,
alpha diversity and a published health index, on 11 held-out case/control studies (at least 20 cases and 10 controls each;
the read-based evaluation below uses 14 studies with at least 20 cases and 20 controls, one of them a largely adolescent
cohort, see `configs/evaluation_settings.md`). Reference-relative features win
modestly: macro-mean AUROC 0.618 vs 0.581 for raw abundances (one-sided Wilcoxon p = 0.021); on the strict baseline the
difference is not significant (p = 0.051). Gene layers alone are at least as good as all layers (0.631). For the read-based
baseline the taxonomy comparison is 0.688 vs 0.664 (10 of 14 studies, p = 0.034), and 0.770 vs 0.715 (8 of 9, p = 0.004)
when restricted to colorectal cancer studies (`reports/reads_calibration_and_evaluation.md`, the run of 2026-10-06 16:12
UTC). Function layers show no advantage in either pipeline.

`reports/assembly_disease_atlas.md` and `reports/reads_disease_atlas.md`: for every case/control study, the share of cases
and controls flagged by the sample-level excess test, the median fraction outside, and per-feature deviations.
`reports/assembly_known_biology.md` tests 43 directions expected from the literature (case versus same-study control
percentile shifts) on the assembly-based pipeline; 14 are confirmed. `reports/assembly_all_gut_analyses_scoring.md` scores
all 8,742 gut analyses of the inventory against the standard assembly-based baseline and tabulates the fraction outside
by health label (reference samples included in the pool are scored in-sample there). `reports/assembly_disease_by_stratum.md`
and `reports/reads_disease_by_stratum.md` repeat the case/control comparison within age, sex, BMI and location strata where
a study has at least 10 cases and 10 controls, and tabulate the balance of cases and controls first.

`reports/example_reports.tsv` lists the individual sample reports used as worked examples; `reports/example_case_counts.tsv` holds every count of the C. difficile worked example (all-layer counts, the 232 common families, and the share of the cohort low or high per family).

## 5. Packaging round trip

The release manifest records, for each tarball, that unpacking it and re-scoring 40 test queries reproduces
`scores.parquet` and `summaries.tsv` exactly. The packaged tool was checked the same way: scoring the test inputs with
the installed package and with the pre-packaging code gives identical `scores.parquet`, `summaries.tsv` and `report.md`
(18,733 score rows for the assembly test, 531 for the read-based test).

## 6. Run time of the scoring step

Measured once with `/usr/bin/time -v` on the README example (`checkgm run-metaphlan` on `examples/synthetic_healthy_adult.txt`,
102 species, against `gut-reads-adult-global-v0.2-lenient`, writing scores, summaries and report): wall time 16.2 s,
105 % of one CPU, 2.8 GB peak resident memory (most of it the NCBI backbone of the bundle), on an AMD EPYC 7742 host
with Python 3.12.14, pandas 3.0.6, numpy 2.5.3 and pyarrow 25.0.0, from a warm file cache. Scoring an assembly-based test set of two MGnify analyses
(18,675 score rows) took about 2 s in the same setting; normalizing a query directory took about 19 s. These are single
measurements, not benchmarks.

## Settings and scripts

`docs/configs/` holds read-only copies of the configuration files the baselines and the evaluation were run with
(`s5.yaml`, `inclusion.yaml`, `splits_rules.yaml`, `s11.yaml`) and the evaluation settings transcribed from the
evaluation script; `scripts/viz/` the scripts that draw every figure of the paper. Both read the project workspace, not
the repository (see their READMEs).

## What is not validated

- HUMAnN 3 layers of the read-based pipeline from raw reads.
- Any body site other than stool, infants, or populations outside the pool's regions (African cohorts fail the fit).
- Individual feature calls as diagnostic statements: they are descriptive; only the sample-level excess test is a
  significance statement.
