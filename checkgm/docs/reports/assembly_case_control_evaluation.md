# S8 — evaluation report (generated 2026-09-29 20:07 UTC)

## Leakage check
No held-out study appears in its own fold reference (lenient and strict).

## Case/control LOSO AUROC — lenient (primary) reference

| feature_set                      |   macro_mean |   macro_median |   pooled_auroc |   within_study_cv_median (ceiling) |
|:---------------------------------|-------------:|---------------:|---------------:|-----------------------------------:|
| alpha_diversity                  |        0.515 |          0.52  |          0.501 |                              0.547 |
| health_index                     |        0.479 |          0.474 |          0.483 |                              0.586 |
| raw_clr                          |        0.581 |          0.591 |          0.537 |                              0.687 |
| reference_relative               |        0.618 |          0.639 |          0.572 |                              0.69  |
| reference_relative_genes_only    |        0.631 |          0.677 |          0.578 |                              0.697 |
| reference_relative_taxonomy_only |        0.585 |          0.581 |          0.539 |                              0.705 |


## Paired comparison across the 11 held-out studies (one-sided Wilcoxon, reference_relative greater)

| comparison                                            |   mean_diff |   wins |   losses |   wilcoxon_p |
|:------------------------------------------------------|------------:|-------:|---------:|-------------:|
| reference_relative − alpha_diversity                  |       0.103 |      9 |        2 |        0.012 |
| reference_relative − health_index                     |       0.139 |     10 |        1 |        0.003 |
| reference_relative − raw_clr                          |       0.037 |      9 |        2 |        0.021 |
| reference_relative − reference_relative_genes_only    |      -0.013 |      2 |        9 |        0.991 |
| reference_relative − reference_relative_taxonomy_only |       0.033 |      6 |        5 |        0.16  |

Per-study intervals: `report_loso.md`.

## Case/control LOSO AUROC — strict (sensitivity) reference

| feature_set                      |   macro_mean |   macro_median |   pooled_auroc |   within_study_cv_median (ceiling) |
|:---------------------------------|-------------:|---------------:|---------------:|-----------------------------------:|
| alpha_diversity                  |        0.515 |          0.52  |          0.501 |                              0.547 |
| health_index                     |        0.479 |          0.474 |          0.483 |                              0.586 |
| raw_clr                          |        0.581 |          0.591 |          0.537 |                              0.687 |
| reference_relative               |        0.615 |          0.643 |          0.57  |                              0.676 |
| reference_relative_genes_only    |        0.618 |          0.639 |          0.574 |                              0.703 |
| reference_relative_taxonomy_only |        0.591 |          0.575 |          0.541 |                              0.704 |


## Paired comparison across the 11 held-out studies (one-sided Wilcoxon, reference_relative greater)

| comparison                                            |   mean_diff |   wins |   losses |   wilcoxon_p |
|:------------------------------------------------------|------------:|-------:|---------:|-------------:|
| reference_relative − alpha_diversity                  |       0.101 |     10 |        1 |        0.005 |
| reference_relative − health_index                     |       0.137 |     10 |        1 |        0.001 |
| reference_relative − raw_clr                          |       0.035 |      7 |        4 |        0.051 |
| reference_relative − reference_relative_genes_only    |      -0.003 |      4 |        7 |        0.711 |
| reference_relative − reference_relative_taxonomy_only |       0.024 |      6 |        5 |        0.35  |

Per-study intervals: `report_loso_strict.md`.

## Known biology (fixed expectation list, cases vs same-study controls)

14 of 43 expectation × study tests agree in direction at p < 0.05. By condition:

| condition          |   tests |   confirmed |
|:-------------------|--------:|------------:|
| C. difficile       |       4 |           3 |
| Colorectal         |       4 |           0 |
| Hypertension       |       2 |           0 |
| Inflammatory bowel |      18 |           3 |
| Kidney stones      |       1 |           0 |
| Parkinson          |      12 |           8 |
| Rheumatoid         |       2 |           0 |

Full table: `report_known_biology.md`.

## Assembly-quality robustness (duplicate analyses of the same sample)
Spearman of per-feature percentiles over jointly assessed features; Jaccard of outside-band deviation calls; Jaccard of expected-but-missing calls. Split by whether the two analyses of the same sample fall in the same quality band.

|                            |   n_pairs |   spearman_med |   spearman_p10 |   jaccard_outside_med |   jaccard_missing_med |
|:---------------------------|----------:|---------------:|---------------:|----------------------:|----------------------:|
| ('ko_eggnog', False)       |       246 |          0.698 |          0.218 |                 0.073 |                 0     |
| ('ko_eggnog', True)        |       645 |          0.898 |          0.724 |                 0.261 |                 0.483 |
| ('ko_kofam', False)        |       246 |          0.632 |          0.172 |                 0.055 |                 0     |
| ('ko_kofam', True)         |       644 |          0.885 |          0.673 |                 0.255 |                 0.459 |
| ('module', False)          |       246 |          0.466 |          0.175 |                 0     |                 0     |
| ('module', True)           |       645 |          0.865 |          0.708 |                 0.25  |                 0.556 |
| ('pfam', False)            |       246 |          0.708 |          0.204 |                 0.08  |                 0     |
| ('pfam', True)             |       645 |          0.914 |          0.728 |                 0.308 |                 0.444 |
| ('taxonomy_family', False) |       246 |          0.316 |          0.01  |                 0     |                 0     |
| ('taxonomy_family', True)  |       629 |          0.797 |          0.597 |                 0.231 |                 0.364 |
| ('taxonomy_genus', False)  |       246 |          0.303 |          0.023 |                 0     |                 0     |
| ('taxonomy_genus', True)   |       638 |          0.814 |          0.625 |                 0.222 |                 0.385 |



## Reference calibration of the scoring bundle (S6, lenient)
Held-out studies: PRJEB39906, PRJEB62833, PRJEB75578; 233 healthy samples. Gate: fraction outside p2.5–p97.5 ≤ 8% per study (hard); < 3% reported as a diagnostic.

| layer           | study_bioproject   |   n_samples |   frac_outside_mean |   frac_below |   frac_above | gate_pass   | below_floor_diagnostic   |
|:----------------|:-------------------|------------:|--------------------:|-------------:|-------------:|:------------|:-------------------------|
| ko_eggnog       | PRJEB39906         |          72 |               0.051 |        0.02  |        0.029 | True        | False                    |
| ko_eggnog       | PRJEB62833         |         133 |               0.13  |        0.052 |        0.077 | False       | False                    |
| ko_eggnog       | PRJEB75578         |          28 |               0.021 |        0.014 |        0.007 | True        | True                     |
| ko_kofam        | PRJEB39906         |          72 |               0.064 |        0.022 |        0.04  | True        | False                    |
| ko_kofam        | PRJEB62833         |         133 |               0.131 |        0.047 |        0.083 | False       | False                    |
| ko_kofam        | PRJEB75578         |          28 |               0.023 |        0.016 |        0.006 | True        | True                     |
| module          | PRJEB39906         |          72 |               0.003 |        0.002 |        0.001 | True        | True                     |
| module          | PRJEB62833         |         133 |               0.009 |        0.008 |        0.002 | True        | True                     |
| module          | PRJEB75578         |          28 |               0.002 |        0.001 |        0     | True        | True                     |
| pfam            | PRJEB39906         |          72 |               0.065 |        0.018 |        0.045 | True        | False                    |
| pfam            | PRJEB62833         |         133 |               0.129 |        0.058 |        0.07  | False       | False                    |
| pfam            | PRJEB75578         |          28 |               0.024 |        0.012 |        0.011 | True        | True                     |
| taxonomy_family | PRJEB39906         |          72 |               0.076 |        0.043 |        0.011 | True        | False                    |
| taxonomy_family | PRJEB62833         |         133 |               0.049 |        0.015 |        0.033 | True        | False                    |
| taxonomy_family | PRJEB75578         |          28 |               0.031 |        0.012 |        0.018 | True        | False                    |
| taxonomy_genus  | PRJEB39906         |          72 |               0.079 |        0.044 |        0.012 | True        | False                    |
| taxonomy_genus  | PRJEB62833         |         133 |               0.054 |        0.014 |        0.039 | True        | False                    |
| taxonomy_genus  | PRJEB75578         |          28 |               0.029 |        0.015 |        0.012 | True        | True                     |


## Ablations
- strict vs lenient reference: see the two LOSO sections above.
- taxonomy-only vs gene-only vs combined reference-relative features: in the lenient LOSO summary.
- calibration on/off: not applicable (assembly query on assembly reference; calibration = none).
- Tier A only: not applicable (S4 reliability tiers not built; cMD read-based layer deferred).
