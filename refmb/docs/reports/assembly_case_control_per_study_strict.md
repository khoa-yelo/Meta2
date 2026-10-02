# S8 — leave-one-study-out case/control AUROC (11 studies, 2,294 in-stratum samples, 1,407 cases)

Same random forest (500 trees, balanced class weights, seed 20261020), same folds, for every feature set. Each study's features were computed against a reference built without that study.

## Per held-out study

|                                                        |   alpha_diversity |   health_index |   raw_clr |   reference_relative |   reference_relative_genes_only |   reference_relative_taxonomy_only |
|:-------------------------------------------------------|------------------:|---------------:|----------:|---------------------:|--------------------------------:|-----------------------------------:|
| ('PRJEB25962', 'Inflammatory bowel disease', 379, 143) |             0.549 |          0.491 |     0.458 |                0.554 |                           0.565 |                              0.47  |
| ('PRJEB26165', 'C. difficile infection', 70, 56)       |             0.564 |          0.365 |     0.596 |                0.749 |                           0.821 |                              0.474 |
| ('PRJEB29015', 'Rheumatoid arthritis', 232, 135)       |             0.458 |          0.479 |     0.551 |                0.573 |                           0.56  |                              0.575 |
| ('PRJEB41218', 'Kidney stones', 56, 31)                |             0.52  |          0.632 |     0.652 |                0.643 |                           0.685 |                              0.68  |
| ('PRJEB43973', 'Inflammatory bowel disease', 162, 142) |             0.556 |          0.458 |     0.77  |                0.75  |                           0.767 |                              0.778 |
| ('PRJEB48911', 'Seasonal allergies', 391, 289)         |             0.492 |          0.497 |     0.541 |                0.55  |                           0.558 |                              0.521 |
| ('PRJEB52172', "Parkinson's disease", 74, 37)          |             0.536 |          0.438 |     0.476 |                0.436 |                           0.369 |                              0.61  |
| ('PRJEB61255', "Parkinson's disease", 692, 468)        |             0.45  |          0.474 |     0.422 |                0.479 |                           0.485 |                              0.508 |
| ('PRJEB62821', 'Hemodialysis', 40, 24)                 |             0.549 |          0.299 |     0.591 |                0.648 |                           0.685 |                              0.648 |
| ('PRJEB62825', 'Hypertension', 116, 60)                |             0.506 |          0.459 |     0.693 |                0.676 |                           0.667 |                              0.665 |
| ('PRJEB67980', 'Colorectal neoplasia', 82, 22)         |             0.482 |          0.673 |     0.638 |                0.711 |                           0.639 |                              0.575 |


## Summary

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


## Per-study 95% bootstrap CI (300 resamples of the test set)

| study      | alpha_diversity    | health_index       | raw_clr            | reference_relative   | reference_relative_genes_only   | reference_relative_taxonomy_only   |
|:-----------|:-------------------|:-------------------|:-------------------|:---------------------|:--------------------------------|:-----------------------------------|
| PRJEB25962 | 0.549 [0.49, 0.60] | 0.491 [0.44, 0.55] | 0.458 [0.39, 0.52] | 0.554 [0.49, 0.62]   | 0.565 [0.51, 0.63]              | 0.470 [0.41, 0.53]                 |
| PRJEB26165 | 0.564 [0.42, 0.69] | 0.365 [0.17, 0.58] | 0.596 [0.41, 0.77] | 0.749 [0.61, 0.86]   | 0.821 [0.71, 0.91]              | 0.474 [0.30, 0.63]                 |
| PRJEB29015 | 0.458 [0.38, 0.54] | 0.479 [0.41, 0.55] | 0.551 [0.48, 0.64] | 0.573 [0.50, 0.65]   | 0.560 [0.49, 0.64]              | 0.575 [0.49, 0.66]                 |
| PRJEB41218 | 0.520 [0.35, 0.67] | 0.632 [0.48, 0.76] | 0.652 [0.49, 0.78] | 0.643 [0.47, 0.79]   | 0.685 [0.53, 0.83]              | 0.680 [0.53, 0.81]                 |
| PRJEB43973 | 0.556 [0.39, 0.70] | 0.458 [0.32, 0.59] | 0.770 [0.65, 0.88] | 0.750 [0.63, 0.85]   | 0.767 [0.65, 0.86]              | 0.778 [0.68, 0.87]                 |
| PRJEB48911 | 0.492 [0.43, 0.55] | 0.497 [0.43, 0.56] | 0.541 [0.47, 0.60] | 0.550 [0.48, 0.62]   | 0.558 [0.49, 0.62]              | 0.521 [0.46, 0.58]                 |
| PRJEB52172 | 0.536 [0.41, 0.66] | 0.438 [0.32, 0.58] | 0.476 [0.35, 0.63] | 0.436 [0.30, 0.59]   | 0.369 [0.26, 0.51]              | 0.610 [0.49, 0.73]                 |
| PRJEB61255 | 0.450 [0.41, 0.50] | 0.474 [0.43, 0.52] | 0.422 [0.37, 0.47] | 0.479 [0.43, 0.53]   | 0.485 [0.43, 0.53]              | 0.508 [0.46, 0.56]                 |
| PRJEB62821 | 0.549 [0.36, 0.72] | 0.299 [0.13, 0.47] | 0.591 [0.39, 0.77] | 0.648 [0.47, 0.81]   | 0.685 [0.51, 0.84]              | 0.648 [0.42, 0.83]                 |
| PRJEB62825 | 0.506 [0.39, 0.59] | 0.459 [0.35, 0.55] | 0.693 [0.60, 0.78] | 0.676 [0.57, 0.77]   | 0.667 [0.57, 0.76]              | 0.665 [0.57, 0.77]                 |
| PRJEB67980 | 0.482 [0.34, 0.62] | 0.673 [0.56, 0.81] | 0.638 [0.49, 0.79] | 0.711 [0.58, 0.83]   | 0.639 [0.51, 0.77]              | 0.575 [0.44, 0.72]                 |


## Per condition group (≥ 2 studies)

|                                   |   alpha_diversity |   health_index |   raw_clr |   reference_relative |   reference_relative_genes_only |   reference_relative_taxonomy_only |
|:----------------------------------|------------------:|---------------:|----------:|---------------------:|--------------------------------:|-----------------------------------:|
| ('Inflammatory bowel disease', 2) |             0.552 |          0.475 |     0.614 |                0.652 |                           0.666 |                              0.624 |
| ("Parkinson's disease", 2)        |             0.493 |          0.456 |     0.449 |                0.457 |                           0.427 |                              0.559 |

