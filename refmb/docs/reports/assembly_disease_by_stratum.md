# S14 — disease deviations by age, sex, BMI and location, pipeline A — 2026-09-30 00:28 UTC

10 case/control studies (>= 15 cases and >= 15 controls), 2,539 samples. A study x level cell is analysed when it holds >= 10 cases and >= 10 controls: 27 cells in 7 studies.

## 1. Balance of cases and controls

| study      | condition                  | region                                                             |   n_case |   n_control |   age_known |   age_case |   age_control |   age_p |   bmi_known |   bmi_case |   bmi_control |   bmi_p |   sex_known |   female_case |   female_control |   sex_p |
|:-----------|:---------------------------|:-------------------------------------------------------------------|---------:|------------:|------------:|-----------:|--------------:|--------:|------------:|-----------:|--------------:|--------:|------------:|--------------:|-----------------:|--------:|
| PRJEB25962 | Inflammatory bowel disease | Europe                                                             |      148 |         246 |        0.96 |       42.5 |          54   |   0     |        0.96 |     23.309 |        27.574 |   0     |        0.96 |         0.65  |            0.521 |   0.019 |
| PRJEB29015 | Rheumatoid arthritis       | E. Asia                                                            |      302 |         195 |        0    |      nan   |         nan   | nan     |        0    |    nan     |       nan     | nan     |        0.46 |         0.805 |            0.698 |   0.088 |
| PRJEB41218 | Kidney stones              | N. America                                                         |       31 |          25 |        1    |       62   |          58   |   0.184 |        0    |    nan     |       nan     | nan     |        1    |         0.387 |            0.36  |   1     |
| PRJEB43973 | Inflammatory bowel disease | N. America                                                         |      143 |          20 |        0    |      nan   |         nan   | nan     |        0    |    nan     |       nan     | nan     |        0    |         0     |            0     | nan     |
| PRJEB48911 | Seasonal allergies         | E. Asia/Europe/N. America/Oceania/S. America/S. Asia/SE Asia/other |      314 |         110 |        0.91 |       49   |          48.5 |   0.804 |        0    |    nan     |       nan     | nan     |        0.89 |         0.561 |            0.444 |   0.06  |
| PRJEB52172 | Parkinson's disease        | E. Asia                                                            |       37 |          37 |        0    |      nan   |         nan   | nan     |        0    |    nan     |       nan     | nan     |        0    |         0     |            0     | nan     |
| PRJEB61255 | Parkinson's disease        | N. America                                                         |      468 |         224 |        1    |       69   |          67   |   0     |        0    |    nan     |       nan     | nan     |        1    |         0.368 |            0.692 |   0     |
| PRJEB62821 | Hemodialysis               |                                                                    |       24 |          16 |        0    |      nan   |         nan   | nan     |        0    |    nan     |       nan     | nan     |        0    |         0     |            0     | nan     |
| PRJEB62825 | Hypertension               | E. Asia                                                            |       61 |          56 |        0.99 |       56   |          55   |   0.588 |        0    |    nan     |       nan     | nan     |        0.99 |         0.417 |            0.464 |   0.742 |
| PRJEB67980 | Colorectal neoplasia       | Europe                                                             |       22 |          60 |        1    |       62.5 |          59   |   0.005 |        1    |     25.5   |        24.6   |   0.683 |        1    |         0.5   |            0.467 |   0.985 |

## 2. Cells analysed

|                                              |   ('age_bin', '18-39') |   ('age_bin', '40-64') |   ('age_bin', '65+') |   ('bmi_class', '18.5-25') |   ('bmi_class', '25-30') |   ('bmi_class', '>=30') |   ('sex', 'female') |   ('sex', 'male') |
|:---------------------------------------------|-----------------------:|-----------------------:|---------------------:|---------------------------:|-------------------------:|------------------------:|--------------------:|------------------:|
| ('PRJEB25962', 'Inflammatory bowel disease') |                     54 |                     86 |                    0 |                         84 |                       38 |                      14 |                  93 |                50 |
| ('PRJEB29015', 'Rheumatoid arthritis')       |                      0 |                      0 |                    0 |                          0 |                        0 |                       0 |                 107 |                26 |
| ('PRJEB41218', 'Kidney stones')              |                      0 |                     15 |                    0 |                          0 |                        0 |                       0 |                   0 |                19 |
| ('PRJEB48911', 'Seasonal allergies')         |                    116 |                     83 |                   85 |                          0 |                        0 |                       0 |                 156 |               122 |
| ('PRJEB61255', "Parkinson's disease")        |                      0 |                    124 |                  342 |                          0 |                        0 |                       0 |                 172 |               296 |
| ('PRJEB62825', 'Hypertension')               |                      0 |                     47 |                   13 |                          0 |                        0 |                       0 |                  25 |                35 |
| ('PRJEB67980', 'Colorectal neoplasia')       |                      0 |                     13 |                    0 |                          0 |                        0 |                       0 |                  11 |                11 |

(number of cases in each analysed cell; 0 = cell not analysed)

## 3. Shifts adjusted for the variable vs crude shifts (atlas-significant features: q <= 0.05 and |shift| >= 10)

| study      | condition                  | variable   |   levels_used |   features |   spearman_crude_vs_adjusted |   atlas_significant |   same_sign |   retained_half_magnitude |   median_ratio_adjusted_to_crude |
|:-----------|:---------------------------|:-----------|--------------:|-----------:|-----------------------------:|--------------------:|------------:|--------------------------:|---------------------------------:|
| PRJEB25962 | Inflammatory bowel disease | age_bin    |             2 |      12395 |                        0.926 |                2992 |       1     |                     0.993 |                            0.918 |
| PRJEB25962 | Inflammatory bowel disease | bmi_class  |             3 |      11981 |                        0.89  |                2980 |       1     |                     0.991 |                            0.928 |
| PRJEB25962 | Inflammatory bowel disease | sex        |             2 |      12592 |                        0.938 |                2998 |       1     |                     0.999 |                            0.994 |
| PRJEB29015 | Rheumatoid arthritis       | sex        |             2 |      12949 |                        0.885 |                   3 |       1     |                     1     |                            0.981 |
| PRJEB41218 | Kidney stones              | age_bin    |             1 |       9180 |                        0.727 |                   0 |     nan     |                   nan     |                          nan     |
| PRJEB41218 | Kidney stones              | sex        |             1 |       9391 |                        0.81  |                   0 |     nan     |                   nan     |                          nan     |
| PRJEB48911 | Seasonal allergies         | age_bin    |             3 |      11564 |                        0.819 |                   0 |     nan     |                   nan     |                          nan     |
| PRJEB48911 | Seasonal allergies         | sex        |             2 |      12032 |                        0.856 |                   0 |     nan     |                   nan     |                          nan     |
| PRJEB61255 | Parkinson's disease        | age_bin    |             2 |      15366 |                        0.929 |                1956 |       1     |                     1     |                            0.98  |
| PRJEB61255 | Parkinson's disease        | sex        |             2 |      15156 |                        0.891 |                1958 |       0.999 |                     0.995 |                            1.015 |
| PRJEB62825 | Hypertension               | age_bin    |             2 |      11053 |                        0.881 |                   1 |       1     |                     1     |                            1.182 |
| PRJEB62825 | Hypertension               | sex        |             2 |      10302 |                        0.898 |                   1 |       1     |                     1     |                            0.985 |
| PRJEB67980 | Colorectal neoplasia       | age_bin    |             1 |       7927 |                        0.808 |                 201 |       1     |                     0.91  |                            0.819 |
| PRJEB67980 | Colorectal neoplasia       | sex        |             2 |       7755 |                        0.914 |                 198 |       1     |                     1     |                            0.967 |

## 4. Agreement of shifts between levels within a study (medians over studies and layers)

| condition                  | variable   | level_a   | level_b   |   studies |   spearman_all |   same_sign_significant |
|:---------------------------|:-----------|:----------|:----------|----------:|---------------:|------------------------:|
| Colorectal neoplasia       | sex        | female    | male      |         1 |          0.165 |                   0.992 |
| Hypertension               | age_bin    | 40-64     | 65+       |         1 |         -0.002 |                 nan     |
| Hypertension               | sex        | female    | male      |         1 |          0.152 |                   1     |
| Inflammatory bowel disease | age_bin    | 18-39     | 40-64     |         1 |          0.211 |                   0.862 |
| Inflammatory bowel disease | bmi_class  | 18.5-25   | 25-30     |         1 |          0.206 |                   0.846 |
| Inflammatory bowel disease | bmi_class  | 18.5-25   | >=30      |         1 |          0.442 |                   0.942 |
| Inflammatory bowel disease | bmi_class  | 25-30     | >=30      |         1 |          0.447 |                   0.794 |
| Inflammatory bowel disease | sex        | female    | male      |         1 |          0.469 |                   0.992 |
| Parkinson's disease        | age_bin    | 40-64     | 65+       |         1 |          0.261 |                   0.934 |
| Parkinson's disease        | sex        | female    | male      |         1 |          0.523 |                   0.998 |
| Rheumatoid arthritis       | sex        | female    | male      |         1 |          0.001 |                   1     |
| Seasonal allergies         | age_bin    | 18-39     | 40-64     |         1 |         -0.118 |                 nan     |
| Seasonal allergies         | age_bin    | 18-39     | 65+       |         1 |          0.039 |                 nan     |
| Seasonal allergies         | age_bin    | 40-64     | 65+       |         1 |          0.005 |                 nan     |
| Seasonal allergies         | sex        | female    | male      |         1 |          0.12  |                 nan     |

## 5. Features consistent in the atlas: do they hold within each level and region?

| condition                  | variable   | level      |   features |   studies |   share_same_direction |   median_abs_shift |
|:---------------------------|:-----------|:-----------|-----------:|----------:|-----------------------:|-------------------:|
| Inflammatory bowel disease | age_bin    | 18-39      |       1041 |         1 |                  0.894 |              9.958 |
| Inflammatory bowel disease | age_bin    | 40-64      |       1042 |         1 |                  1     |             17.685 |
| Inflammatory bowel disease | bmi_class  | 18.5-25    |       1042 |         1 |                  0.995 |             17.439 |
| Inflammatory bowel disease | bmi_class  | 25-30      |       1040 |         1 |                  0.897 |             12.084 |
| Inflammatory bowel disease | bmi_class  | >=30       |       1031 |         1 |                  0.95  |             18.289 |
| Inflammatory bowel disease | region     | Europe     |       1042 |         1 |                  1     |             16.562 |
| Inflammatory bowel disease | region     | N. America |       1042 |         1 |                  1     |             23.287 |
| Inflammatory bowel disease | sex        | female     |       1042 |         1 |                  1     |             17.599 |
| Inflammatory bowel disease | sex        | male       |       1041 |         1 |                  0.992 |             15.141 |
| Parkinson's disease        | age_bin    | 40-64      |         25 |         1 |                  0.96  |              9.833 |
| Parkinson's disease        | age_bin    | 65+        |         25 |         1 |                  1     |             16.107 |
| Parkinson's disease        | region     | E. Asia    |         25 |         1 |                  1     |             24.387 |
| Parkinson's disease        | region     | N. America |         25 |         1 |                  1     |             13.432 |
| Parkinson's disease        | sex        | female     |         25 |         1 |                  1     |             11.26  |
| Parkinson's disease        | sex        | male       |         25 |         1 |                  1     |             14.923 |

(features = atlas-consistent features assessable in the level; share_same_direction = share whose median shift over studies has the atlas direction)

## 6. Sample level: median fraction outside, cases vs controls, by level (medians over studies)

| condition                  | layer           | variable   | level   |   studies |   case |   control |
|:---------------------------|:----------------|:-----------|:--------|----------:|-------:|----------:|
| Colorectal neoplasia       | ko_eggnog       | age_bin    | 40-64   |         1 |  0.043 |     0.027 |
| Colorectal neoplasia       | ko_eggnog       | sex        | female  |         1 |  0.038 |     0.026 |
| Colorectal neoplasia       | ko_eggnog       | sex        | male    |         1 |  0.071 |     0.024 |
| Colorectal neoplasia       | ko_kofam        | age_bin    | 40-64   |         1 |  0.042 |     0.026 |
| Colorectal neoplasia       | ko_kofam        | sex        | female  |         1 |  0.042 |     0.024 |
| Colorectal neoplasia       | ko_kofam        | sex        | male    |         1 |  0.058 |     0.025 |
| Colorectal neoplasia       | module          | age_bin    | 40-64   |         1 |  0     |     0     |
| Colorectal neoplasia       | module          | sex        | female  |         1 |  0.005 |     0     |
| Colorectal neoplasia       | module          | sex        | male    |         1 |  0     |     0     |
| Colorectal neoplasia       | pfam            | age_bin    | 40-64   |         1 |  0.052 |     0.029 |
| Colorectal neoplasia       | pfam            | sex        | female  |         1 |  0.04  |     0.028 |
| Colorectal neoplasia       | pfam            | sex        | male    |         1 |  0.064 |     0.025 |
| Colorectal neoplasia       | taxonomy_family | age_bin    | 40-64   |         1 |  0.034 |     0.031 |
| Colorectal neoplasia       | taxonomy_family | sex        | female  |         1 |  0.036 |     0.03  |
| Colorectal neoplasia       | taxonomy_family | sex        | male    |         1 |  0.039 |     0.034 |
| Colorectal neoplasia       | taxonomy_genus  | age_bin    | 40-64   |         1 |  0.039 |     0.039 |
| Colorectal neoplasia       | taxonomy_genus  | sex        | female  |         1 |  0.051 |     0.037 |
| Colorectal neoplasia       | taxonomy_genus  | sex        | male    |         1 |  0.035 |     0.038 |
| Hypertension               | ko_eggnog       | age_bin    | 40-64   |         1 |  0.03  |     0.022 |
| Hypertension               | ko_eggnog       | age_bin    | 65+     |         1 |  0.032 |     0.018 |
| Hypertension               | ko_eggnog       | sex        | female  |         1 |  0.035 |     0.021 |
| Hypertension               | ko_eggnog       | sex        | male    |         1 |  0.023 |     0.023 |
| Hypertension               | ko_kofam        | age_bin    | 40-64   |         1 |  0.025 |     0.024 |
| Hypertension               | ko_kofam        | age_bin    | 65+     |         1 |  0.037 |     0.018 |
| Hypertension               | ko_kofam        | sex        | female  |         1 |  0.037 |     0.022 |
| Hypertension               | ko_kofam        | sex        | male    |         1 |  0.02  |     0.026 |
| Hypertension               | module          | age_bin    | 40-64   |         1 |  0     |     0     |
| Hypertension               | module          | age_bin    | 65+     |         1 |  0     |     0.002 |
| Hypertension               | module          | sex        | female  |         1 |  0     |     0     |
| Hypertension               | module          | sex        | male    |         1 |  0     |     0     |
| Hypertension               | pfam            | age_bin    | 40-64   |         1 |  0.028 |     0.02  |
| Hypertension               | pfam            | age_bin    | 65+     |         1 |  0.03  |     0.02  |
| Hypertension               | pfam            | sex        | female  |         1 |  0.04  |     0.019 |
| Hypertension               | pfam            | sex        | male    |         1 |  0.022 |     0.026 |
| Hypertension               | taxonomy_family | age_bin    | 40-64   |         1 |  0.043 |     0.033 |
| Hypertension               | taxonomy_family | age_bin    | 65+     |         1 |  0.041 |     0.031 |
| Hypertension               | taxonomy_family | sex        | female  |         1 |  0.052 |     0.038 |
| Hypertension               | taxonomy_family | sex        | male    |         1 |  0.042 |     0.026 |
| Hypertension               | taxonomy_genus  | age_bin    | 40-64   |         1 |  0.045 |     0.032 |
| Hypertension               | taxonomy_genus  | age_bin    | 65+     |         1 |  0.046 |     0.032 |
| Hypertension               | taxonomy_genus  | sex        | female  |         1 |  0.049 |     0.031 |
| Hypertension               | taxonomy_genus  | sex        | male    |         1 |  0.043 |     0.035 |
| Inflammatory bowel disease | ko_eggnog       | age_bin    | 18-39   |         1 |  0.027 |     0.022 |
| Inflammatory bowel disease | ko_eggnog       | age_bin    | 40-64   |         1 |  0.021 |     0.028 |
| Inflammatory bowel disease | ko_eggnog       | bmi_class  | 18.5-25 |         1 |  0.025 |     0.026 |
| Inflammatory bowel disease | ko_eggnog       | bmi_class  | 25-30   |         1 |  0.02  |     0.026 |
| Inflammatory bowel disease | ko_eggnog       | bmi_class  | >=30    |         1 |  0.018 |     0.026 |
| Inflammatory bowel disease | ko_eggnog       | sex        | female  |         1 |  0.023 |     0.025 |
| Inflammatory bowel disease | ko_eggnog       | sex        | male    |         1 |  0.022 |     0.029 |
| Inflammatory bowel disease | ko_kofam        | age_bin    | 18-39   |         1 |  0.027 |     0.022 |
| Inflammatory bowel disease | ko_kofam        | age_bin    | 40-64   |         1 |  0.023 |     0.034 |
| Inflammatory bowel disease | ko_kofam        | bmi_class  | 18.5-25 |         1 |  0.024 |     0.028 |
| Inflammatory bowel disease | ko_kofam        | bmi_class  | 25-30   |         1 |  0.023 |     0.031 |
| Inflammatory bowel disease | ko_kofam        | bmi_class  | >=30    |         1 |  0.017 |     0.034 |
| Inflammatory bowel disease | ko_kofam        | sex        | female  |         1 |  0.022 |     0.031 |
| Inflammatory bowel disease | ko_kofam        | sex        | male    |         1 |  0.026 |     0.031 |
| Inflammatory bowel disease | module          | age_bin    | 18-39   |         1 |  0     |     0     |
| Inflammatory bowel disease | module          | age_bin    | 40-64   |         1 |  0     |     0     |
| Inflammatory bowel disease | module          | bmi_class  | 18.5-25 |         1 |  0     |     0     |
| Inflammatory bowel disease | module          | bmi_class  | 25-30   |         1 |  0     |     0     |
| Inflammatory bowel disease | module          | bmi_class  | >=30    |         1 |  0     |     0     |
| Inflammatory bowel disease | module          | sex        | female  |         1 |  0     |     0     |
| Inflammatory bowel disease | module          | sex        | male    |         1 |  0     |     0     |
| Inflammatory bowel disease | pfam            | age_bin    | 18-39   |         1 |  0.027 |     0.025 |
| Inflammatory bowel disease | pfam            | age_bin    | 40-64   |         1 |  0.022 |     0.03  |
| Inflammatory bowel disease | pfam            | bmi_class  | 18.5-25 |         1 |  0.024 |     0.029 |
| Inflammatory bowel disease | pfam            | bmi_class  | 25-30   |         1 |  0.024 |     0.031 |
| Inflammatory bowel disease | pfam            | bmi_class  | >=30    |         1 |  0.019 |     0.03  |
| Inflammatory bowel disease | pfam            | sex        | female  |         1 |  0.022 |     0.029 |
| Inflammatory bowel disease | pfam            | sex        | male    |         1 |  0.024 |     0.03  |
| Inflammatory bowel disease | taxonomy_family | age_bin    | 18-39   |         1 |  0.039 |     0.036 |
| Inflammatory bowel disease | taxonomy_family | age_bin    | 40-64   |         1 |  0.036 |     0.036 |
| Inflammatory bowel disease | taxonomy_family | bmi_class  | 18.5-25 |         1 |  0.039 |     0.034 |
| Inflammatory bowel disease | taxonomy_family | bmi_class  | 25-30   |         1 |  0.036 |     0.038 |
| Inflammatory bowel disease | taxonomy_family | bmi_class  | >=30    |         1 |  0.045 |     0.039 |
| Inflammatory bowel disease | taxonomy_family | sex        | female  |         1 |  0.039 |     0.034 |
| Inflammatory bowel disease | taxonomy_family | sex        | male    |         1 |  0.037 |     0.038 |
| Inflammatory bowel disease | taxonomy_genus  | age_bin    | 18-39   |         1 |  0.043 |     0.036 |
| Inflammatory bowel disease | taxonomy_genus  | age_bin    | 40-64   |         1 |  0.036 |     0.04  |
| Inflammatory bowel disease | taxonomy_genus  | bmi_class  | 18.5-25 |         1 |  0.037 |     0.04  |
| Inflammatory bowel disease | taxonomy_genus  | bmi_class  | 25-30   |         1 |  0.041 |     0.04  |
| Inflammatory bowel disease | taxonomy_genus  | bmi_class  | >=30    |         1 |  0.043 |     0.039 |
| Inflammatory bowel disease | taxonomy_genus  | sex        | female  |         1 |  0.039 |     0.04  |
| Inflammatory bowel disease | taxonomy_genus  | sex        | male    |         1 |  0.038 |     0.039 |
| Kidney stones              | ko_eggnog       | age_bin    | 40-64   |         1 |  0.032 |     0.024 |
| Kidney stones              | ko_eggnog       | sex        | male    |         1 |  0.047 |     0.021 |
| Kidney stones              | ko_kofam        | age_bin    | 40-64   |         1 |  0.032 |     0.019 |
| Kidney stones              | ko_kofam        | sex        | male    |         1 |  0.053 |     0.018 |
| Kidney stones              | module          | age_bin    | 40-64   |         1 |  0     |     0     |
| Kidney stones              | module          | sex        | male    |         1 |  0.004 |     0     |
| Kidney stones              | pfam            | age_bin    | 40-64   |         1 |  0.032 |     0.023 |
| Kidney stones              | pfam            | sex        | male    |         1 |  0.054 |     0.024 |
| Kidney stones              | taxonomy_family | age_bin    | 40-64   |         1 |  0.061 |     0.056 |
| Kidney stones              | taxonomy_family | sex        | male    |         1 |  0.055 |     0.05  |
| Kidney stones              | taxonomy_genus  | age_bin    | 40-64   |         1 |  0.059 |     0.046 |
| Kidney stones              | taxonomy_genus  | sex        | male    |         1 |  0.048 |     0.044 |
| Parkinson's disease        | ko_eggnog       | age_bin    | 40-64   |         1 |  0.038 |     0.031 |
| Parkinson's disease        | ko_eggnog       | age_bin    | 65+     |         1 |  0.044 |     0.035 |
| Parkinson's disease        | ko_eggnog       | sex        | female  |         1 |  0.053 |     0.033 |
| Parkinson's disease        | ko_eggnog       | sex        | male    |         1 |  0.037 |     0.031 |
| Parkinson's disease        | ko_kofam        | age_bin    | 40-64   |         1 |  0.036 |     0.032 |
| Parkinson's disease        | ko_kofam        | age_bin    | 65+     |         1 |  0.046 |     0.036 |
| Parkinson's disease        | ko_kofam        | sex        | female  |         1 |  0.051 |     0.034 |
| Parkinson's disease        | ko_kofam        | sex        | male    |         1 |  0.039 |     0.035 |
| Parkinson's disease        | module          | age_bin    | 40-64   |         1 |  0     |     0     |
| Parkinson's disease        | module          | age_bin    | 65+     |         1 |  0     |     0     |
| Parkinson's disease        | module          | sex        | female  |         1 |  0     |     0     |
| Parkinson's disease        | module          | sex        | male    |         1 |  0     |     0     |
| Parkinson's disease        | pfam            | age_bin    | 40-64   |         1 |  0.041 |     0.034 |
| Parkinson's disease        | pfam            | age_bin    | 65+     |         1 |  0.048 |     0.035 |
| Parkinson's disease        | pfam            | sex        | female  |         1 |  0.058 |     0.035 |
| Parkinson's disease        | pfam            | sex        | male    |         1 |  0.042 |     0.035 |
| Parkinson's disease        | taxonomy_family | age_bin    | 40-64   |         1 |  0.05  |     0.04  |
| Parkinson's disease        | taxonomy_family | age_bin    | 65+     |         1 |  0.054 |     0.048 |
| Parkinson's disease        | taxonomy_family | sex        | female  |         1 |  0.061 |     0.049 |
| Parkinson's disease        | taxonomy_family | sex        | male    |         1 |  0.05  |     0.044 |
| Parkinson's disease        | taxonomy_genus  | age_bin    | 40-64   |         1 |  0.062 |     0.056 |
| Parkinson's disease        | taxonomy_genus  | age_bin    | 65+     |         1 |  0.062 |     0.057 |
| Parkinson's disease        | taxonomy_genus  | sex        | female  |         1 |  0.069 |     0.057 |
| Parkinson's disease        | taxonomy_genus  | sex        | male    |         1 |  0.058 |     0.055 |
| Rheumatoid arthritis       | ko_eggnog       | sex        | female  |         1 |  0.032 |     0.04  |
| Rheumatoid arthritis       | ko_eggnog       | sex        | male    |         1 |  0.037 |     0.046 |
| Rheumatoid arthritis       | ko_kofam        | sex        | female  |         1 |  0.032 |     0.038 |
| Rheumatoid arthritis       | ko_kofam        | sex        | male    |         1 |  0.038 |     0.047 |
| Rheumatoid arthritis       | module          | sex        | female  |         1 |  0     |     0     |
| Rheumatoid arthritis       | module          | sex        | male    |         1 |  0     |     0     |
| Rheumatoid arthritis       | pfam            | sex        | female  |         1 |  0.03  |     0.035 |
| Rheumatoid arthritis       | pfam            | sex        | male    |         1 |  0.036 |     0.043 |
| Rheumatoid arthritis       | taxonomy_family | sex        | female  |         1 |  0.052 |     0.061 |
| Rheumatoid arthritis       | taxonomy_family | sex        | male    |         1 |  0.055 |     0.057 |
| Rheumatoid arthritis       | taxonomy_genus  | sex        | female  |         1 |  0.055 |     0.06  |
| Rheumatoid arthritis       | taxonomy_genus  | sex        | male    |         1 |  0.065 |     0.064 |
| Seasonal allergies         | ko_eggnog       | age_bin    | 18-39   |         1 |  0.02  |     0.018 |
| Seasonal allergies         | ko_eggnog       | age_bin    | 40-64   |         1 |  0.024 |     0.023 |
| Seasonal allergies         | ko_eggnog       | age_bin    | 65+     |         1 |  0.017 |     0.019 |
| Seasonal allergies         | ko_eggnog       | sex        | female  |         1 |  0.022 |     0.019 |
| Seasonal allergies         | ko_eggnog       | sex        | male    |         1 |  0.022 |     0.02  |
| Seasonal allergies         | ko_kofam        | age_bin    | 18-39   |         1 |  0.021 |     0.022 |
| Seasonal allergies         | ko_kofam        | age_bin    | 40-64   |         1 |  0.027 |     0.021 |
| Seasonal allergies         | ko_kofam        | age_bin    | 65+     |         1 |  0.019 |     0.017 |
| Seasonal allergies         | ko_kofam        | sex        | female  |         1 |  0.023 |     0.02  |
| Seasonal allergies         | ko_kofam        | sex        | male    |         1 |  0.022 |     0.022 |
| Seasonal allergies         | module          | age_bin    | 18-39   |         1 |  0     |     0     |
| Seasonal allergies         | module          | age_bin    | 40-64   |         1 |  0     |     0.005 |
| Seasonal allergies         | module          | age_bin    | 65+     |         1 |  0     |     0     |
| Seasonal allergies         | module          | sex        | female  |         1 |  0     |     0     |
| Seasonal allergies         | module          | sex        | male    |         1 |  0     |     0     |
| Seasonal allergies         | pfam            | age_bin    | 18-39   |         1 |  0.019 |     0.019 |
| Seasonal allergies         | pfam            | age_bin    | 40-64   |         1 |  0.026 |     0.02  |
| Seasonal allergies         | pfam            | age_bin    | 65+     |         1 |  0.018 |     0.02  |
| Seasonal allergies         | pfam            | sex        | female  |         1 |  0.023 |     0.019 |
| Seasonal allergies         | pfam            | sex        | male    |         1 |  0.021 |     0.02  |
| Seasonal allergies         | taxonomy_family | age_bin    | 18-39   |         1 |  0.031 |     0.03  |
| Seasonal allergies         | taxonomy_family | age_bin    | 40-64   |         1 |  0.036 |     0.045 |
| Seasonal allergies         | taxonomy_family | age_bin    | 65+     |         1 |  0.037 |     0.044 |
| Seasonal allergies         | taxonomy_family | sex        | female  |         1 |  0.034 |     0.034 |
| Seasonal allergies         | taxonomy_family | sex        | male    |         1 |  0.037 |     0.04  |
| Seasonal allergies         | taxonomy_genus  | age_bin    | 18-39   |         1 |  0.03  |     0.028 |
| Seasonal allergies         | taxonomy_genus  | age_bin    | 40-64   |         1 |  0.033 |     0.038 |
| Seasonal allergies         | taxonomy_genus  | age_bin    | 65+     |         1 |  0.041 |     0.04  |
| Seasonal allergies         | taxonomy_genus  | sex        | female  |         1 |  0.033 |     0.03  |
| Seasonal allergies         | taxonomy_genus  | sex        | male    |         1 |  0.038 |     0.038 |

Wall 132s.