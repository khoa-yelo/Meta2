# S14 — disease deviations by age, sex, BMI and location, pipeline B — 2026-09-30 00:27 UTC

15 case/control studies (>= 15 cases and >= 15 controls), 1,437 samples. A study x level cell is analysed when it holds >= 10 cases and >= 10 controls: 63 cells in 15 studies.

## 1. Balance of cases and controls

| study                  | condition    | region       |   n_case |   n_control |   age_known |   age_case |   age_control |   age_p |   bmi_known |   bmi_case |   bmi_control |   bmi_p |   sex_known |   female_case |   female_control |   sex_p |
|:-----------------------|:-------------|:-------------|---------:|------------:|------------:|-----------:|--------------:|--------:|------------:|-----------:|--------------:|--------:|------------:|--------------:|-----------------:|--------:|
| FengQ_2015             | CRC          | Europe       |       46 |          61 |           1 |       68.5 |          68   |   0.616 |           1 |     26.895 |        29     |   0.122 |           1 |         0.391 |            0.41  |   1     |
| GuptaA_2019            | CRC          | S. Asia      |       30 |          30 |           1 |       59   |          40   |   0     |           1 |     20.08  |        22.03  |   0.008 |           1 |         0.4   |            0.6   |   0.197 |
| HMP_2019               | IBD          | N. America   |       62 |          20 |           1 |       16.5 |          16.5 |   0.983 |           1 |     21.25  |        23.25  |   0.479 |           1 |         0.5   |            0.5   |   1     |
| HanniganGD_2017        | CRC          | N. America   |       26 |          28 |           1 |       59   |          54   |   0.228 |           1 |     28.719 |        25.18  |   0.05  |           1 |         0.269 |            0.607 |   0.026 |
| KarlssonFH_2013        | T2D          | Europe/other |       53 |          43 |           1 |       70   |          70   |   0.564 |           1 |     28     |        25.7   |   0.006 |           1 |         1     |            1     | nan     |
| NagySzakalD_2017       | ME/CFS       | N. America   |       50 |          50 |           1 |       51   |          51   | nan     |           1 |     25.745 |        24.388 |   0.171 |           1 |         0.82  |            0.82  |   1     |
| RubelMA_2020           | STH          | Africa       |       77 |          82 |           1 |       37   |          45   |   0.262 |           1 |     21.6   |        23.3   |   0.001 |           1 |         0.532 |            0.61  |   0.41  |
| SankaranarayananK_2015 | T2D          | N. America   |       19 |          18 |           1 |       55   |          47   |   0.049 |           1 |     35.1   |        32.5   |   0.362 |           1 |         0.579 |            0.667 |   0.833 |
| ThomasAM_2019_a        | CRC          | Europe       |       28 |          22 |           1 |       70.5 |          67.5 |   0.159 |           1 |     24.5   |        25     |   0.992 |           1 |         0.214 |            0.409 |   0.238 |
| ThomasAM_2019_b        | CRC          | Europe       |       31 |          27 |           1 |       59   |          58   |   0.864 |           1 |     27.441 |        24.163 |   0.029 |           1 |         0.258 |            0.407 |   0.353 |
| VogtmannE_2016         | CRC          | N. America   |       49 |          52 |           1 |       64   |          63   |   0.962 |           1 |     24.325 |        24.014 |   0.62  |           1 |         0.265 |            0.288 |   0.97  |
| WirbelJ_2018           | CRC          | Europe       |       60 |          65 |           1 |       63   |          58   |   0.002 |           1 |     26.05  |        24.6   |   0.045 |           1 |         0.4   |            0.431 |   0.867 |
| YuJ_2015               | CRC          | E. Asia      |       74 |          53 |           1 |       67   |          64   |   0.012 |           1 |     23.85  |        23     |   0.152 |           1 |         0.365 |            0.377 |   1     |
| ZellerG_2014           | CRC          | Europe       |       51 |          59 |           1 |       67   |          63   |   0.009 |           1 |     25     |        24     |   0.487 |           1 |         0.451 |            0.542 |   0.444 |
| ZhuF_2020              | schizofrenia | E. Asia      |       90 |          81 |           1 |       28   |          27   |   0.059 |           1 |     20.28  |        21.6   |   0.017 |           1 |         0.489 |            0.494 |   1     |

## 2. Cells analysed

|                                   |   ('age_bin', '18-39') |   ('age_bin', '40-64') |   ('age_bin', '65+') |   ('bmi_class', '18.5-25') |   ('bmi_class', '25-30') |   ('bmi_class', '<18.5') |   ('bmi_class', '>=30') |   ('sex', 'female') |   ('sex', 'male') |
|:----------------------------------|-----------------------:|-----------------------:|---------------------:|---------------------------:|-------------------------:|-------------------------:|------------------------:|--------------------:|------------------:|
| ('FengQ_2015', 'CRC')             |                      0 |                     17 |                   29 |                         15 |                       24 |                        0 |                       0 |                  18 |                28 |
| ('GuptaA_2019', 'CRC')            |                      0 |                     21 |                    0 |                         26 |                        0 |                        0 |                       0 |                  12 |                18 |
| ('HMP_2019', 'IBD')               |                      0 |                      0 |                    0 |                          0 |                        0 |                        0 |                       0 |                  29 |                31 |
| ('HanniganGD_2017', 'CRC')        |                      0 |                     18 |                    0 |                          0 |                        0 |                        0 |                       0 |                   0 |                18 |
| ('KarlssonFH_2013', 'T2D')        |                      0 |                      0 |                   53 |                         13 |                       23 |                        0 |                       0 |                  53 |                 0 |
| ('NagySzakalD_2017', 'ME/CFS')    |                      0 |                     50 |                    0 |                         21 |                       15 |                        0 |                       0 |                  41 |                 0 |
| ('RubelMA_2020', 'STH')           |                     43 |                     26 |                    0 |                         63 |                        0 |                        0 |                       0 |                  41 |                36 |
| ('SankaranarayananK_2015', 'T2D') |                      0 |                      0 |                    0 |                          0 |                        0 |                        0 |                      14 |                  11 |                 0 |
| ('ThomasAM_2019_a', 'CRC')        |                      0 |                      0 |                   20 |                         14 |                        0 |                        0 |                       0 |                   0 |                22 |
| ('ThomasAM_2019_b', 'CRC')        |                      0 |                     21 |                    0 |                          0 |                        0 |                        0 |                       0 |                   0 |                23 |
| ('VogtmannE_2016', 'CRC')         |                      0 |                     22 |                   24 |                         25 |                       14 |                        0 |                       0 |                  13 |                36 |
| ('WirbelJ_2018', 'CRC')           |                      0 |                     30 |                   28 |                         20 |                       26 |                        0 |                       0 |                  24 |                36 |
| ('YuJ_2015', 'CRC')               |                      0 |                     34 |                   39 |                         44 |                        0 |                        0 |                       0 |                  27 |                47 |
| ('ZellerG_2014', 'CRC')           |                      0 |                     22 |                   29 |                         22 |                       17 |                        0 |                       0 |                  23 |                28 |
| ('ZhuF_2020', 'schizofrenia')     |                     68 |                     13 |                    0 |                         67 |                        0 |                       18 |                       0 |                  44 |                46 |

(number of cases in each analysed cell; 0 = cell not analysed)

## 3. Shifts adjusted for the variable vs crude shifts (atlas-significant features: q <= 0.05 and |shift| >= 10)

| study                  | condition    | variable   |   levels_used |   features |   spearman_crude_vs_adjusted |   atlas_significant |   same_sign |   retained_half_magnitude |   median_ratio_adjusted_to_crude |
|:-----------------------|:-------------|:-----------|--------------:|-----------:|-----------------------------:|--------------------:|------------:|--------------------------:|---------------------------------:|
| FengQ_2015             | CRC          | age_bin    |             2 |       4520 |                        0.883 |                 132 |       1     |                     1     |                            0.978 |
| FengQ_2015             | CRC          | bmi_class  |             2 |       3466 |                        0.647 |                 126 |       1     |                     0.944 |                            0.948 |
| FengQ_2015             | CRC          | sex        |             2 |       4232 |                        0.887 |                 132 |       1     |                     1     |                            0.975 |
| GuptaA_2019            | CRC          | age_bin    |             1 |       1308 |                        0.919 |                 317 |       0.997 |                     0.943 |                            0.988 |
| GuptaA_2019            | CRC          | bmi_class  |             1 |       2207 |                        0.951 |                 433 |       1     |                     0.998 |                            0.999 |
| GuptaA_2019            | CRC          | sex        |             2 |       1349 |                        0.94  |                 288 |       1     |                     0.993 |                            0.953 |
| HMP_2019               | IBD          | sex        |             2 |         25 |                        0.647 |                   0 |     nan     |                   nan     |                          nan     |
| HanniganGD_2017        | CRC          | age_bin    |             1 |       1691 |                        0.765 |                   0 |     nan     |                   nan     |                          nan     |
| HanniganGD_2017        | CRC          | sex        |             1 |       1076 |                        0.631 |                   0 |     nan     |                   nan     |                          nan     |
| KarlssonFH_2013        | T2D          | age_bin    |             1 |       4313 |                        1     |                   2 |       1     |                     1     |                            1     |
| KarlssonFH_2013        | T2D          | bmi_class  |             2 |       2707 |                        0.67  |                   2 |       1     |                     1     |                            1.047 |
| KarlssonFH_2013        | T2D          | sex        |             1 |       4313 |                        1     |                   2 |       1     |                     1     |                            1     |
| NagySzakalD_2017       | ME/CFS       | age_bin    |             1 |       4212 |                        1     |                   8 |       1     |                     1     |                            1     |
| NagySzakalD_2017       | ME/CFS       | bmi_class  |             2 |       2591 |                        0.696 |                   8 |       1     |                     0.875 |                            1.015 |
| NagySzakalD_2017       | ME/CFS       | sex        |             1 |       3937 |                        0.881 |                   8 |       1     |                     1     |                            0.965 |
| RubelMA_2020           | STH          | age_bin    |             2 |       3817 |                        0.926 |                 748 |       1     |                     0.984 |                            1.047 |
| RubelMA_2020           | STH          | bmi_class  |             1 |       3878 |                        0.9   |                 738 |       0.999 |                     0.988 |                            1.074 |
| RubelMA_2020           | STH          | sex        |             2 |       4011 |                        0.945 |                 757 |       1     |                     0.995 |                            0.994 |
| SankaranarayananK_2015 | T2D          | bmi_class  |             1 |       2119 |                        0.866 |                   1 |       1     |                     1     |                            1.006 |
| SankaranarayananK_2015 | T2D          | sex        |             1 |       1856 |                        0.687 |                   2 |       1     |                     1     |                            0.796 |
| ThomasAM_2019_a        | CRC          | age_bin    |             1 |       2097 |                        0.825 |                 126 |       1     |                     0.968 |                            0.831 |
| ThomasAM_2019_a        | CRC          | bmi_class  |             1 |       1239 |                        0.668 |                 101 |       0.97  |                     0.772 |                            0.749 |
| ThomasAM_2019_a        | CRC          | sex        |             1 |       2114 |                        0.854 |                 111 |       1     |                     0.982 |                            0.903 |
| ThomasAM_2019_b        | CRC          | age_bin    |             1 |       2752 |                        0.771 |                  10 |       1     |                     0.9   |                            0.856 |
| ThomasAM_2019_b        | CRC          | sex        |             1 |       2658 |                        0.793 |                  10 |       1     |                     1     |                            1.18  |
| VogtmannE_2016         | CRC          | age_bin    |             2 |       3180 |                        0.823 |                   0 |     nan     |                   nan     |                          nan     |
| VogtmannE_2016         | CRC          | bmi_class  |             2 |       3411 |                        0.78  |                   2 |       1     |                     1     |                            0.969 |
| VogtmannE_2016         | CRC          | sex        |             2 |       3634 |                        0.869 |                   2 |       1     |                     1     |                            0.726 |
| WirbelJ_2018           | CRC          | age_bin    |             2 |       3104 |                        0.819 |                  19 |       1     |                     0.947 |                            0.844 |
| WirbelJ_2018           | CRC          | bmi_class  |             2 |       2985 |                        0.828 |                  19 |       1     |                     1     |                            0.819 |
| WirbelJ_2018           | CRC          | sex        |             2 |       3397 |                        0.902 |                  20 |       1     |                     1     |                            1.005 |
| YuJ_2015               | CRC          | age_bin    |             2 |       4573 |                        0.832 |                  15 |       1     |                     1     |                            0.98  |
| YuJ_2015               | CRC          | bmi_class  |             1 |       4674 |                        0.764 |                  16 |       1     |                     0.938 |                            1.065 |
| YuJ_2015               | CRC          | sex        |             2 |       4566 |                        0.861 |                  16 |       1     |                     1     |                            0.983 |
| ZellerG_2014           | CRC          | age_bin    |             2 |       4330 |                        0.894 |                 232 |       1     |                     0.996 |                            0.957 |
| ZellerG_2014           | CRC          | bmi_class  |             2 |       4082 |                        0.872 |                 229 |       1     |                     1     |                            0.989 |
| ZellerG_2014           | CRC          | sex        |             2 |       4273 |                        0.904 |                 232 |       1     |                     1     |                            1.005 |
| ZhuF_2020              | schizofrenia | age_bin    |             2 |       4664 |                        0.855 |                   2 |       1     |                     1     |                            0.842 |
| ZhuF_2020              | schizofrenia | bmi_class  |             2 |       4723 |                        0.85  |                   2 |       1     |                     1     |                            0.851 |
| ZhuF_2020              | schizofrenia | sex        |             2 |       4628 |                        0.913 |                   2 |       1     |                     1     |                            1     |

## 4. Agreement of shifts between levels within a study (medians over studies and layers)

| condition    | variable   | level_a   | level_b   |   studies |   spearman_all |   same_sign_significant |
|:-------------|:-----------|:----------|:----------|----------:|---------------:|------------------------:|
| CRC          | age_bin    | 40-64     | 65+       |         5 |          0.193 |                       1 |
| CRC          | bmi_class  | 18.5-25   | 25-30     |         4 |          0.166 |                       1 |
| CRC          | sex        | female    | male      |         6 |          0.332 |                       1 |
| ME/CFS       | bmi_class  | 18.5-25   | 25-30     |         1 |          0.131 |                       1 |
| STH          | age_bin    | 18-39     | 40-64     |         1 |          0.324 |                       1 |
| STH          | sex        | female    | male      |         1 |          0.458 |                       1 |
| T2D          | bmi_class  | 18.5-25   | 25-30     |         1 |          0.05  |                       1 |
| schizofrenia | age_bin    | 18-39     | 40-64     |         1 |         -0.002 |                       1 |
| schizofrenia | bmi_class  | <18.5     | 18.5-25   |         1 |          0.02  |                       1 |
| schizofrenia | sex        | female    | male      |         1 |         -0.104 |                       1 |

## 5. Features consistent in the atlas: do they hold within each level and region?

| condition   | variable   | level      |   features |   studies |   share_same_direction |   median_abs_shift |
|:------------|:-----------|:-----------|-----------:|----------:|-----------------------:|-------------------:|
| CRC         | age_bin    | 40-64      |        203 |         8 |                  0.966 |             11.501 |
| CRC         | age_bin    | 65+        |        203 |         6 |                  0.946 |             13.09  |
| CRC         | bmi_class  | 18.5-25    |        203 |         7 |                  0.995 |             13.704 |
| CRC         | bmi_class  | 25-30      |        202 |         4 |                  0.95  |              8.594 |
| CRC         | region     | E. Asia    |        202 |         1 |                  0.738 |              5.142 |
| CRC         | region     | Europe     |        203 |         5 |                  0.995 |             14.074 |
| CRC         | region     | N. America |        203 |         2 |                  0.818 |              5.151 |
| CRC         | region     | S. Asia    |        191 |         1 |                  0.969 |             25.468 |
| CRC         | sex        | female     |        202 |         6 |                  0.941 |             11.201 |
| CRC         | sex        | male       |        203 |         9 |                  0.995 |             11.937 |

(features = atlas-consistent features assessable in the level; share_same_direction = share whose median shift over studies has the atlas direction)

## 6. Sample level: median fraction outside, cases vs controls, by level (medians over studies)

| condition    | layer            | variable   | level   |   studies |   case |   control |
|:-------------|:-----------------|:-----------|:--------|----------:|-------:|----------:|
| CRC          | ko_humann        | age_bin    | 40-64   |         8 |  0.038 |     0.038 |
| CRC          | ko_humann        | age_bin    | 65+     |         6 |  0.041 |     0.034 |
| CRC          | ko_humann        | bmi_class  | 18.5-25 |         7 |  0.045 |     0.04  |
| CRC          | ko_humann        | bmi_class  | 25-30   |         4 |  0.042 |     0.053 |
| CRC          | ko_humann        | sex        | female  |         6 |  0.039 |     0.034 |
| CRC          | ko_humann        | sex        | male    |         9 |  0.041 |     0.035 |
| CRC          | pathway_humann   | age_bin    | 40-64   |         8 |  0.024 |     0.027 |
| CRC          | pathway_humann   | age_bin    | 65+     |         6 |  0.024 |     0.026 |
| CRC          | pathway_humann   | bmi_class  | 18.5-25 |         7 |  0.029 |     0.027 |
| CRC          | pathway_humann   | bmi_class  | 25-30   |         4 |  0.022 |     0.038 |
| CRC          | pathway_humann   | sex        | female  |         6 |  0.022 |     0.027 |
| CRC          | pathway_humann   | sex        | male    |         9 |  0.032 |     0.026 |
| CRC          | taxonomy_family  | age_bin    | 40-64   |         8 |  0.059 |     0.051 |
| CRC          | taxonomy_family  | age_bin    | 65+     |         6 |  0.069 |     0.058 |
| CRC          | taxonomy_family  | bmi_class  | 18.5-25 |         7 |  0.061 |     0.065 |
| CRC          | taxonomy_family  | bmi_class  | 25-30   |         4 |  0.048 |     0.04  |
| CRC          | taxonomy_family  | sex        | female  |         6 |  0.078 |     0.048 |
| CRC          | taxonomy_family  | sex        | male    |         9 |  0.05  |     0.065 |
| CRC          | taxonomy_genus   | age_bin    | 40-64   |         8 |  0.058 |     0.06  |
| CRC          | taxonomy_genus   | age_bin    | 65+     |         6 |  0.064 |     0.066 |
| CRC          | taxonomy_genus   | bmi_class  | 18.5-25 |         7 |  0.059 |     0.066 |
| CRC          | taxonomy_genus   | bmi_class  | 25-30   |         4 |  0.066 |     0.054 |
| CRC          | taxonomy_genus   | sex        | female  |         6 |  0.074 |     0.061 |
| CRC          | taxonomy_genus   | sex        | male    |         9 |  0.058 |     0.065 |
| CRC          | taxonomy_species | age_bin    | 40-64   |         8 |  0.061 |     0.053 |
| CRC          | taxonomy_species | age_bin    | 65+     |         6 |  0.074 |     0.066 |
| CRC          | taxonomy_species | bmi_class  | 18.5-25 |         7 |  0.071 |     0.064 |
| CRC          | taxonomy_species | bmi_class  | 25-30   |         4 |  0.073 |     0.073 |
| CRC          | taxonomy_species | sex        | female  |         6 |  0.069 |     0.058 |
| CRC          | taxonomy_species | sex        | male    |         9 |  0.065 |     0.054 |
| IBD          | taxonomy_family  | sex        | female  |         1 |  0.077 |     0.089 |
| IBD          | taxonomy_family  | sex        | male    |         1 |  0.062 |     0.109 |
| IBD          | taxonomy_genus   | sex        | female  |         1 |  0.093 |     0.113 |
| IBD          | taxonomy_genus   | sex        | male    |         1 |  0.086 |     0.074 |
| IBD          | taxonomy_species | sex        | female  |         1 |  0.075 |     0.082 |
| IBD          | taxonomy_species | sex        | male    |         1 |  0.07  |     0.06  |
| ME/CFS       | ko_humann        | age_bin    | 40-64   |         1 |  0.064 |     0.045 |
| ME/CFS       | ko_humann        | bmi_class  | 18.5-25 |         1 |  0.062 |     0.049 |
| ME/CFS       | ko_humann        | bmi_class  | 25-30   |         1 |  0.064 |     0.06  |
| ME/CFS       | ko_humann        | sex        | female  |         1 |  0.068 |     0.043 |
| ME/CFS       | pathway_humann   | age_bin    | 40-64   |         1 |  0.077 |     0.044 |
| ME/CFS       | pathway_humann   | bmi_class  | 18.5-25 |         1 |  0.086 |     0.057 |
| ME/CFS       | pathway_humann   | bmi_class  | 25-30   |         1 |  0.047 |     0.045 |
| ME/CFS       | pathway_humann   | sex        | female  |         1 |  0.086 |     0.035 |
| ME/CFS       | taxonomy_family  | age_bin    | 40-64   |         1 |  0.083 |     0.049 |
| ME/CFS       | taxonomy_family  | bmi_class  | 18.5-25 |         1 |  0.05  |     0.056 |
| ME/CFS       | taxonomy_family  | bmi_class  | 25-30   |         1 |  0.08  |     0.04  |
| ME/CFS       | taxonomy_family  | sex        | female  |         1 |  0.087 |     0.045 |
| ME/CFS       | taxonomy_genus   | age_bin    | 40-64   |         1 |  0.096 |     0.064 |
| ME/CFS       | taxonomy_genus   | bmi_class  | 18.5-25 |         1 |  0.117 |     0.066 |
| ME/CFS       | taxonomy_genus   | bmi_class  | 25-30   |         1 |  0.093 |     0.042 |
| ME/CFS       | taxonomy_genus   | sex        | female  |         1 |  0.105 |     0.058 |
| ME/CFS       | taxonomy_species | age_bin    | 40-64   |         1 |  0.087 |     0.079 |
| ME/CFS       | taxonomy_species | bmi_class  | 18.5-25 |         1 |  0.083 |     0.08  |
| ME/CFS       | taxonomy_species | bmi_class  | 25-30   |         1 |  0.091 |     0.07  |
| ME/CFS       | taxonomy_species | sex        | female  |         1 |  0.091 |     0.076 |
| STH          | ko_humann        | age_bin    | 18-39   |         1 |  0.044 |     0.025 |
| STH          | ko_humann        | age_bin    | 40-64   |         1 |  0.044 |     0.032 |
| STH          | ko_humann        | bmi_class  | 18.5-25 |         1 |  0.047 |     0.025 |
| STH          | ko_humann        | sex        | female  |         1 |  0.035 |     0.028 |
| STH          | ko_humann        | sex        | male    |         1 |  0.049 |     0.032 |
| STH          | pathway_humann   | age_bin    | 18-39   |         1 |  0.023 |     0.026 |
| STH          | pathway_humann   | age_bin    | 40-64   |         1 |  0.031 |     0.03  |
| STH          | pathway_humann   | bmi_class  | 18.5-25 |         1 |  0.026 |     0.027 |
| STH          | pathway_humann   | sex        | female  |         1 |  0.026 |     0.023 |
| STH          | pathway_humann   | sex        | male    |         1 |  0.026 |     0.037 |
| STH          | taxonomy_family  | age_bin    | 18-39   |         1 |  0     |     0     |
| STH          | taxonomy_family  | age_bin    | 40-64   |         1 |  0.057 |     0.043 |
| STH          | taxonomy_family  | bmi_class  | 18.5-25 |         1 |  0     |     0.038 |
| STH          | taxonomy_family  | sex        | female  |         1 |  0     |     0.043 |
| STH          | taxonomy_family  | sex        | male    |         1 |  0     |     0     |
| STH          | taxonomy_genus   | age_bin    | 18-39   |         1 |  0.054 |     0.029 |
| STH          | taxonomy_genus   | age_bin    | 40-64   |         1 |  0.055 |     0.051 |
| STH          | taxonomy_genus   | bmi_class  | 18.5-25 |         1 |  0.053 |     0.05  |
| STH          | taxonomy_genus   | sex        | female  |         1 |  0.05  |     0.047 |
| STH          | taxonomy_genus   | sex        | male    |         1 |  0.056 |     0.053 |
| STH          | taxonomy_species | age_bin    | 18-39   |         1 |  0.042 |     0.041 |
| STH          | taxonomy_species | age_bin    | 40-64   |         1 |  0.044 |     0.039 |
| STH          | taxonomy_species | bmi_class  | 18.5-25 |         1 |  0.043 |     0.044 |
| STH          | taxonomy_species | sex        | female  |         1 |  0.038 |     0.039 |
| STH          | taxonomy_species | sex        | male    |         1 |  0.045 |     0.05  |
| T2D          | ko_humann        | age_bin    | 65+     |         1 |  0.022 |     0.018 |
| T2D          | ko_humann        | bmi_class  | 18.5-25 |         1 |  0.016 |     0.018 |
| T2D          | ko_humann        | bmi_class  | 25-30   |         1 |  0.039 |     0.019 |
| T2D          | ko_humann        | bmi_class  | >=30    |         1 |  0.04  |     0.041 |
| T2D          | ko_humann        | sex        | female  |         2 |  0.027 |     0.026 |
| T2D          | pathway_humann   | age_bin    | 65+     |         1 |  0.019 |     0.017 |
| T2D          | pathway_humann   | bmi_class  | 18.5-25 |         1 |  0.015 |     0.017 |
| T2D          | pathway_humann   | bmi_class  | 25-30   |         1 |  0.022 |     0.011 |
| T2D          | pathway_humann   | bmi_class  | >=30    |         1 |  0.032 |     0.021 |
| T2D          | pathway_humann   | sex        | female  |         2 |  0.019 |     0.02  |
| T2D          | taxonomy_family  | age_bin    | 65+     |         1 |  0.036 |     0     |
| T2D          | taxonomy_family  | bmi_class  | 18.5-25 |         1 |  0.043 |     0     |
| T2D          | taxonomy_family  | bmi_class  | 25-30   |         1 |  0.036 |     0     |
| T2D          | taxonomy_family  | bmi_class  | >=30    |         1 |  0.047 |     0.08  |
| T2D          | taxonomy_family  | sex        | female  |         2 |  0.036 |     0.023 |
| T2D          | taxonomy_genus   | age_bin    | 65+     |         1 |  0.033 |     0.021 |
| T2D          | taxonomy_genus   | bmi_class  | 18.5-25 |         1 |  0.024 |     0.027 |
| T2D          | taxonomy_genus   | bmi_class  | 25-30   |         1 |  0.048 |     0.021 |
| T2D          | taxonomy_genus   | bmi_class  | >=30    |         1 |  0.045 |     0.071 |
| T2D          | taxonomy_genus   | sex        | female  |         2 |  0.039 |     0.036 |
| T2D          | taxonomy_species | age_bin    | 65+     |         1 |  0.045 |     0.029 |
| T2D          | taxonomy_species | bmi_class  | 18.5-25 |         1 |  0.052 |     0.021 |
| T2D          | taxonomy_species | bmi_class  | 25-30   |         1 |  0.05  |     0.03  |
| T2D          | taxonomy_species | bmi_class  | >=30    |         1 |  0.06  |     0.05  |
| T2D          | taxonomy_species | sex        | female  |         2 |  0.047 |     0.037 |
| schizofrenia | ko_humann        | age_bin    | 18-39   |         1 |  0.054 |     0.048 |
| schizofrenia | ko_humann        | age_bin    | 40-64   |         1 |  0.051 |     0.07  |
| schizofrenia | ko_humann        | bmi_class  | 18.5-25 |         1 |  0.056 |     0.052 |
| schizofrenia | ko_humann        | bmi_class  | <18.5   |         1 |  0.053 |     0.044 |
| schizofrenia | ko_humann        | sex        | female  |         1 |  0.054 |     0.047 |
| schizofrenia | ko_humann        | sex        | male    |         1 |  0.054 |     0.061 |
| schizofrenia | pathway_humann   | age_bin    | 18-39   |         1 |  0.03  |     0.041 |
| schizofrenia | pathway_humann   | age_bin    | 40-64   |         1 |  0.043 |     0.034 |
| schizofrenia | pathway_humann   | bmi_class  | 18.5-25 |         1 |  0.033 |     0.041 |
| schizofrenia | pathway_humann   | bmi_class  | <18.5   |         1 |  0.048 |     0.023 |
| schizofrenia | pathway_humann   | sex        | female  |         1 |  0.033 |     0.027 |
| schizofrenia | pathway_humann   | sex        | male    |         1 |  0.034 |     0.049 |
| schizofrenia | taxonomy_family  | age_bin    | 18-39   |         1 |  0.043 |     0.077 |
| schizofrenia | taxonomy_family  | age_bin    | 40-64   |         1 |  0.069 |     0.044 |
| schizofrenia | taxonomy_family  | bmi_class  | 18.5-25 |         1 |  0.045 |     0.08  |
| schizofrenia | taxonomy_family  | bmi_class  | <18.5   |         1 |  0.05  |     0.034 |
| schizofrenia | taxonomy_family  | sex        | female  |         1 |  0.04  |     0.04  |
| schizofrenia | taxonomy_family  | sex        | male    |         1 |  0.064 |     0.083 |
| schizofrenia | taxonomy_genus   | age_bin    | 18-39   |         1 |  0.063 |     0.057 |
| schizofrenia | taxonomy_genus   | age_bin    | 40-64   |         1 |  0.043 |     0.051 |
| schizofrenia | taxonomy_genus   | bmi_class  | 18.5-25 |         1 |  0.059 |     0.059 |
| schizofrenia | taxonomy_genus   | bmi_class  | <18.5   |         1 |  0.061 |     0.04  |
| schizofrenia | taxonomy_genus   | sex        | female  |         1 |  0.062 |     0.047 |
| schizofrenia | taxonomy_genus   | sex        | male    |         1 |  0.055 |     0.062 |
| schizofrenia | taxonomy_species | age_bin    | 18-39   |         1 |  0.067 |     0.053 |
| schizofrenia | taxonomy_species | age_bin    | 40-64   |         1 |  0.049 |     0.064 |
| schizofrenia | taxonomy_species | bmi_class  | 18.5-25 |         1 |  0.067 |     0.056 |
| schizofrenia | taxonomy_species | bmi_class  | <18.5   |         1 |  0.057 |     0.06  |
| schizofrenia | taxonomy_species | sex        | female  |         1 |  0.061 |     0.055 |
| schizofrenia | taxonomy_species | sex        | male    |         1 |  0.067 |     0.067 |

Wall 72s.