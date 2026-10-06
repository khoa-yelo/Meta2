# S12.3 — disease deviations relative to the reference, pipeline B

15 case/control studies with >= 15 cases and >= 15 controls.

## Sample-level excess test: share flagged (q <= 0.05) and median fraction outside, cases vs controls

| study                  | condition    | layer            |   excess_positive_control |   excess_positive_case |   median_frac_outside_control |   median_frac_outside_case |
|:-----------------------|:-------------|:-----------------|--------------------------:|-----------------------:|------------------------------:|---------------------------:|
| FengQ_2015             | CRC          | ko_humann        |                     0.55  |                  0.478 |                         0.072 |                      0.053 |
| FengQ_2015             | CRC          | pathway_humann   |                     0.433 |                  0.304 |                         0.065 |                      0.034 |
| FengQ_2015             | CRC          | taxonomy_family  |                     0.197 |                  0.043 |                         0.094 |                      0.085 |
| FengQ_2015             | CRC          | taxonomy_genus   |                     0.459 |                  0.152 |                         0.127 |                      0.075 |
| FengQ_2015             | CRC          | taxonomy_species |                     0.508 |                  0.217 |                         0.125 |                      0.08  |
| GuptaA_2019            | CRC          | ko_humann        |                     0.633 |                  0.533 |                         0.073 |                      0.064 |
| GuptaA_2019            | CRC          | pathway_humann   |                     0.5   |                  0.367 |                         0.073 |                      0.054 |
| GuptaA_2019            | CRC          | taxonomy_family  |                     0.033 |                  0.033 |                         0.08  |                      0.095 |
| GuptaA_2019            | CRC          | taxonomy_genus   |                     0     |                  0     |                         0.077 |                      0.073 |
| GuptaA_2019            | CRC          | taxonomy_species |                     0     |                  0     |                         0.032 |                      0.064 |
| HMP_2019               | IBD          | ko_humann        |                     0.444 |                  0.532 |                         0.032 |                      0.07  |
| HMP_2019               | IBD          | pathway_humann   |                     0.444 |                  0.323 |                         0.053 |                      0.059 |
| HMP_2019               | IBD          | taxonomy_family  |                     0.1   |                  0.016 |                         0.109 |                      0.067 |
| HMP_2019               | IBD          | taxonomy_genus   |                     0.1   |                  0.113 |                         0.091 |                      0.09  |
| HMP_2019               | IBD          | taxonomy_species |                     0.1   |                  0.145 |                         0.073 |                      0.071 |
| HanniganGD_2017        | CRC          | ko_humann        |                     0.429 |                  0.36  |                         0.052 |                      0.052 |
| HanniganGD_2017        | CRC          | pathway_humann   |                     0.357 |                  0.64  |                         0.071 |                      0.093 |
| HanniganGD_2017        | CRC          | taxonomy_family  |                     0     |                  0     |                         0.056 |                      0.05  |
| HanniganGD_2017        | CRC          | taxonomy_genus   |                     0.071 |                  0     |                         0.072 |                      0.065 |
| HanniganGD_2017        | CRC          | taxonomy_species |                     0     |                  0     |                         0.048 |                      0.043 |
| KarlssonFH_2013        | T2D          | ko_humann        |                     0.167 |                  0.226 |                         0.018 |                      0.022 |
| KarlssonFH_2013        | T2D          | pathway_humann   |                     0.095 |                  0.075 |                         0.017 |                      0.019 |
| KarlssonFH_2013        | T2D          | taxonomy_family  |                     0     |                  0.019 |                         0     |                      0.036 |
| KarlssonFH_2013        | T2D          | taxonomy_genus   |                     0     |                  0.075 |                         0.021 |                      0.033 |
| KarlssonFH_2013        | T2D          | taxonomy_species |                     0     |                  0.057 |                         0.029 |                      0.045 |
| NagySzakalD_2017       | ME/CFS       | ko_humann        |                     0.408 |                  0.64  |                         0.045 |                      0.064 |
| NagySzakalD_2017       | ME/CFS       | pathway_humann   |                     0.327 |                  0.44  |                         0.044 |                      0.077 |
| NagySzakalD_2017       | ME/CFS       | taxonomy_family  |                     0.04  |                  0.1   |                         0.049 |                      0.083 |
| NagySzakalD_2017       | ME/CFS       | taxonomy_genus   |                     0.06  |                  0.22  |                         0.064 |                      0.096 |
| NagySzakalD_2017       | ME/CFS       | taxonomy_species |                     0.16  |                  0.34  |                         0.079 |                      0.087 |
| RubelMA_2020           | STH          | ko_humann        |                     0.143 |                  0.364 |                         0.03  |                      0.044 |
| RubelMA_2020           | STH          | pathway_humann   |                     0.169 |                  0.143 |                         0.027 |                      0.026 |
| RubelMA_2020           | STH          | taxonomy_family  |                     0.024 |                  0     |                         0.039 |                      0     |
| RubelMA_2020           | STH          | taxonomy_genus   |                     0     |                  0     |                         0.048 |                      0.053 |
| RubelMA_2020           | STH          | taxonomy_species |                     0.024 |                  0     |                         0.043 |                      0.043 |
| SankaranarayananK_2015 | T2D          | ko_humann        |                     0.222 |                  0.421 |                         0.04  |                      0.041 |
| SankaranarayananK_2015 | T2D          | pathway_humann   |                     0.111 |                  0.263 |                         0.024 |                      0.034 |
| SankaranarayananK_2015 | T2D          | taxonomy_family  |                     0.056 |                  0     |                         0.057 |                      0.042 |
| SankaranarayananK_2015 | T2D          | taxonomy_genus   |                     0.222 |                  0.211 |                         0.056 |                      0.054 |
| SankaranarayananK_2015 | T2D          | taxonomy_species |                     0.278 |                  0.211 |                         0.048 |                      0.074 |
| ThomasAM_2019_a        | CRC          | ko_humann        |                     0.409 |                  0.393 |                         0.043 |                      0.043 |
| ThomasAM_2019_a        | CRC          | pathway_humann   |                     0.273 |                  0.429 |                         0.035 |                      0.052 |
| ThomasAM_2019_a        | CRC          | taxonomy_family  |                     0.136 |                  0.143 |                         0.077 |                      0.098 |
| ThomasAM_2019_a        | CRC          | taxonomy_genus   |                     0.182 |                  0.214 |                         0.103 |                      0.109 |
| ThomasAM_2019_a        | CRC          | taxonomy_species |                     0.318 |                  0.357 |                         0.082 |                      0.094 |
| ThomasAM_2019_b        | CRC          | ko_humann        |                     0.148 |                  0.194 |                         0.016 |                      0.031 |
| ThomasAM_2019_b        | CRC          | pathway_humann   |                     0.074 |                  0.032 |                         0.014 |                      0.014 |
| ThomasAM_2019_b        | CRC          | taxonomy_family  |                     0     |                  0     |                         0.032 |                      0.034 |
| ThomasAM_2019_b        | CRC          | taxonomy_genus   |                     0.037 |                  0.065 |                         0.033 |                      0.045 |
| ThomasAM_2019_b        | CRC          | taxonomy_species |                     0.111 |                  0.065 |                         0.038 |                      0.047 |
| VogtmannE_2016         | CRC          | ko_humann        |                     0.265 |                  0.388 |                         0.022 |                      0.035 |
| VogtmannE_2016         | CRC          | pathway_humann   |                     0.122 |                  0.122 |                         0.024 |                      0.024 |
| VogtmannE_2016         | CRC          | taxonomy_family  |                     0.019 |                  0.041 |                         0.046 |                      0.042 |
| VogtmannE_2016         | CRC          | taxonomy_genus   |                     0.096 |                  0.163 |                         0.066 |                      0.062 |
| VogtmannE_2016         | CRC          | taxonomy_species |                     0.212 |                  0.286 |                         0.069 |                      0.07  |
| WirbelJ_2018           | CRC          | ko_humann        |                     0.344 |                  0.333 |                         0.031 |                      0.033 |
| WirbelJ_2018           | CRC          | pathway_humann   |                     0.188 |                  0.05  |                         0.025 |                      0.027 |
| WirbelJ_2018           | CRC          | taxonomy_family  |                     0.031 |                  0.017 |                         0.037 |                      0.046 |
| WirbelJ_2018           | CRC          | taxonomy_genus   |                     0.062 |                  0.083 |                         0.038 |                      0.049 |
| WirbelJ_2018           | CRC          | taxonomy_species |                     0.077 |                  0.15  |                         0.034 |                      0.062 |
| YuJ_2015               | CRC          | ko_humann        |                     0.396 |                  0.324 |                         0.036 |                      0.042 |
| YuJ_2015               | CRC          | pathway_humann   |                     0.208 |                  0.162 |                         0.022 |                      0.021 |
| YuJ_2015               | CRC          | taxonomy_family  |                     0.019 |                  0.027 |                         0.065 |                      0.074 |
| YuJ_2015               | CRC          | taxonomy_genus   |                     0.113 |                  0.189 |                         0.068 |                      0.061 |
| YuJ_2015               | CRC          | taxonomy_species |                     0.113 |                  0.216 |                         0.068 |                      0.078 |
| ZellerG_2014           | CRC          | ko_humann        |                     0.22  |                  0.216 |                         0.035 |                      0.033 |
| ZellerG_2014           | CRC          | pathway_humann   |                     0.119 |                  0.137 |                         0.026 |                      0.016 |
| ZellerG_2014           | CRC          | taxonomy_family  |                     0     |                  0.02  |                         0.037 |                      0.056 |
| ZellerG_2014           | CRC          | taxonomy_genus   |                     0.085 |                  0.098 |                         0.053 |                      0.069 |
| ZellerG_2014           | CRC          | taxonomy_species |                     0.102 |                  0.255 |                         0.055 |                      0.067 |
| ZhuF_2020              | schizofrenia | ko_humann        |                     0.462 |                  0.433 |                         0.052 |                      0.054 |
| ZhuF_2020              | schizofrenia | pathway_humann   |                     0.25  |                  0.189 |                         0.041 |                      0.034 |
| ZhuF_2020              | schizofrenia | taxonomy_family  |                     0     |                  0.011 |                         0.069 |                      0.047 |
| ZhuF_2020              | schizofrenia | taxonomy_genus   |                     0.074 |                  0.067 |                         0.054 |                      0.059 |
| ZhuF_2020              | schizofrenia | taxonomy_species |                     0.086 |                  0.133 |                         0.059 |                      0.066 |

## Feature-level shifts (case vs same-study control reference percentile; q <= 0.05 and |shift| >= 10 points)

| study                  | condition    | layer            |   features |   significant |   up_in_cases |   down_in_cases |
|:-----------------------|:-------------|:-----------------|-----------:|--------------:|--------------:|----------------:|
| FengQ_2015             | CRC          | ko_humann        |       4072 |            59 |            33 |              23 |
| FengQ_2015             | CRC          | pathway_humann   |        387 |            45 |             9 |              35 |
| FengQ_2015             | CRC          | taxonomy_family  |         37 |             7 |             0 |               7 |
| FengQ_2015             | CRC          | taxonomy_genus   |         86 |            14 |             0 |              14 |
| FengQ_2015             | CRC          | taxonomy_species |        154 |            12 |             0 |              12 |
| GuptaA_2019            | CRC          | ko_humann        |       2293 |           334 |            22 |             312 |
| GuptaA_2019            | CRC          | pathway_humann   |        315 |           114 |            36 |              75 |
| GuptaA_2019            | CRC          | taxonomy_family  |         15 |             4 |             3 |               1 |
| GuptaA_2019            | CRC          | taxonomy_genus   |         24 |            10 |             5 |               3 |
| GuptaA_2019            | CRC          | taxonomy_species |         29 |             8 |             4 |               3 |
| HMP_2019               | IBD          | taxonomy_family  |         18 |             0 |             0 |               0 |
| HMP_2019               | IBD          | taxonomy_genus   |         33 |             0 |             0 |               0 |
| HMP_2019               | IBD          | taxonomy_species |         47 |             0 |             0 |               0 |
| HanniganGD_2017        | CRC          | ko_humann        |       1705 |             0 |             0 |               0 |
| HanniganGD_2017        | CRC          | pathway_humann   |        255 |             0 |             0 |               0 |
| HanniganGD_2017        | CRC          | taxonomy_family  |         16 |             0 |             0 |               0 |
| HanniganGD_2017        | CRC          | taxonomy_genus   |         29 |             0 |             0 |               0 |
| HanniganGD_2017        | CRC          | taxonomy_species |         37 |             0 |             0 |               0 |
| KarlssonFH_2013        | T2D          | ko_humann        |       3718 |             0 |             0 |               0 |
| KarlssonFH_2013        | T2D          | pathway_humann   |        372 |             0 |             0 |               0 |
| KarlssonFH_2013        | T2D          | taxonomy_family  |         32 |             0 |             0 |               0 |
| KarlssonFH_2013        | T2D          | taxonomy_genus   |         68 |             0 |             0 |               0 |
| KarlssonFH_2013        | T2D          | taxonomy_species |        123 |             2 |             2 |               0 |
| NagySzakalD_2017       | ME/CFS       | ko_humann        |       3606 |             0 |             0 |               0 |
| NagySzakalD_2017       | ME/CFS       | pathway_humann   |        368 |             0 |             0 |               0 |
| NagySzakalD_2017       | ME/CFS       | taxonomy_family  |         30 |             0 |             0 |               0 |
| NagySzakalD_2017       | ME/CFS       | taxonomy_genus   |         71 |             5 |             4 |               1 |
| NagySzakalD_2017       | ME/CFS       | taxonomy_species |        137 |             3 |             2 |               1 |
| RubelMA_2020           | STH          | ko_humann        |       3705 |           804 |           393 |             347 |
| RubelMA_2020           | STH          | pathway_humann   |        367 |            17 |             4 |              10 |
| RubelMA_2020           | STH          | taxonomy_family  |         27 |             1 |             0 |               1 |
| RubelMA_2020           | STH          | taxonomy_genus   |         49 |             4 |             0 |               3 |
| RubelMA_2020           | STH          | taxonomy_species |         87 |             3 |             0 |               1 |
| SankaranarayananK_2015 | T2D          | ko_humann        |       2092 |             0 |             0 |               0 |
| SankaranarayananK_2015 | T2D          | pathway_humann   |        290 |             0 |             0 |               0 |
| SankaranarayananK_2015 | T2D          | taxonomy_family  |         22 |             0 |             0 |               0 |
| SankaranarayananK_2015 | T2D          | taxonomy_genus   |         42 |             0 |             0 |               0 |
| SankaranarayananK_2015 | T2D          | taxonomy_species |         59 |             2 |             0 |               2 |
| ThomasAM_2019_a        | CRC          | ko_humann        |       2548 |             0 |             0 |               0 |
| ThomasAM_2019_a        | CRC          | pathway_humann   |        324 |           134 |            21 |             113 |
| ThomasAM_2019_a        | CRC          | taxonomy_family  |         21 |             0 |             0 |               0 |
| ThomasAM_2019_a        | CRC          | taxonomy_genus   |         44 |             1 |             0 |               1 |
| ThomasAM_2019_a        | CRC          | taxonomy_species |         61 |             0 |             0 |               0 |
| ThomasAM_2019_b        | CRC          | ko_humann        |       2779 |             0 |             0 |               0 |
| ThomasAM_2019_b        | CRC          | pathway_humann   |        342 |             0 |             0 |               0 |
| ThomasAM_2019_b        | CRC          | taxonomy_family  |         32 |             9 |             0 |               9 |
| ThomasAM_2019_b        | CRC          | taxonomy_genus   |         62 |             1 |             0 |               1 |
| ThomasAM_2019_b        | CRC          | taxonomy_species |        111 |             0 |             0 |               0 |
| VogtmannE_2016         | CRC          | ko_humann        |       3587 |             0 |             0 |               0 |
| VogtmannE_2016         | CRC          | pathway_humann   |        373 |             0 |             0 |               0 |
| VogtmannE_2016         | CRC          | taxonomy_family  |         40 |             0 |             0 |               0 |
| VogtmannE_2016         | CRC          | taxonomy_genus   |         89 |             1 |             0 |               1 |
| VogtmannE_2016         | CRC          | taxonomy_species |        169 |             1 |             0 |               1 |
| WirbelJ_2018           | CRC          | ko_humann        |       3848 |             7 |             6 |               1 |
| WirbelJ_2018           | CRC          | pathway_humann   |        379 |            10 |             1 |               9 |
| WirbelJ_2018           | CRC          | taxonomy_family  |         37 |             0 |             0 |               0 |
| WirbelJ_2018           | CRC          | taxonomy_genus   |         79 |             2 |             1 |               1 |
| WirbelJ_2018           | CRC          | taxonomy_species |        153 |             1 |             1 |               0 |
| YuJ_2015               | CRC          | ko_humann        |       4085 |             0 |             0 |               0 |
| YuJ_2015               | CRC          | pathway_humann   |        407 |             0 |             0 |               0 |
| YuJ_2015               | CRC          | taxonomy_family  |         38 |             3 |             1 |               2 |
| YuJ_2015               | CRC          | taxonomy_genus   |         87 |             7 |             2 |               5 |
| YuJ_2015               | CRC          | taxonomy_species |        164 |             6 |             0 |               6 |
| ZellerG_2014           | CRC          | ko_humann        |       4121 |           188 |            43 |             145 |
| ZellerG_2014           | CRC          | pathway_humann   |        392 |            27 |             2 |              25 |
| ZellerG_2014           | CRC          | taxonomy_family  |         39 |             6 |             1 |               5 |
| ZellerG_2014           | CRC          | taxonomy_genus   |         92 |            10 |             0 |              10 |
| ZellerG_2014           | CRC          | taxonomy_species |        177 |             1 |             0 |               1 |
| ZhuF_2020              | schizofrenia | ko_humann        |       4142 |             0 |             0 |               0 |
| ZhuF_2020              | schizofrenia | pathway_humann   |        412 |             0 |             0 |               0 |
| ZhuF_2020              | schizofrenia | taxonomy_family  |         34 |             0 |             0 |               0 |
| ZhuF_2020              | schizofrenia | taxonomy_genus   |         78 |             2 |             1 |               0 |
| ZhuF_2020              | schizofrenia | taxonomy_species |        145 |             1 |             1 |               0 |

## CRC: features shifted in the same direction in >= 2 of 9 studies

| condition   | layer           | feature_id        |   n_studies |   n_up |   n_down |   median_shift | consistent   |
|:------------|:----------------|:------------------|------------:|-------:|---------:|---------------:|:-------------|
| CRC         | pathway_humann  | GLYCOGENSYNTH-PWY |           4 |      0 |        4 |         -22.54 | down         |
| CRC         | pathway_humann  | TRPSYN-PWY        |           4 |      0 |        4 |         -25.27 | down         |
| CRC         | pathway_humann  | PWY66-422         |           4 |      0 |        4 |         -26.71 | down         |
| CRC         | taxonomy_family | 186803            |           4 |      0 |        4 |         -34.87 | down         |
| CRC         | pathway_humann  | ARGSYNBSUB-PWY    |           3 |      0 |        3 |         -14.35 | down         |
| CRC         | pathway_humann  | ARGSYN-PWY        |           3 |      0 |        3 |         -14.39 | down         |
| CRC         | pathway_humann  | PWY-7400          |           3 |      0 |        3 |         -14.94 | down         |
| CRC         | pathway_humann  | GLUTORN-PWY       |           3 |      0 |        3 |         -17.86 | down         |
| CRC         | pathway_humann  | HSERMETANA-PWY    |           3 |      0 |        3 |         -18.08 | down         |
| CRC         | pathway_humann  | PWY0-1061         |           3 |      0 |        3 |         -19.41 | down         |
| CRC         | ko_humann       | K03488            |           3 |      0 |        3 |         -22.09 | down         |
| CRC         | ko_humann       | K06958            |           3 |      0 |        3 |         -22.63 | down         |
| CRC         | taxonomy_family | 1643826           |           3 |      0 |        3 |         -22.85 | down         |
| CRC         | pathway_humann  | NONMEVIPP-PWY     |           3 |      0 |        3 |         -23.76 | down         |
| CRC         | pathway_humann  | PWY-6549          |           3 |      0 |        3 |         -23.79 | down         |
| CRC         | pathway_humann  | PWY-7115          |           3 |      0 |        3 |         -23.8  | down         |
| CRC         | ko_humann       | K09157            |           3 |      0 |        3 |         -24.13 | down         |
| CRC         | taxonomy_genus  | 1407607           |           3 |      0 |        3 |         -24.57 | down         |
| CRC         | pathway_humann  | PWY-6527          |           3 |      0 |        3 |         -24.71 | down         |
| CRC         | taxonomy_genus  | 28050             |           3 |      0 |        3 |         -24.76 | down         |
| CRC         | pathway_humann  | PWY-6317          |           3 |      0 |        3 |         -27    | down         |
| CRC         | pathway_humann  | PWY-6270          |           3 |      0 |        3 |         -27.78 | down         |
| CRC         | pathway_humann  | NONOXIPENT-PWY    |           3 |      0 |        3 |         -28.18 | down         |
| CRC         | taxonomy_family | 31953             |           3 |      0 |        3 |         -28.96 | down         |
| CRC         | pathway_humann  | PWY-7560          |           3 |      0 |        3 |         -29.14 | down         |

## T2D: features shifted in the same direction in >= 2 of 2 studies

none
