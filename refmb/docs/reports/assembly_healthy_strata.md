# S12.2 — healthy baseline across strata, pipeline A

Out-of-sample scores (per-study LOSO for pool studies; base bundle for held-out studies): 2,174 healthy samples, 29 studies.

## By age

| age_bin   | layer           |   n |   studies |   median_frac_outside |   share_over_gate |   excess_positive |   core_dist_pct |   vs_pooled_pts |
|:----------|:----------------|----:|----------:|----------------------:|------------------:|------------------:|----------------:|----------------:|
| 18-39     | ko_eggnog       | 367 |        11 |                 0.04  |             0.253 |             0.063 |          92.221 |            -0.5 |
| 18-39     | ko_kofam        | 367 |        11 |                 0.04  |             0.294 |             0.074 |          92.221 |            -0.6 |
| 18-39     | module          | 367 |        11 |                 0     |             0.005 |             0     |          92.221 |             0   |
| 18-39     | pfam            | 367 |        11 |                 0.043 |             0.302 |             0.074 |          92.221 |            -0.5 |
| 18-39     | taxonomy_family | 367 |        11 |                 0.047 |             0.24  |             0.03  |          92.221 |            -0.5 |
| 18-39     | taxonomy_genus  | 367 |        11 |                 0.047 |             0.193 |             0.052 |          92.221 |            -0.6 |
| 40-64     | ko_eggnog       | 486 |        12 |                 0.029 |             0.146 |           nan     |          88.218 |            -1.6 |
| 40-64     | ko_kofam        | 486 |        12 |                 0.031 |             0.191 |           nan     |          88.218 |            -1.5 |
| 40-64     | module          | 486 |        12 |                 0     |             0.002 |           nan     |          88.218 |             0   |
| 40-64     | pfam            | 486 |        12 |                 0.031 |             0.181 |           nan     |          88.218 |            -1.7 |
| 40-64     | taxonomy_family | 486 |        12 |                 0.04  |             0.113 |           nan     |          88.218 |            -1.2 |
| 40-64     | taxonomy_genus  | 486 |        12 |                 0.042 |             0.095 |           nan     |          88.218 |            -1.1 |
| 65+       | ko_eggnog       | 219 |         8 |                 0.026 |             0.142 |           nan     |          86.786 |            -1.9 |
| 65+       | ko_kofam        | 219 |         8 |                 0.028 |             0.142 |           nan     |          86.786 |            -1.8 |
| 65+       | module          | 219 |         8 |                 0     |             0.005 |           nan     |          86.786 |             0   |
| 65+       | pfam            | 219 |         8 |                 0.03  |             0.183 |           nan     |          86.786 |            -1.8 |
| 65+       | taxonomy_family | 219 |         8 |                 0.044 |             0.096 |           nan     |          86.786 |            -0.8 |
| 65+       | taxonomy_genus  | 219 |         8 |                 0.049 |             0.164 |           nan     |          86.786 |            -0.4 |

## By region

| region     | layer           |   n |   studies |   median_frac_outside |   share_over_gate |   excess_positive |   core_dist_pct |   vs_pooled_pts |
|:-----------|:----------------|----:|----------:|----------------------:|------------------:|------------------:|----------------:|----------------:|
| Africa     | ko_eggnog       | 372 |         2 |                 0.092 |             0.567 |           nan     |          96.256 |             4.7 |
| Africa     | ko_kofam        | 372 |         2 |                 0.085 |             0.546 |           nan     |          96.256 |             3.9 |
| Africa     | module          | 372 |         2 |                 0.009 |             0.022 |           nan     |          96.256 |             0.9 |
| Africa     | pfam            | 372 |         2 |                 0.1   |             0.624 |           nan     |          96.256 |             5.2 |
| Africa     | taxonomy_family | 372 |         2 |                 0.198 |             0.903 |           nan     |          96.256 |            14.6 |
| Africa     | taxonomy_genus  | 372 |         2 |                 0.226 |             0.898 |           nan     |          96.256 |            17.3 |
| E. Asia    | ko_eggnog       | 304 |         6 |                 0.051 |             0.362 |           nan     |          86.882 |             0.6 |
| E. Asia    | ko_kofam        | 304 |         6 |                 0.049 |             0.362 |           nan     |          86.882 |             0.3 |
| E. Asia    | module          | 304 |         6 |                 0     |             0     |           nan     |          86.882 |             0   |
| E. Asia    | pfam            | 304 |         6 |                 0.051 |             0.359 |           nan     |          86.882 |             0.3 |
| E. Asia    | taxonomy_family | 304 |         6 |                 0.056 |             0.24  |           nan     |          86.882 |             0.4 |
| E. Asia    | taxonomy_genus  | 304 |         6 |                 0.055 |             0.164 |           nan     |          86.882 |             0.2 |
| Europe     | ko_eggnog       | 440 |        10 |                 0.028 |             0.157 |           nan     |          90.18  |            -1.7 |
| Europe     | ko_kofam        | 440 |        10 |                 0.03  |             0.193 |           nan     |          90.18  |            -1.6 |
| Europe     | module          | 440 |        10 |                 0     |             0.009 |           nan     |          90.18  |             0   |
| Europe     | pfam            | 440 |        10 |                 0.029 |             0.182 |           nan     |          90.18  |            -1.9 |
| Europe     | taxonomy_family | 440 |        10 |                 0.037 |             0.111 |           nan     |          90.18  |            -1.5 |
| Europe     | taxonomy_genus  | 440 |        10 |                 0.039 |             0.066 |           nan     |          90.18  |            -1.4 |
| N. America | ko_eggnog       | 750 |        10 |                 0.03  |             0.161 |             0.032 |          87.683 |            -1.5 |
| N. America | ko_kofam        | 750 |        10 |                 0.03  |             0.171 |             0.037 |          87.683 |            -1.6 |
| N. America | module          | 750 |        10 |                 0     |             0.001 |             0     |          87.683 |             0   |
| N. America | pfam            | 750 |        10 |                 0.031 |             0.184 |             0.04  |          87.683 |            -1.7 |
| N. America | taxonomy_family | 750 |        10 |                 0.044 |             0.163 |             0.015 |          87.683 |            -0.8 |
| N. America | taxonomy_genus  | 750 |        10 |                 0.047 |             0.195 |             0.025 |          87.683 |            -0.6 |
| S. Asia    | ko_eggnog       |  88 |         1 |                 0.06  |             0.295 |           nan     |          95.035 |             1.5 |
| S. Asia    | ko_kofam        |  88 |         1 |                 0.066 |             0.364 |           nan     |          95.035 |             2   |
| S. Asia    | module          |  88 |         1 |                 0     |             0     |           nan     |          95.035 |             0   |
| S. Asia    | pfam            |  88 |         1 |                 0.064 |             0.386 |           nan     |          95.035 |             1.6 |
| S. Asia    | taxonomy_family |  88 |         1 |                 0.061 |             0.352 |           nan     |          95.035 |             0.9 |
| S. Asia    | taxonomy_genus  |  88 |         1 |                 0.062 |             0.386 |           nan     |          95.035 |             0.9 |
| W. Asia    | ko_eggnog       |  59 |         1 |                 0.045 |             0.254 |           nan     |          97.715 |            -0   |
| W. Asia    | ko_kofam        |  59 |         1 |                 0.047 |             0.254 |           nan     |          97.715 |             0.1 |
| W. Asia    | module          |  59 |         1 |                 0.017 |             0.237 |           nan     |          97.715 |             1.7 |
| W. Asia    | pfam            |  59 |         1 |                 0.052 |             0.305 |           nan     |          97.715 |             0.4 |
| W. Asia    | taxonomy_family |  59 |         1 |                 0.053 |             0.322 |           nan     |          97.715 |             0.1 |
| W. Asia    | taxonomy_genus  |  59 |         1 |                 0.053 |             0.373 |           nan     |          97.715 |             0   |
| other      | ko_eggnog       |   6 |         2 |                 0.061 |             0.5   |           nan     |          92.306 |             1.6 |
| other      | ko_kofam        |   6 |         2 |                 0.041 |             0.333 |           nan     |          92.306 |            -0.5 |
| other      | module          |   6 |         2 |                 0.007 |             0.167 |           nan     |          92.306 |             0.7 |
| other      | pfam            |   6 |         2 |                 0.055 |             0.5   |           nan     |          92.306 |             0.7 |
| other      | taxonomy_family |   6 |         2 |                 0.055 |             0.333 |           nan     |          92.306 |             0.3 |
| other      | taxonomy_genus  |   6 |         2 |                 0.061 |             0.333 |           nan     |          92.306 |             0.8 |

## By Westernized lifestyle

| westernized   | layer           |   n |   studies |   median_frac_outside |   share_over_gate |   excess_positive |   core_dist_pct |   vs_pooled_pts |
|:--------------|:----------------|----:|----------:|----------------------:|------------------:|------------------:|----------------:|----------------:|
| no            | ko_eggnog       |  22 |         1 |                 0.049 |             0.182 |           nan     |          94.503 |             0.4 |
| no            | ko_kofam        |  22 |         1 |                 0.071 |             0.455 |           nan     |          94.503 |             2.5 |
| no            | module          |  22 |         1 |                 0.006 |             0     |           nan     |          94.503 |             0.6 |
| no            | pfam            |  22 |         1 |                 0.062 |             0.364 |           nan     |          94.503 |             1.4 |
| no            | taxonomy_family |  22 |         1 |                 0.052 |             0.182 |           nan     |          94.503 |             0   |
| no            | taxonomy_genus  |  22 |         1 |                 0.052 |             0.091 |           nan     |          94.503 |            -0.1 |
| yes           | ko_eggnog       | 627 |        10 |                 0.036 |             0.238 |             0.037 |          90.07  |            -0.9 |
| yes           | ko_kofam        | 627 |        10 |                 0.04  |             0.274 |             0.043 |          90.07  |            -0.6 |
| yes           | module          | 627 |        10 |                 0     |             0.006 |             0     |          90.07  |             0   |
| yes           | pfam            | 627 |        10 |                 0.04  |             0.281 |             0.043 |          90.07  |            -0.8 |
| yes           | taxonomy_family | 627 |        10 |                 0.043 |             0.188 |             0.018 |          90.07  |            -0.9 |
| yes           | taxonomy_genus  | 627 |        10 |                 0.043 |             0.148 |             0.03  |          90.07  |            -1   |

## By sex

| sex    | layer           |   n |   studies |   median_frac_outside |   share_over_gate |   excess_positive |   core_dist_pct |   vs_pooled_pts |
|:-------|:----------------|----:|----------:|----------------------:|------------------:|------------------:|----------------:|----------------:|
| female | ko_eggnog       | 715 |        16 |                 0.036 |             0.218 |             0.014 |          88.17  |            -0.9 |
| female | ko_kofam        | 715 |        16 |                 0.037 |             0.245 |             0.02  |          88.17  |            -0.9 |
| female | module          | 715 |        16 |                 0     |             0.001 |             0     |          88.17  |             0   |
| female | pfam            | 715 |        16 |                 0.038 |             0.257 |             0.02  |          88.17  |            -1   |
| female | taxonomy_family | 715 |        16 |                 0.048 |             0.227 |             0.004 |          88.17  |            -0.4 |
| female | taxonomy_genus  | 715 |        16 |                 0.051 |             0.232 |             0.013 |          88.17  |            -0.2 |
| male   | ko_eggnog       | 662 |        17 |                 0.037 |             0.254 |             0.021 |          90.533 |            -0.8 |
| male   | ko_kofam        | 662 |        17 |                 0.041 |             0.27  |             0.021 |          90.533 |            -0.5 |
| male   | module          | 662 |        17 |                 0     |             0.009 |             0     |          90.533 |             0   |
| male   | pfam            | 662 |        17 |                 0.041 |             0.29  |             0.024 |          90.533 |            -0.7 |
| male   | taxonomy_family | 662 |        17 |                 0.049 |             0.304 |             0.012 |          90.533 |            -0.3 |
| male   | taxonomy_genus  | 662 |        17 |                 0.05  |             0.26  |             0.015 |          90.533 |            -0.3 |

## By BMI class

| bmi_class   | layer           |   n |   studies |   median_frac_outside |   share_over_gate |   excess_positive |   core_dist_pct |   vs_pooled_pts |
|:------------|:----------------|----:|----------:|----------------------:|------------------:|------------------:|----------------:|----------------:|
| 18.5-25     | ko_eggnog       | 316 |         6 |                 0.043 |             0.288 |             0.057 |          90.616 |            -0.2 |
| 18.5-25     | ko_kofam        | 316 |         6 |                 0.047 |             0.332 |             0.07  |          90.616 |             0.1 |
| 18.5-25     | module          | 316 |         6 |                 0     |             0.003 |             0     |          90.616 |             0   |
| 18.5-25     | pfam            | 316 |         6 |                 0.046 |             0.332 |             0.07  |          90.616 |            -0.2 |
| 18.5-25     | taxonomy_family | 316 |         6 |                 0.047 |             0.237 |             0.025 |          90.616 |            -0.5 |
| 18.5-25     | taxonomy_genus  | 316 |         6 |                 0.047 |             0.187 |             0.047 |          90.616 |            -0.6 |
| 25-30       | ko_eggnog       | 111 |         6 |                 0.033 |             0.243 |             0.045 |          90.557 |            -1.2 |
| 25-30       | ko_kofam        | 111 |         6 |                 0.04  |             0.261 |             0.045 |          90.557 |            -0.6 |
| 25-30       | module          | 111 |         6 |                 0     |             0.009 |             0     |          90.557 |             0   |
| 25-30       | pfam            | 111 |         6 |                 0.041 |             0.288 |             0.045 |          90.557 |            -0.7 |
| 25-30       | taxonomy_family | 111 |         6 |                 0.04  |             0.207 |             0.027 |          90.557 |            -1.2 |
| 25-30       | taxonomy_genus  | 111 |         6 |                 0.044 |             0.198 |             0.036 |          90.557 |            -0.9 |
| <18.5       | ko_eggnog       |  25 |         5 |                 0.044 |             0.24  |           nan     |          87.959 |            -0.1 |
| <18.5       | ko_kofam        |  25 |         5 |                 0.037 |             0.28  |           nan     |          87.959 |            -0.9 |
| <18.5       | module          |  25 |         5 |                 0     |             0     |           nan     |          87.959 |             0   |
| <18.5       | pfam            |  25 |         5 |                 0.047 |             0.28  |           nan     |          87.959 |            -0.1 |
| <18.5       | taxonomy_family |  25 |         5 |                 0.047 |             0.2   |           nan     |          87.959 |            -0.5 |
| <18.5       | taxonomy_genus  |  25 |         5 |                 0.051 |             0.08  |           nan     |          87.959 |            -0.2 |
| >=30        | ko_eggnog       | 114 |         5 |                 0.027 |             0.132 |           nan     |          90.235 |            -1.8 |
| >=30        | ko_kofam        | 114 |         5 |                 0.034 |             0.193 |           nan     |          90.235 |            -1.2 |
| >=30        | module          | 114 |         5 |                 0     |             0.009 |           nan     |          90.235 |             0   |
| >=30        | pfam            | 114 |         5 |                 0.03  |             0.193 |           nan     |          90.235 |            -1.8 |
| >=30        | taxonomy_family | 114 |         5 |                 0.038 |             0.096 |           nan     |          90.235 |            -1.4 |
| >=30        | taxonomy_genus  | 114 |         5 |                 0.039 |             0.079 |           nan     |          90.235 |            -1.4 |

## Stratum-specific reference rule (>= +3 points on >= 2 layers, n >= 120, >= 3 studies)

- no stratum fires the rule

## Age 65+ vs 18-39 (219 vs 367 samples): features with q <= 0.05 and |shift| >= 15 percentile points

| layer | features tested | significant | shifted up in 65+ | shifted down |
|---|---:|---:|---:|---:|
| ko_eggnog | 4631 | 2674 | 1419 | 357 |
| ko_kofam | 2829 | 1655 | 874 | 200 |
| module | 238 | 0 | 0 | 0 |
| pfam | 4575 | 2592 | 1270 | 389 |
| taxonomy_family | 408 | 102 | 68 | 15 |
| taxonomy_genus | 1080 | 306 | 202 | 52 |

Top 15 by p:
| layer           | feature_id   |   n_a |   n_b |   median_pct_a |   median_pct_b |   shift |   q |
|:----------------|:-------------|------:|------:|---------------:|---------------:|--------:|----:|
| taxonomy_genus  | 946234       |   219 |   357 |         90.074 |         35.206 |  54.868 |   0 |
| taxonomy_family | 3082771      |   218 |   347 |         82.505 |         30.484 |  52.02  |   0 |
| taxonomy_genus  | 1017280      |   218 |   349 |         82.968 |         34.428 |  48.54  |   0 |
| taxonomy_genus  | 1924093      |   218 |   343 |         82.024 |         31.744 |  50.28  |   0 |
| taxonomy_family | 1643826      |   218 |   345 |         77.658 |         26.075 |  51.583 |   0 |
| taxonomy_genus  | 1506553      |   219 |   359 |         79.121 |         27.381 |  51.74  |   0 |
| taxonomy_genus  | 572511       |   219 |   364 |         84.741 |         23.262 |  61.479 |   0 |
| taxonomy_genus  | 244127       |   219 |   352 |         77.204 |         37.276 |  39.928 |   0 |
| taxonomy_genus  | 1649459      |   219 |   349 |         76.518 |         29.926 |  46.593 |   0 |
| taxonomy_genus  | 2719313      |   219 |   357 |         80.096 |         34.046 |  46.051 |   0 |
| taxonomy_genus  | 1924105      |   210 |   250 |         83.581 |         33.348 |  50.233 |   0 |
| pfam            | PF13936      |   219 |   362 |         74.824 |         32.928 |  41.896 |   0 |
| taxonomy_family | 186822       |   218 |   349 |         73.562 |         33.147 |  40.416 |   0 |
| taxonomy_family | 1300         |   218 |   360 |         75.234 |         28.513 |  46.721 |   0 |
| taxonomy_genus  | 44249        |   217 |   342 |         74.866 |         32.317 |  42.549 |   0 |

## Non-Western vs Western lifestyle: too few samples (22 vs 627)

## Female vs male (715 vs 662 samples): features with q <= 0.05 and |shift| >= 15 percentile points

| layer | features tested | significant | shifted up in female | shifted down |
|---|---:|---:|---:|---:|
| ko_eggnog | 5768 | 560 | 13 | 6 |
| ko_kofam | 4289 | 250 | 7 | 3 |
| module | 261 | 0 | 0 | 0 |
| pfam | 5753 | 640 | 9 | 10 |
| taxonomy_family | 584 | 16 | 1 | 0 |
| taxonomy_genus | 1687 | 37 | 5 | 0 |

Top 15 by p:
| layer           | feature_id   |   n_a |   n_b |   median_pct_a |   median_pct_b |   shift |     q |
|:----------------|:-------------|------:|------:|---------------:|---------------:|--------:|------:|
| ko_eggnog       | K10547       |   376 |   310 |         52.827 |         33.645 |  19.182 | 0     |
| pfam            | PF12642      |   549 |   450 |         54.579 |         39.548 |  15.032 | 0     |
| ko_kofam        | K13922       |   225 |   171 |         61.545 |         43.528 |  18.017 | 0.002 |
| pfam            | PF13682      |   167 |   159 |         36.395 |         63.49  | -27.094 | 0.002 |
| ko_kofam        | K10672       |   264 |   183 |         58.785 |         40.525 |  18.26  | 0.003 |
| pfam            | PF11329      |   188 |   145 |         54.518 |         34.839 |  19.679 | 0.002 |
| ko_eggnog       | K10826       |   182 |   111 |         62.163 |         46.587 |  15.576 | 0.003 |
| taxonomy_genus  | 7164         |   101 |    94 |         50.289 |         23.244 |  27.044 | 0.006 |
| taxonomy_genus  | 221065       |   115 |    97 |         43.145 |         15.312 |  27.833 | 0.009 |
| ko_eggnog       | K07469       |   326 |   264 |         64.994 |         49.832 |  15.162 | 0.004 |
| taxonomy_family | 415002       |    95 |    93 |         40.348 |         20.705 |  19.643 | 0.013 |
| taxonomy_genus  | 140458       |   165 |   132 |         51.287 |         27.932 |  23.355 | 0.016 |
| ko_kofam        | K02008       |   262 |   196 |         58.64  |         42.798 |  15.842 | 0.015 |
| ko_kofam        | K12065       |    26 |    26 |         54.793 |         23.689 |  31.104 | 0.016 |
| ko_eggnog       | K18353       |   267 |   191 |         49.724 |         33.24  |  16.484 | 0.012 |

## BMI >= 30 vs 18.5-25 (114 vs 316 samples): features with q <= 0.05 and |shift| >= 15 percentile points

| layer | features tested | significant | shifted up in >=30 | shifted down |
|---|---:|---:|---:|---:|
| ko_eggnog | 3852 | 1548 | 56 | 744 |
| ko_kofam | 2321 | 1138 | 42 | 624 |
| module | 203 | 0 | 0 | 0 |
| pfam | 4112 | 1874 | 90 | 981 |
| taxonomy_family | 233 | 11 | 2 | 5 |
| taxonomy_genus | 571 | 10 | 1 | 3 |

Top 15 by p:
| layer    | feature_id   |   n_a |   n_b |   median_pct_a |   median_pct_b |   shift |   q |
|:---------|:-------------|------:|------:|---------------:|---------------:|--------:|----:|
| pfam     | PF08436      |   114 |   316 |         14.949 |         56.935 | -41.986 |   0 |
| ko_kofam | K00384       |   114 |   313 |         12.3   |         45.353 | -33.052 |   0 |
| pfam     | PF02875      |   114 |   316 |         12.119 |         47.237 | -35.118 |   0 |
| pfam     | PF13288      |   114 |   316 |         18.322 |         56.339 | -38.017 |   0 |
| ko_kofam | K00099       |   114 |   316 |         15.759 |         54.759 | -39     |   0 |
| pfam     | PF01071      |   114 |   316 |          7.071 |         44.227 | -37.156 |   0 |
| pfam     | PF02670      |   114 |   316 |         16.239 |         55.527 | -39.288 |   0 |
| ko_kofam | K03526       |   114 |   316 |         13.995 |         54.53  | -40.535 |   0 |
| ko_kofam | K01945       |   114 |   316 |          7.823 |         47.058 | -39.234 |   0 |
| ko_kofam | K00800       |   114 |   316 |          9.361 |         50.283 | -40.922 |   0 |
| pfam     | PF08245      |   114 |   316 |         10.36  |         45.691 | -35.332 |   0 |
| pfam     | PF00679      |   114 |   316 |          9.441 |         48.541 | -39.101 |   0 |
| pfam     | PF00291      |   114 |   316 |         15.245 |         50.724 | -35.479 |   0 |
| pfam     | PF06574      |   114 |   316 |         15.931 |         53.84  | -37.91  |   0 |
| pfam     | PF00441      |   114 |   316 |         19.819 |         55.377 | -35.558 |   0 |
