# S12.2 — healthy baseline across strata, pipeline B

Out-of-sample scores (per-study LOSO for pool studies; base bundle for held-out studies): 6,811 healthy samples, 22 studies.

## By age

| age_bin   | layer            |    n |   studies |   median_frac_outside |   share_over_gate |   excess_positive |   core_dist_pct |   vs_pooled_pts |
|:----------|:-----------------|-----:|----------:|----------------------:|------------------:|------------------:|----------------:|----------------:|
| 18-39     | ko_humann        | 2117 |        17 |                 0.021 |             0.148 |             0.203 |          50.662 |            -0.5 |
| 18-39     | pathway_humann   | 2117 |        17 |                 0.024 |             0.165 |             0.155 |          50.662 |            -0.2 |
| 18-39     | taxonomy_family  | 2257 |        18 |                 0.036 |             0.204 |             0.019 |          50.688 |            -0.1 |
| 18-39     | taxonomy_genus   | 2257 |        18 |                 0.036 |             0.205 |             0.051 |          50.688 |            -0.2 |
| 18-39     | taxonomy_species | 2257 |        18 |                 0.038 |             0.18  |             0.065 |          50.688 |            -0.3 |
| 40-64     | ko_humann        | 2803 |        17 |                 0.028 |             0.168 |             0.251 |          51.748 |             0.2 |
| 40-64     | pathway_humann   | 2803 |        17 |                 0.029 |             0.158 |             0.145 |          51.748 |             0.3 |
| 40-64     | taxonomy_family  | 3061 |        18 |                 0.038 |             0.22  |             0.021 |          51.663 |             0.1 |
| 40-64     | taxonomy_genus   | 3061 |        18 |                 0.038 |             0.21  |             0.062 |          51.663 |             0   |
| 40-64     | taxonomy_species | 3061 |        18 |                 0.042 |             0.193 |             0.082 |          51.663 |             0.1 |
| 65+       | ko_humann        |  490 |        16 |                 0.033 |             0.204 |             0.271 |          53.318 |             0.7 |
| 65+       | pathway_humann   |  490 |        16 |                 0.025 |             0.159 |             0.145 |          53.318 |            -0.1 |
| 65+       | taxonomy_family  |  537 |        17 |                 0.038 |             0.242 |             0.013 |          51.748 |             0.1 |
| 65+       | taxonomy_genus   |  537 |        17 |                 0.042 |             0.235 |             0.088 |          51.748 |             0.4 |
| 65+       | taxonomy_species |  537 |        17 |                 0.05  |             0.264 |             0.14  |          51.748 |             0.9 |

## By region

| region     | layer            |    n |   studies |   median_frac_outside |   share_over_gate |   excess_positive |   core_dist_pct |   vs_pooled_pts |
|:-----------|:-----------------|-----:|----------:|----------------------:|------------------:|------------------:|----------------:|----------------:|
| Africa     | ko_humann        |  110 |         1 |                 0.102 |             0.636 |             0.764 |          94.948 |             7.6 |
| Africa     | pathway_humann   |  110 |         1 |                 0.053 |             0.391 |             0.382 |          94.948 |             2.7 |
| Africa     | taxonomy_family  |  110 |         1 |                 0.083 |             0.518 |             0.127 |          94.948 |             4.6 |
| Africa     | taxonomy_genus   |  110 |         1 |                 0.082 |             0.518 |             0.191 |          94.948 |             4.4 |
| Africa     | taxonomy_species |  110 |         1 |                 0.062 |             0.327 |             0.127 |          94.948 |             2.1 |
| C. Asia    | ko_humann        |   83 |         1 |                 0.031 |             0.181 |             0.265 |          43.502 |             0.5 |
| C. Asia    | pathway_humann   |   83 |         1 |                 0.021 |             0.157 |             0.157 |          43.502 |            -0.5 |
| C. Asia    | taxonomy_family  |   83 |         1 |                 0.04  |             0.241 |             0.024 |          43.502 |             0.3 |
| C. Asia    | taxonomy_genus   |   83 |         1 |                 0.04  |             0.217 |             0.024 |          43.502 |             0.2 |
| C. Asia    | taxonomy_species |   83 |         1 |                 0.047 |             0.169 |             0.036 |          43.502 |             0.6 |
| E. Asia    | ko_humann        |  859 |         6 |                 0.039 |             0.254 |             0.352 |          53.824 |             1.3 |
| E. Asia    | pathway_humann   |  859 |         6 |                 0.028 |             0.191 |             0.185 |          53.824 |             0.2 |
| E. Asia    | taxonomy_family  |  859 |         6 |                 0.043 |             0.313 |             0.019 |          53.824 |             0.6 |
| E. Asia    | taxonomy_genus   |  859 |         6 |                 0.059 |             0.357 |             0.129 |          53.824 |             2.1 |
| E. Asia    | taxonomy_species |  859 |         6 |                 0.064 |             0.363 |             0.185 |          53.824 |             2.3 |
| Europe     | ko_humann        | 3695 |        10 |                 0.026 |             0.144 |             0.211 |          51.726 |             0   |
| Europe     | pathway_humann   | 3695 |        10 |                 0.027 |             0.157 |             0.145 |          51.726 |             0.1 |
| Europe     | taxonomy_family  | 4326 |        11 |                 0.036 |             0.204 |             0.025 |          51.004 |            -0.1 |
| Europe     | taxonomy_genus   | 4326 |        11 |                 0.036 |             0.195 |             0.06  |          51.004 |            -0.2 |
| Europe     | taxonomy_species | 4326 |        11 |                 0.038 |             0.188 |             0.077 |          51.004 |            -0.3 |
| N. America | ko_humann        |  305 |         2 |                 0.011 |             0.066 |             0.092 |          54.72  |            -1.5 |
| N. America | pathway_humann   |  305 |         2 |                 0.02  |             0.108 |             0.098 |          54.72  |            -0.6 |
| N. America | taxonomy_family  |  305 |         2 |                 0     |             0.108 |             0.02  |          54.72  |            -3.7 |
| N. America | taxonomy_genus   |  305 |         2 |                 0.023 |             0.118 |             0.013 |          54.72  |            -1.5 |
| N. America | taxonomy_species |  305 |         2 |                 0.032 |             0.102 |             0.023 |          54.72  |            -0.9 |
| Oceania    | ko_humann        |  152 |         1 |                 0.047 |             0.224 |             0.355 |          83.002 |             2.1 |
| Oceania    | pathway_humann   |  152 |         1 |                 0.032 |             0.145 |             0.138 |          83.002 |             0.6 |
| Oceania    | taxonomy_family  |  152 |         1 |                 0.043 |             0.27  |             0     |          83.002 |             0.6 |
| Oceania    | taxonomy_genus   |  152 |         1 |                 0.058 |             0.171 |             0.007 |          83.002 |             2   |
| Oceania    | taxonomy_species |  152 |         1 |                 0.053 |             0.164 |             0.02  |          83.002 |             1.2 |
| S. Asia    | ko_humann        |   88 |         1 |                 0.065 |             0.409 |             0.557 |          94.06  |             3.9 |
| S. Asia    | pathway_humann   |   88 |         1 |                 0.08  |             0.5   |             0.477 |          94.06  |             5.4 |
| S. Asia    | taxonomy_family  |   88 |         1 |                 0.077 |             0.466 |             0     |          94.06  |             4   |
| S. Asia    | taxonomy_genus   |   88 |         1 |                 0.071 |             0.409 |             0     |          94.06  |             3.3 |
| S. Asia    | taxonomy_species |   88 |         1 |                 0.038 |             0.216 |             0     |          94.06  |            -0.3 |
| W. Asia    | ko_humann        |  888 |         1 |                 0.013 |             0.087 |             0.114 |          30.093 |            -1.3 |
| W. Asia    | pathway_humann   |  888 |         1 |                 0.018 |             0.092 |             0.078 |          30.093 |            -0.8 |
| W. Asia    | taxonomy_family  |  888 |         1 |                 0.036 |             0.145 |             0     |          30.093 |            -0.1 |
| W. Asia    | taxonomy_genus   |  888 |         1 |                 0.031 |             0.088 |             0.001 |          30.093 |            -0.7 |
| W. Asia    | taxonomy_species |  888 |         1 |                 0.036 |             0.059 |             0.009 |          30.093 |            -0.5 |

## By Westernized lifestyle

| westernized   | layer            |    n |   studies |   median_frac_outside |   share_over_gate |   excess_positive |   core_dist_pct |   vs_pooled_pts |
|:--------------|:-----------------|-----:|----------:|----------------------:|------------------:|------------------:|----------------:|----------------:|
| no            | ko_humann        |  326 |         3 |                 0.061 |             0.38  |             0.518 |          84.406 |             3.5 |
| no            | pathway_humann   |  326 |         3 |                 0.038 |             0.255 |             0.248 |          84.406 |             1.2 |
| no            | taxonomy_family  |  326 |         3 |                 0.062 |             0.374 |             0.049 |          84.406 |             2.5 |
| no            | taxonomy_genus   |  326 |         3 |                 0.068 |             0.347 |             0.101 |          84.406 |             3   |
| no            | taxonomy_species |  326 |         3 |                 0.06  |             0.27  |             0.083 |          84.406 |             1.9 |
| yes           | ko_humann        | 5854 |        19 |                 0.024 |             0.15  |             0.214 |          48.655 |            -0.2 |
| yes           | pathway_humann   | 5854 |        19 |                 0.025 |             0.153 |             0.142 |          48.655 |            -0.1 |
| yes           | taxonomy_family  | 6485 |        20 |                 0.037 |             0.208 |             0.02  |          48.334 |            -0   |
| yes           | taxonomy_genus   | 6485 |        20 |                 0.036 |             0.199 |             0.057 |          48.334 |            -0.2 |
| yes           | taxonomy_species | 6485 |        20 |                 0.04  |             0.187 |             0.077 |          48.334 |            -0.1 |

## By sex

| sex    | layer            |    n |   studies |   median_frac_outside |   share_over_gate |   excess_positive |   core_dist_pct |   vs_pooled_pts |
|:-------|:-----------------|-----:|----------:|----------------------:|------------------:|------------------:|----------------:|----------------:|
| female | ko_humann        | 3388 |        18 |                 0.026 |             0.158 |             0.229 |          49.663 |             0   |
| female | pathway_humann   | 3388 |        18 |                 0.026 |             0.149 |             0.138 |          49.663 |            -0   |
| female | taxonomy_family  | 3696 |        19 |                 0.037 |             0.209 |             0.018 |          49.591 |            -0   |
| female | taxonomy_genus   | 3696 |        19 |                 0.036 |             0.195 |             0.051 |          49.591 |            -0.2 |
| female | taxonomy_species | 3696 |        19 |                 0.041 |             0.189 |             0.072 |          49.591 |             0   |
| male   | ko_humann        | 2369 |        17 |                 0.024 |             0.16  |             0.225 |          53.908 |            -0.2 |
| male   | pathway_humann   | 2369 |        17 |                 0.026 |             0.17  |             0.158 |          53.908 |            -0   |
| male   | taxonomy_family  | 2506 |        18 |                 0.037 |             0.213 |             0.022 |          53.821 |            -0   |
| male   | taxonomy_genus   | 2506 |        18 |                 0.038 |             0.22  |             0.067 |          53.821 |             0   |
| male   | taxonomy_species | 2506 |        18 |                 0.039 |             0.19  |             0.085 |          53.821 |            -0.2 |

## By BMI class

| bmi_class   | layer            |    n |   studies |   median_frac_outside |   share_over_gate |   excess_positive |   core_dist_pct |   vs_pooled_pts |
|:------------|:-----------------|-----:|----------:|----------------------:|------------------:|------------------:|----------------:|----------------:|
| 18.5-25     | ko_humann        | 2863 |        16 |                 0.025 |             0.159 |             0.228 |          51.391 |            -0.1 |
| 18.5-25     | pathway_humann   | 2863 |        16 |                 0.024 |             0.146 |             0.138 |          51.391 |            -0.2 |
| 18.5-25     | taxonomy_family  | 3032 |        17 |                 0.037 |             0.213 |             0.02  |          51.139 |            -0   |
| 18.5-25     | taxonomy_genus   | 3032 |        17 |                 0.038 |             0.216 |             0.064 |          51.139 |             0   |
| 18.5-25     | taxonomy_species | 3032 |        17 |                 0.042 |             0.202 |             0.083 |          51.139 |             0.1 |
| 25-30       | ko_humann        | 1463 |        16 |                 0.027 |             0.17  |             0.24  |          49.245 |             0.1 |
| 25-30       | pathway_humann   | 1463 |        16 |                 0.029 |             0.176 |             0.16  |          49.245 |             0.3 |
| 25-30       | taxonomy_family  | 1480 |        17 |                 0.038 |             0.244 |             0.029 |          49.182 |             0.1 |
| 25-30       | taxonomy_genus   | 1480 |        17 |                 0.038 |             0.236 |             0.074 |          49.182 |             0   |
| 25-30       | taxonomy_species | 1480 |        17 |                 0.042 |             0.22  |             0.103 |          49.182 |             0.1 |
| <18.5       | ko_humann        |  880 |        20 |                 0.029 |             0.175 |             0.244 |          53.814 |             0.3 |
| <18.5       | pathway_humann   |  880 |        20 |                 0.025 |             0.153 |             0.15  |          53.814 |            -0.1 |
| <18.5       | taxonomy_family  |  882 |        21 |                 0.036 |             0.198 |             0.02  |          53.757 |            -0.1 |
| <18.5       | taxonomy_genus   |  882 |        21 |                 0.041 |             0.204 |             0.056 |          53.757 |             0.3 |
| <18.5       | taxonomy_species |  882 |        21 |                 0.042 |             0.194 |             0.06  |          53.757 |             0.1 |
| >=30        | ko_humann        |  899 |        14 |                 0.022 |             0.155 |             0.217 |          45.816 |            -0.4 |
| >=30        | pathway_humann   |  899 |        14 |                 0.031 |             0.181 |             0.16  |          45.816 |             0.5 |
| >=30        | taxonomy_family  | 1155 |        15 |                 0.037 |             0.222 |             0.023 |          47.353 |            -0   |
| >=30        | taxonomy_genus   | 1155 |        15 |                 0.035 |             0.173 |             0.041 |          47.353 |            -0.3 |
| >=30        | taxonomy_species | 1155 |        15 |                 0.039 |             0.154 |             0.06  |          47.353 |            -0.2 |

## Stratum-specific reference rule (>= +3 points on >= 2 layers, n >= 120, >= 3 studies)

- Westernized lifestyle = no: +3.0 to +3.5 points on 2 layers (n=326, 3 studies) -> candidate stratum-specific reference

## Age 65+ vs 18-39 (537 vs 2257 samples): features with q <= 0.05 and |shift| >= 15 percentile points

| layer | features tested | significant | shifted up in 65+ | shifted down |
|---|---:|---:|---:|---:|
| ko_humann | 4507 | 2469 | 193 | 391 |
| pathway_humann | 423 | 292 | 22 | 41 |
| taxonomy_family | 49 | 14 | 1 | 1 |
| taxonomy_genus | 119 | 43 | 2 | 5 |
| taxonomy_species | 253 | 69 | 14 | 5 |

Top 15 by p:
| layer          | feature_id       |   n_a |   n_b |   median_pct_a |   median_pct_b |   shift |   q |
|:---------------|:-----------------|------:|------:|---------------:|---------------:|--------:|----:|
| ko_humann      | K07321           |   472 |  2005 |         48.711 |         79.183 | -30.472 |   0 |
| ko_humann      | K07469           |   469 |  2019 |         48.056 |         77.949 | -29.893 |   0 |
| ko_humann      | K10546           |   452 |  1980 |         52.616 |         77.624 | -25.008 |   0 |
| pathway_humann | METH-ACETATE-PWY |   474 |  2030 |         46.905 |         79.149 | -32.243 |   0 |
| ko_humann      | K10118           |   476 |  2046 |         49.129 |         73.83  | -24.701 |   0 |
| ko_humann      | K13571           |   391 |  1733 |         40.848 |         67.124 | -26.276 |   0 |
| ko_humann      | K03484           |   479 |  2084 |         47.631 |         73.78  | -26.15  |   0 |
| ko_humann      | K01305           |   472 |  2034 |         51.439 |         76.654 | -25.215 |   0 |
| ko_humann      | K05813           |   440 |  1975 |         47.36  |         76.292 | -28.932 |   0 |
| ko_humann      | K10117           |   489 |  2112 |         45.482 |         68.084 | -22.602 |   0 |
| ko_humann      | K13787           |   326 |  1557 |         37.932 |         64.997 | -27.065 |   0 |
| ko_humann      | K07816           |   485 |  2088 |         46.186 |         71.5   | -25.314 |   0 |
| ko_humann      | K16925           |   367 |  1645 |         41.05  |         64.961 | -23.911 |   0 |
| ko_humann      | K00002           |   423 |  1905 |         50.978 |         77.368 | -26.39  |   0 |
| ko_humann      | K12510           |   484 |  2066 |         45.582 |         67.471 | -21.889 |   0 |

## Non-Western vs Western lifestyle (326 vs 6485 samples): features with q <= 0.05 and |shift| >= 15 percentile points

| layer | features tested | significant | shifted up in no | shifted down |
|---|---:|---:|---:|---:|
| ko_humann | 4423 | 3850 | 1137 | 1879 |
| pathway_humann | 420 | 362 | 71 | 161 |
| taxonomy_family | 44 | 34 | 11 | 13 |
| taxonomy_genus | 98 | 75 | 17 | 38 |
| taxonomy_species | 188 | 127 | 27 | 64 |

Top 15 by p:
| layer           | feature_id       |   n_a |   n_b |   median_pct_a |   median_pct_b |   shift |   q |
|:----------------|:-----------------|------:|------:|---------------:|---------------:|--------:|----:|
| pathway_humann  | TEICHOICACID-PWY |   317 |  5799 |         14.359 |         63.737 | -49.378 |   0 |
| ko_humann       | K10947           |   324 |  5832 |         11.261 |         61.049 | -49.788 |   0 |
| ko_humann       | K06871           |   324 |  5833 |         10.329 |         57.415 | -47.086 |   0 |
| ko_humann       | K17103           |   278 |  5787 |          7.36  |         61.353 | -53.993 |   0 |
| pathway_humann  | PWY-4984         |   300 |  5782 |          7.174 |         55.113 | -47.939 |   0 |
| ko_humann       | K01206           |   309 |  5809 |          9.943 |         57.497 | -47.554 |   0 |
| ko_humann       | K07098           |   317 |  5822 |         10.226 |         58.651 | -48.425 |   0 |
| ko_humann       | K02897           |   323 |  5827 |          8.783 |         52.843 | -44.061 |   0 |
| taxonomy_family | 815              |   312 |  6454 |          9.536 |         51.623 | -42.087 |   0 |
| ko_humann       | K02118           |   325 |  5837 |         13.361 |         59.194 | -45.833 |   0 |
| ko_humann       | K01914           |   325 |  5842 |         13.364 |         58.047 | -44.683 |   0 |
| ko_humann       | K04477           |   323 |  5819 |         14.853 |         63.688 | -48.836 |   0 |
| ko_humann       | K12267           |   310 |  5785 |          8.951 |         56.187 | -47.236 |   0 |
| pathway_humann  | CITRULBIO-PWY    |   321 |  5809 |         11.605 |         55.517 | -43.912 |   0 |
| pathway_humann  | PWY-6892         |   325 |  5848 |         16.374 |         64.141 | -47.768 |   0 |

## Female vs male (3696 vs 2506 samples): features with q <= 0.05 and |shift| >= 15 percentile points

| layer | features tested | significant | shifted up in female | shifted down |
|---|---:|---:|---:|---:|
| ko_humann | 4518 | 2733 | 2 | 9 |
| pathway_humann | 423 | 204 | 0 | 0 |
| taxonomy_family | 49 | 25 | 0 | 1 |
| taxonomy_genus | 119 | 45 | 0 | 1 |
| taxonomy_species | 254 | 92 | 0 | 0 |

Top 15 by p:
| layer           | feature_id   |   n_a |   n_b |   median_pct_a |   median_pct_b |   shift |     q |
|:----------------|:-------------|------:|------:|---------------:|---------------:|--------:|------:|
| ko_humann       | K09823       |   911 |   581 |         30.723 |         46.062 | -15.339 | 0     |
| ko_humann       | K10680       |   810 |   529 |         31.189 |         48.116 | -16.927 | 0     |
| ko_humann       | K13821       |  1491 |   848 |         80.287 |         63.947 |  16.34  | 0     |
| ko_humann       | K03676       |   836 |   488 |         29.354 |         44.616 | -15.262 | 0     |
| ko_humann       | K21814       |   641 |   280 |         77.753 |         61.287 |  16.466 | 0     |
| ko_humann       | K19230       |   572 |   385 |         33.784 |         48.894 | -15.11  | 0     |
| ko_humann       | K00844       |   466 |   299 |         30.016 |         46.931 | -16.915 | 0     |
| ko_humann       | K13630       |   419 |   323 |         33.197 |         49.355 | -16.157 | 0     |
| ko_humann       | K07067       |   203 |   192 |         36.623 |         57.384 | -20.761 | 0     |
| taxonomy_genus  | 848          |   362 |   278 |         37.707 |         53.529 | -15.822 | 0     |
| taxonomy_family | 203492       |   371 |   286 |         37.716 |         53.938 | -16.223 | 0     |
| ko_humann       | K00372       |   216 |   252 |         41.476 |         56.963 | -15.487 | 0.001 |
| ko_humann       | K03716       |   158 |   146 |         36.096 |         54.567 | -18.471 | 0.003 |

## BMI >= 30 vs 18.5-25 (1155 vs 3032 samples): features with q <= 0.05 and |shift| >= 15 percentile points

| layer | features tested | significant | shifted up in >=30 | shifted down |
|---|---:|---:|---:|---:|
| ko_humann | 4510 | 1448 | 27 | 23 |
| pathway_humann | 421 | 249 | 6 | 1 |
| taxonomy_family | 49 | 26 | 1 | 2 |
| taxonomy_genus | 119 | 59 | 4 | 7 |
| taxonomy_species | 254 | 99 | 7 | 11 |

Top 15 by p:
| layer            | feature_id     |   n_a |   n_b |   median_pct_a |   median_pct_b |   shift |   q |
|:-----------------|:---------------|------:|------:|---------------:|---------------:|--------:|----:|
| taxonomy_species | 216816         |   923 |  2507 |         42.267 |         61.871 | -19.603 |   0 |
| pathway_humann   | PHOSLIPSYN-PWY |   897 |  2858 |         59.13  |         41.778 |  17.352 |   0 |
| pathway_humann   | PWY4FS-7       |   897 |  2858 |         58.157 |         38.949 |  19.208 |   0 |
| pathway_humann   | PWY4FS-8       |   897 |  2858 |         58.157 |         38.949 |  19.208 |   0 |
| taxonomy_species | 1680           |   861 |  2167 |         43.166 |         59.828 | -16.662 |   0 |
| taxonomy_species | 40520          |  1061 |  2766 |         51.534 |         67.914 | -16.38  |   0 |
| pathway_humann   | PWY-5973       |   898 |  2863 |         55.791 |         40.732 |  15.059 |   0 |
| taxonomy_genus   | 572511         |  1150 |  2979 |         50.672 |         65.745 | -15.073 |   0 |
| taxonomy_species | 116085         |   963 |  2522 |         52.736 |         72.831 | -20.095 |   0 |
| taxonomy_genus   | 3570277        |   963 |  2522 |         55.697 |         72.748 | -17.051 |   0 |
| ko_humann        | K01494         |   760 |  2659 |         42.629 |         58.126 | -15.497 |   0 |
| taxonomy_genus   | 2172           |   386 |  1199 |         46.079 |         63.42  | -17.34  |   0 |
| taxonomy_species | 2173           |   386 |  1192 |         45.222 |         63.899 | -18.677 |   0 |
| taxonomy_family  | 2159           |   393 |  1220 |         45.43  |         64.942 | -19.512 |   0 |
| taxonomy_species | 1262824        |   527 |  1023 |         60.411 |         44.282 |  16.129 |   0 |
