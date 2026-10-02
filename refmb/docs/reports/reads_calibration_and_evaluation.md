# Pipeline B v0.2 (read-based, MetaPhlAn 3 + HUMAnN 3 / cMD3) — calibration and evaluation — 2026-09-30 00:04 UTC

Fraction of assessed features outside the reference 2.5–97.5 band; expectation about 0.05, gate 0.08.

## Held-out healthy studies (never in the pool)

| layer            | study         |   median |   mean |   p90 |   share_over_gate |   n |
|:-----------------|:--------------|---------:|-------:|------:|------------------:|----:|
| ko_humann        | CosteaPI_2017 |    0.028 |  0.043 | 0.103 |             0.178 | 101 |
| ko_humann        | XieH_2016     |    0.036 |  0.055 | 0.117 |             0.22  | 173 |
| ko_humann        | YeZ_2018      |    0.054 |  0.108 | 0.253 |             0.372 |  43 |
| pathway_humann   | CosteaPI_2017 |    0.024 |  0.052 | 0.121 |             0.149 | 101 |
| pathway_humann   | XieH_2016     |    0.024 |  0.044 | 0.111 |             0.133 | 173 |
| pathway_humann   | YeZ_2018      |    0.053 |  0.107 | 0.319 |             0.395 |  43 |
| taxonomy_family  | CosteaPI_2017 |    0.04  |  0.062 | 0.148 |             0.248 | 101 |
| taxonomy_family  | XieH_2016     |    0.061 |  0.071 | 0.148 |             0.353 | 173 |
| taxonomy_family  | YeZ_2018      |    0.1   |  0.088 | 0.2   |             0.535 |  43 |
| taxonomy_genus   | CosteaPI_2017 |    0.045 |  0.053 | 0.094 |             0.198 | 101 |
| taxonomy_genus   | XieH_2016     |    0.057 |  0.068 | 0.118 |             0.283 | 173 |
| taxonomy_genus   | YeZ_2018      |    0.103 |  0.099 | 0.154 |             0.651 |  43 |
| taxonomy_species | CosteaPI_2017 |    0.045 |  0.052 | 0.091 |             0.139 | 101 |
| taxonomy_species | XieH_2016     |    0.065 |  0.072 | 0.115 |             0.329 | 173 |
| taxonomy_species | YeZ_2018      |    0.075 |  0.083 | 0.13  |             0.488 |  43 |

## In-pool leave-one-study-out (reference rebuilt without the study)

| layer            |   studies |   median_of_study_medians |   max_study_median |   studies_over_gate |
|:-----------------|----------:|--------------------------:|-------------------:|--------------------:|
| ko_humann        |        18 |                     0.03  |              0.102 |                   1 |
| pathway_humann   |        18 |                     0.026 |              0.08  |                   1 |
| taxonomy_family  |        19 |                     0.039 |              0.083 |                   2 |
| taxonomy_genus   |        19 |                     0.047 |              0.091 |                   2 |
| taxonomy_species |        19 |                     0.043 |              0.087 |                   1 |

Per study: `pipelineB_loso_all.tsv`.

## Case/control, cross-study leave-one-study-out (macro-mean AUROC; one-sided Wilcoxon, reference-relative greater)

| feature_set   | training    |   studies |   reference_relative |   raw_clr |   wins |   losses |   wilcoxon_p |
|:--------------|:------------|----------:|---------------------:|----------:|-------:|---------:|-------------:|
| all           | all_studies |        14 |                0.687 |     0.684 |      8 |        6 |        0.292 |
| all           | crc_only    |         9 |                0.758 |     0.752 |      7 |        2 |        0.082 |
| function      | all_studies |        14 |                0.676 |     0.679 |      5 |        9 |        0.787 |
| function      | crc_only    |         9 |                0.739 |     0.734 |      6 |        3 |        0.285 |
| taxonomy      | all_studies |        14 |                0.69  |     0.664 |      9 |        5 |        0.052 |
| taxonomy      | crc_only    |         9 |                0.769 |     0.715 |      8 |        1 |        0.004 |

## Per study (training on all other studies)

|                                |   ('all', 'raw_clr') |   ('all', 'reference_relative') |   ('function', 'raw_clr') |   ('function', 'reference_relative') |   ('taxonomy', 'raw_clr') |   ('taxonomy', 'reference_relative') |
|:-------------------------------|---------------------:|--------------------------------:|--------------------------:|-------------------------------------:|--------------------------:|-------------------------------------:|
| ('FengQ_2015', 'CRC')          |                0.674 |                           0.709 |                     0.701 |                                0.69  |                     0.738 |                                0.796 |
| ('GuptaA_2019', 'CRC')         |                0.766 |                           0.805 |                     0.791 |                                0.767 |                     0.622 |                                0.824 |
| ('HMP_2019', 'IBD')            |                0.745 |                           0.671 |                     0.727 |                                0.664 |                     0.503 |                                0.511 |
| ('HanniganGD_2017', 'CRC')     |                0.569 |                           0.565 |                     0.631 |                                0.621 |                     0.458 |                                0.512 |
| ('KarlssonFH_2013', 'T2D')     |                0.574 |                           0.554 |                     0.534 |                                0.537 |                     0.661 |                                0.65  |
| ('NagySzakalD_2017', 'ME/CFS') |                0.741 |                           0.752 |                     0.731 |                                0.706 |                     0.699 |                                0.727 |
| ('RubelMA_2020', 'STH')        |                0.5   |                           0.453 |                     0.479 |                                0.472 |                     0.481 |                                0.439 |
| ('ThomasAM_2019_a', 'CRC')     |                0.772 |                           0.769 |                     0.754 |                                0.75  |                     0.7   |                                0.728 |
| ('ThomasAM_2019_b', 'CRC')     |                0.768 |                           0.783 |                     0.716 |                                0.782 |                     0.658 |                                0.711 |
| ('VogtmannE_2016', 'CRC')      |                0.656 |                           0.64  |                     0.621 |                                0.627 |                     0.739 |                                0.715 |
| ('WirbelJ_2018', 'CRC')        |                0.765 |                           0.795 |                     0.75  |                                0.768 |                     0.833 |                                0.83  |
| ('YuJ_2015', 'CRC')            |                0.755 |                           0.777 |                     0.75  |                                0.728 |                     0.813 |                                0.838 |
| ('ZellerG_2014', 'CRC')        |                0.799 |                           0.839 |                     0.805 |                                0.84  |                     0.775 |                                0.796 |
| ('ZhuF_2020', 'schizofrenia')  |                0.494 |                           0.51  |                     0.509 |                                0.508 |                     0.613 |                                0.587 |

Studies: 14 with >= 20 cases and >= 20 controls. Function sets use the samples that have HUMAnN tables. Wall 215s.