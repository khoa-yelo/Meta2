# S6 calibration check — gut-assembly-adult-global-v0.8-lenient

Variant: detection floor 0.05, completeness correction off, excluded bands none, adults only True.

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


By quality band (batch vs banding diagnostic):

| layer           | quality_band   |   frac_outside |
|:----------------|:---------------|---------------:|
| ko_eggnog       | high           |          0.12  |
| ko_eggnog       | low            |          0.058 |
| ko_eggnog       | medium         |          0.14  |
| ko_kofam        | high           |          0.132 |
| ko_kofam        | low            |          0.064 |
| ko_kofam        | medium         |          0.138 |
| module          | high           |          0.007 |
| module          | low            |          0.004 |
| module          | medium         |          0.012 |
| pfam            | high           |          0.129 |
| pfam            | low            |          0.066 |
| pfam            | medium         |          0.135 |
| taxonomy_family | high           |          0.047 |
| taxonomy_family | low            |          0.063 |
| taxonomy_family | medium         |          0.044 |
| taxonomy_genus  | high           |          0.053 |
| taxonomy_genus  | low            |          0.065 |
| taxonomy_genus  | medium         |          0.05  |


Top excess features per layer (mapping diagnostic):

| feature_id   |   n |   frac_outside | layer           |
|:-------------|----:|---------------:|:----------------|
| PF06629      |  11 |          0.636 | pfam            |
| K10555       |  10 |          0.6   | ko_eggnog       |
| PF14175      |  10 |          0.6   | pfam            |
| K21883       |  10 |          0.6   | ko_eggnog       |
| PF06029      |  14 |          0.571 | pfam            |
| K08352       | 146 |          0.548 | ko_eggnog       |
| K16786       | 233 |          0.536 | ko_eggnog       |
| K17810       | 155 |          0.535 | ko_kofam        |
| PF09989      | 231 |          0.519 | pfam            |
| PF12892      | 185 |          0.519 | pfam            |
| PF12401      | 134 |          0.515 | pfam            |
| K03220       | 151 |          0.51  | ko_eggnog       |
| K00620       | 218 |          0.509 | ko_kofam        |
| K07739       | 124 |          0.508 | ko_kofam        |
| K01317       | 187 |          0.508 | ko_eggnog       |
| K06442       | 225 |          0.507 | ko_eggnog       |
| K09157       | 228 |          0.504 | ko_kofam        |
| 2005386      | 135 |          0.504 | taxonomy_genus  |
| 644652       | 161 |          0.503 | taxonomy_genus  |
| K13613       |  10 |          0.5   | ko_eggnog       |
| K13640       | 134 |          0.5   | ko_kofam        |
| K05539       |  10 |          0.5   | ko_kofam        |
| K16787       | 233 |          0.498 | ko_eggnog       |
| K16254       | 127 |          0.496 | ko_eggnog       |
| PF01910      | 142 |          0.493 | pfam            |
| 395331       |  30 |          0.467 | taxonomy_genus  |
| PF12029      | 133 |          0.466 | pfam            |
| PF06050      | 231 |          0.459 | pfam            |
| PF08353      | 192 |          0.458 | pfam            |
| K06442       | 219 |          0.457 | ko_kofam        |
| K18817       |  11 |          0.455 | ko_kofam        |
| K08151       |  11 |          0.455 | ko_kofam        |
| K03328       | 226 |          0.447 | ko_kofam        |
| 1926677      | 155 |          0.374 | taxonomy_genus  |
| 1091         |  38 |          0.368 | taxonomy_genus  |
| 1473205      | 140 |          0.364 | taxonomy_genus  |
| 1392997      |  11 |          0.364 | taxonomy_genus  |
| 1198         |  17 |          0.353 | taxonomy_genus  |
| 3046411      |  17 |          0.353 | taxonomy_family |
| 82800        |  40 |          0.35  | taxonomy_genus  |
| 3151707      | 183 |          0.344 | taxonomy_family |
| 1427376      | 134 |          0.336 | taxonomy_genus  |
| 1392996      |  13 |          0.308 | taxonomy_family |
| 3031659      |  34 |          0.294 | taxonomy_family |
| 85031        |  14 |          0.286 | taxonomy_family |
| 84107        | 221 |          0.276 | taxonomy_family |
| 1892254      |  28 |          0.25  | taxonomy_family |
| 3031626      |  23 |          0.217 | taxonomy_family |
| 1643826      | 219 |          0.21  | taxonomy_family |
| 89377        |  20 |          0.2   | taxonomy_family |
| M00866       | 229 |          0.083 | module          |
| M00064       | 231 |          0.082 | module          |
| M00060       | 229 |          0.074 | module          |
| M00620       | 233 |          0.06  | module          |
| M00567       | 186 |          0.054 | module          |
| M00572       | 233 |          0.052 | module          |
| M00144       | 233 |          0.047 | module          |
| M00357       | 233 |          0.047 | module          |
| M00374       | 233 |          0.039 | module          |
| M00167       | 233 |          0.039 | module          |


**Gate: FAIL — STOP, diagnose, do not widen intervals**
