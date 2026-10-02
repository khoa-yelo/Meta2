# S7 scores — gut, bundle gut-assembly-adult-global-v0.8-lenient

in 8,742 = scored 8,742 + rejected 0  (OK)

## Median fraction of assessed features outside p2.5–p97.5 (raw), by layer and health label

| layer           |   (blank) |   Diseased |   Healthy |   Other |   Treatment |   cMD3 case |   cMD3 control |
|:----------------|----------:|-----------:|----------:|--------:|------------:|------------:|---------------:|
| ko_eggnog       |     0.033 |      0.055 |     0.04  |   0.036 |       0.044 |       0.058 |          0.052 |
| ko_kofam        |     0.035 |      0.055 |     0.04  |   0.038 |       0.051 |       0.055 |          0.05  |
| module          |     0     |      0     |     0     |   0     |       0.005 |       0     |          0     |
| pfam            |     0.035 |      0.055 |     0.043 |   0.039 |       0.046 |       0.065 |          0.058 |
| taxonomy_family |     0.051 |      0.051 |     0.051 |   0.045 |       0.058 |       0.078 |          0.08  |
| taxonomy_genus  |     0.052 |      0.054 |     0.051 |   0.045 |       0.066 |       0.069 |          0.065 |

## Deviation calls per sample after FDR (median), by layer

| layer           |   n_assessed |   n_low |   n_high |   n_missing |   n_not_assessable |   expected_false_positives_at_alpha |
|:----------------|-------------:|--------:|---------:|------------:|-------------------:|------------------------------------:|
| ko_eggnog       |         3384 |       0 |        0 |           2 |             1637   |                               169.2 |
| ko_kofam        |         2107 |       0 |        0 |           2 |             1085   |                               105.4 |
| module          |          206 |       0 |        0 |           1 |                0   |                                10.3 |
| pfam            |         3554 |       0 |        0 |           3 |             1457.5 |                               177.7 |
| taxonomy_family |          160 |       0 |        0 |           6 |                0   |                                 8   |
| taxonomy_genus  |          405 |       0 |        0 |          17 |                0   |                                20.2 |

## Feature-level FDR
With ~1,200 reference samples the smallest attainable two-sided p per feature is ~1/1,221, so BH across ~3,000 features per sample cannot mark individual features significant; `call_fdr` is therefore almost always `within_after_fdr`. Deviation calls are reported as outside-the-95%-band with the expected false-positive count (0.05 × n_assessed), and significance is assessed at the sample level (binomial excess test, BH across samples within layer).

## Samples with a significant excess of outside features (q ≤ 0.05), by layer and health label

| layer           |   (blank) |   Diseased |   Healthy |   Other |   Treatment |   cMD3 case |   cMD3 control |
|:----------------|----------:|-----------:|----------:|--------:|------------:|------------:|---------------:|
| ko_eggnog       |     0.315 |      0.479 |     0.391 |   0.298 |       0.394 |       0.5   |          0.467 |
| ko_kofam        |     0.282 |      0.475 |     0.375 |   0.294 |       0.394 |       0.46  |          0.419 |
| module          |     0     |      0.005 |     0.003 |   0.021 |       0     |       0.004 |          0.022 |
| pfam            |     0.321 |      0.488 |     0.407 |   0.337 |       0.424 |       0.54  |          0.504 |
| taxonomy_family |     0.114 |      0.108 |     0.212 |   0.082 |       0.121 |       0.25  |          0.265 |
| taxonomy_genus  |     0.204 |      0.213 |     0.304 |   0.17  |       0.273 |       0.339 |          0.319 |

## In-reference vs query samples (eggNOG KO, median frac outside raw)

| in_reference_pool   |   frac_outside_raw |
|:--------------------|-------------------:|
| False               |              0.051 |
| True                |              0.032 |

## Landscape core-distance percentile by health label (family layer)

| group        |   core_distance_pct |
|:-------------|--------------------:|
| (blank)      |                42.4 |
| Diseased     |                44.2 |
| Healthy      |                55.1 |
| Other        |                59.2 |
| Treatment    |                23.8 |
| cMD3 case    |                77.4 |
| cMD3 control |                81   |
