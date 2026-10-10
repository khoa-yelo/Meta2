# Does the reference range transfer across the baseline's own populations? (read pipeline)

Reference adults only, each scored against a baseline rebuilt without its own study. 6810 samples, 22 studies, 17 countries with at least 40 samples. 5% outside is expected by construction.

Statistic: mean per-sample fraction outside within a study, then the median over a country's studies. The per-sample median is not usable here: with a median of 24 assessed families a 5% expectation is about one feature, so over half the samples of several countries have exactly zero outside and the median is 0.0 regardless of fit.

| country   |   studies |    n |   feats |   family |   genus |   species |   family_spread |   genus_spread |   species_spread |
|:----------|----------:|-----:|--------:|---------:|--------:|----------:|----------------:|---------------:|-----------------:|
| MDG       |         1 |  110 |    24   |      9.9 |     9.1 |       7.3 |           nan   |          nan   |            nan   |
| MNG       |         1 |  109 |    25   |      9.7 |    10.2 |       9.3 |           nan   |          nan   |            nan   |
| IND       |         1 |   88 |    16   |      8.6 |     7.4 |       4.3 |           nan   |          nan   |            nan   |
| KAZ       |         1 |   83 |    28   |      6.3 |     5.5 |       5.4 |           nan   |          nan   |            nan   |
| JPN       |         1 |  250 |    24   |      6.2 |     7.9 |       9   |           nan   |          nan   |            nan   |
| CHN       |         4 |  500 |    24.5 |      6.1 |     7.3 |       7.7 |             3.2 |            3.6 |              1.7 |
| IRL       |         1 |  117 |    22   |      5.7 |     4.2 |       3.8 |           nan   |          nan   |            nan   |
| FJI       |         1 |  152 |    31   |      5.6 |     6   |       5.8 |           nan   |          nan   |            nan   |
| NLD       |         2 | 1604 |    24.5 |      5.6 |     6.4 |       5.8 |             3.4 |            4   |              3.7 |
| SWE       |         1 |   99 |    27   |      5.2 |     5   |       5.5 |           nan   |          nan   |            nan   |
| ESP       |         1 |   59 |    27   |      4.7 |     5   |       4.8 |           nan   |          nan   |            nan   |
| DNK       |         3 |  641 |    25   |      4.4 |     3.6 |       3.9 |             4   |            3   |              3.3 |
| DEU       |         2 |  290 |    24   |      4.4 |     4.1 |       4.1 |             2.2 |            1.1 |              0.5 |
| ISR       |         1 |  888 |    26   |      3.9 |     3.6 |       4   |           nan   |          nan   |            nan   |
| GBR       |         3 | 1343 |    25   |      3.7 |     3.8 |       4.6 |             4.1 |            3.5 |              3.7 |
| USA       |         2 |  304 |    21.5 |      3.6 |     3.8 |       4.4 |             0.5 |            0.7 |              1   |
| FRA       |         1 |  173 |    24   |      3.2 |     2.5 |       3.1 |           nan   |          nan   |            nan   |

## Is it the population or the study?

Country and study coincide wherever a country contributes one study. Where a country contributes two or more, the spread between its own studies measures what a protocol alone can do:

| country   |   family |   genus |   species |
|:----------|---------:|--------:|----------:|
| CHN       |      3.2 |     3.6 |       1.7 |
| DEU       |      2.2 |     1.1 |       0.5 |
| DNK       |      4   |     3   |       3.3 |
| GBR       |      4.1 |     3.5 |       3.7 |
| NLD       |      3.4 |     4   |       3.7 |
| USA       |      0.5 |     0.7 |       1   |

Judged on standard deviations rather than ranges, which depend on how many groups are compared:

| layer   |   sd_between_countries |   sd_within_country |   industrialised_median |   top3_excess_in_within_SDs |
|:--------|-----------------------:|--------------------:|------------------------:|----------------------------:|
| family  |                   2.04 |                1.85 |                     3.9 |                         3   |
| genus   |                   2.15 |                1.82 |                     4   |                         2.8 |
| species |                   1.82 |                1.7  |                     4   |                         2.7 |

Study-to-study variation within a country is almost as large as variation between countries (1.85 against 2.04 points on families), so most of the spread in the table above cannot be attributed to population rather than protocol. The exception is the top of it: Madagascar, Mongolia and India sit about 3 within-country standard deviations above the industrialised median, which protocol alone is not observed to do.


## Does the excess replicate across layers?

Excess over the industrialised median (family 4.4%, genus 3.9%, species 4.2%), in points:

| country   |   assessed_families |   family |   genus |   species | replicates        |
|:----------|--------------------:|---------:|--------:|----------:|:------------------|
| MDG       |                24   |      5.5 |     5.2 |       3   | all three layers  |
| MNG       |                25   |      5.3 |     6.2 |       5.1 | all three layers  |
| IND       |                16   |      4.2 |     3.5 |       0   | one layer only    |
| KAZ       |                28   |      1.9 |     1.6 |       1.2 | one layer only    |
| JPN       |                24   |      1.8 |     4   |       4.8 | genus and species |
| CHN       |                24.5 |      1.7 |     3.4 |       3.4 | genus and species |
| IRL       |                22   |      1.3 |     0.3 |      -0.5 | one layer only    |
| FJI       |                31   |      1.2 |     2.1 |       1.5 | one layer only    |

Only Madagascar and Mongolia are elevated on all three layers. India is elevated on families alone, on 16 assessed features where one feature is six points, and sits at the industrialised median on species; it is not counted. Japan and China run the other way, at the industrialised level on families and two to five points above it on genus and species, which is the more surprising result of the two: a reference range can transfer at one rank and not at another.
