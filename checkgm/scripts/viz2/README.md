# Paper figures (v2)

Figures 1-5 of the two-column paper. Two steps, so the numbers and the drawing stay separate:

1. `prep.py` reads the results tables and the scored bundles and writes tidy CSVs to `figures/v2/data/`.
   Every number in a figure comes from a file named there; nothing is typed by hand.
   Needs the project Python env and `CHECKGM_PROJECT` pointing at the project root.
2. `figures.R` draws them with ggplot2 + patchwork into `figures/v2/figN.pdf` (and a 300 dpi PNG preview).
   One theme and one palette throughout: blue below the healthy range (light blue expected but missing),
   orange above it, grey within; violet for checkGM, greys for comparators; teal/ochre for the two pipelines.

```sh
CHECKGM_PROJECT=$PWD python scripts/viz2/prep.py
Rscript scripts/viz2/figures.R            # or: Rscript scripts/viz2/figures.R fig3 fig4
```

The paper itself is `paper/main.tex` (IEEEtran conference class, builds with tectonic).
