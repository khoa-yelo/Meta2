# Paper figures (v2)

These are the scripts that draw figures 1-5 of the manuscript in `paper/`, and they are the only figure set the current
`main.tex` includes; the Python set in `scripts/viz/`, sixteen figures numbered up to 18, drew the superseded drafts and
is kept for provenance. The work is split in two so that the numbers and the drawing stay separate:

1. `prep.py` reads the results tables and the scored bundles and writes tidy CSVs to `figures/v2/data/`.
   Every number in a figure comes from a file named there; nothing is typed by hand.
2. `figures.R` draws them with ggplot2 + patchwork into `figures/v2/figN.pdf` (and a 300 dpi PNG preview).
   One theme and one palette throughout: blue below the healthy range (light blue expected but missing),
   orange above it, grey within; violet for checkGM, greys for comparators; teal/ochre for the two pipelines.

```sh
CHECKGM_PROJECT=$PWD REFMB_PROJECT=$PWD python scripts/viz2/prep.py
Rscript scripts/viz2/figures.R            # or: Rscript scripts/viz2/figures.R fig3 fig4
```

Both scripts expect the project workspace rather than this repository, and neither is portable as it stands. The
project root is hard-coded at the top of each (`P`), so they only run where that tree actually lives. `prep.py`
additionally imports two helpers from `scripts/viz/`, which take the root from the environment instead and fail with a
plain error when it names anything that is not a project directory. Running them elsewhere therefore means editing `P`
as well as exporting the variable. Which variable is read depends on the copy, the rename of the tool being unfinished
here: `tokens.py` as shipped reads `CHECKGM_PROJECT`, while `prep.py` still exports the older `REFMB_PROJECT`, which
is the name the workspace's own `tokens.py` reads. Setting both costs nothing and is the safe course.

The Python side needs the project environment (pandas, numpy, pyarrow). The R side loads ggplot2, patchwork, dplyr,
tidyr, readr and scales, and needs two more that no `library()` call names. `ragg` is called as `ragg::agg_png` by the
`save_fig` helper that every figure goes through, and the same helper writes its PDF with `cairo_pdf`, so a
cairo-enabled R is required too. Both are satisfied in the project's R environment, which is why nothing has ever
complained about them.

The manuscript itself is `paper/main.tex` (IEEEtran conference class, builds with `tectonic main.tex`).
