# Figure scripts of the paper

These are the scripts that produce every figure of the manuscript and the HTML figure summary, copied verbatim from
the project workspace so that the plotting code, palette and lint rules are inspectable. They are **not** part of the
installable package and `pip install checkGM` does not install them.

| file | makes |
|---|---|
| `tokens.py` | palette, typography and the matplotlib print profile shared by every figure |
| `lint.py` | checks a rendered PDF/SVG pair: exact IEEE column width, 6 pt font floor, embedded fonts, palette-only colours, deterministic SVG |
| `make_figures.py` | renders every figure twice (determinism), runs `lint.py` on each, writes `LINT.json` and the figure README |
| `make_summary.py` | the stand-alone HTML summary of all figures with captions |
| `fig1_schematic.py` ... `fig14_sample_card.py` | one figure each; the manuscript figure numbers differ from the file numbers (`fig3` and `fig6` were dropped) |
| `tokens_validation.txt` | output of the palette contrast validator for `tokens.py` |

## Data dependency

The scripts read the evaluation results, baseline bundles and fonts of the project workspace, not the repository:
the project root is taken from the environment variable `CHECKGM_PROJECT` (in `tokens.py` and `make_figures.py`; the
individual figure scripts are being moved to the same variable and may still hold the development path at the time of
this copy). Under that root they expect `results/` (the validation tables, of which `docs/reports/` holds verbatim
copies), `work/` (per-sample score tables, not distributed), `refs/` (the bundles, distributed separately) and
`resources/fonts/`. Several scripts also import `VER`, `VTAG` and `bundle_A` from the workspace's own `checkgm/paths.py` (reached through
`sys.path.insert(0, ...)` in `tokens.py`); the installed package's `paths.py` defines none of these (0.8.2 removed its unused
`bundle_A`/`bundle_B`), so the scripts need the workspace copy of `paths.py`, not the packaged one. Re-running them therefore needs the full workspace; the copies here document how each figure was
drawn and which numbers it shows. Where a figure prints a number, the number is computed from the result tables listed
in `docs/validation.md`.
