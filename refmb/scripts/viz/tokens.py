"""Visual tokens shared by every figure (plan §9.2) and matplotlib print-profile setup (plan §9.1).

Palette choices follow the dataviz method: colour is assigned by job (diverging for deviation direction,
sequential for density, a fixed-order categorical set for the four feature sets), validated with the
skill's validator (see tokens_validation.txt next to this file), text in ink tokens only, texture as the
CVD/greyscale backup, no gradients.
"""
import os
import matplotlib
matplotlib.use("Agg")
from matplotlib import font_manager, rcParams

# Project root: REFMB_PROJECT in the environment; when unset, the directory two levels above this file (the scripts live in
# <project>/scripts/viz/). Every figure script reads paths from here, so the scripts can be copied into the public
# repository unchanged; there the fallback is not a project directory and the check below fails with a plain message.
_DEFAULT_P = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
P = os.environ.get("REFMB_PROJECT") or _DEFAULT_P
if not os.path.isdir(P) or not os.path.isdir(os.path.join(P, "resources")) or not os.path.isdir(os.path.join(P, "refs")):
    raise FileNotFoundError(
        f"refmb project root not found: REFMB_PROJECT is {os.environ.get('REFMB_PROJECT')!r} and the fallback {_DEFAULT_P!r} is not a project directory "
        "(a project root holds resources/, refs/, work/ and results/). Set REFMB_PROJECT=/path/to/project before running the figure scripts.")
# curatedMetagenomicData 3 sample metadata (used by fig2 panel c only); lives beside the project in the development checkout
CMD3_META = os.environ.get("REFMB_CMD3_META", os.path.join(os.path.dirname(P), "cmd3", "cmd3_release_sampleMetadata.csv"))

FONT_DIR = os.path.join(P, "resources", "fonts")
for f in sorted(os.listdir(FONT_DIR)):
    if f.endswith(".ttf"):
        font_manager.fontManager.addfont(os.path.join(FONT_DIR, f))
FONT = "DejaVu Sans"

# --- deviation scale: 5 diverging steps, blue (low) -> neutral -> orange (high); lightness differs as well as hue
DEV = {"low_deep": "#1c5cab", "low_light": "#5598e7", "neutral": "#a9a7a1", "high_light": "#e07f3b", "high_deep": "#b8460f"}
DEV_ORDER = ["low_deep", "low_light", "neutral", "high_light", "high_deep"]
# categorical, fixed order (purple, green, magenta, yellow), used only for the four feature sets of the classification benchmark
# (fig8); direct-labelled. Blue and orange are deliberately absent: they carry the deviation direction. Validated with the dataviz
# validator (tokens_validation.txt): all checks pass within the set; the purple also clears the normal-vision floor against the
# deviation blue (ΔE 16.5), which it never shares a figure with.
CAT = ["#8a3ab9", "#008300", "#e87ba4", "#eda100"]
CAT_LABELS = {"reference_relative": "reference-relative", "raw_clr": "raw CLR", "alpha_diversity": "alpha diversity", "health_index": "GMHI (genus approx.)"}
# sequential (density, hexbin): one blue hue light -> dark
SEQ = ["#86b6ef", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]   # starts at ramp step 250 so the lightest bin clears 2:1 on white
# ink & chrome (print)
INK = {"primary": "#0b0b0b", "secondary": "#52514e", "muted": "#898781", "grid": "#e1e0d9", "axis": "#c3c2b7", "surface": "#ffffff"}
# categorical state marks (plan §9.2): expected_but_missing = hollow ring; not_assessable = hatch
MISSING_RING = {"facecolor": "none", "edgecolor": INK["secondary"], "linewidth": 0.6}
HATCH = "////"

# ---------------------------------------------------------------------------------------------------------------------
# Colour conventions (one meaning per colour across the whole figure set; figures/README.md repeats this table)
#   DEV blue / orange ........ direction of a deviation from the healthy range ONLY: blue = below (low), orange = above
#                              (high); deep step = beyond the 1st / 99th percentile or a shift >= 20 points, light step =
#                              the milder class; DEV neutral grey = within the range. Used in fig1 (report rows), fig5,
#                              fig10, fig12, fig13, fig14 (orange diamond = above the healthy 90th percentile), and for
#                              the "cohort above the pre-set 8 % limit" mark in fig4 (an exceedance, i.e. a high call).
#   CAT (4 slots) ............ the four feature sets of the classification benchmark only (fig8): percentiles (purple), raw
#                              CLR (green), alpha diversity (magenta), GMHI approximation (yellow). Fixed order, direct-labelled;
#                              no blue or orange, so a feature set is never mistaken for a deviation direction.
#   SEQ blues ................ an ordinal magnitude: sample density per hex (fig7) and the assembly-size / read-depth
#                              class of a study's samples (fig2 a, d). Never a category.
#   INK greys ................ everything else: the two pipelines (filled light-grey box / solid bar = assembly pipeline,
#                              outlined white box / outlined bar = read pipeline), layers (marker shape, one hue:
#                              circle = family, square = KO), read layout (circle = paired, triangle = single-end),
#                              medians, quantile bars, reference lines. Pipelines, layers and layouts are therefore
#                              never carried by hue, so blue and orange keep their single meaning.
# ---------------------------------------------------------------------------------------------------------------------
PIPE = {  # the pipeline pair, by fill and outline (not hue)
    "assembly": {"facecolor": INK["grid"], "edgecolor": INK["primary"], "linestyle": "-"},
    "read": {"facecolor": INK["surface"], "edgecolor": INK["primary"], "linestyle": (0, (3, 1.5))},
}
BAR_ASSEMBLY = {"color": INK["secondary"], "linewidth": 0}                                   # filled bar
BAR_READ = {"facecolor": INK["surface"], "edgecolor": INK["primary"], "linewidth": 0.7}      # outlined bar
LAYER_MARK = {"family": "o", "KO": "s"}                                                       # layers by shape, one hue
LAYOUT_MARK = {"paired": "o", "single": "^"}                                                  # read layout by shape

# --- print profile (IEEE): widths in inches; body 7 pt, floor 6 pt
SINGLE_IN = 88.9 / 25.4
DOUBLE_IN = 181.9 / 25.4
PT_BODY, PT_MIN, PT_TITLE = 7, 6, 8


def print_profile():
    rcParams.update({
        "font.family": FONT, "font.size": PT_BODY, "axes.titlesize": PT_TITLE, "axes.labelsize": PT_BODY,
        "xtick.labelsize": PT_MIN, "ytick.labelsize": PT_MIN, "legend.fontsize": PT_MIN,
        "axes.edgecolor": INK["axis"], "axes.linewidth": 0.5, "xtick.color": INK["muted"], "ytick.color": INK["muted"],
        "axes.labelcolor": INK["secondary"], "text.color": INK["primary"], "grid.color": INK["grid"], "grid.linewidth": 0.4,
        "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 1.0, "patch.linewidth": 0.5,
        "pdf.fonttype": 42, "svg.fonttype": "none", "svg.hashsalt": "refmb", "figure.dpi": 300, "savefig.dpi": 300,
        "figure.facecolor": INK["surface"], "axes.facecolor": INK["surface"], "savefig.facecolor": INK["surface"],
        "path.simplify": False,
    })


def save(fig, stem, out_dir):
    """Deterministic SVG + PDF at exact width; metadata that would change between runs is removed."""
    os.makedirs(out_dir, exist_ok=True)
    meta_svg = {"Date": None, "Creator": "refmb viz", "Title": os.path.basename(stem)}
    meta_pdf = {"CreationDate": None, "ModDate": None, "Creator": "refmb viz", "Producer": "matplotlib", "Title": os.path.basename(stem)}
    fig.savefig(os.path.join(out_dir, stem + ".svg"), format="svg", metadata=meta_svg, bbox_inches=None)
    fig.savefig(os.path.join(out_dir, stem + ".pdf"), format="pdf", metadata=meta_pdf, bbox_inches=None)
    os.makedirs(os.path.join(out_dir, "_preview"), exist_ok=True)   # PNG preview for visual inspection only; not a deliverable
    fig.savefig(os.path.join(out_dir, "_preview", stem + ".png"), format="png", dpi=200, bbox_inches=None)
