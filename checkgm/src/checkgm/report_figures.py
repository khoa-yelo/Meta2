"""checkgm.report_figures — the per-sample laboratory-report figure that `assess` emits.

One figure per sample, read like a clinical laboratory report: each row is one feature, the grey bar is the reference
range (2.5th-97.5th percentile of the baseline), the darker inner bar its middle half (25th-75th), the tick the
reference median, the thin line the 1st-99th percentile, and the dot is this sample, coloured by its call. The design is
the paper's figure 17 (scripts/viz/fig17_range_report.py in the research workspace); this is the shippable form of it,
reading nothing but an assess output directory and the bundle that directory names.

matplotlib is optional: it is imported inside the drawing code and `available()` answers False rather than raising when
it is absent, so `assess` keeps working in an install without it.
"""
from __future__ import annotations

import contextlib
import importlib.util
import math
import os
import re
import warnings

import numpy as np
import pandas as pd

from checkgm.paths import resolve_bundle
from checkgm.score import LAYER_SPEC, Bundle

# Visual tokens of the paper's figure set (scripts/viz/tokens.py), copied rather than imported: that module is research
# code which raises unless the project workspace is present.
DOUBLE_IN = 181.9 / 25.4                       # IEEE double-column width in inches
PT_BODY, PT_MIN, PT_TITLE = 7, 6, 8
DEV = {"low_deep": "#1c5cab", "low_light": "#5598e7", "neutral": "#a9a7a1", "high_light": "#e07f3b", "high_deep": "#b8460f"}
INK = {"primary": "#0b0b0b", "secondary": "#52514e", "muted": "#898781", "grid": "#e1e0d9", "axis": "#c3c2b7", "surface": "#ffffff"}

FORMATS = ("png", "pdf", "svg")
ROW_BUDGET = 20          # rows a reader can take in at a glance; the per-panel cap is this divided by the panels drawn
ROW_IN = 0.22            # vertical inches per feature row, so the figure grows with the number of rows rather than squeezing them
HEAD_IN, TITLE_IN, FOOT_IN = 0.62, 0.30, 0.62   # inches for the figure heading, a panel heading, and a panel's axis with its caption
LEFT, WIDTH, MARGIN = 0.27, 0.43, 0.02          # the axes leave the left quarter to feature names and the right third to the readings

# The two panels, each the first of its preference list that the sample has scored rows on: one taxonomic layer (named,
# well populated, the paper's choice) and one functional layer. Showing every layer of an assembly bundle would leave
# three rows per layer, which is worse than showing two layers properly.
TAXON_LAYERS = ("taxonomy_family", "taxonomy_genus", "taxonomy_species")
FUNCTION_LAYERS = ("ko_eggnog", "ko_kofam", "ko_humann", "pfam", "module", "pathway_humann")
LAYER_LABEL = {"taxonomy_family": "bacterial families", "taxonomy_genus": "bacterial genera", "taxonomy_species": "bacterial species",
               "ko_eggnog": "gene families (KEGG orthologs, eggNOG)", "ko_kofam": "gene families (KEGG orthologs, KOfam)",
               "ko_humann": "gene families (KEGG orthologs, HUMAnN)", "pfam": "protein families (Pfam)",
               "module": "KEGG modules", "pathway_humann": "pathways (MetaCyc, HUMAnN)"}

# How a layer's measurement becomes a distance from the typical reference sample. The value column is score.LAYER_SPEC's,
# so a layer added there is handled here as soon as its unit is one of these three.
MODE_OF_VALUE = {"clr": "log2_diff",            # centred log-ratio: a difference in natural-log units, shown in log2
                 "copies_per_genome": "log2_ratio",   # a positive rate: shown as a log2 fold difference from the median
                 "completeness": "linear"}      # a bounded fraction: a fold difference would be meaningless, so plot the difference
X_LABEL = {"log2_diff": "difference from the typical reference sample (fold difference in the centred log-ratio)",
           "log2_ratio": "fold difference from the reference median (copies per genome)",
           "linear": "difference from the reference median (completeness)"}

# The percentiles the row is drawn from. score.score() indexes p2_5 and p97_5 directly, so those two and the median are
# present in any bundle that could have produced the input; the rest degrade gracefully (see _quantiles).
PCOL = {"lo99": "p1", "lo": "p2_5", "q1": "p25", "med": "p50", "q3": "p75", "hi": "p97_5", "hi99": "p99"}
ESSENTIAL = ("lo", "med", "hi")


def available() -> bool:
    """Whether figures can be drawn here, i.e. whether matplotlib is importable. Never raises."""
    try:
        return importlib.util.find_spec("matplotlib") is not None
    except (ImportError, ValueError):   # a broken or shadowed install: treat as absent, the caller only wants to skip
        return False


@contextlib.contextmanager
def _quiet():
    """Swallow library warnings and numpy floating-point notices for the duration.

    The one-line `checkgm: ...` stderr contract is tested, and drawing is a rich enough source of third-party warnings
    (font caches, unwritable MPLCONFIGDIR, empty slices) that none of them may be allowed to leak out of here."""
    with warnings.catch_warnings(), np.errstate(all="ignore"):
        warnings.simplefilter("ignore")
        yield


def _pyplot():
    """pyplot and the print profile, imported on first use. The profile is returned rather than installed: rcParams is
    process-global and permanent, so updating it would restyle every later plot of a process that embeds checkgm. The
    caller applies it with plt.rc_context. An absent matplotlib raises ValueError, not ModuleNotFoundError, because
    only the former reaches cli.main's one-line "checkgm: " handler; available() exists so the CLI can skip first."""
    try:
        import matplotlib
    except ImportError:
        raise ValueError("drawing figures needs matplotlib, which is not installed; "
                         "install it with 'pip install checkgm[figures]' or pass --no-figures") from None
    matplotlib.use("Agg", force=False)   # headless by default; force=False leaves an embedding process's own backend alone
    import matplotlib.pyplot as plt
    style = {
        "font.size": PT_BODY, "axes.titlesize": PT_TITLE, "axes.labelsize": PT_MIN, "xtick.labelsize": PT_MIN,
        "ytick.labelsize": PT_MIN, "axes.edgecolor": INK["axis"], "axes.linewidth": 0.5, "xtick.color": INK["muted"],
        "ytick.color": INK["muted"], "axes.labelcolor": INK["secondary"], "text.color": INK["primary"],
        "grid.color": INK["grid"], "grid.linewidth": 0.4, "axes.spines.top": False, "axes.spines.right": False,
        "figure.facecolor": INK["surface"], "axes.facecolor": INK["surface"], "savefig.facecolor": INK["surface"],
        "savefig.dpi": 300, "pdf.fonttype": 42, "svg.fonttype": "none",
    }
    return plt, style


# ------------------------------------------------------------------------------------------------- reading the input

def _read_assessed(assessed_dir: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """The scores and summaries tables of an assess output directory."""
    sp, mp = os.path.join(assessed_dir, "scores.parquet"), os.path.join(assessed_dir, "summaries.tsv")
    if not os.path.isdir(assessed_dir):
        raise FileNotFoundError(f"assess output directory not found: {assessed_dir}")
    if not os.path.isfile(sp):
        raise FileNotFoundError(f"{assessed_dir} holds no scores.parquet, so it is not an assess output directory "
                                f"(pass the directory 'checkgm assess' wrote its report to)")
    scores = pd.read_parquet(sp)
    try:   # an assess run that rejected every sample writes a summaries.tsv of one newline, which read_csv rejects
        summ = pd.read_csv(mp, sep="\t") if os.path.isfile(mp) else pd.DataFrame()
    except pd.errors.EmptyDataError:
        summ = pd.DataFrame()
    return scores, summ


def _bundle_of(assessed_dir: str, summ: pd.DataFrame) -> Bundle:
    """The bundle the samples were assessed against, named by summaries.tsv and looked up like any bundle argument.

    The figure needs the reference percentiles themselves, which only the bundle holds; the assess output records the
    bundle's id, so the id is resolved through the ordinary CHECKGM_BUNDLES search and the failure is the familiar one."""
    ids = sorted({str(b) for b in summ.get("bundle_id", pd.Series(dtype=object)).dropna().unique()})
    if not ids:
        raise ValueError(f"{assessed_dir} does not name the bundle it was assessed against (summaries.tsv has no bundle_id), "
                         f"so the reference ranges cannot be read; re-run 'checkgm assess' with this release")
    if len(ids) > 1:
        raise ValueError(f"{assessed_dir} mixes {len(ids)} bundles ({', '.join(ids)}); figures are drawn against one "
                         f"reference, so assess each bundle into its own directory")
    try:
        return Bundle(resolve_bundle(ids[0]))
    except FileNotFoundError as e:
        raise FileNotFoundError(f"the figures need the bundle the samples were assessed against, and {e}") from None


# TODO: gene, protein and module rows are labelled by identifier (K01470, PF00001, M00001). A bundle ships no name table
# for them -- the paper's figure took KO names from a project-only KEGG map file -- so no name can be shown without
# inventing one. If a later bundle format carries one, label these rows with it as the taxonomy rows are labelled.
def _taxon_names(bundle: Bundle, taxids) -> dict:
    """{taxid: scientific name} from the bundle's NCBI backbone, or {} when the bundle ships none.

    The same lookup as the text report's, repeated here rather than imported from cli, which imports this module."""
    p = os.path.join(bundle.path, "normalization", "backbone", "ncbi_taxonomy.parquet")
    want = {str(t) for t in taxids}
    if not want or not os.path.isfile(p):
        return {}
    tax = pd.read_parquet(p, columns=["taxid", "name"])
    tax = tax[tax["taxid"].astype(str).isin(want)]
    return dict(zip(tax["taxid"].astype(str), tax["name"]))


# ------------------------------------------------------------------------------------------------- feature selection

def _prevalence(ref: pd.DataFrame, band) -> pd.Series:
    """Reference prevalence per feature: within the sample's own quality band when the bundle reports it, else the
    highest band prevalence (the same quantity the expected-but-missing rule is defined on)."""
    col = f"prevalence_band_{band}"
    if col in ref.columns:
        return ref[col].astype(float)
    cols = [c for c in ref.columns if c.startswith("prevalence_band_")]
    return ref[cols].astype(float).max(axis=1) if cols else pd.Series(0.0, index=ref.index)


def _quotas(cap: int) -> list[tuple[str, int]]:
    """Per-pool row quotas for a panel of `cap` rows: 30% most-deviant low, 30% high, 20% missing, 20% context."""
    low = max(1, cap * 3 // 10)
    miss = max(1, (cap - 2 * low) // 2)
    return [("low", low), ("high", low), ("missing", miss), ("within", max(0, cap - 2 * low - miss))]


def _outside_distance(g: pd.DataFrame, ref: pd.DataFrame) -> pd.Series:
    """How far each value sits outside its reference range, in the layer's own units, for ranking the extremes.

    The percentile cannot do this job past the grid: interp_percentile clamps to 0.5 / 99.5, so the worst twenty
    features tie at one number. Distance beyond the nearer reference bound keeps separating them. Returns 0 for a
    feature inside the range and NaN where a bound is missing, both of which sort last.
    """
    if "value" not in g.columns or not {"p2_5", "p97_5"} <= set(ref.columns):
        return pd.Series(np.nan, index=g.index)   # the bounds live in the bundle, not in scores.parquet
    v = pd.to_numeric(g["value"], errors="coerce")
    lo = pd.to_numeric(ref["p2_5"].reindex(g.index), errors="coerce")
    hi = pd.to_numeric(ref["p97_5"].reindex(g.index), errors="coerce")
    return pd.concat([lo - v, v - hi], axis=1).max(axis=1).clip(lower=0)


def _deviant(sub: pd.DataFrame, d: pd.Series, low_first: bool) -> list[str]:
    """One directional pool, worst first. Not written with DataFrame.assign: assigning a full-length Series to an
    empty frame does not align down to the empty index, it adopts the Series' index, so a sample with no call in this
    direction would come back holding every feature of the layer."""
    if sub.empty:
        return []
    k = pd.DataFrame({"_d": d.reindex(sub.index), "_p": pd.to_numeric(sub["percentile"], errors="coerce")}, index=sub.index)
    return list(k.sort_values(["_d", "_p"], ascending=[False, low_first], kind="mergesort").index)


def _select(g: pd.DataFrame, ref: pd.DataFrame, prev: pd.Series, cap: int) -> list[str]:
    """Which features of one layer the panel shows, by the rule stated in write_figures's docstring.

    g is the sample's score rows for the layer indexed by feature_id; prev is reference prevalence. Pools are filled to
    their quota in a fixed order and unused quota is passed down the same order, so a sample with no high calls shows
    more low ones instead of a panel of context rows. Sorting is stable over a feature_id-sorted frame, so equal
    percentiles and equal prevalences resolve the same way on every run."""
    gs = g.sort_index()
    pr = prev.reindex(gs.index)
    # the two directional pools rank by distance outside the reference range, not by percentile:
    # interp_percentile clamps past the grid edges to 0.5 / 99.5, so every badly deviant feature ties at one number
    # and the panel would fill by feature_id instead of by deviation. The percentile only breaks remaining ties.
    dist = _outside_distance(gs, ref)
    pools = {
        "low": _deviant(gs[gs["call"] == "low"], dist, True),
        "high": _deviant(gs[gs["call"] == "high"], dist, False),
        "missing": list(pr[gs["call"] == "expected_but_missing"].sort_values(ascending=False, kind="mergesort").index),
        "within": list(pr[gs["call"] == "within"].sort_values(ascending=False, kind="mergesort").index),
    }
    quotas = _quotas(cap)
    take = {name: pools[name][:q] for name, q in quotas}
    left = cap - sum(len(v) for v in take.values())
    for name, _ in quotas:
        if left <= 0:
            break
        extra = pools[name][len(take[name]):][:left]
        take[name] += extra
        left -= len(extra)
    return [f for name, _ in quotas for f in take[name]]


def _quantiles(ref_row: pd.Series, band, banded: bool) -> dict | None:
    """The reference percentiles of one feature, band-specific when the layer is scored within bands and the band's own
    columns are complete; None when even the range edges and median are missing (the feature is not drawable).

    All-or-nothing per row: mixing a banded median with a pooled range edge would draw a bar that no reference sample
    ever had. The 1st/99th whisker and the middle-half bar are optional and simply not drawn when the bundle's
    percentile grid lacks them."""
    for suffix in ((f"_band_{band}", "") if banded and band else ("",)):
        q = {}
        for key, col in PCOL.items():
            c = col + suffix
            v = ref_row.get(c, np.nan)
            q[key] = float(v) if pd.notna(v) else None
        if all(q[k] is not None for k in ESSENTIAL):
            return q
    return None


def _dev(mode: str, v, med: float):
    """One value as a distance from the reference median, in the layer's plotting unit."""
    if v is None or not np.isfinite(v):
        return float("nan")
    if mode == "log2_diff":
        return (v - med) / math.log(2)
    if mode == "log2_ratio":
        return math.log2(v / med) if v > 0 and med > 0 else float("nan")
    return v - med


def _panel_rows(g: pd.DataFrame, ref: pd.DataFrame, bundle: Bundle, layer: str, band, cap: int, names: dict) -> list[dict]:
    """The drawable rows of one panel, ordered the way a report reads: absent first, then by deviation, lowest to highest."""
    mode = MODE_OF_VALUE.get(LAYER_SPEC[layer][0])
    if mode is None:
        return []
    banded = layer in bundle.band_layers
    prev = _prevalence(ref, band)
    rows = []
    for fid in _select(g, ref, prev, cap):
        if fid not in ref.index:
            continue
        q = _quantiles(ref.loc[fid], band, banded)
        if q is None or (mode == "log2_ratio" and q["med"] <= 0):
            continue
        r = g.loc[fid]
        d = {k: _dev(mode, q[k], q["med"]) if q[k] is not None else None for k in PCOL}
        label = names.get(str(fid), str(fid))
        rows.append(dict(fid=str(fid), label=label, call=str(r["call"]), pct=float(r["percentile"]) if pd.notna(r["percentile"]) else None,
                         x=_dev(mode, r["value"] if pd.notna(r["value"]) else None, q["med"]), prev=float(prev.get(fid, np.nan)), **d))
    rows.sort(key=lambda r: (0 if r["call"] == "expected_but_missing" else 1,
                             -r["prev"] if r["call"] == "expected_but_missing" else (r["x"] if np.isfinite(r["x"]) else 0.0), r["fid"]))
    return rows


# ------------------------------------------------------------------------------------------------- drawing

def _colour(call: str, pct):
    """The paper's deviation palette: blue below the reference range, light blue expected but missing, orange above, grey
    within. The deep step marks a value beyond the 1st/99th percentile, as in the paper's figures."""
    if call == "expected_but_missing":
        return DEV["low_light"]
    if call == "low":
        return DEV["low_deep"] if pct is not None and pct < 1 else DEV["low_light"]
    if call == "high":
        return DEV["high_deep"] if pct is not None and pct > 99 else DEV["high_light"]
    return DEV["neutral"]


def _fold_text(mode: str, x: float) -> str:
    if not np.isfinite(x):
        return ""
    if mode == "linear":
        return "at the reference median" if abs(x) < 0.005 else f"{x:+.2f} vs the median"
    if abs(x) < 0.05:   # under 3.5%: "1.0x higher" would read as a finding where there is none
        return "at the reference median"
    f = 2.0 ** x
    if f >= 1:
        return (f"{f:.1f}x higher" if f < 10 else f"{f:,.0f}x higher")
    g = 1 / f
    return (f"{g:.1f}x lower" if g < 10 else f"{g:,.0f}x lower")


def _ordinal(p) -> str:
    if p is None or not np.isfinite(p):
        return ""
    if p < 1:
        return "below the 1st percentile"
    if p > 99:
        return "above the 99th percentile"
    n = int(round(p)); suf = "th" if 10 <= n % 100 <= 20 else {1: "st", 2: "nd", 3: "rd"}.get(n % 10, "th")
    return f"{n}{suf} percentile"


def _tick_label(mode: str, t: float) -> str:
    if mode == "linear":
        return f"{t:+g}" if t else "typical"
    if t == 0:
        return "typical"
    f = 2.0 ** t
    return f"{f:g}x" if f >= 1 else f"1/{1 / f:g}"


def _limits(rows: list[dict]) -> tuple[float, float]:
    """Axis limits from the reference bars, widened to hold the sample's dots but never by more than the reference span:
    a sample a thousand-fold off would otherwise shrink every reference bar to a line, so such dots are clipped to the
    edge and drawn as arrows instead."""
    los = [r["lo99"] if r["lo99"] is not None else r["lo"] for r in rows]
    his = [r["hi99"] if r["hi99"] is not None else r["hi"] for r in rows]
    lo, hi = min(los), max(his)
    span = max(hi - lo, 1e-6)
    xs = [r["x"] for r in rows if np.isfinite(r["x"])]
    if xs:
        lo = max(lo - span, min(lo, min(xs)))
        hi = min(hi + span, max(hi, max(xs)))
    pad = 0.08 * max(hi - lo, 1e-6)
    return lo - pad, hi + pad


def _draw_panel(plt, ax, rows: list[dict], mode: str, title: str, italic: bool, title_x: float = 0.0):
    from matplotlib.ticker import MaxNLocator
    xl = _limits(rows)
    inset = 0.02 * (xl[1] - xl[0])   # keeps an off-scale marker inside the axes rather than on the spine
    ys = []
    for y, r in enumerate(rows):
        lo99 = max(r["lo99"] if r["lo99"] is not None else r["lo"], xl[0])
        hi99 = min(r["hi99"] if r["hi99"] is not None else r["hi"], xl[1])
        ax.plot([lo99, hi99], [y, y], color=INK["axis"], lw=0.6, zorder=1, solid_capstyle="butt")
        x0, x1 = max(r["lo"], xl[0]), min(r["hi"], xl[1])
        ax.add_patch(plt.Rectangle((x0, y - 0.3), x1 - x0, 0.6, facecolor=INK["grid"], edgecolor="none", zorder=2))
        if r["q1"] is not None and r["q3"] is not None:
            ax.add_patch(plt.Rectangle((r["q1"], y - 0.3), r["q3"] - r["q1"], 0.6, facecolor=INK["axis"], edgecolor="none", zorder=3))
        ax.plot([0, 0], [y - 0.34, y + 0.34], color=INK["secondary"], lw=0.8, zorder=4)
        if r["call"] == "expected_but_missing" or not np.isfinite(r["x"]):
            ax.scatter([xl[0] + inset], [y], s=22, facecolor="none", edgecolor=_colour("expected_but_missing", None), linewidth=0.9, zorder=6)
            txt = f"not detected; carried by {r['prev'] * 100:.0f}% of the reference" if np.isfinite(r["prev"]) else "not detected"
        else:
            x = float(np.clip(r["x"], xl[0] + inset, xl[1] - inset))
            marker = "<" if r["x"] < xl[0] + inset else (">" if r["x"] > xl[1] - inset else "o")
            ax.scatter([x], [y], s=30, marker=marker, color=_colour(r["call"], r["pct"]), edgecolor=INK["surface"], linewidth=0.5, zorder=6)
            txt = " / ".join(t for t in (_fold_text(mode, r["x"]), _ordinal(r["pct"])) if t)
        ax.text(1.02, y, txt, va="center", ha="left", fontsize=PT_MIN, color=INK["secondary"],
                transform=ax.get_yaxis_transform(), clip_on=False)
        ys.append(y)
    ax.set_yticks(ys)
    ax.set_yticklabels([r["label"] for r in rows], color=INK["primary"], fontstyle="italic" if italic else "normal")
    ax.tick_params(axis="y", length=0, labelsize=PT_MIN)
    ax.set_ylim(len(rows) - 0.4, -0.7)
    ax.set_xlim(*xl)
    if mode != "linear":   # a tick at 0.5 log2 units would read "1.41x"; whole doublings are what the axis is for
        ax.xaxis.set_major_locator(MaxNLocator(nbins=8, integer=True))
    ticks = [t for t in ax.get_xticks() if xl[0] <= t <= xl[1]]
    ax.set_xticks(ticks)
    ax.set_xticklabels([_tick_label(mode, t) for t in ticks])
    ax.grid(axis="x", color=INK["grid"], lw=0.4, zorder=0)
    ax.set_axisbelow(True)
    for sp in ("left", "top", "right"):
        ax.spines[sp].set_visible(False)
    ax.set_title(title, loc="left", x=title_x, fontsize=PT_TITLE, pad=4)


def _int(v) -> int:
    """A summaries.tsv count as an int; a missing or non-finite entry counts as zero (an unscored layer writes NaN)."""
    try:
        f = float(v)
    except (TypeError, ValueError):
        return 0
    return int(f) if np.isfinite(f) else 0


def _panel_title(layer: str, srow) -> str:
    """"<layer>: n of m assessed features outside the reference range", the counts taken from summaries.tsv."""
    head = LAYER_LABEL.get(layer, layer.replace("_", " "))
    if srow is None:
        return head
    n_out = _int(srow.get("n_low_raw")) + _int(srow.get("n_high_raw"))
    n_ass = _int(srow.get("n_assessed"))
    tail = f"{n_out:,} of {n_ass:,} assessed features outside the reference range" if n_ass else "no assessable features"
    miss = _int(srow.get("n_missing"))
    return f"{head}: {tail}" + (f"; {miss:,} expected but missing" if miss else "")


def _figure(plt, sample: str, panels: list[tuple], bundle: Bundle, band) -> "object":
    """One figure, one panel per layer, laid out in inches rather than by gridspec fractions.

    Every row must get the same height whatever a panel holds, and the furniture around a panel (its heading, the axis
    and its caption) is a fixed number of inches, not a fraction of a figure whose height varies with the row count;
    placing the axes directly is the only way to hold both, and it keeps the margins from swallowing the captions."""
    body = [(len(rows) + 1.1) * ROW_IN for _, rows, _ in panels]   # +1.1: the half row of padding above and below (set_ylim)
    h = HEAD_IN + sum(b + TITLE_IN + FOOT_IN for b in body)
    fig = plt.figure(figsize=(DOUBLE_IN, h))
    y = h - HEAD_IN
    for i, ((layer, rows, title), bh) in enumerate(zip(panels, body)):
        bottom = y - TITLE_IN - bh
        ax = fig.add_axes([LEFT, bottom / h, WIDTH, bh / h])
        mode = MODE_OF_VALUE[LAYER_SPEC[layer][0]]
        # the heading starts at the figure's left edge, not the axes', or a count of features runs off the page
        _draw_panel(plt, ax, rows, mode, title, italic=layer.startswith("taxonomy_"), title_x=-(LEFT - MARGIN) / WIDTH)
        # each panel names its own unit (the layers do not share one); the last one also carries the key to the bars
        ax.set_xlabel(X_LABEL[mode] + ("; grey bar = reference range (2.5th-97.5th percentile), darker = middle half"
                                       if i == len(panels) - 1 else ""), fontsize=PT_MIN)
        y = bottom - FOOT_IN
    # two heading lines rather than one: a bundle id plus the pool size runs off a 7-inch page on one line
    man = bundle.manifest
    sub = " · ".join([str(man.get("bundle_id", os.path.basename(bundle.path))),
                      f"{man.get('n_samples', 0):,} reference samples, {man.get('n_studies', 0)} studies"]
                     + ([f"quality band {band}"] if band else []))
    fig.text(MARGIN, 1 - 0.20 / h, sample, ha="left", va="top", fontsize=PT_TITLE, color=INK["primary"])
    fig.text(MARGIN, 1 - 0.36 / h, sub, ha="left", va="top", fontsize=PT_MIN, color=INK["secondary"])
    return fig


# ------------------------------------------------------------------------------------------------- the entry point

def write_figures(assessed_dir, out_dir, samples=None, fmt="png") -> list:
    """Draw the laboratory-report figure for each assessed sample and return the paths written, in the order written.

    `assessed_dir` is an assess output directory (scores.parquet, summaries.tsv); the reference percentiles come from the
    bundle whose id that directory records, resolved like any --bundle argument (a path, or an id under
    CHECKGM_BUNDLES). `samples` restricts and orders the output; the default is every sample with scored rows, sorted.
    `fmt` is one of png, pdf, svg.

    Which features a figure shows. A sample carries hundreds of assessable features and a legible figure holds about
    twenty, so the rows are picked by a fixed rule rather than by hand: two panels, one taxonomic layer and one
    functional layer (the first of each that the sample has rows on), and in each panel ten rows filled in this order --
    the three features with the lowest percentile among those called low, the three with the highest percentile among
    those called high, the two most prevalent reference features that are expected but missing, and the two most
    prevalent features whose value falls within the range, as context for what a normal row looks like. Quota a pool
    cannot fill passes to the pools after it in that same order, so a sample with nothing high shows more low features
    rather than more context. The selection is therefore symmetric in direction, always shows some unremarkable
    features, and is reproducible: a reader can check it, and no choice of rows can flatter or dramatize a sample.

    A sample that cannot be drawn -- it was rejected before scoring, or no layer of it has a drawable feature -- is not
    an error: it is simply absent from the returned list, and the caller reports it.
    """
    fmt = str(fmt).lower().lstrip(".")
    if fmt not in FORMATS:
        raise ValueError(f"unsupported figure format '{fmt}': use one of {', '.join(FORMATS)}")
    scores, summ = _read_assessed(assessed_dir)
    if scores.empty:
        return []
    scores = scores.astype({"analysis_id": str, "feature_id": str})
    bundle = _bundle_of(assessed_dir, summ)
    have = list(dict.fromkeys(scores["analysis_id"]))
    if samples is None:
        want = sorted(have)
    else:
        want = [str(s) for s in samples]
        unknown = [s for s in want if s not in set(have)]
        if unknown:
            shown = ", ".join(sorted(have)[:8]) + (f", ... ({len(have)} in all)" if len(have) > 8 else "")
            raise ValueError(f"sample(s) not in {assessed_dir}: {', '.join(unknown)}; the directory holds {shown}")
    os.makedirs(out_dir, exist_ok=True)
    band_of, srow_of = {}, {}
    if not summ.empty:
        S = summ.copy(); S["analysis_id"] = S["analysis_id"].astype(str)
        if "quality_band" in S.columns:
            band_of = {s: g["quality_band"].dropna().iloc[0] for s, g in S.groupby("analysis_id") if g["quality_band"].notna().any()}
        srow_of = {(s, l): r for s, l, r in zip(S["analysis_id"], S["layer"], S.to_dict("records"))}

    written = []
    with _quiet():
        plt, _style = _pyplot()
        for sid in want:
            sc = scores[scores["analysis_id"] == sid]
            band = band_of.get(sid) or (sc["band"].dropna().iloc[0] if "band" in sc.columns and sc["band"].notna().any() else None)
            # restricted to the bundle's own layers: a directory assessed against another build of the same id would
            # otherwise reach for a feature table that is not there, mid-drawing
            scorable = {l for l, g in sc.groupby("layer") if l in LAYER_SPEC and (g["call"] != "not_assessable").any()} & set(bundle.layers)
            chosen = [next((l for l in pref if l in scorable), None) for pref in (TAXON_LAYERS, FUNCTION_LAYERS)]
            chosen = [l for l in chosen if l]
            if not chosen:
                continue
            cap = max(4, ROW_BUDGET // len(chosen))
            panels = []
            for layer in chosen:
                g = sc[(sc["layer"] == layer) & (sc["call"] != "not_assessable")].set_index("feature_id")
                ref = bundle.features(layer)
                if ref.index.dtype != g.index.dtype:   # never rewrite the Bundle's cached table in place
                    ref = ref.rename(index=str)
                rows = _panel_rows(g, ref, bundle, layer, band, cap, _taxon_names(bundle, g.index) if layer.startswith("taxonomy_") else {})
                if rows:
                    panels.append((layer, rows, _panel_title(layer, srow_of.get((sid, layer)))))
            if not panels:
                continue
            with plt.rc_context(_style):   # scoped: rcParams is global, and checkgm may be embedded in a process that plots
                fig = _figure(plt, sid, panels, bundle, band)
                path = os.path.join(out_dir, f"{re.sub(r'[^A-Za-z0-9._-]', '_', sid)}_range_report.{fmt}")
                fig.savefig(path, format=fmt, bbox_inches=None)
                plt.close(fig)
            written.append(path)
    return written
