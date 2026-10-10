#!/usr/bin/env Rscript
# Figures 1-5 of the two-column paper, drawn with ggplot2 + patchwork from the tidy tables that scripts/viz2/prep.py writes to
# figures/v2/data. One theme and one palette throughout:
#   blue = below the healthy range (light blue = expected but missing), orange = above it, grey = within;
#   violet = checkGM, greys = comparators; teal / ochre = assembly / read pipeline (schematics and curation only).
# Usage: Rscript scripts/viz2/figures.R [fig1 fig2 ...]   (default: all). Output: figures/v2/figN.pdf (+ .png preview, 300 dpi).
suppressPackageStartupMessages({library(ggplot2); library(patchwork); library(dplyr); library(tidyr); library(readr); library(scales)})
options(readr.show_col_types = FALSE, dplyr.summarise.inform = FALSE)
# the project root comes from the environment, as the README documents, and falls back to the tree this was written in
P <- Filter(nzchar, c(Sys.getenv("CHECKGM_PROJECT"), Sys.getenv("REFMB_PROJECT"), "/home/classes/bios/270/khoa/meta2/project"))[1]
D <- file.path(P, "figures/v2/data"); OUT <- file.path(P, "figures/v2")
FONT <- "Nimbus Sans"
SINGLE <- 3.5; DOUBLE <- 7.16                                       # IEEE column and page widths (in)

INK <- "#1A1A1A"; INK2 <- "#4D4D4D"; MUTED <- "#7F7F7F"; AXIS <- "#A6A6A6"; GRID <- "#ECECEC"
BAND <- "#E3DFD5"; IQR <- "#C6C0B0"; MED <- "#5C574C"
LOW <- "#2E6DB4"; LOW_L <- "#A8C3E5"; HIGH <- "#D1612B"; WITHIN <- "#9E9E9E"
ACC <- "#5E3C99"; GREY_D <- "#4A4A4A"; GREY_M <- "#8C8C8C"; GREY_L <- "#BDBDBD"
PA <- "#2C7C84"; PA_BG <- "#E3EFF0"; PB <- "#A8761E"; PB_BG <- "#F5ECDA"; NEUT_BG <- "#F3F1EC"
pt <- function(x) x / .pt                                            # font size in pt -> ggplot text size (mm)
MINUS <- "−"

theme_checkgm <- function(base = 7) {
  theme_minimal(base_size = base, base_family = FONT) %+replace% theme(
    text = element_text(colour = INK, family = FONT, size = base),
    axis.text = element_text(size = base - 0.5, colour = INK2), axis.title = element_text(size = base, colour = INK),
    axis.title.x = element_text(margin = margin(t = 3)), axis.title.y = element_text(margin = margin(r = 3), angle = 90),
    axis.line = element_line(colour = AXIS, linewidth = 0.3), axis.ticks = element_line(colour = AXIS, linewidth = 0.3),
    axis.ticks.length = unit(1.6, "pt"),
    panel.grid.major = element_line(colour = GRID, linewidth = 0.25), panel.grid.minor = element_blank(),
    strip.text = element_text(size = base, face = "bold", hjust = 0, colour = INK, margin = margin(2, 0, 2, 0)),
    legend.text = element_text(size = base - 0.5), legend.title = element_text(size = base - 0.5, colour = INK2),
    legend.key.size = unit(7, "pt"), legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(3, "pt"),
    plot.title = element_text(size = base, face = "bold", hjust = 0, margin = margin(0, 0, 3, 0)),
    plot.subtitle = element_text(size = base - 0.5, hjust = 0, colour = INK2, margin = margin(0, 0, 3, 0)),
    plot.tag = element_text(size = base + 2, face = "bold"),
    plot.margin = margin(2, 4, 2, 2))
}
theme_set(theme_checkgm())
tags <- function() plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(size = 9, face = "bold", family = FONT, hjust = 0, vjust = 1))

save_fig <- function(p, name, w, h) {
  ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h, units = "in", device = cairo_pdf)
  ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, units = "in", dpi = 300, device = ragg::agg_png, bg = "white")
  message("wrote ", name)
}
jnum <- function(path) {                                             # flat JSON of numbers -> named list (no jsonlite needed)
  j <- paste(readLines(path, warn = FALSE), collapse = ""); kv <- regmatches(j, gregexpr('"[^"]+": *[-0-9.eE]+', j))[[1]]
  setNames(as.list(as.numeric(sub('.*: *', '', kv))), gsub('"|:.*', '', kv))
}
fold_lab <- function(x) {                                            # log2 difference -> "139x lower"
  f <- 2^abs(x); s <- ifelse(f >= 10, comma(round(f)), sprintf("%.1f", f)); paste0(s, "× ", ifelse(x < 0, "lower", "higher"))
}
pct_lab <- function(p) ifelse(p <= 1, "below the 1st percentile", ifelse(p >= 99, "above the 99th percentile", paste(ordinal(round(p)), "percentile")))
log2_axis <- function(br) ifelse(br == 0, "typical", ifelse(br > 0, paste0(2^br, "×"), paste0("1/", 2^-br)))
# discrete-axis labeller that italicises the taxon names in `ital`, so families and genera are set the same way in the
# figures as in the text; everything else (gene and pathway names) stays roman
ital_lab <- function(ital) function(br) as.expression(lapply(br, function(n) if (n %in% ital) bquote(italic(.(n))) else bquote(.(n))))

# ---- schematic helpers. The schematics (figs 1a, 2a) are laid out in a unit that renders as the same length in x and in
# y, which coord_cartesian does not give: pinning the patchwork row to H inches leaves the row H + 0.152 in tall, so the
# y unit came out 1.17x the x unit and neither the round corners nor the text block inside a box was what it claimed to
# be. coord_fixed is used instead, and the row is pinned to the inches the canvas then needs. Two numbers are measured in
# the rendered PDF rather than derived: patchwork insets the assembled figure by MARG on every side, which leaves PTU
# points for each unit of a DOUBLE-wide canvas, and gives the pinned row that inset back at top and bottom. Only the
# point-to-unit conversions below depend on them, so a few per cent would move a text block by a quarter of a point, and
# pinning slightly tall merely centres the canvas in its row whereas pinning short would also narrow it.
MARG <- 5.48; PTU <- 70.40
U <- function(p) p / PTU                                             # real points -> canvas units
pin_row <- function(H) unit(c(H * PTU / 72 - 2 * MARG / 72 + 0.007, 1), c("in", "null"))
rr <- function(x0, x1, y0, y1, r = 0.035, n = 7) {
  a <- seq(0, pi / 2, length.out = n)
  data.frame(x = c(x1 - r + r * cos(a), x0 + r - r * sin(a), x0 + r - r * cos(a), x1 - r + r * sin(a)),
             y = c(y1 - r + r * sin(a), y1 - r + r * cos(a), y0 + r - r * sin(a), y0 + r - r * cos(a)))
}
box <- function(x0, x1, y0, y1, title, sub = NULL, fill = "white", col = INK2, lty = "solid", ts = 6.8, ss = 6.2, hj = 0.5) {
  xc <- if (hj == 0) x0 + 0.06 else (x0 + x1) / 2; yc <- (y0 + y1) / 2
  L <- list(geom_polygon(data = rr(x0, x1, y0, y1), aes(x, y), fill = fill, colour = col, linewidth = 0.35, linetype = lty, inherit.aes = FALSE))
  ns <- if (is.null(sub)) 0 else length(strsplit(sub, "\n")[[1]])
  # the block is measured in points, where the type is: a line box is one font size tall and grid sets successive lines
  # 1.2 sizes apart, so the title's box plus its leading is 1.35 sizes. Centring the block on the box then needs only the
  # point-to-unit conversion, and both gaps come out equal instead of 2.4 pt above the text and 5.0 pt below it.
  th <- ts * 1.35; blk <- if (ns > 0) th + ss + (ns - 1) * ss * 1.2 else ts
  L <- c(L, annotate("text", xc, yc + U(blk) / 2, label = title, size = pt(ts), fontface = "bold", family = FONT, colour = INK, vjust = 1, hjust = hj, lineheight = 0.95))
  if (ns > 0) L <- c(L, annotate("text", xc, yc + U(blk) / 2 - U(th), label = sub, size = pt(ss), family = FONT, colour = INK2, vjust = 1, hjust = hj, lineheight = 1.0))
  L
}
arr <- function(x0, y0, x1, y1, col = INK2) annotate("segment", x = x0, y = y0, xend = x1, yend = y1, colour = col, linewidth = 0.4,
                                                      arrow = arrow(length = unit(3.2, "pt"), type = "closed", angle = 25))
seg <- function(x0, y0, x1, y1, col = INK2) annotate("segment", x = x0, y = y0, xend = x1, yend = y1, colour = col, linewidth = 0.4)
canvas <- function(w, h) ggplot() + coord_fixed(ratio = 1, xlim = c(0, w), ylim = c(0, h), expand = FALSE, clip = "off") + theme_void(base_family = FONT) +
  theme(plot.margin = margin(0, 0, 0, 0), plot.tag = element_text(size = 9, face = "bold"), plot.tag.position = c(0.004, 0.985))
range_row <- function(x0, x1, y, lo, q1, q3, hi, med, s, col) list(   # mini healthy-range glyph for the schematics
  annotate("rect", xmin = lo, xmax = hi, ymin = y - 0.035, ymax = y + 0.035, fill = BAND),
  annotate("rect", xmin = q1, xmax = q3, ymin = y - 0.035, ymax = y + 0.035, fill = IQR),
  annotate("segment", x = med, xend = med, y = y - 0.05, yend = y + 0.05, colour = MED, linewidth = 0.45),
  annotate("point", x = s, y = y, colour = "white", fill = col, shape = 21, size = 1.9, stroke = 0.3))

# ======================================================================================== Figure 1: overview + one report
fig1 <- function() {
  W <- DOUBLE; H <- 1.02
  reads <- data.frame(x = 0.70 + c(0, .05, .02, .07, .01, .05, .03), y = seq(0.36, 0.72, length.out = 7), len = c(.25, .2, .27, .18, .24, .21, .26))
  a <- canvas(W, H) +
    geom_polygon(data = rr(0.26, 0.48, 0.32, 0.70, r = 0.08), aes(x, y), fill = "#CDB49A", colour = INK2, linewidth = 0.35) +
    annotate("rect", xmin = 0.23, xmax = 0.51, ymin = 0.68, ymax = 0.77, fill = INK2) +
    geom_segment(data = reads, aes(x = x, xend = x + len, y = y, yend = y), colour = c(MUTED, GREY_L)[rep(1:2, 4)][1:7], linewidth = 0.9, lineend = "round") +
    annotate("text", 0.6, 0.17, label = "gut metagenome\n(raw reads)", size = pt(6.5), family = FONT, colour = INK, lineheight = 0.95) +
    seg(1.03, 0.52, 1.12, 0.52) + seg(1.12, 0.275, 1.12, 0.755) + arr(1.12, 0.755, 1.24, 0.755) + arr(1.12, 0.275, 1.24, 0.275) +
    box(1.25, 2.95, 0.56, 0.95, "assembly pipeline", "MGnify v5 workflow: assemble,\ncall genes, annotate", PA_BG, PA) +
    box(1.25, 2.95, 0.08, 0.47, "read pipeline", "curatedMetagenomicData 3 workflow:\nMetaPhlAn 3 and HUMAnN 3", PB_BG, PB) +
    arr(2.97, 0.755, 3.19, 0.755) + arr(2.97, 0.275, 3.19, 0.275) +
    box(3.2, 4.88, 0.56, 0.95, "scored against the reference", "1,941 adults, 26 studies;\nstatistics only, no samples", NEUT_BG, PA) +
    box(3.2, 4.88, 0.08, 0.47, "scored against the reference", "6,494 adults, 19 studies;\nstatistics only, no samples", NEUT_BG, PB) +
    seg(4.9, 0.755, 4.99, 0.755) + seg(4.9, 0.275, 4.99, 0.275) + seg(4.99, 0.275, 4.99, 0.755) + arr(4.99, 0.515, 5.11, 0.515) +
    # The schematic used to end in a thumbnail of the report, which panel (c) draws properly and Fig. 3 draws for a
    # whole cohort: three renderings of one idea inside one figure. The last step is the pair of outputs instead, which
    # is also what the flow was missing, the deviation score having appeared nowhere in it.
    box(5.12, 7.13, 0.56, 0.95, "a report for every feature",
        "percentile, and a call of within,\nbelow or above the range (c)", NEUT_BG, INK2) +
    box(5.12, 7.13, 0.08, 0.47, "a deviation score",
        "one number per sample,\nweighted over its percentiles", NEUT_BG, ACC)

  # b: how one feature is placed
  m <- jnum(file.path(D, "f1_density_meta.json")); dz <- read_csv(file.path(D, "f1_density.csv"))
  dd <- density(dz$x, adjust = 1.1, n = 1024); dd <- data.frame(x = dd$x, y = dd$y / max(dd$y)); dd <- dd[dd$x > -10.5 & dd$x < 5.5, ]
  inb <- dd[dd$x >= m[["2_5"]] & dd$x <= m[["97_5"]], ]
  b <- ggplot(dd, aes(x, y)) + geom_area(data = inb, fill = BAND) + geom_line(colour = INK2, linewidth = 0.4) +
    # the percentile ticks sit a little below the top of the plotting window, which keeps them clear of the second line of
    # the annotation that labels the sample; the two shared a column and their glyphs interleaved
    annotate("segment", x = c(m[["2_5"]], m[["97_5"]]), xend = c(m[["2_5"]], m[["97_5"]]), y = 0, yend = 1.04, colour = MUTED, linewidth = 0.3, linetype = "22") +
    annotate("text", x = c(m[["2_5"]], m[["97_5"]]), y = 1.07, label = c("2.5th", "97.5th"), size = pt(6), family = FONT, colour = INK2, vjust = 0) +
    annotate("text", x = -0.6, y = 0.3, label = "reference\nrange", size = pt(6.2), family = FONT, colour = INK2, lineheight = 0.95) +
    annotate("segment", x = m$sample, xend = m$sample, y = 0, yend = 1.33, colour = LOW, linewidth = 0.6) +
    annotate("point", x = m$sample, y = 0, colour = "white", fill = LOW, shape = 21, size = 2.2, stroke = 0.3) +
    annotate("text", x = m$sample + 0.3, y = 1.6, label = paste0("this sample: ", fold_lab(m$sample), "\n(", pct_lab(m$pct), ")"), hjust = 0, vjust = 1,
             size = pt(6.2), family = FONT, colour = LOW, lineheight = 0.95) +
    scale_x_continuous(breaks = c(-8, -4, 0, 4), labels = log2_axis, limits = c(-10.5, 5.5), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1.62), expand = c(0, 0)) +
    # The axis names its unit. For a family the plotted quantity is a difference of centred log-ratios (prep.py writes
    # (clr - p50)/log 2, and refmb/normalize.py defines clr = log p - mean log p over the basis), so 2^x is the ratio of
    # CLR-normalised values, not the ratio of relative abundances: this sample's Oscillospiraceae is 364x lower on that
    # scale and 59x lower on proportion_mapped against the pool median. The band, the ticks and the percentile all live in
    # the same CLR space, which is the only space the shipped reference stores for taxa, so the unit is named rather than
    # the number changed.
    # two lines: the one-line form is 170 pt wide against the 136 pt this panel can give an axis title, and the leading
    # word was being clipped at the figure's left edge
    labs(x = "relative to the reference median\n(centred log-ratio)", y = NULL,
         title = expression(bold("placing one feature: ") * bolditalic("Oscillospiraceae"))) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(), panel.grid = element_blank(), panel.grid.major = element_blank())

  # c: the report for the same sample (bars clipped to the axis range)
  XL <- c(-14, 8.5); cl <- function(v) pmin(pmax(v, XL[1]), XL[2])
  r <- read_csv(file.path(D, "f1_report.csv")) %>% mutate(call = ifelse(is.na(x), "absent", call))
  r$label <- sub("^lysine racemase.*", "lysine racemase (K20707)", r$label); r$label <- sub("^ECF sigma.*", "ECF sigma factor (K03088)", r$label)
  r <- r %>% mutate(txt = ifelse(call == "absent", "absent; carried by most reference adults", paste0(fold_lab(x), " · ", pct_lab(pct))),
                    layer = factor(ifelse(layer == "Families", "taxa", "genes"), c("taxa", "genes")),
                    across(c(lo99, lo, q1, q3, hi, hi99), cl))
  fam <- grepl("^taxa", r$layer) & r$call != "absent"; gen <- grepl("^genes", r$layer)
  ordr <- c(r$label[fam][order(r$x[fam])], r$label[r$call == "absent"], r$label[gen][order(r$x[gen])])
  r$label <- factor(r$label, rev(ordr))
  c_ <- ggplot(r, aes(y = label)) +
    geom_rect(data = filter(r, grepl("^genes", layer)) %>% distinct(layer), inherit.aes = FALSE,
              aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), fill = "#F6F4EF") +
    geom_segment(aes(x = lo99, xend = hi99, yend = label), colour = AXIS, linewidth = 0.3) +
    geom_tile(aes(x = (lo + hi) / 2, width = hi - lo), height = 0.62, fill = BAND) +
    geom_tile(aes(x = (q1 + q3) / 2, width = q3 - q1), height = 0.62, fill = IQR) +
    geom_tile(aes(x = 0), width = 0.14, height = 0.78, fill = MED) +
    geom_point(data = filter(r, call != "absent"), aes(x = x, fill = call), shape = 21, colour = "white", size = 2.1, stroke = 0.3) +
    geom_point(data = filter(r, call == "absent"), aes(x = XL[1] + 0.35), shape = 21, colour = LOW_L, fill = "white", size = 2.1, stroke = 0.9) +
    geom_text(aes(x = XL[2] + 0.4, label = txt), hjust = 0, size = pt(6.2), family = FONT, colour = INK2) +
    facet_grid(layer ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_fill_manual(values = c(low = LOW, high = HIGH, within = WITHIN), guide = "none") +
    scale_y_discrete(labels = ital_lab(as.character(r$label[r$layer == "taxa"]))) +
    scale_x_continuous(breaks = c(-12, -8, -4, 0, 4, 8), labels = log2_axis, expand = c(0, 0)) +
    coord_cartesian(xlim = XL, clip = "off") +
    # The two facets share one axis but not one unit: the taxa rows are differences of centred log-ratios (fig17_range_report
    # rows(), x = (clr - p50)/log 2), the gene rows are plain log2 ratios of copies per genome (x = log2(value/median)).
    # Only the taxa half needs naming, so the parenthesis qualifies it alone.
    labs(x = "relative to the reference median; taxa in centred log-ratio", y = NULL,
         title = expression(bold("report for one ") * bolditalic("C. difficile") * bold(" infection patient (excerpt)"))) +
    theme(panel.grid.major.y = element_blank(), axis.text.y = element_text(colour = INK, size = 6.5),
          # the two-row gene facet gives its rotated label a band shorter than the word, and "genes" was being clipped to
          # "jene"; strip.clip = "off" lets the point or so of overflow fall in the panel spacing
          strip.placement = "outside", strip.clip = "off",
          strip.text.y.left = element_text(angle = 90, hjust = 0.5, size = 6.2, face = "bold", colour = INK2),
          panel.spacing.y = unit(5, "pt"),
          plot.margin = margin(2, 122, 2, 2))
  # the schematic row is pinned to the inches its canvas needs rather than given a share of the figure: a relative height
  # rescales the coordinate system but not the pt-sized text, and the box borders then cut through their own last line
  p <- wrap_elements(full = a) / ((b | c_) + plot_layout(widths = c(1, 1.32))) + plot_layout(heights = pin_row(H)) +
    plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(size = 9, face = "bold", family = FONT))
  save_fig(p, "fig1", DOUBLE, 2.84)
}

# ======================================================================================== Figure 2: pipelines, curation, checks
LAYERS <- c("Family", "Genus", "Species", "KO (eggNOG)", "KO (KOfam)", "Pfam", "KEGG module", "KO", "Pathway")
fig2 <- function() {
  W <- DOUBLE; H <- 1.06
  a <- canvas(W, H) +
    annotate("text", 0.16, 0.86, label = "assembly\n(MGnify v5)", hjust = 0, vjust = 1, size = pt(6.8), fontface = "bold", colour = PA, family = FONT, lineheight = 0.95) +
    annotate("text", 0.16, 0.36, label = "read\n(cMD3)", hjust = 0, vjust = 1, size = pt(6.8), fontface = "bold", colour = PB, family = FONT, lineheight = 0.95) +
    arr(0.7, 0.77, 0.81, 0.77, PA) + arr(0.7, 0.27, 0.81, 0.27, PB) +
    box(0.82, 2.14, 0.58, 0.96, "assemble", "metaSPAdes 3.15.3 (paired)\nMEGAHIT 1.2.9 (single-end)", PA_BG, PA) +
    box(2.27, 3.31, 0.58, 0.96, "call genes", "Prodigal 2.6.3,\nFragGeneScan", PA_BG, PA) +
    box(3.44, 5.27, 0.58, 0.96, "annotate", "DIAMOND vs UniRef90 (taxa);\neggNOG-mapper, KOfam, InterProScan", PA_BG, PA) +
    box(0.82, 2.14, 0.08, 0.46, "profile taxa", "MetaPhlAn 3.0.14\n(markers mpa_v30)", PB_BG, PB) +
    box(2.27, 3.31, 0.08, 0.46, "profile genes", "HUMAnN 3\n(KOs, MetaCyc pathways)", PB_BG, PB) +
    arr(2.15, 0.77, 2.26, 0.77, PA) + arr(3.32, 0.77, 3.43, 0.77, PA) + arr(2.15, 0.27, 2.26, 0.27, PB) +
    arr(5.28, 0.77, 5.39, 0.77, PA) + arr(3.32, 0.27, 5.39, 0.27, PB) +
    box(5.4, 6.38, 0.58, 0.96, "features", "families, genera,\nKOs, Pfam, modules", PA_BG, PA) +
    box(5.4, 6.38, 0.08, 0.46, "features", "families, genera,\nspecies, KOs, pathways", PB_BG, PB) +
    seg(6.39, 0.77, 6.45, 0.77) + seg(6.39, 0.27, 6.45, 0.27) + seg(6.45, 0.27, 6.45, 0.77) + arr(6.45, 0.52, 6.51, 0.52) +
    box(6.52, 7.14, 0.08, 0.96, "score", "normalize;\nplace each\nfeature in\nthe range", NEUT_BG, INK2, ts = 6.6, ss = 6.0)

  # b: curation, from public resource to baseline and evaluation sets
  fu <- read_csv(file.path(D, "f2_funnel.csv"))
  ev <- read_tsv(file.path(P, "results/s8/loso_auroc.tsv")) %>% filter(study != "POOLED", feature_set == "reference_relative")
  # the 15 case/control studies of the read pipeline; the 14 with >= 20 cases and >= 20 controls enter the benchmark
  # (SankaranarayananK_2015, 19 cases / 18 controls, does not), but all 15 are counted here as curated samples
  cc <- read_tsv(file.path(P, "results/s14/disease_strata_B_balance.tsv"))
  stopifnot(sum(ev$n_test) == 2294, sum(cc$n_case + cc$n_control) == 1437)
  extra <- tibble(pipeline = c("Assembly (MGnify v5)", "Assembly (MGnify v5)", "Read (cMD3)", "Read (cMD3)"), step = c(6, 7, 6, 7),
                  label = rep(c("held-out healthy adults", "case/control studies"), 2), n = c(233, sum(ev$n_test), 317, sum(cc$n_case + cc$n_control)),
                  studies = c(3, nrow(ev), 3, nrow(cc)))
  fu <- bind_rows(fu, extra) %>% mutate(role = case_when(step == 5 ~ "baseline", step >= 6 ~ "set aside", TRUE ~ "screen"),
                                        label = recode(label, "pass assembly QC" = "pass assembly QC", "gut, one per sample" = "gut, one per sample",
                                                       "adults, profile resolvable" = "adults, resolvable profile", "held-out healthy adults" = "held out, healthy"),
                                        lab = ifelse(is.na(studies), comma(n), paste0(comma(n), " (", studies, ")")),
                                        key = paste(pipeline, step), pal = paste(pipeline, role),
                                        pipeline = recode(pipeline, "Assembly (MGnify v5)" = "assembly pipeline (MGnify v5)", "Read (cMD3)" = "read pipeline (cMD3)"))
  fu$key <- factor(fu$key, rev(fu$key[order(fu$pipeline, fu$step)]))
  fills <- c("Assembly (MGnify v5) screen" = PA_BG, "Assembly (MGnify v5) baseline" = PA, "Assembly (MGnify v5) set aside" = "#9CC5C9",
             "Read (cMD3) screen" = PB_BG, "Read (cMD3) baseline" = PB, "Read (cMD3) set aside" = "#DCC08C")
  # the bottom two rows of each facet are the sets held out of the baseline rather than screened out on the way to it, so
  # a rule separates them from the funnel proper: they are siblings of the baseline and one of them is the longer bar
  b <- ggplot(fu, aes(y = key, x = n, fill = pal)) + geom_col(width = 0.72, colour = NA) +
    geom_hline(yintercept = 2.5, colour = AXIS, linewidth = 0.25, linetype = "22") +
    geom_text(aes(label = lab), hjust = -0.08, size = pt(6), family = FONT, colour = INK2) +
    facet_wrap(~pipeline, ncol = 1, scales = "free_y") +
    scale_fill_manual(values = fills, guide = "none") +
    # thirteen rows of 6 pt counts need 5.7 pt of pitch before the parentheses of consecutive labels touch wherever two
    # bars end at a similar x, which they do five times over in the assembly funnel. The pitch is bought inside the panel
    # rather than from the figure's height: a tighter discrete expansion, closer facets and leaner strips.
    scale_y_discrete(labels = setNames(fu$label, fu$key), expand = expansion(add = 0.45)) +
    scale_x_continuous(labels = function(x) ifelse(x == 0, "0", paste0(x / 1000, "k")), limits = c(0, 33500), breaks = c(0, 10000, 20000), expand = c(0, 0)) +
    labs(x = "samples (studies)", y = NULL, title = "curation") +
    theme(panel.grid.major.y = element_blank(), axis.text.y = element_text(size = 6.2, colour = INK),
          strip.text = element_text(size = 6.5, face = "plain", colour = INK2, margin = margin(1, 0, 1, 0)),
          panel.spacing.y = unit(3, "pt"))

  # c: calibration, share of features outside the range in healthy samples the baseline never saw
  ca <- read_csv(file.path(D, "f2_calibration.csv")) %>% mutate(frac = 100 * frac, layer = factor(layer, rev(LAYERS))) %>% filter(!is.na(layer)) %>%
    mutate(kind = factor(kind, c("held-out cohort", "baseline study left out")),
           pipeline = recode(pipeline, "Assembly (MGnify v5)" = "assembly pipeline", "Read (cMD3)" = "read pipeline"))
  med <- ca %>% group_by(pipeline, layer) %>% summarise(m = median(frac))
  # The Hadza point used to carry a "Hadza study (Fig. 5)" label, which the caption names anyway. A row of this panel is
  # not quite 6 pt, every row is occupied from 3 to 20 % and the family row is the topmost one, so the only y at which a
  # 6 pt line clears both its own circle and the rest of the row lies outside the panel; buying the room from the scale
  # would have left a tenth of each facet blank. The lone point at 20 % against a next-highest 13 % is unmistakable.
  c_ <- ggplot(ca, aes(x = frac, y = layer)) + geom_vline(xintercept = 5, colour = MUTED, linewidth = 0.35, linetype = "22") +
    geom_point(aes(shape = kind), position = position_jitter(height = 0.18, width = 0, seed = 1), size = 0.9, colour = INK2, stroke = 0.35, alpha = 0.8) +
    geom_point(data = med, aes(x = m), shape = 124, size = 3.2, colour = MED) +
    facet_wrap(~pipeline, ncol = 1, scales = "free_y") + scale_shape_manual(values = c(16, 1), name = NULL) +
    scale_y_discrete(expand = expansion(add = c(0.45, 0.7))) +
    scale_x_continuous(limits = c(0, 25), breaks = c(0, 5, 10, 15, 20, 25), labels = function(x) paste0(x, "%"), expand = c(0, 0.3)) +
    labs(x = "features outside the reference range", y = NULL, title = "calibration (5 % expected)") +
    theme(legend.position = "bottom", legend.text = element_text(size = 6), panel.grid.major.y = element_blank(),
          axis.text.y = element_text(size = 6.2, colour = INK2), strip.text = element_text(size = 6.5, face = "plain", colour = INK2))

  # d: the pipelines reproduce the source percentiles from raw reads
  fi <- read_csv(file.path(D, "f2_fidelity.csv")) %>% mutate(layer = factor(layer, rev(LAYERS))) %>% filter(!is.na(layer)) %>%
    mutate(pipeline = recode(pipeline, "Assembly (MGnify v5)" = "assembly pipeline, 61 samples", "Read (cMD3)" = "read pipeline, 34 samples"))
  fm <- fi %>% group_by(pipeline, layer) %>% summarise(m = median(rho), n = n())
  d <- ggplot(fi, aes(x = rho, y = layer)) +
    geom_point(position = position_jitter(height = 0.18, width = 0, seed = 1), size = 0.8, colour = INK2, alpha = 0.55, stroke = 0) +
    geom_point(data = fm, aes(x = m), shape = 124, size = 3.2, colour = MED) +
    geom_text(data = fm, aes(x = 0.405, label = sprintf("%.2f", m)), hjust = 0, size = pt(6), family = FONT, colour = MED) +
    facet_wrap(~pipeline, ncol = 1, scales = "free_y") +
    scale_x_continuous(limits = c(0.4, 1.0), breaks = c(0.4, 0.6, 0.8, 1.0), expand = c(0, 0.01)) +
    labs(x = "Spearman ρ with source percentiles", y = NULL, title = "reproduced from raw reads") +
    theme(panel.grid.major.y = element_blank(), axis.text.y = element_text(size = 6.2, colour = INK2),
          strip.text = element_text(size = 6.5, face = "plain", colour = INK2))
  p <- wrap_elements(full = a) / (b | c_ | d) + plot_layout(heights = pin_row(H)) +   # see fig1: inches, not a share
    plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(size = 9, face = "bold", family = FONT))
  save_fig(p, "fig2", DOUBLE, 3.37)
}

# ======================================================================================== Figure 3: disease benchmark
fig3 <- function() {
  au <- read_csv(file.path(D, "f3_auroc.csv"))
  # prep.py now names the three assembly feature sets at source (reference_relative = taxa and genes together, plus the
  # taxonomy-only and genes-only ablations from results/s8), so no relabelling is needed here
  # row labels drop the tool name: the violet bars are checkGM output, the greys are the comparators
  au$method <- recode(au$method, "GMHI (published)" = "GMHI", "GMWI2 (published)" = "GMWI2",
                      "Alpha diversity" = "alpha diversity",
                      "checkGM deviation score" = "deviation score",
                      "raw abundances, same model" = "raw abundances, same model")
  # "14 of 15" rather than "14": fifteen case/control studies are curated in Fig. 2b, and the one with fewer than 20
  # cases and 20 controls (SankaranarayananK_2015) is held out of the benchmark, which the panel titles otherwise hide
  lv <- c("read pipeline\n14 of 15 studies", "read pipeline\n4 studies GMWI2 never saw")
  au$pipeline <- factor(recode(au$pipeline, "Read pipeline (14 studies)" = lv[1],
                               "Read pipeline, 4 studies GMWI2 never saw" = lv[2]), lv)
  OURS <- c("deviation score")
  stopifnot(all(OURS %in% au$method))   # a rename must not silently drop the violet grouping
  s <- au %>% group_by(pipeline, method) %>% summarise(m = mean(auroc), n = n()) %>% ungroup() %>%
    mutate(lab = ifelse(method == "GMWI2" & pipeline == lv[1], "GMWI2\u2020", method), key = paste(pipeline, method),
           grp = ifelse(method %in% OURS, "checkGM", ifelse(method == "GMWI2", "GMWI2", "other")))
  # facet_wrap gives every panel the same height, so the four-method assembly panel would draw its bars half again as
  # thick as the nine-method read panels for the same quantity on the same axis. Blank rows pad the short panel to nine,
  # and are carried by geom_blank alone so that nothing is drawn and no layer sees a missing value.
  s <- s %>% arrange(pipeline, m) %>% mutate(pipeline = as.character(pipeline))
  npad <- max(table(s$pipeline)) - table(s$pipeline)
  pad <- bind_rows(lapply(names(npad)[npad > 0], function(pl)
    tibble(pipeline = pl, key = paste(pl, "pad", seq_len(npad[[pl]])), lab = "")))
  s <- bind_rows(s, pad) %>% mutate(pipeline = factor(pipeline, lv)) %>% arrange(pipeline, is.na(m), m) %>% mutate(key = factor(key, key))
  sr <- filter(s, !is.na(m))
  au <- au %>% mutate(key = factor(paste(pipeline, method), levels(s$key)))
  g <- ggplot(sr, aes(y = key)) + geom_blank(data = s, aes(y = key)) + geom_vline(xintercept = 0.5, colour = MUTED, linewidth = 0.35) +
    geom_tile(aes(x = (0.5 + m) / 2, width = abs(m - 0.5), fill = grp), height = 0.72) +
    geom_point(data = au, aes(x = auroc), position = position_jitter(height = 0.17, width = 0, seed = 2), size = 0.6, colour = INK, alpha = 0.5, stroke = 0) +
    geom_text(aes(x = 1.08, label = sprintf("%.2f", m), fontface = ifelse(grp == "checkGM", "bold", "plain")), hjust = 1, size = pt(6.4), family = FONT, colour = INK) +
    facet_wrap(~pipeline, nrow = 1, scales = "free_y") +
    scale_fill_manual(values = c(checkGM = ACC, GMWI2 = GREY_D, other = GREY_L), name = NULL,
                      labels = c(checkGM = "checkGM output", GMWI2 = "GMWI2", other = "other comparators"), breaks = c("checkGM", "GMWI2", "other")) +
    scale_y_discrete(labels = setNames(s$lab, s$key)) +
    # the lower limit reaches the smallest study AUROC (alpha diversity on GuptaA_2019, 0.17): clipping the axis would
    # drop a dot while its mean bar still counted it. The upper limit leaves the mean labels a column of their own,
    # clear of the largest study AUROC (0.908, the health score on one read-pipeline study), whose dot the label grazed
    scale_x_continuous(limits = c(0.15, 1.09), breaks = c(0.2, 0.4, 0.6, 0.8, 1.0), expand = c(0, 0)) +
    labs(x = "AUROC in studies left out of training (bar = mean, dots = studies)", y = NULL) +
    guides(fill = guide_legend(keywidth = unit(6, "pt"), keyheight = unit(6, "pt"))) +
    theme(panel.grid.major.y = element_blank(), strip.text = element_text(size = 6.5, face = "bold", colour = INK, hjust = 0.5, margin = margin(1, 0, 3, 0)),
          axis.text.y = element_text(size = 6.5, colour = INK), panel.spacing.x = unit(9, "pt"), legend.position = "bottom",
          legend.margin = margin(-2, 0, 0, 0))
  save_fig(g, "fig3", DOUBLE, 1.95)
}

# ---- "share outside the range" bars shared by Figs 4a and 5b: every bar grows rightward from zero and the panel is split by
# direction, so no half of the axis is left empty. The bars are stacked shares, and the share the text quotes for a row is the
# one in that row's own direction, so that segment is drawn first, flush against the axis, and its far end is where the
# comparator marker for the same quantity sits. `d` therefore carries a `dir` column naming the row's direction.
CALLS <- c("below the range" = LOW, "expected but missing" = LOW_L, "above the range" = HIGH)
call_bars <- function(d) bind_rows(d %>% transmute(name, side, dir, v = low, what = "below the range"),
                                   d %>% transmute(name, side, dir, v = missing, what = "expected but missing"),
                                   d %>% transmute(name, side, dir, v = high, what = "above the range")) %>%
  filter(v > 0) %>% mutate(what = factor(what, names(CALLS)),
                           ord = as.integer(what) - 10 * (dir == "above" & what == "above the range")) %>% group_by(name) %>%
  arrange(ord, .by_group = TRUE) %>% mutate(xmax = 100 * cumsum(v), xmin = xmax - 100 * v, x = (xmin + xmax) / 2, w = xmax - xmin) %>%
  ungroup() %>% select(-ord)

# ======================================================================================== Figure 4: known disease associations
cohort_lab <- function(s) sub("^([A-Z][a-z]+)[A-Z]+_([0-9]{4})_?([ab]?)$", "\\1 \\2\\3", s)
# ============================================= Figure 3: the report view carried into a cohort and across nine cohorts
# Both panels are Fig. 1c's geometry on a percentile axis. In percentile space the reference range is 2.5--97.5 for
# every feature, so one band serves every row and a point beyond it is outside the range by construction: the thing the
# earlier version of this figure asserted with shares and shifts without ever drawing the range they departed from.
ref_band <- function(ymax, label = TRUE) {
  list(annotate("rect", xmin = 2.5, xmax = 97.5, ymin = -Inf, ymax = Inf, fill = BAND, colour = NA),
       annotate("rect", xmin = 25, xmax = 75, ymin = -Inf, ymax = Inf, fill = IQR, colour = NA),
       annotate("segment", x = 50, xend = 50, y = -Inf, yend = Inf, colour = MED, linewidth = 0.3))
}
CALL_F <- c(low = LOW, high = HIGH, within = WITHIN, expected_but_missing = LOW_L)
CALL_L <- c(low = "below the range", high = "above the range", within = "within",
            expected_but_missing = "expected but absent")

fig4 <- function() {
  # (a) one cohort, one dot per patient per feature
  cd <- read_csv(file.path(D, "f3_range_cdi.csv"), show_col_types = FALSE) %>%
    mutate(block = factor(block, c("Bacterial families", "Gene families")))
  # not_assessable is a feature whose coverage is too thin to place at all, which the method leaves unassessed rather
  # than calling; it is neither a finding nor an absence, so it is dropped rather than drawn as either
  n_na <- sum(cd$call == "not_assessable")
  cd <- filter(cd, call != "not_assessable") %>% mutate(call = factor(call, names(CALL_F)))
  stopifnot(!any(is.na(cd$call)))
  if (n_na > 0) message("  fig4a: ", n_na, " not-assessable cells dropped (too thin to place)")
  ordr <- cd %>% filter(group == "C. difficile cases") %>% group_by(block, name) %>%
    summarise(m = median(pct, na.rm = TRUE), .groups = "drop") %>% arrange(block, m) %>% pull(name)
  cases <- cd %>% filter(group == "C. difficile cases") %>% mutate(name = factor(name, unique(ordr)))
  # an absent feature has no percentile; it is drawn as an open ring at the floor so that shape, not colour alone,
  # carries the distinction, and so that the row's absences are visible rather than silently dropped
  miss <- cases %>% filter(is.na(pct)) %>% count(block, name, name = "k")
  # labelled by the side the feature was selected for, so the number states why the row is in the panel
  n_cases <- dplyr::n_distinct(cases$sample)
  share <- cases %>% group_by(block, name, side) %>%
    summarise(lost = sum(call %in% c("low", "expected_but_missing")) / n_cases,
              below = sum(call == "low") / n_cases,
              gain = sum(call == "high") / n_cases, .groups = "drop") %>%
    # a lost row is lost either by falling below the range or by being absent; the label names whichever it is, so
    # that the number beside a row is the finding that put it there
    mutate(gone = lost - below,
           lab = ifelse(side == "gained", sprintf("%.0f%% above", 100 * gain),
                        ifelse(gone > below, sprintf("%.0f%% absent", 100 * gone), sprintf("%.0f%% below", 100 * below))),
           col = ifelse(side == "gained", HIGH, LOW))
  ctrl <- cd %>% filter(group == "controls") %>% group_by(block, name) %>%
    summarise(m = median(pct, na.rm = TRUE), .groups = "drop") %>% mutate(name = factor(name, unique(ordr)))
  a <- ggplot(cases, aes(y = name)) + ref_band() +
    geom_point(data = filter(cases, !is.na(pct), call == "within"), aes(x = pct), shape = 16, size = 0.62,
               colour = WITHIN, alpha = 0.5, position = position_jitter(height = 0.2, width = 0, seed = 4)) +
    geom_point(data = filter(cases, !is.na(pct), call != "within"), aes(x = pct, fill = call), shape = 21,
               size = 1.7, stroke = 0.22, colour = "white",
               position = position_jitter(height = 0.2, width = 0, seed = 4)) +
    geom_point(data = miss, aes(x = 1.2, y = name), shape = 21, size = 1.5, stroke = 0.4,
               colour = LOW, fill = "white", inherit.aes = FALSE) +
    geom_text(data = miss, aes(x = 4.6, y = name, label = k), size = pt(5), family = FONT,
              colour = LOW, hjust = 0, inherit.aes = FALSE) +
    geom_point(data = ctrl, aes(x = m, y = name), shape = 23, size = 1.6, stroke = 0.4,
               colour = INK, fill = "white", inherit.aes = FALSE) +
    geom_text(data = share, aes(x = 103, y = name, label = lab, colour = col), hjust = 0, size = pt(5.6),
              family = FONT, inherit.aes = FALSE) +
    scale_colour_identity() +
    scale_fill_manual(values = CALL_F, name = NULL, breaks = c("low", "high"),
                      labels = CALL_L[c("low", "high")]) +
    # room to the right of the range for the per-row label; the band and the breaks still end at 97.5
    scale_x_continuous("percentile of the reference population", limits = c(0, 152),
                       breaks = c(2.5, 25, 50, 75, 97.5), labels = c("2.5", "25", "50", "75", "97.5"),
                       expand = expansion(mult = c(0.018, 0))) +
    facet_grid(block ~ ., scales = "free_y", space = "free_y", switch = "y") +
    labs(title = expression(bold("(a) One cohort: 56 "*italic("C. difficile")*" patients, one dot per patient")),
         subtitle = "diamond, the 14 controls at their median; ring and count, patients in whom the feature is absent") +
    guides(fill = guide_legend(override.aes = list(size = 2.2, stroke = 0.25))) +
    theme(legend.position = "bottom", panel.grid.major.y = element_blank(),
          axis.title.y = element_blank(), strip.placement = "outside")

  # (b) nine cohorts, one dot per cohort median
  cr <- read_csv(file.path(D, "f3_range_crc.csv"), show_col_types = FALSE) %>%
    mutate(block = factor(block, c("Taxa", "Function")))
  ordr2 <- cr %>% group_by(block, name) %>% summarise(m = median(median_pct_case, na.rm = TRUE), .groups = "drop") %>%
    arrange(block, m) %>% pull(name)
  cr <- cr %>% mutate(name = factor(name, unique(ordr2)),
                      call = ifelse(median_pct_case < 2.5, "low", ifelse(median_pct_case > 97.5, "high", "within")),
                      call = factor(call, names(CALL_F)))
  # the shift is the quantity the text quotes and the panel's reason for existing, so it is stated per row and the
  # segment is coloured by its direction; a grey segment leaves the reader to measure it off the axis
  shift <- cr %>% group_by(block, name) %>%
    summarise(m = dplyr::first(median_shift), k = dplyr::first(n_studies), .groups = "drop") %>%
    mutate(lab = sprintf("%s%.0f in %d", ifelse(m < 0, MINUS, "+"), abs(m), k), col = ifelse(m < 0, LOW, HIGH))
  cr <- left_join(cr, select(shift, block, name, dir = m), by = c("block", "name"))
  b <- ggplot(cr, aes(y = name)) + ref_band() +
    geom_segment(aes(x = median_pct_control, xend = median_pct_case, yend = name,
                     colour = ifelse(dir < 0, LOW, HIGH)), linewidth = 0.32, alpha = 0.55) +
    geom_point(aes(x = median_pct_control), shape = 23, size = 1.1, stroke = 0.28, colour = MUTED, fill = "white") +
    geom_point(aes(x = median_pct_case, colour = ifelse(dir < 0, LOW, HIGH)), shape = 16, size = 1.45) +
    geom_text(data = shift, aes(x = 103, y = name, label = lab, colour = col), hjust = 0, size = pt(5.6),
              family = FONT, inherit.aes = FALSE) +
    scale_colour_identity() +
    scale_x_continuous("median percentile of the reference population", limits = c(0, 152),
                       breaks = c(2.5, 25, 50, 75, 97.5), labels = c("2.5", "25", "50", "75", "97.5"),
                       expand = expansion(mult = c(0.018, 0))) +
    facet_grid(block ~ ., scales = "free_y", space = "free_y", switch = "y") +
    labs(title = "(b) Nine colorectal cancer cohorts: each cohort's median, cases against its own controls",
         subtitle = "diamond, a cohort's controls; dot, its cases; label, the median shift and the cohorts it holds in") +
    theme(panel.grid.major.y = element_blank(), axis.title.y = element_blank(), strip.placement = "outside")

  p <- (a / b) + plot_layout(heights = c(1.0, 0.92)) + tags()
  save_fig(p, "fig4", DOUBLE, 6.0)
}

# ======================================================================================== Figure 5: a population outside the range
fig5 <- function() {
  fr <- read_csv(file.path(D, "f5_frac.csv")) %>% filter(layer %in% c("Family", "KO (eggNOG)")) %>% mutate(pct = 100 * frac_outside_raw)
  n <- fr %>% distinct(sample, group) %>% count(group); nn <- setNames(n$n, n$group)
  st <- c("Other baseline studies" = paste0("other healthy adults\n(28 studies, n = ", comma(nn[["Other baseline studies"]]), ")"),
          "US participants, same study" = paste0("US participants of the\nHadza study (n = ", nn[["US participants, same study"]], ")"),
          "Hadza (Tanzania)" = paste0("Hadza, Tanzania\n(n = ", nn[["Hadza (Tanzania)"]], ")"))
  fr$g <- factor(st[fr$group], rev(st)); fr$layer <- factor(ifelse(fr$layer == "Family", "bacterial families", "gene families (KOs)"), c("bacterial families", "gene families (KOs)"))
  md <- fr %>% group_by(layer, g) %>% summarise(m = median(pct))
  a <- ggplot(fr, aes(x = pct, y = g)) + geom_vline(xintercept = 5, colour = MUTED, linewidth = 0.35, linetype = "22") +
    geom_boxplot(aes(fill = g), width = 0.62, outlier.size = 0.2, outlier.colour = GREY_L, outlier.stroke = 0, linewidth = 0.3, colour = INK2, median.colour = INK, median.linewidth = 0.7) +
    geom_text(data = md, aes(x = 64.5, label = paste0(sprintf("%.1f", m), "%")), hjust = 0, size = pt(6.4), family = FONT, colour = INK) +
    facet_wrap(~layer, ncol = 1) + scale_fill_manual(values = setNames(c("#8A8274", "#C6C0B0", BAND), rev(st)), guide = "none") +
    scale_x_continuous(breaks = c(0, 20, 40, 60), labels = function(x) paste0(x, "%"), expand = c(0, 0)) +
    # the window reaches the largest sample (61.2 %): a narrower one let the boxplot's own outliers escape into the right
    # margin, where they had no axis under them and read as stray ink
    coord_cartesian(xlim = c(0, 63), clip = "off") +
    # the dashed rule is named in the title, as it is in Fig. 2c, rather than labelled inside the panel: at the only y
    # that cleared the top box the label fell outside the panel altogether, into the facet strip's band, where it read as
    # a third item in the strip's own header row instead of as a note on the rule
    labs(x = "features outside the reference range, per sample", y = NULL, title = "share of each sample outside the range (5 % expected)") +
    theme(panel.grid.major.y = element_blank(), axis.text.y = element_text(colour = INK, size = 6.4, lineheight = 0.95),
          strip.text = element_text(size = 6.5, face = "plain", colour = INK2), plot.margin = margin(2, 26, 2, 2))
  hdr <- data.frame(layer = factor("bacterial families", levels(fr$layer)), g = factor(st[[3]], levels(fr$g)))
  a <- a + geom_text(data = hdr, inherit.aes = FALSE, aes(x = 64.5, y = 3.62, label = "median"), hjust = 0, size = pt(6), colour = MUTED, family = FONT)

  # b: genera and gene families together, so the panel shows both layers of the output. Bars grow from zero; colour carries the
  # direction, so no half of the axis is left empty. Genera are italic, gene families roman.
  rd <- function(f, lay, up, dn) {
    d <- read_csv(file.path(D, f)) %>% mutate(name = sub(" \\(ex .*\\)$", "", name), layer = lay, feature_id = as.character(feature_id))
    h <- filter(d, group == "Hadza (Tanzania)") %>% mutate(score = high - low - missing)
    keep <- c(h %>% filter(direction == "above") %>% slice_max(score, n = up) %>% pull(name),
              h %>% filter(direction == "below") %>% slice_min(score, n = dn) %>% pull(name))
    filter(d, name %in% keep)
  }
  ge <- bind_rows(rd("f5_genera.csv", "genera", 5, 3), rd("f5_ko.csv", "gene families", 4, 3)) %>%
    mutate(name = sub(" \\[EC:[^]]*\\]$", "", name), layer = factor(layer, c("genera", "gene families")))
  hz <- filter(ge, group == "Hadza (Tanzania)") %>% mutate(score = high - low - missing)
  ordr <- hz %>% arrange(desc(layer), score) %>% pull(name)
  ge$name <- factor(ge$name, ordr); hz$name <- factor(hz$name, ordr)
  ital <- hz$name[hz$layer == "genera"]
  bars <- call_bars(hz %>% mutate(side = layer, dir = direction)) %>% mutate(name = factor(name, ordr), layer = side)
  cmp <- ge %>% filter(group != "Hadza (Tanzania)") %>% mutate(x = ifelse(direction == "above", high, low + missing),
           who = factor(ifelse(group == "Other baseline studies", "other healthy adults", "US participants, same study"),
                        c("US participants, same study", "other healthy adults")))
  # both comparators are drawn on the row itself. A row is 7.4 pt of a figure this size and a legible marker is close to
  # 5 pt, so nudging the pair apart cost more than it bought: the two still overlapped, each hung out of its own bar, and
  # on the four rows where both sit at the same x the four marks of two rows formed a chain whose rows could not be told
  # apart. They are nested instead, the broad baseline an open circle and this study's own US participants a small solid
  # diamond, which stays readable both where the two coincide and where they are a few points apart.
  cmp_pt <- function(w, sz, fl, st) geom_point(data = filter(cmp, who == w), aes(x = 100 * x, shape = who),
                                               size = sz, colour = INK, fill = fl, stroke = st)
  # Panel (b) is Fig. 1c's geometry again: the reference band, and every Hadza sample's percentile against it, so that
  # the genera driving the excess are read off the same axis as the shift itself rather than as a share bar.
  hz <- read_csv(file.path(D, "f4_range_hadza.csv"), show_col_types = FALSE) %>%
    filter(call != "not_assessable") %>%
    mutate(block = factor(block, c("Common gut genera", "Gene families")), call = factor(call, names(CALL_F)))
  stopifnot(!any(is.na(hz$call)))
  hz_ord <- hz %>% filter(group == "Hadza (Tanzania)") %>% group_by(block, name) %>%
    summarise(m = median(pct, na.rm = TRUE), .groups = "drop") %>% arrange(block, m) %>% pull(name)
  hz <- mutate(hz, name = factor(name, unique(hz_ord)))
  hadza <- filter(hz, group == "Hadza (Tanzania)")
  hmiss <- hadza %>% filter(is.na(pct)) %>% count(block, name, name = "k")
  n_hadza <- dplyr::n_distinct(hadza$sample)
  hshare <- hadza %>% group_by(block, name) %>%
    summarise(lost = sum(call %in% c("low", "expected_but_missing")) / n_hadza,
              below = sum(call == "low") / n_hadza,
              gain = sum(call == "high") / n_hadza, .groups = "drop") %>%
    mutate(gone = lost - below,
           lab = ifelse(gain > lost, sprintf("%.0f%% above", 100 * gain),
                        ifelse(gone > below, sprintf("%.0f%% absent", 100 * gone), sprintf("%.0f%% below", 100 * below))),
           col = ifelse(gain > lost, HIGH, LOW))
  b <- ggplot(hadza, aes(y = name)) + ref_band() +
    geom_point(data = filter(hadza, !is.na(pct), call == "within"), aes(x = pct), shape = 16, size = 0.5,
               colour = WITHIN, alpha = 0.45, position = position_jitter(height = 0.22, width = 0, seed = 7)) +
    geom_point(data = filter(hadza, !is.na(pct), call != "within"), aes(x = pct, fill = call), shape = 21,
               size = 1.35, stroke = 0.18, colour = "white",
               position = position_jitter(height = 0.22, width = 0, seed = 7)) +
    geom_point(data = hmiss, aes(x = 1.2, y = name), shape = 21, size = 1.4, stroke = 0.4,
               colour = LOW, fill = "white", inherit.aes = FALSE) +
    geom_text(data = hmiss, aes(x = 4.6, y = name, label = k), size = pt(5), family = FONT,
              colour = LOW, hjust = 0, inherit.aes = FALSE) +
    geom_text(data = hshare, aes(x = 103, y = name, label = lab, colour = col), hjust = 0, size = pt(5.6),
              family = FONT, inherit.aes = FALSE) +
    scale_colour_identity() +
    scale_fill_manual(values = CALL_F, name = NULL, breaks = c("low", "within", "high"),
                      labels = CALL_L[c("low", "within", "high")]) +
    scale_y_discrete(labels = ital_lab(unique(as.character(filter(hz, block == "Common gut genera")$name)))) +
    scale_x_continuous("percentile of the reference population", limits = c(0, 155),
                       breaks = c(2.5, 25, 50, 75, 97.5), labels = c("2.5", "25", "50", "75", "97.5"),
                       expand = expansion(mult = c(0.018, 0))) +
    facet_grid(block ~ ., scales = "free_y", space = "free_y") +
    labs(y = NULL, title = "the genera and gene families behind the shift",
         subtitle = "one dot per Hadza sample; ring and count, samples lacking it") +
    guides(fill = guide_legend(override.aes = list(size = 2.2, stroke = 0.2))) +
    theme(panel.grid.major.y = element_blank(), axis.text.y = element_text(colour = INK, size = 6.5),
          legend.position = "bottom",
          strip.text.y = element_text(angle = -90, size = 6.2, face = "plain", colour = INK2, hjust = 0.5),
          panel.spacing.y = unit(5, "pt"))
  # (c) the same question asked of every population the baseline is built from, not just the one that fails loudest
  cc <- read_csv(file.path(D, "f5_countries.csv"), show_col_types = FALSE) %>%
    mutate(layer = factor(layer, c("family", "genus", "species")),
           rep = replicates == "all three layers")
  ord <- cc %>% filter(layer == "family") %>% arrange(pct) %>% pull(country)
  cc <- mutate(cc, country = factor(country, ord))
  spans <- cc %>% filter(layer == "family", studies >= 2) %>% distinct(country, study_lo, study_hi)
  cl <- ggplot(cc, aes(y = country)) +
    annotate("rect", xmin = -Inf, xmax = 5, ymin = -Inf, ymax = Inf, fill = BAND, alpha = 0.45) +
    annotate("segment", x = 5, xend = 5, y = -Inf, yend = Inf, colour = MED, linewidth = 0.3, linetype = "22") +
    # the span of a country's own study means: what a protocol alone is observed to do, and the yardstick the
    # elevations have to clear before they can be read as a property of the population
    geom_segment(data = spans, aes(x = study_lo, xend = study_hi, y = country, yend = country), colour = AXIS, linewidth = 1.1,
                 lineend = "round", inherit.aes = FALSE) +
    geom_point(aes(x = pct, shape = layer, colour = ifelse(rep, INK, MUTED), size = rep), stroke = 0.45, fill = "white") +
    scale_colour_identity() + scale_size_manual(values = c(`TRUE` = 1.5, `FALSE` = 1.05), guide = "none") +
    scale_shape_manual(values = c(family = 16, genus = 1, species = 4), name = NULL) +
    scale_x_continuous("outside the range, per sample", labels = function(x) paste0(x, "%"),
                       breaks = c(0, 5, 10), limits = c(0, 11.2), expand = expansion(0)) +
    labs(y = NULL, title = "every population the baseline is built from",
         subtitle = "grey bar spans that country's own studies") +
    guides(shape = guide_legend(override.aes = list(size = 1.4, colour = INK))) +
    theme(panel.grid.major.y = element_blank(), legend.position = "bottom",
          axis.text.y = element_text(colour = INK, size = 6.2))

  p <- ((a | cl) + plot_layout(widths = c(1, 0.92))) / b + plot_layout(heights = c(0.82, 1)) + plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(size = 9, face = "bold", family = FONT))
  save_fig(p, "fig5", DOUBLE, 5.1)
}

args <- commandArgs(trailingOnly = TRUE); if (length(args) == 0) args <- paste0("fig", 1:5)
for (f in args) get(f)()
