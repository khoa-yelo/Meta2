#!/usr/bin/env Rscript
# Figures 1-5 of the two-column paper, drawn with ggplot2 + patchwork from the tidy tables that scripts/viz2/prep.py writes to
# figures/v2/data. One theme and one palette throughout:
#   blue = below the healthy range (light blue = expected but missing), orange = above it, grey = within;
#   violet = checkgm, greys = comparators; teal / ochre = assembly / read pipeline (schematics and curation only).
# Usage: Rscript scripts/viz2/figures.R [fig1 fig2 ...]   (default: all). Output: figures/v2/figN.pdf (+ .png preview, 300 dpi).
suppressPackageStartupMessages({library(ggplot2); library(patchwork); library(dplyr); library(tidyr); library(readr); library(scales)})
options(readr.show_col_types = FALSE, dplyr.summarise.inform = FALSE)
P <- "/home/classes/bios/270/khoa/meta2/project"; D <- file.path(P, "figures/v2/data"); OUT <- file.path(P, "figures/v2")
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

# ---- schematic helpers (coordinates in inches of the panel, so circles stay round)
rr <- function(x0, x1, y0, y1, r = 0.035, n = 7) {
  a <- seq(0, pi / 2, length.out = n)
  data.frame(x = c(x1 - r + r * cos(a), x0 + r - r * sin(a), x0 + r - r * cos(a), x1 - r + r * sin(a)),
             y = c(y1 - r + r * sin(a), y1 - r + r * cos(a), y0 + r - r * sin(a), y0 + r - r * cos(a)))
}
box <- function(x0, x1, y0, y1, title, sub = NULL, fill = "white", col = INK2, lty = "solid", ts = 6.8, ss = 6.2, hj = 0.5) {
  xc <- if (hj == 0) x0 + 0.06 else (x0 + x1) / 2; yc <- (y0 + y1) / 2
  L <- list(geom_polygon(data = rr(x0, x1, y0, y1), aes(x, y), fill = fill, colour = col, linewidth = 0.35, linetype = lty, inherit.aes = FALSE))
  ns <- if (is.null(sub)) 0 else length(strsplit(sub, "\n")[[1]])
  th <- ts / 72 * 1.05; sh <- ss / 72 * 1.12; H <- th + if (ns > 0) 0.03 + ns * sh else 0
  L <- c(L, annotate("text", xc, yc + H / 2, label = title, size = pt(ts), fontface = "bold", family = FONT, colour = INK, vjust = 1, hjust = hj, lineheight = 0.95))
  if (ns > 0) L <- c(L, annotate("text", xc, yc + H / 2 - th - 0.03, label = sub, size = pt(ss), family = FONT, colour = INK2, vjust = 1, hjust = hj, lineheight = 1.0))
  L
}
arr <- function(x0, y0, x1, y1, col = INK2) annotate("segment", x = x0, y = y0, xend = x1, yend = y1, colour = col, linewidth = 0.4,
                                                      arrow = arrow(length = unit(3.2, "pt"), type = "closed", angle = 25))
seg <- function(x0, y0, x1, y1, col = INK2) annotate("segment", x = x0, y = y0, xend = x1, yend = y1, colour = col, linewidth = 0.4)
canvas <- function(w, h) ggplot() + coord_cartesian(xlim = c(0, w), ylim = c(0, h), expand = FALSE, clip = "off") + theme_void(base_family = FONT) +
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
    box(3.2, 4.88, 0.56, 0.95, "scored against healthy adults", "1,941 adults, 26 studies;\nstatistics only, no samples", NEUT_BG, PA) +
    box(3.2, 4.88, 0.08, 0.47, "scored against healthy adults", "6,494 adults, 19 studies;\nstatistics only, no samples", NEUT_BG, PB) +
    seg(4.9, 0.755, 4.99, 0.755) + seg(4.9, 0.275, 4.99, 0.275) + seg(4.99, 0.275, 4.99, 0.755) + arr(4.99, 0.515, 5.11, 0.515) +
    geom_polygon(data = rr(5.12, 7.13, 0.08, 0.95), aes(x, y), fill = "white", colour = INK2, linewidth = 0.35) +
    annotate("text", 5.21, 0.925, label = "report for every feature", size = pt(6.8), fontface = "bold", family = FONT, colour = INK, hjust = 0, vjust = 1)
  rows <- data.frame(y = c(0.665, 0.535, 0.27), lab = c("Oscillospiraceae", "Streptococcaceae", "lysine racemase"), s = c(6.03, 6.33, 6.74),
                     col = c(LOW, WITHIN, HIGH), call = c("low", "within", "high"))
  a <- a + annotate("segment", x = 5.21, xend = 7.04, y = 0.425, yend = 0.425, colour = GRID, linewidth = 0.5) +
    annotate("text", 5.21, 0.765, label = "taxa", hjust = 0, size = pt(5.8), family = FONT, colour = MUTED, fontface = "italic") +
    annotate("text", 5.21, 0.36, label = "genes and pathways", hjust = 0, size = pt(5.8), family = FONT, colour = MUTED, fontface = "italic")
  for (i in 1:3) a <- a + range_row(6.1, 6.65, rows$y[i], 6.1, 6.27, 6.47, 6.65, 6.37, rows$s[i], rows$col[i]) +
    annotate("text", 5.21, rows$y[i], label = rows$lab[i], hjust = 0, size = pt(6.2), family = FONT, colour = INK2) +
    annotate("text", 6.83, rows$y[i], label = rows$call[i], hjust = 0, size = pt(6.2), family = FONT, colour = rows$col[i], fontface = "bold")
  a <- a + annotate("text", 6.375, 0.125, label = "healthy range", size = pt(5.8), family = FONT, colour = MUTED, vjust = 0)

  # b: how one feature is placed
  m <- jnum(file.path(D, "f1_density_meta.json")); dz <- read_csv(file.path(D, "f1_density.csv"))
  dd <- density(dz$x, adjust = 1.1, n = 1024); dd <- data.frame(x = dd$x, y = dd$y / max(dd$y)); dd <- dd[dd$x > -10.5 & dd$x < 5.5, ]
  inb <- dd[dd$x >= m[["2_5"]] & dd$x <= m[["97_5"]], ]
  b <- ggplot(dd, aes(x, y)) + geom_area(data = inb, fill = BAND) + geom_line(colour = INK2, linewidth = 0.4) +
    annotate("segment", x = c(m[["2_5"]], m[["97_5"]]), xend = c(m[["2_5"]], m[["97_5"]]), y = 0, yend = 1.08, colour = MUTED, linewidth = 0.3, linetype = "22") +
    annotate("text", x = c(m[["2_5"]], m[["97_5"]]), y = 1.11, label = c("2.5th", "97.5th"), size = pt(6), family = FONT, colour = INK2, vjust = 0) +
    annotate("text", x = -0.6, y = 0.3, label = "healthy\nrange", size = pt(6.2), family = FONT, colour = INK2, lineheight = 0.95) +
    annotate("segment", x = m$sample, xend = m$sample, y = 0, yend = 1.33, colour = LOW, linewidth = 0.6) +
    annotate("point", x = m$sample, y = 0, colour = "white", fill = LOW, shape = 21, size = 2.2, stroke = 0.3) +
    annotate("text", x = m$sample + 0.3, y = 1.6, label = paste0("this sample: ", fold_lab(m$sample), "\n(", pct_lab(m$pct), ")"), hjust = 0, vjust = 1,
             size = pt(6.2), family = FONT, colour = LOW, lineheight = 0.95) +
    scale_x_continuous(breaks = c(-8, -4, 0, 4), labels = log2_axis, limits = c(-10.5, 5.5), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1.62), expand = c(0, 0)) +
    labs(x = "relative to the typical healthy adult", y = NULL, title = "placing one feature: Oscillospiraceae") +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(), panel.grid = element_blank(), panel.grid.major = element_blank())

  # c: the report for the same sample (bars clipped to the axis range)
  XL <- c(-14, 8.5); cl <- function(v) pmin(pmax(v, XL[1]), XL[2])
  r <- read_csv(file.path(D, "f1_report.csv")) %>% mutate(call = ifelse(is.na(x), "absent", call))
  r$label <- sub("^lysine racemase.*", "lysine racemase (K20707)", r$label); r$label <- sub("^ECF sigma.*", "ECF sigma factor (K03088)", r$label)
  r <- r %>% mutate(txt = ifelse(call == "absent", "absent; carried by most healthy adults", paste0(fold_lab(x), " · ", pct_lab(pct))),
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
    geom_point(data = filter(r, call == "absent"), aes(x = XL[1] + 0.35), shape = 21, colour = LOW, fill = "white", size = 1.9, stroke = 0.6) +
    geom_text(aes(x = XL[2] + 0.4, label = txt), hjust = 0, size = pt(6.2), family = FONT, colour = INK2) +
    facet_grid(layer ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_fill_manual(values = c(low = LOW, high = HIGH, within = WITHIN), guide = "none") +
    scale_x_continuous(breaks = c(-12, -8, -4, 0, 4, 8), labels = log2_axis, expand = c(0, 0)) +
    coord_cartesian(xlim = XL, clip = "off") +
    labs(x = "relative to the typical healthy adult (median)", y = NULL,
         title = expression(bold("report for one ") * bolditalic("C. difficile") * bold(" infection patient (excerpt)"))) +
    theme(panel.grid.major.y = element_blank(), axis.text.y = element_text(colour = INK, size = 6.5),
          strip.placement = "outside", strip.text.y.left = element_text(angle = 90, hjust = 0.5, size = 6.2, face = "bold", colour = INK2),
          panel.spacing.y = unit(5, "pt"),
          plot.margin = margin(2, 122, 2, 2))
  p <- wrap_elements(full = a) / ((b | c_) + plot_layout(widths = c(1, 1.32))) + plot_layout(heights = c(1.02, 1.50)) + plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(size = 9, face = "bold", family = FONT))
  save_fig(p, "fig1", DOUBLE, 2.50)
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
  cc <- read_tsv(file.path(P, "results/s14/disease_strata_B_balance.tsv"))   # the 15 studies with >= 15 cases and 15 controls (14 enter the benchmark)
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
  b <- ggplot(fu, aes(y = key, x = n, fill = pal)) + geom_col(width = 0.72, colour = NA) +
    geom_text(aes(label = lab), hjust = -0.08, size = pt(6), family = FONT, colour = INK2) +
    facet_wrap(~pipeline, ncol = 1, scales = "free_y") +
    scale_fill_manual(values = fills, guide = "none") + scale_y_discrete(labels = setNames(fu$label, fu$key)) +
    scale_x_continuous(labels = function(x) ifelse(x == 0, "0", paste0(x / 1000, "k")), limits = c(0, 33500), breaks = c(0, 10000, 20000), expand = c(0, 0)) +
    labs(x = "samples (studies)", y = NULL, title = "curation") +
    theme(panel.grid.major.y = element_blank(), axis.text.y = element_text(size = 6.2, colour = INK), strip.text = element_text(size = 6.5, face = "plain", colour = INK2),
          panel.spacing.y = unit(7, "pt"))

  # c: calibration, share of features outside the range in healthy samples the baseline never saw
  ca <- read_csv(file.path(D, "f2_calibration.csv")) %>% mutate(frac = 100 * frac, layer = factor(layer, rev(LAYERS))) %>% filter(!is.na(layer)) %>%
    mutate(kind = factor(kind, c("held-out cohort", "baseline study left out")),
           pipeline = recode(pipeline, "Assembly (MGnify v5)" = "assembly pipeline", "Read (cMD3)" = "read pipeline"))
  med <- ca %>% group_by(pipeline, layer) %>% summarise(m = median(frac))
  hz <- ca %>% filter(unit == "PRJEB49206", layer == "Family")
  c_ <- ggplot(ca, aes(x = frac, y = layer)) + geom_vline(xintercept = 5, colour = MUTED, linewidth = 0.35, linetype = "22") +
    geom_point(aes(shape = kind), position = position_jitter(height = 0.18, width = 0, seed = 1), size = 0.9, colour = INK2, stroke = 0.35, alpha = 0.8) +
    geom_point(data = med, aes(x = m), shape = 124, size = 3.2, colour = ACC) +
    geom_text(data = hz, aes(label = "Hadza study (Fig. 5)"), nudge_y = 0.62, size = pt(6), family = FONT, colour = INK2, hjust = 0.75) +
    facet_wrap(~pipeline, ncol = 1, scales = "free_y") + scale_shape_manual(values = c(16, 1), name = NULL) +
    scale_y_discrete(expand = expansion(add = c(0.6, 1.0))) +
    scale_x_continuous(limits = c(0, 25), breaks = c(0, 5, 10, 15, 20, 25), labels = function(x) paste0(x, "%"), expand = c(0, 0.3)) +
    labs(x = "features outside the healthy range", y = NULL, title = "calibration (5 % expected)") +
    theme(legend.position = "bottom", legend.text = element_text(size = 6), panel.grid.major.y = element_blank(), strip.text = element_text(size = 6.5, face = "plain", colour = INK2))

  # d: the pipelines reproduce the source percentiles from raw reads
  fi <- read_csv(file.path(D, "f2_fidelity.csv")) %>% mutate(layer = factor(layer, rev(LAYERS))) %>% filter(!is.na(layer)) %>%
    mutate(pipeline = recode(pipeline, "Assembly (MGnify v5)" = "assembly pipeline, 61 samples", "Read (cMD3)" = "read pipeline, 34 samples"))
  fm <- fi %>% group_by(pipeline, layer) %>% summarise(m = median(rho), n = n())
  d <- ggplot(fi, aes(x = rho, y = layer)) +
    geom_point(position = position_jitter(height = 0.18, width = 0, seed = 1), size = 0.8, colour = INK2, alpha = 0.55, stroke = 0) +
    geom_point(data = fm, aes(x = m), shape = 124, size = 3.2, colour = ACC) +
    geom_text(data = fm, aes(x = 0.405, label = sprintf("%.2f", m)), hjust = 0, size = pt(6), family = FONT, colour = ACC) +
    facet_wrap(~pipeline, ncol = 1, scales = "free_y") +
    scale_x_continuous(limits = c(0.4, 1.0), breaks = c(0.4, 0.6, 0.8, 1.0), expand = c(0, 0.01)) +
    labs(x = "Spearman ρ with source percentiles", y = NULL, title = "reproduced from raw reads") +
    theme(panel.grid.major.y = element_blank(), strip.text = element_text(size = 6.5, face = "plain", colour = INK2))
  p <- wrap_elements(full = a) / (b | c_ | d) + plot_layout(heights = c(1.06, 2.45)) + plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(size = 9, face = "bold", family = FONT))
  save_fig(p, "fig2", DOUBLE, 2.9)
}

# ======================================================================================== Figure 3: disease benchmark
fig3 <- function() {
  au <- read_csv(file.path(D, "f3_auroc.csv"))
  au$method <- recode(au$method, "GMHI (genus approx.)" = "GMHI, genus approx.", "GMHI (published)" = "GMHI", "GMWI2 (published)" = "GMWI2",
                      "Alpha diversity" = "alpha diversity", "Raw abundances" = "raw abundances",
                      "checkgm percentiles" = "checkgm percentiles, taxa", "checkgm percentiles, genes only" = "checkgm percentiles, genes",
                      "checkgm percentiles, taxa + genes" = "checkgm percentiles, taxa + genes")
  lv <- c("assembly pipeline, 11 studies", "read pipeline, 14 studies", "read pipeline, the 4 GMWI2 never saw")
  au$pipeline <- factor(recode(au$pipeline, "Assembly pipeline (11 studies)" = lv[1], "Read pipeline (14 studies)" = lv[2], "Read pipeline, 4 studies GMWI2 never saw" = lv[3]), lv)
  s <- au %>% group_by(pipeline, method) %>% summarise(m = mean(auroc), n = n()) %>% ungroup() %>%
    mutate(lab = ifelse(method == "GMWI2" & pipeline == lv[2], "GMWI2†", method), key = paste(pipeline, method),
           grp = ifelse(grepl("^checkgm", method), "checkgm", ifelse(method == "GMWI2", "GMWI2", "other")))
  s <- s %>% arrange(pipeline, m) %>% mutate(key = factor(key, key))
  au <- au %>% mutate(key = factor(paste(pipeline, method), levels(s$key)))
  g <- ggplot(s, aes(y = key)) + geom_vline(xintercept = 0.5, colour = MUTED, linewidth = 0.35) +
    geom_tile(aes(x = (0.5 + m) / 2, width = abs(m - 0.5), fill = grp), height = 0.7) +
    geom_point(data = au, aes(x = auroc), position = position_jitter(height = 0.17, width = 0, seed = 2), size = 0.6, colour = INK, alpha = 0.5, stroke = 0) +
    geom_text(aes(x = 1.075, label = sprintf("%.2f", m), fontface = ifelse(grp == "checkgm", "bold", "plain")), hjust = 1, size = pt(6.2), family = FONT, colour = INK) +
    facet_wrap(~pipeline, ncol = 1, scales = "free_y") + scale_fill_manual(values = c(checkgm = ACC, GMWI2 = GREY_D, other = GREY_L), guide = "none") +
    scale_y_discrete(labels = setNames(s$lab, s$key)) +
    scale_x_continuous(limits = c(0.25, 1.08), breaks = seq(0.3, 1.0, 0.1), labels = c("0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1.0"), expand = c(0, 0)) +
    labs(x = "AUROC in studies left out of training", y = NULL) +
    theme(panel.grid.major.y = element_blank(), strip.text = element_text(size = 6.5, face = "bold", colour = INK, margin = margin(3, 0, 2, 0)),
          axis.text.y = element_text(size = 6.5, colour = INK), panel.spacing.y = unit(4, "pt"))
  gt <- ggplot_gtable(ggplot_build(g))                                     # facet heights proportional to the number of bars
  rows <- gt$layout$t[grepl("panel", gt$layout$name)]; gt$heights[rows] <- unit(as.numeric(table(s$pipeline)), "null")
  ggsave(file.path(OUT, "fig3.pdf"), gt, width = SINGLE, height = 3.1, units = "in", device = cairo_pdf)
  ggsave(file.path(OUT, "fig3.png"), gt, width = SINGLE, height = 3.1, units = "in", dpi = 300, device = ragg::agg_png, bg = "white"); message("wrote fig3")
}

# ---- "share outside the range" bars shared by Figs 4a and 5b: every bar grows rightward from zero and the panel is split by
# direction, so no half of the axis is left empty
CALLS <- c("below the range" = LOW, "expected but missing" = LOW_L, "above the range" = HIGH)
call_bars <- function(d) bind_rows(d %>% transmute(name, side, v = low, what = "below the range"),
                                   d %>% transmute(name, side, v = missing, what = "expected but missing"),
                                   d %>% transmute(name, side, v = high, what = "above the range")) %>%
  filter(v > 0) %>% mutate(what = factor(what, names(CALLS))) %>% group_by(name) %>%
  arrange(what, .by_group = TRUE) %>% mutate(xmax = 100 * cumsum(v), xmin = xmax - 100 * v, x = (xmin + xmax) / 2, w = xmax - xmin) %>% ungroup()

# ======================================================================================== Figure 4: known disease associations
cohort_lab <- function(s) sub("^([A-Z][a-z]+)[A-Z]+_([0-9]{4})_?([ab]?)$", "\\1 \\2\\3", s)
fig4 <- function() {
  cd <- read_csv(file.path(D, "f4_cdi_calls.csv")) %>% rename(missing = expected_but_missing) %>%
    mutate(side = factor(ifelse(side == "lost", "lost", "gained"), c("lost", "gained")))
  cs <- filter(cd, group == "cases"); ct <- filter(cd, group == "controls")
  ordr <- c(cs %>% filter(side == "lost") %>% arrange(low + missing) %>% pull(name), cs %>% filter(side == "gained") %>% arrange(high) %>% pull(name))
  bars <- call_bars(cs) %>% mutate(name = factor(name, ordr))
  mk <- ct %>% transmute(name = factor(name, ordr), side, x = ifelse(side == "gained", high, low + missing))
  a <- ggplot(bars, aes(y = name)) +
    geom_tile(aes(x = x, width = w, fill = what), height = 0.68) +
    geom_point(data = mk, aes(x = 100 * x, shape = "controls (n = 14)"), size = 1.5, colour = INK, fill = "white", stroke = 0.45) +
    facet_grid(side ~ ., scales = "free_y", space = "free_y") +
    scale_fill_manual(values = CALLS, name = NULL) + scale_shape_manual(values = c("controls (n = 14)" = 23), name = NULL) +
    scale_x_continuous(limits = c(0, 68), breaks = c(0, 20, 40, 60), labels = function(x) paste0(x, "%"), expand = c(0, 0)) +
    labs(x = "patients outside the healthy range (cases, n = 56)", y = NULL,
         title = expression(bolditalic("C. difficile") * bold(" infection, one cohort (assembly pipeline)"))) +
    guides(fill = guide_legend(order = 1, nrow = 1, keywidth = unit(6, "pt"), keyheight = unit(6, "pt")), shape = guide_legend(order = 2)) +
    theme(legend.position = "bottom", legend.box = "horizontal", legend.box.just = "top", legend.spacing.x = unit(5, "pt"),
          legend.justification = "center", strip.text.y = element_text(angle = -90, size = 6.2, face = "plain", colour = INK2, hjust = 0.5),
          panel.spacing.y = unit(5, "pt"), axis.text.y = element_text(colour = INK, size = 6.5), panel.grid.major.y = element_blank())

  short <- function(x) {                                    # keep the long MetaCyc and KO names readable at 6 pt
    x <- sub(" \\[EC:[^]]*\\]$", "", x)
    x <- sub("^two-component system, chemotaxis family, protein-glutamate.*", "chemotaxis methylesterase (CheB)", x)
    x <- sub("^energy-coupling factor transport system.*", "energy-coupling factor transporter", x)
    x <- sub("^ATP phosphoribosyltransferase regulatory subunit$", "ATP phosphoribosyltransferase (reg.)", x)
    x <- sub("^thiamin formation from pyrithiamine and oxythiamine \\(yeast\\)$", "thiamin formation", x)
    x <- sub("^myo-, chiro- and scillo-inositol degradation$", "inositol degradation", x)
    x <- sub("^pyruvate fermentation to acetate and lactate II$", "pyruvate fermentation to acetate", x)
    x <- sub("^pentose phosphate pathway \\(non-oxidative branch\\)$", "pentose phosphate (non-oxidative)", x)
    x <- sub("^methylerythritol phosphate pathway II$", "methylerythritol phosphate pathway", x)
    x
  }
  keep_n <- c(Family = 4, Genus = 5)                        # room for the function rows; the full lists are in the atlas tables
  cr <- read_csv(file.path(D, "f4_crc.csv")) %>% filter(rank != "Species")
  topt <- cr %>% distinct(rank, name, median_shift) %>% group_by(rank) %>% slice_min(median_shift, n = 5) %>% ungroup() %>%
    filter(rank != "Family" | name %in% (cr %>% distinct(rank, name, median_shift) %>% filter(rank == "Family") %>% slice_min(median_shift, n = 4) %>% pull(name)))
  cr <- filter(cr, name %in% topt$name)
  fn <- read_csv(file.path(D, "f4_crc_function.csv")) %>% mutate(name = short(name)) %>%
    group_by(rank) %>% filter(name %in% (distinct(., name, median_shift) %>% slice_min(median_shift, n = 4) %>% pull(name))) %>% ungroup()
  cr <- bind_rows(cr %>% mutate(rank = ifelse(rank == "Family", "families", "genera")),
                  fn %>% mutate(rank = ifelse(rank == "gene families", "genes", rank))) %>% mutate(cohort = cohort_lab(study))
  ital <- cr$name[cr$rank == "genera"]
  ordt <- cr %>% distinct(rank, name, median_shift) %>% mutate(rank = factor(rank, c("families", "genera", "pathways", "genes"))) %>% arrange(rank, desc(median_shift))
  cr$name <- factor(cr$name, ordt$name); cr$rank <- factor(cr$rank, c("families", "genera", "pathways", "genes"))
  cr$cohort <- factor(cr$cohort, sort(unique(cr$cohort)))
  ylab_it <- function(br) as.expression(lapply(br, function(n) if (n %in% ital) bquote(italic(.(n))) else bquote(.(n))))
  sg <- filter(cr, sig)
  b <- ggplot(cr, aes(x = cohort, y = name, fill = shift)) + geom_tile(colour = "white", linewidth = 0.5) +
    geom_point(data = sg, size = 0.75, colour = ifelse(abs(sg$shift) > 30, "white", INK)) +
    facet_grid(rank ~ ., scales = "free_y", space = "free_y") +
    scale_fill_gradient2(low = LOW, mid = "#F2F2F2", high = HIGH, midpoint = 0, limits = c(-50, 50), oob = squish, breaks = c(-50, -25, 0, 25, 50),
                         labels = c(paste0(MINUS, "50"), paste0(MINUS, "25"), "0", "25", "50"), name = "cases − controls\n(percentile points)") +
    scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0), labels = ylab_it) +
    labs(x = NULL, y = NULL, title = "colorectal cancer, nine cohorts (read pipeline)") +
    guides(fill = guide_colourbar(barwidth = unit(4.5, "pt"), barheight = unit(60, "pt"), ticks.colour = "white")) +
    theme(axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1, size = 6.2, colour = INK), panel.grid.major = element_blank(), axis.line = element_blank(),
          axis.ticks = element_blank(), strip.text.y = element_text(angle = -90, size = 6.2, face = "plain", colour = INK2, hjust = 0.5),
          axis.text.y = element_text(colour = INK, size = 6.2), legend.title = element_text(size = 6, colour = INK2), legend.position = "right",
          panel.spacing.y = unit(2.5, "pt"))
  p <- (a | b) + plot_layout(widths = c(0.86, 1.0)) + plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(size = 9, face = "bold", family = FONT))
  save_fig(p, "fig4", DOUBLE, 2.76)
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
    geom_text(data = md, aes(x = 46.5, label = paste0(sprintf("%.1f", m), "%")), hjust = 0, size = pt(6.4), family = FONT, colour = INK) +
    facet_wrap(~layer, ncol = 1) + scale_fill_manual(values = setNames(c("#8A8274", "#C6C0B0", BAND), rev(st)), guide = "none") +
    scale_x_continuous(breaks = c(0, 10, 20, 30, 40), labels = function(x) paste0(x, "%"), expand = c(0, 0)) +
    coord_cartesian(xlim = c(0, 45), clip = "off") +
    labs(x = "features outside the healthy range, per sample", y = NULL, title = "share of each sample outside the range") +
    theme(panel.grid.major.y = element_blank(), axis.text.y = element_text(colour = INK, size = 6.4, lineheight = 0.95),
          strip.text = element_text(size = 6.5, face = "plain", colour = INK2), plot.margin = margin(2, 22, 2, 2))
  a <- a + geom_text(data = data.frame(layer = factor("bacterial families", levels(fr$layer)), g = factor(st[[3]], levels(fr$g))), inherit.aes = FALSE,
                     aes(x = 46.5, y = 3.62, label = "median"), hjust = 0, size = pt(6), colour = MUTED, family = FONT)

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
  bars <- call_bars(hz %>% mutate(side = layer)) %>% mutate(name = factor(name, ordr), layer = side)
  cmp <- ge %>% filter(group != "Hadza (Tanzania)") %>% mutate(x = ifelse(direction == "above", high, low + missing),
           who = factor(ifelse(group == "Other baseline studies", "other healthy adults", "US participants, same study"),
                        c("US participants, same study", "other healthy adults")))
  ylab_it <- function(br) as.expression(lapply(br, function(n) if (n %in% ital) bquote(italic(.(n))) else bquote(.(n))))
  b <- ggplot(bars, aes(y = name)) +
    geom_tile(aes(x = x, width = w, fill = what), height = 0.68) +
    geom_point(data = cmp, aes(x = 100 * x, shape = who), size = 1.5, colour = INK, fill = "white", stroke = 0.45) +
    facet_grid(layer ~ ., scales = "free_y", space = "free_y") +
    scale_fill_manual(values = CALLS, name = "Hadza samples") +
    scale_shape_manual(values = c("US participants, same study" = 23, "other healthy adults" = 21), name = NULL) +
    scale_y_discrete(labels = ylab_it) +
    scale_x_continuous(limits = c(0, 78), breaks = c(0, 20, 40, 60), labels = function(x) paste0(x, "%"), expand = c(0, 0)) +
    labs(x = "samples outside the healthy range", y = NULL, title = "the genera and gene families behind the shift") +
    guides(fill = guide_legend(order = 1, ncol = 1, keywidth = unit(6, "pt"), keyheight = unit(6, "pt")), shape = guide_legend(order = 2, ncol = 1)) +
    theme(panel.grid.major.y = element_blank(), axis.text.y = element_text(colour = INK, size = 6.5), legend.position = "right",
          strip.text.y = element_text(angle = -90, size = 6.2, face = "plain", colour = INK2, hjust = 0.5), panel.spacing.y = unit(5, "pt"))
  p <- (a | b) + plot_layout(widths = c(1, 1.42)) + plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(size = 9, face = "bold", family = FONT))
  save_fig(p, "fig5", DOUBLE, 2.30)
}

args <- commandArgs(trailingOnly = TRUE); if (length(args) == 0) args <- paste0("fig", 1:5)
for (f in args) get(f)()
