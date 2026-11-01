# =============================================================================
# Fig. 7  Adjusted R2 of ordinary least squares regressions of observed
#         full-community species richness on each dominance metric alone, at
#         three spatial scales for grazing (top row), fencing (middle row), and
#         all plots (bottom row).
#         a, dominance metrics from field species abundance (OB_Top*)
#         b, dominance metrics from the deep learning species maps (DL_Top*)
#         c, dominance metrics from the spectral clusters (SC_Top*, k-means on
#            z-scored reflectance)
#
# Steps
#   1. For each source, subset and scale, regress observed SR on each of the
#      five dominance metrics (Top1, Top2, Top3, Top12, Top123) and record
#      adjusted R2 and the P value of the slope.
#   2. Plot adjusted R2 against scale, one line per metric; adjusted R2 below
#      zero is drawn as an open triangle on the baseline. Asterisks mark
#      significant slopes and are placed automatically next to their points,
#      avoiding points, lines and other asterisks.
#   3. Export the single panels and the assembled figure.
#
# Input  : analysis_table.csv (columns Scale, Management, Observed_SR,
#          OB_Top*, DL_Top*, SC_Top*)
# Output : output/Fig07_statistics.csv            adjusted R2 and P of every fit
#          output/Fig07_<a|b|c>_<subset>.tiff     single panels (6.2 x 5.6 in)
#          output/Fig07.tiff                      assembled figure
#
# Usage  : Rscript Fig07_single_metrics.R [data_file] [output_dir]
# Default: data_file  = "../data/simulated/analysis_table.csv"
#          output_dir = "output"
# =============================================================================

source("functions.R")
suppressPackageStartupMessages(library(patchwork))

args      <- commandArgs(trailingOnly = TRUE)
data_file <- if (length(args) >= 1) args[1] else file.path(DATA_DIR, "analysis_table.csv")
out_dir   <- if (length(args) >= 2) args[2] else "output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

dat <- read_table(data_file)

SOURCES <- list(a = list(prefix = "OB", title = "Observed-Top"),
                b = list(prefix = "DL", title = "Deep Learning-Top"),
                c = list(prefix = "SC", title = "Spectral Clustering-Top"))
SUBSETS    <- c("Grazing", "Fencing", "Overall")          # rows, top to bottom
TOP_LEVELS <- c("Top1", "Top2", "Top3", "Top12", "Top123")
SCALE_LABELS <- c("0.5m" = "0.5m×0.5m", "1m" = "1m×1m", "1.5m" = "1.5m×1.5m")

# ---- 1. statistics ------------------------------------------------------------------
res <- do.call(rbind, lapply(names(SOURCES), function(k) {
  do.call(rbind, lapply(SUBSETS, function(s) {
    do.call(rbind, lapply(SCALES, function(sc) {
      d <- dat[dat$Scale == sc & (s == "Overall" | dat$Management == s), ]
      do.call(rbind, lapply(TOP_LEVELS, function(t) {
        v <- paste0(SOURCES[[k]]$prefix, "_", t)
        f <- ols_fit(y = d$Observed_SR, x = d[[v]])
        data.frame(Panel = k, Source = SOURCES[[k]]$prefix, Management = s, Scale = sc,
                   Metric = t, n = f[["n"]], adj_R2 = f[["adj_R2"]], P = f[["P"]])
      }))
    }))
  }))
}))
res$Star <- vapply(res$P, p_to_stars, character(1))
write.csv(res, file.path(out_dir, "Fig07_statistics.csv"), row.names = FALSE)

# ---- 2. panels --------------------------------------------------------------------------
BASE_SIZE  <- 25
TOP_COLS   <- c(Top1 = "#D16E3D", Top2 = "#C0E3DD", Top3 = "#A19700",
                Top12 = "#98B764", Top123 = "#FDD478")
STAR_SIZE  <- 29                                  # asterisk font size (pt)
POINT_SIZE <- 5.85                                # filled points
pd <- position_dodge(width = 0.1)                 # separates the five lines slightly

# y tick interval per panel as in the published figure; another interval is
# chosen when it would give fewer than three or more than six ticks
step_published <- function(src, sub) {
  if (src == "DL") return(0.2)
  if (src == "OB") return(0.1)
  if (sub == "Fencing") return(0.01)
  0.05
}
y_step <- function(src, sub, y_top) {
  s <- step_published(src, sub)
  n <- floor(y_top / s) + 1
  if (n >= 3 && n <= 6) return(s)
  cand <- c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5)
  cand[which(floor(y_top / cand) + 1 <= 6)][1]
}

# Asterisks: each asterisk of a significant slope is placed at the first of a
# series of candidate positions around its point (above, below, right, left,
# diagonals, then further away) at which it overlaps no point, no line, no
# other asterisk and stays inside the panel. Positions are computed in inches
# on the panel; panel_in gives the panel size in the assembled figure.
place_stars <- function(p, panel_in, star_pt = STAR_SIZE) {
  b   <- ggplot_build(p)
  pp  <- b$layout$panel_params[[1]]
  xr  <- if (!is.null(pp$x.range)) pp$x.range else pp$x$continuous_range
  yr  <- if (!is.null(pp$y.range)) pp$y.range else pp$y$continuous_range
  W <- panel_in[1]; H <- panel_in[2]
  to_in  <- function(x, y) cbind((x - xr[1]) / diff(xr) * W, (y - yr[1]) / diff(yr) * H)
  lines  <- b$data[[2]]                             # geom_line
  pts    <- do.call(rbind, lapply(b$data[3:4], function(z)          # points and triangles
    if (nrow(z) > 0) z[, c("x", "y", "colour")] else NULL))
  P  <- to_in(pts$x, pts$y)
  L  <- lapply(split(lines, lines$group), function(g) { g <- g[order(g$x), ]; to_in(g$x, g$y) })
  # glyph geometry of the asterisks (calibrated on bold Arial)
  em     <- star_pt / 72
  glyph_h <- 0.185 / 0.5 * em                       # height of the asterisk glyph
  glyph_w <- 0.193 / 0.5 * em                       # width per asterisk
  rise    <- 0.0792 / 0.5 * em                      # glyph centre above the text anchor
  r_pt   <- 0.09                                    # point radius + margin
  r_line <- 0.045                                   # half line width + margin
  # number of conflicts of a label box centred at (cx, cy)
  conflicts <- function(cx, cy, w, h, boxes) {
    x0 <- cx - w / 2; x1 <- cx + w / 2; y0 <- cy - h / 2; y1 <- cy + h / 2
    n <- 0
    if (x0 < 0.03 || x1 > W - 0.03 || y0 < 0.03 || y1 > H - 0.03) n <- n + 10
    near <- function(px, py, r) {
      dx <- pmax(x0 - px, 0, px - x1); dy <- pmax(y0 - py, 0, py - y1)
      sum(dx^2 + dy^2 < r^2)
    }
    n <- n + 2 * near(P[, 1], P[, 2], r_pt)
    t <- seq(0, 1, length.out = 40)
    for (g in L) for (j in seq_len(nrow(g) - 1))
      n <- n + (near(g[j, 1] + t * (g[j + 1, 1] - g[j, 1]), g[j, 2] + t * (g[j + 1, 2] - g[j, 2]), r_line) > 0)
    gap <- 0.03                                     # minimum gap between asterisks (in)
    for (bx in boxes) n <- n + 3 * (x0 < bx[2] + gap && x1 > bx[1] - gap && y0 < bx[4] + gap && y1 > bx[3] - gap)
    n
  }
  stars <- p$data[p$data$star_lab != "" & !p$data$below_zero, ]
  if (nrow(stars) == 0) return(NULL)
  # match each significant metric to its dodged point
  key_pts <- paste(round(pts$x), pts$colour)
  stars$k <- match(paste(as.integer(stars$Scale), TOP_COLS[as.character(stars$Metric)]), key_pts)
  stars <- stars[order(-stars$y_plot), ]            # upper points first
  boxes <- list(); out <- list()
  for (i in seq_len(nrow(stars))) {
    c0 <- P[stars$k[i], ]
    w  <- glyph_w * nchar(stars$star_lab[i]); h <- glyph_h
    dv <- r_pt + 0.03 + h / 2; dh <- r_pt + 0.03 + w / 2
    # candidates in order of preference: above, below, right, left, diagonals,
    # then the same directions further away
    dirs <- rbind(c(0, 1), c(0, -1), c(1, 0), c(-1, 0), c(1, 1), c(1, -1), c(-1, 1), c(-1, -1),
                  c(0.5, 1), c(0.5, -1), c(-0.5, 1), c(-0.5, -1))
    cand <- do.call(rbind, lapply(c(0, 0.12, 0.25), function(e)
      cbind(dirs[, 1] * (dh + e), dirs[, 2] * (dv + e))))
    nc <- vapply(seq_len(nrow(cand)), function(j)
      conflicts(c0[1] + cand[j, 1], c0[2] + cand[j, 2], w, h, boxes), numeric(1))
    j  <- which(nc == min(nc))[1]
    cx <- c0[1] + cand[j, 1]; cy <- c0[2] + cand[j, 2]
    boxes[[length(boxes) + 1]] <- c(cx - w / 2, cx + w / 2, cy - h / 2, cy + h / 2)
    out[[length(out) + 1]] <- data.frame(
      x = xr[1] + cx / W * diff(xr), y = yr[1] + (cy - rise) / H * diff(yr),
      label = stars$star_lab[i], Metric = stars$Metric[i])
  }
  do.call(rbind, out)
}

r2_lines_panel <- function(k, s, x_labels, y_title, top_strip, right_strip, panel_in) {
  d <- res[res$Panel == k & res$Management == s, ]
  d$Scale  <- factor(SCALE_LABELS[d$Scale], levels = SCALE_LABELS)
  d$Metric <- factor(d$Metric, levels = TOP_LEVELS)
  d$below_zero <- d$adj_R2 < 0
  d$y_plot   <- pmax(d$adj_R2, 0)
  d$star_lab <- ifelse(d$Star %in% c("*", "**", "***"), d$Star, "")
  d$Column <- SOURCES[[k]]$title; d$Row <- s
  ymax  <- max(d$y_plot, na.rm = TRUE); if (!is.finite(ymax) || ymax <= 0) ymax <- 0.1
  y_top <- ymax * 1.15
  step  <- y_step(SOURCES[[k]]$prefix, s, y_top)
  strip_el <- element_rect(fill = STRIP_FILL, colour = STRIP_LINE, linewidth = 0.8)
  p <- ggplot(d, aes(x = Scale, y = y_plot, group = Metric, color = Metric)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.5) +
    geom_line(linewidth = 2.0, position = pd) +
    geom_point(data = ~ subset(.x, !below_zero), size = POINT_SIZE, position = pd) +
    geom_point(data = ~ subset(.x, below_zero), shape = 25, fill = "white",
               size = 5, stroke = 1.5, position = pd, show.legend = FALSE) +
    scale_color_manual(values = TOP_COLS, breaks = TOP_LEVELS, name = NULL) +
    scale_y_continuous(limits = c(0, y_top), breaks = seq(0, y_top, by = step),
                       expand = expansion(mult = c(0.05, 0.05)),
                       labels = scales::number_format(accuracy = 0.01),
                       oob = scales::oob_keep) +           # keeps asterisks placed in the axis margins
    facet_grid(Row ~ Column) +
    labs(x = if (x_labels) "Scale" else NULL, y = if (y_title) "R²" else NULL) +
    theme_nature(BASE_SIZE) +
    theme(legend.position = "none",
          axis.text.x = if (x_labels) element_text(size = BASE_SIZE - 2, face = "bold") else element_blank(),
          strip.text.x = if (top_strip) element_text(size = STRIP_SIZE, face = "bold", margin = ggplot2::margin(10, 8, 10, 8)) else element_blank(),
          strip.background.x = if (top_strip) strip_el else element_blank(),
          strip.text.y = if (right_strip) element_text(size = STRIP_SIZE, face = "bold", angle = -90, margin = ggplot2::margin(8, 10, 8, 10)) else element_blank(),
          strip.background.y = if (right_strip) strip_el else element_blank(),
          plot.background = element_blank(),
          plot.margin = ggplot2::margin(3, 6, 3, 6))
  st <- place_stars(p, panel_in)
  if (!is.null(st))
    p <- p + geom_text(data = st, aes(x = x, y = y, label = label, color = Metric), inherit.aes = FALSE,
                       size = STAR_SIZE / .pt, fontface = "bold", family = "Arial", show.legend = FALSE)
  p
}

STRIP_LINE <- "#969696"
STRIP_SIZE <- 27
FIG_W <- 18.0; FIG_H <- 16.2                       # assembled figure (in)
PANEL_IN <- c(FIG_W * 2139 / 8000, FIG_H * 2030 / 7207)   # panel size in the assembled figure

panels <- list()
for (s in SUBSETS) for (k in names(SOURCES)) {
  p <- r2_lines_panel(k, s, x_labels = s == "Overall", y_title = k == "a",
                      top_strip = s == "Grazing", right_strip = k == "c", panel_in = PANEL_IN)
  if (s == "Grazing") p <- p + labs(tag = k)
  panels[[length(panels) + 1]] <- p
  save_tiff(p, file.path(out_dir, sprintf("Fig07_%s_%s.tiff", k, s)), 6.2, 5.6)
}

# ---- 3. assembled figure -------------------------------------------------------------
LEGEND_SIZE <- 25
LETTER_SIZE <- 40
fig <- wrap_plots(panels, ncol = 3) +
  plot_layout(guides = "collect") &
  guides(color = guide_legend(nrow = 1, override.aes = list(linewidth = 2.0, size = 5))) &
  theme(legend.position = "bottom",
        legend.text = element_text(size = LEGEND_SIZE), legend.background = element_blank(),
        legend.key.width = unit(60, "pt"),
        plot.tag = element_text(size = LETTER_SIZE, face = "bold", family = "Arial"),
        plot.tag.position = c(0.01, 0.99))
fig[[1]] <- fig[[1]] + theme(plot.tag.position = c(0.02, 0.99))
save_tiff(fig, file.path(out_dir, "Fig07.tiff"), FIG_W, FIG_H)

cat("Fig. 7 written to", out_dir, "\n")
