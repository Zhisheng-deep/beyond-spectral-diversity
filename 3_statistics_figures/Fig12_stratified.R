# =============================================================================
# Fig. 12  Relationships between observed full-community species richness (SR)
#          and its two direct remotely sensed predictors in the low and high
#          strata of four ecological factors, at three spatial scales.
#          a, SR predicted by deep learning (DL_SR)
#          b, spectral species richness from spectral clustering (SSR_km_z,
#             k-means on z-scored reflectance)
#
# Steps
#   1. For each predictor, scale and factor (Height_CV, Rare_Species_Ratio,
#      Intra_Species_SDiv and Turnover, which is undefined at 0.5 m), split all
#      observations at the median of the factor: the observations are ordered
#      by the factor and the lower floor(n/2) form the low stratum, the others
#      the high stratum.
#   2. Within each stratum, regress observed SR on the predictor with ordinary
#      least squares and record adjusted R2, the slope with its standard error
#      and the P value of the slope.
#   3. Plot observed SR against the predictor with the least squares lines and
#      their 95% confidence bands of the two strata, one panel per factor and
#      scale, and report adjusted R2 (R2 < 0 when below zero) and slope +/-
#      standard error of each stratum. Export the single panels and the two
#      assembled figures.
#
# Input  : analysis_table.csv (columns Scale, Observed_SR, DL_SR, SSR_km_z,
#          Height_CV, Rare_Species_Ratio, Intra_Species_SDiv, Turnover)
# Output : output/Fig12_statistics.csv                 statistics of every fit
#          output/Fig12_<a|b>_<scale>_<factor>.tiff     single panels (4.65 x 4.2 in)
#          output/Fig12a.tiff, Fig12b.tiff              assembled figures
#
# Usage  : Rscript Fig12_stratified.R [data_file] [output_dir]
# Default: data_file  = "../data/simulated/analysis_table.csv"
#          output_dir = "output"
# =============================================================================

source("functions.R")

args      <- commandArgs(trailingOnly = TRUE)
data_file <- if (length(args) >= 1) args[1] else file.path(DATA_DIR, "analysis_table.csv")
out_dir   <- if (length(args) >= 2) args[2] else "output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

dat <- read_table(data_file)

PREDICTORS <- list(a = list(var = "DL_SR",    x_title = "Predicted Species Richness", header = "Deep Learning"),
                   b = list(var = "SSR_km_z", x_title = "Spectral Species Richness",  header = "Spectral Clustering"))
FACTORS <- c(Height_CV = "Height CV", Rare_Species_Ratio = "Rare Species Ratio",
             Intra_Species_SDiv = "Intra-species SDiv", Turnover = "Turnover")
FACTORS_BY_SCALE <- list("0.5m" = names(FACTORS)[1:3], "1m" = names(FACTORS), "1.5m" = names(FACTORS))
STRATA <- c("Low", "High")

# ---- 1-2. strata and statistics -------------------------------------------------------
stratify <- function(k, sc, f) {
  v <- PREDICTORS[[k]]$var
  d <- dat[dat$Scale == sc, ]
  d <- d[complete.cases(d[, c("Observed_SR", v, f)]), ]
  d <- d[order(d[[f]]), ]
  d$Stratum <- factor(ifelse(seq_len(nrow(d)) <= floor(nrow(d) / 2), "Low", "High"), levels = STRATA)
  data.frame(x = d[[v]], y = d$Observed_SR, Stratum = d$Stratum)
}

stats <- list(); strata_data <- list()
for (k in names(PREDICTORS)) for (sc in SCALES) for (f in FACTORS_BY_SCALE[[sc]]) {
  d <- stratify(k, sc, f)
  strata_data[[paste(k, sc, f)]] <- d
  for (s in STRATA) {
    ds <- d[d$Stratum == s, ]
    cf <- summary(lm(y ~ x, data = ds))
    stats[[length(stats) + 1]] <- data.frame(
      Panel = k, Predictor = PREDICTORS[[k]]$var, Scale = sc, Factor = f, Stratum = s, n = nrow(ds),
      adj_R2 = cf$adj.r.squared, slope = cf$coefficients[2, 1], slope_se = cf$coefficients[2, 2],
      P = cf$coefficients[2, 4])
  }
}
stats <- do.call(rbind, stats)
write.csv(stats, file.path(out_dir, "Fig12_statistics.csv"), row.names = FALSE)

# ---- 3. panels --------------------------------------------------------------------------
BASE_SIZE  <- 20.5
STRATA_COLS <- c(Low = "#78AB30", High = "#7D2E8F")
POINT_SIZE  <- 2.75
PANEL_W <- 4.65; PANEL_H <- 4.2      # single panels (inches)

stratum_panel <- function(k, sc, f, x_title = TRUE, y_title = TRUE) {
  ggplot(strata_data[[paste(k, sc, f)]], aes(x = x, y = y, color = Stratum)) +
    geom_point(alpha = 0.70, size = POINT_SIZE) +
    geom_smooth(aes(fill = Stratum), method = "lm", formula = y ~ x, se = TRUE, linewidth = 1.0, alpha = 0.23) +
    scale_color_manual(values = STRATA_COLS, guide = "none") +
    scale_fill_manual(values = STRATA_COLS, guide = "none") +
    scale_x_continuous(breaks = integer_breaks, expand = expansion(mult = c(0.07, 0.10))) +
    scale_y_continuous(breaks = integer_breaks, expand = expansion(mult = c(0.07, 0.28))) +
    labs(x = if (x_title) PREDICTORS[[k]]$x_title else NULL,
         y = if (y_title) "Observed Species Richness" else NULL) +
    theme_nature(BASE_SIZE) +
    theme(axis.text = element_text(size = BASE_SIZE - 2, face = "bold", colour = "black"),
          plot.title = element_blank(),
          plot.margin = ggplot2::margin(t = 25, r = 30, b = 8, l = 8))
}

for (k in names(PREDICTORS)) for (sc in SCALES) for (f in FACTORS_BY_SCALE[[sc]])
  save_tiff(stratum_panel(k, sc, f), file.path(out_dir, sprintf("Fig12_%s_%s_%s.tiff", k, sc, f)), PANEL_W, PANEL_H)

# adjusted R2 and slope +/- standard error of the two strata; the slopes start
# at the same position in both lines
label_text <- function(k, sc, f) {
  st <- stats[stats$Panel == k & stats$Scale == sc & stats$Factor == f, ]
  st <- st[match(STRATA, st$Stratum), ]
  r2 <- ifelse(st$adj_R2 < 0, "R\u00b2<0",
               paste0("R\u00b2=", sprintf("%.2f", st$adj_R2),
                      vapply(st$P, function(p) { s <- p_to_stars(p); if (s == "NS") "" else s }, character(1))))
  list(r2 = paste0(r2, ","), slope = paste0("Slope=", sprintf("%.2f", st$slope), "\u00b1", sprintf("%.2f", st$slope_se)))
}

# ---- assembled figures ------------------------------------------------------------------
# Positions and sizes below are in inches from the top-left corner; panels are
# listed by row (0.5 m: Height CV, Rare Species Ratio, Intra-species SDiv; 1 m and
# 1.5 m: the same and Turnover).
LAYOUT <- list(
  a = list(W = 16.9627, H = 13.1754, band = c(0.0603, 0.0496, 16.9367, 0.6951),
           letter = c(0.125, 0.489), header_base = 0.5142,
           x_title_x = c(2.5024, 6.5711, 10.6375, 14.7038), x_title_y = 13.0217,
           panel_x = c(0.7045, 4.7708, 8.8300, 0.7045, 4.7850, 8.8442, 12.8869, 0.7045, 4.8110, 8.8395, 12.9058),
           panel_y = c(1.3924, 1.3924, 1.3924, 5.3524, 5.3524, 5.3524, 5.3524, 9.3147, 9.3029, 9.3029, 9.3029),
           label_x = c(0.7861, 4.8501, 8.9117, 0.7861, 4.8643, 8.9235, 12.9662, 0.8074, 4.9021, 8.9235, 13.0229),
           label_y = c(1.6549, 1.6549, 1.6549, 5.6172, 5.6172, 5.6172, 5.6172, 9.5653, 9.5653, 9.5653, 9.5653),
           label_gap = c(rep(0.2577, 3), rep(0.2553, 4), rep(0.2979, 4)),
           title_x = c(2.4481, 6.5274, 10.5334, 14.5560), title_y = c(1.1868, 1.1868, 1.1868, 5.1444),
           strips = rbind(c(12.2770, 1.4019, 12.8586, 4.5250), c(16.3362, 5.3619, 16.9202, 8.4731),
                          c(16.3362, 9.2911, 16.9296, 12.4212)),
           legend_x = c(14.1707, 14.8350), legend_dx = c(0, 0.0118), legend_y = c(2.2424, 3.2684), legend_text_x = 14.9816,
           legend_text_y = c(2.3547, 3.3807)),
  b = list(W = 17.1826, H = 13.2061, band = c(0.0508, 0.0000, 17.1566, 0.6076),
           letter = c(0.157, 0.454), header_base = 0.4244,
           x_title_x = c(2.4445, 6.5629, 10.6824, 14.7995), x_title_y = 13.0808,
           panel_x = c(0.7281, 4.8323, 8.9695, 0.7281, 4.8488, 8.9837, 13.1375, 0.7281, 4.8394, 8.9766, 13.1210),
           panel_y = c(1.3333, 1.3333, 1.3333, 5.3287, 5.3287, 5.3287, 5.3287, 9.3218, 9.3218, 9.3241, 9.3241),
           label_x = c(0.8216, 4.9494, 9.0677, 0.8216, 4.9163, 9.0677, 13.3373, 0.8098, 4.9163, 9.1103, 13.1979),
           label_y = c(1.6053, 1.6053, 1.6053, 5.6219, 5.6148, 5.6219, 5.6148, 9.6173, 9.5866, 9.6173, 9.5890),
           label_gap = c(rep(0.2979, 7), 0.2955, 0.2553, 0.2955, 0.2553),
           title_x = c(2.4611, 6.5321, 10.6788, 14.8173), title_y = c(1.1608, 1.1655, 1.1655, 5.1704),
           strips = rbind(c(12.4141, 1.3263, 12.9933, 4.4659), c(16.5868, 5.3264, 17.1731, 8.4565),
                          c(16.5726, 9.3218, 17.1542, 12.4519)),
           legend_x = c(14.3503, 15.0123), legend_dx = c(0, 0.0118), legend_y = c(2.0320, 3.0580), legend_text_x = 15.1613,
           legend_text_y = c(2.1466, 3.1750)))
LETTER_SIZE <- 29.8; HEADER_SIZE <- 27.7; TITLE_SIZE <- 24; STRIP_SIZE <- 24
LABEL_SIZE  <- 18.54; LEGEND_SIZE <- 24.8
STRIP_LINE  <- "#7F7F7F"; STRIP_LWD <- 1.0                  # strip outline (pt)
LEGEND_LWD  <- 5.1; LEGEND_DOT <- 0.319                     # legend line (pt) and dot diameter (in)
BORDER_IN   <- 0.6 * .pt / 96                               # width of the panel border (in)

gp_text <- function(size, col = "black") grid::gpar(fontfamily = "Arial", fontface = "bold", fontsize = size, col = col)
in_w <- function(u) with_agg(grid::convertWidth(sum(u), "in", valueOnly = TRUE))
in_h <- function(u) with_agg(grid::convertHeight(sum(u), "in", valueOnly = TRUE))

g_ref <- with_agg(ggplotGrob(stratum_panel("a", "1m", "Height_CV")))
pcol  <- unique(g_ref$layout$l[g_ref$layout$name == "panel"])
prow  <- unique(g_ref$layout$t[g_ref$layout$name == "panel"])
AREA_W <- 3.4408; AREA_H <- 3.1128                            # panel area (in)

assemble <- function(k) {
  L <- LAYOUT[[k]]; H <- L$H
  box <- function(x0, y0, x1, y1, fill, col = NA, lwd = 1)
    grid::rectGrob(x = grid::unit(x0, "in"), y = grid::unit(H - y1, "in"), width = grid::unit(x1 - x0, "in"),
                   height = grid::unit(y1 - y0, "in"), just = c("left", "bottom"),
                   gp = grid::gpar(fill = fill, col = col, lwd = lwd * 96 / 72))
  text <- function(txt, x, y, size, hjust = 0.5, col = "black", rot = 0, vjust = 0)
    grid::textGrob(txt, x = grid::unit(x, "in"), y = grid::unit(H - y, "in"), hjust = hjust, vjust = vjust,
                   rot = rot, gp = gp_text(size, col))
  place <- function(p, x, y) {
    g <- with_agg(ggplotGrob(p + theme(plot.background = element_blank())))
    g$widths[pcol] <- grid::unit(AREA_W, "in"); g$heights[prow] <- grid::unit(AREA_H, "in")
    left <- in_w(g$widths[seq_len(pcol - 1)]); top <- in_h(g$heights[seq_len(prow - 1)])
    x <- x + BORDER_IN / 2; y <- y + BORDER_IN / 2
    grid::grobTree(g, vp = grid::viewport(x = grid::unit(x - left, "in"), y = grid::unit(H - y + top, "in"),
                                          width = sum(g$widths), height = sum(g$heights), just = c("left", "top")))
  }
  grobs <- list(box(0, 0, L$W, H, "white"),
                box(L$band[1], L$band[2], L$band[3], L$band[4], STRIP_FILL),
                text(k, L$letter[1], L$letter[2], LETTER_SIZE, hjust = 0),
                text(PREDICTORS[[k]]$header, (L$band[1] + L$band[3]) / 2, L$header_base, HEADER_SIZE))
  i <- 0
  for (r in seq_along(SCALES)) {
    sc <- SCALES[r]; fs <- FACTORS_BY_SCALE[[sc]]
    for (j in seq_along(fs)) {
      i <- i + 1
      p <- stratum_panel(k, sc, fs[j], x_title = FALSE, y_title = j == 1)
      grobs[[length(grobs) + 1]] <- place(p, L$panel_x[i], L$panel_y[i])
      lb <- label_text(k, sc, fs[j])
      wd <- function(s, size) in_w(grid::grobWidth(grid::textGrob(s, gp = gp_text(size))))
      w1 <- max(vapply(lb$r2, wd, numeric(1), size = LABEL_SIZE)); sp <- wd(" ", LABEL_SIZE)
      w2 <- max(vapply(lb$slope, wd, numeric(1), size = LABEL_SIZE))
      # labels wider than the panel are reduced to fit
      room <- L$panel_x[i] + AREA_W + BORDER_IN + 0.015 - L$label_x[i]
      size <- if (w1 + sp + w2 > room) LABEL_SIZE * room / (w1 + sp + w2) else LABEL_SIZE
      w1 <- w1 * size / LABEL_SIZE; sp <- sp * size / LABEL_SIZE
      for (s in 1:2) {
        yy <- L$label_y[i] + (s - 1) * L$label_gap[i]
        grobs <- c(grobs, list(text(lb$r2[s], L$label_x[i], yy, size, hjust = 0, col = STRATA_COLS[[s]]),
                               text(lb$slope[s], L$label_x[i] + w1 + sp, yy, size, hjust = 0, col = STRATA_COLS[[s]])))
      }
    }
    # scale strip to the right of the last panel of the row
    st <- L$strips[r, ]; h <- STRIP_LWD / 72 / 2
    grobs <- c(grobs, list(box(st[1] + h, st[2] + h, st[3] - h, st[4] - h, STRIP_FILL, STRIP_LINE, STRIP_LWD),
                           text(SCALE_TITLES[[sc]], (st[1] + st[3]) / 2, (st[2] + st[4]) / 2, STRIP_SIZE,
                                rot = -90, vjust = 0.5)))
  }
  # x-axis titles below the 1.5 m panels
  for (j in seq_along(FACTORS))
    grobs[[length(grobs) + 1]] <- text(PREDICTORS[[k]]$x_title, L$x_title_x[j], L$x_title_y, BASE_SIZE)
  # factor titles above the 0.5 m panels and above the first Turnover panel
  for (j in seq_along(FACTORS))
    grobs[[length(grobs) + 1]] <- text(FACTORS[[j]], L$title_x[j], L$title_y[j], TITLE_SIZE)
  # legend of the two strata
  for (s in 1:2) {
    yy <- grid::unit(H - L$legend_y[s], "in"); xx <- L$legend_x + L$legend_dx[s]
    grobs <- c(grobs, list(
      grid::segmentsGrob(grid::unit(xx[1], "in"), yy, grid::unit(xx[2], "in"), yy,
                         gp = grid::gpar(col = STRATA_COLS[[s]], lwd = LEGEND_LWD * 96 / 72, lineend = "butt")),
      grid::circleGrob(grid::unit(mean(xx), "in"), yy, r = grid::unit(LEGEND_DOT / 2, "in"),
                       gp = grid::gpar(fill = STRATA_COLS[[s]], col = NA)),
      text(STRATA[s], L$legend_text_x, L$legend_text_y[s], LEGEND_SIZE, hjust = 0)))
  }
  ragg::agg_tiff(file.path(out_dir, sprintf("Fig12%s.tiff", k)), width = L$W, height = H, units = "in",
                 res = 600, compression = "lzw")
  grid::grid.draw(grid::gTree(children = do.call(grid::gList, grobs)))
  invisible(dev.off())
}
for (k in names(PREDICTORS)) assemble(k)

cat("Fig. 12 written to", out_dir, "\n")
