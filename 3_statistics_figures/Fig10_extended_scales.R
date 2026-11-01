# =============================================================================
# Fig. 10  Adjusted R2 of ordinary least squares regressions of observed
#          full-community species richness (SR) on each of four remotely sensed
#          predictors across seven spatial scales for grazing, fencing, and all
#          plots, with 95% cluster-bootstrap intervals.
#          a, SR predicted by deep learning (DL_SR)
#          b, Top123 from the deep learning species maps (DL_Top123)
#          c, spectral species richness (SSR_km_z, k-means on z-scored
#             reflectance)
#          d, Top123 from the spectral clusters (SC_Top123)
#
# Steps
#   1. For each predictor, scale and subset (grazing, fencing, all plots), fit
#      observed SR on the predictor with ordinary least squares (adjusted R2).
#      Scales: the three scales of the sampling design (0.5, 1 and 1.5 m), the
#      two intermediate windows (0.5 x 1 and 1 x 1.5 m) and the pooled plots
#      formed from all combinations of two and three plots within a management
#      unit. SC has no pooled plots.
#   2. Cluster bootstrap (2000 replicates): plots are resampled within
#      management units with all their quadrats, subplots and windows, and the
#      pooled plots are reconstructed from the resampled plots; coefficients are
#      refitted. 95% percentile intervals of adjusted R2.
#   3. Plot adjusted R2 against scale with the intervals as shaded bands;
#      filled points for the three scales of the sampling design, open points
#      for the four additional scales; adjusted R2 below zero is drawn as an
#      open triangle on the baseline and interval bounds below zero are cut at
#      the baseline. Export the single panels and the assembled figure.
#
# Input  : analysis_table.csv   (columns Scale, Management, Unit, Plot15,
#                                Observed_SR, DL_SR, DL_Top123, SSR_km_z,
#                                SC_Top123; scales 0.5m, 1m, 1.5m)
#          multiscale_table.csv (the same columns and member_plots; scales
#                                0.5x1m, 1x1.5m, pool2, pool3)
# Output : output/Fig10_statistics.csv     adjusted R2 and bootstrap interval of
#                                          every model
#          output/Fig10_<predictor>.tiff   single panels (7.6 x 5.6 in)
#          output/Fig10.tiff               assembled figure
#
# Usage  : Rscript Fig10_extended_scales.R [data_file] [multiscale_file] [output_dir]
# Default: data_file       = "../data/simulated/analysis_table.csv"
#          multiscale_file = "../data/simulated/multiscale_table.csv"
#          output_dir      = "output"
# =============================================================================

source("functions.R")

args       <- commandArgs(trailingOnly = TRUE)
data_file  <- if (length(args) >= 1) args[1] else file.path(DATA_DIR, "analysis_table.csv")
multi_file <- if (length(args) >= 2) args[2] else file.path(DATA_DIR, "multiscale_table.csv")
out_dir    <- if (length(args) >= 3) args[3] else "output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

B    <- 2000
SEED <- 20260908
SUBSETS    <- c("Overall", "Grazing", "Fencing")
PREDICTORS <- c(DL_SR = "DL_SR", DL_Top123 = "DL_Top123", SSR = "SSR_km_z", SC_Top123 = "SC_Top123")
EXTENTS    <- c("0.5m", "0.5x1m", "1m", "1x1.5m", "1.5m", "pool2", "pool3")
DESIGNED   <- c("0.5m", "1m", "1.5m")

# ---- data: the three design scales and the additional scales in one table ------------
cols <- c("Scale", "Management", "Unit", "Plot15", "member_plots", "Observed_SR", PREDICTORS)
main <- read_table(data_file)
main$member_plots <- main$Plot15                    # each row belongs to one plot
multi <- read_table(multi_file)
for (v in setdiff(cols, names(multi))) multi[[v]] <- NA
dat <- rbind(main[, cols], multi[multi$Scale %in% EXTENTS, cols])
dat$Unit <- as.character(dat$Unit)
dat$Plot15 <- as.character(dat$Plot15); dat$member_plots <- as.character(dat$member_plots)

# ---- 1. full-data models ---------------------------------------------------------------
models <- expand.grid(Scale = EXTENTS, Management = SUBSETS, Predictor = names(PREDICTORS),
                      stringsAsFactors = FALSE)
specs <- vector("list", nrow(models))
point <- vector("list", nrow(models))
for (i in seq_len(nrow(models))) {
  sc <- models$Scale[i]; mg <- models$Management[i]; v <- PREDICTORS[[models$Predictor[i]]]
  rows <- which(dat$Scale == sc & (mg == "Overall" | dat$Management == mg))
  d <- dat[rows, ]
  if (all(is.na(d[[v]]))) next                       # SC at the pooled plots
  fit <- lm(as.formula(paste("Observed_SR ~", v)), data = d)
  specs[[i]] <- list(rows = rows, X = model.matrix(formula(fit)[-2], data = d), y = d$Observed_SR)
  point[[i]] <- data.frame(Predictor = models$Predictor[i], Scale = sc, Management = mg,
                           n = nrow(d), n_units = length(unique(d$Unit)),
                           adjR2 = summary(fit)$adj.r.squared)
}

# ---- 2. cluster bootstrap --------------------------------------------------------------
boot <- cluster_bootstrap(dat, specs, B = B, seed = SEED)
res <- do.call(rbind, lapply(seq_len(nrow(models)), function(i) {
  if (is.null(point[[i]])) return(NULL)
  q <- boot_ci(boot[, i])
  cbind(point[[i]], boot_lower = q[["lower"]], boot_upper = q[["upper"]], n_boot_valid = q[["n_valid"]])
}))
write.csv(res, file.path(out_dir, "Fig10_statistics.csv"), row.names = FALSE)

# ---- 3. figure ---------------------------------------------------------------------------
BASE_SIZE    <- 25
MGMT3_COLS   <- c(Fencing = "#56B4E9", Grazing = "#E69F00", Overall = "black")      # lines and points
FILL3_COLS   <- c(Fencing = "#A9DBF5", Grazing = "#F5D48F", Overall = "#B8B8B8")    # bands
RIBBON_ALPHA <- 0.45
SCALE_LABELS <- c("0.5", "0.5\u00d71", "1", "1\u00d71.5", "1.5", "2\nplots", "3\nplots")
PRED_TITLES  <- c(DL_SR = "Deep Learning-SR", DL_Top123 = "Deep Learning-Top123",
                  SSR = "Spectral Clustering-SSR", SC_Top123 = "Spectral Clustering-Top123")
DODGE        <- 0.10                      # x offset per regime: Fencing -1, Grazing 0, Overall +1
MGMT_SHIFT   <- c(Fencing = -1, Grazing = 0, Overall = 1)
FILLED_SIZE  <- 5.85                      # filled points (design scales)

d_plot <- res
d_plot$Management_label <- factor(d_plot$Management, levels = c("Fencing", "Grazing", "Overall"))
d_plot$xnum  <- match(d_plot$Scale, EXTENTS)
d_plot$x_d   <- d_plot$xnum + MGMT_SHIFT[as.character(d_plot$Management_label)] * DODGE
d_plot$Layer <- ifelse(d_plot$Scale %in% DESIGNED, "Designed scale", "Additional scale")
d_plot$below_zero <- d_plot$adjR2 < 0
d_plot$y_plot  <- pmax(d_plot$adjR2, 0)           # negative adjusted R2 drawn on the baseline
d_plot$lo_plot <- pmax(d_plot$boot_lower, 0)      # interval bounds below zero cut at the baseline
d_plot$hi_plot <- d_plot$boot_upper

# y range per panel: upper limit = largest interval bound (rounded to two
# decimals) rounded up to 0.1; breaks every 0.2 above 0.6, otherwise every 0.1
y_axis_for <- function(df) {
  y_max <- ceiling(round(max(df$hi_plot, na.rm = TRUE), 2) * 10) / 10
  step  <- if (y_max > 0.6) 0.2 else 0.1
  list(limits = c(0, y_max), breaks = seq(0, y_max, by = step))
}

# layers: band -> zero line -> line separating single and pooled plots -> line
# -> points (filled: design scales, open: additional scales) -> triangles
r2_panel <- function(df) {
  ya <- y_axis_for(df)
  ggplot(df, aes(x = xnum, y = y_plot, group = Management_label, color = Management_label)) +
    geom_ribbon(aes(x = xnum, ymin = lo_plot, ymax = hi_plot, fill = Management_label),
                alpha = RIBBON_ALPHA, colour = NA, show.legend = FALSE) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.5) +
    geom_vline(xintercept = 5.5, linetype = "dotted", colour = "grey50", linewidth = 0.8) +
    geom_line(aes(x = x_d), linewidth = 2.0) +
    geom_point(data = ~ subset(.x, !below_zero & Layer == "Designed scale"), aes(x = x_d), size = FILLED_SIZE) +
    geom_point(data = ~ subset(.x, !below_zero & Layer == "Additional scale"), aes(x = x_d),
               shape = 21, fill = "white", size = 5, stroke = 1.5) +
    geom_point(data = ~ subset(.x, below_zero), aes(x = x_d), shape = 25, fill = "white",
               size = 5, stroke = 1.5, show.legend = FALSE) +
    scale_color_manual(values = MGMT3_COLS, breaks = c("Fencing", "Grazing", "Overall"),
                       name = "Management") +
    scale_fill_manual(values = FILL3_COLS, guide = "none") +
    scale_x_continuous(breaks = seq_along(EXTENTS), labels = SCALE_LABELS,
                       expand = expansion(add = 0.3)) +
    scale_y_continuous(breaks = ya$breaks, expand = expansion(mult = c(0.05, 0.05))) +
    coord_cartesian(xlim = c(1, 7), ylim = ya$limits) +
    labs(x = "Scale (m or pooled plots)", y = "R\u00b2") +
    theme_nature(BASE_SIZE) +
    theme(axis.text.x = element_text(lineheight = 0.8),
          plot.margin = ggplot2::margin(1, 12, 1, 1),
          legend.position = "none")
}

PANEL_W <- 7.6; PANEL_H <- 5.6     # single panel (inches)
panels <- list()
for (k in names(PREDICTORS)) {
  panels[[k]] <- r2_panel(d_plot[d_plot$Predictor == k, ])
  save_tiff(panels[[k]], file.path(out_dir, sprintf("Fig10_%s.tiff", k)), PANEL_W, PANEL_H)
}

# ---- assembled figure ----------------------------------------------------------------
# As in the published figure, the single panels are placed on a 2 x 2 grid
# without rescaling: x-axis labels only in the bottom row, y-axis title only in
# the left column, tick marks on every panel. Positions and sizes below are in
# inches from the top-left corner.
FIG_W <- 15.12; FIG_H <- 12.87
COL_X <- c(1.0037, 8.5114)           # left edge of the panel area of each column
ROW_Y <- c(0.7882, 6.4637)           # top edge of the panel area of each row
HEAD_X <- list(c(0.9945, 7.4947), c(8.5046, 15.0049))   # panel headers, left and right
HEAD_Y <- list(c(0.1177, 0.7810), c(5.7906, 6.4539))    # panel headers, top and bottom
HEAD_BASE <- c(0.5932, 6.2661)       # baseline of the header text
LETTER_X  <- list(c(0.0877, 7.6334), c(0.0853, 7.6357))
LETTER_Y  <- c(0.4113, 6.0985)
STRIP_LINE  <- "#7F7F7F"
STRIP_LWD   <- 0.95                  # outline of the panel headers (pt)
HEADER_SIZE <- 28.0                  # header text (pt)
LETTER_SIZE <- 29.8

# panel area of a single panel at its export size
g_full <- with_agg(ggplotGrob(panels[[1]]))
pcol   <- unique(g_full$layout$l[g_full$layout$name == "panel"])
prow   <- unique(g_full$layout$t[g_full$layout$name == "panel"])
in_w   <- function(u) with_agg(grid::convertWidth(sum(u), "in", valueOnly = TRUE))
in_h   <- function(u) with_agg(grid::convertHeight(sum(u), "in", valueOnly = TRUE))
AREA_W <- PANEL_W - in_w(g_full$widths[-pcol])
AREA_H <- PANEL_H - in_h(g_full$heights[-prow])

# grob of a panel with its panel area fixed, placed with that area at (x, y)
place_panel <- function(p, x, y) {
  g <- with_agg(ggplotGrob(p + theme(plot.background = element_blank())))
  g$widths[pcol] <- grid::unit(AREA_W, "in"); g$heights[prow] <- grid::unit(AREA_H, "in")
  left <- in_w(g$widths[seq_len(pcol - 1)]); top <- in_h(g$heights[seq_len(prow - 1)])
  grid::grobTree(g, vp = grid::viewport(x = grid::unit(x - left, "in"), y = grid::unit(FIG_H - y + top, "in"),
                                        width = sum(g$widths), height = sum(g$heights),
                                        just = c("left", "top")))
}
gp_text <- function(size) grid::gpar(fontfamily = "Arial", fontface = "bold", fontsize = size)
box <- function(x0, y0, x1, y1, fill, col = NA, lwd = 1)
  grid::rectGrob(x = grid::unit(x0, "in"), y = grid::unit(FIG_H - y1, "in"),
                 width = grid::unit(x1 - x0, "in"), height = grid::unit(y1 - y0, "in"),
                 just = c("left", "bottom"), gp = grid::gpar(fill = fill, col = col, lwd = lwd * 96 / 72))
label <- function(text, x, y, size, hjust = 0.5, vjust = 0)
  grid::textGrob(text, x = grid::unit(x, "in"), y = grid::unit(FIG_H - y, "in"),
                 hjust = hjust, vjust = vjust, gp = gp_text(size))

grobs <- list(box(0, 0, FIG_W, FIG_H, "white"))
layout <- matrix(names(PREDICTORS), nrow = 2, byrow = TRUE)     # a b / c d
for (r in 1:2) for (j in 1:2) {
  k <- layout[r, j]
  p <- panels[[k]]
  if (j == 2) p <- p + theme(axis.title.y = element_blank())
  if (r == 1) p <- p + theme(axis.text.x = element_blank(), axis.title.x = element_blank())
  hx <- HEAD_X[[j]]; hy <- HEAD_Y[[r]]
  grobs <- c(grobs, list(place_panel(p, COL_X[j], ROW_Y[r]),
                         box(hx[1], hy[1], hx[2], hy[2], STRIP_FILL, STRIP_LINE, STRIP_LWD),
                         label(PRED_TITLES[[k]], mean(hx), HEAD_BASE[r], HEADER_SIZE),
                         label(letters[(r - 1) * 2 + j], LETTER_X[[r]][j], LETTER_Y[r], LETTER_SIZE, hjust = 0)))
}

# legend: management regimes, keys with line and point
legend_plot <- panels[[1]] +
  guides(color = guide_legend(override.aes = list(shape = 16, fill = NA, stroke = 0, size = 5))) +
  theme(legend.position = "bottom", legend.background = element_blank(),
        legend.key.width = unit(0.237, "in"))   # line visible on both sides of the point
legend <- get_legend_grob(legend_plot)
LEGEND_X <- 7.884; LEGEND_Y <- 12.612   # centre of the legend
grobs[[length(grobs) + 1]] <- grid::grobTree(legend, vp = grid::viewport(
  x = grid::unit(LEGEND_X, "in"), y = grid::unit(FIG_H - LEGEND_Y, "in"),
  width = sum(legend$widths), height = sum(legend$heights)))

fig <- grid::gTree(children = do.call(grid::gList, grobs))
ragg::agg_tiff(file.path(out_dir, "Fig10.tiff"), width = FIG_W, height = FIG_H, units = "in",
               res = 600, compression = "lzw")
grid::grid.draw(fig)
invisible(dev.off())

cat("Fig. 10 written to", out_dir, "\n")
