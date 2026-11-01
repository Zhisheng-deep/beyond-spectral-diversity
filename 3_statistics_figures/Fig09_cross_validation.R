# =============================================================================
# Fig. 9  Robustness of the four remotely sensed predictors of observed
#         full-community species richness (SR) at three spatial scales for
#         grazing, fencing, and all plots.
#         a, SR predicted by deep learning (DL_SR)
#         b, spectral species richness (SSR_km_z)
#         c, dominance metrics from the deep learning species maps (DL_Top)
#         d, dominance metrics from the spectral clusters (SC_Top, k-means on
#            z-scored reflectance)
#
# Steps
#   1. All observations. For each predictor, scale and subset, fit observed SR
#      on the predictor with ordinary least squares (adjusted R2). DL_Top and
#      SC_Top enter as the five dominance metrics (Top1, Top2, Top12, Top3,
#      Top123) screened for collinearity (pairwise Spearman |rho| <= 0.8,
#      VIF <= 5). Leave-one-unit-out cross-validation: each management unit
#      is held out in turn, the model (including the screening) is refitted on
#      the remaining units and predicts the held-out unit; cross-validated
#      R2 = 1 - sum((obs - pred)^2) / sum((obs - mean(obs))^2).
#   2. One observation per plot. The same analysis on data sets with one
#      observation per 1.5 m plot: each of the nine 0.5 m quadrats (Window 1-9),
#      each of the four 1 m subplots (Window 1-4), and the 1.5 m plot itself.
#      Adjusted R2 uses the predictors retained in step 1; cross-validation
#      repeats the screening within each training set. Range of adjusted R2 and
#      median cross-validated R2 over these data sets.
#   3. Plot the four series against scale; values below zero are drawn as open
#      triangles on the baseline and the range is cut at zero. Export the
#      single panels and the assembled figure.
#
# Input  : analysis_table.csv (columns Scale, Management, Unit, Window,
#          Observed_SR, DL_SR, SSR_km_z, DL_Top*, SC_Top*)
# Output : output/Fig09_statistics.csv        retained predictors, adjusted and
#                                             cross-validated R2 of every series
#          output/Fig09_<predictor>_<subset>.tiff   single panels (4.65 x 4.2 in)
#          output/Fig09.tiff                  assembled figure
#
# Usage  : Rscript Fig09_cross_validation.R [data_file] [output_dir]
# Default: data_file  = "../data/simulated/analysis_table.csv"
#          output_dir = "output"
# =============================================================================

source("functions.R")
suppressPackageStartupMessages(library(tidyr))

args      <- commandArgs(trailingOnly = TRUE)
data_file <- if (length(args) >= 1) args[1] else file.path(DATA_DIR, "analysis_table.csv")
out_dir   <- if (length(args) >= 2) args[2] else "output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

RESP    <- "Observed_SR"
SUBSETS <- c("Overall", "Grazing", "Fencing")
PREDS   <- c("DL_SR", "SSR", "DL_Top", "SC_Top")

dat <- as.data.frame(read_table(data_file))
dat$Unit <- as.character(dat$Unit)

# ---- model of each predictor -----------------------------------------------------
# Returns the OLS model fitted on d and the retained predictors
fit_predictor <- function(d, pred) {
  if (pred %in% c("DL_SR", "SSR")) {
    v <- if (pred == "DL_SR") "DL_SR" else "SSR_km_z"
    return(list(model = lm(as.formula(paste(RESP, "~", v)), data = d), vars = v))
  }
  cand <- top_vars(if (pred == "DL_Top") "DL" else "SC")
  # metrics constant in d cannot enter the regression
  cand <- cand[vapply(cand, function(v) isTRUE(stats::sd(d[[v]], na.rm = TRUE) > 0), logical(1))]
  s <- screen_fit(d, RESP, cand)
  list(model = s$model, vars = s$final_vars)
}

# Cross-validated R2 of out-of-unit predictions (NA when fewer than three)
r2_cv <- function(obs, pred) {
  ok <- !is.na(pred)
  if (sum(ok) < 3) return(NA_real_)
  1 - sum((obs[ok] - pred[ok])^2) / sum((obs[ok] - mean(obs[ok]))^2)
}

# Leave-one-unit-out cross-validation with the model refitted on each training set
louo_r2 <- function(d, pred) {
  p <- rep(NA_real_, nrow(d))
  for (u in unique(d$Unit)) {
    test <- d$Unit == u
    f <- fit_predictor(d[!test, ], pred)
    if (!is.null(f$model)) p[test] <- predict(f$model, newdata = d[test, ])
  }
  r2_cv(d[[RESP]], p)
}

# ---- 1-2. statistics -------------------------------------------------------------
res <- list()
for (sc in SCALES) for (mg in SUBSETS) {
  d_all <- dat[dat$Scale == sc & (mg == "Overall" | dat$Management == mg), ]
  # data sets with one observation per plot
  one_obs <- if (sc == "1.5m") list(d_all) else
    lapply(sort(unique(d_all$Window)), function(w) d_all[d_all$Window == w, ])
  for (pred in PREDS) {
    full  <- fit_predictor(d_all, pred)
    form  <- if (is.null(full$model)) NULL else formula(full$model)
    adj_one <- sapply(one_obs, function(d)
      if (is.null(form)) NA_real_ else summary(lm(form, data = d))$adj.r.squared)
    cv_one  <- sapply(one_obs, louo_r2, pred = pred)
    res[[length(res) + 1]] <- data.frame(
      Predictor = pred, Scale = sc, Management = mg,
      n = nrow(d_all), n_units = length(unique(d_all$Unit)),
      retained = paste(full$vars, collapse = ";"),
      adjR2_all = if (is.null(form)) NA_real_ else summary(full$model)$adj.r.squared,
      cvR2_all  = louo_r2(d_all, pred),
      n_one_obs_sets = length(one_obs), n_one_obs = nrow(one_obs[[1]]),
      adjR2_one_min = min(adj_one, na.rm = TRUE), adjR2_one_max = max(adj_one, na.rm = TRUE),
      cvR2_one_median = median(cv_one, na.rm = TRUE))
  }
}
res <- do.call(rbind, res)
write.csv(res, file.path(out_dir, "Fig09_statistics.csv"), row.names = FALSE)

# ---- 3. figure -----------------------------------------------------------------------
BASE_SIZE    <- 20.5
PRED_TITLES  <- c(DL_SR = "Deep Learning-SR", SSR = "Spectral Clustering-SSR",
                  DL_Top = "Deep Learning-Top", SC_Top = "Spectral Clustering-Top")
SCALE_LABELS <- c("0.5m" = "0.5m\u00d70.5m", "1m" = "1m\u00d71m", "1.5m" = "1.5m\u00d71.5m")
MGMT3_COLS   <- c(Grazing = "#E69F00", Fencing = "#56B4E9", Overall = "black")
RIBBON_COL   <- "#B8B8B8"
RIBBON_ALPHA <- 0.45
SERIES       <- c("adj_all", "cv_all", "cv_one")
SERIES_SHAPE <- c(adj_all = 16, cv_all = 21, cv_one = 23)
SERIES_LTY   <- c(adj_all = "solid", cv_all = "solid", cv_one = "22")
SERIES_LWD   <- c(adj_all = 2.0, cv_all = 2.0, cv_one = 1.2)
pd <- position_dodge(width = 0.20)       # lines and points share the same dodge

d_plot <- res
d_plot$Scale <- factor(SCALE_LABELS[d_plot$Scale], levels = SCALE_LABELS)
ribbon_df <- transform(d_plot, lo_plot = pmax(adjR2_one_min, 0), hi_plot = pmax(adjR2_one_max, 0))
pts <- d_plot %>%
  select(Predictor, Scale, Management, adj_all = adjR2_all, cv_all = cvR2_all, cv_one = cvR2_one_median) %>%
  pivot_longer(all_of(SERIES), names_to = "Series", values_to = "value") %>%
  mutate(Series = factor(Series, levels = SERIES),
         below_zero = value < 0,          # drawn as an open triangle on the baseline
         y_plot = pmax(value, 0))

# y axis shared by all panels: upper limit = largest value (rounded to two
# decimals) rounded up to 0.1; breaks every 0.2 above 0.6, otherwise every 0.1
y_max   <- ceiling(round(max(c(pts$y_plot, ribbon_df$hi_plot), na.rm = TRUE), 2) * 10) / 10
Y_AXIS  <- list(limits = c(0, y_max), breaks = seq(0, y_max, by = if (y_max > 0.6) 0.2 else 0.1))

# layers: range -> zero line -> lines -> points -> triangles
r2_panel <- function(pred, mg) {
  p_sub <- pts[pts$Predictor == pred & pts$Management == mg, ]
  r_sub <- ribbon_df[ribbon_df$Predictor == pred & ribbon_df$Management == mg, ]
  ggplot(p_sub, aes(x = Scale)) +
    geom_ribbon(data = r_sub, aes(x = Scale, ymin = lo_plot, ymax = hi_plot, group = Management),
                fill = RIBBON_COL, alpha = RIBBON_ALPHA, colour = NA, inherit.aes = FALSE) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.5) +
    geom_line(aes(y = y_plot, colour = Management, group = Series,
                  linetype = Series, linewidth = Series), position = pd) +
    geom_point(data = ~ subset(.x, !below_zero),
               aes(y = y_plot, colour = Management, group = Series, shape = Series),
               fill = "white", size = 5, stroke = 1.5, position = pd) +
    geom_point(data = ~ subset(.x, below_zero),
               aes(y = y_plot, colour = Management, group = Series),
               shape = 25, fill = "white", size = 5, stroke = 1.5, position = pd) +
    scale_colour_manual(values = MGMT3_COLS) +
    scale_shape_manual(values = SERIES_SHAPE) +
    scale_linetype_manual(values = SERIES_LTY) +
    scale_linewidth_manual(values = SERIES_LWD) +
    scale_y_continuous(breaks = Y_AXIS$breaks, expand = expansion(mult = c(0.05, 0.05))) +
    coord_cartesian(ylim = Y_AXIS$limits) +
    labs(x = "Scale", y = "R\u00b2") +
    theme_nature(BASE_SIZE) +
    theme(legend.position = "none")
}

PANEL_W <- 4.65; PANEL_H <- 4.2     # single panel (inches)
panels <- list()
for (pred in PREDS) for (mg in SUBSETS[c(2, 3, 1)]) {
  p <- r2_panel(pred, mg)
  panels[[paste(pred, mg)]] <- p
  save_tiff(p, file.path(out_dir, sprintf("Fig09_%s_%s.tiff", pred, mg)), PANEL_W, PANEL_H)
}

# ---- assembled figure ----------------------------------------------------------------
# As in the published figure, the single panels are placed on a grid without
# rescaling: y-axis labels only in the left column, x-axis labels only in the
# bottom row, tick marks on every panel. Positions and sizes below are in inches
# from the top-left corner.
FIG_W <- 18.21; FIG_H <- 13.96
COL_X <- 0.8436 + (0:3) * 4.2533        # left edge of the panel area of each column
ROW_Y <- 0.7673 + (0:2) * 3.8984        # top edge of the panel area of each row
STRIP_LINE  <- "#7F7F7F"
STRIP_LWD   <- 0.95                     # outline of the panel headers (pt)
HEADER_SIZE <- 22.7                     # panel headers and legend text (pt)
LETTER_SIZE <- 29.8

# panel area of a single panel at its export size
g_full  <- with_agg(ggplotGrob(panels[[1]]))
pcol    <- unique(g_full$layout$l[g_full$layout$name == "panel"])
prow    <- unique(g_full$layout$t[g_full$layout$name == "panel"])
in_w    <- function(u) with_agg(grid::convertWidth(sum(u), "in", valueOnly = TRUE))
in_h    <- function(u) with_agg(grid::convertHeight(sum(u), "in", valueOnly = TRUE))
AREA_W  <- PANEL_W - in_w(g_full$widths[-pcol])
AREA_H  <- PANEL_H - in_h(g_full$heights[-prow])

# grob of a panel with its panel area fixed, placed with that area at (x, y)
place_panel <- function(p, x, y) {
  g <- with_agg(ggplotGrob(p + theme(plot.background = element_blank())))
  g$widths[pcol] <- grid::unit(AREA_W, "in"); g$heights[prow] <- grid::unit(AREA_H, "in")
  left <- in_w(g$widths[seq_len(pcol - 1)]); top <- in_h(g$heights[seq_len(prow - 1)])
  grid::grobTree(g, vp = grid::viewport(x = grid::unit(x - left, "in"), y = grid::unit(FIG_H - y + top, "in"),
                                        width = sum(g$widths), height = sum(g$heights),
                                        just = c("left", "top")))
}
gp_text <- function(size, col = "black") grid::gpar(fontfamily = "Arial", fontface = "bold", fontsize = size, col = col)
box <- function(x0, y0, x1, y1, fill, col = NA, lwd = 1)
  grid::rectGrob(x = grid::unit(x0, "in"), y = grid::unit(FIG_H - y1, "in"),
                 width = grid::unit(x1 - x0, "in"), height = grid::unit(y1 - y0, "in"),
                 just = c("left", "bottom"), gp = grid::gpar(fill = fill, col = col, lwd = lwd * 96 / 72))
label <- function(text, x, y, size, hjust = 0.5, vjust = 0.5, rot = 0)
  grid::textGrob(text, x = grid::unit(x, "in"), y = grid::unit(FIG_H - y, "in"), hjust = hjust,
                 vjust = vjust, rot = rot, gp = gp_text(size))

grobs <- list(box(0, 0, FIG_W, FIG_H, "white"))
for (i in seq_along(PREDS)) for (j in 1:3) {
  mg <- SUBSETS[c(2, 3, 1)][j]
  p  <- panels[[paste(PREDS[i], mg)]]
  if (i > 1) p <- p + theme(axis.text.y = element_blank(), axis.title.y = element_blank())
  if (j < 3) p <- p + theme(axis.text.x = element_blank(), axis.title.x = element_blank())
  grobs[[length(grobs) + 1]] <- place_panel(p, COL_X[i], ROW_Y[j])
}
# column headers and panel letters
for (i in seq_along(PREDS)) {
  x0 <- COL_X[i] - 0.0090; x1 <- COL_X[i] + AREA_W - 0.0135
  grobs <- c(grobs, list(box(x0, 0.0880, x1, 0.7518, STRIP_FILL, STRIP_LINE, STRIP_LWD),
                         label(PRED_TITLES[PREDS[i]], (x0 + x1) / 2, 0.5389, HEADER_SIZE, vjust = 0),
                         label(letters[i], c(0.3868, 4.7541, 9.0312, 13.2797)[i], 0.3319,
                               LETTER_SIZE, hjust = 0, vjust = 0)))
}
# row headers
for (j in 1:3) {
  y0 <- ROW_Y[j] - 0.0100; y1 <- ROW_Y[j] + AREA_H - 0.0140
  grobs <- c(grobs, list(box(17.4794, y0, 18.1416, y1, STRIP_FILL, STRIP_LINE, STRIP_LWD),
                         label(SUBSETS[c(2, 3, 1)][j], 17.6940, (y0 + y1) / 2, HEADER_SIZE, vjust = 0, rot = -90)))
}
# legend as in the published figure: line with filled circle, line with open
# circle, open diamond, shaded range
key_pts <- function(x, y, pch, size_pt, lwd_pt, fill = "black")
  grid::pointsGrob(grid::unit(x, "in"), grid::unit(FIG_H - y, "in"), pch = pch,
                   gp = grid::gpar(col = "black", fill = fill, fontsize = size_pt, lwd = lwd_pt * 96 / 72))
key_line <- function(x, y, w, lwd_pt)
  grid::segmentsGrob(grid::unit(x - w / 2, "in"), grid::unit(FIG_H - y, "in"),
                     grid::unit(x + w / 2, "in"), grid::unit(FIG_H - y, "in"),
                     gp = grid::gpar(col = "black", lwd = lwd_pt * 96 / 72, lineend = "butt"))
# legend text with a superscript 2 (two thirds of the text size, raised by 0.3 em)
sup_label <- function(pre, post, x, y, size) {
  x0 <- grid::unit(x, "in"); y0 <- grid::unit(FIG_H - y, "in")
  t1 <- grid::textGrob(pre, x0, y0, hjust = 0, vjust = 0, gp = gp_text(size))
  x2 <- x0 + grid::grobWidth(t1)
  t2 <- grid::textGrob("2", x2, y0 + grid::unit(0.3 * size, "pt"), hjust = 0, vjust = 0, gp = gp_text(size * 0.67))
  t3 <- grid::textGrob(post, x2 + grid::grobWidth(t2), y0, hjust = 0, vjust = 0, gp = gp_text(size))
  grid::gTree(children = grid::gList(t1, t2, t3))
}
LEG <- list(
  list(pre = "All Observations: Adjusted R", post = "", tx = 4.7304, ty = 13.2227),
  list(pre = "All Observations: Cross-validated R", post = "", tx = 9.8831, ty = 13.2227),
  list(pre = "One Observation Per Plot: Cross-validated R", post = " (Median)", tx = 1.3648, ty = 13.8484),
  list(pre = "One Observation Per Plot: Adjusted R", post = " (Range)", tx = 10.3646, ty = 13.8532))
for (k in LEG) grobs[[length(grobs) + 1]] <- sup_label(k$pre, k$post, k$tx, k$ty, HEADER_SIZE)
grobs <- c(grobs, list(
  key_line(4.3257, 13.0930, 0.4296, 7.54), key_pts(4.3257, 13.0930, 16, 32.2, 0),
  key_line(9.5189, 13.1013, 0.4272, 6.34), key_pts(9.5189, 13.1013, 21, 31.8, 3.94, fill = "white"),
  key_pts(0.9696, 13.7354, 23, 33.2, 3.51, fill = "white"),
  box(9.8595, 13.5474, 10.2345, 13.9186, grDevices::adjustcolor(RIBBON_COL, RIBBON_ALPHA))))

fig <- grid::gTree(children = do.call(grid::gList, grobs))
ragg::agg_tiff(file.path(out_dir, "Fig09.tiff"), width = FIG_W, height = FIG_H, units = "in",
               res = 600, compression = "lzw")
grid::grid.draw(fig)
invisible(dev.off())

cat("Fig. 9 written to", out_dir, "\n")
