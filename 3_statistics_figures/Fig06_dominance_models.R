# =============================================================================
# Fig. 6  Adjusted R2 of multiple regressions of observed full-community species
#         richness on the dominance metrics retained after collinearity
#         screening, at three spatial scales for grazing, fencing, and all
#         plots, with 95% cluster-bootstrap intervals.
#         a, dominance metrics from field species abundance (OB_Top)
#         b, dominance metrics from the deep learning species maps (DL_Top)
#         c, dominance metrics from the spectral clusters (SC_Top, k-means on
#            z-scored reflectance)
#
# Steps
#   1. For each source, scale and subset (grazing, fencing, all plots), screen
#      the five dominance metrics (Top1, Top2, Top3, Top12, Top123) for
#      collinearity on all observations (pairwise Spearman |rho| <= 0.8,
#      VIF <= 5) and fit observed SR on the retained metrics with ordinary least
#      squares (point estimate of adjusted R2).
#   2. Cluster bootstrap (2000 replicates): plots are resampled within
#      management units, each with all its quadrats and subplots; the retained
#      metrics are held fixed and the coefficients refitted. 95% percentile
#      intervals of adjusted R2.
#   3. Plot adjusted R2 against scale with the intervals as shaded bands;
#      adjusted R2 below zero is drawn as an open triangle on the baseline and
#      interval bounds below zero are cut at the baseline. Export the single
#      panels and the assembled figure.
#
# Input  : analysis_table.csv (columns Scale, Management, Unit, Plot15,
#          Observed_SR, OB_Top*, DL_Top*, SC_Top*)
# Output : output/Fig06_statistics.csv       retained metrics, adjusted R2 and
#                                            bootstrap interval of every model
#          output/Fig06_<a|b|c>.tiff         single panels (6.2 x 5.6 in)
#          output/Fig06.tiff                 assembled figure
#
# Usage  : Rscript Fig06_dominance_models.R [data_file] [output_dir]
# Default: data_file  = "../data/simulated/analysis_table.csv"
#          output_dir = "output"
# Runtime: about one minute
# =============================================================================

source("functions.R")
suppressPackageStartupMessages(library(patchwork))

args      <- commandArgs(trailingOnly = TRUE)
data_file <- if (length(args) >= 1) args[1] else file.path(DATA_DIR, "analysis_table.csv")
out_dir   <- if (length(args) >= 2) args[2] else "output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

B    <- 2000
SEED <- 20260908
SUBSETS <- c("Overall", "Grazing", "Fencing")
SOURCES <- list(a = list(prefix = "OB", title = "Observed-Top"),
                b = list(prefix = "DL", title = "Deep Learning-Top"),
                c = list(prefix = "SC", title = "Spectral Clustering-Top"))

dat <- read_table(data_file)
dat$member_plots <- dat$Plot15          # each row belongs to one plot
dat <- as.data.frame(dat)

# ---- 1. screening and full-data models ----------------------------------------------
models <- expand.grid(Scale = SCALES, Management = SUBSETS, Panel = names(SOURCES),
                      stringsAsFactors = FALSE)
specs <- vector("list", nrow(models))
point <- vector("list", nrow(models))
for (i in seq_len(nrow(models))) {
  sc <- models$Scale[i]; mg <- models$Management[i]; k <- models$Panel[i]
  rows <- which(dat$Scale == sc & (mg == "Overall" | dat$Management == mg))
  d <- dat[rows, ]
  s <- screen_fit(d, "Observed_SR", top_vars(SOURCES[[k]]$prefix))
  if (is.null(s$model)) next
  specs[[i]] <- list(rows = rows, X = model.matrix(formula(s$model)[-2], data = d), y = d$Observed_SR)
  point[[i]] <- data.frame(Panel = k, Source = SOURCES[[k]]$prefix, Scale = sc, Management = mg,
                           n = nrow(d), retained = paste(s$final_vars, collapse = ";"),
                           adjR2 = summary(s$model)$adj.r.squared)
}

# ---- 2. cluster bootstrap -----------------------------------------------------------
boot <- cluster_bootstrap(dat, specs, B = B, seed = SEED)
res <- do.call(rbind, lapply(seq_len(nrow(models)), function(i) {
  if (is.null(point[[i]])) return(NULL)
  q <- boot_ci(boot[, i])
  cbind(point[[i]], boot_lower = q[["lower"]], boot_upper = q[["upper"]], n_boot_valid = q[["n_valid"]])
}))
write.csv(res, file.path(out_dir, "Fig06_statistics.csv"), row.names = FALSE)

# ---- 3. figure -----------------------------------------------------------------------
BASE_SIZE    <- 25
MGMT3_COLS   <- c(Fencing = "#56B4E9", Grazing = "#E69F00", Overall = "black")      # lines and points
FILL3_COLS   <- c(Fencing = "#A9DBF5", Grazing = "#F5D48F", Overall = "#B8B8B8")    # bands
RIBBON_ALPHA <- 0.45
POINT_SIZE   <- 5.85                       # filled points
SCALE_LABELS <- c("0.5m" = "0.5m×0.5m", "1m" = "1m×1m", "1.5m" = "1.5m×1.5m")

d_plot <- res
d_plot$Scale <- factor(SCALE_LABELS[d_plot$Scale], levels = SCALE_LABELS)
d_plot$Management_label <- factor(d_plot$Management, levels = c("Fencing", "Grazing", "Overall"))
d_plot$below_zero <- d_plot$adjR2 < 0
d_plot$y_plot  <- pmax(d_plot$adjR2, 0)          # negative adjusted R2 drawn on the baseline
d_plot$lo_plot <- pmax(d_plot$boot_lower, 0)     # interval bounds below zero cut at the baseline
d_plot$hi_plot <- d_plot$boot_upper

pd <- position_dodge(width = 0.10)               # separates points and lines, not the bands

# y range per panel: upper limit = largest interval bound (rounded to two
# decimals) rounded up to 0.1; breaks every 0.2 above 0.6, otherwise every 0.1
y_axis_for <- function(df) {
  y_max <- ceiling(round(max(df$hi_plot, na.rm = TRUE), 2) * 10) / 10
  step  <- if (y_max > 0.6) 0.2 else 0.1
  list(limits = c(0, y_max), breaks = seq(0, y_max, by = step))
}

# layers: band -> zero line -> line -> point -> triangle
r2_panel <- function(df, title, y_title) {
  ya <- y_axis_for(df)
  df$Title <- title
  ggplot(df, aes(x = Scale, y = y_plot, group = Management_label, color = Management_label)) +
    geom_ribbon(aes(ymin = lo_plot, ymax = hi_plot, fill = Management_label),
                alpha = RIBBON_ALPHA, colour = NA, show.legend = FALSE) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.5) +
    geom_line(linewidth = 2.0, position = pd) +
    geom_point(data = ~ subset(.x, !below_zero), size = POINT_SIZE, position = pd) +
    geom_point(data = ~ subset(.x, below_zero), shape = 25, fill = "white",
               size = 5, stroke = 1.5, position = pd, show.legend = FALSE) +
    scale_color_manual(values = MGMT3_COLS, breaks = c("Fencing", "Grazing", "Overall"),
                       name = "Management") +
    scale_fill_manual(values = FILL3_COLS, guide = "none") +
    scale_y_continuous(breaks = ya$breaks, expand = expansion(mult = c(0.05, 0.05))) +
    coord_cartesian(ylim = ya$limits) +
    facet_grid(. ~ Title) +
    labs(x = "Scale", y = if (y_title) "R²" else NULL) +
    theme_nature(BASE_SIZE) +
    theme(legend.position = "none",
          strip.background = element_rect(fill = STRIP_FILL, colour = STRIP_LINE, linewidth = 0.8),
          strip.text = element_text(size = STRIP_SIZE, face = "bold", margin = ggplot2::margin(10, 8, 10, 8)),
          plot.margin = ggplot2::margin(4, 10, 4, 4))
}

STRIP_LINE <- "#969696"
STRIP_SIZE <- 23.5
panels <- lapply(names(SOURCES), function(k)
  r2_panel(d_plot[d_plot$Panel == k, ], SOURCES[[k]]$title, y_title = k == "a"))
names(panels) <- names(SOURCES)
for (k in names(panels)) save_tiff(panels[[k]], file.path(out_dir, sprintf("Fig06_%s.tiff", k)), 6.2, 5.6)

# assembled figure: panel letters, one legend below
LEGEND_SIZE <- 25
LETTER_SIZE <- 36
FIG_W <- 18.6; FIG_H <- 7.0   # inches
fig <- wrap_plots(panels, nrow = 1) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a") &
  guides(color = guide_legend(override.aes = list(linewidth = 2.0, size = 5))) &
  theme(legend.position = "bottom",
        legend.text = element_text(size = LEGEND_SIZE), legend.title = element_text(size = LEGEND_SIZE),
        legend.background = element_blank(),
        legend.key.width = unit(24, "pt"),       # wide enough for the line to show beside the point
        plot.tag = element_text(size = LETTER_SIZE, face = "bold", family = "Arial"),
        plot.tag.position = c(0.01, 0.99))
fig[[1]] <- fig[[1]] + theme(plot.tag.position = c(0.06, 0.99))   # letter a beside the y-axis title
save_tiff(fig, file.path(out_dir, "Fig06.tiff"), FIG_W, FIG_H)

cat("Fig. 6 written to", out_dir, "\n")
