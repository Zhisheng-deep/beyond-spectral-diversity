# =============================================================================
# Fig. 4  Relationships between remotely sensed and observed full-community
#         species richness (SR) at three spatial scales under grazing and
#         fencing.
#         a, SR predicted by deep learning (DL_SR)
#         b, spectral species richness (SSR_km_z, k-means on z-scored
#            reflectance)
#
# Steps
#   1. For each predictor and scale (0.5 m, 1 m, 1.5 m), regress observed SR
#      on the predictor with ordinary least squares for grazing, fencing and
#      all observations, and record adjusted R2 and the P value of the slope.
#   2. Plot the predictor against observed SR with least squares lines and
#      their 95% confidence bands, one panel per predictor and scale.
#   3. Export the single panels and the assembled figure.
#
# Input  : analysis_table.csv (columns Scale, Management, Observed_SR, DL_SR,
#          SSR_km_z)
# Output : output/Fig04_statistics.csv     adjusted R2 and P of every fit
#          output/Fig04_<a|b>_<scale>.tiff  single panels (6.2 x 5.6 in)
#          output/Fig04.tiff                assembled figure
#
# Usage  : Rscript Fig04_direct_SR.R [data_file] [output_dir]
# Default: data_file  = "../data/simulated/analysis_table.csv"
#          output_dir = "output"
# =============================================================================

source("functions.R")

args      <- commandArgs(trailingOnly = TRUE)
data_file <- if (length(args) >= 1) args[1] else file.path(DATA_DIR, "analysis_table.csv")
out_dir   <- if (length(args) >= 2) args[2] else "output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

data <- read_table(data_file)

BASE_SIZE  <- 20.5   # base font size of the panels
ANNOT_SIZE <- 6.5    # size of the R2 labels (mm)

# Placement of the R2 labels per scale (adjust to the data):
# labels  : corner, "bottom" or "top" right
# label_x : right end of the labels as a fraction of the panel width
PREDICTORS <- list(
  a = list(var = "DL_SR",    y_title = "Predicted Species Richness",
           header = "Deep Learning-Taxonomic Species",
           labels  = c("0.5m" = "bottom", "1m" = "bottom", "1.5m" = "bottom"),
           label_x = c("0.5m" = 1, "1m" = 1, "1.5m" = 1)),
  b = list(var = "SSR_km_z", y_title = "Spectral Species Richness",
           header = "Spectral Clustering-Spectral Species",
           labels  = c("0.5m" = "top", "1m" = "top", "1.5m" = "top"),
           label_x = c("0.5m" = 1, "1m" = 1, "1.5m" = 0.83))
)

# ---- 1. statistics ---------------------------------------------------------------
stats <- do.call(rbind, lapply(names(PREDICTORS), function(k) {
  v <- PREDICTORS[[k]]$var
  do.call(rbind, lapply(SCALES, function(sc) {
    d <- data[data$Scale == sc, ]
    groups <- list(Overall = d, Grazing = d[d$Management == "Grazing", ],
                   Fencing = d[d$Management == "Fencing", ])
    do.call(rbind, lapply(names(groups), function(g) {
      f <- ols_fit(y = groups[[g]]$Observed_SR, x = groups[[g]][[v]])
      data.frame(Panel = k, Predictor = v, Scale = sc, Group = g,
                 n = f[["n"]], adj_R2 = f[["adj_R2"]], P = f[["P"]])
    }))
  }))
}))
write.csv(stats, file.path(out_dir, "Fig04_statistics.csv"), row.names = FALSE)

# ---- 2. panels -------------------------------------------------------------------
# The three R2 labels are listed from top to bottom as Fencing, Grazing, Overall.
sr_panel <- function(k, sc, keep_y_title, labels, label_x = 1) {
  v  <- PREDICTORS[[k]]$var
  d  <- data[data$Scale == sc, ]
  st <- stats[stats$Panel == k & stats$Scale == sc, ]

  p <- ggplot(d, aes(x = Observed_SR, y = .data[[v]], color = Management)) +
    geom_point(alpha = 0.7, size = 2.8) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, aes(fill = Management), alpha = 0.25) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "black", fill = "grey70",
                alpha = 0.35, linetype = "dashed") +
    scale_color_manual(values = MGMT_COLS, breaks = c("Fencing", "Grazing")) +
    scale_fill_manual(values = MGMT_COLS, breaks = c("Fencing", "Grazing")) +
    scale_x_continuous(breaks = integer_breaks) +
    scale_y_continuous(breaks = integer_breaks) +
    labs(title = SCALE_TITLES[[sc]], x = "Observed Species Richness",
         y = if (keep_y_title) PREDICTORS[[k]]$y_title else NULL) +
    theme_nature(BASE_SIZE) +
    guides(color = "none", fill = "none")

  # R2 labels, from top to bottom Fencing, Grazing, Overall; vertical offsets
  # in text heights from the panel edge
  grp   <- c("Fencing", "Grazing", "Overall")
  text  <- sapply(grp, function(g) paste0(g, ": R\u00b2 = ", sprintf("%.2f", st$adj_R2[st$Group == g])))
  stars <- sapply(grp, function(g) { s <- p_to_stars(st$P[st$Group == g]); if (s == "NS") "" else s })
  vj    <- if (labels == "bottom") c(-4.25, -2.4, -0.55) else c(2.6, 4.45, 6.3)
  p + annotation_custom(r2_labels_grob(text, stars,
                                       colour = c(MGMT_COLS[["Fencing"]], MGMT_COLS[["Grazing"]], "black"),
                                       corner = labels, vjust = vj, x = label_x, size_mm = ANNOT_SIZE))
}

panels <- list()
for (k in names(PREDICTORS)) for (sc in SCALES) {
  panels[[k]][[sc]] <- sr_panel(k, sc, keep_y_title = sc == "0.5m",
                                labels = PREDICTORS[[k]]$labels[[sc]],
                                label_x = PREDICTORS[[k]]$label_x[[sc]])
}

# ---- 3. export -------------------------------------------------------------------
PANEL_W <- 6.2; PANEL_H <- 5.6       # single panel (in)
HEADER_H <- 0.64; LEGEND_H <- 0.48   # header and legend rows of the assembled figure (in)

for (k in names(panels)) for (sc in SCALES) {
  save_tiff(panels[[k]][[sc]], file.path(out_dir, sprintf("Fig04_%s_%s.tiff", k, sc)),
            PANEL_W, PANEL_H)
}

# legend taken from one panel drawn with its legend
LEGEND_SIZE <- 25
legend_plot <- panels$a[["0.5m"]] +
  guides(color = guide_legend(title = "Management", override.aes = list(
           colour = light_colour(MGMT_COLS[c("Fencing", "Grazing")]), linewidth = 2, size = 4)),
         fill  = guide_legend(title = "Management")) +
  theme(legend.text = element_text(size = LEGEND_SIZE), legend.title = element_text(size = LEGEND_SIZE),
        legend.key.width = unit(30, "pt"), legend.key.height = unit(20, "pt"))
legend <- get_legend_grob(legend_plot)

rows <- lapply(names(panels), function(k)
  cowplot::plot_grid(plotlist = panels[[k]], nrow = 1))
fig <- cowplot::plot_grid(
  panel_header("a", PREDICTORS$a$header), rows[[1]],
  panel_header("b", PREDICTORS$b$header), rows[[2]],
  legend,
  ncol = 1, rel_heights = c(HEADER_H, PANEL_H, HEADER_H, PANEL_H, LEGEND_H))
fig <- fig + theme(plot.background = element_rect(fill = "white", colour = NA))
save_tiff(fig, file.path(out_dir, "Fig04.tiff"), 3 * PANEL_W, 2 * (HEADER_H + PANEL_H) + LEGEND_H)

cat("Fig. 4 written to", out_dir, "\n")
