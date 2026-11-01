# =============================================================================
# Fig. 5  Relationships between remotely sensed and observed species richness
#         (SR) of the sunlight-exposed species and of the full community at
#         three spatial scales for grazing, fencing, and all plots.
#         a, SR predicted by deep learning (DL_SR)
#         b, spectral species richness (SSR_km_z, k-means on z-scored
#            reflectance)
#
# Steps
#   1. For each predictor, scale and subset (grazing, fencing, all plots),
#      regress observed sunlight-exposed SR and observed full-community SR on
#      the predictor with ordinary least squares, recording adjusted R2 and the
#      P value of the slope, and calculate the mean absolute error (MAE)
#      between the predictor and each observed SR.
#   2. Plot the predictor against both observed SR with least squares lines,
#      their 95% confidence bands and the 1:1 line, one panel per scale and
#      subset. Tick intervals are shared by the panels of a scale (x axis) and
#      of a subset (y axis).
#   3. Export the single panels and the assembled figure.
#
# Input  : analysis_table.csv (columns Scale, Management, Observed_SR,
#          Observed_Visible_SR, DL_SR, SSR_km_z)
# Output : output/Fig05_statistics.csv                adjusted R2, P and MAE
#          output/Fig05_<a|b>_<scale>_<subset>.tiff    single panels (6.2 x 5.6 in)
#          output/Fig05a.tiff, Fig05b.tiff             assembled figure, a and b
#                                                      (legend below b)
#
# Usage  : Rscript Fig05_visible_full_SR.R [data_file] [output_dir]
# Default: data_file  = "../data/simulated/analysis_table.csv"
#          output_dir = "output"
# =============================================================================

source("functions.R")
suppressPackageStartupMessages(library(patchwork))

args      <- commandArgs(trailingOnly = TRUE)
data_file <- if (length(args) >= 1) args[1] else file.path(DATA_DIR, "analysis_table.csv")
out_dir   <- if (length(args) >= 2) args[2] else "output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

data <- read_table(data_file)

BASE_SIZE  <- 25
ANNOT_SIZE <- 8.1    # size of the R2 and MAE labels (mm)
SR_COLS    <- c(Visible = "#56B4E9", Full = "#E69F00")   # sunlight-exposed, full community
SUBSETS    <- c("Grazing", "Fencing", "Overall")
OBSERVED   <- c(Full = "Observed_SR", Visible = "Observed_Visible_SR")

# Labels of R2 and MAE in each panel:
#   style  "stack" : R2 and MAE on separate lines (panel a)
#          "inline": R2 and MAE on one line per observed SR (panel b)
#   corner "bottomright" or "topleft" (adjust to the data)
PREDICTORS <- list(
  a = list(var = "DL_SR",    y_title = "Predicted Species Richness", style = "stack",
           corner = "topleft"),
  b = list(var = "SSR_km_z", y_title = "Spectral Species Richness",  style = "inline",
           corner = "bottomright")
)

subset_data <- function(d, s) if (s == "Overall") d else d[d$Management == s, ]

# ---- 1. statistics -------------------------------------------------------------
stats <- do.call(rbind, lapply(names(PREDICTORS), function(k) {
  v <- PREDICTORS[[k]]$var
  do.call(rbind, lapply(SCALES, function(sc) {
    do.call(rbind, lapply(SUBSETS, function(s) {
      d <- subset_data(data[data$Scale == sc, ], s)
      do.call(rbind, lapply(names(OBSERVED), function(o) {
        f <- ols_fit(y = d[[OBSERVED[[o]]]], x = d[[v]])
        data.frame(Panel = k, Predictor = v, Scale = sc, Subset = s,
                   Observed = OBSERVED[[o]], n = f[["n"]], adj_R2 = f[["adj_R2"]], P = f[["P"]],
                   MAE = mean(abs(d[[v]] - d[[OBSERVED[[o]]]])))
      }))
    }))
  }))
}))
write.csv(stats, file.path(out_dir, "Fig05_statistics.csv"), row.names = FALSE)

# ---- 2. panels -----------------------------------------------------------------
# Integer tick interval giving about five ticks over the range of x
choose_tick_by <- function(x, target_n = 5) {
  x <- x[is.finite(x)]
  span <- diff(range(x))
  if (span <= 0) return(1)
  candidates <- c(1, 2, 3, 4, 5, 6, 8, 10, 12, 15, 20, 25, 30)
  candidates[which.min(abs(candidates - span / target_n))]
}
make_breaks_by <- function(by) {
  force(by)
  function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0) return(NULL)
    seq(floor(min(x) / by) * by, ceiling(max(x) / by) * by, by = by)
  }
}

# R2 and MAE labels of one panel
sr_labels_grob <- function(st, style, corner, size_mm) {
  gp <- function(col) grid::gpar(fontfamily = "Arial", fontface = "bold",
                                 fontsize = size_mm * .pt, col = col)
  r2  <- function(o) {
    s <- p_to_stars(st$P[st$Observed == OBSERVED[[o]]])
    paste0("R\u00b2 = ", sprintf("%.2f", st$adj_R2[st$Observed == OBSERVED[[o]]]),
           if (s == "NS") "" else paste0(" (", s, ")"))
  }
  mae <- function(o) paste0("MAE = ", sprintf("%.1f", st$MAE[st$Observed == OBSERVED[[o]]]))
  width <- function(txt) do.call(max, lapply(txt, function(t) grid::grobWidth(grid::textGrob(t, gp = gp("black")))))
  pad <- grid::unit(4, "pt")
  top <- corner == "topleft"
  y   <- grid::unit(if (top) 1 else 0, "npc")
  # vj: line offsets from the panel edge in text heights, listed from the top line
  if (style == "stack") {
    txt <- c(r2("Full"), mae("Full"), r2("Visible"), mae("Visible"))
    col <- SR_COLS[c("Full", "Full", "Visible", "Visible")]
    x0  <- if (top) grid::unit(12, "pt") else grid::unit(1, "npc") - pad - width(txt)
    vj  <- if (top) c(1.3, 3.16, 5.03, 6.9) else c(-6.9, -5.03, -3.16, -1.3)
    grobs <- lapply(1:4, function(i) grid::textGrob(txt[i], x = x0, y = y,
                                                    hjust = 0, vjust = vj[i], gp = gp(col[[i]])))
  } else {
    left  <- c(paste0(r2("Full"), ";"), paste0(r2("Visible"), ";"))
    right <- c(mae("Full"), mae("Visible"))
    col   <- SR_COLS[c("Full", "Visible")]
    gap   <- grid::unit(10, "pt")
    col1  <- width("R\u00b2 = 0.00 (***);")                     # width of the R2 column
    x0    <- if (top) grid::unit(12, "pt") else grid::unit(1, "npc") - pad - width(right) - gap - col1
    x1    <- x0 + col1 + gap                                     # start of the MAE column
    vj    <- if (top) c(1.6, 3.45) else c(-3.45, -1.6)
    grobs <- c(lapply(1:2, function(i) grid::textGrob(left[i], x = x0, y = y,
                                                      hjust = 0, vjust = vj[i], gp = gp(col[[i]]))),
               lapply(1:2, function(i) grid::textGrob(right[i], x = x1, y = y,
                                                      hjust = 0, vjust = vj[i], gp = gp(col[[i]]))))
  }
  grid::gTree(children = do.call(grid::gList, grobs))
}

sr_panel <- function(k, sc, s, x_by, y_by, x_title, y_title, top_strip, right_strip) {
  v <- PREDICTORS[[k]]$var
  d <- subset_data(data[data$Scale == sc, ], s)
  # long format, the two observed SR of each sampling unit in consecutive rows
  long <- data.frame(
    Predicted = rep(d[[v]], each = 2),
    Observed  = as.vector(rbind(d$Observed_SR, d$Observed_Visible_SR)),
    SR_type   = factor(rep(c("Full", "Visible"), nrow(d)), levels = c("Full", "Visible")),
    Column    = SCALE_TITLES[[sc]], Row = s)
  st <- stats[stats$Panel == k & stats$Scale == sc & stats$Subset == s, ]

  ggplot(long, aes(x = Observed, y = Predicted, color = SR_type)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth = 1.1) +
    geom_point(alpha = 0.7, size = 2.8) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, aes(fill = SR_type), alpha = 0.25) +
    scale_color_manual(values = SR_COLS) +
    scale_fill_manual(values = SR_COLS) +
    scale_x_continuous(breaks = make_breaks_by(x_by), expand = expansion(mult = c(0.04, 0.10))) +
    scale_y_continuous(breaks = make_breaks_by(y_by), expand = expansion(mult = c(0.06, 0.08))) +
    facet_grid(Row ~ Column) +
    labs(x = if (x_title) "Observed Species Richness" else NULL,
         y = if (y_title) PREDICTORS[[k]]$y_title else NULL) +
    annotation_custom(sr_labels_grob(st, PREDICTORS[[k]]$style, PREDICTORS[[k]]$corner, ANNOT_SIZE)) +
    theme_nature(BASE_SIZE) +
    theme(legend.position = "none",
          strip.background = element_rect(fill = STRIP_FILL, colour = STRIP_LINE, linewidth = 0.8),
          strip.text = element_text(size = STRIP_SIZE, face = "bold", margin = ggplot2::margin(8, 8, 8, 8)),
          strip.text.y = element_text(angle = -90),
          strip.text.x = if (top_strip) element_text(size = STRIP_SIZE, face = "bold") else element_blank(),
          strip.background.x = if (top_strip) element_rect(fill = STRIP_FILL, colour = STRIP_LINE, linewidth = 0.8) else element_blank(),
          strip.text.y.right = if (right_strip) element_text(size = STRIP_SIZE, face = "bold", angle = -90) else element_blank(),
          strip.background.y = if (right_strip) element_rect(fill = STRIP_FILL, colour = STRIP_LINE, linewidth = 0.8) else element_blank(),
          plot.background = element_blank(),     # tick labels at the edge stay visible
          plot.margin = ggplot2::margin(PANEL_MARGIN, PANEL_MARGIN, PANEL_MARGIN, PANEL_MARGIN))
}

STRIP_LINE   <- "#969696"
STRIP_SIZE   <- 32
PANEL_MARGIN <- 6

blocks <- list()
for (k in names(PREDICTORS)) {
  v <- PREDICTORS[[k]]$var
  # tick intervals: x shared within a scale, y shared within a subset
  x_by <- sapply(SCALES, function(sc) {
    d <- data[data$Scale == sc, ]; choose_tick_by(c(d$Observed_SR, d$Observed_Visible_SR)) })
  y_by <- sapply(SUBSETS, function(s) choose_tick_by(subset_data(data, s)[[v]]))
  panels <- list()
  for (s in SUBSETS) for (sc in SCALES) {
    p <- sr_panel(k, sc, s, x_by[[sc]], y_by[[s]],
                  x_title = s == "Overall", y_title = sc == "0.5m",
                  top_strip = s == "Grazing", right_strip = sc == "1.5m")
    panels[[length(panels) + 1]] <- p
    save_tiff(p, file.path(out_dir, sprintf("Fig05_%s_%s_%s.tiff", k, sc, s)), 6.2, 5.6)
  }
  blocks[[k]] <- wrap_plots(panels, ncol = 3)
}

# ---- 3. assembled figure ---------------------------------------------------------
LEGEND_SIZE <- 27
SR_LABELS   <- c(Visible = "Sunlight-exposed Species Richness", Full = "Full-community Species Richness")
legend_data <- data.frame(x = rep(1:3, 2), y = rep(1:3, 2),
                          SR_type = factor(rep(c("Visible", "Full"), each = 3), levels = c("Visible", "Full")))
key_aes <- list(colour = light_colour(SR_COLS[c("Visible", "Full")]), linewidth = 2,
                size = 5, linetype = "solid")
legend_plot <- ggplot(legend_data, aes(x, y, color = SR_type, fill = SR_type)) +
  geom_point(alpha = 0.7, size = 2.8) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, alpha = 0.25) +
  geom_line(aes(linetype = "1:1 Line", group = SR_type), color = "black", linewidth = 1.1) +
  scale_color_manual(values = SR_COLS, labels = SR_LABELS, name = NULL) +
  scale_fill_manual(values = SR_COLS, labels = SR_LABELS, name = NULL) +
  scale_linetype_manual(values = c("1:1 Line" = "51"), name = NULL) +
  guides(color = guide_legend(order = 1, override.aes = key_aes),
         fill  = guide_legend(order = 1),
         linetype = guide_legend(order = 2, keywidth = unit(60, "pt"),
                                 override.aes = list(linewidth = 2.2, fill = NA))) +
  theme_nature(BASE_SIZE) +
  theme(legend.position = "bottom", legend.text = element_text(size = LEGEND_SIZE),
        legend.key.width = unit(55, "pt"), legend.key.height = unit(36, "pt"),
        legend.spacing.x = unit(6, "pt"), legend.box.spacing = unit(0, "pt"))
legend <- get_legend_grob(legend_plot)

block_with_letter <- function(block, letter) {
  cowplot::ggdraw(with_agg(patchwork::patchworkGrob(block))) +
    cowplot::draw_label(letter, x = 0.006, y = 0.997, hjust = 0, vjust = 1,
                        fontface = "bold", fontfamily = "Arial", size = LETTER_SIZE)
}
LETTER_SIZE <- 44
BLOCK_W <- 17.5; BLOCK_H <- 16.84; LEGEND_H <- 0.85   # inches

# a and b are separate files; the legend is placed below b
fig_a <- block_with_letter(blocks$a, "a") + theme(plot.background = element_rect(fill = "white", colour = NA))
save_tiff(fig_a, file.path(out_dir, "Fig05a.tiff"), BLOCK_W, BLOCK_H)
fig_b <- cowplot::plot_grid(block_with_letter(blocks$b, "b"), legend, ncol = 1,
                            rel_heights = c(BLOCK_H, LEGEND_H)) +
  theme(plot.background = element_rect(fill = "white", colour = NA))
save_tiff(fig_b, file.path(out_dir, "Fig05b.tiff"), BLOCK_W, BLOCK_H + LEGEND_H)

cat("Fig. 5 written to", out_dir, "\n")
