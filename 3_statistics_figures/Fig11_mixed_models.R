# =============================================================================
# Fig. 11  Contributions of management regime, canopy structure, dominance,
#          intra-species spectral diversity (SDiv), and compositional dynamics
#          to observed full-community species richness (SR) in linear mixed
#          models with the dominance metrics from deep learning (DL) or spectral
#          clustering (SC), at three spatial scales on all plots.
#          a, sequential decomposition of the marginal R2 (R2m)
#          b, hierarchical partitioning of the R2m of the full model and
#             standardized coefficients of the fixed effects
#
# Steps
#   1. For each scale and approach (DL_Top* or SC_Top*), screen the candidate
#      variables jointly for collinearity (pairwise Spearman |rho| <= 0.8,
#      VIF <= 5) and standardize the retained variables within scale. Groups:
#      management regime; canopy structure (FVC, Height_Mean, Height_CV);
#      dominance (Top1, Top2, Top3, Top12, Top123); intra-species SDiv
#      (Intra_Species_SDiv); compositional dynamics (Rare_Species_Ratio and
#      Turnover, which is undefined at 0.5 m).
#   2. Fit linear mixed models of log observed SR by maximum likelihood with
#      random intercepts for management unit and plot within unit at 0.5 m and
#      1 m and for management unit at 1.5 m. Add the groups in the order above;
#      R2m and conditional R2 (R2c) of each model, and a likelihood ratio test
#      of each added group against the preceding model.
#   3. Partition the R2m of the full model hierarchically among the groups
#      (glmm.hp) and test the standardized coefficients of its fixed effects
#      with Satterthwaite degrees of freedom (lmerTest; 95% confidence
#      intervals).
#   4. Figures: a, stacked increments of R2m with asterisks from the likelihood
#      ratio tests; b, stacked percentages of R2m and coefficients with 95%
#      confidence intervals. Export the single panels and the two assembled
#      figures.
#
# Input  : analysis_table.csv (columns Scale, Management, Unit, Plot15,
#          Observed_SR, FVC, Height_Mean, Height_CV, DL_Top*, SC_Top*,
#          Intra_Species_SDiv, Rare_Species_Ratio, Turnover)
# Output : output/Fig11_sequential.csv     R2m, R2c, increments and likelihood
#                                          ratio tests of the sequential models
#          output/Fig11_partition.csv      hierarchical partitioning of R2m
#          output/Fig11_coefficients.csv   standardized coefficients of the full
#                                          models
#          output/Fig11a_<scale>.tiff, Fig11b_<DL|SC>_<scale>.tiff
#                                          single panels (7.5 x 5 in)
#          output/Fig11a.tiff, Fig11b.tiff assembled figures
#
# Usage  : Rscript Fig11_mixed_models.R [data_file] [output_dir]
# Default: data_file  = "../data/simulated/analysis_table.csv"
#          output_dir = "output"
# =============================================================================

source("functions.R")
suppressPackageStartupMessages({
  library(lme4); library(lmerTest); library(MuMIn); library(glmm.hp); library(patchwork)
})

args      <- commandArgs(trailingOnly = TRUE)
data_file <- if (length(args) >= 1) args[1] else file.path(DATA_DIR, "analysis_table.csv")
out_dir   <- if (length(args) >= 2) args[2] else "output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

RESP     <- "Observed_SR"
OFFSET   <- 0.01
PATHWAYS <- c("DL", "SC")
BLOCKS <- function(path, sc) list(
  Structure = c("FVC", "Height_Mean", "Height_CV"),
  Dominance = top_vars(path),
  SDiv      = "Intra_Species_SDiv",
  Context   = if (sc == "0.5m") "Rare_Species_Ratio" else c("Turnover", "Rare_Species_Ratio"))
BLOCK_ORDER <- c("Management", "Structure", "Dominance", "SDiv", "Context")

dat <- read_table(data_file)
dat$Management <- relevel(factor(dat$Management), ref = "Grazing")
dat$Unit <- factor(dat$Unit); dat$Plot15 <- factor(dat$Plot15)

r2 <- function(m) {
  r <- suppressWarnings(MuMIn::r.squaredGLMM(m))
  c(R2m = unname(r[1, "R2m"]), R2c = unname(r[1, "R2c"]))
}

# ---- 1-3. statistics ---------------------------------------------------------------
seq_out <- list(); hp_out <- list(); coef_out <- list()
for (sc in SCALES) for (path in PATHWAYS) {
  blocks <- BLOCKS(path, sc)
  cand <- unlist(blocks, use.names = FALSE)
  d <- dat[dat$Scale == sc, ]
  d <- d[complete.cases(d[, c(RESP, "Management", "Unit", "Plot15", cand)]), ]

  # joint screening; retained variables are standardized and returned to their groups
  s <- screen_fit(d, RESP, cand)
  kept <- lapply(blocks, function(v) intersect(v, s$final_vars))
  for (v in s$final_vars) d[[paste0("z_", v)]] <- as.numeric(scale(d[[v]]))
  d$logSR <- log(d[[RESP]] + OFFSET)
  zk <- lapply(kept, function(v) if (length(v)) paste0("z_", v) else character(0))
  re <- if (sc == "1.5m") "(1 | Unit)" else "(1 | Unit/Plot15)"

  # sequential models
  steps <- list(M0 = "Management",
                M1 = c("Management", zk$Structure),
                M2 = c("Management", zk$Structure, zk$Dominance),
                M3 = c("Management", zk$Structure, zk$Dominance, zk$SDiv),
                M4 = c("Management", zk$Structure, zk$Dominance, zk$SDiv, zk$Context))
  prev <- NULL; prev_r2 <- c(R2m = 0, R2c = 0)
  for (k in seq_along(steps)) {
    m <- lme4::lmer(as.formula(paste("logSR ~", paste(steps[[k]], collapse = " + "), "+", re)),
                    data = d, REML = FALSE, control = lmerControl(check.conv.singular = "ignore"))
    rr <- r2(m)
    added <- if (k == 1) "Management" else setdiff(steps[[k]], steps[[k - 1]])
    lrt <- if (is.null(prev) || length(added) == 0) c(NA, NA, NA) else {
      a <- anova(prev, m); c(a$Chisq[2], a$Df[2], a$`Pr(>Chisq)`[2]) }
    seq_out[[length(seq_out) + 1]] <- data.frame(
      Scale = sc, Pathway = path, step = names(steps)[k], block = BLOCK_ORDER[k], n = nrow(d),
      terms_added = paste(sub("^z_", "", added), collapse = ";"),
      R2m = rr[["R2m"]], R2c = rr[["R2c"]], dR2m = rr[["R2m"]] - prev_r2[["R2m"]],
      LRT_chisq = lrt[1], LRT_df = lrt[2], LRT_P = lrt[3])
    prev <- m; prev_r2 <- rr
  }

  # hierarchical partitioning of the R2m of the full model among the groups
  grp <- list(Management = "Management", Structure = zk$Structure, Dominance = zk$Dominance,
              SDiv = zk$SDiv, Context = zk$Context)
  grp <- grp[lengths(grp) > 0]
  hp  <- as.data.frame(glmm.hp(m, grp)$hierarchical.partitioning)
  ci  <- grep("^I$|^Individual", names(hp), value = TRUE)[1]
  cp  <- grep("perc", names(hp), value = TRUE)[1]
  hp_out[[length(hp_out) + 1]] <- data.frame(Scale = sc, Pathway = path, block = rownames(hp),
                                             I = hp[[ci]], I_perc = hp[[cp]],
                                             R2m = rr[["R2m"]], R2c = rr[["R2c"]])

  # standardized coefficients of the full model (Satterthwaite degrees of freedom)
  mt <- lmerTest::lmer(formula(m), data = d, REML = FALSE, control = lmerControl(check.conv.singular = "ignore"))
  ct <- summary(mt)$coefficients
  ct <- ct[rownames(ct) != "(Intercept)", , drop = FALSE]
  vars <- ifelse(grepl("^Management", rownames(ct)), "Management", sub("^z_", "", rownames(ct)))
  coef_out[[length(coef_out) + 1]] <- data.frame(
    Scale = sc, Pathway = path, variable = vars,
    block = vapply(vars, function(v) { b <- names(kept)[vapply(kept, function(x) v %in% x, logical(1))]
                                       if (length(b) == 0) "Management" else b }, character(1)),
    estimate = ct[, "Estimate"], std.error = ct[, "Std. Error"], df = ct[, "df"],
    conf.low  = ct[, "Estimate"] - qt(0.975, ct[, "df"]) * ct[, "Std. Error"],
    conf.high = ct[, "Estimate"] + qt(0.975, ct[, "df"]) * ct[, "Std. Error"],
    p.value = ct[, "Pr(>|t|)"], row.names = NULL)
}
seq_out <- do.call(rbind, seq_out); hp_out <- do.call(rbind, hp_out); coef_out <- do.call(rbind, coef_out)
write.csv(seq_out,  file.path(out_dir, "Fig11_sequential.csv"),   row.names = FALSE)
write.csv(hp_out,   file.path(out_dir, "Fig11_partition.csv"),    row.names = FALSE)
write.csv(coef_out, file.path(out_dir, "Fig11_coefficients.csv"), row.names = FALSE)

# ---- 4. figures ----------------------------------------------------------------------
BASE_SIZE <- 22
GROUPS <- c(Management = "Management", Context = "Compositional Dynamics", Structure = "Canopy Structure",
            SDiv = "Intra-species SDiv", Dominance = "Dominance of Top Species")
GROUP_COLS <- c("Management" = "#D16E3D", "Compositional Dynamics" = "#C0E3DD",
                "Canopy Structure" = "#A19700", "Intra-species SDiv" = "#98B76A",
                "Dominance of Top Species" = "#FDD47B")
SEQ_ORDER <- GROUPS[c("Management", "Structure", "Dominance", "SDiv", "Context")]   # bottom to top in a
HP_ORDER  <- unname(GROUPS)                                                       # bottom to top in b
PATH_LABELS <- c(DL = "Deep Learning", SC = "Spectral Clustering")

theme_fig11 <- function(base_size = BASE_SIZE) {
  theme_bw(base_size = base_size) +
    theme(
      text = element_text(family = "Arial", face = "bold", colour = "black"),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.25),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      axis.line  = element_line(color = "black", linewidth = 0.8),
      axis.ticks = element_line(color = "black", linewidth = 0.8),
      axis.text  = element_text(size = base_size - 2, face = "bold", colour = "grey30"),
      axis.title.x = element_text(margin = ggplot2::margin(t = 0.4, unit = "cm"), size = base_size,
                                  face = "bold", colour = "black"),
      axis.title.y = element_text(margin = ggplot2::margin(r = 0.4, unit = "cm"), size = base_size,
                                  face = "bold", colour = "black"),
      legend.position = "none",
      plot.title = element_blank(),
      strip.text = element_text(size = base_size, face = "bold", colour = "black"),
      plot.margin = ggplot2::margin(0.5, 0.5, 0.5, 0.5, unit = "cm"))
}

# ---- a: sequential increments of R2m --------------------------------------------------
seq_plot <- seq_out
seq_plot$Block   <- factor(GROUPS[seq_plot$block], levels = SEQ_ORDER)
seq_plot$Pathway <- factor(PATH_LABELS[seq_plot$Pathway], levels = PATH_LABELS)
seq_plot$star    <- vapply(seq_plot$LRT_P, function(p) { s <- p_to_stars(p); if (s == "NS") "" else s }, character(1))

STAR_SIZE <- 21.2          # asterisks (pt), centred on their segment
STAR_RISE <- -0.0101       # vertical offset of the asterisk glyph (data units)
star_positions <- function(df) {
  pos <- df[df$dR2m > 0, ]
  pos <- pos[order(pos$Pathway, pos$Block), ]
  pos$top  <- ave(pos$dR2m, pos$Pathway, FUN = cumsum)
  pos$ymid <- pos$top - pos$dR2m / 2
  pos[pos$star != "", ]
}

seq_panel <- function(sc, stars = FALSE) {
  df <- seq_plot[seq_plot$Scale == sc, ]
  p <- ggplot(df, aes(x = Pathway, y = dR2m, fill = Block)) +
    geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.5) +
    geom_col(width = 0.4, colour = "black", linewidth = 0.3, position = position_stack(reverse = TRUE))
  if (stars) {
    st <- star_positions(df)
    p <- p + geom_text(data = st, aes(x = Pathway, y = ymid + STAR_RISE, label = star), inherit.aes = FALSE,
                       size = STAR_SIZE / .pt, fontface = "bold", family = "Arial")
  }
  p + scale_fill_manual(values = GROUP_COLS, name = NULL, drop = FALSE) +
    scale_y_continuous(breaks = seq(0, 0.8, 0.2), expand = expansion(mult = c(0, 0.02))) +
    coord_cartesian(ylim = c(-0.02, 0.85)) +
    labs(x = NULL, y = "Marginal R\u00b2 (sequential)") +
    theme_fig11() +
    theme(axis.text.x = element_text(size = BASE_SIZE, colour = "black"))
}

PANEL_W <- 7.5; PANEL_H <- 5     # single panels (inches)
for (sc in SCALES) save_tiff(seq_panel(sc), file.path(out_dir, sprintf("Fig11a_%s.tiff", sc)), PANEL_W, PANEL_H)

# ---- b: hierarchical partitioning and coefficients -----------------------------------
# percentages of R2m
hp_plot <- merge(expand.grid(Scale = SCALES, Pathway = PATHWAYS, block = names(GROUPS), stringsAsFactors = FALSE),
                 hp_out[, c("Scale", "Pathway", "block", "I_perc")], all.x = TRUE)
hp_plot$I_perc[is.na(hp_plot$I_perc)] <- 0
hp_plot <- do.call(rbind, lapply(split(hp_plot, paste(hp_plot$Scale, hp_plot$Pathway)), function(g) {
  if (any(g$I_perc < 0)) g$I_perc <- pmax(g$I_perc, 0)
  g$perc <- if (sum(g$I_perc) > 0) g$I_perc / sum(g$I_perc) * 100 else 0
  # integer percentages for the labels, summing to 100
  fl <- floor(g$perc); r <- 100 - sum(fl)
  if (r > 0) fl[order(g$perc - fl, decreasing = TRUE)[seq_len(r)]] <- fl[order(g$perc - fl, decreasing = TRUE)[seq_len(r)]] + 1
  g$label <- fl
  g
}))
hp_plot$Group <- factor(GROUPS[hp_plot$block], levels = HP_ORDER)

VAR_ORDER <- c("Management", "Turnover", "Rare_Species_Ratio", "FVC", "Height_Mean", "Height_CV",
               "Intra_Species_SDiv", "Top1", "Top2", "Top3", "Top12", "Top123")        # bottom to top
VAR_LABELS <- c(Rare_Species_Ratio = "Rare Species Ratio", Intra_Species_SDiv = "Intra-species SDiv",
                Height_CV = "Height CV", Height_Mean = "Height Mean")
coef_plot <- coef_out
coef_plot$term  <- sub("^(DL|SC)_", "", coef_plot$variable)
coef_plot$Group <- factor(GROUPS[coef_plot$block], levels = HP_ORDER)
coef_plot$sig   <- vapply(coef_plot$p.value, function(p) { s <- p_to_stars(p); if (s == "NS") "" else s }, character(1))
coef_plot$label <- paste0(ifelse(coef_plot$term %in% names(VAR_LABELS), VAR_LABELS[coef_plot$term], coef_plot$term),
                          " ", coef_plot$sig)

NUM_SIZE  <- 21.7          # percentages in the bars (pt)
NUM_WHITE <- c("Management", "Canopy Structure", "Intra-species SDiv")
NUM_ABOVE <- 5.2           # a top segment thinner than this (percent) is labelled above the bar
NUM_GAP   <- 1.3           # gap between the bar and a label above it (percent)
number_positions <- function(df) {
  df <- df[order(df$Group), ]
  df$top  <- cumsum(df$perc); df$ymid <- df$top - df$perc / 2
  df <- df[df$label > 0, ]
  top_seg <- which.max(df$top)
  df$above  <- seq_len(nrow(df)) == top_seg & df$perc < NUM_ABOVE
  df$y      <- ifelse(df$above, 100 + NUM_GAP, pmax(df$ymid, NUM_ABOVE / 2))   # within the bar
  df$vjust  <- ifelse(df$above, 0, 0.5)
  df$colour <- ifelse(!df$above & as.character(df$Group) %in% NUM_WHITE, "white", "black")
  df
}

hp_panel <- function(sc, path, numbers = FALSE) {
  bars <- hp_plot[hp_plot$Scale == sc & hp_plot$Pathway == path, ]
  cf <- coef_plot[coef_plot$Scale == sc & coef_plot$Pathway == path, ]
  cf <- cf[order(factor(cf$term, levels = VAR_ORDER)), ]
  cf$label <- factor(cf$label, levels = unique(cf$label))
  p_bar <- ggplot(bars, aes(x = 0.75, y = perc, fill = Group)) +
    geom_bar(stat = "identity", width = 0.4, color = "black", linewidth = 0.3,
             position = position_stack(reverse = TRUE))
  if (numbers) {
    nb <- number_positions(bars)
    p_bar <- p_bar + geom_text(data = nb, aes(x = 0.75, y = y, label = label, vjust = vjust, colour = colour),
                               inherit.aes = FALSE, size = NUM_SIZE / .pt, family = "Arial") +
      scale_colour_identity()
  }
  p_bar <- p_bar +
    scale_fill_manual(values = GROUP_COLS, drop = FALSE) +
    scale_y_continuous(breaks = seq(0, 100, by = 20), expand = expansion(mult = c(0, 0))) +
    scale_x_continuous(breaks = NULL, expand = expansion(mult = c(0, 0))) +
    coord_cartesian(xlim = c(0.45, 1.05), ylim = c(0, 100), clip = if (numbers) "off" else "on") +
    labs(x = NULL, y = if (sc == "0.5m") "Percentage Contribution (%)" else "") +
    theme_fig11() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.line.x = element_blank(), axis.title.x = element_blank(),
          axis.text.y = element_text(size = BASE_SIZE - 2, colour = "grey30"),
          axis.title.y = element_text(colour = "black"),
          panel.border = element_blank(),
          panel.background = element_rect(fill = "transparent", color = NA),
          plot.background  = element_rect(fill = "transparent", color = NA),
          plot.margin = ggplot2::margin(2, 0, 2, 0, "mm"),
          axis.title.y.left = element_text(margin = ggplot2::margin(t = -2, r = -4.5, b = -2, l = -2)))
  p_coef <- ggplot(cf, aes(x = estimate, y = label, color = Group)) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 2.0, color = "grey40") +
    geom_point(size = 7, show.legend = FALSE) +
    geom_errorbar(aes(xmin = conf.low, xmax = conf.high), width = 0, orientation = "y", linewidth = 2.5,
                  show.legend = FALSE) +
    scale_color_manual(values = GROUP_COLS, drop = FALSE) +
    scale_y_discrete(position = "right", expand = expansion(add = c(0.5, 0.5))) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 4)) +
    labs(x = "Estimate", y = NULL) +
    theme_fig11() +
    theme(axis.text.y = element_text(size = BASE_SIZE - 2, colour = "black", face = "bold",
                                     hjust = 0, margin = ggplot2::margin(l = 6)),
          axis.line.y = element_blank(), axis.ticks.y = element_blank(),
          panel.border = element_blank(),
          axis.text.x = element_text(size = BASE_SIZE - 2, colour = "grey30"),
          plot.margin = ggplot2::margin(2, 0, 2, 0, "mm"),
          axis.title.x = element_text(margin = ggplot2::margin(t = 2)),
          panel.spacing = unit(1, "mm"))
  p_bar + p_coef + plot_layout(ncol = 2, widths = c(0.6, 2.5)) &
    theme(plot.margin = ggplot2::margin(2, 0, 0, 0.5, "mm"))
}

for (path in PATHWAYS) for (sc in SCALES)
  save_tiff(hp_panel(sc, path), file.path(out_dir, sprintf("Fig11b_%s_%s.tiff", path, sc)), PANEL_W, PANEL_H)

# ---- assembled figures -----------------------------------------------------------------
# As in the published figures, the single panels are placed without rescaling.
# Positions and sizes below are in inches from the top-left corner.
gp_text <- function(size, col = "black", face = "bold") grid::gpar(fontfamily = "Arial", fontface = face,
                                                                 fontsize = size, col = col)
in_w <- function(u) with_agg(grid::convertWidth(sum(u), "in", valueOnly = TRUE))
in_h <- function(u) with_agg(grid::convertHeight(sum(u), "in", valueOnly = TRUE))
draw_tiff <- function(grobs, file, w, h) {
  ragg::agg_tiff(file, width = w, height = h, units = "in", res = 600, compression = "lzw")
  grid::grid.draw(grid::gTree(children = do.call(grid::gList, grobs)))
  invisible(dev.off())
}
make_box <- function(fig_h) function(x0, y0, x1, y1, fill, col = NA, lwd = 1)
  grid::rectGrob(x = grid::unit(x0, "in"), y = grid::unit(fig_h - y1, "in"),
                 width = grid::unit(x1 - x0, "in"), height = grid::unit(y1 - y0, "in"),
                 just = c("left", "bottom"), gp = grid::gpar(fill = fill, col = col, lwd = lwd * 96 / 72))
make_label <- function(fig_h) function(text, x, y, size, hjust = 0.5, col = "black", face = "bold")
  grid::textGrob(text, x = grid::unit(x, "in"), y = grid::unit(fig_h - y, "in"), hjust = hjust,
                 vjust = 0, gp = gp_text(size, col, face))

# a: three panels with headers; y-axis title only on the left panel
FIG_A_W <- 21.29; FIG_A_H <- 5.43
A_COL_X  <- c(1.0157, 8.0580, 15.0843)     # left edge of the panel area of each panel
A_ROW_Y  <- c(0.7600, 0.7693, 0.7480)      # top edge of the panel area of each panel
A_HEAD_X <- list(c(1.0032, 7.1474), c(8.0469, 14.1911), c(15.0718, 21.2161))
A_HEAD_Y <- c(0.0905, 0.7517)
A_HEAD_BASE  <- 0.5774
HEADER_SIZE_A <- 29.9
LETTER_SIZE   <- 29.8
STRIP_LINE <- "#7F7F7F"; STRIP_LWD <- 0.95

box_a <- make_box(FIG_A_H); label_a <- make_label(FIG_A_H)
g_full <- with_agg(ggplotGrob(seq_panel(SCALES[1])))
pcol <- unique(g_full$layout$l[g_full$layout$name == "panel"])
prow <- unique(g_full$layout$t[g_full$layout$name == "panel"])
AREA_W <- PANEL_W - in_w(g_full$widths[-pcol]); AREA_H <- PANEL_H - in_h(g_full$heights[-prow])
place_panel <- function(p, x, y, fig_h) {
  g <- with_agg(ggplotGrob(p + theme(plot.background = element_blank())))
  g$widths[pcol] <- grid::unit(AREA_W, "in"); g$heights[prow] <- grid::unit(AREA_H, "in")
  left <- in_w(g$widths[seq_len(pcol - 1)]); top <- in_h(g$heights[seq_len(prow - 1)])
  grid::grobTree(g, vp = grid::viewport(x = grid::unit(x - left, "in"), y = grid::unit(fig_h - y + top, "in"),
                                        width = sum(g$widths), height = sum(g$heights), just = c("left", "top")))
}
grobs_a <- list(box_a(0, 0, FIG_A_W, FIG_A_H, "white"),
                label_a("a", 0.0559, 0.2608, LETTER_SIZE, hjust = 0))
for (i in seq_along(SCALES)) {
  p <- seq_panel(SCALES[i], stars = TRUE)
  if (i > 1) p <- p + theme(axis.title.y = element_blank())
  hx <- A_HEAD_X[[i]]
  grobs_a <- c(grobs_a, list(place_panel(p, A_COL_X[i], A_ROW_Y[i], FIG_A_H),
                             box_a(hx[1], A_HEAD_Y[1], hx[2], A_HEAD_Y[2], STRIP_FILL, STRIP_LINE, STRIP_LWD),
                             label_a(SCALE_TITLES[[SCALES[i]]], mean(hx), A_HEAD_BASE, HEADER_SIZE_A)))
}
draw_tiff(grobs_a, file.path(out_dir, "Fig11a.tiff"), FIG_A_W, FIG_A_H)

# b: two rows of three panels under the headers of the two approaches
FIG_B_W <- 22.13; FIG_B_H <- 14.09
B_X <- list(DL = c(0.0115, 7.2868, 14.5648), SC = c(0.0115, 7.3034, 14.6118))   # left edge of each panel
B_Y <- list(DL = c(1.9282, 1.9213, 1.9144), SC = c(8.5077, 8.5105, 8.5022))     # top edge of each panel
B_BAND_Y <- list(DL = c(0.0083, 0.6197), SC = c(7.0124, 7.6210))
B_BAND_X <- c(0.0083, 22.1079)
B_BAND_BASE <- c(DL = 0.4689, SC = 7.4647)
B_SCALE_X <- c(2.8299, 9.9834, 17.4523); B_SCALE_BASE <- 1.0802
B_R2_X <- list(DL = c(0.2711, 7.5408, 14.8547), SC = c(0.2573, 7.5546, 14.9045))
B_R2_BASE <- list(DL = c(1.6653, 1.6625, 1.6708), SC = c(8.1909, 8.1909, 8.1909))
BAND_SIZE <- 29.9; SCALE_SIZE <- 28.0; R2_SIZE <- 24.2
SCALE_COL <- "#4472C4"; R2_COL <- "#7030A0"

box_b <- make_box(FIG_B_H); label_b <- make_label(FIG_B_H)
place_canvas <- function(p, x, y) {
  g <- with_agg(patchwork::patchworkGrob(p & theme(plot.background = element_blank())))
  grid::grobTree(g, vp = grid::viewport(x = grid::unit(x, "in"), y = grid::unit(FIG_B_H - y, "in"),
                                        width = grid::unit(PANEL_W, "in"), height = grid::unit(PANEL_H, "in"),
                                        just = c("left", "top")))
}
# R2m and R2c of the full model, with sub- and superscripts as in the published figure
r2_label <- function(R2m, R2c, x, y) {
  sz <- R2_SIZE; x0 <- grid::unit(x, "in"); y0 <- grid::unit(FIG_B_H - y, "in")
  parts <- list(list("R", 1, 0), list("2", 0.67, 0.30), list("m", 0.67, -0.25),
                list(sprintf("= %.2f, R", R2m), 1, 0), list("2", 0.67, 0.30), list("c ", 0.67, -0.25),
                list(sprintf("= %.2f", R2c), 1, 0))
  out <- list(); xx <- x0
  for (pt in parts) {
    tg <- grid::textGrob(pt[[1]], x = xx, y = y0 + grid::unit(pt[[3]] * sz, "pt"), hjust = 0, vjust = 0,
                         gp = gp_text(sz * pt[[2]], R2_COL))
    out[[length(out) + 1]] <- tg; xx <- xx + grid::grobWidth(tg)
  }
  grid::gTree(children = do.call(grid::gList, out))
}
r2_full <- seq_out[seq_out$step == "M4", ]
grobs_b <- list(box_b(0, 0, FIG_B_W, FIG_B_H, "white"))
for (path in PATHWAYS) {
  grobs_b <- c(grobs_b, list(box_b(B_BAND_X[1], B_BAND_Y[[path]][1], B_BAND_X[2], B_BAND_Y[[path]][2], STRIP_FILL),
                             label_b(PATH_LABELS[[path]], FIG_B_W / 2, B_BAND_BASE[[path]], BAND_SIZE)))
  for (i in seq_along(SCALES)) {
    rr <- r2_full[r2_full$Scale == SCALES[i] & r2_full$Pathway == path, ]
    grobs_b <- c(grobs_b, list(place_canvas(hp_panel(SCALES[i], path, numbers = TRUE), B_X[[path]][i], B_Y[[path]][i]),
                               r2_label(rr$R2m, rr$R2c, B_R2_X[[path]][i], B_R2_BASE[[path]][i])))
    if (path == "DL") grobs_b[[length(grobs_b) + 1]] <-
      label_b(SCALE_TITLES[[SCALES[i]]], B_SCALE_X[i], B_SCALE_BASE, SCALE_SIZE, col = SCALE_COL)
  }
}
grobs_b[[length(grobs_b) + 1]] <- label_b("b", 0.0719, 0.4371, LETTER_SIZE, hjust = 0)

# legend: the five groups
LEG_SQ_X  <- c(1.3887, 4.1217, 8.7469, 12.2739, 15.9696)    # left edge of each square
LEG_SQ_Y  <- c(13.6985, 13.7123, 13.7012, 13.7123, 13.6985) # top edge of each square
LEG_SQ    <- 0.3098                                         # side of the squares
LEG_TXT_X <- c(1.8174, 4.5372, 9.1922, 12.7546, 16.4292)
LEG_TXT_Y <- c(13.9917, 14.0000, 13.9945, 13.9917, 13.9779)
LEG_TXT   <- c("Management", "Compositional Dynamics", "Canopy Structure", "Intra-species SDiv",
               "Dominance of Top Species")
for (k in seq_along(LEG_TXT))
  grobs_b <- c(grobs_b, list(box_b(LEG_SQ_X[k], LEG_SQ_Y[k], LEG_SQ_X[k] + LEG_SQ, LEG_SQ_Y[k] + LEG_SQ,
                                   GROUP_COLS[[LEG_TXT[k]]]),
                             label_b(LEG_TXT[k], LEG_TXT_X[k], LEG_TXT_Y[k], R2_SIZE, hjust = 0)))
draw_tiff(grobs_b, file.path(out_dir, "Fig11b.tiff"), FIG_B_W, FIG_B_H)

cat("Fig. 11 written to", out_dir, "\n")
