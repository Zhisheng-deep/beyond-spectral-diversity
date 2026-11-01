# =============================================================================
# Settings and functions shared by the figure scripts in this folder
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
  library(ragg)
  library(dplyr)
})

# All figures are drawn with the ragg device, which renders text, lines and
# points identically on Windows, macOS and Linux. The same device measures text
# when panels are assembled.
cowplot::set_null_device("agg")

# ---- data -------------------------------------------------------------------------
DATA_DIR <- "../data/simulated"

# Read a data table and write management regimes as "Grazing" / "Fencing"
read_table <- function(file) {
  d <- read.csv(file, check.names = FALSE, fileEncoding = "UTF-8-BOM")
  d$Management <- ifelse(d$Management == "grazing", "Grazing",
                         ifelse(d$Management == "fencing", "Fencing", d$Management))
  d
}

# ---- figure style ------------------------------------------------------------------
MGMT_COLS    <- c(Grazing = "#E69F00", Fencing = "#56B4E9")
SCALES       <- c("0.5m", "1m", "1.5m")
SCALE_TITLES <- c("0.5m" = "0.5 m \u00d7 0.5 m", "1m" = "1 m \u00d7 1 m",
                  "1.5m" = "1.5 m \u00d7 1.5 m")
STRIP_FILL   <- "#F2F2F2"   # background of the panel headers

# Colour mixed with white (light tone for legend keys); amount = share of the colour
light_colour <- function(col, amount = 0.45) {
  m <- grDevices::col2rgb(col) / 255
  grDevices::rgb(t(1 - amount * (1 - m)))
}

theme_nature <- function(base_size) {
  theme_bw(base_size = base_size) +
    theme(
      text = element_text(family = "Arial", face = "bold"),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.25),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
      axis.line = element_line(color = "black", linewidth = 0.6),
      axis.ticks = element_line(color = "black", linewidth = 0.6),
      axis.text  = element_text(size = base_size - 2, face = "bold"),
      axis.title.x = element_text(margin = ggplot2::margin(t = 0.5), size = base_size, face = "bold"),
      axis.title.y = element_text(margin = ggplot2::margin(r = 1), size = base_size, face = "bold"),
      legend.position = "bottom",
      legend.key = element_blank(),
      legend.background = element_blank(),
      legend.text  = element_text(size = base_size, face = "bold"),
      legend.title = element_text(size = base_size, face = "bold"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = base_size + 4,
                                margin = ggplot2::margin(b = 2)),
      strip.text = element_text(size = base_size, face = "bold"),
      plot.margin = ggplot2::margin(1, 1, 1, 1)
    )
}

# Default axis breaks, keeping only whole numbers (counts such as species richness)
integer_breaks <- function(limits) {
  b <- scales::extended_breaks()(limits)
  b[abs(b - round(b)) < 1e-8]
}

# ---- statistics ------------------------------------------------------------------------
p_to_stars <- function(p) {
  if (is.na(p)) return("NS")
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  "NS"
}

# Adjusted R2 and P value of the slope of an ordinary least squares fit y ~ x
ols_fit <- function(y, x) {
  s <- summary(lm(y ~ x))
  c(n = length(y), adj_R2 = s$adj.r.squared, P = s$coefficients[2, 4])
}

# Stacked R2 labels in a panel corner. Labels with asterisks end at x; a label
# without asterisks (P >= 0.05) ends where its value lines up with the values
# of the other labels.
#   text    : label text up to the value, e.g. "Fencing: R\u00b2 = 0.38"
#   stars   : asterisks of each label ("" when not significant)
#   colour  : colour of each label
#   corner  : "bottom" or "top"; vjust: offsets of the labels in text heights
#   x       : right end of the labels as a fraction of the panel width
#   size_mm : text size (mm)
r2_labels_grob <- function(text, stars, colour, corner, vjust, x = 1, size_mm = 6.5) {
  gp <- function(col) grid::gpar(fontfamily = "Arial", fontface = "bold",
                                 fontsize = size_mm * .pt, col = col)
  suffix <- ifelse(stars == "", "", paste0(" (", stars, ")"))
  shift  <- if (all(suffix == "")) grid::unit(0, "pt") else
    do.call(max, lapply(suffix[suffix != ""], function(s) grid::grobWidth(grid::textGrob(s, gp = gp("black")))))
  y <- grid::unit(if (corner == "bottom") 0 else 1, "npc")
  grobs <- lapply(seq_along(text), function(i) {
    xi <- if (suffix[i] == "") grid::unit(x, "npc") - shift else grid::unit(x, "npc")
    grid::textGrob(paste0(text[i], suffix[i]), x = xi, y = y, hjust = 1,
                   vjust = vjust[i], gp = gp(colour[i]))
  })
  grid::gTree(children = do.call(grid::gList, grobs))
}

# ---- assembly of a multi-panel figure ---------------------------------------------------
# Grey header across the figure with the panel letter on the left and the
# title in the centre
panel_header <- function(letter, title, letter_size = 32, title_size = 28) {
  ggplot() +
    annotate("text", x = 0.004, y = 0.5, label = letter, hjust = 0, fontface = "bold",
             family = "Arial", size = letter_size / .pt) +
    annotate("text", x = 0.5, y = 0.5, label = title, fontface = "bold",
             family = "Arial", size = title_size / .pt) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    theme_void() +
    theme(plot.background = element_rect(fill = STRIP_FILL, colour = NA),
          plot.margin = ggplot2::margin(0, 0, 0, 0))
}

# Evaluate an expression that measures text (e.g. converts plots to grobs)
# with the ragg device, which also renders the figures
with_agg <- function(expr) {
  ragg::agg_png(tempfile(fileext = ".png"))
  on.exit(dev.off())
  expr
}

# Legend of a ggplot as a grob (works across ggplot2 versions)
get_legend_grob <- function(plot) {
  g   <- with_agg(ggplotGrob(plot))
  idx <- grep("^guide-box", g$layout$name)
  idx <- idx[!vapply(g$grobs[idx], inherits, logical(1), "zeroGrob")]
  g$grobs[[idx[1]]]
}

# Save a figure as LZW-compressed TIFF at 600 dpi
save_tiff <- function(plot, file, width, height) {
  ggsave(file, plot, device = ragg::agg_tiff, width = width, height = height,
         units = "in", dpi = 600, compression = "lzw")
}

# ---- collinearity screening of a candidate set of predictors ---------------------
# Removes, at each step, the predictor less correlated with the response (|Spearman
# rho|) until all pairwise |rho| <= cor_cutoff, then until all variance inflation
# factors <= vif_cutoff. Returns the retained predictors and the OLS model.
screen_fit <- function(df_sub, resp, pred_vars, cor_cutoff = 0.8, vif_cutoff = 5) {
  empty <- list(final_vars = character(0), model = NULL)
  vars_available <- pred_vars[pred_vars %in% names(df_sub)]
  vars_available <- vars_available[sapply(vars_available, function(v) !all(is.na(df_sub[[v]])))]
  if (nrow(df_sub) < 5 || length(vars_available) == 0) return(empty)
  df_use <- df_sub %>% dplyr::select(all_of(c(resp, vars_available))) %>% tidyr::drop_na()
  if (nrow(df_use) < 5) return(empty)
  y <- df_use[[resp]]; sel_vars <- vars_available
  sp <- function(a, b) stats::cor(a, b, method = "spearman", use = "pairwise.complete.obs")
  if (length(sel_vars) > 1) repeat {
    cor_abs <- abs(stats::cor(df_use[, sel_vars, drop = FALSE], method = "spearman", use = "pairwise.complete.obs"))
    cor_abs[lower.tri(cor_abs, diag = TRUE)] <- 0
    max_cor <- max(cor_abs, na.rm = TRUE)
    if (is.na(max_cor) || max_cor <= cor_cutoff || length(sel_vars) <= 1) break
    idx <- which(cor_abs == max_cor, arr.ind = TRUE)[1, ]
    v1 <- rownames(cor_abs)[idx[1]]; v2 <- colnames(cor_abs)[idx[2]]
    drop_var <- if (abs(sp(df_use[[v1]], y)) >= abs(sp(df_use[[v2]], y))) v2 else v1
    sel_vars <- setdiff(sel_vars, drop_var)
    if (length(sel_vars) <= 1) break
  }
  if (length(sel_vars) == 0) return(empty)
  df_model <- df_use[, c(resp, sel_vars), drop = FALSE]
  while (nrow(df_model) <= length(sel_vars) + 1 && length(sel_vars) > 1) {
    cr <- sapply(sel_vars, function(v) sp(df_model[[v]], df_model[[resp]]))
    sel_vars <- setdiff(sel_vars, names(which.min(abs(cr)))); df_model <- df_use[, c(resp, sel_vars), drop = FALSE]
  }
  if (nrow(df_model) <= length(sel_vars) + 1) return(empty)
  if (length(sel_vars) > 1) repeat {
    lm_now <- lm(as.formula(paste(resp, "~", paste(sel_vars, collapse = " + "))), data = df_model)
    ali <- alias(lm_now)$Complete
    if (!is.null(ali)) {
      av <- rownames(ali); cr <- sapply(av, function(v) sp(df_model[[v]], df_model[[resp]]))
      sel_vars <- setdiff(sel_vars, names(which.min(abs(cr)))); if (length(sel_vars) == 0) break
      df_model <- df_use[, c(resp, sel_vars), drop = FALSE]; if (length(sel_vars) == 1) break; next
    }
    vv <- car::vif(lm_now); if (is.matrix(vv)) vv <- diag(vv)
    if (is.na(max(vv, na.rm = TRUE)) || max(vv, na.rm = TRUE) <= vif_cutoff || length(sel_vars) == 1) break
    high <- names(vv)[vv > vif_cutoff]; cr <- sapply(sel_vars, function(v) sp(df_model[[v]], df_model[[resp]]))
    sel_vars <- setdiff(sel_vars, names(sort(abs(cr[high]), decreasing = FALSE))[1]); if (length(sel_vars) == 0) break
    df_model <- df_use[, c(resp, sel_vars), drop = FALSE]; if (length(sel_vars) == 1) break
  }
  if (length(sel_vars) == 0) return(empty)
  list(final_vars = sel_vars, model = lm(as.formula(paste(resp, "~", paste(sel_vars, collapse = " + "))), data = df_model))
}

# The five dominance metrics of a source ("DL", "SC" or "OB")
top_vars <- function(prefix) paste0(prefix, "_Top", c("1", "2", "12", "3", "123"))

# ---- cluster bootstrap of adjusted R2 ------------------------------------------------
# Adjusted R2 of an OLS fit with rows replicated w times (exact frequency weights)
adjR2_w <- function(X, y, w) {
  idx <- rep.int(seq_along(w), w); Xb <- X[idx, , drop = FALSE]; yb <- y[idx]
  n <- length(yb); p <- ncol(Xb)
  if (n <= p + 1 || length(unique(yb)) < 2) return(NA_real_)
  fit <- .lm.fit(Xb, yb); rss <- sum(fit$residuals^2); tss <- sum((yb - mean(yb))^2)
  if (tss == 0) return(NA_real_)
  1 - (rss / (n - p)) / (tss / (n - 1))
}

# Resampling unit: the 1.5 m plot (Plot15), stratified by management unit: within
# each unit, as many plots as the unit has are drawn with replacement. A drawn
# plot enters with all its rows (quadrats, subplots, windows). A row lists its
# plots in member_plots (";"-separated; pooled plots list several) and enters as
# many times as the product of its members' draw counts, i.e. not at all if a
# member was not drawn. The retained predictors of each model are held fixed and
# the coefficients refitted.
#   dat   : data with columns Scale, Unit, Plot15 and member_plots
#   specs : list of models, each list(rows, X, y) on the rows of dat
#   B     : number of replicates; seed: random seed
# Returns a B x length(specs) matrix of adjusted R2.
cluster_bootstrap <- function(dat, specs, B = 2000, seed = 20260908) {
  frame   <- dat %>% dplyr::filter(Scale == "1.5m") %>% dplyr::distinct(Unit, Plot15)
  units   <- split(frame$Plot15, as.character(frame$Unit))
  members <- strsplit(dat$member_plots, ";", fixed = TRUE)
  set.seed(seed)
  boot <- matrix(NA_real_, nrow = B, ncol = length(specs))
  for (b in seq_len(B)) {
    cnt <- integer(0)
    for (u in names(units)) {
      pl <- units[[u]]; dr <- sample(pl, length(pl), replace = TRUE)
      cnt[pl] <- tabulate(match(dr, pl), length(pl))
    }
    w_all <- vapply(members, function(m) prod(cnt[m]), numeric(1))
    for (i in seq_along(specs)) {
      s <- specs[[i]]; if (is.null(s)) next
      w <- w_all[s$rows]; if (sum(w > 0) < 5) next
      boot[b, i] <- adjR2_w(s$X, s$y, as.integer(w))
    }
  }
  boot
}

# 95% percentile interval (NA when fewer than 100 valid replicates)
boot_ci <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 100) return(c(lower = NA, upper = NA, n_valid = length(x)))
  c(lower = unname(quantile(x, 0.025)), upper = unname(quantile(x, 0.975)), n_valid = length(x))
}
