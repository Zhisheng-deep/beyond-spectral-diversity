# =============================================================================
# Fig. 8  Spearman correlations of observed species richness (outer three
#         rings) and observed Shannon-Wiener diversity (inner three rings) with
#         the predicted diversity of the deep learning species maps (DL) and of
#         the spectral clusters (SC), the dominance metrics (DL-Top*, SC-Top*)
#         and three conventional spectral diversity indices (CV, CHA, Rao's Q),
#         at three spatial scales.
#
# Steps
#   1. For each ring (response x scale) and each of the 15 predictors, compute
#      Spearman's rho and its P value; adjust the 15 P values of each ring with
#      the Benjamini-Hochberg false discovery rate. For the DL and SC columns
#      the predictor is predicted species richness for the species-richness
#      rings and predicted Shannon-Wiener diversity for the Shannon-Wiener
#      rings.
#   2. Draw the fan: one cell per ring and predictor, filled by rho, with rho
#      printed in the cell and asterisks from the adjusted P value
#      (* <= 0.05, ** <= 0.01, *** <= 0.001).
#   3. Add the ring symbols, the colour bar and the symbol legend, and export
#      the figure.
#
# Input  : analysis_table.csv (columns Scale, Observed_SR, Observed_Shannon,
#          DL_SR, SSR_km_z, DL_Shannon, SC_Shannon, DL_Top*, SC_Top*, CV, CHA,
#          Rao'S Q)
# Output : output/Fig08_statistics.csv   rho, P and adjusted P of every pair
#          output/Fig08.tiff             figure (6.53 x 4.24 in, 600 dpi)
#
# Usage  : Rscript Fig08_correlation_fan.R [data_file] [output_dir]
# Default: data_file  = "../data/simulated/analysis_table.csv"
#          output_dir = "output"
# =============================================================================

source("functions.R")

args      <- commandArgs(trailingOnly = TRUE)
data_file <- if (length(args) >= 1) args[1] else file.path(DATA_DIR, "analysis_table.csv")
out_dir   <- if (length(args) >= 2) args[2] else "output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

dat <- read_table(data_file)

# predictors from left to right; column used for the species-richness rings
# (x_SR) and for the Shannon-Wiener rings (x_SW)
FEATURES <- data.frame(
  label = c("DL", "SC",
            "DL-Top1", "DL-Top2", "DL-Top12", "DL-Top3", "DL-Top123",
            "SC-Top1", "SC-Top2", "SC-Top12", "SC-Top3", "SC-Top123",
            "CV", "CHA", "Rao's Q"),
  x_SR  = c("DL_SR", "SSR_km_z",
            "DL_Top1", "DL_Top2", "DL_Top12", "DL_Top3", "DL_Top123",
            "SC_Top1", "SC_Top2", "SC_Top12", "SC_Top3", "SC_Top123",
            "CV", "CHA", "Rao'S Q"),
  stringsAsFactors = FALSE)
FEATURES$x_SW <- FEATURES$x_SR
FEATURES$x_SW[1:2] <- c("DL_Shannon", "SC_Shannon")

# rings from outside to inside
RINGS <- data.frame(Response = rep(c("Observed_SR", "Observed_Shannon"), each = 3),
                    Scale    = rep(SCALES, 2),
                    stringsAsFactors = FALSE)

stars_fdr <- function(q) {
  if (is.na(q)) return("")
  if (q <= 0.001) return("***")
  if (q <= 0.01)  return("**")
  if (q <= 0.05)  return("*")
  ""
}

# ---- 1. statistics ----------------------------------------------------------------
res <- do.call(rbind, lapply(seq_len(nrow(RINGS)), function(r) {
  d    <- dat[dat$Scale == RINGS$Scale[r], ]
  xcol <- if (RINGS$Response[r] == "Observed_SR") FEATURES$x_SR else FEATURES$x_SW
  out  <- do.call(rbind, lapply(seq_len(nrow(FEATURES)), function(k) {
    x  <- d[[xcol[k]]]
    y  <- d[[RINGS$Response[r]]]
    ok <- complete.cases(x, y)
    if (sum(ok) < 3) {
      rho <- NA_real_; p <- 1
    } else {
      ct  <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman"))
      rho <- unname(ct$estimate); p <- ct$p.value
    }
    data.frame(Ring = r, Response = RINGS$Response[r], Scale = RINGS$Scale[r],
               Feature = k, Predictor = FEATURES$label[k], Column = xcol[k],
               n = sum(ok), rho = rho, P = p, stringsAsFactors = FALSE)
  }))
  out$P_FDR <- p.adjust(out$P, method = "BH")          # within each ring
  out
}))
res$Stars <- vapply(res$P_FDR, stars_fdr, character(1))
write.csv(res, file.path(out_dir, "Fig08_statistics.csv"), row.names = FALSE)

# ---- 2. fan -------------------------------------------------------------------------
# Plot units: a 7.6-in-wide plot with x limits +/-14.2 (and R's 4% axis extension)
UNIT_IN <- 7.6 / (2 * 14.2 * 1.08)                      # inches per plot unit
FIG_W   <- 6.53; FIG_H <- 4.24                          # inches
X_LEFT  <- -13.14; Y_TOP <- 13.34                       # plot units at the top-left corner

pal <- colorRampPalette(c("#588157", "#FFFFFF", "#F9A620"))(200)
rho_to_col <- function(x) {
  if (is.na(x)) return(NA_character_)
  pal[max(1, min(200, round((x + 1) / 2 * 199) + 1))]
}

N_feat  <- nrow(FEATURES)
N_ring  <- nrow(RINGS)
fan_rad <- 160 * pi / 180
theta_centers <- seq(pi / 2 + fan_rad / 2, pi / 2 - fan_rad / 2, length.out = N_feat)
half_span     <- fan_rad / N_feat / 2
bar_height    <- 1.10
ring_gap      <- 0.12
radii         <- 11.8 - (0:(N_ring - 1)) * (bar_height + ring_gap)
label_radius  <- max(radii) + 1.25

ragg::agg_tiff(file.path(out_dir, "Fig08.tiff"), width = FIG_W, height = FIG_H,
               units = "in", res = 600, compression = "lzw")
par(family = "Arial", font = 2, mar = c(0, 0, 0, 0))
plot.new()
plot.window(xlim = c(X_LEFT, X_LEFT + FIG_W / UNIT_IN),
            ylim = c(Y_TOP - FIG_H / UNIT_IN, Y_TOP), xaxs = "i", yaxs = "i")

for (i in seq_len(nrow(res))) {
  k   <- res$Feature[i]
  r_c <- radii[res$Ring[i]]
  th  <- theta_centers[k]
  ang <- c(th - half_span, th - half_span, th + half_span, th + half_span)
  rs  <- c(r_c - bar_height / 2, r_c + bar_height / 2, r_c + bar_height / 2, r_c - bar_height / 2)
  if (is.na(res$rho[i])) next
  polygon(rs * cos(ang), rs * sin(ang), col = rho_to_col(res$rho[i]), border = "white", lwd = 0.7)
  xt  <- r_c * cos(th); yt <- r_c * sin(th)
  srt <- th * 180 / pi - 90
  lab <- sprintf("%.2f", res$rho[i])
  if (lab == "-0.00") lab <- "0.00"
  text(xt, yt, lab, cex = 0.65, srt = srt, xpd = TRUE, font = 2)
  if (res$Stars[i] != "") {                    # outwards by 0.35 and to the left by 0.27
    text(xt + 0.35 * cos(th) + 0.27 * sin(th), yt + 0.35 * sin(th) - 0.27 * cos(th),
         res$Stars[i], cex = 0.55, srt = srt, xpd = TRUE, font = 2)
  }
}
for (k in seq_len(N_feat)) {
  th <- theta_centers[k]
  text(label_radius * cos(th), label_radius * sin(th), FEATURES$label[k],
       cex = 0.65, srt = th * 180 / pi - 90, xpd = TRUE, font = 2)
}

# ---- 3. ring symbols, colour bar and symbol legend -----------------------------------
# black symbols; (x, y) is the centre of the bounding box, w its width
draw_symbol <- function(type, x, y, w) {
  ring <- function(n, r, rot) { a <- rot + 2 * pi * (0:(n - 1)) / n; cbind(r * cos(a), r * sin(a)) }
  v <- switch(type,
    triangle = ring(3, 1, pi / 2),
    diamond  = ring(4, 1, pi / 2),
    star     = { o <- ring(5, 1, pi / 2); i <- ring(5, 0.382, pi / 2 + pi / 5)
                 rbind(o, i)[c(rbind(1:5, 6:10)), ] },
    pentagon = ring(5, 1, pi / 2),
    hexagon  = ring(6, 1, 0),
    circle   = ring(100, 1, 0))
  h  <- if (type == "star") w else w * diff(range(v[, 2])) / diff(range(v[, 1]))
  v[, 1] <- (v[, 1] - mean(range(v[, 1]))) / diff(range(v[, 1])) * w + x
  v[, 2] <- (v[, 2] - mean(range(v[, 2]))) / diff(range(v[, 2])) * h + y
  polygon(v[, 1], v[, 2], col = "black", border = NA)
}
SYMBOLS   <- c("triangle", "diamond", "star", "pentagon", "hexagon", "circle")   # rings 1-6
SYMBOL_W  <- c(triangle = 0.56, diamond = 0.554, star = 0.58,
               pentagon = 0.577, hexagon = 0.577, circle = 0.478)              # plot units
TEXT_CEX  <- 0.75                                     # colour-bar title and legend text

# ring symbols: below the lower edge of the DL column, 0.48 units beyond the
# edge and 0.044 units outwards along the radius
th_edge <- theta_centers[1] + half_span
for (r in seq_len(N_ring)) {
  rr <- radii[r] + 0.044
  draw_symbol(SYMBOLS[r], rr * cos(th_edge) - 0.48 * sin(th_edge),
              rr * sin(th_edge) + 0.48 * cos(th_edge), SYMBOL_W[SYMBOLS[r]])
}

# colour bar: drawn as in a 6.2 x 1.2 in legend plot (x limits +/-1.1, y limits
# 0.6-1.05, margins 2.8/2/0.8/2 lines), reduced to 0.92 of that size
LEG_SCALE <- 0.92
LEG_X_IN  <- 5.4  / (2.2  * 1.08)                     # inches per legend x unit
LEG_Y_IN  <- 0.48 / (0.45 * 1.08)                     # inches per legend y unit
BAR_LEFT  <- -8.369; BAR_BOTTOM <- -1.128             # plot units
lx <- function(x) BAR_LEFT + (x + 1) * LEG_X_IN * LEG_SCALE / UNIT_IN
ly <- function(y) BAR_BOTTOM + (y - 0.86) * LEG_Y_IN * LEG_SCALE / UNIT_IN
rasterImage(matrix(pal, nrow = 1), lx(-1), ly(0.86), lx(1), ly(1.01),
            interpolate = FALSE)                       # 200 colour steps, no seams
rect(lx(-1), ly(0.86), lx(1), ly(1.01), border = "black", lwd = 1.3 * LEG_SCALE)
ticks <- seq(-1, 1, by = 0.25)
for (x in ticks) {
  len <- if (x %in% c(-1, 0, 1)) 0.10 else 0.07
  segments(lx(x), ly(0.86), lx(x), ly(0.86 - len), lwd = 1.2 * LEG_SCALE)
}
text(lx(ticks), ly(0.86 - 0.14), labels = ticks, cex = 0.8 * LEG_SCALE, font = 2)
text(lx(0), ly(1.01) + 0.51, "Spearman Correlation Coefficient", cex = TEXT_CEX, font = 2)

# symbol legend
LEG_ROW_Y   <- c(-2.561, -3.464)                        # row centres (plot units)
LEG_HEAD_X  <- -6.78                                    # row titles, left aligned
LEG_SYM_X   <- c(-2.17, 1.22, 3.77)                     # symbol centres
LEG_LAB_X   <- c(-1.77, 1.60, 4.17)                     # labels, left aligned
LEG_HEADS   <- c("Species Richness", "Shannon-Wiener")
LEG_LABELS  <- paste0(c("0.5m", "1m", "1.5m"), "\u00d7", c("0.5m", "1m", "1.5m"))
for (row in 1:2) {
  text(LEG_HEAD_X, LEG_ROW_Y[row], LEG_HEADS[row], adj = c(0, 0.5), cex = TEXT_CEX, font = 2)
  for (j in 1:3) {
    s <- SYMBOLS[(row - 1) * 3 + j]
    draw_symbol(s, LEG_SYM_X[j], LEG_ROW_Y[row], SYMBOL_W[s])
    text(LEG_LAB_X[j], LEG_ROW_Y[row], LEG_LABELS[j], adj = c(0, 0.5), cex = TEXT_CEX, font = 2)
  }
}
invisible(dev.off())
