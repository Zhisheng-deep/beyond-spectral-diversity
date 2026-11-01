# =============================================================================
# Predicted species richness and dominance metrics from DK-CAT species maps
#
# Input : one single-band species map (GeoTIFF) per 1.5 m x 1.5 m plot, as
#         produced by the DK-CAT inference code. Each pixel is assigned a
#         species, stored as an integer code (codes and species are listed in
#         ../data/species_code.csv). Pixels removed by the radiometric
#         filtering (soil and shaded pixels) are stored as NoData and are
#         excluded.
#         File names start with "<unit>_<plot>", e.g. "49_15_species_map.tif".
# Output: one CSV row per sampling unit at three scales,
#           0.5m : nine 0.5 m x 0.5 m quadrats         (Window 1-9)
#           1m   : four overlapping 1 m x 1 m subplots (Window 1-4)
#           1.5m : the whole plot                     (Window 0)
#         with predicted species richness (DL_SR) and the five dominance
#         metrics (DL_Top1, DL_Top2, DL_Top3, DL_Top12, DL_Top123).
#
# Definitions
#   Canopy share of a species = number of valid pixels assigned to it divided
#   by the number of valid pixels in the sampling unit. Species are ranked by
#   canopy share; Top1, Top2 and Top3 are the shares of the first-, second-
#   and third-ranked species, Top12 = Top1 + Top2 and Top123 = Top1 + Top2 +
#   Top3. A rank that is not occupied (fewer than three species in the unit)
#   has a share of 0, e.g. with two species Top3 = 0 and Top12 = Top123 = 1.
#   DL_SR = number of distinct species among the valid pixels.
#
# Sampling units (pixel size 2 cm)
#   0.5m : 3 x 3 windows of 0.5 m x 0.5 m laid from the upper-left corner of
#          the plot image, numbered row by row (1-3 top row, 4-6 middle,
#          7-9 bottom).
#   1m   : 1 m x 1 m windows anchored at the four corners of the plot image,
#          1 = upper left, 2 = upper right, 3 = lower left, 4 = lower right.
#   1.5m : all pixels of the plot image.
#   Each 0.5 m and 1 m window is defined by its geographic extent and cut
#   from the plot image with terra::crop(snap = "out").
#
# Usage   : Rscript dominance_metrics.R [input_dir] [output_csv]
# Defaults: input_dir  = "../data/species_map"
#           output_csv = "output/dominance_metrics.csv"
# =============================================================================

suppressPackageStartupMessages(library(terra))

args       <- commandArgs(trailingOnly = TRUE)
input_dir  <- if (length(args) >= 1) args[1] else "../data/species_map"
output_csv <- if (length(args) >= 2) args[2] else "output/dominance_metrics.csv"

QUAD_PX <- 25L   # 0.5 m at 2 cm
SUB_PX  <- 50L   # 1 m at 2 cm
GRID_PX <- 75L   # 1.5 m at 2 cm

# ---- metrics for the pixels of one sampling unit -----------------------------
unit_metrics <- function(v) {
  v <- v[is.finite(v) & v != 0]
  if (length(v) == 0) {
    warning("A sampling unit has no valid pixels; its metrics are set to NA.")
    return(c(DL_SR = 0, DL_Top1 = NA, DL_Top2 = NA, DL_Top3 = NA,
             DL_Top12 = NA, DL_Top123 = NA))
  }
  share <- sort(as.numeric(table(v)) / length(v), decreasing = TRUE)
  top   <- c(share, 0, 0, 0)[1:3]              # unoccupied ranks have a share of 0
  c(DL_SR     = length(share),
    DL_Top1   = top[1],
    DL_Top2   = top[2],
    DL_Top3   = top[3],
    DL_Top12  = top[1] + top[2],
    DL_Top123 = top[1] + top[2] + top[3])
}

# ---- geographic extent of a window given by pixel rows r0-r1 and columns c0-c1
# (1-based; rows counted from the top, columns from the left of the plot image)
window_extent <- function(r, r0, r1, c0, c1) {
  ext(xmin(r) + (c0 - 1) * xres(r), xmin(r) + c1 * xres(r),
      ymax(r) - r1 * yres(r),       ymax(r) - (r0 - 1) * yres(r))
}

# ---- sampling units of one plot ---------------------------------------------
unit_windows <- function(nr, nc) {
  s    <- seq(1L, by = QUAD_PX, length.out = 3L)
  quad <- expand.grid(c0 = s, r0 = s)                    # row by row
  quad <- data.frame(Scale = "0.5m", Window = 1:9,
                     r0 = quad$r0, r1 = quad$r0 + QUAD_PX - 1L,
                     c0 = quad$c0, c1 = quad$c0 + QUAD_PX - 1L)
  sub  <- data.frame(Scale = "1m", Window = 1:4,
                     r0 = c(1L, 1L, nr - SUB_PX + 1L, nr - SUB_PX + 1L),
                     c0 = c(1L, nc - SUB_PX + 1L, 1L, nc - SUB_PX + 1L))
  sub$r1 <- sub$r0 + SUB_PX - 1L
  sub$c1 <- sub$c0 + SUB_PX - 1L
  plot <- data.frame(Scale = "1.5m", Window = 0L, r0 = 1L, r1 = nr, c0 = 1L, c1 = nc)
  rbind(quad, sub[, names(quad)], plot)
}

# ---- main --------------------------------------------------------------------
files <- list.files(input_dir, pattern = "\\.tif$", full.names = TRUE, ignore.case = TRUE)
if (length(files) == 0) stop("No .tif files found in ", input_dir)

rows <- list()
for (f in files) {
  id <- regmatches(basename(f), regexec("^([0-9]+)_([0-9]+)", basename(f)))[[1]]
  if (length(id) == 0) {
    warning("File name does not start with <unit>_<plot>; skipped: ", basename(f))
    next
  }
  r  <- rast(f)[[1]]
  nr <- nrow(r); nc <- ncol(r)
  if (nr < GRID_PX || nc < GRID_PX) {
    warning("Plot image smaller than 75 x 75 pixels; skipped: ", basename(f))
    next
  }
  w <- unit_windows(nr, nc)
  for (k in seq_len(nrow(w))) {
    unit <- if (w$Scale[k] == "1.5m") r else
      crop(r, window_extent(r, w$r0[k], w$r1[k], w$c0[k], w$c1[k]), snap = "out")
    met <- unit_metrics(values(unit, mat = FALSE))
    rows[[length(rows) + 1L]] <- data.frame(
      Unit = as.integer(id[2]), Plot = as.integer(id[3]),
      Plot15 = paste(id[2], id[3], sep = "_"),
      Scale = w$Scale[k], Window = w$Window[k],
      t(met), check.names = FALSE)
  }
}
if (length(rows) == 0) stop("No plot was processed.")

res <- do.call(rbind, rows)
res$DL_SR <- as.integer(res$DL_SR)
res <- res[order(match(res$Scale, c("0.5m", "1m", "1.5m")), res$Unit, res$Plot, res$Window), ]

dir.create(dirname(output_csv), showWarnings = FALSE, recursive = TRUE)
write.csv(res, output_csv, row.names = FALSE)
cat("Processed", length(unique(res$Plot15)), "plot(s);", nrow(res), "rows written to", output_csv, "\n")
