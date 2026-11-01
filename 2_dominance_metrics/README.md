# 2. Dominance metrics from the species maps

This folder contains the R script that derives predicted species richness (DL-SR) and the five dominance metrics (Top1, Top2, Top3, Top12 and Top123) from the DK-CAT species maps at the three spatial scales of the study (Section 3.1.5 of the paper).

## Input

One single-band GeoTIFF species map per 1.5 m × 1.5 m plot, as written by the inference code in `../1_classification`. The species maps of the two example plots are provided in `../data/species_map`.

- Each pixel is assigned a species. As GeoTIFF pixel values are numeric, each species is stored as an integer code; the codes and the species they stand for are listed in `../data/species_code.csv`.
- Pixels removed by the radiometric filtering (soil and shaded pixels) are stored as NoData and are excluded.
- File names start with `<unit>_<plot>`, where `<unit>` is the management unit and `<plot>` the plot within it, for example `49_15_species_map.tif`.

## Sampling units

Pixel size is 2 cm.

| Scale | Sampling unit | `Window` |
| --- | --- | --- |
| 0.5m | nine 0.5 m × 0.5 m quadrats laid in a 3 × 3 grid from the upper-left corner of the plot image | 1–9, row by row from the upper left |
| 1m | four overlapping 1 m × 1 m subplots anchored at the corners of the plot image | 1 upper left, 2 upper right, 3 lower left, 4 lower right |
| 1.5m | the whole plot image | 0 |

Each 0.5 m and 1 m window is defined by its geographic extent and cut from the plot image with `terra::crop(snap = "out")`. The species maps must therefore keep their original georeferencing.

## Metrics

- **Canopy share** of a species is the number of valid pixels assigned to it divided by the number of valid pixels in the sampling unit.
- **DL_SR** is the number of distinct species among the valid pixels.
- Species are ranked by canopy share. **DL_Top1**, **DL_Top2** and **DL_Top3** are the shares of the first-, second- and third-ranked species, **DL_Top12** = Top1 + Top2 and **DL_Top123** = Top1 + Top2 + Top3.
- A rank that is not occupied has a share of 0. With two species, Top3 = 0 and Top12 = Top123 = 1; with one species, Top2 = Top3 = 0 and Top1 = Top12 = Top123 = 1.

## Usage

Run from this folder:

```bash
Rscript dominance_metrics.R ../data/species_map output/dominance_metrics.csv
```

Both arguments are optional; the values shown are the defaults. In RStudio, set the working directory to this folder and source the script.

## Output

`output/dominance_metrics.csv`, one row per sampling unit:

| Column | Description |
| --- | --- |
| `Unit`, `Plot`, `Plot15` | management unit, plot within the unit, and `<unit>_<plot>` |
| `Scale` | `0.5m`, `1m` or `1.5m` |
| `Window` | sampling unit within the plot (see the table above) |
| `DL_SR` | predicted species richness |
| `DL_Top1` … `DL_Top123` | dominance metrics |

The same column names are used in the data tables read by the scripts in `../3_statistics_figures`.

## Requirements

R 4.4.2 with the *terra* package (version 1.8.5).
