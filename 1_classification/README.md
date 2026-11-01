# 1. Species classification

This folder contains the implementation of the two-stage species classification framework of the paper (Section 3.1.2) and its trained weights. DK-CAT, a transformer with a knowledge-guided cross-attention module, extracts spectral features, and a CatBoost classifier assigns each pixel a species. The code generates species maps from hyperspectral images.

## Input

One multi-band GeoTIFF per 1.5 m × 1.5 m plot, processed as described in Section 2.3 of the paper: 215 bands from 450 to 895 nm, brightness-normalized, at 2 cm pixel size, with the pixels removed by the radiometric filtering (soil and shaded pixels) stored as NoData. The images of the two example plots are provided in `../data/hyperspectral`. File names end with `_hyperspectral.tif`, for example `49_15_hyperspectral.tif`.

## Model

| Component | Description |
| --- | --- |
| `model.py` | DK-CAT feature extraction and CatBoost classification |
| `predict.py` | applies the model to every image of the input folder and writes the species maps |
| `weights/encoder.pt` | trained DK-CAT weights |
| `weights/classifier.cbm` | trained CatBoost classifier |
| `weights/model.json` | band wavelengths, the 23 species in the order of the model outputs and their codes |

The knowledge stream of DK-CAT encodes the spectral regions of pigment absorption in the visible (450–500 nm), the red edge (700–750 nm) and the near-infrared (> 750 nm), together with spectral derivatives and statistical moments. Each pixel is assigned to the species with the highest predicted probability (Section 3.1.3).

## Usage

Run from this folder:

    python predict.py ../data/hyperspectral output

Both arguments are optional; the values shown are the defaults.

## Output

One single-band GeoTIFF per input image, `output/<unit>_<plot>_species_map.tif`, on the grid and coordinate system of the input. Each pixel is assigned a species, stored as an integer code; the codes and the species they stand for are listed in `../data/species_code.csv`. Pixels that are NoData in the input are NoData in the species map. For the two example plots, the output is identical to the species maps in `../data/species_map`, which are the input of `../2_dominance_metrics`.

## Requirements

Python 3.9 with the packages listed in `requirements.txt` (*numpy* 1.26.4, *scipy* 1.13.1, *torch* 2.5.1, *catboost* 1.2.8, *rasterio* 1.4.3):

    pip install -r requirements.txt
