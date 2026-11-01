"""Species maps from the 215-band hyperspectral images of the plots.

Usage  : python predict.py [input_dir] [output_dir] [--device cpu]
Default: input_dir  = ../data/hyperspectral
         output_dir = output
Every file <unit>_<plot>_hyperspectral.tif in input_dir is written to
output_dir as <unit>_<plot>_species_map.tif.
"""

import argparse
from pathlib import Path

import numpy as np
import rasterio

from model import SpeciesModel

HERE = Path(__file__).resolve().parent


def species_map(model, input_file, output_file):
    with rasterio.open(input_file) as src:
        if src.count != model.num_bands:
            raise ValueError(f"{input_file.name}: expected {model.num_bands} bands, found {src.count}")
        cube = src.read(out_dtype="float32")
        profile = src.profile.copy()
    # pixels removed by the radiometric filtering are NoData in the input
    valid = np.isfinite(cube).all(axis=0)
    out = np.full(valid.shape, np.nan, dtype=np.float32)
    if valid.any():
        out[valid] = model.predict(cube[:, valid].T)
    profile.update(count=1, dtype="float32", nodata=np.nan, interleave="band")
    with rasterio.open(output_file, "w", **profile) as dst:
        dst.write(out, 1)
    return int(valid.sum())


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("input_dir", nargs="?", default=str(HERE.parent / "data" / "hyperspectral"))
    parser.add_argument("output_dir", nargs="?", default="output")
    parser.add_argument("--device", default="cpu", help="cpu or cuda")
    args = parser.parse_args()

    model = SpeciesModel(HERE / "weights", device=args.device)
    files = sorted(Path(args.input_dir).glob("*_hyperspectral.tif"))
    if not files:
        raise SystemExit(f"No *_hyperspectral.tif files found in {args.input_dir}")
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    for f in files:
        target = out_dir / f.name.replace("_hyperspectral.tif", "_species_map.tif")
        n = species_map(model, f, target)
        print(f"{f.name}: {n} vegetation pixels classified -> {target}")


if __name__ == "__main__":
    main()
