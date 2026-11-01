# Beyond Spectral Diversity: Spectral detection of dominant species decodes grassland plant diversity across management regimes

Code and example data for the paper of the same title.

The repository follows the workflow of the paper in three parts. Each part has its own README with its input, output, usage and software requirements.

| Folder | Content | Language |
| --- | --- | --- |
| `1_classification/` | DK-CAT species classification framework with its trained weights; generates species maps from hyperspectral images | Python |
| `2_dominance_metrics/` | predicted species richness and the five dominance metrics derived from the species maps at three spatial scales | R |
| `3_statistics_figures/` | one script for each main-text figure (Fig. 4–12), combining the statistical analysis and the figure | R |
| `data/` | example data and simulated data (see `data/README.md`) | — |

## Data

- Hyperspectral images (`data/hyperspectral/`) and species maps (`data/species_map/`) of two example plots of the study, one under grazing and one under fencing. The images are the input of part 1, and the species maps the input of part 2.
- Simulated data (`data/simulated/`) with the sampling design of the study, the input of part 3. They demonstrate the analytical workflow and reproduce the layout of the figures.

## Running the workflow

    cd 1_classification
    python predict.py                       # species maps of the example plots
    cd ../2_dominance_metrics
    Rscript dominance_metrics.R             # dominance metrics from the species maps
    cd ../3_statistics_figures
    Rscript Fig04_direct_SR.R               # and the other figure scripts

## Requirements

Part 1: Python 3.9 with the packages listed in `1_classification/requirements.txt`. Parts 2 and 3: R 4.4.2 with the packages listed in their READMEs.

## License

Code: MIT License (see `LICENSE`). Data: CC BY 4.0.
