# Microclimate connectivity scripts

This repository contains R scripts used to calculate microclimate connectivity in the Amazon for the manuscript "Amazon microclimate connectivity collapses under high-emission scenarios". The scripts assume that the required public climate, topographic, and protected area raster/vector files have already been prepared following the processing workflow described in the manuscript. No local machine-specific paths are used.

## Repository structure

```
amazon-microclimate-connectivity/
└── Amazon_microclimate_connectivity/
    ├── data/
    │   ├── example_data/          # Small spatial subset for testing
    │   │   ├── micro_t_test_data.tif
    │   │   ├── micro_245_test_data.tif
    │   │   ├── micro_585_test_data.tif
    │   │   └── river_test_data.gpkg
    │   └── outputs/               # Example outputs
    │       ├── connectivity_polygons.csv
    │       └── r_roles.tif
    └── scripts/
        ├── 01_calculate_climate_connectivity_amazon.R
        ├── 02_climate_connectivity_statistics.R
        ├── 03_patch_role_analysis.R
        └── 04_sensitivity_analysis.R
```

## Scripts

| Script | Description |
| --- | --- |
| `01_calculate_climate_connectivity_amazon.R` | Calculates microclimate (or macroclimate) connectivity following McGuire et al. (2016). |
| `02_climate_connectivity_statistics.R` | Statistical analyses for microclimate connectivity. |
| `03_patch_role_analysis.R` | Assigns a role to each temperature patch (origin, neutral, intermediate, destination) based on the diversity of simulated species movements.&#32; |
| `04_sensitivity_analysis.R` | Evaluates sensitivity of microclimate connectivity to different warming magnitudes between present conditions and the SSP5-8.5 scenario. |

## System requirements

All analyses are implemented in **R**. The following R packages are required:

| Package | Version tested | Used in |
| --- | --- | --- |
| terra | ≥ 1.7 | Scripts 01, 02, 03, 04 |
| sf | ≥ 1.0 | Scripts 01, 03, 04 |
| igraph | ≥ 1.5 | Script 01 |
| dplyr | ≥ 1.1 | Scripts 01, 02, 03, 04 |
| mgcv | ≥ 1.9 | Script 02 |
| tidyr | ≥ 1.3 | Script 02 |
| effectsize | ≥ 0.8 | Script 02 |
| ape | ≥ 5.7 | Script 02 |
| data.table | ≥ 1.14 | Script 03 |
| tidyverse | ≥ 2.0 | Script 03 |
| ggplot2 | ≥ 3.4 | Scripts 03, 04 |
| patchwork | ≥ 1.2 | Script 03 |

The scripts were developed and tested with R 4.5.2 on a Windows 10/11 desktop environment. They should also run on macOS and Linux.

No specialized hardware is required. For the full Amazon-wide analysis, a machine with at least 16 GB RAM is recommended due to the size of the rasters. The example/demo dataset runs on any standard desktop or laptop computer.

## Installation

1. Clone or download this repository:

   ```
   git clone https://github.com/irisstarck/amazon-microclimate-connectivity.git
   ```
2. Install required R packages by running the following in R:

   ```
   install.packages(c(
     "terra", "sf", "igraph", "dplyr",
     "mgcv", "tidyr", "effectsize", "ape",
     "data.table", "tidyverse", "ggplot2", "patchwork"
   ))
   ```

Typical installation time on a standard desktop computer is approximately 5–15 minutes, depending on internet speed and whether geospatial dependencies are already available.

## Demo

A small spatial subset of the study area is included in `Amazon_microclimate_connectivity/data/example_data/` to test the climate connectivity workflow without access to the full Amazon-wide datasets.

### Input files

| File | Description |
| --- | --- |
| `micro_t_test_data.tif` | Present-day microclimate temperature raster (2013–2022 mean) |
| `micro_245_test_data.tif` | Future microclimate temperature raster (SSP2-4.5, 2093–2100 mean) |
| `micro_585_test_data.tif` | Future microclimate temperature raster (SSP5-8.5, 2093–2100 mean) |
| `river_test_data.gpkg` | River polygon layer (subset of GRWL) |

### Instructions to run the demo

1. Open `Amazon_microclimate_connectivity/scripts/01_calculate_climate_connectivity_amazon.R` in R or RStudio.
2. Set the `data.dir` and `result.folder` path variables at the top of the script to point to the `example_data/` and `outputs/` directories, respectively.
3. Run the script from top to bottom. Note: all lines that save outputs are initially commented out.

### Expected output

Running the demo produces:

- Temperature patch polygons (vector layer)
- Microclimate connectivity rasters for SSP2-4.5 and SSP5-8.5
- A connectivity summary CSV (`connectivity_polygons.csv`)
- A connectivity shapefile (`connectivity_polygons.shp`)

Results from the example data will differ from those reported in the manuscript, which used the complete Amazon-wide datasets.

## Running the scripts

All scripts are located in `Amazon_microclimate_connectivity/scripts/`. Scripts should be run in numerical order.

1. `01_calculate_climate_connectivity_amazon.R`&#32;
2. `02_climate_connectivity_statistics.R`&#32;
3. `03_patch_role_analysis.R`&#32;
4. `04_sensitivity_analysis.R`&#32;

**Before running any script:**

1. Set the `data.dir` and `result.folder` path variables to match your local directory structure.
2. Ensure all raster and vector inputs are projected to EPSG:5880, cover the same spatial extent, and have matching resolution and alignment.
3. Uncomment output saving lines to write results to disk.

Large raster datasets and restricted field-derived intermediate products are not included directly in this repository because of file size and data-sharing constraints. Publicly available input datasets are described in the manuscript Data Availability section.

## Expected run time

Run times for the full dataset will vary with available memory, raster sizes and temperature class number. Creating temperature patches in script 01 can take several hours. The other scripts take approximately 30-60 min.
