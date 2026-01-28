# Code for "Abattoir registrations of liver fluke in Danish cattle reveal transmission hotspots but no evidence of climatic effects"

Authors: Tiem van der Deure, Matthew Denwood, Stig Milan Thamsborg & Anna-Sofie Stensgaard

Article is in preparation.

This repository contains code and data for analyzing *Fasciola hepatica* (liver fluke) prevalence in Danish cattle using Bayesian geostatistical models.

Access to the Danish cattle database and slaughter database was provided for the purpose of this study by SEGES. Data is therefore not available.

## Repository Setup

### Environment Setup

This repository is structured as a Julia package. Dependencies are specified in `Project.toml` and `Manifest.toml`. To set up:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

### Directory Structure

- **`src/`** - Julia source code - defines functions used in the scripts
  - `FasciolaDK.jl` - Main module
  - `climate.jl` - Climate data download processing
  - `plots.jl` - Utilities used for figures
  - `registry.jl` - Reads in cattle registry data

- **`scripts/`** - Analysis and data processing scripts
  - `1_data_wrangling.jl` - Weather and cattle data is read in, filtered, and written to .csv files used by the other scripts
  - `2_inla.R` - Runs INLA models in R and writes resulting posteriors to .csv files
  - `3_figuredata.jl`  & `3_figures.jl` - Figure and Table generation
  - `4_estimate_cor.jl` - Supplementary analysis - estimates correlation between liver fluke and other liver disease within herds

- **`images/`** - Generated figures and tables


