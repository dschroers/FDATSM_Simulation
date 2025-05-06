# FDATSM Simulation Study

This repository contains the code and results for the Monte Carlo simulation study accompanying the paper "**[Your Paper Title]**", which introduces the `FDATSM` method.

## Description

The simulation study evaluates the performance of the `FDATSM` method under controlled settings. It consists of 500 Monte Carlo iterations, which took approximately 15 hours to run on macOS Big Sur with R 4.2.1.

The core method is implemented in a separate R package: [`FDATSM`](https://github.com/yourusername/FDATSM)  
(Replace this link with the actual URL of your package.)

## File Structure

- `Simulation_Study.R`: Main script including all functions, parameter setup, and simulation code.
- `results/simulation_results.csv`: Output from the simulation.
- `results/simulation_results.rds`: RDS version of the results (smaller, faster to load).
- `session_info.txt`: R session and package version info for reproducibility.

## Instructions

To rerun the simulation (not recommended unless necessary, as it takes ~15 hours):

```r
source("Simulation_Study.R")
```

To inspect the results without rerunning:

```r
sim_results <- read.csv("results/simulation_results.csv")
# or
sim_results <- readRDS("results/simulation_results.rds")
```

## Requirements

- **R version**: 4.2.1
- **Operating system**: macOS Big Sur 11.7.4
- **Key packages**:
  - `FDATSM` (your method package)
  - `pracma`
  - `matrixcalc`
  - `npreg`
  - `MASS`
  - ... (see `session_info.txt` for full list and versions)

## License

This code is released under the MIT License.  
See the `LICENSE` file for details.
