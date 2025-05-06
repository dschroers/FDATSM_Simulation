# FDATSM Simulation

## Description

## Description

## Description

This repository contains the code for the simulation study presented in the article ["Dynamically Consistent Analysis of Realized Covariations in Term Structure Models"](https://arxiv.org/abs/2406.19412) by Dennis Schroers. The core methodology is implemented and documented in the [`FDATSM`](https://github.com/dschroers/FDATSM) R package.

The simulation study investigates the finite-sample performance of the truncated covariation estimator for difference returns, as outlined in the article. It aims to provide insights into the robustness of the estimator against jumps and its ability to recover the dimension of the original quadratic covariation.


## Files

- `Simulation Study.R`: The main R script for replicating the simulation study.
- `results.csv`: A CSV file containing the results of the simulation (500 Monte Carlo iterations).

## Installation

To use the code, you must have the following R packages installed:

- `MASS`
- `matrixcalc`
- `pracma`
- `npreg`
- `splines`


You can install these packages using:

```r
install.packages(c("MASS", "matrixcalc", "pracma", "npreg", "splines"))
```


## Usage

To run the simulation study, execute the Simulation Study.R script in RStudio or from the command line:

```r
source("Simulation Study.R")
```

The simulation involves 500 Monte Carlo iterations, which may take a significant amount of time. To reduce the simulation time it is possible to use an ex-post smoothing which will yield similar, but not exactly matching results.

The simulation results will be saved as a CSV file (results.csv) in the results folder.
